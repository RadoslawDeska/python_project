"""
Headless-compatible Integration and Fitting classes.
No GUI dependencies - can be used anywhere.
"""

from math import factorial
import numpy as np
from numpy.typing import NDArray
from lmfit import Minimizer, Parameters
from lmfit.minimizer import MinimizerResult
from scipy.special import hyp2f1
from typing import Optional, Tuple


class Integration:
    """
    Integrates the electric field according to Sheik-Bahae procedure.
    Completely decoupled from GUI - all parameters passed explicitly.
    """

    def __init__(
        self,
        beta: float,
        n2: float,
        DPhi0: float,
        positions: NDArray,
        d0: float,
        aperture_radius: float,
        wavelength: float,
        beamwaist: float,
        n_components: int,
        integration_steps: int,
        stype: str = "CA",
    ):
        # Data range
        self.z = positions  # evenly spaced

        # Apertures and distances
        self.d0 = d0  # distance from z=0 to aperture plane [m]
        self.ra = aperture_radius  # aperture radius [m]

        # Beam properties
        self.lda = wavelength  # wavelength [m]
        self.w0 = beamwaist  # beam waist [m]

        # Sample properties
        self.n2 = n2  # non-linear refractive index [m²/W]

        # Calculate transmittance properly
        try:
            self.T = beta * self.lda / self.n2 if self.n2 != 0 else 0
        except (ZeroDivisionError, TypeError):
            self.T = 0

        self.DPhi0 = DPhi0  # on-axis phase shift [rad]

        # Integration parameters
        self.mm = n_components
        self.ir = integration_steps
        self.stype = stype

        # Derived values (calculated in derive)
        self.closed_sum: NDArray[np.complex128] = np.zeros_like(self.z, dtype=np.complex128)
        self.open_sum: NDArray[np.complex128] = np.zeros_like(self.z, dtype=np.complex128)
        self.Tznorm: np.ndarray = np.zeros_like(self.z, dtype=np.complex128)
        self.Tz: np.ndarray = np.zeros_like(self.z, dtype=np.complex128)
        self.z0: float = 0.0

        self.E: NDArray[np.complex128] = np.zeros_like(self.z, dtype=np.complex128)

        # Internal beam decomposition variables (initialized for mypy)
        self.dr: float = 0.0
        self.wa: float = 0.0
        self.wm0: NDArray = np.array([])
        self.dm: NDArray = np.array([])
        self.wm: NDArray = np.array([])
        self.tm: NDArray = np.array([])
        self.Rm: NDArray = np.array([])
        self.wz: NDArray[np.float64] = np.zeros_like(self.z)
        self.Rz: NDArray[np.float64] = np.zeros_like(self.z)
        self.d: NDArray[np.float64] = np.zeros_like(self.z)
        self.g: NDArray[np.float64] = np.zeros_like(self.z)
        
        self.derive(self.DPhi0, self.w0, self.d0, self.ra, stype)

    def derive(
        self, DPhi0: float, w0: float, d0: float, ra: float, stype: str
    ) -> None:
        """Calculate derived parameters and perform integration."""
        # Beam properties
        self.k: float = 2 * np.pi / self.lda  # wave vector
        self.z0 = 0.5 * self.k * w0**2  # Rayleigh range
        self.wa = w0 * np.sqrt(
            1 + d0**2 / self.z0**2
        )  # beam radius at aperture

        # Aperture radius
        self.ra = ra

        # Sample properties
        self.Dphi0 = DPhi0 / (1 + self.z**2 / self.z0**2)

        # Additional derived parameters
        self.wz = w0 * np.sqrt(1 + self.z**2 / self.z0**2)
        self.Rz = self.z + self.z0**2 / self.z
        self.d = d0 - self.z
        self.g = 1 + self.d / self.Rz

        if stype == "CA":
            self.bigproduct()
        elif stype == "OA":
            self.calculate_Tz_for_OA(model="2PA")

    def calculate_Tz_for_OA(self, model: str = "2PA") -> NDArray:
        """Calculate transmittance for open aperture measurement."""
        self.Tz = np.zeros_like(self.z, dtype=float)

        if self.T != 0:
            Psi1 = self.T
            Psi2 = self.T * 8
        else:
            Psi1 = 0
            Psi2 = 0

        if model == "2PA":
            psi1 = Psi1 / (1 + self.z**2 / self.z0**2)
            self.Tz = hyp2f1(1, 1, 2, -psi1)
        elif model == "3PA":
            psi2 = Psi2 / (1 + self.z**2 / self.z0**2)
            self.Tz = hyp2f1(1 / 2, 1 / 2, 3 / 2, -((psi2) ** 2))
        elif model == "2PA+3PA":
            if Psi1 != 0:
                psi1 = Psi1 / (1 + self.z**2 / self.z0**2)
                psi2 = Psi2 / (1 + self.z**2 / self.z0**2)
                f_coupling = 1 + psi1 * (
                    0.339 * np.sin(0.498 * psi2) - 0.029
                ) / (1 + 0.966 * psi1 * psi2**-0.718)
                self.Tz = (
                    hyp2f1(1, 1, 2, -psi1)
                    * hyp2f1(1 / 2, 1 / 2, 3 / 2, -((psi2) ** 2))
                    * f_coupling
                )
            else:
                self.Tz = np.ones_like(self.z)
        else:  # RSA, SA, etc.
            self.Tz = np.ones_like(self.z)

        # Normalize
        self.Tznorm = (
            2
            * self.Tz
            / (np.average(self.Tz[0:10]) + np.average(self.Tz[-10:]))
        )
        return self.Tznorm

    def calculate_fm(self):
        """Calculate Gaussian decomposition coefficients."""
        self.fm = [
            (1j * self.Dphi0) ** m / factorial(m) * self.product[m]
            for m in range(self.mm)
        ]
        return self.fm

    def bigproduct(self) -> None:
        """Calculate big product for CA and get normalized transmittance."""
        self.product = []
        for m in range(self.mm):
            if m == 0:
                self.product.append(1)
            else:
                prod = np.cumprod(
                    [
                        1 + 1j * (j - 1 / 2) / (2 * np.pi) * self.T
                        for j in range(1, m + 1)
                    ]
                )[-1]
                self.product.append(prod)

        self.fm = self.calculate_fm()
        Tzo = self.open()
        Tzc = self.closed()
        self.Tznorm = Tzc / Tzo

    def open(self) -> NDArray[np.complex128]:
        """Calculate open aperture transmittance."""
        self.dr = 3 * self.wa / self.ir
        self.open_sum = self.bigsum()
        return self.open_sum

    def closed(self) -> NDArray[np.complex128]:
        """Calculate closed aperture transmittance."""
        self.dr = self.ra / self.ir
        self.closed_sum = self.bigsum()
        return self.closed_sum

    def bigsum(self) -> NDArray:
        """Perform radial integration for transmittance calculation."""
        self.Tz = np.zeros_like(self.z, dtype=np.complex128)

        for rr in range(self.ir):
            self.E = np.zeros_like(self.z, dtype=np.complex128)

            for m in range(self.mm):
                self.wm0 = self.wz / np.sqrt(2 * m + 1)
                self.dm = 0.5 * self.k * self.wm0**2
                self.wm = self.wm0 * np.sqrt(self.g**2 + self.d**2 / self.dm**2)
                self.tm = np.arctan(self.g / (self.d / self.dm))
                self.Rm = self.d / (
                    1 - self.g / (self.g**2 + self.d**2 / self.dm**2)
                )

                self.E += (
                    self.fm[m]
                    * np.exp(1j * self.tm)
                    * self.wm0
                    / self.wm
                    / self.wz
                    * np.exp(
                        (-1 / self.wm**2 + 1j * np.pi / self.lda / self.Rm)
                        * (rr * self.dr) ** 2
                    )
                )

            self.Tz += np.abs(self.E) ** 2 * rr * self.dr

        self.Tznorm = (
            2
            * self.Tz
            / (np.average(self.Tz[0:10]) + np.average(self.Tz[-10:]))
        )
        return self.Tznorm


class Fitting:
    """
    Handles manual and automatic fitting of Z-scan curves.
    No GUI dependencies - all parameters passed explicitly.
    """

    def __init__(
        self,
        integration: Integration,
        amplitude: float,
        beamwaist: float,
        zero_level: float,
        centerpoint: float,
        nop: int,
        y_data: NDArray,
    ):
        self.integration = integration
        self.amplitude = amplitude
        self.beamwaist = beamwaist
        self.zero_level = zero_level
        self.centerpoint = centerpoint
        self.nop = nop
        self.ydata = y_data
        # FIX: Initialize as float instead of None
        self.z_range: float = 0.0

    def manual(
        self,
        zero_level: float,
        centerpoint: float,
        amplitude: float,
        beamwaist: float,
        z_range: float,
        d0: float,
        ra: float,
        stype: str = "CA",
    ) -> NDArray:
        """
        Calculate fitted curve for given parameters (manual fitting).

        Args:
            zero_level: Baseline transmittance level
            centerpoint: Center position of curve [data points]
            amplitude: Phase shift (DPhi0) or transmittance (T)
            beamwaist: Beam waist radius [m]
            z_range: Total z-scan range [m]
            d0: Distance from focus to aperture [m]
            ra: Aperture radius [m]
            stype: "CA" or "OA"

        Returns:
            Fitted transmittance curve
        """
        self.z_range = z_range

        self.integration.z = np.array(
            [
                self.z_range * zz / self.nop - self.z_range / 2 - centerpoint
                for zz in range(self.nop)
            ]
        )

        # Derive with new parameters
        if stype == "CA":
            self.integration.derive(amplitude, beamwaist, d0, ra, stype)

            assert self.integration.closed_sum is not None
            assert self.integration.open_sum is not None

            cas = self.integration.closed_sum
            oas = self.integration.open_sum
            result = (cas / oas) + (zero_level - 1)
        elif stype == "OA":
            self.integration.derive(amplitude, beamwaist, d0, ra, stype)

            assert self.integration.Tznorm is not None

            oas = self.integration.Tznorm
            result = oas + (zero_level - 1)
        else:
            raise ValueError(f"Unknown stype: {stype}")

        return np.asarray(result, dtype=float)

    def fcn2min(
        self,
        params: Parameters,
        d0: float,
        ra: float,
        stype: str,
        weights: Optional[NDArray] = None,
    ) -> NDArray:
        """Function to minimize during automatic fitting."""
        if weights is None:
            weights = np.ones(self.nop)

        vals = params.valuesdict()

        # Extract parameters
        zero_level = vals["Zero"]
        centerpoint = vals["Center"]
        beamwaist = vals["Beamwaist"]
        z_range = vals["Zrange"]

        if stype == "CA":
            amplitude = vals["DPhi0"]
        else:
            amplitude = vals["T"]

        # Calculate fitted curve
        ynew = self.manual(
            zero_level,
            centerpoint,
            amplitude,
            beamwaist,
            z_range,
            d0,
            ra,
            stype,
        )

        # Calculate weighted residuals
        weights = np.asarray(weights, dtype=float).ravel()
        ynew = np.asarray(ynew, dtype=float).ravel()
        ydata = np.asarray(self.ydata, dtype=float).ravel()

        residuals = ynew - ydata
        residuals = np.where(np.isfinite(residuals), residuals, 1e6)

        return np.asarray(weights * residuals**2, dtype=float)

    def automatic(
        self,
        z_range: float,
        d0: float,
        ra: float,
        stype: str,
        vary_beamwaist: bool = True,
        vary_centerpoint: bool = True,
        max_iterations: int = 1000,
        weights: Optional[NDArray] = None,
    ) -> Tuple[object, NDArray]:
        """
        Automatically fit curve using least-squares optimization.

        Args:
            z_range: Total z-scan range [m]
            d0: Distance from focus to aperture [m]
            ra: Aperture radius [m]
            stype: "CA" or "OA"
            vary_beamwaist: Allow beamwaist to vary
            vary_centerpoint: Allow center point to vary
            max_iterations: Maximum optimization iterations
            weights: Optional weight array for data points

        Returns:
            Tuple of (fitted_parameters, result_curve)
        """
        if weights is None:
            weights = np.ones(self.nop)

        self.z_range = z_range

        # Set up fitting parameters
        params = Parameters()
        params.add("Zero", value=self.zero_level, min=0.75, max=1.25)
        params.add(
            "Center",
            value=self.centerpoint,
            min=-50,
            max=50,
            vary=vary_centerpoint,
        )

        if stype == "CA":
            params.add("DPhi0", value=self.amplitude, min=-2, max=2)
            params.add(
                "Beamwaist",
                value=self.beamwaist,
                min=15e-6,
                max=150e-6,
                vary=vary_beamwaist,
            )
        elif stype == "OA":
            params.add("T", value=self.amplitude, min=-2, max=2)
            params.add(
                "Beamwaist",
                value=self.beamwaist,
                min=15e-6,
                max=150e-6,
                vary=False,
            )

        params.add("Zrange", value=z_range, vary=False)

        # Perform fitting
        fitter = Minimizer(
            self.fcn2min,
            params,
            fcn_args=(d0, ra, stype, weights),
        )

        result: MinimizerResult = fitter.minimize(
            method="least_squares", max_nfev=max_iterations
        )

        # Generate final curve with fitted parameters
        params = getattr(result, 'params')
        vals = params.valuesdict()
        zero_level = vals["Zero"]
        centerpoint = vals["Center"]
        beamwaist = vals["Beamwaist"]

        if stype == "CA":
            amplitude = vals["DPhi0"]
        else:
            amplitude = vals["T"]

        result_curve = self.manual(
            zero_level,
            centerpoint,
            amplitude,
            beamwaist,
            z_range,
            d0,
            ra,
            stype,
        )

        return result, result_curve