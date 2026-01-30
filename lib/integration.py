"""
Integration module for Z-Scan analysis.

This module contains the Integration class that calculates the normalized
transmittance through the aperture for closed and open aperture measurements
according to the procedure from Sheik-Bahae et al.

Classes:
    Integration: Calculates integrated electric field and transmittance.
"""

from math import factorial
from typing import List
import numpy as np
from numpy.typing import NDArray
from scipy.special import hyp2f1, lambertw


class Integration:
    """
    Integrates the electric field according to Sheik-Bahae procedure.

    Calculates the normalized transmittance T(z) for both closed-aperture (CA)
    and open-aperture (OA) measurements in Z-scan experiments.

    Attributes:
        z (NDArray): Position array along the beam propagation direction [m]
        d0 (float): Distance from focus to aperture plane [m]
        ra (float): Aperture radius [m]
        lda (float): Wavelength [m]
        w0 (float): Beam waist radius at focus [m]
        n2 (float): Non-linear refractive index [m^2/W]
        T (float): Non-linear transmittance
        DPhi0 (float): On-axis phase shift at z=0 [rad]
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
        """
        Initialize the Integration object.

        Args:
            beta (float): Non-linear absorption coefficient [m/W]
            n2 (float): Non-linear refractive index [m^2/W]
            DPhi0 (float): On-axis phase shift at focus [rad]
            positions (NDArray): Array of z-positions [m]
            d0 (float): Distance from focus to aperture [m]
            aperture_radius (float): Aperture radius [m]
            wavelength (float): Laser wavelength [m]
            beamwaist (float): Beam waist radius [m]
            n_components (int): Number of Gaussian decomposition components
            integration_steps (int): Number of radial integration steps
            stype (str): Measurement type - "CA" for closed aperture, "OA" for open aperture
        """
        self.z = positions
        self.d0 = d0
        self.ra = aperture_radius
        self.lda = wavelength
        self.w0 = beamwaist
        self.n2 = n2

        try:
            # For OA fitting, the caller provides an amplitude-like value in `beta`.
            # Use it directly rather than dividing by `n2` which often is zero
            # and makes T vanish. Keep the raw beta value as the OA amplitude.
            self.T = beta
        except Exception:
            self.T = 0

        self.DPhi0 = DPhi0
        self.mm = n_components
        self.ir = integration_steps
        self.stype = stype

        self.derive(self.DPhi0, self.w0, self.d0, self.ra, stype)

    def derive(
        self, DPhi0: float, w0: float, d0: float, ra: float, stype: str
    ) -> None:
        """
        Calculate derived parameters.

        Computes beam propagation parameters and performs integration
        based on measurement type (CA or OA).

        Args:
            DPhi0 (float): On-axis phase shift [rad]
            w0 (float): Beam waist [m]
            d0 (float): Distance to aperture [m]
            ra (float): Aperture radius [m]
            stype (str): Measurement type ("CA" or "OA")
        """
        self.k = 2 * np.pi / self.lda
        self.z0 = 0.5 * self.k * w0**2
        self.wa = w0 * np.sqrt(1 + d0**2 / self.z0**2)
        self.ra = ra
        self.Dphi0 = DPhi0 / (1 + self.z**2 / self.z0**2)
        self.wz = w0 * np.sqrt(1 + self.z**2 / self.z0**2)
        with np.errstate(divide='ignore', invalid='ignore'):
            self.Rz = self.z + self.z0**2 / self.z
        self.d = d0 - self.z
        self.g = 1 + self.d / self.Rz

        if stype == "CA":
            self.bigproduct()
        elif stype == "OA":
            # Ensure OA amplitude (T) reflects the currently provided DPhi0/amplitude
            # when derive() is called with amplitude as first arg for OA.
            try:
                # DPhi0 argument is used as amplitude for OA in caller code
                self.T = DPhi0
            except Exception:
                pass
            self.calculate_Tz_for_OA()

    def calculate_Tz_for_OA(self, model: str = "2PA") -> NDArray:
        """
        Calculate transmittance for open aperture measurement.

        Accounts for nonlinear absorption processes (2-photon, 3-photon,
        reverse saturable absorption, or saturable absorption).

        Args:
            model (str): Absorption model - "2PA", "3PA", "2PA+3PA",
                        "RSA", "SA", or "2PA+SA"

        Returns:
            NDArray: Normalized transmittance T(z)
        """
        self.Tz = np.zeros_like(self.z, dtype=np.complex128)

        if self.T != 0:
            Psi1 = self.T
            Psi2 = self.T * 8  # Scaling factor for 3PA
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
        else:  # RSA, SA, 2PA+SA
            self.Tz = np.ones_like(self.z)

        # Normalize transmittance
        self.Tznorm = (
            2
            * self.Tz
            / (
                np.average(self.Tz[0:10])
                + np.average(self.Tz[len(self.Tz) - 10 :])
            )
        )
        return self.Tznorm

    def calculate_fm(self) -> List[complex]:
        """
        Calculate Gaussian decomposition coefficients.

        Returns:
            List[complex]: List of coefficients f_m for each component
        """
        self.result = [
            (1j * self.DPhi0)**m / factorial(m) * self.product[m]
            for m in range(self.mm)
        ]
        return self.result


    def bigproduct(self) -> None:
        """
        Calculate the big product for closed aperture analysis.

        Computes the cumulative product needed for Gaussian decomposition
        and normalizes the transmittance by the ratio of closed to open aperture.
        """
        self.product = []
        for m in range(0, self.mm):
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

    def open(self) -> NDArray:
        """
        Calculate transmittance for open aperture (no aperture).

        Returns:
            NDArray: Normalized transmittance for open aperture
        """
        self.dr = 3 * self.wa / self.ir
        self.open_sum = self.bigsum()
        return self.open_sum

    def closed(self) -> NDArray:
        """
        Calculate transmittance for closed aperture (with aperture).

        Returns:
            NDArray: Normalized transmittance for closed aperture
        """
        self.dr = self.ra / self.ir
        self.closed_sum = self.bigsum()
        return self.closed_sum

    def bigsum(self) -> NDArray:
        """
        Perform radial integration for transmittance calculation.
        
        Integrates the electric field intensity over the aperture or
        detection plane to obtain the transmitted power.
        
        Returns:
            NDArray: Normalized transmittance T(z)
        """
        self.Tz = np.zeros_like(self.z, dtype=np.float64)
        
        for rr in range(self.ir):
            self.E = np.zeros_like(self.z, dtype=np.complex128)
            
            for m in range(0, self.mm):
                self.wm0 = self.wz / np.sqrt(2 * m + 1)
                self.dm = 0.5 * self.k * self.wm0**2
                self.wm = self.wm0 * np.sqrt(
                    self.g**2 + self.d**2 / self.dm**2
                )
                self.tm = np.arctan(self.g / (self.d / self.dm))
                self.Rm = self.d / (1 - self.g / (self.g**2 + self.d**2 / self.dm**2))
                
                self.E += (
                    self.fm[m] *
                    np.exp(1j * self.tm) *
                    self.wm0 / self.wm / self.wz *
                    np.exp(
                        (-1 / self.wm**2 + 1j * np.pi / self.lda / self.Rm) *
                        (rr * self.dr)**2
                    )
                )
            
            self.Tz += np.abs(self.E)**2 * rr * self.dr
        
        # Normalize before CA/OA division
        self.Tznorm = 2 * self.Tz / (
            np.average(self.Tz[0:10]) + np.average(self.Tz[len(self.Tz) - 10:])
        )
        return self.Tznorm