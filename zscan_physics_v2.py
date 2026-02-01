"""
Physics-based Z-scan fitting engine
Independent of UI - can be tested separately
Provides fit_ca and fit_oa functions
"""

from dataclasses import dataclass
import numpy as np
from typing import Optional

try:
    from lib.integration_headless import Integration, Fitting
    from lib.config import N_COMPONENTS, INTEGRATION_STEPS
except ImportError:
    print("WARNING: lib.integration_headless not found")
    print("Using placeholder values")
    N_COMPONENTS = 50
    INTEGRATION_STEPS = 100


@dataclass
class FittingParams:
    """Container for fitting parameters"""

    amplitude: float  # DPhi0 for CA, T for OA
    beamwaist: float  # meters
    zero_level: float  # 0.75 - 1.25
    centerpoint: float  # data points
    d0: float  # aperture distance (meters)
    ra: float  # aperture radius (meters)

    def get_amplitude_for_physics(self) -> float:
        """Return absolute value for physics calculations"""
        return abs(self.amplitude)


@dataclass
class FittingResult:
    """Result from physics fitting"""

    y_fit: np.ndarray
    params: FittingParams
    r_squared: float
    chi_squared: float
    n2: Optional[float] = None


class ClosedAperturePhysics:
    """Physics for closed aperture (refraction) measurements"""

    @staticmethod
    def estimate_n2(wavelength_nm: float) -> float:
        """Estimate n2 for silica at given wavelength"""
        wavelength = wavelength_nm * 1e-9
        n2 = 2.8203e-20 - 3e-27 / wavelength + 2e-33 / (wavelength**2)
        return n2

    @staticmethod
    def infer_initial_params(
        ca: np.ndarray,
        position_centered_mm: np.ndarray,
        wavelength_nm: float,
    ) -> FittingParams:
        """
        Infer reasonable initial parameters from data geometry
        Uses peak-to-valley heuristics
        """
        pos_m = position_centered_mm * 1e-3  # Convert to meters

        # Find peak and valley
        max_idx = np.argmax(ca)
        min_idx = np.argmin(ca)

        z_pv = abs(pos_m[max_idx] - pos_m[min_idx])
        pv_amp = ca[max_idx] - ca[min_idx]

        # Estimate Rayleigh range and get beam waist
        wavelength = wavelength_nm * 1e-9
        z0 = z_pv / 1.7
        beamwaist = np.sqrt(z0 * wavelength / np.pi)

        # Estimate DPhi0 using proper formula
        # ΔT_p-v ≈ 0.405(1-S)^0.25 |ΔΦ₀| for |ΔΦ₀| ≤ π
        r_a = 0.001  # 1 mm aperture
        r_w = beamwaist  # beam waist radius
        S = 1 - (r_a / r_w) ** 2
        S = np.clip(S, 0, 1)

        coefficient = 0.405 * (1 - S) ** 0.25
        dphi0_magnitude = pv_amp / coefficient

        # ================================================================
        # CRITICAL: Determine sign from peak/valley POSITIONS, not just amplitude
        # ================================================================
        # For positive DPhi0: valley before focus (negative z), peak after (positive z)
        # For negative DPhi0: peak before focus (negative z), valley after (positive z)

        peak_pos = pos_m[max_idx]
        valley_pos = pos_m[min_idx]

        print("\n[infer_initial_params]")
        print(f"  Peak at z={peak_pos*1e3:.1f}mm (index {max_idx})")
        print(f"  Valley at z={valley_pos*1e3:.1f}mm (index {min_idx})")

        if valley_pos < 0 and peak_pos > 0:
            # Valley before focus, peak after → positive DPhi0 ✓
            sign = +1
            print("  → Valley BEFORE focus, peak AFTER → DPhi0 is POSITIVE")
        elif peak_pos < 0 and valley_pos > 0:
            # Peak before focus, valley after → negative DPhi0 ✓
            sign = -1
            print("  → Peak BEFORE focus, valley AFTER → DPhi0 is NEGATIVE")
        elif valley_pos < peak_pos:
            # Valley comes first along position axis
            sign = +1
            print("  → Valley position < peak position → DPhi0 is POSITIVE")
        else:
            # Peak comes first along position axis
            sign = -1
            print("  → Peak position < valley position → DPhi0 is NEGATIVE")

        DPhi0 = sign * dphi0_magnitude
        print(f"  Calculated DPhi0: {DPhi0:+.6f}\n")

        return FittingParams(
            amplitude=DPhi0,
            beamwaist=beamwaist,
            zero_level=1.0,
            centerpoint=0.0,
            d0=0.26,  # 260 mm
            ra=0.001,  # 1 mm aperture
        )

    @staticmethod
    def fit_ca_manual(
        ca: np.ndarray,
        position_centered_mm: np.ndarray,
        wavelength_nm: float,
        params: FittingParams,
    ) -> Optional[FittingResult]:
        """
        Generate CA curve using physics model and provided arguments.

        Args:
            ca: CA data
            position_centered_mm: Position array centered at focal point
            wavelength_nm: Wavelength in nanometers
            params: Fitting parameters

        Returns:
            FittingResult with fit curve and metrics, or None if failed
        """
        try:
            wavelength = wavelength_nm * 1e-9
            z_range = (
                np.max(position_centered_mm) - np.min(position_centered_mm)
            ) * 1e-3
            positions = position_centered_mm * 1e-3  # Convert to meters

            # Calculate n2
            n2 = ClosedAperturePhysics.estimate_n2(wavelength_nm)

            amplitude_mag = (
                params.get_amplitude_for_physics()
            )  # Get positive value

            # Create integration object
            integration = Integration(  # type: ignore
                beta=0,
                n2=n2,
                DPhi0=amplitude_mag,
                positions=positions,
                d0=params.d0,
                aperture_radius=params.ra,
                wavelength=wavelength,
                beamwaist=params.beamwaist,
                n_components=N_COMPONENTS,
                integration_steps=INTEGRATION_STEPS,
                stype="CA",
            )
            
            # Create fitter
            fitter = Fitting(  # type: ignore
                integration=integration,
                amplitude=params.amplitude,
                beamwaist=params.beamwaist,
                zero_level=params.zero_level,
                centerpoint=params.centerpoint,
                nop=len(ca),
                y_data=ca,
            )

            # Generate fit
            y_fit = fitter.manual(
                zero_level=params.zero_level,
                centerpoint=params.centerpoint,
                amplitude=params.amplitude,
                beamwaist=params.beamwaist,
                z_range=z_range,
                d0=params.d0,
                ra=params.ra,
                stype="CA",
            )

            y_fit = np.asarray(y_fit, dtype=float)

            # Calculate metrics
            residuals = y_fit - ca
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((ca - np.mean(ca)) ** 2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
            chi_squared = ss_res / len(ca)

            return FittingResult(
                y_fit=y_fit,
                params=params,
                r_squared=r_squared,
                chi_squared=chi_squared,
                n2=n2,
            )

        except Exception as e:
            print(f"Fitting error: {e}")
            import traceback

            traceback.print_exc()
            return None

    @staticmethod
    def fit_ca_automatic(
        ca: np.ndarray,
        position_centered_mm: np.ndarray,
        wavelength_nm: float,
        params: FittingParams,
    ) -> Optional[FittingResult]:
        """
        Optimize fit to CA data using physics model

        Args:
            ca: CA data
            position_centered_mm: Position array centered at focal point
            wavelength_nm: Wavelength in nanometers
            params: Fitting parameters

        Returns:
            FittingResult with fit curve and metrics, or None if failed
        """
        try:
            wavelength = wavelength_nm * 1e-9
            z_range = (
                np.max(position_centered_mm) - np.min(position_centered_mm)
            ) * 1e-3
            positions = position_centered_mm * 1e-3  # Convert to meters

            # Calculate n2
            n2 = ClosedAperturePhysics.estimate_n2(wavelength_nm)

            amplitude_mag = (
                params.get_amplitude_for_physics()
            )  # Get positive value

            # Create integration object
            integration = Integration(  # type: ignore
                beta=0,
                n2=n2,
                DPhi0=amplitude_mag,
                positions=positions,
                d0=params.d0,
                aperture_radius=params.ra,
                wavelength=wavelength,
                beamwaist=params.beamwaist,
                n_components=N_COMPONENTS,
                integration_steps=INTEGRATION_STEPS,
                stype="CA",
            )

            # Create fitter
            fitter = Fitting(  # type: ignore
                integration=integration,
                amplitude=params.amplitude,
                beamwaist=params.beamwaist,
                zero_level=params.zero_level,
                centerpoint=params.centerpoint,
                nop=len(ca),
                y_data=ca,
            )

            # Generate fit
            fit_params, y_fit = fitter.automatic(
                z_range=z_range,
                d0=params.d0,
                ra=params.ra,
                stype="CA",
                vary_beamwaist=True,
                vary_centerpoint=True
            )

            # Calculate metrics
            residuals = y_fit - ca
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((ca - np.mean(ca)) ** 2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
            chi_squared = ss_res / len(ca)

            return FittingResult(
                y_fit=y_fit,
                params=params,
                r_squared=r_squared,
                chi_squared=chi_squared,
                n2=n2,
            )

        except Exception as e:
            print(f"Fitting error: {e}")
            import traceback

            traceback.print_exc()
            return None


class OpenAperturePhysics:
    """Physics for open aperture (absorption) measurements"""

    @staticmethod
    def fit_oa_manual(
        oa: np.ndarray,
        position_centered_mm: np.ndarray,
        wavelength_nm: float,
        beamwaist_m: float,
        beta: float,
        zero_level: float = 1.0,
        centerpoint: float = 0.0,
        d0_m: float = 0.26,
        ra_m: float = 0.001,
        absorption_model: str = "2PA",
    ) -> Optional[FittingResult]:
        """
        Generate OA curve using physics model and provided arguments.

        Args:
            oa: Open aperture (absorption) data
            position_centered_mm: Position array centered at focal point
            wavelength_nm: Wavelength in nanometers
            beamwaist_m: Beam waist (from silica CA reference)
            beta: Absorption coefficient
            zero_level: Baseline transmittance
            centerpoint: Position offset
            d0_m: Aperture distance (meters)
            ra_m: Aperture radius (meters)
            absorption_model: '2PA', '3PA', '2PA+3PA', etc.

        Returns:
            FittingResult with fit curve and metrics, or None if failed
        """
        try:
            wavelength = wavelength_nm * 1e-9
            z_range = (
                np.max(position_centered_mm) - np.min(position_centered_mm)
            ) * 1e-3
            positions = position_centered_mm * 1e-3  # Convert to meters

            # Create integration object for OA
            integration = Integration(  # type: ignore
                beta=beta,
                n2=0,  # No refraction for OA
                DPhi0=0,
                positions=positions,
                d0=d0_m,
                aperture_radius=ra_m,
                wavelength=wavelength,
                beamwaist=beamwaist_m,
                n_components=N_COMPONENTS,
                integration_steps=INTEGRATION_STEPS,
                stype="OA",
            )

            # Create fitter
            fitter = Fitting(  # type: ignore
                integration=integration,
                amplitude=beta,
                beamwaist=beamwaist_m,
                zero_level=zero_level,
                centerpoint=centerpoint,
                nop=len(oa),
                y_data=oa,
            )

            # Generate fit
            y_fit = fitter.manual(
                zero_level=zero_level,
                centerpoint=centerpoint,
                amplitude=beta,
                beamwaist=beamwaist_m,
                z_range=z_range,
                d0=d0_m,
                ra=ra_m,
                stype="OA",
            )

            y_fit = np.asarray(y_fit, dtype=float)

            # Calculate metrics
            residuals = y_fit - oa
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((oa - np.mean(oa)) ** 2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
            chi_squared = ss_res / len(oa)

            return FittingResult(
                y_fit=y_fit,
                params=FittingParams(
                    amplitude=beta,
                    beamwaist=beamwaist_m,
                    zero_level=zero_level,
                    centerpoint=centerpoint,
                    d0=d0_m,
                    ra=ra_m,
                ),
                r_squared=r_squared,
                chi_squared=chi_squared,
                n2=0,  # Not applicable for OA
            )

        except Exception as e:
            print(f"OA Fitting error: {e}")
            import traceback

            traceback.print_exc()
            return None

    @staticmethod
    def fit_oa_automatic(
        oa: np.ndarray,
        position_centered_mm: np.ndarray,
        wavelength_nm: float,
        beamwaist_m: float,
        beta: float,
        zero_level: float = 1.0,
        centerpoint: float = 0.0,
        d0_m: float = 0.26,
        ra_m: float = 0.001,
        absorption_model: str = "2PA",
    ) -> Optional[FittingResult]:
        """
        Fit OA data using absorption model

        Args:
            oa: Open aperture (absorption) data
            position_centered_mm: Position array centered at focal point
            wavelength_nm: Wavelength in nanometers
            beamwaist_m: Beam waist (from silica CA reference)
            beta: Absorption coefficient
            zero_level: Baseline transmittance
            centerpoint: Position offset
            d0_m: Aperture distance (meters)
            ra_m: Aperture radius (meters)
            absorption_model: '2PA', '3PA', '2PA+3PA', etc.

        Returns:
            FittingResult with fit curve and metrics, or None if failed
        """
        try:
            wavelength = wavelength_nm * 1e-9
            z_range = (
                np.max(position_centered_mm) - np.min(position_centered_mm)
            ) * 1e-3
            positions = position_centered_mm * 1e-3  # Convert to meters

            # Create integration object for OA
            integration = Integration(  # type: ignore
                beta=beta,
                n2=0,  # No refraction for OA
                DPhi0=0,
                positions=positions,
                d0=d0_m,
                aperture_radius=ra_m,
                wavelength=wavelength,
                beamwaist=beamwaist_m,
                n_components=N_COMPONENTS,
                integration_steps=INTEGRATION_STEPS,
                stype="OA",
            )

            # Create fitter
            fitter = Fitting(  # type: ignore
                integration=integration,
                amplitude=beta,
                beamwaist=beamwaist_m,
                zero_level=zero_level,
                centerpoint=centerpoint,
                nop=len(oa),
                y_data=oa,
            )

            # Generate fit
            y_fit = fitter.manual(
                zero_level=zero_level,
                centerpoint=centerpoint,
                amplitude=beta,
                beamwaist=beamwaist_m,
                z_range=z_range,
                d0=d0_m,
                ra=ra_m,
                stype="OA",
            )

            y_fit = np.asarray(y_fit, dtype=float)

            # Calculate metrics
            residuals = y_fit - oa
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((oa - np.mean(oa)) ** 2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
            chi_squared = ss_res / len(oa)

            return FittingResult(
                y_fit=y_fit,
                params=FittingParams(
                    amplitude=beta,
                    beamwaist=beamwaist_m,
                    zero_level=zero_level,
                    centerpoint=centerpoint,
                    d0=d0_m,
                    ra=ra_m,
                ),
                r_squared=r_squared,
                chi_squared=chi_squared,
                n2=0,  # Not applicable for OA
            )

        except Exception as e:
            print(f"OA Fitting error: {e}")
            import traceback

            traceback.print_exc()
            return None


if __name__ == "__main__":
    # Test with your silica data
    from zscan_data_parser_v2 import ZScanFileParser, ZScanProcessor

    print("=" * 70)
    print("PHYSICS ENGINE TEST")
    print("=" * 70)

    filepath = "data/2021_07_22__10_10__silica_0-0_1600-0_2.txt"

    raw = ZScanFileParser.parse(filepath)
    if not raw:
        print(f"✗ Could not parse {filepath}")
        exit(1)

    processed = ZScanProcessor.normalize(raw)
    print(f"✓ Loaded: {processed.sample_code}")
    print(f"  Points: {len(processed.ca)}")
    print(f"  λ: {processed.wavelength_nm} nm")
    print()

    # Infer initial parameters
    initial_params = ClosedAperturePhysics.infer_initial_params(
        processed.ca,
        processed.position_centered,
        processed.wavelength_nm,
    )

    print("Initial parameters (inferred from data):")
    print(f"  DPhi0: {initial_params.amplitude:.4f} rad")
    print(f"  w0: {initial_params.beamwaist * 1e6:.2f} µm")
    print(f"  zero_level: {initial_params.zero_level:.4f}")
    print(f"  d0: {initial_params.d0 * 1000:.1f} mm")
    print(f"  ra: {initial_params.ra * 1000:.2f} mm")
    print()

    # Try fitting
    print("Running physics fit...")
    result = ClosedAperturePhysics.fit_ca_manual(
        processed.ca,
        processed.position_centered,
        processed.wavelength_nm,
        initial_params,
    )

    if result:
        print("✓ Fit successful!")
        print(f"  R²: {result.r_squared:.6f}")
        print(f"  χ²: {result.chi_squared:.6f}")
        print(f"  n2: {result.n2:.3e} m²/W")
        print()
        print("Fitted parameters:")
        print(f"  DPhi0: {result.params.amplitude:.4f} rad")
        print(f"  w0: {result.params.beamwaist * 1e6:.2f} µm")
        print(f"  zero_level: {result.params.zero_level:.4f}")
    else:
        print("✗ Fit failed")
