"""
Scientific rounding module for error propagation.

This module provides functions for rounding numerical values and their
uncertainties according to standard scientific notation practices.
"""

from decimal import Decimal, ROUND_UP, ROUND_HALF_UP
from typing import Any, Optional, Tuple
import math

def error_rounding(
    value: Any,
    error: Any
) -> Tuple[float, Optional[float], Optional[int]]:
    """
    Scientific rounding with Decimal. Always returns a float for the value (0.0 on error for `value`)
    """

    # --- Normalize value to Decimal ---
    try:
        value_dec = Decimal(str(value))
    except Exception:
        # fallback: try float
        try:
            value_dec = Decimal(str(float(value)))
        except Exception:
            value_dec = Decimal("0")

    # --- Normalize error to Decimal ---
    try:
        # Handle invalid error early
        if error is None:
            raise Exception
        err_float = float(error)
        if math.isnan(err_float) or err_float <= 0:
            raise Exception
        error_dec = Decimal(str(err_float))
    except Exception:
        # Invalid error → return value unchanged, no error, no precision
        return float(value_dec), None, None

    # --- Determine order of magnitude ---
    error_mag = error_dec.adjusted()

    # Scale error to [1, 10)
    scaled = error_dec.scaleb(-error_mag)

    # --- Round error to 1–2 significant digits ---
    if scaled < Decimal("1.5"):
        rounded_scaled = scaled.quantize(Decimal("0.01"), rounding=ROUND_UP)
        sig_digits = 2
    else:
        rounded_scaled = scaled.quantize(Decimal("0.1"), rounding=ROUND_UP)
        sig_digits = 1

    rounded_error_dec = rounded_scaled.scaleb(error_mag)

    # --- Determine decimal places ---
    decimal_places = max(0, -(error_mag) + (sig_digits - 1))

    # --- Round value to match error precision ---
    quant = Decimal(1).scaleb(-decimal_places)
    rounded_value_dec = value_dec.quantize(quant, rounding=ROUND_HALF_UP)

    return float(rounded_value_dec), float(rounded_error_dec), decimal_places
