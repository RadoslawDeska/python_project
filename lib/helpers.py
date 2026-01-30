from typing import Any

def as_float(x: Any, default: float = 0.0) -> float:
    try:
        if x is None:
            return default
        return float(x)
    except Exception:
        return default


def ensure_float(x: float | None, default: float = 0.0) -> float:
    return default if x is None else x
