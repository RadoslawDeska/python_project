from decimal import *
from unicodedata import decimal
import numpy as np
import traceback
import logging

# Round value according to uncertainty
def error_rounding(value, error):
    decimal_places = 15 # This precision is needed for the values in zscan1.py program
    precision = Decimal(10) ** -decimal_places

    try:
        val_sign, val_digits, val_exp = Decimal(value).quantize(precision).as_tuple()
        err_sign, err_digits, err_exp = Decimal(error).quantize(precision).as_tuple()
    except Exception as e:
        logging.error(traceback.format_exc())
        print(f"{value=}\n{error=}\n{precision=}")
        return float(value), float(error), decimal_places

    value = Decimal((val_sign, val_digits, val_exp))
    error = Decimal((err_sign, err_digits, err_exp))

    #print(f'Taking {value} and {error}.')

    try:
        error_mag = np.floor(error.log10())
    except Exception as e:
        logging.error(traceback.format_exc())
        print(f"{value=}\n{error=}")
        return float(value), float(error), decimal_places
    
    #print(f'Magnitude of the error is {error_mag}.')
    
    try:
        error_amp = Decimal(float(error)*10**-(error_mag)).quantize(precision)
    except Exception as e:
        logging.error(traceback.format_exc())
        print(f"{value=}\n{error=}\n{error_mag=}")
    #print(f'Error amplitude is {error_amp}.')

    try:
        rounded_error_amp = error_amp.quantize(0, rounding=ROUND_UP)
    except Exception as e:
        logging.error(traceback.format_exc())
        print(f"{value=}\n{error=}\n{error_mag=}\n{error_amp=}")
    #print(f'Rounded error amplitude is equal to {rounded_error_amp}.')

    try:
        if Decimal((rounded_error_amp-error_amp)/error_amp) > 0.1:
            rounded_error_amp = error_amp.quantize(Decimal('.1'), rounding=ROUND_UP)
            #print(f'It is more than 10%. New rounded error amplitude is equal to {rounded_error_amp}.')
            
        else:
            pass
            #print(f'It is less than or equal to 10% of the original value.')
    
    except Exception as e:
        logging.error(traceback.format_exc())
        print(f"{value=}\n{error=}\n{error_mag=}\n{error_amp=}\n{rounded_error_amp=}")

    _, round_err_dig, round_err_exp = rounded_error_amp.as_tuple()
    for _ in range(len(err_digits)-len(round_err_dig)):
            round_err_dig = round_err_dig + (0,)
    
    # Check if rounding didn't cause change in order of magnitude
    error_amp_mag = np.floor(error_amp.log10())
    rounded_error_amp_mag = np.floor(rounded_error_amp.log10())

    if rounded_error_amp_mag != error_amp_mag:
        #print(f"{rounded_error_amp_mag = }")
        #print(f"{error_amp_mag = }")
        
        round_err_dig = round_err_dig + (0,)
        round_err_exp = round_err_exp + 1
        
    rounded_error = Decimal((err_sign, round_err_dig, err_exp)).quantize(Decimal(10) ** Decimal(error_mag+round_err_exp))
    
    #print(val_digits, val_exp)
    #print(err_digits, err_exp)
    #print(round_err_dig, round_err_exp)
    #print(f'Rounded error is equal to {rounded_error}.')
    
    #rounded_error_mag = int(np.floor(np.log10(float(rounded_error))))
    #print(f'Magnitude of the rounded error is {rounded_error_mag}.')

    rounded_value = Decimal(value).quantize(rounded_error, rounding=ROUND_HALF_UP)

    #print(f'Rounded value is equal to {rounded_value}.')

    _,_, exp = rounded_value.as_tuple()

    return float(rounded_value), float(rounded_error), np.abs(exp)