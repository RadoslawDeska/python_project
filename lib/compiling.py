import logging
import os
import traceback
from pathlib import Path  # noqa: F401

from PyQt5 import uic


def compile_GUI(uif):
    """uif is the path to the file with the file name and extension *.ui"""
    try:
        if Path(uif).is_file():
            with open(os.path.join("./lib/", Path(uif).stem + "_gui.py"), "w", encoding="utf-8") as pyf:
                uic.compileUi(uif, pyf)
            return True  # Compilation succeeded
        else:
            raise FileNotFoundError
    except FileNotFoundError:
        print(f"Couldn't compile the UI file {uif}. File not found.")
        logging.error(traceback.format_exc())
        return False  # Compilation failed
