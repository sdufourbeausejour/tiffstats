# -*- coding: utf-8 -*-

# Sophie Dufour-Beauséjour
# s.dufour.beausejour@gmail.com
# Doctorante en sciences de l"eau
# INRS-ETE
# Québec

# Text file writing utils
# This file is imported in main.py
# Contents:
# - line_pre_adder: adds line at the beginning of a file ()

from __future__ import print_function
import fileinput

def line_pre_adder(filename, line_to_prepend):
    """Add a line at the beginning of a text file
    from http://stackoverflow.com/questions/5914627/prepend-line-to-beginning-of-a-file
    """
    f = fileinput.input(filename, inplace=1)
    for xline in f:
        if f.isfirstline():
            print(line_to_prepend.rstrip("\r\n") + "\n" + xline.rstrip("\n"))
        else:
            print(xline.rstrip("\n"))
