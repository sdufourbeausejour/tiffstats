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
# - find_matching_file_list : returns a list of matching files (list)

from __future__ import print_function
import fileinput
import os
import platform
import pprint as pp

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


def find_matching_file_list(directory, exp_file_name, print_list=0):
    """Get list of files in directory matching exp_file_name regex.
    Optionaly print the list"""
    system_name = platform.system()
    file_paths = os.listdir(directory)
    matching_file_paths = []
    for index, file_path in enumerate(file_paths):
        # Try to match file name to regex
        # If matching fails, continue to the top of the loop
        try:
            matched_pattern = exp_file_name.search(file_path).group(0)
            matching_file_paths.append(file_path)
        except:
            continue

        try:
            if print_list:
                print("Found matching files : ")
                pp.pprint(matching_file_paths)
        except:
            print("Found matching files : ")
            pp.pprint(matching_file_paths)

    return matching_file_paths
