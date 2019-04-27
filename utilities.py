#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tools for
"""

# -------- Import libraries -------- #
import pexpect
from itertools import islice

def decimal_to_rgb_gray(decimal, minimum=0, maximum=255):
    """Converts a decimal value (0-1) to a hex RGB grayscale value between a max and min (default is 0-255)"""

    # Calculate the hex value of the decimal
    scale = maximum - minimum
    temp_string = hex(int(decimal*scale + minimum))

    # Remove the 0x from the start of the string
    hex_scale = temp_string.split('x')[1]

    # Test if hex_scale is a single digit
    if len(hex_scale) == 1:
        hex_scale = '0' + hex_scale

    # Repeat the string three times and add a # to the start
    hex_string = '#' + hex_scale*3
    return(hex_string)


def fasta(dna_sequence_file, protein_input, output_file):
    """
    Runs FASTA given a DNA sequence FASTA file, a protein FASTA file, and the name of the file to output to.

    The first should be a DNA sequence file, and the second should be
    a protein sequence (amino acid) file.
    """

    # Download FASTA
    command = 'bash -c "if ! [ -e ./fasta-36.3.8g/bin ]; then wget -q https://github.com/wrpearson/fasta36/releases/download/fasta-v36.3.8g/fasta-36.3.8g-linux64.tar.gz && tar -xzf fasta-36.3.8g-linux64.tar.gz; fi"'
    pexpect.run(command)

    # Run FASTX using variables
    command = 'bash -c ' + '"if [ \':$PATH:\' != *\':./fasta-36.3.8g/bin\'* ]; then PATH=\'./fasta-36.3.8g/bin\'; fi && fastx36 -m 8 -E 0.05 ' + dna_sequence_file + ' ' + protein_input + ' > ' + output_file + '"'
    pexpect.run(command)


def sequence_finder(dna_file, plasmid):
    """Finds corresponding DNA sequence in file given plasmid name"""

    # Variable setup
    line_number = 0
    sequence = ''
    recording = False

    with open(dna_file) as dna_input:
        # Move to line identifying the plasmid being searched for
        for line in islice(dna_input, None):
            line_number += 1
            if recording:
                if '>' not in line:
                    sequence += line.strip()
                else:
                    break

            if plasmid in line:
                recording = True

    sequence = sequence.strip().replace(" ", "")

    return sequence