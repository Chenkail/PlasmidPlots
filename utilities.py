"""Miscellaneous utilities"""

# -------- Import libraries -------- #
import pexpect
from itertools import islice

def decimal_to_rgb_gray(decimal, darkest=0, lightest=255, reverse=False, lower=0, upper=1):
    """
    Converts a decimal value to a hex RGB grayscale value
    between a max and min shade (default is 0-255).

    By default, lower values mean lighter shades of gray.
    The function will return a value based on the location of the input
    between the given lower and upper values for input
    (by default 0 and 1, respectively)
    """

    # Test for input outside of min/max bounds
    if decimal < lower or decimal > upper:
        return "#FFFFFF"

    # Calculate the hex value of the decimal
    # Calculate percentagewise where the input is between the min and max
    scale = upper - lower
    percent = (decimal - lower)/scale

    # If scale is not reversed, flip so that the highest decimal
    # value corresponds with the darkest shade
    if reverse == False:
        percent = 1 - percent

    # Use percentage and apply to location between darkest and lightest shades of gray
    shade_scale = lightest - darkest
    value = percent*shade_scale + darkest

    temp_string = hex(int(value))

    # Remove "0x" from the start of the string
    hex_scale = temp_string.split('x')[1]

    # Pad hex_scale to length 2
    while len(hex_scale) < 2:
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