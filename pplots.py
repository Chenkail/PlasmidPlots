#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plots sequences corresponding to protein families on plasmids"""

# -------- Import libraries -------- #
# Built-in
import argparse
import math
import os
from itertools import islice
from timeit import default_timer as timer

# 3rd party
from Bio import SeqIO
from Bio.SeqUtils import GC, GC_skew
import imagemergetools as imt
from PIL import Image
import matplotlib.pyplot as plt
import numpy as np
import ncbi_tools as ncbi
from plasmidplots import utilities as pputil


# -------- Authorship -------- #
__author__ = "Chenkai Luo and Chris Lausted"
__license__ = "GPLv3"
__version__ = "1.6.0"

# -------- Begin functions -------- #

def read_colors(color_file, subgroup_list_file=None):
    """
    Reads colors from a text file.

    File should include each color on a separate line, with each line
    consisting of the name of the corresponding protein family
    (as listed in the protein family text file) followed by a colon and the
    hex code of the corresponding color. The base color should be prefixed by
    'Base Color'. If no base color is specified, the program will use #C0C0C0
    as the default.

    Example file:
    Base Color:#C0C0C0
    32:#FF0000
    49:#00FF00
    50:#FFFF00
    """

    color_dict = {}
    subgroup_file_dict = {}
    if subgroup_list_file != None:
        with open(subgroup_list_file) as file:
            for line in islice(color_file, None):
                subgroup, subgroup_file = line.strip().split(":")
                subgroup_file_dict[subgroup] = subgroup_file
    else:
        subgroup_file_dict = {'32':'pf32.txt'}

    with open(color_file) as file:
        for line in islice(file, None):
            family, color = line.strip().split(":")

            # Add color to dictionary
            color_dict[family] = color

    if 'Base Color' not in color_dict.keys():
        color_dict['Base Color'] = '#C0C0C0'

    return color_dict, subgroup_file_dict


def gc_content_dict(dna_file, plasmid, window=100):
    """Returns list of GC content percentages for each window"""

    # Pull DNA sequence into a string
    sequence = pputil.sequence_finder(dna_file, plasmid)

    length = len(sequence)
    chunks = math.ceil(length/window)

    gc_data_dict = {}
    for i in range(chunks):
        start = i*window

        if i == chunks - 1:
            end = length

        else:
            end = (i+1) * window

        sequence_chunk = sequence[start:end]

        gc_content = GC(sequence_chunk)/100

        key = str(start + 1) + "-" + str(end)
        gc_data_dict[key] = gc_content

    return gc_data_dict


def gc_skew_dict(dna_file, plasmid, window=100):
    """
    Returns dictionary of locations and GC skew

    Takes dna file, sequence name, and (optionally) window size as arguments
    """

    # Pull DNA sequence into a string
    sequence = pputil.sequence_finder(dna_file, plasmid)

    length = len(sequence)

    # Get list of GC skew values
    skew_list = GC_skew(sequence, window)

    chunks = math.ceil(length/window)

    gc_skew_dict = {}
    for i in range(chunks):
        start = i*window

        if i == chunks - 1:
            end = length

        else:
            end = (i+1) * window

        key = str(start + 1) + "-" + str(end)
        skew = skew_list[i]
        gc_skew_dict[key] = skew

    return gc_skew_dict


def linear_plot(plasmid, data, sequence_color_dict,
                baseline_custom_colors=None,
                output_dir="./plasmidplots_temp/"):
    """Plots a linear plasmid given the plasmid name and data list"""

    # Set constants
    HORIZONTAL_SCALE_CONSTANT = 3/4000
    # Offsets to avoid labels intersecting plot
    LABEL_Y_ADJUST = 0.1
    LABEL_LEFT_X_ADJUST = -1000
    LABEL_RIGHT_X_ADJUST = 500
    LABEL_UP_ADJUST = 0.8
    LABEL_DOWN_ADJUST = -1

    # Set variables
    sequence_name = plasmid
    plasmid_length = data[0]
    data_count = len(data)

    # Scaling
    horizontal_scale = plasmid_length * HORIZONTAL_SCALE_CONSTANT
    plt.figure()
    plt.rcParams['figure.figsize'] = (horizontal_scale, 0.5)
    plt.close('all')

    linear = plt.subplot(111)
    linear.set_facecolor('white')
    background = linear.barh(0, plasmid_length, height=20, color='white')

    if baseline_custom_colors != None:
        # Use dictionary to plot baseline with custom colors
        for location, color in baseline_custom_colors.items():
            first_base = int(location.split('-')[0])
            last_base = int(location.split('-')[1])
            length = last_base - first_base
            baseline = linear.barh(0, length, left=first_base,
                                   height=1, color=color)

    else:
        baseline = linear.barh(0, plasmid_length, height=1,
                               color=sequence_color_dict['Base Color'])

    for index in range(1, data_count):
        sequence = data[index]
        # Get sequence data
        sequence_start = min(sequence[1], sequence[2])
        sequence_end = max(sequence[1], sequence[2])
        sequence_length = sequence_end - sequence_start
        sequence_family = str(sequence[0])

        # Plot gene onto plasmid
        current_color = sequence_color_dict[sequence_family]
        gene_plot = linear.barh(0, sequence_length, left=sequence_start,
                                height=1, color=current_color)

        # Add subgroup labels
        subgroup = sequence[3]
        if subgroup != '':
            label_x = (sequence_start + sequence_end)/2
            label_y = LABEL_DOWN_ADJUST

            # Add line connecting label to plasmid
            linear.annotate(subgroup,
                            xy=(label_x, 0),
                            xytext=(label_x, label_y),
                            xycoords='data',
                            ha='center',
                            va='center',
                            fontsize=9,
                            rotation=90,
                            arrowprops=dict(facecolor='black',
                                            arrowstyle='-',))

    # Add plot labels
    label_y = LABEL_Y_ADJUST
    # Left side label - name of plasmid (second column)
    label_x = LABEL_LEFT_X_ADJUST
    plt.text(label_x, label_y, "0")
    # Right side label - length of plasmid
    label_x = plasmid_length + LABEL_RIGHT_X_ADJUST
    plt.text(label_x, label_y, str(plasmid_length))
    # Top label - plasmid name
    label_x = 0
    label_y = LABEL_UP_ADJUST
    label_text = str(sequence_name) + " - " + str(plasmid_length) + " base pairs"
    plt.text(label_x, label_y, label_text, fontsize=14)

    # Fix axes limits
    xmin, xmax, ymin, ymax = plt.axis('tight')
    plt.ylim(0, 1)

    # Hide background and save plot to image
    plt.axis('off')
    filename = output_dir + "temp.png"
    plt.savefig(filename, bbox_inches='tight')
    plt.close('all')


def circular_plot(plasmid, data, sequence_color_dict,
                  baseline_custom_colors=None,
                  output_dir="./plasmidplots_temp/"):
    """Plots a circular plasmid given the name and the data list"""

    # Set constants
    CHART_BOTTOM = 4
    CHART_THICKNESS = 0.7
    CIRCULAR_SCALE_CONSTANT = 4

    # Set variables
    sequence_name = plasmid
    plasmid_length = data[0]
    data_count = len(data)

    # Sizing
    graph_scale = np.sqrt(plasmid_length/9000) * CIRCULAR_SCALE_CONSTANT
    circle_width = CHART_THICKNESS/graph_scale * CIRCULAR_SCALE_CONSTANT
    circle_bottom = CHART_BOTTOM
    plt.figure()
    plt.rcParams['figure.figsize'] = (graph_scale, graph_scale)
    plt.close('all')

    circle = plt.subplot(111, polar=True)

    # Use data from first protein family sequence
    # for initial minimum/maximum values
    first_data_point = data[1]
    minimum = min(first_data_point[1], first_data_point[2])
    maximum = max(first_data_point[1], first_data_point[2])

    # Center the plot
    for index in range(1, data_count):
        sequence = data[index]
        sequence_start = min(sequence[1], sequence[2])
        sequence_end = max(sequence[1], sequence[2])
        if sequence_start < minimum:
            minimum = sequence_start
        if sequence_end > maximum:
            maximum = sequence_end

    # Center the plot
    center = int(minimum + maximum)/2
    offset_angle = -(center/plasmid_length) * np.radians(360) + np.radians(90)
    circle.set_theta_offset(offset_angle)

    # Plot baseline
    if baseline_custom_colors != None:
        # Use dictionary to plot baseline with custom colors
        for location, color in baseline_custom_colors.items():
            first_base = int(location.split('-')[0])
            last_base = int(location.split('-')[1])
            length = last_base - first_base
            gene_plot_start = ((first_base) / plasmid_length) * np.radians(360)
            gene_plot_width = (length / plasmid_length) * np.radians(360)
            baseline = circle.bar(gene_plot_start,
                          circle_width,
                          width=gene_plot_width,
                          bottom=circle_bottom,
                          align='edge',
                          color=color,)
    else:
        baseline = circle.bar(0,
                          circle_width,
                          width=-np.radians(360),
                          bottom=circle_bottom,
                          align='edge',
                          color=sequence_color_dict['Base Color'],)

    # Plot data
    for index in range(1, data_count):
        sequence = data[index]
        # Get sequence data
        sequence_start = min(sequence[1], sequence[2])
        sequence_end = max(sequence[1], sequence[2])
        sequence_length = sequence_end - sequence_start
        sequence_family = str(sequence[0])

        # Plot gene onto plasmid
        gene_plot_start = ((sequence_start) / plasmid_length) * np.radians(360)
        gene_plot_width = (sequence_length / plasmid_length) * np.radians(360)
        current_color = sequence_color_dict[sequence_family]
        gene_plot = circle.bar(gene_plot_start,
                               circle_width,
                               width=gene_plot_width,
                               bottom=circle_bottom,
                               color=current_color,)

        # Label subgroup
        subgroup = sequence[3]
        if subgroup != '':
            # Get angles, centered on protein segment
            line_angle = gene_plot_start

            # Make sure angle is between 0 and 2π
            line_angle %= (2*np.pi)

            # Rotate label text to point towards center of circle
            label_angle = np.degrees(line_angle + offset_angle)

            # Flip text if on left side of plot so it's not upside down
            if label_angle > 90 and label_angle < 270:
                label_angle += 180

            # Add annotation to plot
            circle.annotate(subgroup,
                            xy=(line_angle, 4.2),
                            xytext=(line_angle, 5.4),
                            xycoords='data',
                            ha='center',
                            va='center',
                            fontsize=12,
                            rotation=label_angle,
                            arrowprops=dict(facecolor='black',
                                            arrowstyle='-',),)

    # Label - plasmid name and length
    label_x = 0.5
    if plasmid_length >= 25000:
        label_y = 0.5
    else:
        label_y = -0.25

    label_text = (str(sequence_name) + "\n"
                    + str(plasmid_length) + " base pairs")

    plt.text(label_x, label_y, label_text, fontsize=24, ha='center',
             transform=circle.transAxes)

    # Hide background and save plot to image
    plt.axis('off')
    filename = output_dir + "temp.png"
    plt.savefig(filename, bbox_inches='tight')
    plt.close('all')


def subgroup_search(sequence, family, subgroup_file_dict):
    """
    Find best subgroup match using FASTA

    Takes a sequence string and the protein family to search
    """

    # Setup
    subgroup_file = subgroup_file_dict[family]
    short_sequence_file = "plasmidplots_temp/short_sequence.txt"

    # Create blank file
    protein_matches_file = 'plasmidplots_temp/protein_matches.txt'
    open(protein_matches_file, 'w').close()

    # Create text file using sequence
    with open(short_sequence_file, 'w') as ss_file:
        ss_file.write(">short_sequence\n")
        ss_file.write(sequence)

    # Run FASTX using variables
    pputil.fasta(short_sequence_file, subgroup_file, protein_matches_file)

    # Find subgroup in FASTA output
    with open(protein_matches_file) as subgroup_text:
        first_line = subgroup_text.readline().strip()
        if "short_sequence" not in first_line:
            return ''
        subgroup = str(first_line.split()[1])
    return subgroup


def sequence_trim(plasmid, start, end, dna_file):
    """
    Returns trimmed sequence string from a plasmid

    Takes plasmid name and file to search through,
    along with starting/ending locations
    """

    # Return sliced portion of sequence
    untrimmed_sequence = pputil.sequence_finder(dna_file, plasmid)
    trimmed_sequence = untrimmed_sequence[start:end - 1]
    return trimmed_sequence


def file_to_dict(filename, dna_sequence_file, replen, subgroup_file_dict):
    """Converts FASTA output file to a dictionary for plotting"""

    data_dict = {}
    line_number = 0
    plasmid_count = 0
    with open(filename) as data:
        for line in islice(data, None):
            line_number += 1
            plasmid = line.split()[0]

            # If sequence does not exist in dictionary,
            # initialize with the plasmid type and length
            if plasmid not in data_dict:
                plasmid_length = replen[plasmid]
                data_dict[plasmid] = [plasmid_length]
                plasmid_count += 1

            # Pull data from relevant columns in FASTA output and create tuple
            family = line.split()[1]
            start = min(int(line.split()[6]), int(line.split()[7]))
            end = max(int(line.split()[6]), int(line.split()[7]))
            if family in subgroup_file_dict.keys():
                protein_sequence = (sequence_trim(plasmid, start,
                                                  end, dna_sequence_file))
                subgroup = subgroup_search(protein_sequence, family,
                                           subgroup_file_dict)
            else:
                subgroup = ''
            pf_tuple = (family, start, end, subgroup)

            # Add tuple to dictionary value for the plasmid
            data_dict[plasmid].append(pf_tuple)

            if line_number % 10 == 0:
                print("FASTA output lines processed: " + str(line_number))

        print("Total lines processed: " + str(line_number))
        print("Plasmid count: " + str(plasmid_count))
    return data_dict


def dict_to_plot(strain, data_dict, sequence_color_dict,
                 circular_plot_columns=5, legend='legend.png',
                 border=True, baseline_colors=None, dna_file=None):
    """Creates a plot for each key in dictionary and generates images"""

    circular_plot_list = []
    linear_plot_list = []
    temp_file = "./plasmidplots_temp/temp.png"
    plot_image_dir = "./plasmidplots_temp/plot_images/"

    if not os.path.exists(plot_image_dir):
        os.makedirs(plot_image_dir)

    # Plot based on plasmid type and sort into two lists
    for plasmid, data in data_dict.items():

        if dna_file != None and baseline_colors != None:
            # Set color scale for plot baselines
            if baseline_colors == 'gc':
                color_scale_dict = gc_content_dict(dna_file, plasmid)

            elif baseline_colors == 'gcskew':
                color_scale_dict = gc_skew_dict(dna_file, plasmid)

            else:
                color_scale_dict = None

            # If color scale exists, generate grayscale values
            if color_scale_dict != None:
                baseline_color_scale = {}

                plasmid_length = data[0]
                first = True

                # Find minimum and maximum values on decimal scale
                scale_min = 0
                scale_max = 0
                first_scale = True
                for decimal_scale in color_scale_dict.values():
                    if first_scale:
                        scale_min = decimal_scale
                        scale_max = decimal_scale
                        first_scale = False
                    else:
                        scale_min = min(scale_min, decimal_scale)
                        scale_max = max(scale_max, decimal_scale)

                for location, decimal_scale in color_scale_dict.items():
                    start, end = location.split('-')
                    start = int(start)
                    end = int(end)

                    # Set window on first key/value pair
                    if first:
                        window = end - start + 1
                        first = False

                    # Increase width of gray bar to prevent artifacts
                    if end + window < plasmid_length:
                        end += window
                    else:
                        end = plasmid_length

                    # Add data for gray bar to dictionary
                    location = str(start) + '-' + str(end)
                    baseline_color_scale[location] = pputil.decimal_to_rgb_gray(decimal_scale,
                                                 lower=scale_min,
                                                 upper=scale_max)

            else:
                baseline_color_scale = None


        else:
            baseline_color_scale = None

        if 'cp' in plasmid.split('_')[1]:
            # Plot
            circular_plot(plasmid, data, sequence_color_dict,
                          baseline_custom_colors=baseline_color_scale)
            image = Image.open(temp_file)
            image.load()
            circular_plot_list.append(image)

        else:
            # Plot
            linear_plot(plasmid, data, sequence_color_dict,
                        baseline_custom_colors=baseline_color_scale)
            image = Image.open(temp_file)
            image.load()
            linear_plot_list.append(image)

    # Combine circular plot images
    circular_plots = imt.image_grid(circular_plot_list, circular_plot_columns)

    if circular_plots == None:
        print("No circular plasmids found for " + strain + ".")
        circular_plot_file = None
    else:
        circular_plot_file = plot_image_dir + strain + "_circular_plots.png"
        circular_plots.save(circular_plot_file)
        circular_plots.close()

    # Combine linear plot images
    linear_plots = imt.image_grid(linear_plot_list, 1)

    if linear_plots == None:
        print("No linear plasmids found for " + strain + ".")
        linear_plot_file = None
    else:
        linear_plot_file = plot_image_dir + strain + "_linear_plots.png"
        linear_plots.save(linear_plot_file)
        linear_plots.close()

    for image in (circular_plot_file, linear_plot_file):
        if image != None:
            imt.append_legend(image, legend, strain)

            if border:
                # Add black border and then white buffer around it
                imt.add_border(image)
                imt.add_border(image, border_color='#FFFFFF',
                               border_dimensions=(5, 5))

    # Return image file names
    return circular_plot_file, linear_plot_file


def strain_sort(data_dict):
    """
    Sorts a dictionary by strain

    The input dictionary should consist of plasmid keys and data for values
    The sorted dictionary with have strains for the keys and plasmid/data
    dictionaries as the values.
    """

    # Separate plasmids into groups
    strain_dict = {}
    for plasmid, data in data_dict.items():
        strain = plasmid.split('_')[0]
        # Create blank entry if none exists
        if strain not in strain_dict.keys():
            strain_dict[strain] = {}

        strain_dict[strain][plasmid] = data

    return strain_dict


# -------- End functions -------- #

# -------- Begin main program -------- #

def main(url_input_file, protein_input, color_file,
         subgroup_list_file, baseline_scale):
    # Set constants
    LEGEND_FONT_SIZE = 48
    TIMER_FORMAT = "%.1f"
    TIME_STRING = " Time: %s seconds"

    # Create folder for temporary files
    temp_dir = "./plasmidplots_temp/"
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    # Get NCBI urls
    url_list = ncbi.url_input(url_input_file)

    # Take filepath input
    # Can just be file name if in content folder (e.g. foo.txt)
    dna_sequence_file = temp_dir + 'replicons.txt'
    if protein_input == None:
        protein_input = 'pfam.txt'
    if color_file == None:
        color_file = 'colors.txt'

    # Start timer for entire program's run time
    program_start = timer()

    # Read data on colors and files for family subgroups
    sequence_color_dict, subgroup_file_dict = read_colors(color_file)

    # Generate legend for plots
    start_time = timer()

    legend_image_file = temp_dir + 'legend.png'
    imt.generate_legend(sequence_color_dict,
                        LEGEND_FONT_SIZE,
                        legend_image_file)

    end_time = timer()
    time = TIMER_FORMAT%(end_time - start_time)
    print("Legend generated." + TIME_STRING%time)

    # Add all plasmid sequences from each url to replicons.txt
    start_time = timer()

    ncbi_id_dict = ncbi.ncbi_scrape(url_list)

    end_time = timer()
    time = TIMER_FORMAT%(end_time - start_time)
    print("Plasmid IDs found." + TIME_STRING%time)

    start_time = timer()

    ncbi.sequence_download(ncbi_id_dict)

    end_time = timer()
    time = TIMER_FORMAT%(end_time - start_time)
    print("Sequences downloaded." + TIME_STRING%time)

    # Run FASTA
    start_time = timer()

    fasta_output = 'output.txt'
    pputil.fasta(dna_sequence_file, protein_input, fasta_output)

    end_time = timer()
    time = TIMER_FORMAT%(end_time - start_time)
    print("Initial FASTA run complete." + TIME_STRING%time)


    # Get sequence length of each record in "replicons.txt"
    start_time = timer()

    sequence_file = temp_dir + 'replicons.txt'
    replen = {}
    with open(sequence_file, 'r') as fna:
        for rec in SeqIO.parse(fna, 'fasta'):
            replen[rec.id] = len(rec.seq)

    end_time = timer()
    time = TIMER_FORMAT%(end_time - start_time)
    print("Plasmid length checking complete." + TIME_STRING%time)

    # Create dictionary from FASTA output
    start_time = timer()

    data_dict = file_to_dict(fasta_output, dna_sequence_file,
                             replen, subgroup_file_dict)

    end_time = timer()
    time = TIMER_FORMAT%(end_time - start_time)
    print("Data dictionary generation complete." + TIME_STRING%time)

    # Separate plasmids into groups
    start_time = timer()

    sorted_dict = strain_sort(data_dict)

    end_time = timer()
    time = TIMER_FORMAT%(end_time - start_time)
    print("Data dictionary sorted." + TIME_STRING%time)

    # Plot all strains
    start_time = timer()

    image_list = []
    for strain, data in sorted_dict.items():
        circular, linear = dict_to_plot(strain, data, sequence_color_dict,
                                        5, border=True,
                                        baseline_colors=baseline_scale,
                                        dna_file=dna_sequence_file,
                                        legend=legend_image_file)
        print("Strain plotted: " + strain)

        if circular != None:
            image_list.append(circular)
        if linear != None:
            image_list.append(linear)

    imt.images_to_pdf(image_list, 'plots.pdf')
    print("PDF generated.")

    end_time = timer()
    time = TIMER_FORMAT%(end_time - start_time)
    print("Plotting complete." + TIME_STRING%time)

    # Clean up temporary files
    start_time = timer()

    #TODO

    end_time = timer()
    time = TIMER_FORMAT%(end_time - start_time)
    print("Cleanup complete." + TIME_STRING%time)


    # Print total run time after input
    program_end = timer()
    time = TIMER_FORMAT%(program_end - program_start)
    print("Program run complete." + TIME_STRING%time)


# -------- Run -------- #
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plots plasmids')
    parser.add_argument("urls",
                        help="Text file with NCBI urls \
                        to download plasmids from")
    parser.add_argument("protein",
                        help="FASTA file with sequences of \
                        protein families to search for")
    parser.add_argument("colors",
                        help="Text file with protein families \
                        and hex value colors")
    parser.add_argument("subgroups",
                        help="Text file with protein families \
                        and names of text files with subfamily sequences")
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s ' + __version__)
    parser.add_argument("-b", "--baseline", type=str,
                        choices=["gc", "gcskew"],
                        help="Data to be used for baseline color scale")

    args = parser.parse_args()

    url_input_file = args.urls
    protein_input_file = args.protein
    color_file = args.colors
    subgroup_list_file = args.subgroups
    baseline_scale = args.baseline

    main(url_input_file, protein_input_file, color_file,
         subgroup_list_file, baseline_scale)