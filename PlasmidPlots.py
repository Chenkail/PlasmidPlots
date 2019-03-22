## Use matplotlib to draw the plasmids in replicons.fna / output.txt.


# -------- Begin functions -------- #

def url_input():
    """
    Takes urls and returns a list.
    
    The urls can be input directly, or in the form of a text file with each url on a separate line.
    """
    
    # Import libraries
    from itertools import islice
    
    url_list = []
    url_count = 0
    next_input = input("Enter an NCBI strain url to grab plasmids from, or a file containing NCBI urls: ").strip()
    while next_input != '':
        # If input is a url, add it to the list
        if 'ncbi.nlm.nih.gov/genome' in next_input:
            if next_input not in url_list:
                url_list.append(next_input)
                url_count += 1
        # Otherwise, assume input is a file and add each line
        else: 
            with open(next_input) as url_file:
                for line in islice(url_file, None):
                    if line not in url_list:
                        url_list.append(line)
                        url_count += 1
        
        
        next_input = input("Enter an NCBI strain url to grab plasmids from, or a file containing NCBI urls: ").strip()
    
    print("URL count: " + str(url_count))
    return url_list

def strain_name_scrape(url):
    """
    Scrapes NCBI for name of a strain given its url
    """
    
    import urllib
    from bs4 import BeautifulSoup
    
    # Get html
    html = urllib.request.urlopen(url)

    # Parse html
    soup = BeautifulSoup(html, 'html.parser')
    
    # Find GCA
    table = soup.find('table', attrs={'class': 'summary'})
    a = table.find_all('a', attrs={'target': '_blank'})
    for item in a:
        line = item.text
        if "GCA_" in line:
            gca_line = line
            break

    split_line = gca_line.split()
    gca_id = split_line[0]
    
    # Find strain given GCA ID
    asm_cleaned = gca_id.split('.')[0] + ".1"
    asm_id = asm_cleaned.replace("GCA", "GCF")
    asm_url = "https://www.ncbi.nlm.nih.gov/assembly/" + asm_id
    asm_html = urllib.request.urlopen(asm_url)
    
    # Parse html
    asm_soup = BeautifulSoup(asm_html, 'html.parser')
    
    # Find strain ID
    dlclass = asm_soup.find('dl', attrs={'class': 'assembly_summary_new margin_t0'})
    dd = dlclass.find_all('dd')
    
    # Search for strain ID in items
    for item in dd:
        line = item.text
        if "Strain:" in line:
            strain_text = line
            break
    
    split_text = strain_text.split()
    strain = split_text[1]
    return strain



def ncbi_scrape(url_list):
    """
    Scrapes tables for lists of IDs and names given urls for each strain
    
    Takes input in the form of a list of NCBI links, e.g.
    https://www.ncbi.nlm.nih.gov/genome/738?genome_assembly_id=335284
    https://www.ncbi.nlm.nih.gov/genome/738?genome_assembly_id=300340
    
    Enter a blank line to stop input.
    """
    
    # Import libraries
    import urllib
    from bs4 import BeautifulSoup
    
    
    # Loop over each url in the list and add data to dictionary
    ncbi_dict = {}
    for url in url_list:
        # Get html
        html = urllib.request.urlopen(url)
        
        # Parse html
        soup = BeautifulSoup(html, 'html.parser')
        
        # Grab table from site
        table = soup.find('table', attrs={'class': 'GenomeList2'})
        
        # Only take body from table
        body = table.find('tbody')
        
        # Get strain for url
        strain = strain_name_scrape(url)
        
        # Get plasmid names and IDs and add to dictionary
        for rows in body.find_all('tr', attrs={'align': 'center'}):
            no_plasmids_found = True
            
            data = rows.find_all('td')
            gene_type = data[0].text
            # Only read data for plasmids, skip the main sequence
            if gene_type == 'Plsm':
                no_plasmids_found = False
                
                # Get name of plasmid
                name = data[1].text
                
                # Clean up plasmid name formatting and add strain
                name = name.strip()
                if ' ' in name:
                    name = name.replace(' ', '_')
                if '_' in name:
                    name = name.split('_')[1]
                name = strain + "_" + name
                
                insdc = data[3].text
                ncbi_dict[name] = insdc
            
        if no_plasmids_found:
            print("Warning: No plasmids found for strain " + strain + " (URL: " + url + ")")
    
    return ncbi_dict


def sequence_download(id_dict):
    """Downloads sequences given a dictionary of names and IDs"""
    
    # Import libraries
    import os
    import pexpect
    import shutil
    
    # Name files
    temp1 = 'temp1.txt'
    temp2 = 'temp2.txt'
    
    # This will be the name of the output file
    sequence_file = 'replicons.txt'
    
    # Delete old file if it exists
    if os.path.isfile(sequence_file):
        os.remove(sequence_file)
    
    # Loop over each plasmid
    for plasmid_name, sequence_id in id_dict.items():
        # Download sequence
        command = 'bash -c ' + '"wget -q -O ' + temp1 + ' "https://www.ncbi.nlm.nih.gov/search/api/sequence/' + sequence_id + '?report=fasta"'
        pexpect.run(command)
        
        # Label with plasmid name
        with open(temp1) as in_file:
            in_file.readline()
            line = ">" + plasmid_name + "\n"
            
            with open(temp2, 'w') as out_file:
                out_file.write(line)
                shutil.copyfileobj(in_file, out_file)
        
        # Add to file
        with open(temp2) as in_file:
            with open(sequence_file, 'a') as out_file:
                for line in in_file:
                    out_file.write(line)
        
        # Clean up temporary files
        temp_files = [temp1,
                      temp2,]

        for file in temp_files:
            if os.path.isfile(file):
                os.remove(file)


def fasta(dna_input, protein_input, output_file):
    """
    Runs FASTA given two filepaths and the name of the file to output to.
    
    The first should be a DNA sequence file, and the second should be 
    a protein sequence (amino acid) file.
    """
    
    # Import libraries
    import pexpect
    
    # Download FASTA
    command = 'bash -c "if ! [ -e /content/fasta-36.3.8g/bin ]; then wget -q https://github.com/wrpearson/fasta36/releases/download/fasta-v36.3.8g/fasta-36.3.8g-linux64.tar.gz && tar -xzf fasta-36.3.8g-linux64.tar.gz; fi"'
    pexpect.run(command)
    
    # Run FASTX using variables
    command = 'bash -c ' + '"if [ \':$PATH:\' != *\':/content/fasta-36.3.8g/bin/fasta36:\'* ]; then PATH=\'$PATH:/content/fasta-36.3.8g/bin\'; fi && fastx36 -m 8 -E 0.05 ' + dna_input + ' ' + protein_input + ' > ' + output_file + '"'
    pexpect.run(command)


def read_colors(file):
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
    
    # Import libraries
    from itertools import islice
    
    color_dict = {}
    subgroup_dict = {}
    
    # Ask user if they would like to define subgroups, otherwise use default
    define_subgroups = False
    define_subgroups_input = input("Would you like to define subgroup files for each protein family? (y/n) ")
    if 'y' in define_subgroups_input:
        define_subgroups = True
    else:
        subgroup_dict = {'32':'pf32.txt'}
    
    with open(file) as color_file:
        for line in islice(color_file, None):
            family, color = line.strip().split(":")
            
            # Add color to dictionary
            color_dict[family] = color
            
            # Prompt user for file corresponding to protein family
            if family not in subgroup_dict.keys():
                subgroup_file_name = ''
                if define_subgroups and family != 'Base Color':
                    subgroup_file_name = input("Enter filepath for Family " + family + " protein sequences: ")
                
                subgroup_dict[family] = subgroup_file_name
    
    if 'Base Color' not in color_dict.keys():
        color_dict['Base Color'] = '#C0C0C0'
    
    return color_dict, subgroup_dict


def pil_grid(images, columns):
    """Takes an array of images and a column count, and generates a grid using the images"""
    
    # Import libraries
    from PIL import Image
    import numpy as np
    
    # Count the images and set row/column count
    image_count = len(images)
    
    # If the number of images is less than the number of columns,
    # set the number of columns to the number of images
    column_count = min(image_count, columns)
    
    # Set row count
    rows = image_count // column_count
    if image_count % column_count != 0:
        rows += 1
    
    # Create arrays of 0s
    left_x, top_y = [0] * column_count, [0] * rows
    
    # Calculate max image width and height for each row/column
    for index, image in enumerate(images):
        horizontal_index, vertical_index = index % column_count, index // column_count
        left_x[horizontal_index] = max(left_x[horizontal_index], image.size[0])
        top_y[vertical_index] = max(top_y[vertical_index], image.size[1])
    
    # Set column width and row height for each row/column
    left_x, top_y = np.cumsum([0] + left_x), np.cumsum([0] + top_y)
    
    # Create blank image
    image_grid = Image.new('RGB', (left_x[-1], top_y[-1]), color='white')
    
    # Add images to grid
    for index, image in enumerate(images):
        column = index % column_count
        row = index // column_count # Note: Max height for last row is irrelevant because there is unlimited space below
        
        if columns == 1:
            x_offset = 0
        else:
            x_offset = (left_x[column + 1] - left_x[column])/2 - image.size[0]/2
        
        x = left_x[column] + x_offset
        y_offset = (top_y[row + 1] - top_y[row])/2 - image.size[1]/2
        if y_offset < 0:
            y_offset = 0
        y = top_y[row] + y_offset
        
        x = int(x)
        y = int(y)
        image_grid.paste(image, (x, y))
    
    return image_grid


def generate_legend(color_dict, font_size, file_name='legend.png'):
    """Creates a legend given a dictionary of colors and names"""
    
    # Import libraries
    import copy
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt
    
    # Remove base color from legend
    legend_colors = copy.deepcopy(color_dict)
    legend_colors.pop('Base Color')
    
    # Create list of patches with colors and names
    patch_list = []
    for label, color in sorted(legend_colors.items()):
        # Create patch and add to list
        legend_key = mpatches.Patch(label=label, color=color)
        patch_list.append(legend_key)
    
    # Plot and save to image
    plt.legend(fontsize=font_size, handles=patch_list)
    plt.axis('off')
    plt.savefig(file_name, bbox_inches='tight')
    plt.close('all')


def append_legend(image, legend_image, legend_size=(300, 300), legend_location='bottom right'):
    """Adds legend to an image"""
    
    from PIL import Image
    
    # Get dimensions of plot image
    old_image = Image.open(image)
    plot_width, plot_height = old_image.size
    
    # Resize legend
    legend = Image.open(legend_image)
    legend.thumbnail(legend_size)
    legend_width, legend_height = legend.size
    
    legend_x = plot_width - legend_width
    legend_y = plot_height
    
    # New image
    new_image = Image.new('RGB', (plot_width, plot_height + legend_height), (255, 255, 255))
    new_image.paste(old_image, (0, 0))
    old_image.close()
    new_image.paste(legend, (legend_x, legend_y), legend)
    legend.close()
    new_image.save(image)


def linear_plot(plasmid, data, sequence_color_dict):
    """Plots a linear plasmid given the plasmid name and data list"""
    
    # Import libraries
    import matplotlib.pyplot as plt
    
    # Set constants
    HORIZONTAL_SCALE_CONSTANT = 1/4000
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
    baseline = linear.barh(0, plasmid_length, height=1, color=sequence_color_dict['Base Color'])
    
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
                            arrowprops=dict(facecolor='black', arrowstyle='-',))
    
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
    plt.text(label_x, label_y, label_text)
    
    # Fix axes limits
    xmin, xmax, ymin, ymax = plt.axis('tight')
    plt.ylim(0, 1)
    
    # Hide background and save plot to image
    plt.axis('off')
    plt.savefig('temp.png', bbox_inches='tight')
    plt.close('all')


def circular_plot(plasmid, data, sequence_color_dict):
    """Plots a circular plasmid given the name and the data list"""
    
    # Import libraries
    import matplotlib.pyplot as plt
    import numpy as np
    
    # Set constants
    CHART_BOTTOM = 4
    CHART_THICKNESS = 1.5
    CIRCULAR_SCALE_CONSTANT = 1.3
    
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
    circle.set_theta_offset(np.radians(90)) # Move origin to top instead of right
    baseline = circle.bar(0, 
                          circle_width, 
                          width=-np.radians(360), 
                          bottom=circle_bottom, 
                          align='edge', 
                          color=sequence_color_dict['Base Color'],)
    
    # Use data from first protein family sequence for initial minimum/maximum values
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
    
    center = (minimum + maximum)/2
    
    # Plot
    for index in range(1, data_count):
        sequence = data[index]
        # Get sequence data
        sequence_start = min(sequence[1], sequence[2])
        sequence_end = max(sequence[1], sequence[2])
        sequence_length = sequence_end - sequence_start
        sequence_family = str(sequence[0])

        # Plot gene onto plasmid
        gene_plot_start = ((center - sequence_start) / plasmid_length) * np.radians(360)
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
            
            # Rotate label to point towards center of circle and
            # make the bottom side of the text face the bottom,
            # no matter which side of the plot it is on
            if np.degrees(line_angle) < 180:
                label_angle = np.degrees(line_angle + np.pi*3/2)
            else:
                label_angle = np.degrees(line_angle + np.pi/2)
            # Add annotation to plot
            circle.annotate(subgroup, 
                            xy=(line_angle, 4), 
                            xytext=(line_angle, 7), 
                            xycoords='data', 
                            ha='center', 
                            va='center', 
                            fontsize=9, 
                            rotation=label_angle, 
                            arrowprops=dict(facecolor='black', arrowstyle='-',),)
    
    # Label - plasmid name and length
    label_x = 0
    label_y = 0
    label_text = str(sequence_name) + "\n" + str(plasmid_length) + " base pairs"
    plt.text(label_x, label_y, label_text, ha='center')
    
    # Hide background and save plot to image
    plt.axis('off')
    plt.savefig('temp.png', bbox_inches='tight')
    plt.close('all')


def subgroup_search(sequence, family, subgroup_dict):
    """Find best subgroup match given a short sequence and the family to search"""
    
    # Import libraries
    import os
    
    # Setup
    subgroup_file = subgroup_dict[family]
    short_sequence_file = "short_sequence.txt"
    
    # Create blank file
    temp_subgroup_text_file = 'subgroup.txt'
    open(temp_subgroup_text_file, 'w').close()
    
    # Create text file using sequence
    with open(short_sequence_file, 'w') as ss_file:
        ss_file.write(">short_sequence\n")
        ss_file.write(sequence)
    
    # Run FASTX using variables
    fasta(short_sequence_file, subgroup_file, temp_subgroup_text_file)
    
    # Find subgroup in FASTA output    
    with open(temp_subgroup_text_file) as subgroup_text:
        first_line = subgroup_text.readline().strip()
        if "short_sequence" not in first_line:
            return ''
        subgroup = str(first_line.split()[1])
    return subgroup
    
    # Clean up temporary files
    temp_files = [short_sequence_file,
                  temp_subgroup_text_file,]
    
    for file in temp_files:
        if os.path.isfile(file):
            os.remove(file)


def short_protein_sequence_search(plasmid, start, end, dna_input):
    """Returns trimmed sequence string given the plasmid name to cut and starting/ending locations"""
    
    # Import libraries
    from itertools import islice
    
    # Variable setup
    line_number = 0
    sequence = ''
    recording = False
    
    with open(dna_input) as dna_file:
        # Move to line identifying the plasmid being searched for
        for line in islice(dna_file, None):
            line_number += 1
            if recording:
                if '>' not in line:
                    sequence += line.strip()
                else:
                    break

            if plasmid in line:
                recording = True
    
    # Return sliced portion of sequence
    untrimmed_sequence = sequence.strip().replace(" ", "")
    trimmed_sequence = untrimmed_sequence[start:end - 1]
    return trimmed_sequence


def file_to_dict(filename, dna_input, replen, subgroup_dict):
    """Converts FASTA output file to a dictionary for plotting"""
    
    from itertools import islice
    
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
            protein_sequence = short_protein_sequence_search(plasmid, start, end, dna_input)
            if family in subgroup_dict.keys():
                subgroup = subgroup_search(protein_sequence, family, subgroup_dict)
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


def dict_to_plot(strain, data_dict, sequence_color_dict, circular_plot_columns, legend='legend.png'):
    """Loops over every key in dictionary and creates a plot for each, then generates images"""
    
    # Import libraries
    import os
    from PIL import Image
    
    circular_plot_list = []
    linear_plot_list = []
    temp_file = 'temp.png'
    for plasmid, data in data_dict.items():
        if ('cp' in plasmid):
            circular_plot(plasmid, data, sequence_color_dict)
            image = Image.open(temp_file)
            image.load()
            circular_plot_list.append(image)
            
        else:
            linear_plot(plasmid, data, sequence_color_dict)
            image = Image.open(temp_file)
            image.load()
            linear_plot_list.append(image)
    
    circular_plots = pil_grid(circular_plot_list, circular_plot_columns)
    circular_plot_file = strain + "_circular_plots.png"
    circular_plots.save(circular_plot_file)
    circular_plots.close()
    
    linear_plots = pil_grid(linear_plot_list, 1)
    linear_plot_file = strain + "_linear_plots.png"
    linear_plots.save(linear_plot_file)
    linear_plots.close()
    
    append_legend(circular_plot_file, legend)
    append_legend(linear_plot_file, legend)
    
    # Clean up temporary files
    temp_files = [temp_file,]
    
    for file in temp_files:
        if os.path.isfile(file):
            os.remove(file)


def strain_sort(data_dict):
    """Returns a dictionary where the keys are the strains and the values are plasmid dictionaries"""
    
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

def main():
    # Import libraries
    from Bio import SeqIO
    from timeit import default_timer as timer
    import os
    
    # Set constants
    LEGEND_FONT_SIZE = 48
    TIMER_FORMAT = "%.1f"
    TIME_STRING = " Time: %s seconds"
    
    # Get NCBI urls
    url_list = url_input()

    # Take filepath input
    # Can just be file name if in content folder (e.g. replicons.fna)
    dna_input = 'replicons.txt'
    protein_input = input("Enter filepath for protein sequences: ").strip()
    if protein_input == '':
        protein_input = 'pfam.txt'
    color_file = input("Enter filepath for colors: ").strip()
    if color_file == '':
        color_file = 'colors.txt'
    
    # Start timer for entire program's run time
    program_start = timer()
    
    # Read data on colors and files for family subgroups
    sequence_color_dict, subgroup_dict = read_colors(color_file)
    
    # Generate legend for plots
    start_time = timer()
    
    legend_image_file = 'legend.png'
    generate_legend(sequence_color_dict, LEGEND_FONT_SIZE, legend_image_file)
    
    end_time = timer()
    time = TIMER_FORMAT%(end_time - start_time)
    print("Legend generated." + TIME_STRING%time)
    
    # Add all plasmid sequences from each url to replicons.fna
    start_time = timer()
    
    ncbi_id_dict = ncbi_scrape(url_list)
    
    
    end_time = timer()
    time = TIMER_FORMAT%(end_time - start_time)
    print("Plasmid IDs found." + TIME_STRING%time)
    
    start_time = timer()
    
    sequence_download(ncbi_id_dict)
    
    end_time = timer()
    time = TIMER_FORMAT%(end_time - start_time)
    print("Sequences downloaded." + TIME_STRING%time)
    
    # Run FASTA
    start_time = timer()
    
    fasta_output = 'output.txt'
    fasta(dna_input, protein_input, fasta_output)
    
    end_time = timer()
    time = TIMER_FORMAT%(end_time - start_time)
    print("Initial FASTA run complete." + TIME_STRING%time)
    
    
    # Get sequence length of each record in "replicons.txt"
    start_time = timer()
    
    sequence_file = 'replicons.txt'
    replen = {}
    with open(sequence_file, 'r') as fna:
        for rec in SeqIO.parse(fna, 'fasta'):
            replen[rec.id] = len(rec.seq)
    
    end_time = timer()
    time = TIMER_FORMAT%(end_time - start_time)
    print("Plasmid length checking complete." + TIME_STRING%time)
    
    # Create dictionary from FASTA output
    start_time = timer()
    
    data_dict = file_to_dict(fasta_output, dna_input, replen, subgroup_dict)
    
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
    
    for strain, data in sorted_dict.items():
        dict_to_plot(strain, data, sequence_color_dict, 5)
        print("Strain plotted: " + strain)
    
    end_time = timer()
    time = TIMER_FORMAT%(end_time - start_time)
    print("Plotting complete." + TIME_STRING%time)
    
    # Clean up temporary files
    start_time = timer()
    
    temp_files = [fasta_output,]
    for file in temp_files:
        if os.path.isfile(file):
            os.remove(file)
    
    end_time = timer()
    time = TIMER_FORMAT%(end_time - start_time)
    print("Cleanup complete." + TIME_STRING%time)
    
    
    # Print total run time after input
    program_end = timer()
    time = TIMER_FORMAT%(program_end - program_start)
    print("Program run complete." + TIME_STRING%time)


# Run main
if __name__ == "__main__":
    main()