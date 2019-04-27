"""Tools for plasmidplots"""

# Import libraries
import os
import pexpect
import shutil
import urllib
from bs4 import BeautifulSoup
from itertools import islice

def url_input(input_file):
    """
    Takes urls and returns a list

    The urls can be input directly, or in the form of a text file with each url on a separate line.
    """

    url_list = []
    url_count = 0
    with open(input_file) as url_file:
        for line in islice(url_file, None):
            if line not in url_list and line.strip() != '':
                url_list.append(line)
                url_count += 1

    print("URL count: " + str(url_count))
    return url_list


def strain_name_scrape(url):
    """
    Scrapes NCBI for name of a strain given its url
    """

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

    temp_dir = "./plasmidplots_temp/"

    # Name files
    temp1 = temp_dir + 'temp1.txt'
    temp2 = temp_dir + 'temp2.txt'

    # This will be the name of the output file
    sequence_file = temp_dir + 'replicons.txt'

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
