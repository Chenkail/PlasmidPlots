"""Tools for plasmidplots"""

# Import libraries
# Built-in
from itertools import islice
import os
import re
import shutil

# 3rd party
from bs4 import BeautifulSoup
import pexpect
from pyvirtualdisplay import Display
from selenium import webdriver
from selenium.webdriver.common.desired_capabilities import DesiredCapabilities
from selenium.webdriver.firefox.options import Options
import urllib


def get_dynamic_content(url, is32bit=False):
    """
    Loads and returns html from a dynamically generated webpage
    """

    if is32bit:
        gecko_bit = "linux32"
    else:
        gecko_bit = "linux64"

    # Download geckodriver
    geckodriver_version = "0.24.0"
    geckodriver_file = ("geckodriver-v" + geckodriver_version + "-"
                        + gecko_bit + ".tar.gz")
    command = ('bash -c "if ! [ -e ./geckodriver ]; then wget -q \
               https://github.com/mozilla/geckodriver/releases/download/v'
               + geckodriver_version + '/' + geckodriver_file
               + ' && tar -xzf ' + geckodriver_file + '; fi"')
    pexpect.run(command)

    display = Display(visible=0, size=(1024, 768))
    display.start()

    cap = DesiredCapabilities().FIREFOX
    cap["marionette"] = True

    options = Options()
    options.set_headless(headless=True)

    browser = webdriver.Firefox(capabilities=cap,
                                executable_path="./geckodriver",
                                firefox_options=options)
    browser.get(url)
    html = browser.page_source

    browser.close()
    display.stop()

    return html


def all_strain_list(url):
    """
    Returns a list of genome urls for every available strain of a given species

    Input should be a url similar to the following:
    https://www.ncbi.nlm.nih.gov/genome/genomes/16651

    Note: will not see any strains beyond the first 100 (NCBI limit per page)
    """

    # Get html and parse
    html = get_dynamic_content(url)
    soup = BeautifulSoup(html, 'html.parser')

    # Find header row
    table = soup.find('table', attrs={'id': 'genomestbl_id'})
    head = table.find('thead')
    header_list = []
    head = head.text.split()
    for item in head:
        header = item.strip()
        header_list.append(header)

    # Get table body rows
    body = table.find('tbody')
    rows = body.find_all('tr')

    url_set = set([])

    for row in rows:
        attributes = {'href': re.compile("/genome/\w*\?genome_assembly_id=")}
        for link in row.find_all('a', attrs=attributes):
            url_end = link.get('href')
            genome_url = "https://www.ncbi.nlm.nih.gov" + url_end
            url_set.add(genome_url)

    return url_set


def url_input(input_file):
    """
    Takes urls and returns a list

    The urls can be input directly, or in the form of a text file
    with each url on a separate line.
    """

    url_list = set([])
    with open(input_file) as url_file:
        for line in islice(url_file, None):
            url = line.strip()
            # Test for empty lines
            if url != '':
                # If url points to strain table for a species,
                # add urls for each strain to list
                if "/genome/genomes/" in line:
                    temp_list = all_strain_list(url)
                    for item in temp_list:
                        url_list.add(item)

                else:
                    url_list.add(line)

    # Count urls in list
    url_count = 0
    for item in url_list:
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

    # Find GCA and gene assembly version
    table = soup.find('table', attrs={'class': 'summary'})
    a = table.find_all('a', attrs={'target': '_blank'})
    for item in a:
        line = item.text
        if "GCA_" in line:
            gca_line = line
            break

    split_line = gca_line.split()
    gca_id = split_line[0]
    asm_id = split_line[1]

    # Find strain given GCA ID
    asm_version = asm_id.split('v')[1]
    gcf = gca_id.split('.')[0] + "." + asm_version
    gcf = gcf.replace("GCA", "GCF")
    assembly_url = "https://www.ncbi.nlm.nih.gov/assembly/" + gcf
    asm_html = urllib.request.urlopen(assembly_url)

    # Parse html
    asm_soup = BeautifulSoup(asm_html, 'html.parser')

    # Find organism name
    dlclass = asm_soup.find('dl',
                            attrs={'class':'assembly_summary_new margin_t0'})
    organism = dlclass.find('dd').text
    organism = organism.split('(')[0].strip()

    # Find strain IDs
    dd = dlclass.find_all('dd')
    for item in dd:
        line = item.text
        if "Strain:" in line:
            strain_text = line
            break

    strain = strain_text.split()[1]

    # Remove strain name and any additional formatting from organism name
    if strain in organism:
        organism = organism.split(strain)[0].strip()

    organism = organism.replace(' ', '-')

    # Add strain to organism name
    name = organism + "--" + strain

    return name


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
            print("Warning: No plasmids found for strain "
                  + strain + " (URL: " + url + ")")

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
        command = ('bash -c ' + '"wget -q -O ' + temp1
                    + ' "https://www.ncbi.nlm.nih.gov/search/api/sequence/'
                    + sequence_id + '?report=fasta"')
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
