def url_input():
    """
    Takes urls and returns a list
    
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
    
    # Import libraries
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