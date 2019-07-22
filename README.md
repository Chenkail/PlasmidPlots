![](https://img.shields.io/github/release/chenkail/plasmidplots.svg)
![](https://img.shields.io/github/license/Chenkail/plasmidplots.svg)
# PlasmidPlots

Borrelia burgdorferi, the bacteria responsible for Lyme disease, contains a large number of plasmids in its genome. These plasmids frequently contain copies of various gene families, some of which are essential for gene partitioning and replication. In *A bacterial genome in flux: the twelve linear and nine circular extrachromosomal DNAs in an infectious isolate of the Lyme disease spirochete Borrelia burgdorferi* [(full paper here)](https://onlinelibrary.wiley.com/doi/full/10.1046/j.1365-2958.2000.01698.x), Casjens et al. plot these genes on the strain B31 as seen below (Figure 8 in the paper):

![](https://wol-prod-cdn.literatumonline.com/cms/attachment/ca139954-6948-4d91-a567-c720e741c7f6/mmi_1698_f8.gif)

This program generates similar plots for any strain of any bacteria and any protein families, and can also show GC content or GC skew on the plots.

## Installing

### Prerequisites

#### Python
 - Python 3
 - [imagemergetools](https://github.com/Chenkail/imagemergetools)
 - matplotlib
 - pyvirtualdisplay
 - Biopython
 - PIL

#### Other
 - Firefox

### Pip
```
pip3 install plasmidplots
```

## Running
The program requires several text files to run, and a folder of example files can be found [here](https://drive.google.com/open?id=1oZpTRwPDGC5YsEs74oZMlLodlZzwWHXg). To run the program, open bash, navigate to the folder containing the files, and run the following command (this will take several minutes, depending on your computer):
```
python3 -m plasmidplots.pplots borrelia_urls.txt pfam.txt colors.txt subgroups.txt
```
To add GC content or GC skew to the plots, add `--baseline gc` or `--baseline gcskew`, respectively.
