# Manual

To launch the program, open the terminal, move to the folder you have downloaded, and write the following command line: `Phyton3 CODENAME.py`.

The script automatically creates the "output" folder inside the same directory. If it already exists, the script raises an error but it continues to work. Therefore, we may suggest to rename or remove previous files before continuing. 

At this point, you are required to give as input the names of the files (including the format) you would like to work on. In case they are not in the same folder, you need to provide the entire file paths. 

If you insert the wrong file names, the code will ask you to re-input the correct names.

Remember to check for NaN values or duplicates, otherwise the program stops, asking you to correct them (then, you need to restart the script from the beginning). 

If everything is correct, the code prints out the number of the total exchanges as follow: 

`... In the community there are n exchanges of m different compunds`.

Meanwhile, it creates three general files describing the entire community:

-  compunds_exchanged.tsv

-  donors_for_compund.tsv

-  receivers_for_compound.tsv

Now you are required to give as input the compound name. You can both use the extended name or the biggID (see for reference: http://bigg.ucsd.edu/universal/metabolites). 

In case the compound is not found in the smetana, the script warns you. We may suggest to check for spelling errors or the correct biggID. If the error persists, the compound is probably not exchanged in the community. 

Alternatively, if the compound is exchanged in the community, you see on screen the number of total exchanges and a browser page opens showing the final network.

The network is interactive, hence you can move as you like all the present nodes. Using ctrl you can select more than one node. We may remember that the color palette is indicative of the phylum level, the diameters of the nodes are proportional to MAGs abundances, and the thickness of the edges is proportional to smetana value. A small info box appears as soon as you position your mouse either on a node or an edge.
 
You can zoom in and zoom out. 

Finally, the script creates two compound specific files: 

-  compundName_exchanges.tsv

-  compundName_species_behavior.tsv

These two files, as well as the network, can be found inside a specific directory named after the compound biggID: '.../CODENAME/outputs/compoundID'

As you finish reading all the files, you can either ask for another compound by writing "Y", or stop the program writing "N".

