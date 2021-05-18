<h1 align='center'> CODENAME </h1>

<p align='center'>
  Michela Palamin, Elisabetta Offer, and Alice Frisinghelli
</p>  
<p align='center'> 
  Institute of Molecular Biology
</p>
<p align='center'>
  University of Padua
</p>

<p align='center'>
  <i>This code was developed as part of a project carried out during the course of Microbial Metagenomics (Molecular Biology master degree) at the University of Padova. The project was supervised by Prof. Stefano Campanaro and Dr. Arianna Basile.</i>
</p>



Contents: 
1) INTRODUCTION
2) REQUIREMENTS
3) INPUT FILES
4) OUTPUT FILES
5) TIPS

______________________________________________________________________________________
## 1) INTRODUCTION

CODENAME (CompOunD Exchanges iN A MicrobiomE) is designed for giving graphical visualizations of metabolite exchanges within complex communities. Flux Balance analysis (FBA) results are not easy to interpret, and certainly not intuitive. CODENAME aims to provide a simple and useful visualization for managing high-throughput data from metagenomic and FBA analyses.
CODENAME directly work on SMETANA output files.

______________________________________________________________________________________
## 2) REQUIREMENTS

To run the script it is necessary to have Phyton v.3.x installed on your laptop. Moreover, you need the following libraries:

-  `Numpy v1.15.1`
-  `Pandas v0.23.4`
-  `NetworkX v2.1`
-  `PyVis v0.1.9`
-  `os`

_______________________________________________________________________________________
## 3) INPUT FILES

CODENAME needs four different input files in order to work: 

- **bigg_compunds_conversion_table_CORRECT.txt**: tabular file allowing the conversion from the biggID to the extended name of the compounds. This file is provided directly and can be found within the downloaded folder. 

The remaining three files must be provided by the user:

- the output file (tsv format) generated by `SMETANA v1.0.0`.
- a tabular file (txt format) which comprises the species (MAGs) abundances of the focus community generated using `checkM v1.0.12`. 
- a tabular file (txt format) reporting the MAGs taxonomy generated using `GTDB-Tk v1.3.0` and then converted to the NCBI taxonomy using `gtdb_to_ncbi_majority_vote.py`.

________________________________________________________________________________________
## 4) OUTPUT FILES

For the entire community:

-  **compounds_exchanged.tsv**: file reporting all compounds exchanged, each reported with the relative number of exchanges and the average probability of exchange (smetana_avg) (file sorted by the number of exchanges in a descending way).
-  **donors_for_compound.tsv**: file reporting all donors for each compound, each reported with the number of receivers and the average probability of donation.
-  **receivers_for_compound.tsv**: file reporting all receivers for each compound, each reported with the number of donors and the average probability of receiving it.

For each compound asked by the user, the script creates a new folder named after the biggID of the compound containing:

-  **compoundName_exchanges.tsv**: file reporting a subset of the smetana output only relative to that compound.
-  **compoundName_species_behaviour.tsv**: file reporting the characteristics of all species involved (relative abundance in %, taxonomy, number of receivers with the relative average probability, number of donors with the relative average probability) and the behaviour of the same depending on donor/receiver probability ratio (file sorted by relative abundances in a descenting way).
-  **compoundName_exchanges.html**: file with the interactive network relative to all exchanges of that compound in the community, where the colour of the nodes depends on the taxonomy, the size of the nodes is proportional to MAGs abundances and the thickness of the edges is proportional to smetana value.

_________________________________________________________________________________________
## 5) TIPS

In order to explore completely the community, we suggest the user to group the species given in the supplementary output file by their taxonomy. This may provide a more general view of the community. To do so, you can run the following code where `biggID` and `compoundName` are general references and must be substituted:

```python
result = pd.read_csv('.\outputs\biggID\compoundName_species_behaviour.tsv', delimiter = "\t", index_col='Species')
result.groupby(['taxonomy', 'behaviour']).mean()
```

Since the compound names in the smetana output are reported as “bigg_IDs” and not extended names, we provide a conversion file created by downloading the conversion file from [bigg site](http://bigg.ucsd.edu/universal/metabolites) and cleaning it for both NaN values and duplicates: **bigg_compunds_conversion_table_CORRECT.txt**.

In case you would like to update it, here we report the code which can be used to create the file:

```python
# download of the bigg compounds conversion file from the web (Internet connection is required)
import requests
import json

response = requests.get('http://bigg.ucsd.edu/api/v2/universal/metabolites')
data = response.text
res = json.loads(data)

# convert it to a Pandas DataFrame
df = pd.DataFrame.from_dict(res['results'])

# check for NA and keep only the first occurency
df.dropna(inplace=True)
# check for duplicates (only at index level) and remove them
df.drop_duplicates(inplace=True)

# rename columns to future managing 
df.rename(columns={'name':'extended name','bigg_id':'biggID'}, inplace=True)
id_conversion_table = df.set_index('extended name')

# save as id_conversion_table_CORRECT.txt file
path = os.getcwd()
id_conversion_table.to_csv(path+'/bigg_compounds_conversion_table_CORRECT.txt', sep = '\t')
```
 
