[![Snakemake](https://img.shields.io/badge/snakemake-≥6.0.2-brightgreen.svg)](https://snakemake.github.io)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)


Snakemake Workflow for Gene Set Enrichment: Pathway and GO Analysis   
=================================================================================

Author: Sherine Awad 

Change the config.yaml file appropriately according to your data. 
Update parameters of Genrich in the config file. Also, change workdir where the reference genome,  etc.



Then run: snakemake -jnumber_of_cores, for example for 5 cores use:

    snakemake -j5 

and for a dry run use: 

    snakemake -j1 -n 


and to print the commands in a dry run use:

    snakemake -j1 -n -p 

For the sake reproducibility, use conda to pull same versions of tools. Snakemake and conda have to be installed in your system:

    snakemake --cores --use-conda


## TO DO 

1. Add more Genome for enrichment 
2. Add more parameters for a genric pipeline 



