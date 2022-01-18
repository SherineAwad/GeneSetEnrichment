configfile: "config.yaml"

OUTNAME = config['OUTNAME']
ORGANISM = config['ORGANISM'] 
print(OUTNAME)
print(ORGANISM)

rule all: 
    input:
         expand("{outname}-KEGG.csv", outname = OUTNAME), 
         expand("{outname}-GO.csv", outname = OUTNAME)

rule pathway:
    input: 
       config['DGE']
    output:
       expand("{outname}-KEGG.csv", outname = OUTNAME) 
    params: 
        organims = ORGANISM,
        outname = OUTNAME 
    shell: 
        """
        Rscript scripts/pathway.R {input} {params.organims} {params.outname}
        """


rule GO:
    input:
       config['DGE']
    output:
       expand("{outname}-GO.csv", outname = OUTNAME)
    params:
        organims = ORGANISM,
        outname = OUTNAME
    shell:
        """
        Rscript scripts/GO.R {input} {params.organims} {params.outname}
        """
