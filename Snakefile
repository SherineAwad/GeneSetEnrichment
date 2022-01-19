configfile: "config.yaml"

OUTNAME = config['OUTNAME']
ORGANISM = config['ORGANISM'] 
COMPARISON = config['COMPARISON']
rule all: 
    input:
         expand("{outname}-{type}-KEGG.csv", outname = OUTNAME, type = COMPARISON), 
         expand("{outname}-{type}-GO.csv", outname = OUTNAME, type = COMPARISON)

rule pathway:
    input: 
       config['DGE']
    output:
       expand("{outname}-{type}-KEGG.csv", outname = OUTNAME, type = COMPARISON), 
    params: 
        organims = ORGANISM,
        outname = OUTNAME, 
        comparison = COMPARISON
    shell: 
        """
        Rscript scripts/pathway.R {input} {params.organims} {params.comparison} {params.outname}
        """


rule GO:
    input:
       config['DGE']
    output:
       expand("{outname}-{type}-GO.csv", outname = OUTNAME, type = COMPARISON),
    params:
        organims = ORGANISM,
        outname = OUTNAME,
        comparison = COMPARISON
    shell:
        """
        Rscript scripts/GO.R {input} {params.organims} {params.comparison} {params.outname}
        """
