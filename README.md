# calculus_study_Lena_2021
This repo contains scripts for reproducing results reported in the calculus study, Lena et al., 2021


~~~
git clone https://github.com/SegataLab/calculus_study_Lena_2021.git
~~~

Dependencies required:
* [numpy (>= 1.x)](https://numpy.org/install/)
* [matplotlib (>= 3.x)](https://matplotlib.org/stable/users/installing.html)
* [seaborn](https://seaborn.pydata.org/installing.html)
* [biopython](https://biopython.org/wiki/Download)
* [trimal (>= 1.4.1)](http://trimal.cgenomics.org/)

## Genome alignment post-processing

~~~Bash
usage: genome_aln_tailoring.py [-h] [-a] [-sa SEMI_AUTOMATED_TAILORING]
                               [-sl SHORT_LIST] [-m MANUAL_TAILORING]
                               [genome_alignment_file]
                               [tailored_genome_alignment_file]

positional arguments:
  genome_alignment_file
                        Input a genome alignment needs to be tailored. [fasta]
                        format
  tailored_genome_alignment_file
                        Output a tailored genome alignment in fasta format

optional arguments:
  -h, --help            show this help message and exit
  -a, --automated_tailoring
                        This option allows automated tailoring for minimizing
                        missing information in the alignment
  -sa SEMI_AUTOMATED_TAILORING, --semi_automated_tailoring SEMI_AUTOMATED_TAILORING
                        This option allows semi-automated tailoring with
                        customized quantile for gap score. e.g. -sa 0.75,0.6
                        [Modern samples with gap score outside uppper 0.75
                        quantile are out, ancient samples with gap score
                        outside upper 0.60 quantile are out.]
  -sl SHORT_LIST, --short_list SHORT_LIST
                        Input a list of taxa to keep for alignment tailoring.
                        Optionally, the limit for missing info in each column
                        can be given and it is delimited by comma, or column
                        tailoring will be performed automatically. e.g.
                        list.txt,0.1
  -m MANUAL_TAILORING, --manual_tailoring MANUAL_TAILORING
                        Missing info for modern, ancient samples and each
                        column is manually. e.g. 0.1,0.2,0.3 [Modern taxa with
                        missing info > 0.1 are removed, and ancient taxa with
                        missing info >0.2 are removed. Afterwards, columns
                        with missing info > 0.3 are removed.]

~~~

## Gene content analysis (based on eggNOG annotation)

~~~Bash
usage: kegg_analysis.py [-h] [-gap_roary [GENE_ABSENCE_PRESENCE_ROARY]]
                        [-p [PANGENOME]] [-eggnog [EGGNOG_ANNOTATION]]
                        [-cd CORE_DENSITY] [-hmap] [-o_fig OUTPUT_FIGURE]
                        [-o_eggnog_subset OUTPUT_EGGNOG_SUBSET]
                        [-o_unique_genes OUTPUT_UNIQUE_GENES] [-s SPECIES]

optional arguments:
  -h, --help            show this help message and exit
  -gap_roary [GENE_ABSENCE_PRESENCE_ROARY], --gene_absence_presence_roary [GENE_ABSENCE_PRESENCE_ROARY]
                        The table of gene absence and presence from Roary.
  -p [PANGENOME], --pangenome [PANGENOME]
                        The pangenome file from Roary.
  -eggnog [EGGNOG_ANNOTATION], --eggnog_annotation [EGGNOG_ANNOTATION]
                        The annotation file from EggNOG.
  -cd CORE_DENSITY, --core_density CORE_DENSITY
                        proportion of isolates a gene must be in to be core.
                        default [1]
  -hmap, --heatmap      plotting the heatmap
  -o_fig OUTPUT_FIGURE, --output_figure OUTPUT_FIGURE
                        Specify the output figure name.
  -o_eggnog_subset OUTPUT_EGGNOG_SUBSET, --output_eggnog_subset OUTPUT_EGGNOG_SUBSET
                        output eggnog file specific for unique genes
  -o_unique_genes OUTPUT_UNIQUE_GENES, --output_unique_genes OUTPUT_UNIQUE_GENES
                        output the unique genes in a file.
  -s SPECIES, --species SPECIES
                        Choose the species you want to check for unique genes.
                        [TS_1, TS_2, Moralis]
  -p, --heatmap      plotting the heatmap
  -o_fig OUTPUT_FIGURE, --output_figure OUTPUT_FIGURE
                        Specify the output figure name.
  -o_eggnog_subset OUTPUT_EGGNOG_SUBSET, --output_eggnog_subset OUTPUT_EGGNOG_SUBSET
                        output eggnog file specific for unique genes
  -o_unique_genes OUTPUT_UNIQUE_GENES, --output_unique_genes OUTPUT_UNIQUE_GENES
                        output the unique genes in a file.
  -s SPECIES, --species SPECIES
                        Choose the species you want to check for unique genes.
                        [TS_1, TS_2, Moralis]
~~~
