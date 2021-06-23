#!/usr/bin/env python

import pandas as pd
import argparse
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def read_args(args):

    parser = argparse.ArgumentParser()


    parser.add_argument('-gap_roary',
                        '--gene_absence_presence_roary',
                        nargs = '?',
                        help = "The table of gene absence and presence from Roary.",
                        type = str)

    parser.add_argument('-p',
                        '--pangenome',
                        nargs = '?',
                        help = 'The pangenome file from Roary.',
                        type = str)

    parser.add_argument('-eggnog',
                        '--eggnog_annotation',
                        nargs = '?',
                        help = 'The annotation file from EggNOG.',
                        type = str)

    parser.add_argument('-cd',
                        '--core_density',
                        help = 'proportion of isolates a gene must be in to be core. default [1]',
                        type = float,
                        default = 1.0)
    
    parser.add_argument('-hmap',
                      '--heatmap',
                      help = 'plotting the heatmap',
                      action = 'store_true')

    parser.add_argument('-o_fig',
                      '--output_figure',
                      help = 'Specify the output figure name.',
                      type = str)

    parser.add_argument('-o_eggnog_subset',
                        '--output_eggnog_subset',
                        help = 'output eggnog file specific for unique genes',
                        type = str)

    parser.add_argument('-o_unique_genes',
                        '--output_unique_genes',
                        help = 'output the unique genes in a file.',
                        type = str)

    parser.add_argument('-s',
                        '--species',
                        help = 'Choose the species you want to check for unique genes. [TS_1, TS_2, Moralis]',
                        type = str)
    
    return vars(parser.parse_args())

def gap_file_parser(gap_file):

    gap_df = pd.read_csv(gap_file, header = 0)

    return gap_df

species_genomes = {'calc_2082.bin.4': 'TS1', 'calc_2086.bin.1': 'TS1',
                   'calc_2084.bin.6': 'TS2', 'calc_2090.bin.4': 'TS2',
                   'calc_2091.bin.4': 'TS2', 'calc_2093.bin.1': 'TS2', 
                   'calc_2094.bin.7': 'TS2', 'calc_2095.bin.8': 'TS2',
                   'calc_2096.bin.3': 'TS2', 'calc_2099.bin.5': 'TS2',
                   'GCA_001639275_SGB720': 'Moralis', 'GCF_000529525.1_SGB720': 'Moralis',
                   'GCF_900289035.1_PRJEB24872_SGB720': 'Moralis', 'GCF_902384065.1_SGB720': 'Moralis', 
                   'calc_2102.bin.3': 'Moralis'}

def pangenome_paser(pangenome_file):

    seqs = SeqIO.parse(pangenome_file, 'fasta')
    seqs_2 = SeqIO.parse(pangenome_file, 'fasta')
    label_link = {}
    
    for seq in seqs:
        labels = seq.description.split(' ')
        label_link[labels[0]] = labels[1]
    for seq in seqs_2:
        labels = seq.description.split(' ')
        label_link[labels[1]] = labels[0]
        
    return label_link

def eggnog_parser(eggnog_annotation):

    eggnog_annotation = pd.read_csv(eggnog_annotation, skiprows = 4, skipfooter = 3, sep = '\t', engine = 'python')

    return eggnog_annotation


def select_ko_pathway_genes(eggnog_df):

    selected_df = eggnog_df.loc[(eggnog_df['KEGG_ko'] != '-') & (eggnog_df['KEGG_Pathway'] != '-')]

    return selected_df    

def extracting_gap_file(gap_file, cut_off = 1.0):
    # n_exclude: excluding genes which are present in N genomes
    gap_file['No. isolates'] = gap_file['No. isolates'].astype(int)
    gap_file['No. sequences'] = gap_file['No. sequences'].astype(int)
    
    n_exclude = cut_off * 15
    gap_file = gap_file.loc[(gap_file['No. isolates'] < n_exclude) & (gap_file['No. sequences'] < n_exclude)]

    subset_df = gap_file[['Gene', 'calc_2082.bin.4', 'calc_2086.bin.1',
                          'calc_2084.bin.6', 'calc_2090.bin.4',
                          'calc_2091.bin.4', 'calc_2093.bin.1',
                          'calc_2094.bin.7', 'calc_2095.bin.8',
                          'calc_2096.bin.3', 'calc_2099.bin.5',
                          'GCA_001639275_SGB720', 'GCF_000529525.1_SGB720',
                          'GCF_900289035.1_PRJEB24872_SGB720', 'GCF_902384065.1_SGB720',
                          'calc_2102.bin.3']]
    return subset_df

def heatmap_maker(df_, opt_fig):

    sns.set(font_scale=1)
    fig, ax = plt.subplots()

    df_heatmap = df_

    fig.set_size_inches(8, 12)

    ax = sns.heatmap(df_heatmap, xticklabels=True, yticklabels=False)
    ax.set(xlabel='', ylabel='KEGG')

    fig.savefig(opt_fig, dpi=300, bbox_inches = 'tight')


def unique_gene_identifier(cleaned_gap_df, genomes_list):
    # cleaned_gap_file: The gene absence and presence file 
    # genomes_list: a list of genomes to search for uniques to the rest
    # eggnog_file: the annotation file from eggnog. 
    
    Moralis_genomes = ['GCA_001639275_SGB720', 'GCF_000529525.1_SGB720',
                      'GCF_900289035.1_PRJEB24872_SGB720', 'GCF_902384065.1_SGB720',
                      'calc_2102.bin.3']
    
    remaining_genomes = [i for i in list(cleaned_gap_df) if i not in genomes_list + ['Gene']]
    unique_genes = []
    print('comparing to M. oralis')
    print(genomes_list)
    print(Moralis_genomes)
    for index, row in cleaned_gap_df.iterrows():
        if (row[genomes_list].isna().sum().sum() == 0) and (row[remaining_genomes].isna().sum().sum() == len(remaining_genomes)):
            unique_genes.append(row['Gene'])

    return unique_genes



if __name__ == '__main__':

    pars = read_args(sys.argv)

    gap_df = gap_file_parser(pars['gene_absence_presence_roary']) 

    gene_group_label = pangenome_paser(pars['pangenome'])
    
    eggnog_annotations = eggnog_parser(pars['eggnog_annotation'])

    genes_ko_pathway = list(select_ko_pathway_genes(eggnog_annotations)['#query']) # Genes identified with KEGG_ko and KEGG_pathway
    
    gap_genes_ko_pathway = [gene_group_label[i] for i in genes_ko_pathway] # Converting gene label to group label
    
    subset_gap_file = extracting_gap_file(gap_df, pars['core_density'])
    abs_pre_in_ko_pathway_df = subset_gap_file[subset_gap_file.Gene.isin(gap_genes_ko_pathway)]
    abs_pre_in_ko_pathway_df_4selecting = subset_gap_file[subset_gap_file.Gene.isin(gap_genes_ko_pathway)]
    

    abs_pre_in_ko_pathway_df = abs_pre_in_ko_pathway_df.drop('Gene', axis = 1)
    binary_gap_file = abs_pre_in_ko_pathway_df.notnull().astype('int')


    if pars['heatmap']:
        heatmap_maker(binary_gap_file, pars['output_figure'])
    else:
    	print("Heatmap generating skipped.")
    
    TS_1_genomes = ['calc_2082.bin.4', 'calc_2086.bin.1']
    TS_2_genomes = ['calc_2084.bin.6', 'calc_2090.bin.4',
                    'calc_2091.bin.4', 'calc_2093.bin.1',
                    'calc_2094.bin.7', 'calc_2095.bin.8',
                    'calc_2096.bin.3', 'calc_2099.bin.5']
    
    Moralis_genomes = ['GCA_001639275_SGB720', 'GCF_000529525.1_SGB720',
                          'GCF_900289035.1_PRJEB24872_SGB720', 'GCF_902384065.1_SGB720',
                          'calc_2102.bin.3']
    
    
    species_genomes_link = {'TS_1': TS_1_genomes, 'TS_2': TS_2_genomes, 'Moralis': Moralis_genomes}
    if pars['output_eggnog_subset'] and pars['species']:
        print("We are searching for unique genes (identified with KEGG KO and KEGG Pathway) to {}".format(pars['species']))
        u_gene = unique_gene_identifier(abs_pre_in_ko_pathway_df_4selecting, species_genomes_link[pars['species']])
        u_gene_group = [gene_group_label[i] for i in u_gene]
        eggnog_subset_tab = eggnog_annotations[eggnog_annotations['#query'].isin(u_gene_group)] 
        eggnog_subset_tab.to_csv(pars['output_eggnog_subset'], sep = '\t', index = False)
    
    else:
        print('Skip unique genes (identified with KEGG KO and KEGG Pathway) searching......')


    if pars['output_unique_genes'] and pars['species']:
        print("We are searching for unique genes to {}".format(pars['species']))

        general_u_gene = unique_gene_identifier(subset_gap_file, species_genomes_link[pars['species']])
        general_u_gene_label = [gene_group_label[i] for i in general_u_gene]
        opt_unique_genes = open(pars['output_unique_genes'], 'w')
        for g in general_u_gene_label:
            opt_unique_genes.write(g + '\n')
        opt_unique_genes.close()

    else:
        print('Skip general unique genes.......')







