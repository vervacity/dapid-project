"""Contains code to select motifs with evidence
"""

import os












def run(args):


    # RNA - first set up a histogram to determine cutoff point of expressed
    # and then output smaller matrix of expressed genes
    prefix = 'dapid'
    filter_rna = ("filter_rna.R {0} {1} {2}").format(args.rna['rna_norm_mat'],
                                                     "{0}/{1}".format(args.folders['results_dir'], prefix),
                                                     4)
    os.system(filter_rna)
    

    quit()
    
    # for motifs
    # First grab metadata and match ENSEMBL ID to it


    # Then use those ENSEMBL IDs to filter the expression matrix



    # At end, have 3 lists: ALL, EARLY, LATE
    # with associated gene ID ensembl, hgnc, motif sequence, motif name


    

    
    

    return None
