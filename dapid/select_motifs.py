"""Contains code to select motifs with evidence
"""

import os
import pandas as pd

from dapid.motifs import PWM
from dapid.motifs import get_consensus_from_weights

def overlap_by_ensembl(hocomoco_metadata_file, homer_metadata_file, expression_file, out_file):
    """Merge all the tables
    """
    hocomoco_metadata = pd.read_table(hocomoco_metadata_file, sep='\t')
    homer_metadata = pd.read_table(homer_metadata_file, sep='\t', header=None)
    homer_metadata.columns = ['homer_consensus', 'Model', 'homer_score']
    rna_file = pd.read_table(expression_file, sep='\t')
    rna_file['ensembl_gene_id'] = rna_file.index
    
    # first merge hocomoco and homer
    all_metadata = hocomoco_metadata.merge(homer_metadata, on='Model')

    # then merge metadata with RNA
    full_df = all_metadata.merge(rna_file, on='ensembl_gene_id')

    # and save out
    full_df.to_csv(out_file, sep='\t')
    
    return None


def filter_by_id(metadata_file, gene_cluster_file, out_file):
    """Goes from cluster file to metadata to filter
    """
    metadata = pd.read_table(metadata_file, sep='\t')
    gene_clusters = pd.read_table(gene_cluster_file, sep='\t', header=None)
    gene_clusters.columns = ['trajectory_cluster', 'ensembl_gene_id']

    # merge
    full_df = metadata.merge(gene_clusters, on='ensembl_gene_id')

    # write
    full_df.to_csv(out_file, sep='\t')
    
    return None


def add_motif_consensus_sequence(metadata_file, pwm_dict, out_file, model_col=2):
    """Add consensus motif
    """
    with open(out_file, 'w') as out:
        with open(metadata_file, 'r') as fp:
            for line in fp:
                if "ensembl" in line:
                    # write out header
                    out.write('{}\tconsensus\n'.format(line.strip()))
                    continue

                pwm_name = line.strip().split('\t')[model_col]
                pwm_sequence = get_consensus_from_weights(pwm_dict[pwm_name].weights)
                out.write('{0}\t{1}\n'.format(line.strip(), pwm_sequence))

    return None


def run(args):

    # RNA - first set up a histogram to determine cutoff point of expressed
    # and then output smaller matrix of expressed genes
    prefix = 'dapid'
    args.rna['rna_norm_filt_mat'] = '{0}/{1}.expressed.mat.gz'.format(args.folders['results_dir'], prefix)
    if not os.path.isfile(args.rna['rna_norm_filt_mat']):
        filter_rna = ("filter_rna.R {0} {1} {2}").format(args.rna['rna_norm_mat'],
                                                         "{0}/{1}".format(args.folders['results_dir'], prefix),
                                                         4)
        os.system('mkdir -p {}'.format(args.folders['results_dir']))
        os.system(filter_rna)
    
    # for motifs
    # First grab metadata and match ENSEMBL ID to it
    os.system('mkdir -p {}'.format(args.folders['annot_dir']))
    args.annotations['hocomoco_homer_file'] = "{0}/{1}".format(
        args.folders['annot_dir'],
        os.path.basename(args.annotations['hocomoco_homer_file_url']))
    if not os.path.isfile(args.annotations['hocomoco_homer_file']):
        get_file = "wget {0} -O {1}".format(args.annotations['hocomoco_homer_file_url'],
                                            args.annotations['hocomoco_homer_file'])
        print get_file
        os.system(get_file)

    # get HOMER simplified data
    args.annotations['hocomoco_homer_simple'] = "{0}/hocomoco_homer_names.txt".format(args.folders['annot_dir'])
    if not os.path.isfile(args.annotations['hocomoco_homer_simple']):
        get_homer_names = ("cat {0} | "
                           "grep '>' | "
                           "awk -F '>' '{{ print $2 }}' > "
                           "{1}").format(args.annotations['hocomoco_homer_file'],
                                         args.annotations['hocomoco_homer_simple'])
        print get_homer_names
        os.system(get_homer_names)

    # Now match hocomoco metadata to expression matrix and only keep complete overlaps
    args.annotations['motifs_w_expression_file'] = "{0}/{1}.motifs_w_expression.txt".format(
        args.folders['results_dir'], prefix)
    if not os.path.isfile(args.annotations['motifs_w_expression_file']):
        overlap_by_ensembl(args.annotations['hocomoco_metadata_file'],
                           args.annotations['hocomoco_homer_simple'],
                           args.rna['rna_norm_filt_mat'],
                           args.annotations['motifs_w_expression_file'])


    # Now do the same where we filter for early and late TF sets
    # go ahead and use trajectory groups
    early_gene_clusters = [1, 6, 10]
    late_gene_clusters = [2, 4, 5]

    # early
    args.rna['early_gene_set'] = '{0}/{1}.progenitor.gene_ids.txt'.format(args.folders['results_dir'],
                                                                          prefix)
    if not os.path.isfile(args.rna['early_gene_set']):
        early_grep = "\"^[{}]\t\"".format('|'.join([str(cluster_num) for cluster_num in early_gene_clusters]))
        get_gene_ids = ("cat {0} | "
                        "grep -P {1} | "
                        "awk -F '.' '{{ print $1 }}' > "
                        "{2}").format(
                            args.rna['rna_trajectory_clusters'],
                            early_grep,
                            args.rna['early_gene_set'])
        print get_gene_ids
        os.system(get_gene_ids)

    # and filter by geneid
    args.annotations['progenitor_motifs_w_expr'] = '{0}/{1}.motifs_w_expression.progenitor.txt'.format(
        args.folders['results_dir'], prefix)
    
    if not os.path.isfile(args.annotations['progenitor_motifs_w_expr']):
        filter_by_id(args.annotations['motifs_w_expression_file'],
                     args.rna['early_gene_set'],
                     args.annotations['progenitor_motifs_w_expr'])

    # late
    args.rna['late_gene_set'] = '{0}/{1}.differentiated.gene_ids.txt'.format(args.folders['results_dir'],
                                                                             prefix)
    if not os.path.isfile(args.rna['late_gene_set']):
        late_grep = "\"^[{}]\t\"".format('|'.join([str(cluster_num) for cluster_num in late_gene_clusters]))
        get_gene_ids = ("cat {0} | "
                        "grep -P {1} | "
                        "awk -F '.' '{{ print $1 }}' > "
                        "{2}").format(
                            args.rna['rna_trajectory_clusters'],
                            late_grep,
                            args.rna['late_gene_set'])
        print get_gene_ids
        os.system(get_gene_ids)
    
    # and filter by geneid
    args.annotations['late_motifs_w_expr'] = '{0}/{1}.motifs_w_expression.differentiated.txt'.format(
        args.folders['results_dir'], prefix)
    
    if not os.path.isfile(args.annotations['late_motifs_w_expr']):
        filter_by_id(args.annotations['motifs_w_expression_file'],
                     args.rna['late_gene_set'],
                     args.annotations['late_motifs_w_expr'])
    

    # extract consensus from PWM file
    pwms = PWM.get_homer_pwms(args.annotations['hocomoco_homer_file'])
    
    # At end, have 3 lists: ALL, EARLY, LATE
    args.annotations['all_final'] = '{0}/{1}.motifs.consensus.expression.all.txt'.format(args.folders['results_dir'], prefix)
    args.annotations['early_final'] = '{0}/{1}.motifs.consensus.expression.progenitor.txt'.format(args.folders['results_dir'], prefix)
    args.annotations['late_final'] = '{0}/{1}.motifs.consensus.expression.differentiated.txt'.format(args.folders['results_dir'], prefix)

    if not os.path.isfile(args.annotations['all_final']):
        add_motif_consensus_sequence(args.annotations['motifs_w_expression_file'], pwms, args.annotations['all_final'], model_col=2)

    if not os.path.isfile(args.annotations['early_final']):
        add_motif_consensus_sequence(args.annotations['progenitor_motifs_w_expr'], pwms, args.annotations['early_final'], model_col=3)

    if not os.path.isfile(args.annotations['late_final']):
        add_motif_consensus_sequence(args.annotations['late_motifs_w_expr'], pwms, args.annotations['late_final'], model_col=3)

    return None
