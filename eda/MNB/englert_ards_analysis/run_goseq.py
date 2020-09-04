##############################################################################
#   Given a TSV file of expression data, run the GSVA algorithm
##############################################################################

from optparse import OptionParser
import json
import numpy as np
import pandas as pd
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.conversion import localconverter
import rpy2.robjects as ro
import sys
from collections import defaultdict

def run_GOseq(de_genes, all_genes, gene_to_gene_sets, gene_lens):
    #gene_set_names = []
    #gene_lists = []
    #for gene_set, genes in gene_set_to_genes.items():
    #    gene_set_names.append(gene_set)
    #    gene_lists.append(genes)
    
    de_genes = set(de_genes)
    de_status = [int(g in de_genes) for g in all_genes]

    with localconverter(ro.default_converter + pandas2ri.converter):
        #r_cts = ro.conversion.py2rpy(df_counts)
        r_de_status = ro.conversion.py2rpy(de_status)
        r_all_genes = ro.conversion.py2rpy(all_genes)
        r_gene_sets = ro.vectors.ListVector(gene_to_gene_sets)
        r_gene_lens = ro.conversion.py2rpy(gene_lens)
    rstring="""
        function(de_status, all_genes, gene_sets, gene_lens) {
            library(goseq)
            gene_lens <- unlist(gene_lens)
            de_status <- unlist(de_status)
            all_genes <- unlist(all_genes)
            names(de_status) <- all_genes
            gene_sets <- lapply(gene_sets, as.character)
            pwf <- nullp(de_status, "hg19", "geneSymbol", bias.data = gene_lens)
            res <- goseq(pwf, "hg19","geneSymbol", gene2cat=gene_sets)
            fdr <- p.adjust(res$over_represented_pvalue, method="BH")
            res$fdr <- fdr
            #res <- res[p.adjust(res$over_represented_pvalue, method="BH")<.05]
            df <- data.frame(res)
            df
        }
    """
    r_func = ro.r(rstring)
    r_res = r_func(r_de_status, r_all_genes, r_gene_sets, r_gene_lens)
    with localconverter(ro.default_converter + pandas2ri.converter):
        df_res = ro.conversion.rpy2py(r_res)
    df_res = df_res.set_index('category')
    return df_res

def _parse_gene_sets(gene_sets_f):
    gene_set_to_genes = {}
    with open(gene_sets_f, 'r') as f:
        for l in f:
            toks = l.split('\t')
            gene_set = toks[0]
            genes = [x.strip() for x in toks[2:]]
            gene_set_to_genes[gene_set] = genes
    return gene_set_to_genes


def main():
    usage = "" # TODO
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--transpose", action="store_true", help="Take transpose of input")
    parser.add_option("-o", "--out_file", help="Output file")
    (options, args) = parser.parse_args()

    with open('../raw_data/Albany_and_Englert/ebseq_ARDS.COVID.v.NO_HSCT/Down.Genes.pp95.txt', 'r') as f:
        de_genes = [l.strip() for l in f]

    df = pd.read_csv('../raw_data/AHNMJYDMXX_rsem/genes.tpm.no_hg.no_C054.tab', sep='\t', index_col=0)
    all_genes = list(df.index)


    len_df = pd.read_csv('effective_gene_length.tsv', sep='\t', index_col=0)
    len_df = len_df.loc[all_genes]
    print(len_df)

    assert set(de_genes) < set(all_genes)

    #gene_set_to_genes = _parse_gene_sets('../gene_sets/h.all.v7.1.symbols.gmt') 
    gene_set_to_genes = _parse_gene_sets('../gene_sets/c5.bp.v7.1.symbols.gmt')
    gene_to_gene_sets = defaultdict(lambda: [])
    for gene_set, genes in gene_set_to_genes.items():
        for gene in genes:
            gene_to_gene_sets[gene].append(gene_set)
    gene_to_gene_sets = dict(gene_to_gene_sets)
    
    run_GOseq(de_genes, all_genes, gene_to_gene_sets, list(len_df['effective_length']))

if __name__ == "__main__":
    main()
