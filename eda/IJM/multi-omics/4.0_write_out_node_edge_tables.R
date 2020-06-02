
data_path = "/Volumes/projects/All_20200428_COVID_plasma_multiomics/Correlation/node_edge_list_kendall_holm.RData"
load(data_path)

#out_path = "~/Desktop/UW_2020/COVID19_multi-omics/COVID-19_Multi-Omics/data/node_edge_list_kendall_holm.txt"
#lapply(ne_kendall, write, out_path, append=TRUE)

# node df
node_df_outpath = "~/Desktop/UW_2020/COVID19_multi-omics/COVID-19_Multi-Omics/data/kendall_nodes.txt"
write.table(ne_kendall$node, node_df_outpath, sep="\t", row.names = F)

# edge df
edge_df_outpath = "~/Desktop/UW_2020/COVID19_multi-omics/COVID-19_Multi-Omics/data/kendall_edges.txt"
write.table(ne_kendall$edge, edge_df_outpath, sep="\t", row.names = F)
