#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )
library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

read_and_process <- function(file_path, file_tag) {
  # read the tsv file
  df <- read_tsv(file_path, show_col_types = FALSE)
  
  if (!"Number_of_Homology_Groups" %in% names(df)) {
    warning(paste("Column 'Number_of_Homology_Groups' does not exist in", file_path))
    df$Number_of_Homology_Groups <- NA
  }
  
  # process the dataframe
  df <- df %>%
    select(TermID, Number_of_Homology_Groups) %>%
    rename(!!sym(file_tag) := Number_of_Homology_Groups)
  
  return(df)
}
# set file paths
file_paths <- c(
  "node11_bottomGO_HGcount_novel_reordered_descrip.tsv",
  "node11_anc3_novel_Hsap_GO_BP_CC_MF_count_HG_only_bottom_10GO_reordered_descrip.tsv",
  "node11_anc2_novel_Hsap_GO_BP_CC_MF_count_HG_only_bottom_10GO_reordered_descrip.tsv",
  "node11_anc1_novel_Hsap_GO_BP_CC_MF_count_HG_only_bottom_10GO_reordered_descrip.tsv",
  "node10_bottomGO_HGcount_novel_reordered_descrip.tsv",
  "node10_anc3_novel_Dmel_GO_BP_CC_MF_count_HG_only_bottom_10GO_reordered_descrip.tsv",
  "node10_anc2_novel_Dmel_GO_BP_CC_MF_count_HG_only_bottom_10GO_reordered_descrip.tsv",
  "node8anc3_node10anc1_Pancrustacea_count_HG_only_bottom_10GO_reordered_descrip.tsv",
  "node9_bottomGO_HGcount_novel_reordered_descrip.tsv",
  "node9_anc2_novel_Anas_GO_BP_CC_MF_count_HG_only_bottom_10GO_reordered_descrip.tsv",
  "node8_bottomGO_HGcount_novel_reordered_descrip.tsv",
  "node8anc3_node10anc1_Pancrustacea_count_HG_only_bottom_10GO_reordered_descrip.tsv",
  "node7anc1_node8anc2_Arthropoda_count_HG_only_bottom_10GO_reordered_descrip.tsv",
  "node6anc3_node8anc1_Lobopodia_count_HG_only_bottom_10GO_reordered_descrip.tsv",
  "node7_bottomGO_HGcount_novel_reordered_descrip.tsv",
  "node7_anc3_novel_Cscu_GO_BP_CC_MF_count_HG_only_bottom_10GO_reordered_descrip.tsv",
  "node7anc1_node8anc2_Arthropoda_count_HG_only_bottom_10GO_reordered_descrip.tsv",
  "node6_bottomGO_HGcount_novel_reordered_descrip.tsv",
  "node6anc3_node8anc1_Lobopodia_count_HG_only_bottom_10GO_reordered_descrip.tsv",
  "node5anc3_node6anc2_Panarthropoda_count_HG_only_bottom_10GO_reordered_descrip.tsv",
  "node4anc3_node6anc1_Cyptovermes_count_HG_only_bottom_10GO_reordered_descrip.tsv",
  "node5_bottomGO_HGcount_novel_reordered_descrip.tsv",
  "node5anc3_node6anc2_Panarthropoda_count_HG_only_bottom_10GO_reordered_descrip.tsv",
  "node4_bottomGO_HGcount_novel_reordered_descrip.tsv",
  "node4anc3_node6anc1_Cyptovermes_count_HG_only_bottom_10GO_reordered_descrip.tsv",
  "node4_anc2_novel_Ppac_GO_BP_CC_MF_count_HG_only_bottom_10GO_reordered_descrip.tsv",
  "node1anc1_node4anc1_Protostomia_count_HG_only_bottom_10GO_reordered_descrip.tsv",
  "node3_bottomGO_HGcount_novel_reordered_descrip.tsv",
  "node3_anc3_novel_Cuni_GO_BP_CC_MF_count_HG_only_bottom_10GO_reordered_descrip.tsv",
  "node3_anc2_novel_Cuni_GO_BP_CC_MF_count_HG_only_bottom_10GO_reordered_descrip.tsv",
  "node3_anc1_novel_Cuni_GO_BP_CC_MF_count_HG_only_bottom_10GO_reordered_descrip.tsv",
  "node2_bottomGO_HGcount_novel_reordered_descrip.tsv",
  "node2_anc3_novel_Eand_GO_BP_CC_MF_count_HG_only_bottom_10GO_reordered_descrip.tsv",
  "node2_anc1_novel_Eand_GO_BP_CC_MF_count_HG_only_bottom_10GO_reordered_descrip.tsv",
  "node1_bottomGO_HGcount_novel_reordered_descrip.tsv",
  "node1_anc3_novel_Rsor_GO_BP_CC_MF_count_HG_only_bottom_10GO_reordered_descrip.tsv",
  "node1_anc2_novel_Rsor_GO_BP_CC_MF_count_HG_only_bottom_10GO_reordered_descrip.tsv",
  "node1anc1_node4anc1_Protostomia_count_HG_only_bottom_10GO_reordered_descrip.tsv"
)

# column naming
file_tags <- c("node11_event_Tetropoda", "node11_ancestor3_Dipnotetrapodomorpha", "node11_ancestor2_Sarcopterygii", "node11_ancestor1_Osteichthyans", "node10_event_Hexapoda", "node10_ancestor3_Allotriocarida", "node10_ancestor2_Altocrustacea", "node10_ancestor1_Pancrustacea", "node9_event_Armadillidium", "node9_ancestor2_Malacostraca", "node8_event_Myriapoda", "node8_ancestor3_Pancrustacea", "node8_ancestor2_Arthropoda", "node8_ancestor1_Lobopodia", "node7_event_Arachnida", "node7_ancestor3_Euchelicerata", "node7_ancestor1_Arthropoda", "node6_event_Onychophora", "node6_ancestor3_Lobopodia", "node6_ancestor2_Panarthropoda", "node6_ancestor1_Cryptovermes", "node5_event_Tardigrada", "node5_ancestor3_Panarthropoda", "node4_event_Nematoda", "node4_ancestor3_Cryptovermes", "node4_ancestor2_Ecdysozoa", "node4_ancestor1_Protostomia", "node3_event_Stylommatophora", "node3_ancestor3_Pulmonata", "node3_ancestor2_Panpulmonata", "node3_ancestor1_Heterobranchia", "node2_event_Clitellata", "node2_ancestor3_Capi-Clite", "node2_ancestor1_Pleistoannelida", "node1_event_Bdelloidea", "node1_ancestor3_Rotifera", "node1_ancestor2_Lophotrochozoa", "node1_ancestor1_Protostomia")

list_dfs <- map2(file_paths, file_tags, read_and_process)

# merge all dataframes by TermID 
merged_df <- reduce(list_dfs, full_join, by = "TermID")


print(merged_df)

# melt the data for ggplot
data_melted <- melt(merged_df, id.vars = 'TermID')

# replace Inf and -Inf with NA
data_melted$value[data_melted$value == Inf] <- NA
data_melted$value[data_melted$value == -Inf] <- NA

# calculate finite min and max values for the scale
valid_values <- data_melted$value[is.finite(data_melted$value)]
scale_min <- min(valid_values)
scale_max <- max(valid_values)


# replace NA values with 0
data_melted$value <- replace(data_melted$value, is.na(data_melted$value), 0)





# Creating the heatmap
##color: purple to green
coul <- colorRampPalette(brewer.pal(9, "PRGn"))(25)

##color: brown to dark green
coul <- colorRampPalette(brewer.pal(9, "BrBG"))(25)

##color: red to blue
coul <- colorRampPalette(brewer.pal(9, "RdBu"))(25)

###log color
ggplot(data_melted, aes(x = TermID, y = variable, fill = value)) + 
  geom_tile() + 
  geom_text(aes(label = as.integer(value)), size = 3, color = "black") +
  scale_fill_gradientn(colors = coul, trans = "log",
                       limits = c(scale_min,scale_max),
                       na.value = "darkgray") + 
  theme_minimal() + 
  theme(
    axis.text.x = element_text(size = 10, angle = 30, vjust = 1, hjust = 1, color = "black"), 
    axis.text.y = element_text(size = 10, color = "black"), 
    axis.title.x = element_text(vjust = 0), 
    axis.title.y = element_text(size = 13, vjust = 2, margin = margin(r = 20, unit = "pt")), 
    legend.text = element_text(size = 10),
#    plot.margin = margin(t = 30, r = 10, b = 10, l = 80, unit = "pt"), 
    plot.title = element_text(size = 15, face = "bold")
  ) +
  coord_fixed(ratio = 0.6) +
  xlab("GO Terms") + 
  ylab("Variables") +
  ggtitle("HG_numbers_heatmap_10bottom_GO") +
  labs(fill = "Number of homology groups")





