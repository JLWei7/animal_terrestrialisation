
rm( list = ls( ) )

library(UpSetR)

setwd("pfams_novel_event/")


# Set the names of 11 sets
file_names <- list.files(pattern = "\\.txt$")
sets <- lapply(file_names, function(x) {
  readLines(x)
})
names(sets) <- gsub(".txt$", "", file_names)


# Convert list of sets to an incidence matrix
incidence_matrix <- fromList(sets)


# Filter sets that appear in more than one file
filtered_matrix <- incidence_matrix[rowSums(incidence_matrix) > 1, ]

diagram_order <- c("Node11_Tetrapoda_Homo_sapiens",
                   "Node10_Hexapoda_Drosophila_melanogaster",
                   "Node9_Armadillidium_Armadillidium_nasatum",
                   "Node8_Myriapoda_Rhysida_immarginata",
                   "Node7_Arachnida_Centruroides_sculpturatus",
                   "Node6_Onychophora_Epiperipatus_broadwayi",
                   "Node5_Tardigrada_Ramazzottius_varieornatusa",
                   "Node4_Nematoda_Pristionchus_pacificus",
                   "Node3_Stylommatophora_Candidula_unifasciata",
                   "Node2_Clitellata_Eisenia_andrei",
                   "Node1_Bdelloida_Rotaria_sordida")

# Create UpSet plot with the filtered matrix
upset(
  filtered_matrix,
  sets = diagram_order,
  keep.order = TRUE,
  order.by = "freq",
  nintersects = 200,
  ylim(0,100)
)



