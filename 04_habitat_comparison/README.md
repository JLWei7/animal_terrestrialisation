# Principal Coordinates Analysis

We compared GO and Pfam compositions associated with novel genes of terrestrial animal clades to capture the function variation, performing both Principal Coordinates Analysis (PCoA) and Principal component analysis(PCA). 

For PCoA (which is described in main text), please find the R scripts in **04_PCoA_and_PCA/scripts/pcoa_goterms_terrestrial_species_ellipses_statistic.R** and **pcoa_pfams_terrestrial_species_ellipses_statistic.R**, and the input files **novel_goterms_terrestrialspecies_df.csv**, **novel_pfams_terrestrial_species_df.csv**. 

For PCA (which is described in supplementary text), please find the R script in **04_principal_component_analysis/scripts/pca_novel_goterms.R**, and the input file **novel_goterms_terrestrialspecies_df.csv**. (This is the script example of GO terms of novel genes, we also generated PCA diagrams for GO terms of novel genes without Bdelloidea, GO terms of ancestral genes, GO terms of ancestral genes without tetrapods. Please see all DataFrames in **04_principal_component_analysis**)

In addition, based on the dataframe of semi-terrestrial and fully terrestrial groups generated above, we did enrichment analysis to infer which biological function enriched in each group.  Please find the R scripts in **04_PCoA_and_PCA/scripts/pcoa_goterms_terrestrial_species_enrich.R** and **pcoa_pfams_terrestrial_species_enrich.R**.



