# Phylogenetic tree inference

Here are the steps of the phylogenetic tree inference. Some of them are so big that cannot be upload here, but are available upon request to the authors.

Based on the tree position of species (guide tree) inferred by previous literature, we built phylogeny of metazoans by the following steps:

1) We identified conserved orthologous genes based on Homo sapiens (943) using the [BUSCO v5.4.7](https://busco.ezlab.org/). These identified conserved protein sequences were aligned by [MAFFT v7.505](https://mafft.cbrc.jp/alignment/software/)

```bash
for fasta_file in busco_genes/*.fasta; do

    # Extract the base name of the file
    base_name=$(basename "${fasta_file}")

    # Construct the output file name
    output_file="./aligned_${base_name}"

    # Apply the mafft command
    mafft --amino --globalpair --maxiterate 1000 "${fasta_file}" > "${output_file}"
done 
```

2. The sequences are trimmed by [trimAl v1.4.rev15](https://vicfero.github.io/trimal/)

```bash
for fasta_file in aligned_busco_genes/*.fasta; do

    # Extract the base name of the file
    base_name=$(basename "${fasta_file}")

    # Construct the output file name
    output_file="./trimmed_${base_name}"
    
    # Run trimal
    trimal -in "${fasta_file}" -out "${output_file}" -automated1
done
```

3. Next, concatenate the trimmed alignments into a single supermatrix using [FASconCAT-G v1.05.1](https://github.com/PatrickKueck/FASconCAT-G)

```
perl FASconCAT-G_v1.05.1.pl -s -p -l
```

4. Build phylogeny tree using IQ-TREE v2.2.2.6101, with C60+G+I model, constraining the guide tree

```
iqtree2 -s FcC_supermatrix.phy -p FcC_supermatrix_partition.txt -m C60+G+I -g SpeciesTree_phylogeny.nwk -B 1000 -pre C60_CT_metazoa
```

Please see the final treefile in Newick format in **02_phylogeny_inference/C60_CT_metazoa.treefile**
