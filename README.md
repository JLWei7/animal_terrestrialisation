# Animal Terrestrialisation
Scripts and datasets for the analysis on convergent genome evolution of animal terrestrialisation

## Genomes download
We compiled 154 genome samplings from published project uploaded in UniProt, NCBI, Ensembl and other resources. After download, we extracted canonical proteins from these whole genomes.

- Genomes from Uniprot: canonical proteins can be dowloaded directly from [Uniprot webpage](https://www.uniprot.org/).
- Genomes from NCBI: we run [Orthofinder](https://github.com/davidemms/OrthoFinder) in-house script [primary_transcript.py](https://davidemms.github.io/orthofinder_tutorials/running-an-example-orthofinder-analysis.html) to retrieve the canonical proteins. (Note there are several genomes cannot be identified).
- Unidentified genomes from NCBI and genomes from Ensembl and other resources: we run [cd-hit](https://github.com/weizhongli/cdhit) to retrieve the canonical proteins.

## Quality check
We use [BUSCO](https://busco.ezlab.org/) to check completeness of the genomes, with database *eukaryota_odb10* for these three unicells genomes, and *metazoa_odb10* for all other metazoan genomes.

## Homology groups (HG) inference
For homology groups inference, we run Orthofinder2

`orthofinder.py -f 154_canonical_genomes/ -t 28 -M msa -S diamond -A mafft -T fasttree`

Please find the HGs infered by our dataset in the folder **00_homology_groups**