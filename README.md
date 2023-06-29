# singlecell_mousemelanoma
Code for scRNAseq of CD45+ FACS sorted cells from 5555_V600-BRAFmut mouse melanoma tumours treated with Ranolazine and anti-PD

Raw and processed data can be downloaded from the following link: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-12315

- 0.read_data.r for loading data into R.
- 1.process_samples.r for filtering low quality cells and ribosomal genes, identifying samples and excluding non CD45 cells.
- 2.annotate_clusters_Javi for annotating clusters.
- 3.plot4paper.r for generating plots for the article
- 4.cellExhaustion_90q.r for characterizing cell exhaustion gene signature
