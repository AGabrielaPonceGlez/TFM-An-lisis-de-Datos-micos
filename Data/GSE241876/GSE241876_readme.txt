RNASeq

Raw reads in the FASTQ format were aligned to the mm10 reference genome using STAR aligner v2.7.6a [1]. The featureCounts read summarization program [2] was used to count mapped reads for each gene in the annotation file (Human GENECODE: Release 42 (GRCh38.p13)). Count files were then merged by an in-house R script and differential gene expression analysis was performed using DESeq2 [3]. DESeq2 normalizes the count data to account for differences in sequencing depth and composition. This is done using the median of ratios method, which involves dividing each gene's raw count by a size factor that accounts for library size and composition. DESeq2 models the variation in gene expression using a negative binomial distribution. The dispersion estimates are used to model the relationship between mean expression and variance, which is critical for accurate detection of differentially expressed genes. 
Raw (raw_count_matrix.csv) and DESeq normalized (DeseqNormalizedCount.csv) matrices are made available along with the FASTQ files used in current study.


1.	A. Dobin, T.R. Gingeras, Mapping RNA-seq Reads with STAR., Curr. Protoc. Bioinformatics. 51 (2015) 11.14.1-11.14.19. https://doi.org/10.1002/0471250953.bi1114s51 

2.	Y. Liao, G.K. Smyth, W. Shi, featureCounts: an efficient general-purpose program for assigning sequence reads to genomic features., Bioinformatics. 30 (2014) 923?930. https://doi.org/10.1093/bioinformatics/btt656 

3.	M.I. Love, W. Huber, S. Anders, Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2, Genome Biol. 15 (2014) 550?550. https://doi.org/10.1186/s13059-014-0550-8

