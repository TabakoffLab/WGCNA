**I make this file for every version of WGCNA

**Updated reconstructed transcriptome from Laura:
/data2/saba/BNLx.SHR/RNA-Seq.Brain.total/reconstruction/reconstruct.total.brain.FINAL.26Aug14.gtf

Parameters Used:

Getting dataset prepped for WGCNA:
1. DABG is considered p-value < 0.0001
2. Need 5% of the arrays to be expressed above background
3. Use the gene-level RNAseq reconstruction
4. Find present probesets that hit an exon in the RNAseq reconstruction (ps to exon/geneID key available from LS)
5. Look for outlier arrays
6. Correlate (based on strain means) the probesets within each gene
7. Collapse down (use collapseRows function) based on correlation > 0.25
8. If collapsed down, summarize measurement using 1st principal component
9. Heritability filter: keep only TC that have herit > 0.25
10. Calculate the strain means
11. Look for outlier strains
CODE = step1_getTranscriptClusters.R & intersectPSandRNAseqTC.txt (BEDTOOLS code)


Performing WGCNA:
1. unsigned network
2. power = 9
3. minimum module size = 5
4. deepSplit = 4
5. correlation type = biweight midcorrelation (bicor)
CODE: step2_RNAseqTC_WGCNA_20141218.R


  
