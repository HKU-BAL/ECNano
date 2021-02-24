# ECNano
A target captured medical exome ONT sequencing with amplicon variant calling workflow

## Introduction

Having a high sequencing error of ONT and limited throughput from a single MinION flowcell ultimately limits its applicability for accurate variant detection. Medical exome sequencing (MES) targets clinically significant exon regions, allowing rapid and comprehensive screening of pathogenic variants. By applying MES with MinION sequencing, the technology can achieve a more uniform capture of the target regions, shorter turnaround time, and lower sequencing cost per sample. 

ECNano is an out-of-the-box workflow comprising (1) a wet-lab protocol for ONT target enrichment sequencing and (2) a downstream variant detection workflow, Clair-ensemble. The ECNano wet-lab protocol was optimized to perform long-read target enrichment and ONT library preparation to stably generate high-quality MES data with adequate coverage. The subsequent variant-calling workflow, Clair-ensemble, adopted a fast RNN-based variant caller, Clair, and was optimized for target enrichment data. ECNano was evaluated on both reference DNA samples and patient samples.