# scGRNom: a computational pipeline of integrative multi-omics analyses for predicting cell-type disease genes and regulatory networks


## Summary

Understanding gene regulatory mechanisms in varous human disesaes is still challenging, especially at the cell-type level. Recent single-cell multi-omics data enable studying cell-type gene regulatory networks (GRNs). For instance,  a bunch of algorithms have been developed to decipher cell-type GRNs from single-cell transcriptomic data [1]. As shown in the flowchart below, Here we introduce **scGRNom** - a new computational pipeline to integrate single cell multi-omics data including cell-type chromatin interactions, cell-type epigenomics and single cell transcriptomics to predict cell-type GRNs. The output networks link transcription factors, distal regulatory elements and target genes at the cell-type level. As a demo, we applied **scGRNom** to single cell data in the human brain (e.g., excitatory and inhibitory neurons, microglia, oligodendrocyte). In addition to predicting gene regulatory networks, the **scGRNom** pipeline is also able to identify cell-type disease genes and regulatory elements (e.g., enhancers, promoters) using the disesae associated SNPs (e.g., from GWAS) and cell-type GRN. Finally, **scGRNom** also works as a general purpose tool for predicting gene regulatory networks (e.g., at the bulk tissue level) from multi-omics and disease genes. 
<p align="left">
  <img width="1000" src="https://github.com/daifengwanglab/scGRNom/blob/master/pipeline.png">
</p>

## Hardward Requirements

The analysis is based on R 4.0. You will only need a standard computer with enough RAM to support the operations. For predicting gene regulatory networks, a *Linux* system with 32 GB RAM and 32GB storage would be enough to support.

## Software Requirements

Users should install the following packages prior to using the scGRNom from an R terminal:
```{r}
install.packages(c('glmnet', 'data.table', 'dplyr', 'parallel', 'doParallel', 'foreach', 'Seurat', 'Rmagic'))
```
Besides, data from packages *BSgenome.Hsapiens.UCSC.hg19*, *TxDb.Hsapiens.UCSC.hg19.knownGene* , *SNPlocs.Hsapiens.dbSNP142.GRCh37* and *JASPAR2018* are also used for our project. These packages can also be installed in Bioconductor.
```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg19', 'TxDb.Hsapiens.UCSC.hg19.knownGene','JASPAR2018','GenomicInteractions','MotifBreakR','S4Vectors',
                       'biomaRt','SNPlocs.Hsapiens.dbSNP142.GRCh37','GenomicRanges','TFBSTools','motifmatchr','GenomicFeatures','GenomeInfoDb','IRanges'))
```

## Functions in the scGRNom pipeline

There are three functions inside the package.

1. **scGRNom_interaction** :
    * main input:
        * hic_interaction: a data frame containing the variables chr1,start1,end1,chr2,start2,end2
          or chr,start1,end1,start2,end2; 
          chr1,ch2 or chr should have the following format 'chrD'(D represents digit 1-22 or X Y);                          
          start1,end1,start2,end2 should be integer type;
          
        * enhancer: a data frame containing chr,start,end for enhancers
        
    * output: a data frame containing gene, gene_chr, promoter_start, promoter_end,               
      enh_chr,enh_start,enh_end
      
The function, scGRNom_interaction inputs the chromatin interaction data (e.g., Hi-C) and predicts all possible interactions between enhancers and promoters in the data or the user-provided list; i.e., ones from Topologically Associating Domains (TADs) in Hi-C data. In addition, the function uses an R package, GenomicInteractions [2] to annotate interacting regions and link them to genes. The genome annotation was from TxDb.Hsapiens.UCSC.hg19.knownGene [3].
 
2. **scGRNom_getTF** :
    * main input: 
        * df: output data frame from **scGRNom_interaction**
    
    * output: a data.table which contains transcription factors for each region
    
Function scGRNom_getTF infers the transcription factor binding sites (TFBS) based on consensus binding site sequences in the enhancers and promoters that potentially interact from the previous step, scGRNom_interaction. It outputs a reference gene regulatory network linking these TF, enhancers and/or promoters of genes. In particular, this function uses TFBSTools [4] to obtain the position weight matrices of the TFBS motifs from the JASPAR database [5] and predict the TFBS locations on the enhancers and promoters via mapping TF motifs. It further links TFs with binding sites on all possible interacting enhancers and promoters, and outputs the reference regulatory network. Furthermore, this function can run on a parallel computing version via an R package, motifmatchr [6] for computational speed-up. 

3. **scGRNom_getNt** :
    * main input: 
        * df: output data.table from **scGRNom_getTF**
        * gexpr: gene expression data in which each row represents a gene and each column represents an observation
        * open_chrom : a list of user-defined chromatin accessible regions
        * extension_bps : an extension length (bp) for overlapping enhancers with chromatin accessible regions

        
    * output:  a data frame containing TG, TF, promoter, enhancer and coef
    
The function, scGRNom_getNt predicts the final gene regulatory network based on the TF-target gene expression relationships in the reference network. The reference gene regulatory network from the previous step provides all possible regulatory relationships (wires) between TF, enhancers, and target genes. However, changes in gene expression may trigger different regulatory wires. To refine our maps and determine the activity status of regulatory wires, this function applies elastic net regression, a machine learning method that has successfully modelled gene regulatory networks in our previous work [7].  In particular, given a gene expression dataset and a reference network from scGRNom_getTF, the function uses the TF expression to predict each target gene expression, and finds the TF with high regression coefficients, indicating an active regulatory influence on the target gene’s expression in the gene expression data. The final gene regulatory network consists of the TF with high elastic net coefficients, target genes and the linked enhancers from their reference network links if any. However, the chromatin interacting regions are broad, so that many TFs likely have binding sites on them.  If the chromatin accessibility information is available (e.g., from scATAC-seq data for a cell type), the function is also able to filter the enhancers based on their chromatin accessibilities and then output the network links only having the enhancers with high accessibility (e.g., overlapped with scATAC-seq peaks). The parameter “open_chrom” inputs a list of user-defined chromatin accessible regions.

4. **scGRNom_disGenes** :
    * main input: 
        * df: output data.table from **scGRNom_getNt**
        * gwas_snps: a list of GWAS SNPs associated with a disease
        * extension_bps : set an extension length (bp)
        
    * output: a data frame containing disease genes
    
For identifying cell-type disease genes and regulatory elements (e.g., enhancers, promoters). This function’s input includes a cell-type gene regulatory network and a list of GWAS SNPs associated with a disease. The function uses an R package, GenomicRanges [8], to overlap these disease SNPs with the enhancers and promoters of the input cell-type gene regulatory network, and then find the ones that interrupt the binding sites of  regulatory TFs (TFBSs) on the enhancers and promoters by motifbreakR[9]. It finally maps the overlapped enhancers or promoters and TFs with interrupted TFBSs onto the input network to find the linked genes and enhancers/promoters as the output cell-type disease genes and regulatory elements.
    
    
## Demo

This demo applies **scGRNom** to predict the gene regulatory network from single cell multi-omics for microglia, an important cell type in the brains. In particular, the microglia's chromatin interactome data and enhancers data are available in [10]. The single-cell gene expression data (UMI) including microglial cells is *DER-22_Single_cell_expression_raw_UMI* at http://resource.psychencode.org/ [11].

**Step 0: data preprocesssing**

```{r}
library(readxl)
interactome_data = read_xlsx("PLAC-seq promoter interactome map.xlsx",
                       sheet = 'Microglia interactome',skip = 2)[,1:6]

enhancers = read_xlsx("PLAC-seq promoter interactome map.xlsx",
                       sheet = 'Microglia enhancers',skip = 2,
                       col_names = c('chr','start','end'))
    
gene_expression = read.table('DER-22_Single_cell_expression_raw_UMI.tsv',
                       header = T,row.names = 1,sep = '\t')
                       
# Remove genes that are expressed in less than 100 cells (based on the gene expression data)
gene_expression = gene_expression[rowSums(gene_expression != 0) > 100,]

# Normalize and scale gene expression by Seurat 4.0 [14] for further removing noises and batches across cell types.
library(Seurat)
gexpr_dropout <- CreateSeuratObject(counts = gene_expression)
gexpr_dropout <- NormalizeData(gexpr_dropout)
gexpr_dropout <- ScaleData(gexpr_dropout, features = rownames(gexpr_dropout))

# Impute the single cell gene expression of all cells to address potential dropout issues using MAGIC method.
library(Rmagic)
gexpr_imputed <- magic(t(data.frame(gexpr_dropout@assays$RNA@data)))

# Select microglia cells
gexpr <- gexpr_imputed$result[grep('Micro',rownames(gexpr_imputed$result)),]

# For each cell type, remove the lowly expressed genes with log10(sum of imputed gene expression levels of the cells of the cell type+1) < 1.
gexpr <- t(gexpr[,log10(colSums(gexpr)+1)> 1])
```
The processed gene expression data for microglia is uploaed in the data folder.

**Step 1: Find chromatin interactions between enhancers and promoters**

We find the interacting enhancers and gene promoters in microglia.

```{r}
df1 <- scGRNom_interaction(interactome_data,enhancers)
# For convenience, we show the relevant outputs of a subset of df1.
head(df1[sample(nrow(df1),20),])
```
    ##        gene gene_chr promoter_start promoter_end enh_chr enh_start   enh_end
    ## 32919 EFCAB5    chr17       28293309     28298309   chr17  28134967  28135555
    ## 69344   SSH2    chr17       28085909     28090909   chr17  28007608  28011467
    ## 73487  SOCS6    chr18       67978927     67983927   chr18  67967584  67974384
    ## 4198   ABCA7    chr19        1056117      1061117   chr19   1179745   1181123
    ## 42653   SSH1    chr12      109248860    109253860   chr12 108981796 108988629
    ## 69901  RUNX1    chr21       36250511     36255511   chr21  36472485  36476546


**Step 2: Infer the transcription factor binding sites (TFBSs) on interacting enhancers and promoters**

We futher infer the TFBSs on the interacting enhancers and promoters from Step 1.

```{r}
df2 <- scGRNom_getTF(df1)
head(df2,1)
```
    ##    gene                 promoter                 enhancer                        promoter_TF
    ## 1: AKT3 chr1:244004085-244009085 chr1:244485607-244489015 IRF2,PPARG,RXRA,RREB1,VDR,TP53,...
    ##                            enhancer_TF
    ## 1: FOXF2,IRF2,NR3C1,NFIC,TLX1,STAT1,...
    
**Step 3: Predict TF-target genes with high expression relationships by Elastic net regression**  

We input the microglial gene expression data from Step 0 to predict the microglial gene regulatory network that links TFs to target genes (TGs) via enhancers and/or promoters. 

```{r}
df3 <- scGRNom_getNt(df = df2, gexpr = gexpr)
head(df3)
```
     ##    TG    TF                 enhancer                 promoter     TFbs      coef      mse 
     ##  AKT3 FOXP1 chr1:244485607-244489015 chr1:244004085-244009085     both 0.4408181 0.0010998
     ##  AKT3 FOXP1 chr1:243649950-243653651 chr1:244004085-244009085 promoter 0.4408181 0.0010998 
     ##  AKT3 FOXP1 chr1:243653917-243655561 chr1:244004085-244009085     both 0.4408181 0.0010998 
     ##  AKT3 FOXP1 chr1:244485607-244489015 chr1:244004387-244009387     both 0.4408181 0.0010998 
     ##  AKT3 FOXP1 chr1:243649950-243653651 chr1:244004387-244009387 promoter 0.4408181 0.0010998 
     ##  AKT3 FOXP1 chr1:243653917-243655561 chr1:244004387-244009387     both 0.4408181 0.0010998

**Step 3 (optional): Predict TF-target genes with high expression relationships by Elastic net regression**

We input the open chromatin regions for microglia by scATAC-seq in [12] to predict the microglial network that only consiste of edges linking high accessible enhancers in microglia. 

```{r}
library(readxl)
chromatin_access_regions <- read_xlsx("open_chromatin_regions.xlsx", sheet = 'Feature Binarization Peaks',skip = 16)

# Select chromation accessible regions of microglia
open_chrom_regions <- chromatin_access_regions[which(chromatin_access_regions$Microglia == 1),]
open_chrom_regions <- data.frame(na.omit(open_chrom_regions[,c('hg38_Chromosome', 'hg38_Start', 'hg38_Stop')]))
head(open_chrom_regions)
```
    ##   hg38_Chromosome hg38_Start hg38_Stop
    ## 1 chr1                816832    817459
    ## 2 chr1                817749    818154
    ## 3 chr1                818160    818369
    ## 4 chr1                906823    907024
    ## 5 chr1               1259974   1260175
    
```{r}
# get links having the enhancers with high accessiblity
df3 <- scGRNom_getNt(df = df2, gexpr = gexpr, open_chrom = open_chrom_regions, extension_bps = 2000)
head(df3)
```
    ##        TG    TF               enhancer               promoter   TFbs        coef      mse   
    ## TMEM170B BACH2 chr6:11212580-11213758 chr6:11535960-11540960 enhancer 0.08476907 0.0002148594 
    ## TMEM170B BACH2 chr6:11524968-11527308 chr6:11535960-11540960 enhancer 0.08476907 0.0002148594 
    ## TMEM170B BATF3 chr6:11205411-11208812 chr6:11535960-11540960 enhancer 0.03291128 0.0002148594 
    ## TMEM170B BATF3 chr6:11217315-11221289 chr6:11535960-11540960 enhancer 0.03291128 0.0002148594 
    ## TMEM170B CEBPB chr6:11212580-11213758 chr6:11535960-11540960 enhancer -0.09745105 0.0002148594 
    ## TMEM170B CLOCK chr6:11212580-11213758 chr6:11535960-11540960 enhancer -0.05542389 0.0002148594 
    
**Step 4: Identify cell-type disease genes and regulatory elements from cell-type GRN and disease associated SNPs**

Last, using microglial gene regulaotry network from Step 3 and the SNPs associated with Alzheimer’s Disease (AD) from the AD GWAS summary statistics [13] (p<5e-5, MAF>0.01), we identify the microglial AD genes. 

```{r}
ad_gwas <- data.frame(read.table("ad_gwas.txt"))[,c('CHR', 'BP', 'SNP')]
head(ad_gwas)
```
    ##                CHR       BP         SNP
    ## 1:5177909_D_I    1  5177909 rs368441845
    ## 1:5177910_C_A    1  5177910  rs60669416
    ## 1:6989416_G_T    1  6989416 rs111677930
    ## 1:7000462_T_A    1  7000462 rs112900196
    ## 1:19294695_I_D   1 19294695 rs199535293
    ## 1:19340155_D_I   1 19340155  rs11445162
    
```{r}
disease_genes <- scGRNom_disGenes(df = df3, gwas_snps = ad_gwas, extension_bps =1000)
head(disease_genes)
```
     ## disease_genes
     ##   ZNF689
     ##   CPEB2
     ##   IL34
     ##   ZNF688
     ##   ABAT
     ##   ERCC1
     ##   NYAP1
     ##   ATMIN
     ##   AP4E1
     
## Reference and source
1. Pratapa A, Jalihal AP, Law JN, Bharadwaj A, Murali TM. Benchmarking algorithms for gene regulatory network inference from single-cell transcriptomic data. Nat Methods. 2020;17:147–54. 
2. Harmston, N., Ing-Simmons, E., Perry, M., Baresic, A., Lenhard, B. GenomicInteractions: R package for handling genomic interaction data [Internet]. 2020. Available from: https://github.com/ComputationalRegulatoryGenomicsICL/GenomicInteractions/
3. Carlson, Marc. TxDb.Hsapiens.UCSC.hg19.knownGene: Annotation package for TxDb object(s) [Internet]. Bioconductor; 2015. Available from: https://bioconductor.org/packages/release/data/annotation/html/TxDb.Hsapiens.UCSC.hg19.knownGene.html
4. Tan G, Lenhard B. TFBSTools: an R/bioconductor package for transcription factor binding site analysis. Bioinforma Oxf Engl. 2016;32:1555–6. 
5. Fornes O, Castro-Mondragon JA, Khan A, van der Lee R, Zhang X, Richmond PA, et al. JASPAR 2020: update of the open-access database of transcription factor binding profiles. Nucleic Acids Res. 2020;48:D87–92. 
6. Schep, Alicia. motifmatchr: Fast Motif Matching in R [Internet]. 2019. Available from: https://www.bioconductor.org/packages/release/bioc/html/motifmatchr.html
7. Wang D, Liu S, Warrell J, Won H, Shi X, Navarro FCP, et al. Comprehensive functional genomic resource and integrative model for the human brain. Science. 2018;362. 
8. Lawrence M, Huber W, Pagès H, Aboyoun P, Carlson M, Gentleman R, et al. Software for computing and annotating genomic ranges. PLoS Comput Biol. 2013;9:e1003118.
9. Coetzee SG, Coetzee GA, Hazelett DJ. motifbreakR: an R/Bioconductor package for predicting variant effects at transcription factor binding sites. Bioinforma Oxf Engl. 2015;31:3847–9. 
10. Nott A, Holtman IR, Coufal NG, Schlachetzki JCM, Yu M, Hu R, et al. Brain cell type-specific enhancer-promoter interactome maps and disease-risk association. Science. 2019;366:1134–9. 
11. Lake BB, Chen S, Sos BC, Fan J, Kaeser GE, Yung YC, et al. Integrative single-cell analysis of transcriptional and epigenetic states in the human adult brain. Nat Biotechnol. 2018;36:70–80. 
12. Corces MR, Shcherbina A, Kundu S, Gloudemans MJ, Frésard L, Granja JM, et al. Single-cell epigenomic analyses implicate candidate causal variants at inherited risk loci for Alzheimer’s and Parkinson’s diseases. Nat Genet. 2020;52:1158–68. 
13. Panagiotou OA, Ioannidis JPA, for the Genome-Wide Significance Project. What should the genome-wide significance threshold be? Empirical replication of borderline genetic associations. Int J Epidemiol. 2012;41:273–86. 
14. Stuart T, Butler A, Hoffman P, Hafemeister C, Papalexi E, Mauck WM, et al. Comprehensive Integration of Single-Cell Data. Cell. 2019;177:1888-1902.e21. 
