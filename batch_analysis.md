### Load packages

First of all we will R packages which we are going to use in this lab:

``` r
library("Seurat")
suppressMessages(require(gplots))
suppressMessages(require(ggplot2))
suppressMessages(require(matrixStats))
suppressMessages(require(mixOmics))
suppressMessages(require(RColorBrewer))
suppressMessages(require(Rtsne))
```

### Load Expression Values and Metadata

Let start with loading the matrix of raw gene expression counts and filtering away genes with median count across all cells below 1, it is very conservative but speeds up computations for this lab. Those genes are lowly expressed genes which should be excluded from the downstream analysis as they might lead to spurious results:

``` r
D <- read.table("data/ILC/ensembl_countvalues_ILC.csv",sep=",",header=T,row.names=1)
library("matrixStats")
D<-D[rowMedians(as.matrix(D))>=1,]
D[1:5,1:5]
```

    ##                 T74_P1_A9_ILC1 T74_P1_B4_NK T74_P1_B7_ILC2 T74_P1_B9_NK
    ## ENSG00000002587             25           44             41           40
    ## ENSG00000003056              0            0              0            1
    ## ENSG00000003402              2           11              2          657
    ## ENSG00000003756              0            1             48            3
    ## ENSG00000004534              1            1           4102            2
    ##                 T74_P1_D10_ILC2
    ## ENSG00000002587               5
    ## ENSG00000003056               1
    ## ENSG00000003402             924
    ## ENSG00000003756               1
    ## ENSG00000004534             672

``` r
dim(D)
```

    ## [1] 1465  648

For the sake of speed and simplicity of this lab we will select only 50% of most varying genes using coefficient of variation as a criterion:

``` r
D_var<-apply(D,1,function(x) sd(x)/mean(x))
D<-D[D_var>quantile(D_var,0.5),]
D[1:5,1:5]
```

    ##                 T74_P1_A9_ILC1 T74_P1_B4_NK T74_P1_B7_ILC2 T74_P1_B9_NK
    ## ENSG00000003056              0            0              0            1
    ## ENSG00000003402              2           11              2          657
    ## ENSG00000003756              0            1             48            3
    ## ENSG00000004534              1            1           4102            2
    ## ENSG00000004897              0            1              2            2
    ##                 T74_P1_D10_ILC2
    ## ENSG00000003056               1
    ## ENSG00000003402             924
    ## ENSG00000003756               1
    ## ENSG00000004534             672
    ## ENSG00000004897               0

``` r
dim(D)
```

    ## [1] 732 648

The rows of the matrix represent Ensembl gene IDs (you can convert them to gene symbols using biomaRt package) from 732 genes, and the columns are IDs from 648 cells from different individuals sequenced at different plates. To see how many individuals and plates we have let us load the meta-information and have a look:

``` r
M <- read.table("data/ILC/Metadata_ILC.csv",sep=",",header=T,row.names=1)
M$Plate<-matrix(unlist(strsplit(as.character(M$Plate),"_")),byrow=TRUE,ncol=2)[,2]
head(M)
```

    ##                 Plate Donor Celltype
    ## T74_P1_A9_ILC1     P1   T74     ILC1
    ## T74_P1_B4_NK       P1   T74       NK
    ## T74_P1_B7_ILC2     P1   T74     ILC2
    ## T74_P1_B9_NK       P1   T74       NK
    ## T74_P1_D10_ILC2    P1   T74     ILC2
    ## T74_P1_E1_ILC3     P1   T74     ILC3

``` r
summary(M)
```

    ##     Plate           Donor     Celltype  
    ##  Length:648         T74: 75   ILC1:127  
    ##  Class :character   T75:273   ILC2:139  
    ##  Mode  :character   T86:300   ILC3:308  
    ##                               NK  : 74

Thus we have cells from 3 individuals with IDs:

``` r
levels(factor(M$Donor))
```

    ## [1] "T74" "T75" "T86"

that were pooled together and sequenced at 4 plates with IDs:

``` r
levels(factor(M$Plate))
```

    ## [1] "P1" "P2" "P3" "P4"

and finally we have 4 cell-types with the following IDs:

``` r
levels(factor(M$Celltype))
```

    ## [1] "ILC1" "ILC2" "ILC3" "NK"

The variables Plate and Donor might potentially reflect technical variation due to batch.

### Checking for Genome-Wide Batch-Effects

Now let us check potential batch-effects in the data set. As we saw previously the cells were pooled from 3 and sequenced on 4 plates. Thoese are potential batches. We need to check how they affect gene expression genome-wide. One way to see it is to plot PCA and tSNE and color cells by batch:

``` r
library("mixOmics")
pca.ILC<-pca(log10(t(D+1)),ncomp=10,center=TRUE,scale=FALSE)
plot(pca.ILC)
```

![](batch_analysis_files/figure-markdown_github/PCA%20and%20tSNE-1.png)

``` r
plotIndiv(pca.ILC,group=factor(M$Plate),ind.names=FALSE,ellipse=FALSE,legend=TRUE,title="PCA PLOT, PLATE EFFECT",cex=1)
```

![](batch_analysis_files/figure-markdown_github/PCA%20and%20tSNE-2.png)

``` r
plotIndiv(pca.ILC,group=factor(M$Donor),ind.names=FALSE,ellipse=FALSE,legend=TRUE,title="PCA PLOT, DONOR EFFECT",cex=1)
```

![](batch_analysis_files/figure-markdown_github/PCA%20and%20tSNE-3.png)

``` r
library("Rtsne")
library("RColorBrewer")
set.seed(1)
tsne.out_expr<-Rtsne(t(log10(D+1)),initial_dims=20,verbose=TRUE,perplexity=30)
```

    ## Read the 648 x 20 data matrix successfully!
    ## Using no_dims = 2, perplexity = 30.000000, and theta = 0.500000
    ## Computing input similarities...
    ## Normalizing input...
    ## Building tree...
    ##  - point 0 of 648
    ## Done in 0.07 seconds (sparsity = 0.186466)!
    ## Learning embedding...
    ## Iteration 50: error is 64.110210 (50 iterations in 0.86 seconds)
    ## Iteration 100: error is 61.119195 (50 iterations in 0.34 seconds)
    ## Iteration 150: error is 60.871258 (50 iterations in 0.17 seconds)
    ## Iteration 200: error is 60.864903 (50 iterations in 0.16 seconds)
    ## Iteration 250: error is 60.867336 (50 iterations in 0.16 seconds)
    ## Iteration 300: error is 1.046757 (50 iterations in 0.17 seconds)
    ## Iteration 350: error is 0.978871 (50 iterations in 0.16 seconds)
    ## Iteration 400: error is 0.962507 (50 iterations in 0.15 seconds)
    ## Iteration 450: error is 0.954225 (50 iterations in 0.15 seconds)
    ## Iteration 500: error is 0.949245 (50 iterations in 0.15 seconds)
    ## Iteration 550: error is 0.944535 (50 iterations in 0.14 seconds)
    ## Iteration 600: error is 0.939772 (50 iterations in 0.14 seconds)
    ## Iteration 650: error is 0.936925 (50 iterations in 0.14 seconds)
    ## Iteration 700: error is 0.933510 (50 iterations in 0.14 seconds)
    ## Iteration 750: error is 0.930883 (50 iterations in 0.14 seconds)
    ## Iteration 800: error is 0.930321 (50 iterations in 0.14 seconds)
    ## Iteration 850: error is 0.929130 (50 iterations in 0.14 seconds)
    ## Iteration 900: error is 0.928918 (50 iterations in 0.14 seconds)
    ## Iteration 950: error is 0.929116 (50 iterations in 0.14 seconds)
    ## Iteration 1000: error is 0.928702 (50 iterations in 0.14 seconds)
    ## Fitting performed in 3.87 seconds.

``` r
palette(brewer.pal(length(levels(factor(M$Plate))),'Dark2'))
plot(tsne.out_expr$Y,main="tSNE PLOT, PLATE EFFECT",col=factor(M$Plate),xlab="tSNE1",ylab="tSNE2")
legend("topleft",levels(factor(M$Plate)),cex=1,fill=brewer.pal(length(levels(factor(M$Plate))),'Dark2'),inset=0.02)
```

![](batch_analysis_files/figure-markdown_github/PCA%20and%20tSNE-4.png)

``` r
palette(brewer.pal(length(levels(factor(M$Donor))),'Dark2'))
plot(tsne.out_expr$Y,main="tSNE PLOT, DONOR EFFECT",col=factor(M$Donor),xlab="tSNE1",ylab="tSNE2")
legend("topleft",levels(factor(M$Donor)),cex=1,fill=brewer.pal(length(levels(factor(M$Donor))),'Dark2'),inset=0.02)
```

![](batch_analysis_files/figure-markdown_github/PCA%20and%20tSNE-5.png)

We can immediately see that there is a slight plate related and a more pronounced donor related batch-effect. To further quantify it let us display how much of variation in each principal component is explained by the batch variables:

``` r
M$Plate<-factor(M$Plate)
M$Donor<-factor(M$Donor)
M$Celltype<-factor(M$Celltype)

pc_adj_r_squared<-matrix(NA,ncol=dim(pca.ILC$x)[2],nrow=dim(M)[2])
for(i in 1:dim(pca.ILC$x)[2])
{
  print(i)
  for(j in 1:dim(M)[2])
  {
    pc_adj_r_squared[j,i]<-summary(lm(pca.ILC$x[,i]~M[,j]))$adj.r.squared
  }
}
```

    ## [1] 1
    ## [1] 2
    ## [1] 3
    ## [1] 4
    ## [1] 5
    ## [1] 6
    ## [1] 7
    ## [1] 8
    ## [1] 9
    ## [1] 10

``` r
pc_adj_r_squared<-as.data.frame(pc_adj_r_squared)
colnames(pc_adj_r_squared)<-colnames(pca.ILC$x)
rownames(pc_adj_r_squared)<-colnames(M)
pc_adj_r_squared
```

    ##                 PC1        PC2        PC3        PC4        PC5
    ## Plate    0.05180021 0.30492744 0.01855468 0.01573704 0.12226500
    ## Donor    0.44214488 0.02107709 0.03756549 0.05522430 0.08987487
    ## Celltype 0.01117032 0.80817822 0.69380232 0.12293283 0.52031296
    ##                   PC6        PC7          PC8         PC9         PC10
    ## Plate     0.101180017 0.05824515  0.004760321 0.029266531 0.0333732685
    ## Donor     0.266018738 0.05135150  0.039968802 0.114280245 0.0081112259
    ## Celltype -0.002524276 0.02196863 -0.004043204 0.002775554 0.0001119801

``` r
library("gplots")
my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
heatmap.2(data.matrix(pc_adj_r_squared),cellnote=round(pc_adj_r_squared,3),notecol="black",density.info="none",trace="none",col = my_palette, margins=c(8,10),dendrogram="row",Colv="NA",scale="row",main="ILC scRNAseq",cexRow=1,cexCol=1)
mtext("Adjusted R^2 of Association between PCs and Phenotypes")
```

![](batch_analysis_files/figure-markdown_github/heatmap%20batch%20effects-1.png)

From the heatmap above it is clear that 44% of PC1 is explained by Donor batch (Celltype contribution is negligible) while 31% of PC2 is explained by Plate batch. In addition, PC2 variation is mainly due to both Celltype and Plate batch, so those two variables are coupled and removing Plate batch could affect the variation due to Celltype which is our phenotype of interest.

### Checking for Per-Gene Batch-Effects

Now let us check batch-effects in the individual genes and figure out genes that are most influenced by batch. Let us check the effect of e.g. Donor on the expression of individual genes. For this purpose we will add a batch variable to the meta information:

``` r
M$batch<-M$Donor
head(M)
```

    ##                 Plate Donor Celltype batch
    ## T74_P1_A9_ILC1     P1   T74     ILC1   T74
    ## T74_P1_B4_NK       P1   T74       NK   T74
    ## T74_P1_B7_ILC2     P1   T74     ILC2   T74
    ## T74_P1_B9_NK       P1   T74       NK   T74
    ## T74_P1_D10_ILC2    P1   T74     ILC2   T74
    ## T74_P1_E1_ILC3     P1   T74     ILC3   T74

Now we will rank all genes by the percentage of variation in their expression explained by the batch factor variable:

``` r
adj_r_squared<-vector()
for(i in 1:dim(D)[1])
{
  adj_r_squared<-append(adj_r_squared,summary(lm(as.numeric(D[i,])~M$batch))$adj.r.squared)
}
adj_r_squared[adj_r_squared<0]<-0
var_expl<-data.frame(genes=rownames(D),var_expl=adj_r_squared)
var_expl<-var_expl[order(-var_expl$var_expl),]
head(var_expl,20)
```

    ##               genes   var_expl
    ## 686 ENSG00000237973 0.57125722
    ## 662 ENSG00000225630 0.44484858
    ## 711 ENSG00000256618 0.34349035
    ## 665 ENSG00000227097 0.23221658
    ## 663 ENSG00000225840 0.18835454
    ## 271 ENSG00000128016 0.15331071
    ## 723 ENSG00000265150 0.10330701
    ## 678 ENSG00000233927 0.09043639
    ## 725 ENSG00000265735 0.07175026
    ## 706 ENSG00000253676 0.07022771
    ## 315 ENSG00000135018 0.06775005
    ## 606 ENSG00000197183 0.06372930
    ## 17  ENSG00000026025 0.06247762
    ## 692 ENSG00000243101 0.05875688
    ## 691 ENSG00000241529 0.05507034
    ## 198 ENSG00000115085 0.05442117
    ## 363 ENSG00000142347 0.05411619
    ## 434 ENSG00000159388 0.05308169
    ## 641 ENSG00000210184 0.05274778
    ## 602 ENSG00000197043 0.05198533

``` r
barplot(var_expl$var_expl[1:7],names=var_expl$genes[1:7],ylab="Variance Explained",main="Top Genes Influenced by Batch",col="darkred",las=1,cex.names=0.7,ylim=c(0,0.6))
```

![](batch_analysis_files/figure-markdown_github/genes%20affected%20by%20batch-1.png)

Thus we conclude that the batch-effect might explain up to ~60% of variation in gene expression of most affected genes.

Let us also check which batch is the most influential:

``` r
my_batches<-levels(M$batch)
my_genes<-as.character(var_expl$genes)
adj_r_squared_per_species<-list()
for(j in 1:length(my_genes))
{
  adj_r_squared_per_batch<-vector()
  for(i in 1:length(my_batches))
  {
    this_batch<-factor(ifelse(as.character(M$batch)==my_batches[i],my_batches[i],paste0("NOT_",my_batches[i])))
    adj_r_squared_per_batch<-append(adj_r_squared_per_batch,summary(lm(as.numeric(D[my_genes[j],])~this_batch))$adj.r.squared)
    adj_r_squared_per_batch[adj_r_squared_per_batch<0]<-0
  }
  adj_r_squared_per_species[[j]]<-adj_r_squared_per_batch
}
batch_matrix<-matrix(unlist(adj_r_squared_per_species),ncol=length(my_batches),byrow=TRUE)
batch_df<-as.data.frame(batch_matrix)
rownames(batch_df)<-my_genes
colnames(batch_df)<-my_batches
batch_df[1:10,1:3]
```

    ##                          T74        T75          T86
    ## ENSG00000237973 0.5716699265 0.06034322 5.576458e-02
    ## ENSG00000225630 0.0411298712 0.44570792 2.790453e-01
    ## ENSG00000256618 0.3444955019 0.03258765 3.626696e-02
    ## ENSG00000227097 0.0110460685 0.16974774 2.306282e-01
    ## ENSG00000225840 0.0009035628 0.15678101 1.797534e-01
    ## ENSG00000128016 0.0195376825 0.09133304 1.544442e-01
    ## ENSG00000265150 0.0995958719 0.02603911 1.785723e-05
    ## ENSG00000233927 0.0002767789 0.08905885 7.171379e-02
    ## ENSG00000265735 0.0552986786 0.03877870 5.575910e-04
    ## ENSG00000253676 0.0000000000 0.05945957 6.767231e-02

``` r
n <- length(my_batches)
library("RColorBrewer")
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
barplot(t(batch_matrix[1:7,]),beside=TRUE,ylab="Variance Explained",names.arg=my_genes[1:7],legend=my_batches,col=col_vector[1:length(my_batches)],cex.names=0.7,xlab="Genes",main="Batch Contribution to Genes")
```

![](batch_analysis_files/figure-markdown_github/most%20influential%20batch-1.png)

Overall, it seems that different genes are affected by contributions from different individuals.

However, what would be variance explained by chance alone? To elucidate this we will perform a number of shuffling of expression vectors for each gene individually and calculate the shuffled variance explained, i.e. the variance explained by chance for each gene. Further, we will plot the noise zone as three standard deviations beyond the mean of shuffled variance explained.

``` r
N_shuffle_ind<-100
ranked_genes<-as.character(var_expl$genes)
shuffle_stat_ind<-list()
for(i in 1:length(ranked_genes))
{
  adj_r_squared_shuffle_ind<-vector()
  for(j in 1:N_shuffle_ind)
  {
    gene_shuffle<-
      as.numeric(D[ranked_genes[i],][sample(1:dim(D)[2],dim(D)[2])])
    adj_r_squared_shuffle_ind<-
      append(adj_r_squared_shuffle_ind,summary(lm(gene_shuffle~M$batch))$adj.r.squared)
  }
  adj_r_squared_shuffle_ind[adj_r_squared_shuffle_ind<0]<-0
  shuffle_stat_ind[[i]]<-adj_r_squared_shuffle_ind
}
shuffle_matrix_ind<-t(matrix(unlist(shuffle_stat_ind),byrow=TRUE,ncol=N_shuffle_ind))
shuffle_matrix_ind[1:5,1:5]
```

    ##            [,1]        [,2]         [,3]        [,4]        [,5]
    ## [1,] 0.00000000 0.000000000 0.0051290157 0.000000000 0.002402295
    ## [2,] 0.00000000 0.000000000 0.0000000000 0.004522597 0.007753193
    ## [3,] 0.00000000 0.000000000 0.0001351077 0.000000000 0.003867364
    ## [4,] 0.00141051 0.003177335 0.0003162582 0.003711627 0.000000000
    ## [5,] 0.00000000 0.000000000 0.0000000000 0.005447060 0.000000000

``` r
library("matrixStats")
noise<-colMeans(shuffle_matrix_ind)+3*colSds(shuffle_matrix_ind)

library("ggplot2")
library("RColorBrewer")
Observed<-data.frame(ranked_genes=ranked_genes,var_expl=var_expl$var_expl)[1:6,]
ByChance<-data.frame(ranked_genes=ranked_genes,var_expl=noise)[1:6,]
ggplot(NULL, aes(ranked_genes,var_expl)) + 
  geom_bar(aes(fill="Observed"),data=Observed,stat='identity') +
  geom_bar(aes(fill="ByChance"),data=ByChance,stat='identity') +
  ggtitle("Observed vs. Resampled Variance Explained by Batch") +
  xlab("Rank of Genes") + ylab("Variance Explained") +
  scale_x_discrete(limits = ranked_genes[1:6], expand = c(0.01,0.01)) +
  scale_y_continuous(expand = c(0.0,0)) + coord_cartesian(ylim=c(0,0.6)) +
  theme(plot.title=element_text(hjust = 0.5)) + 
  scale_fill_manual(name="Legend",values=c(brewer.pal(8,'Dark2')[2],brewer.pal(8,'Dark2')[1]))
```

![](batch_analysis_files/figure-markdown_github/noise%20zone%20calculation-1.png)

Above we displayed just a few most influenced by batch genes (observed) together with the shuffled variance explained (by chance). Here we can see that e.g. ENSG00000237973 seems to be strongly influenced by batch effects since the observed variance explained is beyond three standard deviations from the mean of variance explained by chance.

Now let us display the variance explained by batch for all genes in the ordered way (from largest to lowest) by a curve (observed) together with shuffled variance explained (by chance). Again, mean + 3 standard deviations from the mean is the noise zone boundary:

``` r
plot(var_expl$var_expl~seq(1:length(var_expl$var_expl)),xlab="Rank of Genes",ylab="Variance Explained",col="blue",main="Observed vs. Resampled Variance Explained by Batch",type="l",lwd=2,ylim=c(0,0.2))
lines(noise~seq(1:length(var_expl$var_expl)),col="red")
legend("topright",inset=0.02,c("Observed","ByChance"),col=c("blue","red"),lty=c(1,1))
```

![](batch_analysis_files/figure-markdown_github/unnamed-chunk-7-1.png)

Here we can see that the observed variance explained hits the noie zone for approximately gene \#200 meaning that approximately top 200 genes ordered by their variance explained by batch are significantly influenced by batch, the rest genes are safe to use in the downstream analysis. We can also do a formal statistical test and calculate a p-value of significance of deviation from the noise zone. The p-value represents how many times shuffled variance explained by batch is equal or below the noise zone. We also apply Benjamini-Hochberg correction of th p-values for multiple testing:

``` r
p_res<-vector()
for(i in 1:dim(shuffle_matrix_ind)[2])
{
  p_res<-append(p_res,sum(shuffle_matrix_ind[,i]>=var_expl$var_expl[i])/dim(shuffle_matrix_ind)[1])
}
p_res_BH<-p.adjust(p_res,method="BH")
plot(p_res_BH~seq(1:dim(shuffle_matrix_ind)[2]),type='l',col="darkgreen",xlab="Rank of Genes",ylab="p-value BH",main="Significance of Deviation from Noise")
abline(h=0.05,col="red")
```

![](batch_analysis_files/figure-markdown_github/unnamed-chunk-8-1.png)

Again, we see that the top ca. 200 genes seem to be significanly influenced by the batch-effects. Finally, let us display the genes that are significantly influenced by batch-effects to have a look and memorize them:

``` r
problematic_genes<-data.frame(species=ranked_genes,var_expl_by_batch=var_expl$var_expl,pvalue=p_res,FDR=p_res_BH)
problematic_genes<-problematic_genes[order(problematic_genes$FDR,problematic_genes$pvalue,-problematic_genes$var_expl_by_batch),]
bad_genes<-problematic_genes[problematic_genes$FDR<=0.05,]
good_genes<-problematic_genes[problematic_genes$FDR>0.05,]
```

Thus here are "bad genes" ordered by how strongly they are affected by batch-effects, i.e. the higher in the list the more affected:

``` r
head(bad_genes,50)
```

    ##            species var_expl_by_batch pvalue FDR
    ## 1  ENSG00000237973        0.57125722      0   0
    ## 2  ENSG00000225630        0.44484858      0   0
    ## 3  ENSG00000256618        0.34349035      0   0
    ## 4  ENSG00000227097        0.23221658      0   0
    ## 5  ENSG00000225840        0.18835454      0   0
    ## 6  ENSG00000128016        0.15331071      0   0
    ## 7  ENSG00000265150        0.10330701      0   0
    ## 8  ENSG00000233927        0.09043639      0   0
    ## 9  ENSG00000265735        0.07175026      0   0
    ## 10 ENSG00000253676        0.07022771      0   0
    ## 11 ENSG00000135018        0.06775005      0   0
    ## 12 ENSG00000197183        0.06372930      0   0
    ## 13 ENSG00000026025        0.06247762      0   0
    ## 14 ENSG00000243101        0.05875688      0   0
    ## 15 ENSG00000241529        0.05507034      0   0
    ## 16 ENSG00000115085        0.05442117      0   0
    ## 17 ENSG00000142347        0.05411619      0   0
    ## 18 ENSG00000159388        0.05308169      0   0
    ## 19 ENSG00000210184        0.05274778      0   0
    ## 20 ENSG00000197043        0.05198533      0   0
    ## 21 ENSG00000004534        0.05183540      0   0
    ## 22 ENSG00000226958        0.05041399      0   0
    ## 23 ENSG00000137076        0.04637637      0   0
    ## 24 ENSG00000079332        0.04278522      0   0
    ## 25 ENSG00000167895        0.04254538      0   0
    ## 26 ENSG00000141524        0.04224829      0   0
    ## 27 ENSG00000162231        0.04215924      0   0
    ## 28 ENSG00000243829        0.04189671      0   0
    ## 29 ENSG00000168887        0.04150684      0   0
    ## 30 ENSG00000096060        0.04113025      0   0
    ## 31 ENSG00000196683        0.04056862      0   0
    ## 32 ENSG00000211829        0.04055738      0   0
    ## 33 ENSG00000117984        0.03936293      0   0
    ## 34 ENSG00000183878        0.03906877      0   0
    ## 35 ENSG00000100385        0.03861384      0   0
    ## 36 ENSG00000223001        0.03724164      0   0
    ## 37 ENSG00000107742        0.03718143      0   0
    ## 38 ENSG00000170348        0.03658703      0   0
    ## 39 ENSG00000090104        0.03647631      0   0
    ## 40 ENSG00000139641        0.03626200      0   0
    ## 41 ENSG00000042493        0.03531457      0   0
    ## 42 ENSG00000171223        0.03421048      0   0
    ## 43 ENSG00000235065        0.03372497      0   0
    ## 44 ENSG00000177606        0.03346526      0   0
    ## 45 ENSG00000228847        0.03284347      0   0
    ## 46 ENSG00000204256        0.03279591      0   0
    ## 47 ENSG00000134545        0.03226591      0   0
    ## 48 ENSG00000076928        0.03225244      0   0
    ## 49 ENSG00000222414        0.03194271      0   0
    ## 50 ENSG00000068028        0.03130674      0   0

``` r
dim(bad_genes)[1]
```

    ## [1] 204

And here come genes that are ok to use for the downstream analysis since they are not significantly affected by batch effects:

``` r
head(good_genes,50)
```

    ##             species var_expl_by_batch pvalue        FDR
    ## 100 ENSG00000259865       0.018080726   0.02 0.06337662
    ## 118 ENSG00000135506       0.016423659   0.02 0.06337662
    ## 157 ENSG00000130303       0.012668994   0.02 0.06337662
    ## 159 ENSG00000126353       0.012503202   0.02 0.06337662
    ## 161 ENSG00000227765       0.012175431   0.02 0.06337662
    ## 168 ENSG00000158869       0.011413924   0.02 0.06337662
    ## 184 ENSG00000204482       0.010101335   0.02 0.06337662
    ## 194 ENSG00000111716       0.009557891   0.02 0.06337662
    ## 199 ENSG00000171862       0.009310417   0.02 0.06337662
    ## 202 ENSG00000110031       0.009150156   0.02 0.06337662
    ## 204 ENSG00000100353       0.009034326   0.02 0.06337662
    ## 209 ENSG00000188994       0.008662661   0.02 0.06337662
    ## 212 ENSG00000197622       0.008493330   0.02 0.06337662
    ## 224 ENSG00000148175       0.007909481   0.02 0.06337662
    ## 233 ENSG00000170906       0.007631268   0.02 0.06337662
    ## 234 ENSG00000113387       0.007546467   0.02 0.06337662
    ## 243 ENSG00000185043       0.007268468   0.02 0.06337662
    ## 247 ENSG00000157404       0.007095425   0.02 0.06337662
    ## 248 ENSG00000099783       0.007094136   0.02 0.06337662
    ## 249 ENSG00000127328       0.007025443   0.02 0.06337662
    ## 250 ENSG00000156587       0.006867084   0.02 0.06337662
    ## 251 ENSG00000070770       0.006848395   0.02 0.06337662
    ## 256 ENSG00000105939       0.006466572   0.02 0.06337662
    ## 287 ENSG00000106682       0.005597842   0.02 0.06337662
    ## 294 ENSG00000068308       0.005263436   0.02 0.06337662
    ## 299 ENSG00000123416       0.005208908   0.02 0.06337662
    ## 303 ENSG00000179144       0.004960911   0.02 0.06337662
    ## 116 ENSG00000125835       0.016517078   0.03 0.08749004
    ## 180 ENSG00000198931       0.010388816   0.03 0.08749004
    ## 206 ENSG00000163682       0.008900627   0.03 0.08749004
    ## 208 ENSG00000091317       0.008747865   0.03 0.08749004
    ## 215 ENSG00000103479       0.008295441   0.03 0.08749004
    ## 220 ENSG00000116560       0.008143729   0.03 0.08749004
    ## 221 ENSG00000130770       0.008138738   0.03 0.08749004
    ## 225 ENSG00000167996       0.007883188   0.03 0.08749004
    ## 228 ENSG00000166128       0.007767251   0.03 0.08749004
    ## 235 ENSG00000213684       0.007542630   0.03 0.08749004
    ## 236 ENSG00000152061       0.007471682   0.03 0.08749004
    ## 239 ENSG00000163902       0.007372931   0.03 0.08749004
    ## 260 ENSG00000198830       0.006408042   0.03 0.08749004
    ## 261 ENSG00000120738       0.006375578   0.03 0.08749004
    ## 271 ENSG00000055609       0.006134222   0.03 0.08749004
    ## 277 ENSG00000100592       0.005820364   0.03 0.08749004
    ## 278 ENSG00000165280       0.005813736   0.03 0.08749004
    ## 280 ENSG00000111679       0.005763356   0.03 0.08749004
    ## 283 ENSG00000168653       0.005686616   0.03 0.08749004
    ## 284 ENSG00000166136       0.005653444   0.03 0.08749004
    ## 203 ENSG00000119707       0.009094716   0.04 0.10804428
    ## 211 ENSG00000145287       0.008502212   0.04 0.10804428
    ## 216 ENSG00000153922       0.008295003   0.04 0.10804428

``` r
dim(good_genes)[1]
```

    ## [1] 528

### Correcting for Batch-Effects with Combat

Here we will use Combat in order to remove batch-effects due to different Donors. Combat has an advantage compared to simple linear regression batch-effects removal implemented e.g. in Limma since it uses Bayesian framework and works well even for low numbers of samples/cells. Particular advantage of the Combat batch effect removal is that it uses so-called “Bayesian shrinkage towards the mean”, i.e. correction of each gene uses a “mean” information from other genes which Limma does Frequentist batch-effect correction gene-by-gene independently. We assume that the different Donor and Plate can be viewed as a bacth which can be removed with Combat.

An important step to do before running Combat is to show the phenotype of interst (which is Celltype in our case) to the algorithm. In this way we are saying: do not touch the variation due to the Celltype, correct for the differences in Donor only.

``` r
library("sva")
modcombat = model.matrix(~as.factor(as.character(M$Celltype)), data=M)
combat_edata = ComBat(dat=D, batch=M$batch, mod=modcombat, par.prior=TRUE, prior.plots=TRUE)
```

    ## Found3batches

    ## Adjusting for3covariate(s) or covariate level(s)

    ## Standardizing Data across genes

    ## Fitting L/S model and finding priors

    ## Finding parametric adjustments

    ## Adjusting the Data

![](batch_analysis_files/figure-markdown_github/Combat-1.png)

``` r
D[1:5,1:5]
```

    ##                 T74_P1_A9_ILC1 T74_P1_B4_NK T74_P1_B7_ILC2 T74_P1_B9_NK
    ## ENSG00000003056              0            0              0            1
    ## ENSG00000003402              2           11              2          657
    ## ENSG00000003756              0            1             48            3
    ## ENSG00000004534              1            1           4102            2
    ## ENSG00000004897              0            1              2            2
    ##                 T74_P1_D10_ILC2
    ## ENSG00000003056               1
    ## ENSG00000003402             924
    ## ENSG00000003756               1
    ## ENSG00000004534             672
    ## ENSG00000004897               0

``` r
combat_edata[1:5,1:5]
```

    ##                 T74_P1_A9_ILC1 T74_P1_B4_NK T74_P1_B7_ILC2 T74_P1_B9_NK
    ## ENSG00000003056     -13.666835    -8.867236      -11.32170    -7.770285
    ## ENSG00000003402      36.225949    49.976078       31.28344   653.085592
    ## ENSG00000003756      16.882753    14.389821       73.11964    16.748636
    ## ENSG00000004534      -3.000701     9.716396     2301.19951    10.270803
    ## ENSG00000004897       5.345433     3.075161       12.11363     4.611484
    ##                 T74_P1_D10_ILC2
    ## ENSG00000003056      -10.224751
    ## ENSG00000003402      892.068221
    ## ENSG00000003756       17.687487
    ## ENSG00000004534      399.583837
    ## ENSG00000004897        9.040983

Let us now display a tSNE plot for the Combat corrected data set:

``` r
library("Rtsne")
set.seed(1)
tsne.out_combat<-Rtsne(t(combat_edata),initial_dims=10,verbose=TRUE,perplexity=30)
```

    ## Read the 648 x 10 data matrix successfully!
    ## Using no_dims = 2, perplexity = 30.000000, and theta = 0.500000
    ## Computing input similarities...
    ## Normalizing input...
    ## Building tree...
    ##  - point 0 of 648
    ## Done in 0.06 seconds (sparsity = 0.215435)!
    ## Learning embedding...
    ## Iteration 50: error is 63.721036 (50 iterations in 0.22 seconds)
    ## Iteration 100: error is 61.824577 (50 iterations in 0.23 seconds)
    ## Iteration 150: error is 61.761268 (50 iterations in 0.22 seconds)
    ## Iteration 200: error is 61.772186 (50 iterations in 0.22 seconds)
    ## Iteration 250: error is 61.775711 (50 iterations in 0.23 seconds)
    ## Iteration 300: error is 0.851962 (50 iterations in 0.17 seconds)
    ## Iteration 350: error is 0.776001 (50 iterations in 0.17 seconds)
    ## Iteration 400: error is 0.759151 (50 iterations in 0.18 seconds)
    ## Iteration 450: error is 0.751650 (50 iterations in 0.18 seconds)
    ## Iteration 500: error is 0.747354 (50 iterations in 0.18 seconds)
    ## Iteration 550: error is 0.743511 (50 iterations in 0.18 seconds)
    ## Iteration 600: error is 0.741052 (50 iterations in 0.18 seconds)
    ## Iteration 650: error is 0.739281 (50 iterations in 0.18 seconds)
    ## Iteration 700: error is 0.738501 (50 iterations in 0.18 seconds)
    ## Iteration 750: error is 0.737509 (50 iterations in 0.18 seconds)
    ## Iteration 800: error is 0.736602 (50 iterations in 0.19 seconds)
    ## Iteration 850: error is 0.735737 (50 iterations in 0.18 seconds)
    ## Iteration 900: error is 0.735655 (50 iterations in 0.18 seconds)
    ## Iteration 950: error is 0.735604 (50 iterations in 0.19 seconds)
    ## Iteration 1000: error is 0.735335 (50 iterations in 0.18 seconds)
    ## Fitting performed in 3.83 seconds.

``` r
#BEFORE CORRECTION
plot(tsne.out_expr$Y,main="ILC, tSNE PLOT, BEFORE COMBAT",col=factor(M$Donor),xlab="tSNE1",ylab="tSNE2")
legend("topleft",levels(factor(M$Donor)),cex=1,fill=brewer.pal(length(levels(factor(M$Donor))),'Dark2'),inset=0.02)
```

![](batch_analysis_files/figure-markdown_github/tSNE%20Combat-1.png)

``` r
#AFTER COMBAT CORRECTION
plot(tsne.out_combat$Y,main="ILC, tSNE PLOT, AFTER COMBAT",col=as.factor(M$Donor),xlab="tSNE1",ylab="tSNE2")
legend("topright",levels(as.factor(M$Donor)),cex=1,fill=unique(as.factor(M$Donor)),inset=0.02)
```

![](batch_analysis_files/figure-markdown_github/tSNE%20Combat-2.png)

We see that cells from different Donors seem to be more mixed after Combat correction.

### Correcting for Batch-Effects with MNN

Now we will perform batch-effects correction with Mutial Nearest Neighbours (MNN) algorithm from the "scran" R package developed in the lab of John Marioni, University of Cambridge. The idea of the method is to identify cells between batches having shortest Euclidean distance and hence (presumably!) belonging to the same sub-population (nearest neighbours). Then the difference in expression between the cells is assumed due to technical/batch variation and all other cells are corrected with respect to the identified nearest neighbour cells.

For simplicity we assume that T74 and T75 belong to the first batch while T86 is the second batch.

``` r
library("scran")
D_T74<-subset(D,select=colnames(D)[grepl("T74",colnames(D))])
D_T75<-subset(D,select=colnames(D)[grepl("T75",colnames(D))])
D_T86<-subset(D,select=colnames(D)[grepl("T86",colnames(D))])
mnn_corrected <- mnnCorrect(as.matrix(D_T74), as.matrix(D_T75),as.matrix(D_T86))
mnn_corrected$corrected[[1]][1:5,1:5]
```

    ##                 T74_P1_A9_ILC1 T74_P1_B4_NK T74_P1_B7_ILC2 T74_P1_B9_NK
    ## ENSG00000003056   0.0000000000 0.000000e+00   0.0000000000 4.494933e-05
    ## ENSG00000003402   0.0004876917 3.732144e-04   0.0001286902 2.953171e-02
    ## ENSG00000003756   0.0000000000 3.392858e-05   0.0030885647 1.348480e-04
    ## ENSG00000004534   0.0002438459 3.392858e-05   0.2639435885 8.989867e-05
    ## ENSG00000004897   0.0000000000 3.392858e-05   0.0001286902 8.989867e-05
    ##                 T74_P1_D10_ILC2
    ## ENSG00000003056    0.0001459003
    ## ENSG00000003402    0.1348119107
    ## ENSG00000003756    0.0001459003
    ## ENSG00000004534    0.0980450260
    ## ENSG00000004897    0.0000000000

``` r
mnn_merge<-cbind(mnn_corrected$corrected[[1]],mnn_corrected$corrected[[2]],mnn_corrected$corrected[[3]])
dim(mnn_merge)
```

    ## [1] 732 648

``` r
library("Rtsne")
set.seed(1)
tsne.out_mnn<-Rtsne(t(mnn_merge),initial_dims=10,verbose=TRUE,perplexity=30)
```

    ## Read the 648 x 10 data matrix successfully!
    ## Using no_dims = 2, perplexity = 30.000000, and theta = 0.500000
    ## Computing input similarities...
    ## Normalizing input...
    ## Building tree...
    ##  - point 0 of 648
    ## Done in 0.06 seconds (sparsity = 0.192006)!
    ## Learning embedding...
    ## Iteration 50: error is 63.614562 (50 iterations in 0.38 seconds)
    ## Iteration 100: error is 62.048376 (50 iterations in 0.22 seconds)
    ## Iteration 150: error is 62.042376 (50 iterations in 0.22 seconds)
    ## Iteration 200: error is 62.042389 (50 iterations in 0.21 seconds)
    ## Iteration 250: error is 62.042523 (50 iterations in 0.21 seconds)
    ## Iteration 300: error is 0.914109 (50 iterations in 0.19 seconds)
    ## Iteration 350: error is 0.834585 (50 iterations in 0.18 seconds)
    ## Iteration 400: error is 0.820615 (50 iterations in 0.17 seconds)
    ## Iteration 450: error is 0.811635 (50 iterations in 0.17 seconds)
    ## Iteration 500: error is 0.804897 (50 iterations in 0.17 seconds)
    ## Iteration 550: error is 0.802146 (50 iterations in 0.17 seconds)
    ## Iteration 600: error is 0.799145 (50 iterations in 0.17 seconds)
    ## Iteration 650: error is 0.796975 (50 iterations in 0.18 seconds)
    ## Iteration 700: error is 0.794473 (50 iterations in 0.17 seconds)
    ## Iteration 750: error is 0.792337 (50 iterations in 0.17 seconds)
    ## Iteration 800: error is 0.790972 (50 iterations in 0.17 seconds)
    ## Iteration 850: error is 0.789396 (50 iterations in 0.17 seconds)
    ## Iteration 900: error is 0.788083 (50 iterations in 0.18 seconds)
    ## Iteration 950: error is 0.787897 (50 iterations in 0.17 seconds)
    ## Iteration 1000: error is 0.788348 (50 iterations in 0.17 seconds)
    ## Fitting performed in 3.84 seconds.

``` r
plot(tsne.out_mnn$Y,main="ILC, tSNE PLOT, MNN CORRECTION",col=factor(M$Donor),xlab="tSNE1",ylab="tSNE2")
legend("topleft",levels(factor(M$Donor)),cex=1,fill=brewer.pal(length(levels(factor(M$Donor))),'Dark2'),inset=0.02)
```

![](batch_analysis_files/figure-markdown_github/MNN-1.png)

### Correcting for Batch-Effects with SCMAP

Now we are going to use a method called SCMAP for integration of two scRNAseq data sets which can be view as a correction for batch-effects. The idea of the method is to use one scRNAseq data set as a reference (e.g. a Human Cell Atlas data can be considered as a reference) and project cells from our particular scRNAseq experiment onto the reference data in order to determine the cell types in our scRNAseq data.

For simplicity we will consider cells from T74 and T75 donors as a reference and project cells from T86 onto the reference. Let us start with creating a SingleCellExperiment object for the reference data set.

``` r
library("SingleCellExperiment")
library("scmap")
ref_df<-subset(D,select=colnames(D)[grepl("T74|T75",colnames(D))])
ref_annot<-M[grepl("T74|T75",as.character(M$Donor)),]
ref_annot_scmap<-data.frame(cell_type1=ref_annot$Celltype)
rownames(ref_annot_scmap)<-rownames(ref_annot)
head(ref_annot_scmap)  
```

    ##                 cell_type1
    ## T74_P1_A9_ILC1        ILC1
    ## T74_P1_B4_NK            NK
    ## T74_P1_B7_ILC2        ILC2
    ## T74_P1_B9_NK            NK
    ## T74_P1_D10_ILC2       ILC2
    ## T74_P1_E1_ILC3        ILC3

``` r
sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(ref_df)), colData = ref_annot_scmap)
logcounts(sce) <- log2(normcounts(sce) + 1)
# use gene names as feature symbols
rowData(sce)$feature_symbol <- rownames(sce)
isSpike(sce, "ERCC") <- grepl("^ERCC-", rownames(sce))
# remove features with duplicated names
sce <- sce[!duplicated(rownames(sce)), ]
sce
```

    ## class: SingleCellExperiment 
    ## dim: 732 348 
    ## metadata(0):
    ## assays(2): normcounts logcounts
    ## rownames(732): ENSG00000003056 ENSG00000003402 ...
    ##   ERCC_29.296875:mix1_43.9453125:mix2
    ##   ERCC_468.75:mix1_117.1875:mix2
    ## rowData names(1): feature_symbol
    ## colnames(348): T74_P1_A9_ILC1 T74_P1_B4_NK ... T75_P3_H9_ILC2
    ##   T75_P4_H9_ILC1
    ## colData names(1): cell_type1
    ## reducedDimNames(0):
    ## spikeNames(1): ERCC

Now we will select most informative genes that provide cell clustering in the reference data set:

``` r
sce <- selectFeatures(sce, suppress_plot = FALSE)
```

    ## Warning in linearModel(object, n_features): Your object does not contain
    ## counts() slot. Dropouts were calculated using logcounts() slot...

![](batch_analysis_files/figure-markdown_github/Select%20Features-1.png)

``` r
table(rowData(sce)$scmap_features)
```

    ## 
    ## FALSE  TRUE 
    ##   232   500

Here we index reference data set clustering using medians of genes for each cluster:

``` r
sce <- indexCluster(sce)
head(metadata(sce)$scmap_cluster_index)
```

    ##                     ILC1       NK     ILC2 ILC3
    ## ENSG00000003056 1.000000 0.000000 0.000000    0
    ## ENSG00000003402 1.584963 6.194009 1.000000    1
    ## ENSG00000005955 0.000000 0.000000 0.000000    1
    ## ENSG00000006114 0.000000 1.000000 0.000000    0
    ## ENSG00000008513 3.700440 1.000000 1.000000    0
    ## ENSG00000008952 5.832890 1.000000 2.584963    1

``` r
heatmap(as.matrix(metadata(sce)$scmap_cluster_index))
```

![](batch_analysis_files/figure-markdown_github/Heatmap-1.png)

Lets us now create a SCE object for the test data set tp be projected onto the reference:

``` r
test_df<-subset(D,select=colnames(D)[grepl("T86",colnames(D))])
sce_test <- SingleCellExperiment(assays = list(normcounts = as.matrix(test_df)))
logcounts(sce_test) <- log2(normcounts(sce_test) + 1)
# use gene names as feature symbols
rowData(sce_test)$feature_symbol <- rownames(sce_test)
isSpike(sce_test, "ERCC") <- grepl("^ERCC-", rownames(sce_test))
# remove features with duplicated names
sce_test <- sce_test[!duplicated(rownames(sce_test)), ]
sce_test
```

    ## class: SingleCellExperiment 
    ## dim: 732 300 
    ## metadata(0):
    ## assays(2): normcounts logcounts
    ## rownames(732): ENSG00000003056 ENSG00000003402 ...
    ##   ERCC_29.296875:mix1_43.9453125:mix2
    ##   ERCC_468.75:mix1_117.1875:mix2
    ## rowData names(1): feature_symbol
    ## colnames(300): T86_P1_A1_ILC3 T86_P1_A10_ILC3 ... T86_P4_H5_ILC3
    ##   T86_P4_H8_ILC3
    ## colData names(0):
    ## reducedDimNames(0):
    ## spikeNames(1): ERCC

Now we can project the test data set onto the refernce:

``` r
scmapCluster_results <- scmapCluster(
  projection = sce_test, threshold = 0.3,
  index_list = list(
    reference = metadata(sce)$scmap_cluster_index
  )
)
head(scmapCluster_results$scmap_cluster_labs,20)
```

    ##       reference   
    ##  [1,] "unassigned"
    ##  [2,] "ILC3"      
    ##  [3,] "ILC3"      
    ##  [4,] "unassigned"
    ##  [5,] "unassigned"
    ##  [6,] "unassigned"
    ##  [7,] "ILC3"      
    ##  [8,] "ILC3"      
    ##  [9,] "unassigned"
    ## [10,] "unassigned"
    ## [11,] "unassigned"
    ## [12,] "unassigned"
    ## [13,] "unassigned"
    ## [14,] "NK"        
    ## [15,] "unassigned"
    ## [16,] "unassigned"
    ## [17,] "unassigned"
    ## [18,] "ILC3"      
    ## [19,] "unassigned"
    ## [20,] "unassigned"

``` r
head(scmapCluster_results$scmap_cluster_siml,20)
```

    ##        reference
    ##  [1,] 0.17298459
    ##  [2,] 0.44572133
    ##  [3,] 0.49633172
    ##  [4,] 0.09672277
    ##  [5,]         NA
    ##  [6,] 0.15442758
    ##  [7,] 0.49557526
    ##  [8,] 0.43295443
    ##  [9,] 0.16884769
    ## [10,] 0.13381652
    ## [11,] 0.13938067
    ## [12,] 0.12702673
    ## [13,] 0.16973052
    ## [14,] 0.56767457
    ## [15,] 0.19470757
    ## [16,] 0.18266300
    ## [17,] 0.14917527
    ## [18,] 0.51367106
    ## [19,] 0.17732457
    ## [20,] 0.11778121

``` r
results<-data.frame(CELS=colnames(test_df),ASSIGNED_LABEL=as.vector(scmapCluster_results$scmap_cluster_labs),TRUE_LABEL=matrix(unlist(strsplit(colnames(test_df),"_")),ncol=4,byrow=TRUE)[,4],SIMILARITY=as.vector(scmapCluster_results$scmap_cluster_siml))
head(results,20)
```

    ##               CELS ASSIGNED_LABEL TRUE_LABEL SIMILARITY
    ## 1   T86_P1_A1_ILC3     unassigned       ILC3 0.17298459
    ## 2  T86_P1_A10_ILC3           ILC3       ILC3 0.44572133
    ## 3  T86_P1_A12_ILC3           ILC3       ILC3 0.49633172
    ## 4   T86_P1_A3_ILC3     unassigned       ILC3 0.09672277
    ## 5   T86_P1_A5_ILC3     unassigned       ILC3         NA
    ## 6   T86_P1_A8_ILC3     unassigned       ILC3 0.15442758
    ## 7   T86_P1_B1_ILC3           ILC3       ILC3 0.49557526
    ## 8  T86_P1_B12_ILC3           ILC3       ILC3 0.43295443
    ## 9   T86_P1_B3_ILC3     unassigned       ILC3 0.16884769
    ## 10  T86_P1_B5_ILC3     unassigned       ILC3 0.13381652
    ## 11  T86_P1_B6_ILC3     unassigned       ILC3 0.13938067
    ## 12  T86_P1_B7_ILC3     unassigned       ILC3 0.12702673
    ## 13  T86_P1_B8_ILC3     unassigned       ILC3 0.16973052
    ## 14    T86_P1_B9_NK             NK         NK 0.56767457
    ## 15  T86_P1_C1_ILC3     unassigned       ILC3 0.19470757
    ## 16 T86_P1_C10_ILC3     unassigned       ILC3 0.18266300
    ## 17 T86_P1_C11_ILC3     unassigned       ILC3 0.14917527
    ## 18 T86_P1_C12_ILC3           ILC3       ILC3 0.51367106
    ## 19  T86_P1_C2_ILC3     unassigned       ILC3 0.17732457
    ## 20  T86_P1_C3_ILC3     unassigned       ILC3 0.11778121

We can see that there are 127 "unassigned" cells, i.e. where KNN classifier fails to assign their class. Let as calculate the accuracy of assignment for the assigned cells:

``` r
results_assigned<-results[as.character(results$ASSIGNED_LABEL)!="unassigned",]
head(results_assigned,20)
```

    ##               CELS ASSIGNED_LABEL TRUE_LABEL SIMILARITY
    ## 2  T86_P1_A10_ILC3           ILC3       ILC3  0.4457213
    ## 3  T86_P1_A12_ILC3           ILC3       ILC3  0.4963317
    ## 7   T86_P1_B1_ILC3           ILC3       ILC3  0.4955753
    ## 8  T86_P1_B12_ILC3           ILC3       ILC3  0.4329544
    ## 14    T86_P1_B9_NK             NK         NK  0.5676746
    ## 18 T86_P1_C12_ILC3           ILC3       ILC3  0.5136711
    ## 23  T86_P1_C6_ILC3           ILC3       ILC3  0.4655949
    ## 28 T86_P1_D10_ILC3           ILC2       ILC3  0.3970456
    ## 29   T86_P1_D11_NK             NK         NK  0.5262334
    ## 30 T86_P1_D12_ILC3           ILC3       ILC3  0.5097175
    ## 34  T86_P1_D6_ILC3           ILC3       ILC3  0.4664650
    ## 38 T86_P1_E10_ILC3           ILC3       ILC3  0.4750463
    ## 40 T86_P1_E12_ILC3           ILC3       ILC3  0.4968623
    ## 41  T86_P1_E2_ILC3           ILC3       ILC3  0.4253116
    ## 42  T86_P1_E3_ILC3           ILC3       ILC3  0.4706919
    ## 45    T86_P1_E6_NK             NK         NK  0.5219235
    ## 46    T86_P1_E7_NK             NK         NK  0.5405412
    ## 48  T86_P1_E9_ILC3           ILC3       ILC3  0.4489822
    ## 50 T86_P1_F10_ILC3           ILC3       ILC3  0.4821195
    ## 51 T86_P1_F11_ILC3           ILC3       ILC3  0.4667251

``` r
table(results_assigned$ASSIGNED_LABEL,results_assigned$TRUE_LABEL)
```

    ##             
    ##              ILC1 ILC2 ILC3 NK
    ##   ILC1         65   26    0  0
    ##   ILC2          0    2    2  0
    ##   ILC3          0    0   60  0
    ##   NK            0    0    0 18
    ##   unassigned    0    0    0  0

``` r
sum(as.character(results_assigned$ASSIGNED_LABEL)==as.character(results_assigned$TRUE_LABEL))/dim(results_assigned)[1]
```

    ## [1] 0.8381503

We conclude that the accuracy of assignment is 84% which is not fantastic taking into account that SCMAP failed assignment of almost a half of the cells in the test data set.

### Correcting for Batch Effects with Seurat

Now we will use Seurat “alignment” method which is based on Canonical Correlation Analysis (CCA). CCA extract common sources of variation and facilitate identification of common cell types between two data sets.

``` r
library("Seurat")
ref_df<-subset(D,select=colnames(D)[grepl("T74|T75",colnames(D))])
test_df<-subset(D,select=colnames(D)[grepl("T86",colnames(D))])

# Set up reference object
ref_seurat <- CreateSeuratObject(raw.data = ref_df, project = "ref", min.cells = 5)
ref_seurat@meta.data$stim <- "ref"
ref_seurat <- FilterCells(ref_seurat, subset.names = "nGene", low.thresholds = 100, high.thresholds = Inf)
ref_seurat <- NormalizeData(ref_seurat)
ref_seurat <- ScaleData(ref_seurat, display.progress = T)
```

    ## [1] "Scaling data matrix"
    ## 
      |                                                                       
      |                                                                 |   0%
      |                                                                       
      |=================================================================| 100%

``` r
# Set up T75 object
test_seurat <- CreateSeuratObject(raw.data = test_df, project = "test", min.cells = 5)
test_seurat@meta.data$stim <- "test"
test_seurat <- FilterCells(test_seurat, subset.names = "nGene", low.thresholds = 100, high.thresholds = Inf)
test_seurat <- NormalizeData(test_seurat)
test_seurat <- ScaleData(test_seurat, display.progress = T)
```

    ## [1] "Scaling data matrix"
    ## 
      |                                                                       
      |                                                                 |   0%
      |                                                                       
      |=================================================================| 100%

``` r
# Gene selection for input to CCA
ref_seurat <- FindVariableGenes(ref_seurat, do.plot = F)
test_seurat <- FindVariableGenes(test_seurat, do.plot = F)
g.1 <- head(rownames(ref_seurat@hvg.info), 1000)
g.2 <- head(rownames(test_seurat@hvg.info), 1000)
genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(ref_seurat@scale.data))
genes.use <- intersect(genes.use, rownames(test_seurat@scale.data))

#CCA
seurat_combined <- RunCCA(ref_seurat, test_seurat , genes.use = genes.use, num.cc = 30)
```

    ## [1] "Scaling data matrix"
    ## 
      |                                                                       
      |                                                                 |   0%
      |                                                                       
      |=================================================================| 100%

``` r
p1 <- DimPlot(object = seurat_combined,reduction.use="cca", group.by = "stim", pt.size=0.5, do.return=TRUE)
p2 <- VlnPlot(object = seurat_combined, features.plot = "CC1", group.by = "stim", do.return = TRUE)
plot_grid(p1,p2)
```

![](batch_analysis_files/figure-markdown_github/Seurat-1.png)

``` r
PrintDim(object = seurat_combined, reduction.type = "cca", dims.print = 1:2, genes.print = 10)
```

    ## [1] "CC1"
    ##  [1] "ENSG00000162594" "ENSG00000157404" "ENSG00000110002"
    ##  [4] "ENSG00000104951" "ENSG00000011600" "ENSG00000204482"
    ##  [7] "ENSG00000158869" "ENSG00000102524" "ENSG00000160593"
    ## [10] "ENSG00000198574"
    ## [1] ""
    ##  [1] "ENSG00000142227" "ENSG00000170989" "ENSG00000213684"
    ##  [4] "ENSG00000197956" "ENSG00000111913" "ENSG00000111716"
    ##  [7] "ENSG00000160255" "ENSG00000226958" "ENSG00000141232"
    ## [10] "ENSG00000179144"
    ## [1] ""
    ## [1] ""
    ## [1] "CC2"
    ##  [1] "ENSG00000115523" "ENSG00000158050" "ENSG00000134545"
    ##  [4] "ENSG00000100385" "ENSG00000223336" "ENSG00000011600"
    ##  [7] "ENSG00000223001" "ENSG00000126264" "ENSG00000145287"
    ## [10] "ENSG00000222414"
    ## [1] ""
    ##  [1] "ENSG00000135074" "ENSG00000215788" "ENSG00000197956"
    ##  [4] "ENSG00000112486" "ENSG00000160593" "ENSG00000138468"
    ##  [7] "ENSG00000104951" "ENSG00000135426" "ENSG00000171867"
    ## [10] "ENSG00000167601"
    ## [1] ""
    ## [1] ""

``` r
p3 <- MetageneBicorPlot(seurat_combined, grouping.var = "stim", dims.eval = 1:30, display.progress = FALSE)
```

    ## `geom_smooth()` using method = 'loess'

![](batch_analysis_files/figure-markdown_github/Seurat-2.png)

``` r
DimHeatmap(object = seurat_combined, reduction.type = "cca", cells.use = 500, dim.use = 1:9, do.balanced = TRUE)
```

![](batch_analysis_files/figure-markdown_github/Seurat-3.png)

``` r
seurat_combined <- AlignSubspace(seurat_combined, reduction.type = "cca", grouping.var = "stim", dims.align = 1:20)
```

    ## [1] "Scaling data matrix"
    ## 
      |                                                                       
      |                                                                 |   0%
      |                                                                       
      |=================================================================| 100%
    ## [1] "Scaling data matrix"
    ## 
      |                                                                       
      |                                                                 |   0%
      |                                                                       
      |=================================================================| 100%

``` r
# t-SNE and Clustering
seurat_combined <- RunTSNE(seurat_combined, reduction.use = "cca.aligned", dims.use = 1:20, do.fast = T)
seurat_combined <- FindClusters(seurat_combined, reduction.type = "cca.aligned", resolution = 0.6, dims.use = 1:20,print.output = FALSE)

# Visualization
p1 <- TSNEPlot(seurat_combined, do.return = T, pt.size = 0.5, group.by = "stim")
p2 <- TSNEPlot(seurat_combined, do.label = T, do.return = T, pt.size = 0.5)
plot_grid(p1, p2)
```

![](batch_analysis_files/figure-markdown_github/Seurat-4.png)

Again we see good mixing of cells from the refernce and test data sets.

### Session Information

Finally here is the details on the system on which this document was compiled:

``` r
sessionInfo()
```

    ## R version 3.4.3 (2017-11-30)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 14.04.5 LTS
    ## 
    ## Matrix products: default
    ## BLAS: /usr/lib/libblas/libblas.so.3.0
    ## LAPACK: /usr/lib/lapack/liblapack.so.3.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=sv_SE.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=sv_SE.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=sv_SE.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=sv_SE.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] bindrcpp_0.2               scmap_1.1.5               
    ##  [3] scran_1.6.9                SingleCellExperiment_1.0.0
    ##  [5] SummarizedExperiment_1.8.1 DelayedArray_0.4.1        
    ##  [7] Biobase_2.38.0             GenomicRanges_1.30.1      
    ##  [9] GenomeInfoDb_1.14.0        IRanges_2.12.0            
    ## [11] S4Vectors_0.16.0           BiocGenerics_0.24.0       
    ## [13] sva_3.26.0                 BiocParallel_1.12.0       
    ## [15] genefilter_1.60.0          mgcv_1.8-23               
    ## [17] nlme_3.1-131               Rtsne_0.13                
    ## [19] RColorBrewer_1.1-2         mixOmics_6.3.1            
    ## [21] lattice_0.20-35            MASS_7.3-48               
    ## [23] matrixStats_0.53.1         gplots_3.0.1              
    ## [25] Seurat_2.3.0               Matrix_1.2-11             
    ## [27] cowplot_0.9.2              ggplot2_2.2.1             
    ## [29] rmarkdown_1.9              knitr_1.20                
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] shinydashboard_0.6.1   R.utils_2.6.0          tidyselect_0.2.3      
    ##   [4] RSQLite_2.0            AnnotationDbi_1.40.0   htmlwidgets_1.0       
    ##   [7] grid_3.4.3             trimcluster_0.1-2      ranger_0.9.0          
    ##  [10] munsell_0.4.3          codetools_0.2-15       ica_1.0-1             
    ##  [13] DT_0.4                 statmod_1.4.30         withr_2.1.1.9000      
    ##  [16] colorspace_1.3-2       rstudioapi_0.7         ROCR_1.0-7            
    ##  [19] robustbase_0.92-8      dtw_1.18-1             dimRed_0.1.0          
    ##  [22] labeling_0.3           lars_1.2               tximport_1.6.0        
    ##  [25] GenomeInfoDbData_1.0.0 mnormt_1.5-5           bit64_0.9-7           
    ##  [28] rhdf5_2.22.0           rprojroot_1.3-2        ipred_0.9-6           
    ##  [31] randomForest_4.6-12    diptest_0.75-7         R6_2.2.2              
    ##  [34] ggbeeswarm_0.6.0       VGAM_1.0-5             locfit_1.5-9.1        
    ##  [37] flexmix_2.3-14         DRR_0.0.3              bitops_1.0-6          
    ##  [40] assertthat_0.2.0       SDMTools_1.1-221       scales_0.5.0.9000     
    ##  [43] nnet_7.3-12            beeswarm_0.2.3         gtable_0.2.0          
    ##  [46] ddalpha_1.3.1.1        timeDate_3043.102      rlang_0.2.0.9000      
    ##  [49] CVST_0.2-1             scatterplot3d_0.3-41   RcppRoll_0.2.2        
    ##  [52] splines_3.4.3          lazyeval_0.2.1         ModelMetrics_1.1.0    
    ##  [55] acepack_1.4.1          broom_0.4.3            checkmate_1.8.5       
    ##  [58] rgl_0.99.9             yaml_2.1.18            reshape2_1.4.3        
    ##  [61] crosstalk_1.0.0        backports_1.1.2        httpuv_1.3.6.2        
    ##  [64] Hmisc_4.1-1            caret_6.0-79           tools_3.4.3           
    ##  [67] lava_1.6               psych_1.7.8            proxy_0.4-22          
    ##  [70] dynamicTreeCut_1.63-1  ggridges_0.5.0         Rcpp_0.12.16          
    ##  [73] plyr_1.8.4             progress_1.1.2         zlibbioc_1.24.0       
    ##  [76] base64enc_0.1-3        purrr_0.2.4            RCurl_1.95-4.10       
    ##  [79] prettyunits_1.0.2      rpart_4.1-12           viridis_0.5.0         
    ##  [82] pbapply_1.3-4          zoo_1.8-1              sfsmisc_1.1-1         
    ##  [85] cluster_2.0.6          magrittr_1.5           data.table_1.10.4-3   
    ##  [88] RSpectra_0.12-0        lmtest_0.9-35          RANN_2.5.1            
    ##  [91] mvtnorm_1.0-7          fitdistrplus_1.0-9     mime_0.5              
    ##  [94] evaluate_0.10.1        xtable_1.8-2           XML_3.98-1.1          
    ##  [97] mclust_5.4             gridExtra_2.3          biomaRt_2.34.2        
    ## [100] scater_1.6.3           compiler_3.4.3         ellipse_0.4.1         
    ## [103] tibble_1.4.2           KernSmooth_2.23-15     R.oo_1.21.0           
    ## [106] htmltools_0.3.6        segmented_0.5-3.0      corpcor_1.6.9         
    ## [109] Formula_1.2-2          snow_0.4-2             tidyr_0.8.0           
    ## [112] tclust_1.3-1           lubridate_1.7.2        DBI_0.7               
    ## [115] diffusionMap_1.1-0     fpc_2.1-11             R.methodsS3_1.7.1     
    ## [118] gdata_2.18.0           metap_0.8              bindr_0.1             
    ## [121] gower_0.1.2            igraph_1.1.2           pkgconfig_2.0.1       
    ## [124] sn_1.5-1               numDeriv_2016.8-1      foreign_0.8-69        
    ## [127] recipes_0.1.2          foreach_1.4.4          rARPACK_0.11-0        
    ## [130] annotate_1.56.1        vipor_0.4.5            XVector_0.18.0        
    ## [133] prodlim_1.6.1          stringr_1.3.0          digest_0.6.15         
    ## [136] tsne_0.1-3             htmlTable_1.11.2       edgeR_3.20.8          
    ## [139] googleVis_0.6.2        kernlab_0.9-25         shiny_1.0.5           
    ## [142] gtools_3.5.0           modeltools_0.2-21      rjson_0.2.15          
    ## [145] jsonlite_1.5           viridisLite_0.3.0      limma_3.34.9          
    ## [148] pillar_1.2.0           httr_1.3.1             DEoptimR_1.0-8        
    ## [151] survival_2.41-3        glue_1.2.0             FNN_1.1               
    ## [154] png_0.1-7              prabclus_2.2-6         iterators_1.0.9       
    ## [157] bit_1.1-12             class_7.3-14           stringi_1.1.7         
    ## [160] mixtools_1.1.0         blob_1.1.0             doSNOW_1.0.16         
    ## [163] latticeExtra_0.6-28    caTools_1.17.1         memoise_1.1.0         
    ## [166] dplyr_0.7.4            e1071_1.6-8            irlba_2.3.2           
    ## [169] ape_5.1
