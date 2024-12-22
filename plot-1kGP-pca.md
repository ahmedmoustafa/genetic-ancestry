Plot Principal Component Analysis (PCA) of 1kGP
================
Ahmed Moustafa
22 December, 2024 18:02 +0200

## Loading Libraries

``` r
library(tidyverse)
library(tidylog)
library(Rtsne)
```

## Data Loading: Eigenvalues

``` r
# Load the eigenvalues
eigenval = read_tsv("working/1kGP_pca.eigenval", col_names = c("value"))
```

    ## Rows: 10 Columns: 1
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## dbl (1): value
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
head(eigenval)
```

    ## # A tibble: 6 × 1
    ##    value
    ##    <dbl>
    ## 1 347.  
    ## 2 151.  
    ## 3  42.6 
    ## 4  31.4 
    ## 5   5.56
    ## 6   5.19

## Calculating Variance Explained by Each Principal Component

``` r
variance_explained = tibble(pc = paste0("PC", 1:nrow(eigenval)), var = eigenval$value / sum(eigenval))
head(variance_explained)
```

    ## # A tibble: 6 × 2
    ##   pc        var
    ##   <chr>   <dbl>
    ## 1 PC1   0.581  
    ## 2 PC2   0.253  
    ## 3 PC3   0.0713 
    ## 4 PC4   0.0526 
    ## 5 PC5   0.00930
    ## 6 PC6   0.00868

## Plotting the Variance Explained by Principal Components

``` r
ggplot(variance_explained %>% top_n(3)) + 
  geom_bar(aes(x = pc, y = var), stat = "identity") + 
  theme_light() +
  labs (x = "PC", y = "Variance Explained") +
  scale_y_continuous(labels = scales::percent)
```

    ## Selecting by var
    ## top_n: removed 7 rows (70%), 3 rows remaining

![](plot-1kGP-pca_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

## Data Loading: Eigenvectors

``` r
# Load the eigenvectors
eigenvec = read.table("working/1kGP_pca.eigenvec", sep = " ")
dim(eigenvec)
```

    ## [1] 3202   12

``` r
# Display the first few rows of the PCA data to inspect the structure
head(eigenvec)
```

    ##        V1      V2          V3        V4          V5        V6           V7
    ## 1 HG00096 HG00096 -0.00979116 0.0246988  0.00300170 0.0175045 -0.000177201
    ## 2 HG00097 HG00097 -0.00860299 0.0246735  0.00216549 0.0173607 -0.005536730
    ## 3 HG00099 HG00099 -0.00940134 0.0242093  0.00428509 0.0207030 -0.007276840
    ## 4 HG00100 HG00100 -0.00974342 0.0227575 -0.00062867 0.0176841 -0.008215210
    ## 5 HG00101 HG00101 -0.00949396 0.0236118  0.00404194 0.0183175  0.000795819
    ## 6 HG00102 HG00102 -0.01004650 0.0232639  0.00285737 0.0144682 -0.001366250
    ##            V8          V9          V10          V11         V12
    ## 1 -0.02154430 -0.00710697  4.00740e-03 -0.011880200 -0.00355934
    ## 2  0.00307645 -0.00263790  9.36712e-05  0.007745490  0.00669956
    ## 3 -0.01732370 -0.00154378 -4.40629e-03 -0.004739370  0.00448401
    ## 4 -0.01135090 -0.00653149  4.08700e-04 -0.004569230 -0.00936117
    ## 5 -0.01490100 -0.00882091  7.70109e-03  0.002463610  0.02843140
    ## 6 -0.01884990  0.00368535 -1.34158e-02 -0.000063343 -0.00326615

``` r
pca = tibble(sample = eigenvec$V1, PC1 = eigenvec$V3, PC2 = eigenvec$V4, PC3 = eigenvec$V5)
head(pca)
```

    ## # A tibble: 6 × 4
    ##   sample       PC1    PC2       PC3
    ##   <chr>      <dbl>  <dbl>     <dbl>
    ## 1 HG00096 -0.00979 0.0247  0.00300 
    ## 2 HG00097 -0.00860 0.0247  0.00217 
    ## 3 HG00099 -0.00940 0.0242  0.00429 
    ## 4 HG00100 -0.00974 0.0228 -0.000629
    ## 5 HG00101 -0.00949 0.0236  0.00404 
    ## 6 HG00102 -0.0100  0.0233  0.00286

## Plotting PCA: Visualization of Population Structure

``` r
ggplot(pca) + 
  geom_jitter(aes(x = PC1, y = PC2), alpha = 0.7) + 
  theme_light()
```

![](plot-1kGP-pca_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

## Enhancing the PCA Plot: Color by Population

``` r
samples = read_tsv("data/populations.tsv")
```

    ## Rows: 4978 Columns: 7
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (7): sample, sex, biosample, population, population_name, superpopulatio...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
head(samples)
```

    ## # A tibble: 6 × 7
    ##   sample  sex    biosample   population population_name superpopulation
    ##   <chr>   <chr>  <chr>       <chr>      <chr>           <chr>          
    ## 1 HG00271 male   SAME123417  FIN        Finnish         EUR            
    ## 2 HG00276 female SAME123424  FIN        Finnish         EUR            
    ## 3 HG00288 female SAME1839246 FIN        Finnish         EUR            
    ## 4 HG00290 male   SAME1839057 FIN        Finnish         EUR            
    ## 5 HG00303 male   SAME1840115 FIN        Finnish         EUR            
    ## 6 HG00308 male   SAME124161  FIN        Finnish         EUR            
    ## # ℹ 1 more variable: superpopulation_name <chr>

``` r
pca2 = pca %>% inner_join(samples)
```

    ## Joining with `by = join_by(sample)`
    ## inner_join: added 6 columns (sex, biosample, population, population_name,
    ## superpopulation, …)
    ## > rows only in x ( 0)
    ## > rows only in samples (1,776)
    ## > matched rows 3,202
    ## > =======
    ## > rows total 3,202

``` r
head(pca2)
```

    ## # A tibble: 6 × 10
    ##   sample       PC1    PC2       PC3 sex    biosample  population population_name
    ##   <chr>      <dbl>  <dbl>     <dbl> <chr>  <chr>      <chr>      <chr>          
    ## 1 HG00096 -0.00979 0.0247  0.00300  male   SAME123268 GBR        British        
    ## 2 HG00097 -0.00860 0.0247  0.00217  female SAME123267 GBR        British        
    ## 3 HG00099 -0.00940 0.0242  0.00429  female SAME123271 GBR        British        
    ## 4 HG00100 -0.00974 0.0228 -0.000629 female SAME125154 GBR        British        
    ## 5 HG00101 -0.00949 0.0236  0.00404  male   SAME125153 GBR        British        
    ## 6 HG00102 -0.0100  0.0233  0.00286  female SAME123945 GBR        British        
    ## # ℹ 2 more variables: superpopulation <chr>, superpopulation_name <chr>

``` r
ggplot(pca2) + 
  geom_jitter(aes(x = PC1, y = PC2, color = superpopulation), alpha = 0.7) + 
  theme_light()
```

![](plot-1kGP-pca_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
ggplot(pca2) + 
  geom_jitter(aes(x = PC1, y = PC2, color = population), alpha = 0.7) + 
  theme_light()
```

![](plot-1kGP-pca_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

## Alternative Clustering using t-SNE

``` r
n = 10
df = data.frame (eigenvec[, 3:(n+2)])
colnames(df) = paste0("PC", 1:n)
rownames(df) = eigenvec$V1
head(df)
```

    ##                 PC1       PC2         PC3       PC4          PC5         PC6
    ## HG00096 -0.00979116 0.0246988  0.00300170 0.0175045 -0.000177201 -0.02154430
    ## HG00097 -0.00860299 0.0246735  0.00216549 0.0173607 -0.005536730  0.00307645
    ## HG00099 -0.00940134 0.0242093  0.00428509 0.0207030 -0.007276840 -0.01732370
    ## HG00100 -0.00974342 0.0227575 -0.00062867 0.0176841 -0.008215210 -0.01135090
    ## HG00101 -0.00949396 0.0236118  0.00404194 0.0183175  0.000795819 -0.01490100
    ## HG00102 -0.01004650 0.0232639  0.00285737 0.0144682 -0.001366250 -0.01884990
    ##                 PC7          PC8          PC9        PC10
    ## HG00096 -0.00710697  4.00740e-03 -0.011880200 -0.00355934
    ## HG00097 -0.00263790  9.36712e-05  0.007745490  0.00669956
    ## HG00099 -0.00154378 -4.40629e-03 -0.004739370  0.00448401
    ## HG00100 -0.00653149  4.08700e-04 -0.004569230 -0.00936117
    ## HG00101 -0.00882091  7.70109e-03  0.002463610  0.02843140
    ## HG00102  0.00368535 -1.34158e-02 -0.000063343 -0.00326615

### Perform t-SNE

``` r
# Perform t-SNE
set.seed(123)  # Set a seed for reproducibility
tsne_result <- Rtsne(as.matrix(df), dims = 2, perplexity = 30, theta = 0.5, verbose = TRUE)
```

    ## Performing PCA
    ## Read the 3202 x 10 data matrix successfully!
    ## Using no_dims = 2, perplexity = 30.000000, and theta = 0.500000
    ## Computing input similarities...
    ## Building tree...
    ## Done in 0.14 seconds (sparsity = 0.036737)!
    ## Learning embedding...
    ## Iteration 50: error is 80.026927 (50 iterations in 0.20 seconds)
    ## Iteration 100: error is 64.730357 (50 iterations in 0.20 seconds)
    ## Iteration 150: error is 62.662983 (50 iterations in 0.19 seconds)
    ## Iteration 200: error is 61.832415 (50 iterations in 0.19 seconds)
    ## Iteration 250: error is 61.380744 (50 iterations in 0.19 seconds)
    ## Iteration 300: error is 1.752580 (50 iterations in 0.22 seconds)
    ## Iteration 350: error is 1.456704 (50 iterations in 0.17 seconds)
    ## Iteration 400: error is 1.307920 (50 iterations in 0.16 seconds)
    ## Iteration 450: error is 1.221664 (50 iterations in 0.17 seconds)
    ## Iteration 500: error is 1.173401 (50 iterations in 0.17 seconds)
    ## Iteration 550: error is 1.147737 (50 iterations in 0.17 seconds)
    ## Iteration 600: error is 1.131094 (50 iterations in 0.17 seconds)
    ## Iteration 650: error is 1.116694 (50 iterations in 0.17 seconds)
    ## Iteration 700: error is 1.106174 (50 iterations in 0.17 seconds)
    ## Iteration 750: error is 1.097210 (50 iterations in 0.19 seconds)
    ## Iteration 800: error is 1.088086 (50 iterations in 0.19 seconds)
    ## Iteration 850: error is 1.079247 (50 iterations in 0.17 seconds)
    ## Iteration 900: error is 1.072591 (50 iterations in 0.17 seconds)
    ## Iteration 950: error is 1.068788 (50 iterations in 0.16 seconds)
    ## Iteration 1000: error is 1.065359 (50 iterations in 0.17 seconds)
    ## Fitting performed in 3.60 seconds.

``` r
# The results are stored in tsne_result$Y
head(tsne_result$Y)
```

    ##          [,1]     [,2]
    ## [1,] 24.77744 21.75235
    ## [2,] 20.83018 17.98410
    ## [3,] 22.44079 21.63661
    ## [4,] 25.42296 16.83022
    ## [5,] 17.74731 22.45361
    ## [6,] 20.66063 20.02540

### Visualizing t-SNE Clustering

``` r
# Create a data frame for plotting
tsne_df = data.frame(TSNE1 = tsne_result$Y[, 1], 
                     TSNE2 = tsne_result$Y[, 2], 
                     sample = rownames(df)) %>% inner_join(samples)
```

    ## Joining with `by = join_by(sample)`
    ## inner_join: added 6 columns (sex, biosample, population, population_name,
    ## superpopulation, …)
    ## > rows only in x ( 0)
    ## > rows only in samples (1,776)
    ## > matched rows 3,202
    ## > =======
    ## > rows total 3,202

``` r
head(tsne_df)
```

    ##      TSNE1    TSNE2  sample    sex  biosample population population_name
    ## 1 24.77744 21.75235 HG00096   male SAME123268        GBR         British
    ## 2 20.83018 17.98410 HG00097 female SAME123267        GBR         British
    ## 3 22.44079 21.63661 HG00099 female SAME123271        GBR         British
    ## 4 25.42296 16.83022 HG00100 female SAME125154        GBR         British
    ## 5 17.74731 22.45361 HG00101   male SAME125153        GBR         British
    ## 6 20.66063 20.02540 HG00102 female SAME123945        GBR         British
    ##   superpopulation superpopulation_name
    ## 1             EUR    European Ancestry
    ## 2             EUR    European Ancestry
    ## 3             EUR    European Ancestry
    ## 4             EUR    European Ancestry
    ## 5             EUR    European Ancestry
    ## 6             EUR    European Ancestry

``` r
ggplot(tsne_df, aes(x = TSNE1, y = TSNE2, color = population), alpha = 0.7) +
  geom_point() +
  theme_light()
```

![](plot-1kGP-pca_files/figure-gfm/plot-tsne-by-population-1.png)<!-- -->

``` r
  labs(x = "t-SNE1", y = "t-SNE2")
```

    ## $x
    ## [1] "t-SNE1"
    ## 
    ## $y
    ## [1] "t-SNE2"
    ## 
    ## attr(,"class")
    ## [1] "labels"

``` r
ggplot(tsne_df, aes(x = TSNE1, y = TSNE2, color = superpopulation), alpha = 0.7) +
  geom_point() +
  theme_light()
```

![](plot-1kGP-pca_files/figure-gfm/plot-tsne-by-superpopulation-1.png)<!-- -->

``` r
  labs(x = "t-SNE1", y = "t-SNE2")
```

    ## $x
    ## [1] "t-SNE1"
    ## 
    ## $y
    ## [1] "t-SNE2"
    ## 
    ## attr(,"class")
    ## [1] "labels"

## Session Info

``` r
sessionInfo()
```

    ## R version 4.4.2 (2024-10-31)
    ## Platform: aarch64-apple-darwin24.1.0
    ## Running under: macOS Sequoia 15.2
    ## 
    ## Matrix products: default
    ## BLAS:   /opt/homebrew/Cellar/openblas/0.3.28/lib/libopenblasp-r0.3.28.dylib 
    ## LAPACK: /opt/homebrew/Cellar/r/4.4.2_2/lib/R/lib/libRlapack.dylib;  LAPACK version 3.12.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: Africa/Cairo
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] Rtsne_0.17      tidylog_1.1.0   lubridate_1.9.4 forcats_1.0.0  
    ##  [5] stringr_1.5.1   dplyr_1.1.4     purrr_1.0.2     readr_2.1.5    
    ##  [9] tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.1   tidyverse_2.0.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] utf8_1.2.4        generics_0.1.3    stringi_1.8.4     hms_1.1.3        
    ##  [5] digest_0.6.37     magrittr_2.0.3    evaluate_1.0.1    grid_4.4.2       
    ##  [9] timechange_0.3.0  fastmap_1.2.0     scales_1.3.0      cli_3.6.3        
    ## [13] rlang_1.1.4       crayon_1.5.3      bit64_4.5.2       munsell_0.5.1    
    ## [17] withr_3.0.2       yaml_2.3.10       tools_4.4.2       parallel_4.4.2   
    ## [21] tzdb_0.4.0        colorspace_2.1-1  vctrs_0.6.5       R6_2.5.1         
    ## [25] lifecycle_1.0.4   bit_4.5.0.1       clisymbols_1.2.0  vroom_1.6.5      
    ## [29] pkgconfig_2.0.3   pillar_1.10.0     gtable_0.3.6      glue_1.7.0       
    ## [33] Rcpp_1.0.13       highr_0.11        xfun_0.47         tidyselect_1.2.1 
    ## [37] rstudioapi_0.17.1 knitr_1.48        farver_2.1.2      htmltools_0.5.8.1
    ## [41] rmarkdown_2.28    labeling_0.4.3    compiler_4.4.2
