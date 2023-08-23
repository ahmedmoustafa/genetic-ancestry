Ancestral SNPs Retrieval and Cleaning
================
Ahmed Moustafa
23 August, 2023 09:46 +0300

## Introduction

This notebook is to fetch the details of the 10,000 SNPs defined by
[NCBI’s GRAF](https://github.com/ncbi/graf) as the 10,000 ancestral
SNPs. Once fetched, the data is processed to filter out SNPs and retain
only the biallelic ones based on the information retrieved. The
resultant dataset provides a focused set of SNPs that can be used for
various analyses related to ancestry and population genetics.

### Data

The [rsIDs of the 10,000 SNPs](data/fingerprinting_snps.txt) were
extracted from GRAF’s
[G1000FpGeno](https://github.com/ncbi/graf/tree/master/data) collection
using [PLINK](https://www.cog-genomics.org/plink/)’s `--recode` option.

### References

- [GRAF GitHub](https://github.com/ncbi/graf)
- [GRAF-pop: A Fast Distance-Based Method To Infer Subject Ancestry from
  Multiple Genotype Datasets Without Principal Components
  Analysis](https://pubmed.ncbi.nlm.nih.gov/31151998/)

## Load Required Libraries

The following libraries are needed for SNP data retrieval and
transformation:

``` r
library(rsnps)
library(tidyverse)
library(tidylog)
library(tictoc)
```

## Load Data

Load a list of rsIDs (reference SNP identifiers) from a provided text
file.

``` r
# Load the list of rsIDs from a text file
rsids = readLines("data/fingerprinting_snps.txt")
head(rsids)
```

    ## [1] "rs2887286" "rs6685064" "rs2840528" "rs3890745" "rs1798246" "rs1181875"

``` r
length(rsids)
```

    ## [1] 10000

``` r
batch_size = 1000
batch_number = 1
```

## Define Function to Fetch SNP Data

A function is defined to fetch SNP data for a batch of rsIDs from the
NCBI database. The batch processing is done to avoid overwhelming the
server with too many requests at once.

``` r
# Define a function to fetch SNP data for a batch of rsIDs
fetch_batch_snp_data = function(batch) {
    tryCatch({
        print (paste0("Processing batch number: ", batch_number))
        batch_number <<- batch_number + 1
        tic()
        result = tibble(ncbi_snp_query(batch))
        toc()
        return (result)
    }, error = function(e) {
        message(paste0("Error fetching data: ", e$message))
        NULL
  })
}
```

The rsIDs are split into batches, and the SNP data for each batch is
fetched using the previously defined function.

``` r
# Split rsIDs into batches of 1,000 and fetch the SNP data
snp_data = map_dfr(split(rsids, ceiling(seq_along(rsids) / batch_size)), fetch_batch_snp_data)
```

    ## [1] "Processing batch number: 1"
    ## 244.427 sec elapsed
    ## [1] "Processing batch number: 2"
    ## 268.934 sec elapsed
    ## [1] "Processing batch number: 3"
    ## 269.18 sec elapsed
    ## [1] "Processing batch number: 4"
    ## 281.656 sec elapsed
    ## [1] "Processing batch number: 5"
    ## 281.486 sec elapsed
    ## [1] "Processing batch number: 6"
    ## 276.075 sec elapsed
    ## [1] "Processing batch number: 7"
    ## 285.83 sec elapsed
    ## [1] "Processing batch number: 8"
    ## 292.347 sec elapsed
    ## [1] "Processing batch number: 9"
    ## 284.766 sec elapsed
    ## [1] "Processing batch number: 10"
    ## 275.864 sec elapsed

## Overview of the Fetched SNP Data

Here’s a glimpse of the fetched SNP data:

``` r
glimpse(snp_data)
```

    ## Rows: 10,000
    ## Columns: 16
    ## $ query            <chr> "rs2887286", "rs6685064", "rs2840528", "rs3890745", "…
    ## $ chromosome       <chr> "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1"…
    ## $ bp               <dbl> 1220751, 1275912, 2352457, 2622185, 3164291, 3765267,…
    ## $ class            <chr> "snv", "snv", "snv", "snv", "snv", "snv", "snv", "snv…
    ## $ rsid             <chr> "rs2887286", "rs6685064", "rs2840528", "rs3890745", "…
    ## $ gene             <chr> "SDF4", "", "MORN1/LOC100129534", "MMEL1", "PRDM16/LO…
    ## $ alleles          <chr> "T,C", "C,A,T", "A,G", "T,C", "A,C,G", "T,A,C", "G,A,…
    ## $ ancestral_allele <chr> "T", "C", "A", "T", "A", "T", "G", "C", "C", "C", "A"…
    ## $ variation_allele <chr> "C", "A,T", "G", "C", "C,G", "A,C", "A,C", "A,T", "A,…
    ## $ seqname          <chr> "NC_000001.11", "NC_000001.11", "NC_000001.11", "NC_0…
    ## $ hgvs             <chr> "NC_000001.11:g.1220751T>C", "NC_000001.11:g.1275912C…
    ## $ assembly         <chr> "GRCh38.p14", "GRCh38.p14", "GRCh38.p14", "GRCh38.p14…
    ## $ ref_seq          <chr> "T", NA, "A", "T", NA, "T", "G", NA, "C", "C", "A", "…
    ## $ minor            <chr> "C", NA, "G", "C", NA, "C", "A", NA, "T", "T", "G", "…
    ## $ maf              <dbl> 0.3008832, NA, 0.5262692, 0.4006754, NA, 0.1418150, 0…
    ## $ maf_population   <list> [<data.frame[25 x 4]>], [<data.frame[26 x 4]>], [<da…

## Data Cleaning

Remove unnecessary columns and rename certain columns for clarity.

``` r
# Mainly removing the list of the MAF per population (last column)
snp_data2 = snp_data %>% select(-query, -maf_population) %>% rename(position = bp)
```

    ## select: dropped 2 variables (query, maf_population)

    ## rename: renamed one variable (position)

``` r
glimpse(snp_data2)
```

    ## Rows: 10,000
    ## Columns: 14
    ## $ chromosome       <chr> "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1"…
    ## $ position         <dbl> 1220751, 1275912, 2352457, 2622185, 3164291, 3765267,…
    ## $ class            <chr> "snv", "snv", "snv", "snv", "snv", "snv", "snv", "snv…
    ## $ rsid             <chr> "rs2887286", "rs6685064", "rs2840528", "rs3890745", "…
    ## $ gene             <chr> "SDF4", "", "MORN1/LOC100129534", "MMEL1", "PRDM16/LO…
    ## $ alleles          <chr> "T,C", "C,A,T", "A,G", "T,C", "A,C,G", "T,A,C", "G,A,…
    ## $ ancestral_allele <chr> "T", "C", "A", "T", "A", "T", "G", "C", "C", "C", "A"…
    ## $ variation_allele <chr> "C", "A,T", "G", "C", "C,G", "A,C", "A,C", "A,T", "A,…
    ## $ seqname          <chr> "NC_000001.11", "NC_000001.11", "NC_000001.11", "NC_0…
    ## $ hgvs             <chr> "NC_000001.11:g.1220751T>C", "NC_000001.11:g.1275912C…
    ## $ assembly         <chr> "GRCh38.p14", "GRCh38.p14", "GRCh38.p14", "GRCh38.p14…
    ## $ ref_seq          <chr> "T", NA, "A", "T", NA, "T", "G", NA, "C", "C", "A", "…
    ## $ minor            <chr> "C", NA, "G", "C", NA, "C", "A", NA, "T", "T", "G", "…
    ## $ maf              <dbl> 0.3008832, NA, 0.5262692, 0.4006754, NA, 0.1418150, 0…

## Save Cleaned SNP Data to TSV

The cleaned SNP data is saved to a TSV file for further analysis or
storage.

``` r
write_tsv(snp_data2, "data/fingerprinting_snps.tsv")
```

## Convert SNP Data to BED Format

The SNP data is converted into BED format, which is commonly used in
genomics to represent genomic regions.

``` r
# Create BED format
bed_df = snp_data2 %>%
    transmute(
        chrom = chromosome,
        start = position - 1,
        end = position,
        rsid = rsid)
```

    ## transmute: dropped 13 variables (chromosome, position, class, gene, alleles, …)

    ##            new variable 'chrom' (character) with 22 unique values and 0% NA

    ##            new variable 'start' (double) with 10,000 unique values and 0% NA

    ##            new variable 'end' (double) with 10,000 unique values and 0% NA

``` r
glimpse(bed_df)
```

    ## Rows: 10,000
    ## Columns: 4
    ## $ chrom <chr> "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1",…
    ## $ start <dbl> 1220750, 1275911, 2352456, 2622184, 3164290, 3765266, 3826754, 4…
    ## $ end   <dbl> 1220751, 1275912, 2352457, 2622185, 3164291, 3765267, 3826755, 4…
    ## $ rsid  <chr> "rs2887286", "rs6685064", "rs2840528", "rs3890745", "rs1798246",…

## Save BED Formatted Data

The BED formatted data is saved to a file.

``` r
write_tsv(bed_df, "data/fingerprinting_snps.bed", col_names = FALSE)
```

## Filter for Bi-Allelic SNPs

Filter the SNP data to include only biallelic variants (SNPs with two
alleles).

``` r
snp_data3 = snp_data2 %>% filter(!is.na(minor))
```

    ## filter: removed 1,290 rows (13%), 8,710 rows remaining

``` r
glimpse(snp_data3)
```

    ## Rows: 8,710
    ## Columns: 14
    ## $ chromosome       <chr> "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1"…
    ## $ position         <dbl> 1220751, 2352457, 2622185, 3765267, 3826755, 4304166,…
    ## $ class            <chr> "snv", "snv", "snv", "snv", "snv", "snv", "snv", "snv…
    ## $ rsid             <chr> "rs2887286", "rs2840528", "rs3890745", "rs1181875", "…
    ## $ gene             <chr> "SDF4", "MORN1/LOC100129534", "MMEL1", "CCDC27", "CEP…
    ## $ alleles          <chr> "T,C", "A,G", "T,C", "T,A,C", "G,A,C", "C,A,G,T", "C,…
    ## $ ancestral_allele <chr> "T", "A", "T", "T", "G", "C", "C", "A", "T", "A", "T"…
    ## $ variation_allele <chr> "C", "G", "C", "A,C", "A,C", "A,G,T", "T", "G", "C,G"…
    ## $ seqname          <chr> "NC_000001.11", "NC_000001.11", "NC_000001.11", "NC_0…
    ## $ hgvs             <chr> "NC_000001.11:g.1220751T>C", "NC_000001.11:g.2352457A…
    ## $ assembly         <chr> "GRCh38.p14", "GRCh38.p14", "GRCh38.p14", "GRCh38.p14…
    ## $ ref_seq          <chr> "T", "A", "T", "T", "G", "C", "C", "A", "T", "A", "T"…
    ## $ minor            <chr> "C", "G", "C", "C", "A", "T", "T", "G", "G", "G", "C"…
    ## $ maf              <dbl> 0.3008832, 0.5262692, 0.4006754, 0.1418150, 0.3889897…

## Save Bi-Allelic SNP Data to TSV

The bi-allelic SNP data is saved to a TSV file.

``` r
write_tsv(snp_data3, "data/fingerprinting_snps_biallelic.tsv")
```

## Convert Bi-Allelic SNP Data to BED Format

The bi-allelic SNP data is converted into BED format.

``` r
bed_df2 = snp_data3 %>%
    transmute(
        chrom = chromosome,
        start = position - 1,
        end = position,
        rsid = rsid)
```

    ## transmute: dropped 13 variables (chromosome, position, class, gene, alleles, …)

    ##            new variable 'chrom' (character) with 22 unique values and 0% NA

    ##            new variable 'start' (double) with 8,710 unique values and 0% NA

    ##            new variable 'end' (double) with 8,710 unique values and 0% NA

``` r
glimpse(bed_df2)
```

    ## Rows: 8,710
    ## Columns: 4
    ## $ chrom <chr> "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1",…
    ## $ start <dbl> 1220750, 2352456, 2622184, 3765266, 3826754, 4304165, 4436598, 4…
    ## $ end   <dbl> 1220751, 2352457, 2622185, 3765267, 3826755, 4304166, 4436599, 4…
    ## $ rsid  <chr> "rs2887286", "rs2840528", "rs3890745", "rs1181875", "rs6663840",…

## Save Bi-Allelic BED Formatted Data

The BED formatted data for bi-allelic SNPs is saved to a file.

``` r
write_tsv(bed_df2, "data/fingerprinting_snps_biallelic.bed", col_names = FALSE)
```

## Session Info

``` r
sessionInfo()
```

    ## R version 4.3.1 (2023-06-16)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 22.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: Africa/Cairo
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] tictoc_1.2      tidylog_1.0.2   lubridate_1.9.2 forcats_1.0.0  
    ##  [5] stringr_1.5.0   dplyr_1.1.2     purrr_1.0.2     readr_2.1.4    
    ##  [9] tidyr_1.3.0     tibble_3.2.1    ggplot2_3.4.3   tidyverse_2.0.0
    ## [13] rsnps_0.6.1    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] utf8_1.2.3        generics_0.1.3    stringi_1.7.12    httpcode_0.3.0   
    ##  [5] hms_1.1.3         digest_0.6.33     magrittr_2.0.3    evaluate_0.21    
    ##  [9] grid_4.3.1        timechange_0.2.0  fastmap_1.1.1     jsonlite_1.8.7   
    ## [13] plyr_1.8.8        crul_1.4.0        httr_1.4.7        fansi_1.0.4      
    ## [17] scales_1.2.1      cli_3.6.1         crayon_1.5.2      rlang_1.1.1      
    ## [21] bit64_4.0.5       munsell_0.5.0     withr_2.5.0       yaml_2.3.7       
    ## [25] parallel_4.3.1    tools_4.3.1       tzdb_0.4.0        colorspace_2.1-0 
    ## [29] curl_5.0.2        vctrs_0.6.3       R6_2.5.1          lifecycle_1.0.3  
    ## [33] bit_4.0.5         vroom_1.6.3       clisymbols_1.2.0  pkgconfig_2.0.3  
    ## [37] pillar_1.9.0      gtable_0.3.4      glue_1.6.2        Rcpp_1.0.11      
    ## [41] xfun_0.40         tidyselect_1.2.0  rstudioapi_0.15.0 knitr_1.43       
    ## [45] htmltools_0.5.6   rmarkdown_2.24    compiler_4.3.1
