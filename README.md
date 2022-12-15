# Bioinformatics in Agricultural Crops

#### Abstract
Agriculture is deeply affected by climate change in terms of yield and quality of crops. Bioinformatics has been contributing to the sustainable development of agriculture; it has helped understand the plant’s genetic resources and its response to numerous factors. Due to the complexity of stress conditions, there is a crucial need to discover the mechanisms of genetic response against various stress factors and study the patterns of the stress-responsive genome molecular mechanisms. Keeping these responses in mind, crops that can stand stronger in harsh conditions could be reproduced with bioengineering. Previous works include microarray analysis to determine the cross-responsive genes in plants under influence of biotic and abiotic factors. This study is an effort to perform multiclass classification on the expressed genes. Random Forest algorithm is used to classify differentially expressed genes (DEGs) in cotton influenced by abiotic factors like cold, acid, salinity, and drought using the GEO (Series GSE50770) dataset. The model performance was outlined by a confusion matrix. Future work will be focused on the efficient filtering of the least responsive genes, improving the accuracy of the model, identifying the differentially expressed genes, and deriving relevant conclusions that can contribute to the survival of agricultural crops on the planet.

#
Steps for Windows OS:
1) Download and install R version 4.2.2 or latest from [here](https://cran.r-project.org/bin/windows/base/)
2) Download and install RStudio from [here](https://posit.co/download/rstudio-desktop/)
3) Download and open code.R in RStudio
4) If prompted, install the necessary packages with respect to the libraries mentioned in first few lines of code.R. If not prompted, use the following command to install the packages:
 ``` bash
install.packages('nameofpackage') 
 ```
5) Run code.R 
6) Dataset will be automatically fetched from Gene Expression Omnibus [website](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50770) 
7) Plots can be viewed on bottom right corner of RStudio (Note: Keep adequate margin or space of the plot section to view plots correctly)
