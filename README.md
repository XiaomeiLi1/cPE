# cPE
## Identifying preeclampsia-associated genes using a control theory method
### Xiaomei Li<sup>1</sup>, Lin Liu<sup>1</sup>, Clare Whitehead<sup>2</sup>, Jiuyong Li<sup>1</sup>, Benjamin Thierry<sup>3</sup>,Thuc D. Le<sup>1</sup>, and Marnie Winter<sup>3</sup>

1. UniSA STEM, University of South Australia, Mawson Lakes, SA 5095, Australia
2. Pregnancy Research Centre, Dept of Obstetrics & Gynaecology, University of Melbourne, Royal Women’s Hospital, Melbourne, VIC
3052, Australia
3. ARC Centre of Excellence in Convergent BioNano Science and Technology and Future Industries Institute, University of South
Australia, Adelaide, SA 5095, Australia

Affecting 10 million women and causing 500,000 fetal and 76,000 maternal deaths globally each year, preeclampsia is a serious pregnancy complication which affects 5-8% of all pregnancies. While it is now recognized that preeclampsia causes serious short- and long-term complications for both the mother and their offspring, treatments to slow, stop or prevent preeclampsia development are lacking. Moreover, it is not possible to predict with a high accuracy which women will develop the condition in early pregnancy. Importantly, the triggering insults and precise molecular mechanisms which lead to the heterogenous condition are poorly understood which has hampered development of new treatments, biomarker discovery and early diagnosis. 
While the collection of preeclampsia relevant omic datasets is rapidly increasing and will undoubtedly yield significant insights. Reports of advanced computational approaches necessary to untangle the complexities of preeclampsia are sparse. In this manuscript, we describe for the first time a control theory method to identify preeclampsia associated genes. Firstly, a preeclampsia gene regulatory network was established, and its controllability analyzed. Next, two types of critical preeclampsia associated genes were identified which play important roles in the constructed preeclampsia specific network. Benchmarking against differential expression and hub analysis we demonstrated that the proposed method offers novel insights. Importantly, this control theory approach could aid in a further understanding of preeclampsia molecular mechanisms, identify biomarkers and ultimately contribute to improved preeclampsia diagnosis and treatment.

This repository includes the scripts and data of the proposed method. 

The scripts are in the R folder and include:

- PE.R - Script for identify critical preeclampsia associated genes
- PE_subtype.R - Script for identify critical early onset or late onset preeclampsia associated genes
- ProposedMethod_Functions.R - Script for other support functions

The input data are in the data folder and include:

- GSE75010.rda - preeclampsia gene expression data
- GSE75010_pd.rda - preeclampsia clinical information
- PESNPdb.csv - SNP database from http://bejerano.stanford.edu/pesnpdb
- PPI.xls - human directed PPI network from https://www.flyrnai.org/DirectedPPI/
- Browse Transcription Factors hg19 - resource_browser.csv - transcription factor (TF) list from the FANTOM5 transcriptome catalog database
- gold_standard.rda - a literature-based database (dbPEC) for preeclampsia
- proteinAffectingMut.rda - genes with preeclampsia related SNPs

The output data are in the result folder and include:
 - EoPE - results for early onset preeclampsia
 - LoPE - results for late onset preeclampsia
 - PE - results for preeclampsia
 - subtypes.xlsx - preeclampsia subtype specific critical genes
 
 Notes:
 
 The source code of controllability analysis can be downloaded from https://scholar.harvard.edu/yyl/code. Please compile it and put the libraray in a folder.


