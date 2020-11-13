TSAPA R package
====================

Identification of Tissue-Specific Alternative Polyadenylation sites in plants

About
====================
TSAPA is an open-source R package to identify the tissue-specific alternative polyadenylation sites in plants. This package can identify tissue-specific and constitutive poly(A) sites (tsPAs and csPAs) using the information entropy method. A generic model for rice japonica has been trained and integrated in TSAPA, which can be directly applied to predict tsPAs from any new data of rice japonica without additional inputs. TSAPA also allows users to customize the classification model by selecting desirable features or classifiers for model training on their own data. Given a list of poly(A) sites and the corresponding reference genome sequence, TSAPA is capable of predicting tsPAs for species other than rice japonica. TSAPA can also be used to identify and predict tsPAs in intergenic regions.

Installing TSAPA
=============
Mandatory 
---------

* R (>3.1). [R 3.3.3](https://www.r-project.org/) is recommended.

Required R Packages
---------
* [e1071](https://CRAN.R-project.org/package=e1071), [TCC](http://www.bioconductor.org/packages/release/bioc/html/TCC.html), [seqinr](https://CRAN.R-project.org/package=seqinr ), [stringr](https://CRAN.R-project.org/package=stringr), [Biostrings](http://www.bioconductor.org/packages/release/bioc/html/Biostrings.html), [NuPoP](http://master.bioconductor.org/packages/release/bioc/html/NuPoP.html), [adabag](https://CRAN.R-project.org/package=adabag), [randomForest](https://CRAN.R-project.org/package=randomForest), [Boruta](https://CRAN.R-project.org/package=Boruta), [DMwR](https://CRAN.R-project.org/package=DMwR), [GenomicFeatures](http://www.bioconductor.org/packages/release/bioc/html/GenomicFeatures.html)

Installation
---------
* Install the R package using the following commands on the R console:
```
install.packages("devtools")
library(devtools)
install_github("BMILAB/TSAPA")
library(TSAPA)
```

Using TSAPA
=============
To get started, the user is recommended to use the example dataset which comes with the packages.

Section 1 
---------
* Users can use the “row_data” to calculate the information entropy (H) and adjusted information entropy (modeH).
```
load(system.file("data","row_data.Rdata",package = "TSAPA"))
data <- roku_modeH(row_data)
load(system.file("data","sample_data.Rdata",package = "TSAPA"))
pa_data <- select_tsPA(data,1,1,0.8)
```

Section 2
---------
* Use the annotations file(rice_gff7_japonica.gff3), Users can obtain the standard annotations file.
```
file <- system.file("extdata","rice_gff7_japonica.gff3",package = "TSAPA")
chrLenFile <- system.file("extdata","rice_gff7.chrlen.csv",package = "TSAPA")
gff <- parseGFF(file=file,format='gff3',chrLenFile=chrLenFile,ofilePre=NULL)
```
* Users can map the information of ploy(A) site to the standard annotation file.
```
file <- system.file("extdata","rice_gff7_japonica.gff3",package = "TSAPA")
chrLenFile <- system.file("extdata","rice_gff7.chrlen.csv",package = "TSAPA")
gff <- parseGFF(file=file,format='gff3',chrLenFile=chrLenFile,ofilePre=NULL)
gff <- gff$ftrs
pac <- gff[1:10000,c('chr','strand','ftr_start')]
colnames(pac) <- c('chr','strand','coord')
pac$leaf <- 'leaf'
pac$root <- 'root'
pac2 <- mapPA2GFF(pac,gff)
```
* With the annotation file and the information table of ploy(A) sites, to obtain the APA feature.
```
path <- system.file("extdata","annotation.csv",package = "TSAPA")
annotation <- read.table(path,sep=",",header=T,stringsAsFactors = FALSE)
path1 <- system.file("extdata","tsPA_coordinate.csv",package = "TSAPA")
tsPA <- read.table(path1,sep=",",header=T,stringsAsFactors = FALSE)
apa_feature <- apa_context(annotation,tsPA)
```
* If users provide the annotation file of the genome, then the annotation information of ploy(A) sites can be obtained.
```
genome <- read.table(system.file("extdata","japonica.csv",package = "TSAPA"),sep=",",header=T,stringsAsFactors = FALSE)
pagene <- read.table(system.file("extdata","annotation.csv",package = "TSAPA"),sep=",",header=T,stringsAsFactors = FALSE)
gene_info <- geneinfo(genome,pagene,12429)
```
* Given a genome sequence file and a list of poly(A) sites, user can extracts the sequences of ploy(A) sites.
```
rice <- readDNAStringSet(system.file("extdata","rice_tsetseq.fasta",package = "TSAPA"))
tsPA <- read.table(system.file("extdata","extract_seq_PA.csv",package = "TSAPA"),sep=",",header=T,stringsAsFactors = FALSE)
extract_seq(rice,tsPA,"tsPA.fasta",-100,100)
```
* With the same sequence file, the features of Z-curve, nucleosome occupancy scores, K-grams can be calculated. 
```
zcurve_feature <- Zcurve(system.file("extdata","sample.fasta",package = "TSAPA"))
nupop_feature<-nupop(system.file("extdata","sample.fasta",package = "TSAPA"),100,150,151,200)
gram <- gram(system.file("extdata","sample.fasta",package = "TSAPA"),1)
```
* Using the sequence file and the motif file, the FHMM_feature and pwm_feature can be obtained.
```
FHMM_model <- FHMM(system.file("extdata","6-motif.fasta",package = "TSAPA"))
fhmm_feature <- FHMM_feature(6,system.file("extdata","sample.fasta",package = "TSAPA"),FHMM_model)
pwm_model <- pwm(system.file("extdata","6-motif.fasta",package = "TSAPA"))
pwmfeature <- pwm_feature(6,system.file("extdata","sample.fasta",package = "TSAPA"),pwm_model)
```
* With the file of feature space, users can use feature selection algorithm.
```
load(system.file("data","feature_selection_data.Rdata",package = "TSAPA"))
feature <- feature_selection(feature_selection_data,2)
```
* Users can test the accuracy of individual feature with Adaboost classification model and ranks the individual feature for the predict accuracy.
```
load(system.file("data","features_specific.Rdata",package = "TSAPA"))
load(system.file("data","features_unspecific.Rdata",package = "TSAPA"))
name <- c('Kgram','APAcontext','Conservation_score','FHMM','Zcurve','Nucleosome_Positioning','PWM','Secondary_structure')
featurerank <- feature_rank(features_specific,features_unspecific,name,normalized = T)
```
* To solve the unbalanced classification problem.
```
load(system.file("data","unbalanced_train_data.Rdata",package = "TSAPA"))
newdata<-balanced_data(unbalanced_train_data.Rdata,200,200)
```

Section 3
---------
* For training dataset and testing dataset, users can test the performance of classification model.
```
load(system.file("data","train.Rdata",package = "TSAPA"))
load(system.file("data","test.Rdata",package = "TSAPA"))
modelrank <- model_rank(train,test,normalized = F)
```
* Providing the feature file, users can use the real dataset to customize the classification model. 
```
load(system.file("data","features_specific.Rdata",package = "TSAPA"))
load(system.file("data","features_unspecific.Rdata",package = "TSAPA"))
data<-data.frame(rbind(features_specific,features_unspecific),label=c(rep("+1",nrow(features_specific)),rep("-1",nrow(features_unspecific))))
tem <- sample(2,nrow(data),replace=TRUE,prob=c(0.75,0.25))
train <- data[tem==1,]
test<-data[tem==2,]
model <- classification_model(train,2)
predict <- predict(model,test)
class <- predict$class
reslut <- predict_result(class,test)
reslut
```

Section 4
---------
*  Users can use the generic model to predict the test dataset. 
```
load(system.file("data","generic_model.Rdata",package = "TSAPA"))
load(system.file("data","test_data.Rdata",package = "TSAPA"))
predict1 <- Generic_model(generic_model,test_data)
class <- predict1$class
reslut <- predict_result(class,test_data)
reslut
```

Section 5
---------
*  Users can also predict tsPAs in intergenic regions. 
```
load(system.file("data","intergenic_model.Rdata",package = "TSAPA"))
load(system.file("data","intergenic_test.Rdata",package = "TSAPA"))
predict2<-Intergenic_model(intergenic_model,intergenic_test)
class2 <- predict2$class
reslut <- predict_result(class2,intergenic_test)
reslut
```

Citation
---------
If you are using TSAPA, please cite: [Ji G, Chen M, Ye W, Zhu S, Ye C, Su Y, Peng H and Wu X (2018) TSAPA: identification of tissue-specific alternative polyadenylation sites in plants, Bioinformatics, 34, 2123-2125.](https://academic.oup.com/bioinformatics/article/34/12/2123/4827682)
