# setup
# clear the environment
rm(list=ls())

DATA_DIR <- './data'
IMAGES_DIR <- './images'
OUTPUT_DIR <- './output'

make_dir <- function(d) {
    if (file.exists(d)) unlink(d, recursive=TRUE, force=TRUE)
    dir.create(d)
}
lapply(c(IMAGES_DIR, OUTPUT_DIR),make_dir)


## function that concatenates strings (useful for directory paths)
concat <- function(x1,x2) {
    result <- paste(x1,x2,sep="")
    return(result)
}

## function that checks to see if a package is installed and,if not,installs it
## portions of this code came from http://stackoverflow.com/questions/9341635/how-can-i-check-for-installed-r-packages-before-running-install-packages
load_package <- function(x) {
    if (x %in% rownames(installed.packages())) { 
        print(concat("package already installed: ", x))
    }
    else { 
        install.packages(x, repos="http://mirror.las.iastate.edu/CRAN/") 
    }
    library(x, character.only=TRUE)
}

####################################################
# http://www.ats.ucla.edu/stat/data/discrim.sav

load_package("foreign")
URL <- "http://www.ats.ucla.edu/stat/data/discrim.sav"
data <- read.spss(concat(DATA_DIR,'/discrim.sav'), to.data.frame=TRUE)

# make factors and factor names
# data$JOB <- as.factor(data$JOB)
load_package('plyr')
data$JOB_NAME <- as.factor(data$JOB)
data$JOB_NAME <- revalue(data$JOB_NAME, c("1"="customer service", "2"="mechanic","3"="dispatcher"))
data$JOB_NAME <- as.factor(data$JOB_NAME)
str(data)
data <- data[,-5]
head(data)

# test for homogeneity of covariance matrices
log(det(x=cov(as.matrix(data[data[,4]=="1",1:3]))))
log(det(x=cov(as.matrix(data[data[,4]=="2",1:3]))))
log(det(x=cov(as.matrix(data[data[,4]=="3",1:3]))))
source("BoxMTest.R")
BoxMTest(data[,1:3],data[,5])

# perform cca
load_package("CCA")
X <- as.matrix(data[,1:3])
Y <- as.matrix(as.numeric(data[,4]))
cc1 <- cc(X, Y)
cc1$cor

# Step 1: compute d-dimesnional mean vectors; d is 3
load_package("DiscriMiner")
t(groupMeans(variables=data[,1:3], group=data[,5]))
t(groupStds(variables=data[,1:3], group=data[,5]))
# test correlations
cor.test(data[,1],data[,2])
cor.test(data[,1],data[,3])
cor.test(data[,2],data[,3])

# Step 2: compute d-dimensional scatter matrices
betweenCov(variables=data[,1:3], group=data[,4], div_by_n = FALSE)
withinCov(variables=data[,1:3], group=data[,4])
(Sb <- betweenSS(variables=data[,1:3], group=data[,4]))
(Sw <- withinSS(variables=data[,1:3], group=data[,4]))

# Step 3: calculate eigenvalues and eigenvectors
(eigenvalues <- eigen(solve(Sw) %*% Sb)$values)


# Step 4: Select eigenvalues
plot(eigenvalues, type="b")

# Step 5: transform samples onto new subspace provided by LDA
# linear discriminant analysis
load_package("MASS")
scaled_data <- data
scaled_data[,1:3] <- scale(data[,1:3],
                           center=TRUE, 
                           scale=apply(data[,1:3],2, sd, na.rm=TRUE))
lda <- lda(scaled_data[,1:3], scaled_data[,5], prior=c(1,1,1)/3)
lda$scaling
plda <- predict(lda)

ldahist(data=plda$x[,1], g=data[,5])
ldahist(data=plda$x[,2], g=data[,5])

load_package('dplyr')
grouped_by_job <- group_by(data[,1:3],data[,5])
summarise(grouped_by_job, n=n())

# confusion matrix
(confusion_matrix <- table(plda$class, data[,5]))
round(diag(prop.table(confusion_matrix)), 4)
round(sum(diag(prop.table(confusion_matrix))), 4)
# cross-validated
lda2 <- lda(scaled_data[,5] ~ ., data=scaled_data[,1:3], CV=TRUE)
# confusion matrix
(confusion_matrix <- table(lda2$class, data[,5]))
round(sum(diag(prop.table(confusion_matrix))), 4)
# Calculate Wilk's Lambda
# It also gives us the group centroids
load_package('rrcov')
Wilks.test(plda$x,data[,5])


load_package("ggplot2")
load_package("scales")

class <- scaled_data[,5]
x <- plda$x[,1]
y <- plda$x[,2]
df <- data.frame(class, x, y)
centroids <- aggregate(cbind(x,y)~class,df,mean)
prop.lda = lda$svd^2/sum(lda$svd^2)
ggplot(df,aes(x,y,color=factor(class),shape=factor(class))) +
    geom_point(size=2.5) + 
    geom_point(data=centroids,size=5) +
    labs(x = paste("LD1 (", percent(prop.lda[1]), ")", sep=""),
         y = paste("LD2 (", percent(prop.lda[2]), ")", sep="")) +
    ggtitle("LDA Projection Data") + 
    theme(plot.title = element_text(lineheight=.8, face="bold"))
