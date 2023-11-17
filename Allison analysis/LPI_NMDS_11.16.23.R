## Running an NMDS on the bat data! 

library(tidyverse)
library(dplyr)
library(reshape)

##############################################################
#here reading in the data (using RAI x10s and rounded to develop 
#a site by species matrix accounting for effort.... [not accounting 
#for detection differences though...])

library(reshape)

############################################################
#installing data and converting to a site by species matrix 
#reading in mm_cov2
mm_cov2 <- read.csv("mm_cov2.csv")

#and cov1
mm_cov1 <- read.csv("mm_cov1.csv")

#now subsetting for the data we want (using 100r because some ie mink
# would be all 0s)
colnames(mm_cov1)
mm_nmds <- mm_cov1 %>% dplyr::select(Point, MarshRabbit_RAI100r, deer_RAI100r, 
                                GraySquirrel_RAI100r, Bear_RAI100r, 
                                Bobcat_RAI100r, Hog_RAI100r, Panther_RAI100r,
                                Opossum_RAI100r, SpottedSkunk_RAI100r, 
                                Raccoon_RAI100r, Rodent_RAI100r, Mink_RAI100r)
#already is a site by species matrix 

##### Running NMDS on data #####
#loading necessary packages 
library(raster)
library(cluster)
library(vegan)
library(ade4)
library(ggplot2)

####running the PerMANOVA####
#this is a nonparametric method for testing 
#the hypothesis of no differences between 2 
#or more groups based on rank dissimilarities 
#often paired with NMDS

#transforming data --> wisconsin transformation
#then bray-curtis, then NMDS
#going to not transform and then transform and see if its different

#first going to run a NMDS
#first making sure the data is numeric 
#subsetting for only the mammal RAIs
only_mm <- mm_cov1 %>% dplyr::select(MarshRabbit_RAI100r, deer_RAI100r, 
                                     GraySquirrel_RAI100r, Bear_RAI100r, 
                                     Bobcat_RAI100r, Hog_RAI100r, Panther_RAI100r,
                                     Opossum_RAI100r, SpottedSkunk_RAI100r, 
                                     Raccoon_RAI100r, Rodent_RAI100r, Mink_RAI100r)

#vs rai10
only_mm1 <- mm_cov1 %>% dplyr::select(MarshRabbit_RAI10r, deer_RAI10r, 
                                     GraySquirrel_RAI10r, Bear_RAI10r, 
                                     Bobcat_RAI10r, Hog_RAI10r, Panther_RAI10r,
                                     Opossum_RAI10r, SpottedSkunk_RAI10r, 
                                     Raccoon_RAI10r, Rodent_RAI10r, Mink_RAI10r)



#then running a dist matrix
jmm <-vegdist(only_mm1, "bray") 
jmm

#??metaMDS
nmds_mammal<-metaMDS(jmm,k=2, trace=T) #k =2 bc only interested
#in 2 dimensions

#looking at the stress plot
stressplot(nmds_mammal) 
#nonmetric R2 is 0.986 and linear fit R2 is 0.967

#and what the actual stress is
nmds_mammal
# stress is 0.1193787 which is "fair" I think...
# this is showing how well this represents the original distance matrix 

####NMDS Plotting with ggplot2####

#ggplot
data.scores <- as.data.frame(scores(nmds_mammal))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$rest <- mm_cov1$updated_rest_status  # create a column of rest categories, from the data frame "rest_subset"
data.scores$site <- mm_cov1$Point  #  create a column for site names

head(data.scores)  #look at the data

#plotting the NMDS
ggplot(data.scores, aes(x= NMDS1, y= NMDS2, col=rest, shape=rest)) + #denote groupings by color "col" and shape
  geom_point() +#adds points
  geom_text(aes(label=mm_cov1$Point),hjust=0, vjust=0)+#adds site names
  stat_ellipse() +#adds ellipses
  theme_bw() +
  xlim(-1.2, 1.2)+
  ylim(-.75,.75)+
  labs(title = "NMDS Plot")
#Here we can see the NMDS visualized

#looks pretty much exactly the same between x100 RAI and x10 RAI

###################################################################3
##############also going to look at it by veg category 
###################################################################3

#ggplot
data.scores <- as.data.frame(scores(nmds_mammal))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$veg <- mm_cov1$sum_veg_cat2  # create a column of rest categories, from the data frame "rest_subset"
data.scores$site <- mm_cov1$Point  #  create a column for site names

head(data.scores)  #look at the data

#plotting the NMDS
ggplot(data.scores, aes(x= NMDS1, y= NMDS2, col=veg, shape=veg)) + #denote groupings by color "col" and shape
  geom_point() +#adds points
  geom_text(aes(label=mm_cov1$Point),hjust=0, vjust=0)+#adds site names
  stat_ellipse() +#adds ellipses
  theme_bw() +
  xlim(-1.1, 1.1)+
  ylim(-.75,.75)+
  labs(title = "NMDS Plot")

#### Wisconsin double transformation NMDS ####
#first going to run a NMDS
tmm <- wisconsin(only_mm1)

tjmm<-vegdist(tmm, "bray") 
tjmm

#??metaMDS
tnmds_mm<-metaMDS(tjmm,k=2, trace=T) #k =2 bc only interested
#in 2 dimensions

#looking at the stress plot
stressplot(tnmds_mm)
#and what the actual stress it
tnmds_mm
# stress is .16 which is worse than untransformed

##NMDS Plotting with ggplot2#

#ggplot
tdata.scores <- as.data.frame(scores(tnmdsBat))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
tdata.scores$rest <- bats2$rest_status  # create a column of rest categories, from the data frame "rest_subset"
tdata.scores$site <- bats2$sites  #  create a column for site names

head(tdata.scores)  #look at the data

#plotting the NMDS
ggplot(tdata.scores, aes(x= NMDS1, y= NMDS2, col=rest, shape=rest)) + #denote groupings by color "col" and shape
  geom_point() +#adds points
  geom_text(aes(label=bats2$site),hjust=0, vjust=0)+#adds site names
  stat_ellipse() +#adds ellipses
  theme_bw() +
  xlim(-.75, .75)+
  ylim(-.75,.75)+
  labs(title = "NMDS Plot")
#Here we can see the NMDS visualized
#does not look very different from the untransformed data
#going to move on with the untransformed data

#####Next going to use the adonis function in the vegan package####
#to run the PerMANOVA
#?adonis

set.seed(11)

#making the the restoration subset a matrix
rest_subset <- mm_cov1$updated_rest_status
rest_subset_m <- as.matrix(rest_subset)
str(rest_subset_m)

#running the adonis function
perm_mm <- adonis2(jmm ~ rest_subset_m, permutations=1000)
perm_mm

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 1000
# 
# adonis2(formula = jmm ~ rest_subset_m, permutations = 1000)
#                 Df SumOfSqs      R2      F Pr(>F)
# rest_subset_m  3   0.4289 0.04634 1.2959   0.2348
# Residual      80   8.8260 0.95366              
# Total         83   9.2549 1.00000  

#not significant differences in community among the groups

#looking @ the output table and plot the permuted f-ratios 
hist(perm_mm$f.perms, main = "Histogram of F statistics for Bats" ,xlim=c(0,12))
points(perm_mm$aov.tab$F[1],0, pch=19,col="red", bg="red", cex=2)
#groups appear very different! 

##hmm this part isn't running for some reason, 
# although doesn't really matter anyways since no sig difs. 




#####betadisps function to understand homogeneity of multivariate####
dis <- jbats
groups <- rest_subset_m

#betadisp --> on distance matrix 
MVdisp <- betadisper(dis, groups)
MVdisp

#to look at a boxplot of the distances to the centroid 
boxplot(MVdisp)

#Average distance to median:
#partially_restored      reference_Fak           restored         unrestored 
#0.2683                      0.3109               0.2197             0.2440 
#this helps me understand which are more clustered/less clustered (variance-wise)
#so we can see that the restored is the most clustered, then 
#unrestored, then partially restored and reference has the most variance (most dispersed)..

#F-test
anova(MVdisp)
#F value = 1.7588
#P value = 0.1624 <- not significant, so variances of 
#groups are homogeneous
#which fits one of the assumptions of PerMANOVA (equal variance
#between groups)

#Permutation test for F
pMVdisp <- permutest(MVdisp, permutations = 99, pairwise = TRUE)
pMVdisp
#looking at pairwise homogeneity of variance
#only significant difference in variance pairwise is between
#restored and reference p-value = (0.04)
#but not considered significant..

##simper - shows you which species are responsible for the 
#difference between groups that you observe
?simper
sim <- simper(bats3, groups)
sim

#summary of important species
summary(sim)
#TADBRA <- (brazilian free tailed bats) 
#most abundant species, so not surprising that 
#this sepcies is having the largest impact on 
#differences between groups 


#### Running CART analysis on Veg data####
#in order to understand what vegetation variables 
#are most strongly associated with the restoration categories

#In this case the restoration categories are the response
#variable, and the rest of the veg categories are the explanatory
#variables

#reading in packages 
library(MASS)
library(rpart)
library(ade4)
library(vegan)
library(rattle) 
library(dplyr)

#using veg dataset --> veg_subF

#first going to create a training dataset 
#with half of the dataset - 39 out of 78
set.seed(5)
train_sub <- sample(1:78, 39) 

#here checking the frequency of each rest category in the training data
#to make sure they are relatively proportional to the frequency 
#of the complete data set 
freq<-table(veg_subF$rest_status[train_sub]) # <- change to whatever I 
#save as the master veg file
freq

#next going to specify the CART model
#in my vegetation dataset, restoration category is the categorical response 
#variables and the other measures of veg are the predictor
#variables 
model <- rest_status ~ .

##Running the CART algorithm 
#I want to use a Classification tree because my response variables are 
#categorical (restoration status)
#you will use the rpart function --> here using the subsetted veg_class
library(rpart)
#library(rattle)
veg_rpart_cat <- rpart(model, data = veg_subF[train_sub,], method="class", 
                       control = rpart.control(minsplit = 10))

###Plotting the CART tree and viewing summary 
#going to plot using the function post.rpart
#?post.rpart
#post(veg_rpart_cat, file = "", title = " Veg Classification Tree")
#woww that is ugly...

#plotting fancy plot 
fancyRpartPlot(veg_rpart_cat, sub = "" )  
#?fancyRpartPlot

#Now looking at the node by node summary of the tree
#and the variable importance 
summary(veg_rpart_cat)
#Variable importance!! 
#most important variable is avg cc! 

#?rpart
###Cost-complexity pruning##
#looking at the cross-validation results
printcp(veg_rpart_cat)
#root node error; 23/39 = 0.589

#now going to plot the results to determine the optimal tree
plotcp(veg_rpart_cat)
#from this is looks like canopy coverage (my first node definer) is my most important, and 
#maybe only important variable...

### Pruning the tree ###
#going to prune the the tree using prune.rpart function

#first going to set the cp for the optimal tree size
#this extracts the cp value according to the 1-SE rule 
cp<-printcp(veg_rpart_cat)[2,1]

#and here we are pruning 
veg_prune <- prune(veg_rpart_cat, cp = cp)
print(veg_prune)

#and here plotting the pruned tree 
post(veg_prune, file = "", title = " Veg pruned tree")
fancyRpartPlot(veg_prune, sub = "")
#play with me more

##Classification accuracy##
#now we;re going to compare how accurately both the unpruned and 
#pruned tree classify the testing data

#classification matrix
Ct_unprune<-table(predict(veg_rpart_cat, veg_subF[-train_sub,], type = "class"), veg_subF[-train_sub, "rest_status"])
Ct_prune<-table(predict(veg_prune, veg_subF[-train_sub,], type = "class"), veg_subF[-train_sub, "rest_status"])

#classification accuracy 
class_unprune<-sum(diag(prop.table(Ct_unprune)))
class_prune<-sum(diag(prop.table(Ct_prune)))

class_unprune
#0.487 <- only ~49% accuracy when applied to testing data....not great

class_prune
#0.487 the same!, no improvement with accuracy with pruning..


