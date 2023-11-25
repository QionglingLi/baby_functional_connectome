# install.packages("nlme")
# install.packages("mgcv")
# install.packages("parallel")
# install.packages("foreach")
# install.packages("doParallel")
# install.packages("tcltk")
# install.packages("doSNOW")

library(nlme)
library(mgcv)
# library(parallel)
# library(foreach)
# library(doParallel)
# library(tcltk)
# library(doSNOW)

args <- commandArgs(TRUE)
voxel <- args[1]
voxel <- as.numeric(voxel)

# load data
if (voxel==1) {
Measure <- read.csv("/HeLabData2/qlli/project2/NodalMeasures/FCs/FCs.csv", nrows=1, header = FALSE) 
} else {
Measure <- read.csv("/HeLabData2/qlli/project2/NodalMeasures/FCs/FCs.csv", nrows=1, skip=voxel-1, header = FALSE)
}

Data <- read.csv("/HeLabData2/qlli/project2/NodalMeasures/FCs/basicparas930.csv")
Data$sex <- as.factor(Data$sex)
Data$site <- as.factor(Data$site)
Data$sub_id <- as.factor(Data$sub_id)
scan_age=Data$scan_age
sex=Data$sex
site=Data$site
mFD=Data$mFD
sub_id=Data$sub_id

# set up gam model
Data$Measure <- t(Measure)
mod_gam <- gam(Measure ~ s(scan_age, bs="cs", k=3) + sex + s(sub_id, bs="re") + site + mFD, 
               data=Data, method="REML", na.action="na.omit")

## prediction with real age              
# age_term = predict(mod_gam, type = 'terms')[, 's(scan_age)']
# res_partial = residuals(mod_gam) + age_term
# predvalue <- t(c(voxel,res_partial))

# prediction with generative age
Length=ncol(Measure)
G1_pred <- with(Data, 
                expand.grid(scan_age=seq(min(scan_age), max(scan_age), length = Length)))
G1_pred <-cbind(G1_pred,
                predict(mod_gam,
                        newdata=G1_pred,
                        se.fit = TRUE,
                        type = "link",
                        exclude = c("sex", "s(sub_id)", "mFD", "site")))
age_term=G1_pred$`fit`
age_se=G1_pred$`se.fit`
predvalue <- t(c(voxel,age_term))
predse <- t(c(voxel,age_se))

# pick p&f values for age term
summary_model <- summary(mod_gam)
p_value <- summary_model[["s.table"]]["s(scan_age)","p-value"]
F_value <- summary_model[["s.table"]]["s(scan_age)","F"]
pvalues <- t(c(voxel,p_value))
Fvalues <- t(c(voxel,F_value))

# write p&f values with their corresponding voxel index
write.table( pvalues, file = paste0("/HeLabData2/qlli/project2/NodalMeasures/FCs/outvalues/FCs_p_value_", voxel, ".txt", sep="" ), row.names = FALSE, col.names = FALSE, quote = FALSE )
write.table( Fvalues, file = paste0("/HeLabData2/qlli/project2/NodalMeasures/FCs/outvalues/FCs_F_value_", voxel, ".txt", sep="" ), row.names = FALSE, col.names = FALSE, quote = FALSE )
write.table( predvalue, file = paste0("/HeLabData2/qlli/project2/NodalMeasures/FCs/outvalues/predvalue_", voxel, ".txt", sep="" ), row.names = FALSE, col.names = FALSE, quote = FALSE )
write.table( predse, file = paste0("/HeLabData2/qlli/project2/NodalMeasures/FCs/outvalues/predse_", voxel, ".txt", sep="" ), row.names = FALSE, col.names = FALSE, quote = FALSE )
#write.table( max_age, file = paste0("/HeLabData2/qlli/project2/NodalMeasures/FCs/outvalues/maxage_", voxel, ".txt", sep="" ), row.names = FALSE, col.names = FALSE, quote = FALSE )
#write.table( min_age, file = paste0("/HeLabData2/qlli/project2/NodalMeasures/FCs/outvalues/minage_", voxel, ".txt", sep="" ), row.names = FALSE, col.names = FALSE, quote = FALSE )


