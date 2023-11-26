library(nlme)
library(mgcv)
library(reshape2)
library(RColorBrewer)
library(viridis)
library(ggplot2)

outdir="~/pca/sub_PCA"
Measures <- read.csv("~/pc_bin10.csv", header = FALSE)
Data <- read.csv("~/basicparas930.csv")

Data$sub_id <- as.factor(Data$sub_id)
sub_id=Data$sub_id
scan_age=Data$scan_age
Data$sex <- as.factor(Data$sex)
Data$site <- as.factor(Data$site)
sex=Data$sex
site=Data$site
mFD=Data$mFD
Len=nrow(Data)

cols <- magma(10)[matrix(1:10)]
names(cols) <- c(paste("bin", matrix(1:10), sep=""))

p_values <- c()
F_values <- c()
yy<- c()
for (voxel in 1:10) {
  print(voxel)
  Data$Measure <- as.numeric(t(Measures[voxel,]))
  Measure = Data$Measure
  y=Measure
  
  mod_gam <- gam(y ~ s(scan_age, bs="cs", k=3) + sex + s(sub_id, bs="re") + site + mFD, 
                 data=Data, method="REML", na.action="na.omit")
  summary_model <- summary(mod_gam)
  
  p_values = cbind(p_values,summary_model[["s.table"]]["s(scan_age)","p-value"])
  F_values = cbind(F_values,summary_model[["s.table"]]["s(scan_age)","F"])
  
  
  G1_pred <- with(Data, 
                  expand.grid(scan_age=seq(min(scan_age), max(scan_age), length = Len)))
  G1_pred <-cbind(G1_pred,
                  predict(mod_gam,
                          newdata=G1_pred,
                          se.fit = TRUE,
                          type = "link",
                          exclude = c("sex", "s(sub_id)", "mFD", "site")))
  age_term = G1_pred$`fit`
  yy = cbind(yy,age_term)
  
}
write.table( p_values, file = paste0( outdir,"/seed_p_values.txt", sep="" ), row.names = FALSE, col.names = FALSE, quote = FALSE )
write.table( F_values, file = paste0( outdir,"/seed_F_values.txt", sep="" ), row.names = FALSE, col.names = FALSE, quote = FALSE )


#plot trajectories for all bins
scan_age=seq(min(scan_age), max(scan_age), length = Len)
myplot <- ggplot() +
  geom_line(aes(y=yy[,1], x = scan_age), colour=cols[1]) + 
  labs(y = expression(FCs), x = expression(scan_age))
theme_classic()


for (voxel in 2:10) {
  myplot <- myplot+
    geom_line(aes_(y=yy[,voxel], x = scan_age[order(scan_age)]), colour=cols[voxel]) + 
    labs(y = expression(FCs), x = expression(scan_age))
  theme_classic()
}
myplot <- myplot+
  scale_colour_manual("",values=cols, breaks=names(cols))
