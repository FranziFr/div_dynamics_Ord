library(reshape2) 
library(reshape)
library(plyr)
library(dplyr)
library(tidyr)
library(RMark)
library(ggplot2)
library(ggrepel)

## set working directory to load the PBDB file
setwd("M:/Franeck_and_Liow_2018")

## read in csv download from PBDB
genus <- read.csv("PBDB_Ord_1.csv", sep = ",", header=T)

## filter all genus data
genus <- filter(genus, grepl("genus", accepted_rank))


## model specifications for RMARK - Pradel seniority model
Phi.time <- list(formula=~time)
Gamma.time <- list(formula=~time)
p.time <- list(formula=~time)
L.time <- list(formula=~time)
Phi.time.plate <- list(formula=~time+geoplate)
p.time.plate <- list(formula=~time+geoplate)
Gamma.time.plate <- list(formula=~time+geoplate)
L.time.plate <- list(formula=~time+geoplate)

Phi.Time.plate <- list(formula=~time+geoplate+time*geoplate)
p.Time.plate <- list(formula=~time+geoplate+time*geoplate)
Gamma.Time.plate <- list(formula=~time+geoplate+time*geoplate)
L.Time.plate <- list(formula=~time+geoplate+time*geoplate)

Phi.const <- list(formula=~1)
p.const <- list(formula=~1)
Gamma.const <- list(formula=~1)
L.const <- list(formula=~1)

Phi.plate <- list(formula=~geoplate)
p.plate <- list(formula=~geoplate)
Gamma.plate <- list(formula=~geoplate)
L.plate <- list(formula=~geoplate)

## model specifications for RMARK - POPAN model
Phi.time <- list(formula=~time)
p.time <- list(formula=~time)
pent.time <- list(formula=~time)
N.time <- list(formula=~time)

Phi.const <- list(formula=~1)
p.const <- list(formula=~1)
pent.const <- list(formula=~1)
N.const <- list(formula=~1)



##############
## difference between max and min_ma -> length of time bins 
genus_diff <- cbind(genus[,c("occurrence_no","accepted_name","max_ma","min_ma","geoplate")],
                    "diff"=genus$max_ma-genus$min_ma)


## Figure S1 - duration of time bins
hist(genus_C3_S2$diff,
     breaks = 25,
     xlim = c(0,70),
     main = "Length of time bins (global dataset, PBDB)",
     xlab = "Length of time bin in myr",
     cex.lab =1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
min(genus_C3_S2$diff)
max(genus_C3_S2$diff)
median(genus_C3_S2$diff)
mean(genus_C3_S2$diff)

#########################################################################################################
## Pradel model


reruns.gb=list()
reruns.lambda.gb=list()


for (i in 1:100){
  
  tryCatch({
    
    genus_val <- genus_diff %>% group_by(1:n()) %>% mutate(Ma=ifelse(diff<=12, runif(1,min_ma,max_ma),NA))
    base_SS <- c(509,497,485.4,477.7,470,467.3,458.4,453,445.2,443.8,440.8,438.5)
    name_SS<-c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2")
    genus_SS <- as.data.frame(genus_val %>% group_by(Ma) %>% mutate(SS=cut(Ma, 
                                                                           breaks = base_SS, 
                                                                           labels = rev(name_SS),
                                                                           right=FALSE)))
    genus_SS$SS <- factor(genus_SS$SS,
                          levels = c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2"))
    melt_gen <- melt(genus_SS, id=c("accepted_name", "SS"), na.rm = TRUE)
    melt_gen <- na.omit(melt_gen)
    cast_gen <- cast(melt_gen, accepted_name~SS, length)
    
    ## create input for RMARK
    inp <- as.matrix(cast_gen[,2:12])
    inp <- ifelse(inp >=1, 1, 0)
    inp <- unite(as.data.frame(inp), "ch", c(1:11), sep = "")
    
    proc.Berg_pradsen <- process.data(inp, model= "Pradsen")#, groups = "geoplate")#Pradel's survival and "growth"
    time.global <- mark(proc.Berg_pradsen, model.parameters = list(Phi=Phi.time, p=p.time, Gamma=L.time))#, options="SIMANNEAL")
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  reruns.gb[[i]]=time.global$results$real
  
}



for (i in 1:100){
  
  tryCatch({
    
    genus_val <- genus_diff%>% group_by(1:n()) %>% mutate(Ma=ifelse(diff<=12,runif(1,min_ma,max_ma),NA))
    base_SS <- c(509,497,485.4,477.7,470,467.3,458.4,453,445.2,443.8,440.8,438.5)
    name_SS<-c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2")
    genus_SS <- as.data.frame(genus_val %>% group_by(Ma) %>% mutate(SS=cut(Ma, 
                                                                           breaks = base_SS, 
                                                                           labels = rev(name_SS),
                                                                           right=FALSE)))
    genus_SS$SS <- factor(genus_SS$SS,
                          levels = c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2"))
    melt_gen <- melt(genus_SS, id=c("accepted_name", "SS"), na.rm = TRUE)
    melt_gen <- na.omit(melt_gen)
    cast_gen <- cast(melt_gen, accepted_name~SS, length)
    
    inp <- as.matrix(cast_gen[,2:12])
    inp <- ifelse(inp >=1, 1, 0)
    inp <- unite(as.data.frame(inp), "ch", c(1:11), sep = "")
    
    proc.Berg_pradlambda <- process.data(inp, model= "Pradlambda")#, groups = "geoplate")#Pradel's survival and "growth"
    timeD.global <- mark(proc.Berg_pradlambda, model.parameters = list(Phi=Phi.time, p=p.time, Lambda=L.time))#, options="SIMANNEAL")
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  reruns.lambda.gb[[i]]=timeD.global$results$real
  
}


## extracting results from reruns.gb, averaging over 25 runs of the models
## sample means
mean_results <- rowMeans(do.call(cbind, lapply(reruns.gb, function(x) x$estimate)))

estimates.gb <- do.call(cbind, lapply(reruns.gb, function(x) x$estimate))
lcl.gb <- do.call(cbind, lapply(reruns.gb, function(x) x$lcl))
ucl.gb <- do.call(cbind, lapply(reruns.gb, function(x) x$ucl))


## extracting results from reruns.lambda.gb, averaging over 25 runs of the models
mean_results.lambda <- rowMeans(do.call(cbind, lapply(reruns.lambda.gb, function(x) x$estimate)))

estimates.gb.lb <- do.call(cbind, lapply(reruns.lambda.gb, function(x) x$estimate))
lcl.gb.lb <- do.call(cbind, lapply(reruns.lambda.gb, function(x) x$lcl))
ucl.gb.lb <- do.call(cbind, lapply(reruns.lambda.gb, function(x) x$ucl))


###########################################################################################
## fully time-varying model for global data
## extract from averaged results of the models
## phi = survival probability (from Pradsen parametrization)
## p = sampling probability (from Pradsen parametrization)
## gam = seniority probability (from Pradsen parametrization)
## div = growth rate lambda (from Pradlambda parametrization)

phi <- estimates.gb[2:8,]
phiu <- ucl.gb[2:8,]
phil <- lcl.gb[2:8,]

p <- estimates.gb[13:19,]
pu <- ucl.gb[13:19,]
pl <- lcl.gb[13:19,]

gam <- estimates.gb[23:29,]
gamu <- ucl.gb[23:29,]
gaml <- lcl.gb[23:29,]

div <- estimates.gb.lb[23:29,]
divu <- ucl.gb.lb[23:29,]
divl <- lcl.gb.lb[23:29,]

mean.phi <- mean_results[2:8]
mean.p <- mean_results[13:19]
mean.gam <- mean_results[23:29]
mean.div <- mean_results.lambda[23:29]

## time between midpoints of two neighbouring stageslices ## used to calculate origination and extinction rates
t <- c(9.65,7.7,5.2,5.8,7.15,6.6,4.6)
## time per interval/stageslice ## used to calculate sampling rate (p)
tp <- c(7.7,7.7,2.7,8.9,5.4,7.8,1.4)

## conversion of probabilities into rates 
origprob <- 1-gam
origCIl <- 1-gaml
origCIu <- 1-gamu
Orig_rate <- -log(1-origprob)/t
Orig_rate_CIl <- -log(1-origCIl)/t
Orig_rate_CIu <- -log(1-origCIu)/t

extinctprob <- 1-phi
extinctCIl <- 1-phil
extinctCIu <- 1-phiu
Ext_rate <- -log(1-extinctprob)/t
Ext_rate_CIl <- -log(1-extinctCIl)/t
Ext_rate_CIu <- -log(1-extinctCIu)/t

rate_p <- -log(1-p)/tp
ratel_p <- -log(1-pl)/tp
rateu_p <- -log(1-pu)/tp

## mean values
mean.origprob <- 1-mean.gam
mean.origrate <- -log(1-mean.origprob)/t
mean.extprob <- 1-mean.phi
mean.extrate <- -log(1-mean.extprob)/t
mean.prate <- -log(1-mean.p)/tp


#########################################################################################################
## POPAN model

genus <- read.csv("PBDB_Ord_1.csv", sep = ",", header=T)
genus <- filter(genus, grepl("genus", accepted_rank))
genus_diff <- cbind(genus[,c("occurrence_no","accepted_name","max_ma","min_ma","geoplate")],
                    "diff"=genus$max_ma-genus$min_ma)

reruns.POP.gb=list()

for (i in 1:100){
  
  tryCatch({
    
    genus_val <- genus_diff%>% group_by(1:n()) %>% mutate(Ma=ifelse(diff<=12,runif(1,min_ma,max_ma),NA))
    base_SS <- c(509,497,485.4,477.7,470,467.3,458.4,453,445.2,443.8,440.8,438.5)
    name_SS<-c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2")
    genus_SS <- as.data.frame(genus_val %>% group_by(Ma) %>% mutate(SS=cut(Ma, 
                                                                           breaks = base_SS, 
                                                                           labels = rev(name_SS),
                                                                           right=FALSE)))
    genus_SS$SS <- factor(genus_SS$SS,
                          levels = c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2"))
    melt_gen <- melt(genus_SS, id=c("accepted_name", "SS"), na.rm = TRUE)
    melt_gen <- na.omit(melt_gen)
    cast_gen <- cast(melt_gen, accepted_name~SS, length)
    
    inp <- as.matrix(cast_gen[,2:12])
    inp <- ifelse(inp >=1, 1, 0)
    inp <- unite(as.data.frame(inp), "ch", c(1:11), sep = "")
    
    proc.global  <- process.data(inp, model= "POPAN", time.intervals = c(11.8,9.65,7.7,5.2,5.8,7.15,6.6,4.6,2.2,2.65))#
    popan_Phi_time_p_time_pent_time <- mark(proc.global, model.parameters=list(Phi=Phi.time, p=p.time, pent=pent.time, N=N.const))
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  reruns.POP.gb[[i]]=popan_Phi_time_p_time_pent_time$results$derived$`N Population Size`
  
}

meanN_results_global <- rowMeans(do.call(cbind, lapply(reruns.POP.gb, function(x) x$estimate)))

estimates.gb.N <- do.call(cbind, lapply(reruns.POP.gb, function(x) x$estimate))
lcl.gb.N <- do.call(cbind, lapply(reruns.POP.gb, function(x) x$lcl))
ucl.gb.N <- do.call(cbind, lapply(reruns.POP.gb, function(x) x$ucl))


##########################################################
##########################################################
############## plot ### Figure 1 #########################

Stagebase <-c(485.4,477.7,470,467.3,458.4,453,445.2)
Stagemidpoints <- c(481.55,473.85,468.65,462.85,455.7,449.1,444.5)

##
## timescale function from
## http://simpson-carl.github.io/articles/15/timescales.to.base
tpx <- c("Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Sil")
per <- data.frame(c(485.4,477.7,470,467.3,458.4,453,445.2,443.8))
Ord.Stage <- cbind.data.frame(tpx,"Stage"=per)


tscales.Ord <- function(top, bot, s.bot, ...){
  bases <- Ord.Stage
  cc <- rep(c("grey95","grey97 "),length(bases))
  
  rect(xleft = bases[-8, 2], ybottom = rep(bot,7), xright = bases[-1, 2],
       ytop = rep(top, 7), col = cc, border=NA)
  
  rect(xleft = bases[-8, 2], ybottom = rep(bot,7), xright = bases[-1, 2],
       ytop = rep(s.bot, 7), border = 'grey90')
  
  bt <- (bot+s.bot)/2
  tpl <- bases[,2]+c(diff(bases[,2])/2,0)
  
  text(x=tpl[-7], y=bt[-8], labels = bases[1:7, 1])
}


##########################################################
par(mfrow=c(3,1),
    oma = c(1,0,0,0),
    mai = c(0.4, 0.7, 0.1, 0.7),
    cex = 1.1)


plot(Stagebase-0.2, Orig_rate[,27],
     type = "b", 
     pch=16,
     ylim = c(-0.02,0.2),
     xlim = rev(c(444.18,485.4)),
     axes = F,
     xlab = "",
     ylab = "")
tscales.Ord(0.2, 0, -0.02)

lines(Stagebase-0.2, Orig_rate[,27], type = "b", pch= 16, lwd=0.8)
arrows(x0=Stagebase-0.2, y0=Orig_rate_CIl[,27], x1=Stagebase-0.2, y1=Orig_rate_CIu[,27], length=0.02, lwd = 0.8, angle = 90, code = 3)
text(Stagebase[1]-0.6, Orig_rate[1,27]-0.015,labels="O",lwd=0.8)

lines(Stagebase, Ext_rate[,27], type = "b", pch= 1, lwd=1.2,lty=5)
arrows(x0=Stagebase, y0=Ext_rate_CIl[,27], x1=Stagebase, y1=Ext_rate_CIu[,27],length=0.02, lwd = 1.2, angle = 90, code = 3, lty=5)
text(Stagebase[1]-0.4, Ext_rate[1,27]-0.015,labels="E",lwd=0.8)

legend("topleft", legend="A", bty="n", cex = 1.25)

axis(2, col = 'grey75', line = -0.2, at = seq(0, 0.2, 0.05))

mtext("Orig./Ext. events per myr", side = 2, line = 2, cex = 0.75)

##########################################################
plot(Stagebase, div[,27]-1, type = "b", 
     pch=24,
     lwd=0.8,
     ylim = c(-1.02,2),
     xlim = rev(c(444.18,485.4)),
     axes = F,
     xlab = "",
     ylab = "")

tscales.Ord(2, -0.7, -1.02)
abline(h = 0, col="darkgrey")

lines(Stagebase, div[,27]-1, type = "b", pch=24, lwd=0.8)
arrows(x0=Stagebase, y0=divl[,27]-1, x1=Stagebase, y1=divu[,27]-1, length=0.02, lwd = 0.8, angle = 90, code = 3)

axis(2, col = 'grey75', line = -0.2, at = seq(-1, 2, 0.5))

legend("topleft", legend="B", bty="n", cex = 1.25)

mtext("Net diversification rate", side = 2, line = 2, cex = 0.75)

par(new=TRUE)

plot(Stagemidpoints, estimates.gb.N[3:9,27], type = "b",
     pch=15,
     lty=6,
     lwd=1.2,
     ylim = c(500,2000),
     xlim = rev(c(444.18,485.4)),
     axes = F,
     xlab = "",
     ylab = "")


arrows(x0=Stagemidpoints, y0=lcl.gb.N[3:9,27], x1=Stagemidpoints, y1=ucl.gb.N[3:9,27],
       length=0.02, lwd = 1.2, angle = 90, code = 3 ,lty=6)

axis(4, col = 'grey75', line = -1, at = seq(500, 2000, 500))
mtext("Genus richness", side = 4, line = 2, cex = 0.75)

##########################################################
plot(Stagemidpoints, rate_p[,27], type = "b", 
     ylim = c(-0.1,0.9),
     xlim = rev(c(444.18,485.4)),
     axes = F,
     xlab = "",
     ylab = "")

tscales.Ord(0.9, 0, -0.1)

lines(Stagemidpoints, rate_p[,27], type = "b")
arrows(x0=Stagemidpoints, y0=ratel_p[,27], x1=Stagemidpoints, y1=rateu_p[,27], length=0.02, lwd = 1, angle = 90, code = 3)

legend("topleft", legend="C", bty="n", cex = 1.25)

axis(1, col = 'grey75', line = -0.05, at = seq(445,485,10))
axis(2, col = 'grey75', line = -0.2, at = seq(0, 0.9, 0.15))

mtext("Age (Ma)", side = 1, line = 1.5, cex = 0.75)
mtext("Sampling events per myr", side = 2, line = 2, cex = 0.75)

##########################################################
##########################################################
############## plot ### Figure S2 ########################
########### 100 replicate runs of the model ##############

par(mfrow=c(3,1),
    oma = c(1,0,0,0),
    mai = c(0.4, 0.7, 0.1, 0.7),
    cex = 1.1)

plot(Stagebase-0.2, mean.origrate,
     type = "b", 
     pch=16,
     ylim = c(-0.02,0.2),
     xlim = rev(c(444.18,485.4)),
     axes = F,
     xlab = "",
     ylab = "")
tscales.Ord(0.2, 0, -0.02)

for (i in 1:100){
  lines(Stagebase-0.2, Orig_rate_CIl[,i], col = "grey70")}
for (i in 1:100){
  lines(Stagebase-0.2, Orig_rate_CIu[,i], col = "grey70")
}  


for (i in 1:100){
  lines(Stagebase, Ext_rate_CIl[,i], col = "grey55")}
for (i in 1:100){
  lines(Stagebase, Ext_rate_CIu[,i], col = "grey55")}

lines(Stagebase-0.2, mean.origrate, type = "b", pch= 16, lwd=1.5)
lines(Stagebase, mean.extrate, type = "b", pch= 1, lwd=1.5,lty=5)

text(Stagebase[1]-0.6, Orig_rate[1,27]-0.015,labels="O",lwd=0.8)
text(Stagebase[1]-0.4, Ext_rate[1,27]-0.015,labels="E",lwd=0.8)


legend("topleft", legend="A", bty="n", cex = 1.25)

axis(2, col = 'grey75', line = -0.2, at = seq(0, 0.2, 0.05))

mtext("Orig./Ext. events per myr", side = 2, line = 2, cex = 0.75)

##########################################################
plot(Stagebase, mean.div-1, type = "b", 
     pch=24,
     lwd=0.8,
     ylim = c(-1.02,2),
     xlim = rev(c(444.18,485.4)),
     axes = F,
     xlab = "",
     ylab = "")

tscales.Ord(2, -0.7, -1.02)
abline(h = 0, col="white", lwd = 2)


for (i in 1:100){
  lines(Stagebase, divl[,i]-1, col = "grey70")}
for (i in 1:100){
  lines(Stagebase, divu[,i]-1, col = "grey70")
}  

lines(Stagebase, mean.div-1, type = "b", pch=24, lwd=0.8)

axis(2, col = 'grey75', line = -0.2, at = seq(-1, 2, 0.5))

legend("topleft", legend="B", bty="n", cex = 1.25)

mtext("Net diversification rate", side = 2, line = 2, cex = 0.75)

##########################################################
par(new=TRUE)

plot(Stagemidpoints, meanN_results_global[3:9], type = "b",
     pch=15,
     lty=6,
     lwd=1.2,
     ylim = c(500,2000),
     xlim = rev(c(444.18,485.4)),
     axes = F,
     xlab = "",
     ylab = "")


for (i in 1:100){
  lines(Stagemidpoints, lcl.gb.N[3:9,i], col = "grey55")}
for (i in 1:100){
  lines(Stagemidpoints, ucl.gb.N[3:9,i], col = "grey55")}
lines(Stagemidpoints,meanN_results_global[3:9],type="b")


axis(4, col = 'grey75', line = 0, at = seq(500, 2000, 500))

mtext("Genus richness", side = 4, line = 2, cex = 0.75)

##########################################################
plot(Stagemidpoints, mean.prate, type = "b", 
     ylim = c(-0.1,0.9),
     xlim = rev(c(444.18,485.4)),
     axes = F,
     xlab = "",
     ylab = "")

tscales.Ord(0.9, 0, -0.1)

for (i in 1:100){
  lines(Stagemidpoints, ratel_p[,i], col = "grey55")}
for (i in 1:100){
  lines(Stagemidpoints, rateu_p[,i], col = "grey55")}
  lines(Stagemidpoints,mean.prate, type="b")

legend("topleft", legend="C", bty="n", cex = 1.25)

axis(1, col = 'grey75', line = -0.05)
axis(2, col = 'grey75', line = -0.2, at = seq(0, 0.9, 0.15))

mtext("Age (Ma)", side = 1, line = 1.5, cex = 0.75)
mtext("Sampling events per myr", side = 2, line = 2, cex = 0.75)


#########################################################################################################
#### phyla represented by genus occurrences #############################################################

genusg <- filter(genus, grepl("genus", accepted_rank))
genus_Ord <- filter(genusg, max_ma < 481.55 & max_ma > 443.8 | min_ma < 481.55 & min_ma > 443.8)
genus_phyla <- unique(genus_Ord[,c("accepted_name","phylum")])

n_phyla <- aggregate(list(genus_phyla$accepted_name),
                     list(genus_phyla$phylum),
                     length)

n_phyla <- arrange(n_phyla,desc(n_phyla[,2]))#, n = 10)
n_phyla[,2] <- n_phyla[,2]*100/colSums(n_phyla[2], na.rm = T)
n_phyla$group <- ifelse(n_phyla[,2]>1, 1, 2)
n_phyla1 <- n_phyla[n_phyla$group==1,]
n_phyla1 <- n_phyla1 %>% filter(!Group.1 %in% c("","Problematica"))


colnames(n_phyla1)[1] <- "Phylum"
colnames(n_phyla1)[2] <- "percent"

n_phyla1$Phylum <- factor(n_phyla1$Phylum,
                        levels = c(levels(n_phyla1$Phylum),"Others"))

n_phyla1[nrow(n_phyla1)+1,1] <- c("Others")
n_phyla1[nrow(n_phyla1),2] <- c(100-colSums(n_phyla1[2], na.rm = T))
Phyla_global <- n_phyla1[,1:2]

## pie plot of all phyla that are represented by the Ordovician data in the dataset
ggplot(Phyla_global, aes(x="", y=percent, fill=Phylum)) +
  geom_bar(stat="identity", colour = "white") +
  # geom_text(aes(x= factor(1), y=pos, label = Phylum), size=10) +  # note y = pos
  theme_minimal() +
  scale_fill_grey(start = 0, end = 0.9, na.value = "white")+
  coord_polar("y", start=0) +
  ggtitle("Phyla represented by Genera globally")


#########################################################################################################
## model comparison #####################################################################################

genus_val <- genus_diff %>% group_by(1:n()) %>% mutate(Ma=ifelse(diff<=12, runif(1,min_ma,max_ma),NA))
base_SS <- c(509,497,485.4,477.7,470,467.3,458.4,453,445.2,443.8,440.8,438.5)
name_SS<-c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2")
genus_SS <- as.data.frame(genus_val %>% group_by(Ma) %>% mutate(SS=cut(Ma, 
                                                                       breaks = base_SS, 
                                                                       labels = rev(name_SS),
                                                                       right=FALSE)))
genus_SS$SS <- factor(genus_SS$SS,
                      levels = c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2"))
melt_gen <- melt(genus_SS, id=c("accepted_name", "SS"), na.rm = TRUE)
melt_gen <- na.omit(melt_gen)
cast_gen <- cast(melt_gen, accepted_name~SS, length)

## create input for RMARK
inp <- as.matrix(cast_gen[,2:12])
inp <- ifelse(inp >=1, 1, 0)
inp <- unite(as.data.frame(inp), "ch", c(1:11), sep = "")

##########################################################
##########################################################
## model comparison based on Pradsen model
## Table S1

proc.Berg_pradsen <- process.data(inp, model= "Pradsen")

time.global <- mark(proc.Berg_pradsen, model.parameters = list(Phi=Phi.time, p=p.time, Gamma=L.time))
phi.t.p.1.g.t <- mark(proc.Berg_pradsen, model.parameters = list(Phi=Phi.time, p=p.const, Gamma=L.time))
phi.1.p.t.g.t <- mark(proc.Berg_pradsen, model.parameters = list(Phi=Phi.const, p=p.time, Gamma=L.time))
phi.t.p.t.g.1 <- mark(proc.Berg_pradsen, model.parameters = list(Phi=Phi.time, p=p.time, Gamma=L.const))
phi.1.p.1.g.t <- mark(proc.Berg_pradsen, model.parameters = list(Phi=Phi.const, p=p.const, Gamma=L.time))
phi.1.p.t.g.1 <- mark(proc.Berg_pradsen, model.parameters = list(Phi=Phi.const, p=p.time, Gamma=L.const))
phi.t.p.1.g.1 <- mark(proc.Berg_pradsen, model.parameters = list(Phi=Phi.time, p=p.const, Gamma=L.const))
## compare all models in following table
comparison <- collect.models(type = "Pradsen")
comparison
##########################################################
##########################################################

