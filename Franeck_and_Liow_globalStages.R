library(reshape2) 
library(reshape)
library(plyr)
library(dplyr)
library(tidyr)
library(RMark)
library(ggplot2)
library(stringr)

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
## only considering occurrences that are already assigned to global Ordovician stages.

gS <- genus %>% filter(str_detect(early_interval, "Tremadoc|Tremadocian|Floian|Dapingian|Darriwilian|Sandbian|Katian|Hirnantian")|
                         str_detect(late_interval, "Tremadoc|Tremadocian|Floian|Dapingian|Darriwilian|Sandbian|Katian|Hirnantian"))

## difference between max and min_ma -> length of time bins 
genus_diff <- cbind(gS[,c("occurrence_no","accepted_name","max_ma","min_ma","geoplate")],
                    "diff"=gS$max_ma-gS$min_ma)

#########################################################################################################
## Pradel model


reruns.gb=list()
reruns.lambda.gb=list()


for (i in 1:10){
  
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
    inp <- as.matrix(cast_gen[,2:10])
    inp <- ifelse(inp >=1, 1, 0)
    inp <- unite(as.data.frame(inp), "ch", c(1:9), sep = "")
    
    proc.Berg_pradsen <- process.data(inp, model= "Pradsen")#, groups = "geoplate")#Pradel's survival and "growth"
    time.global <- mark(proc.Berg_pradsen, model.parameters = list(Phi=Phi.time, p=p.time, Gamma=L.time))#, options="SIMANNEAL")
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  reruns.gb[[i]]=time.global$results$real
  
}



for (i in 1:10){
  
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
    
    inp <- as.matrix(cast_gen[,2:10])
    inp <- ifelse(inp >=1, 1, 0)
    inp <- unite(as.data.frame(inp), "ch", c(1:9), sep = "")
    
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

phi <- estimates.gb[1:7,]
phiu <- ucl.gb[1:7,]
phil <- lcl.gb[1:7,]

p <- estimates.gb[9:15,]
pu <- ucl.gb[9:15,]
pl <- lcl.gb[9:15,]

gam <- estimates.gb[18:24,]
gamu <- ucl.gb[18:24,]
gaml <- lcl.gb[18:24,]

div <- estimates.gb.lb[18:24,]
divu <- ucl.gb.lb[18:24,]
divl <- lcl.gb.lb[18:24,]

mean.phi <- mean_results[1:7]
mean.p <- mean_results[9:15]
mean.gam <- mean_results[18:24]
mean.div <- mean_results.lambda[18:24]

## time between midpoints of two neighbouring stageslices ## used to calculate origination and extinction rates
t <- c(7.7,5.2,5.8,7.15,6.6,4.6,2.2)

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

gS <- genus %>% filter(str_detect(early_interval, "Tremadoc|Tremadocian|Floian|Dapingian|Darriwilian|Sandbian|Katian|Hirnantian")|
                         str_detect(late_interval, "Tremadoc|Tremadocian|Floian|Dapingian|Darriwilian|Sandbian|Katian|Hirnantian"))

## difference between max and min_ma -> length of time bins 
genus_diff <- cbind(gS[,c("occurrence_no","accepted_name","max_ma","min_ma","geoplate")],
                    "diff"=gS$max_ma-gS$min_ma)

# genus_diff <- cbind(genus[,c("occurrence_no","accepted_name","max_ma","min_ma","geoplate")],
#                     "diff"=genus$max_ma-genus$min_ma)

reruns.POP.gb=list()

for (i in 1:10){
  
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
    
    inp <- as.matrix(cast_gen[,2:10])
    inp <- ifelse(inp >=1, 1, 0)
    inp <- unite(as.data.frame(inp), "ch", c(1:9), sep = "")
    
    proc.global  <- process.data(inp, model= "POPAN", time.intervals = c(7.7,5.2,5.8,7.15,6.6,4.6,2.2,2.65))#
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
############## plot ### Figure S6 ########################

Stagebase <-c(477.7,470,467.3,458.4,453,445.2,443.8)
Stagemidpoints <- c(473.85,468.65,462.85,455.7,449.1,444.5)
StagemidpointsN <- c(481.55,473.85,468.65,462.85,455.7,449.1,444.5)

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
    oma = c(3,0,0,0),
    mai = c(0, 0.7, 0.1, 0.7),
    cex = 0.85)


plot(Stagebase-0.2, Orig_rate[,2],
     type = "b", 
     pch=16,
     ylim = c(-0.05,0.2),
     xlim = rev(c(444.18,485.4)),
     axes = F,
     xlab = "",
     ylab = "")
tscales.Ord(0.5, -0.001, -0.05)

lines(Stagebase-0.2, Orig_rate[,2], type = "b", pch= 16, lwd=0.8)
arrows(x0=Stagebase-0.2, y0=Orig_rate_CIl[,2], x1=Stagebase-0.2, y1=Orig_rate_CIu[,2], length=0.02, lwd = 0.8, angle = 90, code = 3)
text(Stagebase[1]+1, Orig_rate[1,2]+0.015,labels="O",lwd=0.8, cex = 0.9)

lines(Stagebase, Ext_rate[,2], type = "b", pch= 1, lwd=1.2,lty=5)
arrows(x0=Stagebase, y0=Ext_rate_CIl[,2], x1=Stagebase, y1=Ext_rate_CIu[,2],length=0.02, lwd = 1.2, angle = 90, code = 3, lty=5)
text(Stagebase[1]+1, Ext_rate[1,2]+0.045,labels="E",lwd=0.8,cex = 0.9)

legend("topleft", legend="A", bty="n", cex = 1.25)

axis(2, col = 'grey75', line = -0.2, at = seq(0, 0.2, 0.05))

mtext("Orig./Ext. events\n per myr", side = 2, line = 2, cex = 0.8)

##########################################################
plot(Stagebase, div[,2]-1, type = "b", 
     pch=24,
     lwd=0.8,
     ylim = c(-1.6,1.5),
     xlim = rev(c(444.18,485.4)),
     axes = F,
     xlab = "",
     ylab = "")

tscales.Ord(1.7, -1, -1.6)
abline(h = 0, col="darkgrey")

lines(Stagebase, div[,2]-1, type = "b", pch=24, lwd=0.8)
arrows(x0=Stagebase, y0=divl[,2]-1, x1=Stagebase, y1=divu[,2]-1, length=0.02, lwd = 0.8, angle = 90, code = 3)

axis(2, col = 'grey75', line = -0.2, at = seq(-1, 1.5, 0.5))

legend("topleft", legend="B", bty="n", cex = 1.25)

mtext("Net diversification rate", side = 2, line = 2, cex = 0.8)

par(new=TRUE)

plot(StagemidpointsN, estimates.gb.N[1:7,2], type = "b",
     pch=15,
     lty=6,
     lwd=1.2,
     ylim = c(100,1000),
     xlim = rev(c(444.18,485.4)),
     axes = F,
     xlab = "",
     ylab = "")


arrows(x0=StagemidpointsN, y0=lcl.gb.N[1:7,2], x1=StagemidpointsN, y1=ucl.gb.N[1:7,2],
       length=0.02, lwd = 1.2, angle = 90, code = 3 ,lty=6)

axis(4, col = 'grey75', line = 0, at = seq(100, 1000, 200))
mtext("Genus richness", side = 4, line = 2, cex = 0.8)

##########################################################
plot(Stagemidpoints, rate_p[1:6,2], type = "b", 
     ylim = c(-0.1,0.3),
     xlim = rev(c(444.18,485.4)),
     axes = F,
     xlab = "",
     ylab = "")

tscales.Ord(0.3, -0.01, -0.1)

lines(Stagemidpoints, rate_p[1:6,2], type = "b")
arrows(x0=Stagemidpoints, y0=ratel_p[1:6,2], x1=Stagemidpoints, y1=rateu_p[1:6,2], length=0.02, lwd = 1, angle = 90, code = 3)

legend("topleft", legend="C", bty="n", cex = 1.25)

axis(1, col = 'grey75', line = -0.05, at = seq(445,485,10))
axis(2, col = 'grey75', line = -0.2, at = seq(0, 0.3, 0.1))

mtext("Age (Ma)", side = 1, line = 2, cex = 0.8)
mtext("Sampling events\n per myr", side = 2, line = 2, cex = 0.8)
