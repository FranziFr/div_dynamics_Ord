library(reshape2) 
library(reshape)
library(plyr)
library(dplyr)
library(tidyr)
library(RMark)
library(ggplot2)
library(plotrix)
library(stringr)


## set working directory to load the PBDB file
# setwd("M:/Franeck_and_Liow_2018")
setwd("~/offline/Franeck_and_Liow_2018")

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
Phi.Time.plate <- list(formula=~time+geoplate+time*geoplate)
p.Time.plate <- list(formula=~time+geoplate+time*geoplate)
pent.Time.plate <- list(formula=~time+geoplate+time*geoplate)
N.time <- list(formula=~time)
N.const <- list(formula=~1)
N.plate <- list(formula=~geoplate)

#####################################################################################################
## extract data from Laurentia
gen_Lau <- genus %>% filter(geoplate %in% c("101"))
gen_Lau_diff <- cbind(gen_Lau[,c("occurrence_no","accepted_name","max_ma","min_ma","geoplate","ref_pubyr")],
                      "diff"=gen_Lau$max_ma-gen_Lau$min_ma)

## extract data from Baltica
gen_Bal <- genus %>% filter(str_detect(geoplate,'302')|
                              str_detect(cc, 'SE')|
                              str_detect(cc, 'EE')|
                              str_detect(cc, 'FI')|
                              str_detect(cc, 'NO')|
                              str_detect(state, "St. Petersburg"),
                            !str_detect(state,'Spitsbergen|Svalbard'),
                            !str_detect(geoplate, '311'))

gen_Bal_diff <- cbind(gen_Bal[,c("occurrence_no","accepted_name","max_ma","min_ma","geoplate","ref_pubyr")],
                      "diff"=gen_Bal$max_ma-gen_Bal$min_ma)


## extract all data EXCEPT from Baltica or Laurentia 
gen_cont <- genus %>% filter(!str_detect(geoplate, "101|302"),
                             !str_detect(cc, "SE|EE|FI|NO"),
                             !str_detect(state, "St. Petersburg")|
                             str_detect(state,'Spitsbergen|Svalbard')|
                             str_detect(geoplate, "311"))
                                              

gen_cont_diff <- cbind(gen_cont[,c("occurrence_no","accepted_name","max_ma","min_ma","geoplate","ref_pubyr")],
                       "diff"=gen_cont$max_ma-gen_cont$min_ma)



reruns=list()
reruns.lambda=list()

#########################################################################################################
## Pradel model
## Pradsen parametrization
for (i in 1:100){
  
  tryCatch({
    
    ## random assignment of a specific age betwenn min and max age of the geological unit the fossil was found in, given the unit is longer than 12 myrs
    gen_Lau_val <- gen_Lau_diff%>% group_by(1:n()) %>% mutate(Ma=ifelse(diff<=12,runif(1,min_ma,max_ma),NA))
    
    ## definition of bases of stages and their names 
    base_SS <- c(509,497,485.4,477.7,470,467.3,458.4,453,445.2,443.8,440.8,438.5)
    name_SS<-c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2")
    
    ## adding a new column with names of Stages for the randomly assigned ages
    gen_Lau_val_SS <- as.data.frame(gen_Lau_val %>% group_by(Ma) %>% mutate(SS=cut(Ma, 
                                                                                   breaks = base_SS, 
                                                                                   labels = rev(name_SS),
                                                                                   right=FALSE)))
    gen_Lau_val_SS$SS <- factor(gen_Lau_val_SS$SS,
                                levels = c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2"))
    
    ## creating a matrix/dataframe that consists of the genus names (lines) and stages (columns) with counts of how many observations were done in each stage of each genus
    melt_Lau <- melt(gen_Lau_val_SS, id=c("accepted_name", "SS"), na.rm = TRUE)
    melt_Lau <- na.omit(melt_Lau)
    cast_Lau <- cast(melt_Lau, accepted_name~SS, length)
    #####################################################################################################
    ## Baltica
    ## commands are as above, for Laurentia
    ## 
    gen_Bal_val <- gen_Bal_diff%>% group_by(1:n()) %>% mutate(Ma=ifelse(diff<=12,runif(1,min_ma,max_ma),NA))
    gen_Bal_val_SS <- as.data.frame(gen_Bal_val %>% group_by(Ma) %>% mutate(SS=cut(Ma, 
                                                                                   breaks = base_SS, 
                                                                                   labels = rev(name_SS),
                                                                                   right=FALSE)))
    gen_Bal_val_SS$SS <- factor(gen_Bal_val_SS$SS,
                                levels = c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2"))
    
    melt_Bal <- melt(gen_Bal_val_SS, id=c("accepted_name", "SS"), na.rm = TRUE)
    melt_Bal <- na.omit(melt_Bal)
    cast_Bal <- cast(melt_Bal, accepted_name~SS, length)
    
    
    #####################################################################################################
    #####################################################################################################
    ## the same procedure for data from all other continents, except Baltica and Laurentia
    ## commands are as for Laurentia
    ## note that we include all genera here (not putting constraints on duration of time-bin)
    ## this is because we are not interested in the distribution, but in the names of genera that are overlapping in general 
    gen_cont_val <- gen_cont_diff%>% group_by(1:n()) %>% mutate(Ma=ifelse(diff<=56,runif(1,min_ma,max_ma),NA))
    gen_cont_val_SS <- as.data.frame(gen_cont_val %>% group_by(Ma) %>% mutate(SS=cut(Ma,
                                                                                     breaks = base_SS, 
                                                                                     labels = rev(name_SS),
                                                                                     right=FALSE)))
    gen_cont_val_SS$SS <- factor(gen_cont_val_SS$SS,
                                 levels = c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2"))
    
    melt_cont <- melt(gen_cont_val_SS, id=c("accepted_name", "SS"), na.rm = TRUE)
    melt_cont <- na.omit(melt_cont)
    cast_cont <- cast(melt_cont, accepted_name~SS, length)
    
    ## now we can check for overlapping genera between Baltica and Laurentia, Baltica, and all other data, and Laurentia and all other data
    ## in this order
    overlap <- as.factor(Reduce(intersect, list(cast_Bal$accepted_name, cast_Lau$accepted_name)))
    overlap_Bcont <- as.factor(Reduce(intersect, list(cast_Bal$accepted_name, cast_cont$accepted_name)))
    overlap_Lcont <- as.factor(Reduce(intersect, list(cast_Lau$accepted_name, cast_cont$accepted_name)))
    
    ## from the cast_Laurentia table we filter now all genera that are NOT included in the overlap/overlap_Lcont data.frames
    inp_L <- cast_Lau %>% filter(!accepted_name %in% overlap)
    inp_L <- inp_L  %>% filter(!accepted_name %in% overlap_Lcont)
    
    ## now we substitute all numbers with 1 and all NAs/0s with 0 and create so the presence/absence matrix
    inp_Lmat <- data.frame(inp_L[,1],ifelse(inp_L[,2:12]>=1,1,0))
    
    ## finally, we collapse the matrix into an RMARK input file
    inp_L <- as.matrix(inp_L[,2:12])
    inp_L <- ifelse(inp_L >=1, 1, 0)
    inp_L <- unite(as.data.frame(inp_L), "ch", c(1:11), sep = "")
    inp_L <- cbind(inp_L, geoplate = "NAC", ";")
    
    
    ## we do the same for Baltica
    inp_B <- cast_Bal %>% filter(!accepted_name %in% overlap)
    inp_B <- inp_B  %>% filter(!accepted_name %in% overlap_Bcont)
    inp_Bmat <- data.frame(inp_B[,1],ifelse(inp_B[,2:12]>=1,1,0))
    
    inp_B <- as.matrix(inp_B[,2:12])
    inp_B <- ifelse(inp_B >=1, 1, 0)
    inp_B <- unite(as.data.frame(inp_B), "ch", c(1:11), sep = "")
    inp_B <- cbind(inp_B, geoplate = "BAL", ";")
    
    
    inp <- rbind.data.frame(inp_B, inp_L)
    inp$geoplate <- factor(inp$geoplate,
                           levels = c("BAL", "NAC"))
    
    proc.Berg_pradsen <- process.data(inp, model= "Pradsen", groups = "geoplate")#Pradel's survival and "growth"
    Time.plate <- mark(proc.Berg_pradsen, model.parameters = list(Phi=Phi.Time.plate, p=p.Time.plate, Gamma=L.Time.plate))#, options="SIMANNEAL")
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  reruns[[i]]=Time.plate$results$real
  
}



## Pradlambda parametrization 
## commands to create the input file are as above for the Pradsen parametrization

for (i in 1:100){
  
  tryCatch({
    #####################################################################################################
    ## Laurentia
    gen_Lau_val <- gen_Lau_diff%>% group_by(1:n()) %>% mutate(Ma=ifelse(diff<=12,runif(1,min_ma,max_ma),NA))
    base_SS <- c(509,497,485.4,477.7,470,467.3,458.4,453,445.2,443.8,440.8,438.5)
    name_SS<-c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2")
    gen_Lau_val_SS <- as.data.frame(gen_Lau_val %>% group_by(Ma) %>% mutate(SS=cut(Ma, 
                                                                                   breaks = base_SS, 
                                                                                   labels = rev(name_SS),
                                                                                   right=FALSE)))
    gen_Lau_val_SS$SS <- factor(gen_Lau_val_SS$SS,
                                levels = c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2"))
    melt_Lau <- melt(gen_Lau_val_SS, id=c("accepted_name", "SS"), na.rm = TRUE)
    melt_Lau <- na.omit(melt_Lau)
    cast_Lau <- cast(melt_Lau, accepted_name~SS, length)
    
    #####################################################################################################
    ## Baltica
    
    gen_Bal_val <- gen_Bal_diff%>% group_by(1:n()) %>% mutate(Ma=ifelse(diff<=12,runif(1,min_ma,max_ma),NA))
    gen_Bal_val_SS <- as.data.frame(gen_Bal_val %>% group_by(Ma) %>% mutate(SS=cut(Ma, 
                                                                                   breaks = base_SS, 
                                                                                   labels = rev(name_SS),
                                                                                   right=FALSE)))
    gen_Bal_val_SS$SS <- factor(gen_Bal_val_SS$SS,
                                levels = c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2"))
    
    melt_Bal <- melt(gen_Bal_val_SS, id=c("accepted_name", "SS"), na.rm = TRUE)
    melt_Bal <- na.omit(melt_Bal)
    cast_Bal <- cast(melt_Bal, accepted_name~SS, length)
    
  
    #####################################################################################################
    ## all other data, except data from Baltica or Laurentia
    gen_cont_val <- gen_cont_diff%>% group_by(1:n()) %>% mutate(Ma=ifelse(diff<=56,runif(1,min_ma,max_ma),NA))
    gen_cont_val_SS <- as.data.frame(gen_cont_val %>% group_by(Ma) %>% mutate(SS=cut(Ma,
                                                                                     breaks = base_SS, 
                                                                                     labels = rev(name_SS),
                                                                                     right=FALSE)))
    gen_cont_val_SS$SS <- factor(gen_cont_val_SS$SS,
                                 levels = c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2"))
    
    melt_cont <- melt(gen_cont_val_SS, id=c("accepted_name", "SS"), na.rm = TRUE)
    melt_cont <- na.omit(melt_cont)
    cast_cont <- cast(melt_cont, accepted_name~SS, length)
    
    #####################################################################################################
    overlap <- as.factor(Reduce(intersect, list(cast_Bal$accepted_name, cast_Lau$accepted_name)))
    overlap_Bcont <- as.factor(Reduce(intersect, list(cast_Bal$accepted_name, cast_cont$accepted_name)))
    overlap_Lcont <- as.factor(Reduce(intersect, list(cast_Lau$accepted_name, cast_cont$accepted_name)))
    
    
    inp_L <- cast_Lau %>% filter(!accepted_name %in% overlap)
    inp_L <- inp_L  %>% filter(!accepted_name %in% overlap_Lcont)
    inp_Lmat <- data.frame(inp_L[,1],ifelse(inp_L[,2:12]>=1,1,0))
    
    inp_L <- as.matrix(inp_L[,2:12])
    inp_L <- ifelse(inp_L >=1, 1, 0)
    inp_L <- unite(as.data.frame(inp_L), "ch", c(1:11), sep = "")
    inp_L <- cbind(inp_L, geoplate = "NAC", ";")
    
    inp_B <- cast_Bal %>% filter(!accepted_name %in% overlap)
    inp_B <- inp_B  %>% filter(!accepted_name %in% overlap_Bcont)
    inp_Bmat <- data.frame(inp_B[,1],ifelse(inp_B[,2:12]>=1,1,0))
    
    inp_B <- as.matrix(inp_B[,2:12])
    inp_B <- ifelse(inp_B >=1, 1, 0)
    inp_B <- unite(as.data.frame(inp_B), "ch", c(1:11), sep = "")
    inp_B <- cbind(inp_B, geoplate = "BAL", ";")
    
    
    inp <- rbind.data.frame(inp_B, inp_L)
    inp$geoplate <- factor(inp$geoplate,
                           levels = c("BAL", "NAC"))
    
    proc.Berg_pradlambda <- process.data(inp, model= "Pradlambda", groups = "geoplate")#Pradel's survival and "growth"
    TimeD.plate <- mark(proc.Berg_pradlambda, model.parameters = list(Phi=Phi.Time.plate, p=p.Time.plate, Lambda=L.Time.plate))#, options="SIMANNEAL")
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  reruns.lambda[[i]]=TimeD.plate$results$real
  
}

##################
##################
##################

## now we estimate the mean for all model estimates after 100 runs
## Pradsen parametrization
## mean parameter estimates
mean_results <- rowMeans(do.call(cbind, lapply(reruns, function(x) x$estimate)))

## collecting 100 replicate parameter estimates
estimates.pc <- do.call(cbind, lapply(reruns, function(x) x$estimate))
lcl.pc <- do.call(cbind, lapply(reruns, function(x) x$lcl))
ucl.pc <- do.call(cbind, lapply(reruns, function(x) x$ucl))

## Pradlambda parametrization
mean_results.lambda <- rowMeans(do.call(cbind, lapply(reruns.lambda, function(x) x$estimate)))

estimates.pc.lambda <- do.call(cbind, lapply(reruns.lambda, function(x) x$estimate))
lcl.pc.lambda <- do.call(cbind, lapply(reruns.lambda, function(x) x$lcl))
ucl.pc.lambda <- do.call(cbind, lapply(reruns.lambda, function(x) x$ucl))

estimates.pc[estimates.pc==0|estimates.pc==1] <- NA
lcl.pc[is.na(estimates.pc)] <- NA
ucl.pc[is.na(estimates.pc)] <- NA

estimates.pc[ucl.pc-lcl.pc==0|ucl.pc-lcl.pc==1]<-NA
lcl.pc[is.na(estimates.pc)] <- NA
ucl.pc[is.na(estimates.pc)] <- NA


estimates.pc.lambda[ucl.pc.lambda-lcl.pc.lambda>10] <- NA
estimates.pc.lambda[ucl.pc.lambda-lcl.pc.lambda==0|ucl.pc.lambda-lcl.pc.lambda==1]<-NA
lcl.pc.lambda[is.na(estimates.pc.lambda)] <- NA
ucl.pc.lambda[is.na(estimates.pc.lambda)] <- NA



###########################################################################################
## extracting parameter estimates from the results above
phi_BAL <- estimates.pc[2:8,]
phiu_BAL <- ucl.pc[2:8,]
phil_BAL <- lcl.pc[2:8,]

phi_NAC <- estimates.pc[12:18,]
phiu_NAC <- ucl.pc[12:18,]
phil_NAC <- lcl.pc[12:18,]

p_BAL <- estimates.pc[23:29,]
pu_BAL <- ucl.pc[23:29,]
pl_BAL <- lcl.pc[23:29,]

p_NAC <- estimates.pc[34:40,]
pu_NAC <- ucl.pc[34:40,]
pl_NAC <- lcl.pc[34:40,]

gam_BAL <- estimates.pc[44:50,]
gamu_BAL <- ucl.pc[44:50,]
gaml_BAL <- lcl.pc[44:50,]

gam_NAC <- estimates.pc[54:60,]
gamu_NAC <- ucl.pc[54:60,]
gaml_NAC <- lcl.pc[54:60,]

div_BAL <- estimates.pc.lambda[44:50,]
divu_BAL <- ucl.pc.lambda[44:50,]
divl_BAL <- lcl.pc.lambda[44:50,]

div_NAC <- estimates.pc.lambda[54:60,]
divu_NAC <- ucl.pc.lambda[54:60,]
divl_NAC <- lcl.pc.lambda[54:60,]


mean.phi.B <- mean_results[2:8]
mean.p.B <- mean_results[23:29]
mean.gam.B <- mean_results[44:50]
mean.div.B <- mean_results.lambda[44:50]

mean.phi.L <- mean_results[12:18]
mean.p.L <- mean_results[34:40]
mean.gam.L <- mean_results[54:60]
mean.div.L <- mean_results.lambda[54:60]

## time between midpoints of two neighbouring stageslices 
t <- c(9.65,7.7,5.2,5.8,7.15,6.6,4.6)
## time per interval/stageslice
tp <- c(7.7,7.7,2.7,8.9,5.4,7.8,1.4)


###########################################################################################
## transformation of mean parameter estimates into rates
BAL_origprob <- 1-gam_BAL
BAL_origCIl <- 1-gaml_BAL
BAL_origCIu <- 1-gamu_BAL
BAL_Orig_rate <- -log(1-BAL_origprob)/t
BAL_Orig_rate_CIl <- -log(1-BAL_origCIl)/t
BAL_Orig_rate_CIu <- -log(1-BAL_origCIu)/t

NAC_origprob <- 1-gam_NAC
NAC_origCIl <- 1-gaml_NAC
NAC_origCIu <- 1-gamu_NAC
NAC_Orig_rate <- -log(1-NAC_origprob)/t
NAC_Orig_rate_CIl <- -log(1-NAC_origCIl)/t
NAC_Orig_rate_CIu <- -log(1-NAC_origCIu)/t


BAL_extinctprob <- 1-phi_BAL
BAL_extinctCIl <- 1-phil_BAL
BAL_extinctCIu <- 1-phiu_BAL
BAL_Ext_rate <- -log(1-BAL_extinctprob)/t
BAL_Ext_rate_CIl <- -log(1-BAL_extinctCIl)/t
BAL_Ext_rate_CIu <- -log(1-BAL_extinctCIu)/t

NAC_extinctprob <- 1-phi_NAC
NAC_extinctCIl <- 1-phil_NAC
NAC_extinctCIu <- 1-phiu_NAC
NAC_Ext_rate <- -log(1-NAC_extinctprob)/t
NAC_Ext_rate_CIl <- -log(1-NAC_extinctCIl)/t
NAC_Ext_rate_CIu <- -log(1-NAC_extinctCIu)/t


rate_p_BAL <- -log(1-p_BAL)/tp
ratel_p_BAL <- -log(1-pl_BAL)/tp
rateu_p_BAL <- -log(1-pu_BAL)/tp

rate_p_LAU <- -log(1-p_NAC)/tp
ratel_p_LAU <- -log(1-pl_NAC)/tp
rateu_p_LAU <- -log(1-pu_NAC)/tp


## mean rates
mean.origprob.B <- 1-mean.gam.B
mean.origrate.B <- -log(1-mean.origprob.B)/t
mean.extprob.B <- 1-mean.phi.B
mean.extrate.B <- -log(1-mean.extprob.B)/t
mean.prate.B <- -log(1-mean.p.B)/tp

mean.origprob.L <- 1-mean.gam.L
mean.origrate.L <- -log(1-mean.origprob.L)/t
mean.extprob.L <- 1-mean.phi.L
mean.extrate.L <- -log(1-mean.extprob.L)/t
mean.prate.L <- -log(1-mean.p.L)/tp


#########################################################################################################
## run the following section in order to get time-constant parameter estimates for survivial and seniority
## for the whole Ordovician (only considering Ordovician data)
## plate.model
# 
# genus_Ord <- filter(genus, max_ma < 481.55 & max_ma > 443.8 | min_ma < 481.55 & min_ma >= 443.8)
# 
# 
# gen_Lau <- genus_Ord %>% filter(geoplate %in% c("101"))
# gen_Lau_diff <- cbind(gen_Lau[,c("occurrence_no","accepted_name","max_ma","min_ma","geoplate")],
#                       "diff"=gen_Lau$max_ma-gen_Lau$min_ma)
# 
# # extract data from Baltica
# gen_Bal <- genus_Ord %>% filter(str_detect(geoplate,'302')|
#                                   str_detect(cc, 'SE')|
#                                   str_detect(cc, 'EE')|
#                                   str_detect(cc, 'FI')|
#                                   str_detect(cc, 'NO')|
#                                   str_detect(state, "St. Petersburg"),
#                                 !str_detect(state,'Spitsbergen|Svalbard'),
#                                 !str_detect(geoplate, '311'))
# 
# gen_Bal_diff <- cbind(gen_Bal[,c("occurrence_no","accepted_name","max_ma","min_ma","geoplate")],
#                       "diff"=gen_Bal$max_ma-gen_Bal$min_ma)
# 
# # extract all data EXCEPT from Baltica or Laurentia
# gen_cont <- genus_Ord %>% genus %>% filter(!str_detect(geoplate, "101|302"),
#                                            !str_detect(cc, "SE|EE|FI|NO"),
#                                            !str_detect(state, "St. Petersburg")|
#                                              str_detect(state,'Spitsbergen|Svalbard')|
#                                              str_detect(geoplate, "311"))
# gen_cont_diff <- cbind(gen_cont[,c("occurrence_no","accepted_name","max_ma","min_ma","geoplate","ref_pubyr")],
#                        "diff"=gen_cont$max_ma-gen_cont$min_ma)
# 
# 
# 
# ###############
# ## Laurentia
# 
# gen_Lau_val <- gen_Lau_diff%>% group_by(1:n()) %>% mutate(Ma=ifelse(diff<=12,runif(1,min_ma,max_ma),NA))
# #gen_Lau_val <- gen_Lau %>% group_by(1:n()) %>% mutate(Ma=round(runif(1,oldest,youngest)))
# 
# ###
# base_SS <- c(485.4,477.7,470,467.3,458.4,453,445.2,443.8)
# name_SS<-c("Tr", "Fl","Dp","Dw","Sa","Ka","Hi")
# 
# gen_Lau_val_SS <- as.data.frame(gen_Lau_val %>% group_by(Ma) %>% mutate(SS=cut(Ma,
#                                                                                breaks = base_SS,
#                                                                                labels = rev(name_SS),
#                                                                                right=FALSE)))
# gen_Lau_val_SS$SS <- factor(gen_Lau_val_SS$SS,
#                             levels = c("Tr", "Fl","Dp","Dw","Sa","Ka","Hi"))
# 
# 
# 
# melt_Lau <- melt(gen_Lau_val_SS, id=c("accepted_name", "SS"), na.rm = TRUE)
# melt_Lau <- na.omit(melt_Lau)
# cast_Lau <- cast(melt_Lau, accepted_name~SS, length)
# 
# #####################################################################################################
# ## Baltica
# gen_Bal_val <- gen_Bal_diff%>% group_by(1:n()) %>% mutate(Ma=ifelse(diff<=12,runif(1,min_ma,max_ma),NA))
# gen_Bal_val_SS <- as.data.frame(gen_Bal_val %>% group_by(Ma) %>% mutate(SS=cut(Ma,
#                                                                                breaks = base_SS,
#                                                                                labels = rev(name_SS),
#                                                                                right=FALSE)))
# gen_Bal_val_SS$SS <- factor(gen_Bal_val_SS$SS,
#                             levels = c("Tr", "Fl","Dp","Dw","Sa","Ka","Hi"))
# 
# melt_Bal <- melt(gen_Bal_val_SS, id=c("accepted_name", "SS"), na.rm = TRUE)
# melt_Bal <- na.omit(melt_Bal)
# cast_Bal <- cast(melt_Bal, accepted_name~SS, length)
# #
# #####################################################################################################
# ## all other continents
# gen_cont_val <- gen_cont_diff%>% group_by(1:n()) %>% mutate(Ma=ifelse(diff<=56,runif(1,min_ma,max_ma),NA))
# gen_cont_val_SS <- as.data.frame(gen_cont_val %>% group_by(Ma) %>% mutate(SS=cut(Ma,
#                                                                                  breaks = base_SS,
#                                                                                  labels = rev(name_SS),
#                                                                                  right=FALSE)))
# gen_cont_val_SS$SS <- factor(gen_cont_val_SS$SS,
#                              levels = c("C1","C2","C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2"))
# 
# melt_cont <- melt(gen_cont_val_SS, id=c("accepted_name", "SS"), na.rm = TRUE)
# melt_cont <- na.omit(melt_cont)
# cast_cont <- cast(melt_cont, accepted_name~SS, length)
# 
# #####################################################################################################
# 
# overlap <- as.factor(Reduce(intersect, list(cast_Bal$accepted_name, cast_Lau$accepted_name)))
# overlap_Bcont <- as.factor(Reduce(intersect, list(cast_Bal$accepted_name, cast_cont$accepted_name)))
# overlap_Lcont <- as.factor(Reduce(intersect, list(cast_Lau$accepted_name, cast_cont$accepted_name)))
# 
# 
# inp_L <- cast_Lau %>% filter(!accepted_name %in% overlap)
# inp_L <- inp_L  %>% filter(!accepted_name %in% overlap_Lcont)
# inp_L <- as.matrix(inp_L[,2:8])
# inp_L <- ifelse(inp_L >=1, 1, 0)
# inp_L <- unite(as.data.frame(inp_L), "ch", c(1:7), sep = "")
# inp_L <- cbind(inp_L, geoplate = "NAC", ";")
# 
# inp_B <- cast_Bal %>% filter(!accepted_name %in% overlap)
# inp_B <- inp_B  %>% filter(!accepted_name %in% overlap_Bcont)
# inp_B <- as.matrix(inp_B[,2:8])
# inp_B <- ifelse(inp_B >=1, 1, 0)
# inp_B <- unite(as.data.frame(inp_B), "ch", c(1:7), sep = "")
# inp_B <- cbind(inp_B, geoplate = "BAL", ";")
# 
# 
# inp <- rbind.data.frame(inp_B,inp_L)
# ## we now exclude all data that has only-0 observation histories, because the Pradel model cannot deal with them
# inp <- inp[inp$ch != "0000000",]
# 
# ###########################################################################
# ## we run models as above
# ## and extract the parameter estimates below
# ## in addition, we also include the duration of time.intervals
# ## otherwise Pradel would assume that time-intervals are of equal length
# 
# proc.Berg_pradsen <- process.data(inp, model= "Pradsen", groups = "geoplate",
#                                   time.intervals = c(7.7,5.2,5.8,7.15,6.6,4.6))#Pradel's survival and "growth"
# plate <- mark(proc.Berg_pradsen, model.parameters = list(Phi=Phi.plate, p=p.plate, Gamma= Gamma.plate))
# 
# 
# proc.Berg_pradlambda <- process.data(inp, model= "Pradlambda", groups = "geoplate",
#                                      time.intervals = c(7.7,5.2,5.8,7.15,6.6,4.6))#Pradel's survival and "growth"
# plateD <-mark(proc.Berg_pradlambda, model.parameters = list(Phi=Phi.plate, p=p.plate, Lambda= L.plate))
# 
# 
# phi_BAL_Ord <- plate$results$real$estimate[1]
# phiu_BAL_Ord <- plate$results$real$ucl[1]
# phil_BAL_Ord <- plate$results$real$lcl[1]
# 
# phi_NAC_Ord <- plate$results$real$estimate[2]
# phiu_NAC_Ord <- plate$results$real$ucl[2]
# phil_NAC_Ord <- plate$results$real$lcl[2]
# 
# p_BAL_Ord <- plate$results$real$estimate[3]
# pu_BAL_Ord <- plate$results$real$ucl[3]
# pl_BAL_Ord <- plate$results$real$lcl[3]
# 
# p_NAC_Ord <- plate$results$real$estimate[4]
# pu_NAC_Ord <- plate$results$real$ucl[4]
# pl_NAC_Ord <- plate$results$real$lcl[4]
# 
# gam_BAL_Ord <- plate$results$real$estimate[5]
# gamu_BAL_Ord <- plate$results$real$ucl[5]
# gaml_BAL_Ord <- plate$results$real$lcl[5]
# 
# gam_NAC_Ord <- plate$results$real$estimate[6]
# gamu_NAC_Ord <- plate$results$real$ucl[6]
# gaml_NAC_Ord <- plate$results$real$lcl[6]
# 
# div_BAL_Ord <- plateD$results$real$estimate[5]
# divu_BAL_Ord <- plateD$results$real$ucl[5]
# divl_BAL_Ord <- plateD$results$real$lcl[5]
# 
# div_NAC_Ord <- plateD$results$real$estimate[6]
# divu_NAC_Ord <- plateD$results$real$ucl[6]
# divl_NAC_Ord <- plateD$results$real$lcl[6]
# 
# 
# BAL_origprob_Ord <- 1-gam_BAL_Ord
# BAL_origCIl_Ord <- 1-gaml_BAL_Ord
# BAL_origCIu_Ord <- 1-gamu_BAL_Ord
# 
# NAC_origprob_Ord <- 1-gam_NAC_Ord
# NAC_origCIl_Ord <- 1-gaml_NAC_Ord
# NAC_origCIu_Ord <- 1-gamu_NAC_Ord
# 
# BAL_extinctprob_Ord <- 1-phi_BAL_Ord
# BAL_extinctCIl_Ord <- 1-phil_BAL_Ord
# BAL_extinctCIu_Ord <- 1-phiu_BAL_Ord
# 
# NAC_extinctprob_Ord <- 1-phi_NAC_Ord
# NAC_extinctCIl_Ord <- 1-phil_NAC_Ord
# NAC_extinctCIu_Ord <- 1-phiu_NAC_Ord

#########################################################################################################
## POPAN MODEL FOR GETTING TOTAL N (superpopulation size) FOR BALTICA AND LAURENTIA DURING THE ORDOVICIAN
# 
# reruns.POP_ORD = list()
# for (i in 1:1){
# 
#   tryCatch({
# genus_Ord <- filter(genus, max_ma < 481.55 & max_ma > 443.8 | min_ma < 481.55 & min_ma >= 443.8)
# gen_Lau <- genus_Ord %>% filter(geoplate %in% c("101"))
# gen_Lau_diff <- cbind(gen_Lau[,c("occurrence_no","accepted_name","max_ma","min_ma","geoplate")],
#                       "diff"=gen_Lau$max_ma-gen_Lau$min_ma)
# 
# 
# 
# gen_Bal <- genus_Ord %>% filter(str_detect(geoplate,'302')|
#                                   str_detect(cc, 'SE')|
#                                   str_detect(cc, 'EE')|
#                                   str_detect(cc, 'FI')|
#                                   str_detect(cc, 'NO')|
#                                   str_detect(state, "St. Petersburg"),
#                                 !str_detect(state,'Spitsbergen|Svalbard'),
#                                 !str_detect(geoplate, '311'))
# 
# gen_Bal_diff <- cbind(gen_Bal[,c("occurrence_no","accepted_name","max_ma","min_ma","geoplate")],
#                       "diff"=gen_Bal$max_ma-gen_Bal$min_ma)
# 
# 
# gen_cont <- genus_Ord %>% filter(!str_detect(geoplate, "101|302"),
#                                  !str_detect(cc, "RU|SE|EE|FI|NO")|
#                                    str_detect(state,'Spitsbergen'))
# 
# gen_cont_diff <- cbind(gen_cont[,c("occurrence_no","accepted_name","max_ma","min_ma","geoplate","ref_pubyr")],
# "diff"=gen_cont$max_ma-gen_cont$min_ma)
# 
# 
# 
#     ###############
#     ## Laurentia
# 
#     gen_Lau_val <- gen_Lau_diff%>% group_by(1:n()) %>% mutate(Ma=ifelse(diff<=12,runif(1,min_ma,max_ma),NA))
#     #gen_Lau_val <- gen_Lau %>% group_by(1:n()) %>% mutate(Ma=round(runif(1,oldest,youngest)))
# 
#     ###
#     base_SS <- c(485.4,477.7,470,467.3,458.4,453,445.2,443.8)
#     name_SS<-c("Tr", "Fl","Dp","Dw","Sa","Ka","Hi")
# 
#     gen_Lau_val_SS <- as.data.frame(gen_Lau_val %>% group_by(Ma) %>% mutate(SS=cut(Ma,
#                                                                                    breaks = base_SS,
#                                                                                    labels = rev(name_SS),
#                                                                                    right=FALSE)))
#     gen_Lau_val_SS$SS <- factor(gen_Lau_val_SS$SS,
#                                 levels = c("Tr", "Fl","Dp","Dw","Sa","Ka","Hi"))
# 
# 
# 
#     melt_Lau <- melt(gen_Lau_val_SS, id=c("accepted_name", "SS"), na.rm = TRUE)
#     melt_Lau <- na.omit(melt_Lau)
#     cast_Lau <- cast(melt_Lau, accepted_name~SS, length)
# 
#     #####################################################################################################
#     ## Baltica
#     gen_Bal_val <- gen_Bal_diff%>% group_by(1:n()) %>% mutate(Ma=ifelse(diff<=12,runif(1,min_ma,max_ma),NA))
#     gen_Bal_val_SS <- as.data.frame(gen_Bal_val %>% group_by(Ma) %>% mutate(SS=cut(Ma,
#                                                                                    breaks = base_SS,
#                                                                                    labels = rev(name_SS),
#                                                                                    right=FALSE)))
#     gen_Bal_val_SS$SS <- factor(gen_Bal_val_SS$SS,
#                                 levels = c("Tr", "Fl","Dp","Dw","Sa","Ka","Hi"))
# 
#     melt_Bal <- melt(gen_Bal_val_SS, id=c("accepted_name", "SS"), na.rm = TRUE)
#     melt_Bal <- na.omit(melt_Bal)
#     cast_Bal <- cast(melt_Bal, accepted_name~SS, length)
# 
#     #####################################################################################################
#     ## all other continents
#     gen_cont_val <- gen_cont_diff%>% group_by(1:n()) %>% mutate(Ma=ifelse(diff<=56,runif(1,min_ma,max_ma),NA))
#     gen_cont_val_SS <- as.data.frame(gen_cont_val %>% group_by(Ma) %>% mutate(SS=cut(Ma,
#                                                                                      breaks = base_SS,
#                                                                                      labels = rev(name_SS),
#                                                                                      right=FALSE)))
#     gen_cont_val_SS$SS <- factor(gen_cont_val_SS$SS,
#                                  levels = c("C1","C2","C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2"))
# 
#     melt_cont <- melt(gen_cont_val_SS, id=c("accepted_name", "SS"), na.rm = TRUE)
#     melt_cont <- na.omit(melt_cont)
#     cast_cont <- cast(melt_cont, accepted_name~SS, length)
# 
#     #####################################################################################################
# 
#     overlap <- as.factor(Reduce(intersect, list(cast_Bal$accepted_name, cast_Lau$accepted_name)))
#     overlap_Bcont <- as.factor(Reduce(intersect, list(cast_Bal$accepted_name, cast_cont$accepted_name)))
#     overlap_Lcont <- as.factor(Reduce(intersect, list(cast_Lau$accepted_name, cast_cont$accepted_name)))
# 
# 
#     inp_L <- cast_Lau %>% filter(!accepted_name %in% overlap)
#     inp_L <- inp_L  %>% filter(!accepted_name %in% overlap_Lcont)
#     inp_L <- as.matrix(inp_L[,2:8])
#     inp_L <- ifelse(inp_L >=1, 1, 0)
#     inp_L <- unite(as.data.frame(inp_L), "ch", c(1:7), sep = "")
#     inp_L <- cbind(inp_L, geoplate = "NAC", ";")
# 
#     inp_B <- cast_Bal %>% filter(!accepted_name %in% overlap)
#     inp_B <- inp_B  %>% filter(!accepted_name %in% overlap_Bcont)
#     inp_B <- as.matrix(inp_B[,2:8])
#     inp_B <- ifelse(inp_B >=1, 1, 0)
#     inp_B <- unite(as.data.frame(inp_B), "ch", c(1:7), sep = "")
#     inp_B <- cbind(inp_B, geoplate = "BAL", ";")
# 
# 
#     inp <- rbind.data.frame(inp_B,inp_L)
#     inp <- inp[inp$ch != "0000000",]
# 
#     proc.POPAN_ORD <- process.data(inp, model= "POPAN", groups = "geoplate",
#                                    time.intervals = c(7.7,5.2,5.8,7.15,6.6,4.6)) ## time.intervals from one sampling occasion (midpoint) to the next midpoint
#     POP.Time.plate_ORD <- mark(proc.POPAN_ORD,
#                                model.parameters=list(Phi=Phi.Time.plate, p=p.Time.plate, pent=pent.Time.plate, N=N.plate))#, options="SIMANNEAL")
# 
#   }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
# 
#   reruns.POP_ORD[[i]]=POP.Time.plate_ORD$results$real
# 
# }
# 
# ##################
# 
# mean_results_POP.ORD <- rowMeans(do.call(cbind, lapply(reruns.POP_ORD, function(x) x$estimate)))
# mean_lcl_POP.ORD <- rowMeans(do.call(cbind, lapply(reruns.POP_ORD, function(x) x$lcl)))
# mean_ucl_POP.ORD <- rowMeans(do.call(cbind, lapply(reruns.POP_ORD, function(x) x$ucl)))
# 
# POP.ORD_results <- data.frame("estimate"=mean_results_POP.ORD, "ucl"=mean_lcl_POP.ORD, "lcl"=mean_ucl_POP.ORD)
# 
# N_ORD_BAL <- POP.ORD_results[27,1]
# N_ORD_BAL_lcI <- POP.ORD_results[27,2]
# N_ORD_BAL_ucI <- POP.ORD_results[27,3]
# 
# N_ORD_LAU <- POP.ORD_results[28,1]
# N_ORD_LAU_lcI <- POP.ORD_results[28,2]
# N_ORD_LAU_ucI <- POP.ORD_results[28,3]
# 
# 
#########################################################################################################
#########################################################################################################
######### plot figure 2 #################################################################################
Stagebase <-c(485.4,477.7,470,467.3,458.4,453,445.2)
Stagemidpoints <- c(481.55,473.85,468.65,462.85,455.7,449.1,444.5)


tpx <- c("Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Sil")
per <- data.frame(c(485.4,477.7,470,467.3,458.4,453,445.2,443.8))
Ord.Stage <- cbind.data.frame(tpx,"Stage"=per)

##
## timescale function from
## http://simpson-carl.github.io/articles/15/timescales.to.base

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


##############################################################################
# par(mfrow=c(2,2), mar = c(3,1,0.1,1), oma = c(0,2,0,0))
par(mfrow=c(2,2), mar = c(0,3,0.1,1), oma = c(4,2,0,0))
##############################################################################
plot(Stagebase-0.2, BAL_Orig_rate[,2], type = "b", 
     pch=19,
     ylim = c(-0.08,0.65),
     xlim = rev(c(444.18,485)),
     axes = F,
     xlab = "",
     ylab = "")

tscales.Ord(0.65, -0.01, -0.08)

lines(Stagebase-0.3, BAL_Orig_rate[,2], type = "b", pch=19, lwd=0.8)
lines(Stagebase[2]-1.5, BAL_Orig_rate[2,2]+0.05, type = "b", pch="B")
arrows(x0=Stagebase-0.3, y0=BAL_Orig_rate_CIl[,2], x1=Stagebase-0.3, y1=BAL_Orig_rate_CIu[,2], length=0.02, lwd = 0.8, angle = 90, code = 3)
lines(Stagebase, NAC_Orig_rate[,2], type = "b", lty = 2, pch=17,lwd=1.2)
lines(Stagebase[1]-1, NAC_Orig_rate[1,2]+0.025, type = "b", pch="L")
arrows(x0=Stagebase, y0=NAC_Orig_rate_CIl[,2], x1=Stagebase, y1=NAC_Orig_rate_CIu[,2], length=0.02, lwd = 1.2, lty=2, angle = 90, code = 3)

legend("topleft", legend="A", bty="n", cex = 1.25)

axis(2, col = 'grey75', line = -0.2, at = seq(0, 0.65, 0.2))

mtext("Origination events per myr", side = 2, line = 2, cex = 1)


##############################################################################
## Extinction Rate
plot(Stagebase-0.2, BAL_Ext_rate[,2], type = "b",
     pch=19,
     ylim = c(-0.08,0.6),
     xlim = rev(c(444.18,485.4)),
     axes = F,
     xlab = "",
     ylab = "")

tscales.Ord(0.6, -0.01, -0.08)

lines(Stagebase-0.3, BAL_Ext_rate[,2], type = "b", pch=19, lwd=0.8)
arrows(x0=Stagebase-0.3, y0=BAL_Ext_rate_CIl[,2], x1=Stagebase-0.3, y1=BAL_Ext_rate_CIu[,2], length=0.02, lwd = 0.8, angle = 90, code = 3)
lines(Stagebase, NAC_Ext_rate[,2], type = "b", lty = 2, lwd=1.2, pch=17)
arrows(x0=Stagebase, y0=NAC_Ext_rate_CIl[,2], x1=Stagebase, y1=NAC_Ext_rate_CIu[,2], length=0.02, lwd = 1.2,lty=2,  angle = 90, code = 3)

legend("topleft", legend="B", bty="n", cex = 1.25)

axis(2, col = 'grey75', line = -0.2, at = seq(0, 0.6, 0.2))

mtext("Extinction events per myr", side = 2, line = 2, cex = 1)


##############################################################################
## net diversification 
plot(Stagebase-0.2, div_BAL[,2]-1, type = "b", 
     ylim = c(-1.8,7),
     xlim = rev(c(444.18,485.4)),
     axes = F,
     xlab = "",
     ylab = "")

tscales.Ord(7, -1.1, -1.8)
abline(h = 0, col="darkgrey")

lines(Stagebase-0.2, div_BAL[,2]-1, type = "b", lwd=0.8, pch=19)
arrows(x0=Stagebase-0.2, y0=divl_BAL[,2]-1, x1=Stagebase-0.2, y1=divu_BAL[,2]-1, length=0.02, lwd = 0.8, angle = 90, code = 3)
lines(Stagebase, div_NAC[,2]-1, type = "b", lty = 2, pch=17, lwd=0.8)
arrows(x0=Stagebase, y0=divl_NAC[,2]-1, x1=Stagebase, y1=divu_NAC[,2]-1, length=0.02, lwd = 1.2, lty=2,  angle = 90, code = 3)

legend("topleft", legend="C", bty="n", cex = 1.25)

axis(1, col = 'grey75', line = 0.15, at = seq(445,485,10))
axis(2, col = 'grey75', line = -0.2, at = seq(-1, 7, 1))

mtext("Age (Ma)", side = 1, line = 2.5, cex = 1)
mtext("Net diversification rate", side = 2, line = 2, cex = 1)

##############################################################################
## Sampling rate
plot(Stagemidpoints-0.1, rate_p_BAL[,2], type = 'b', 
     xlim = rev(c(444.18,485.4)),
     ylim = c(-0.1,0.8),
     axes = F,
     xlab = "", ylab = "")

tscales.Ord(0.8, -0.01, -0.1)

lines(Stagemidpoints-0.1, rate_p_BAL[,2], type = "b", pch=19,lwd=0.8)
arrows(x0=Stagemidpoints-0.1, y0=ratel_p_BAL[,2], x1=Stagemidpoints-0.1, y1=rateu_p_BAL[,2], length=0.02, lwd = 0.8, angle = 90, code = 3)
lines(Stagemidpoints+0.1, rate_p_LAU[,2], type = "b" , lty = 2, lwd=1.2, pch=17)
arrows(x0=Stagemidpoints+0.1, y0=ratel_p_LAU[,2], x1=Stagemidpoints+0.1, y1=rateu_p_LAU[,2], length=0.02, lwd = 1.2, lty=2, angle = 90, code = 3)

legend("topleft", legend="D", bty="n", cex = 1.25)

axis(1, col = 'grey75', line = 0.15, at = seq(445,485,10))
axis(2, col = 'grey75', line = -0.15, at = seq(0, 0.8, 0.2))

mtext("Age (Ma)", side = 1, line = 2.5, cex = 1)
mtext("Sampling events per myr", side = 2, line = 2, cex = 1)




################################################################################################
##### Figure S7 ## replicate of 100 runs #######################################################


par(mfrow=c(2,2), mar = c(1,2,0.1,1), oma = c(2,2,0,0))

plot(Stagebase-0.2, mean.origrate.B, type = "b", 
     pch=19,
     ylim = c(-0.12,1),
     xlim = rev(c(444.18,485)),
     axes = F,
     xlab = "",
     ylab = "")

tscales.Ord(1, -0.01, -0.12)

for (i in 1:100){
  lines(Stagebase-0.2, BAL_Orig_rate_CIl[,i], col = "grey70")}
for (i in 1:100){
  lines(Stagebase-0.2, BAL_Orig_rate_CIu[,i], col = "grey70")
}  


for (i in 1:100){
  lines(Stagebase, NAC_Orig_rate_CIl[,i], col = "grey55")}
for (i in 1:100){
  lines(Stagebase, NAC_Orig_rate_CIu[,i], col = "grey55")}



lines(Stagebase-0.3, mean.origrate.B, type = "b", pch=19, lwd=0.8)
lines(Stagebase, mean.origrate.L, type = "b", lty = 2, pch=17,lwd=1.2)

lines(Stagebase[1]-1, mean.origrate.B[1]-0.02, type = "b", pch="B")
lines(Stagebase[1]-1, mean.origrate.L[1]-0.02, type = "b", pch="L")

legend("topleft", legend="A", bty="n", cex = 1.25)

axis(2, col = 'grey75', line = -0.2, at = seq(0, 1, 0.2))

mtext("Origination events per myr", side = 2, line = 2, cex = 0.75)


##############################################################################
## Extinction Rate
plot(Stagebase-0.2, mean.extrate.B, type = "b",
     pch=19,
     ylim = c(-0.12,1),
     xlim = rev(c(444.18,485.4)),
     axes = F,
     xlab = "",
     ylab = "")

tscales.Ord(1, -0.01, -0.12)

for (i in 1:100){
  lines(Stagebase-0.2, BAL_Ext_rate_CIl[,i], col = "grey70")}
for (i in 1:100){
  lines(Stagebase-0.2, BAL_Ext_rate_CIu[,i], col = "grey70")
}  


for (i in 1:100){
  lines(Stagebase, NAC_Ext_rate_CIl[,i], col = "grey55")}
for (i in 1:100){
  lines(Stagebase, NAC_Ext_rate_CIu[,i], col = "grey55")}

lines(Stagebase-0.3, mean.extrate.B, type = "b", pch=19, lwd=0.8)
lines(Stagebase, mean.extrate.L, type = "b", lty = 2, lwd=1.2, pch=17)


legend("topleft", legend="B", bty="n", cex = 1.25)

axis(2, col = 'grey75', line = -0.2, at = seq(0, 1, 0.2))

mtext("Extinction events per myr", side = 2, line = 2, cex = 0.75)


##############################################################################
## net diversification 
plot(Stagebase-0.2, mean.div.B-1, type = "b", 
     ylim = c(-2,10),
     xlim = rev(c(444.18,485.4)),
     axes = F,
     xlab = "",
     ylab = "")

tscales.Ord(10, -1.1, -2)
abline(h = 0, col="darkgrey")


for (i in 1:100){
  lines(Stagebase-0.2, divl_BAL[,i]-1, col = "grey70")}
for (i in 1:100){
  lines(Stagebase-0.2, divu_BAL[,i]-1, col = "grey70")
}  


for (i in 1:100){
  lines(Stagebase, divl_NAC[,i]-1, col = "grey55")}
for (i in 1:100){
  lines(Stagebase, divu_NAC[,i]-1, col = "grey55")}

lines(Stagebase-0.2, mean.div.B-1, type = "b", lwd=0.8, pch=19)
lines(Stagebase, mean.div.L-1, type = "b", lty = 2, pch=17, lwd=0.8)

legend("topleft", legend="C", bty="n", cex = 1.25)

axis(1, col = 'grey75', line = 0.15, at = seq(445,485,10))
axis(2, col = 'grey75', line = -0.2, at = seq(-2, 10, 2))

mtext("Age (Ma)", side = 1, line = 2, cex = 0.75)
mtext("Net diversification rate", side = 2, line = 2, cex = 0.75)

##############################################################################
## Sampling rate
plot(Stagemidpoints-0.1, mean.prate.B, type = 'b', 
     xlim = rev(c(444.18,485.4)),
     ylim = c(-0.15,1),
     axes = F,
     xlab = "", ylab = "")

tscales.Ord(1, -0.01, -0.15)

for (i in 1:100){
  lines(Stagemidpoints-0.2, ratel_p_BAL[,i], col = "grey70")}
for (i in 1:100){
  lines(Stagemidpoints-0.2, rateu_p_BAL[,i], col = "grey70")
}  


for (i in 1:100){
  lines(Stagemidpoints, ratel_p_LAU[,i], col = "grey55")}
for (i in 1:100){
  lines(Stagemidpoints, rateu_p_LAU[,i], col = "grey55")}


lines(Stagemidpoints-0.1, mean.prate.B, type = "b", pch=19,lwd=0.8)
lines(Stagemidpoints+0.1, mean.prate.L, type = "b" , lty = 2, lwd=1.2, pch=17)


legend("topleft", legend="D", bty="n", cex = 1.25)

axis(1, col = 'grey75', line = 0.15, at = seq(445,485,10))
axis(2, col = 'grey75', line = -0.15, at = seq(0, 1, 0.2))

mtext("Age (Ma)", side = 1, line = 2, cex = 0.75)
mtext("Sampling events per myr", side = 2, line = 2, cex = 0.75)



#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
## POPAN model
## processing of the data to get the input file  as before

reruns.popan=list()

for (i in 1:100){
  
  gen_Lau <- genus %>% filter(geoplate %in% c("101"))
  gen_Lau_diff <- cbind(gen_Lau[,c("occurrence_no","accepted_name","max_ma","min_ma","geoplate","ref_pubyr")],
                        "diff"=gen_Lau$max_ma-gen_Lau$min_ma)
  pubyear_L <- unique(gen_Lau[,c("ref_pubyr","reference_no")])
  
  
  gen_Bal <- genus %>% filter(str_detect(geoplate,'302')|
                                str_detect(cc, 'SE')|
                                str_detect(cc, 'EE')|
                                str_detect(cc, 'FI')|
                                str_detect(cc, 'NO')|
                                str_detect(state, "St. Petersburg"),
                              !str_detect(state,'Spitsbergen|Svalbard'),
                              !str_detect(geoplate, '311'))
  
  gen_Bal_diff <- cbind(gen_Bal[,c("occurrence_no","accepted_name","max_ma","min_ma","geoplate","ref_pubyr")],
                        "diff"=gen_Bal$max_ma-gen_Bal$min_ma)
  pubyear_B <- unique(gen_Bal[,c("ref_pubyr","reference_no")])
  
  
  #####################################################################################################
  ## Laurentia
  
  gen_Lau_val <- gen_Lau_diff%>% group_by(1:n()) %>% mutate(Ma=ifelse(diff<=12,runif(1,min_ma,max_ma),NA))
  base_SS <- c(509,497,485.4,477.7,470,467.3,458.4,453,445.2,443.8,440.8,438.5)
  name_SS<-c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2")
  gen_Lau_val_SS <- as.data.frame(gen_Lau_val %>% group_by(Ma) %>% mutate(SS=cut(Ma, 
                                                                                 breaks = base_SS, 
                                                                                 labels = rev(name_SS),
                                                                                 right=FALSE)))
  gen_Lau_val_SS$SS <- factor(gen_Lau_val_SS$SS,
                              levels = c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2"))
  melt_Lau <- melt(gen_Lau_val_SS, id=c("accepted_name", "SS"), na.rm = TRUE)
  melt_Lau <- na.omit(melt_Lau)
  cast_Lau <- cast(melt_Lau, accepted_name~SS, length)
  
  #####################################################################################################
  ## Baltica
  
  gen_Bal_val <- gen_Bal_diff%>% group_by(1:n()) %>% mutate(Ma=ifelse(diff<=12,runif(1,min_ma,max_ma),NA))
  gen_Bal_val_SS <- as.data.frame(gen_Bal_val %>% group_by(Ma) %>% mutate(SS=cut(Ma, 
                                                                                 breaks = base_SS, 
                                                                                 labels = rev(name_SS),
                                                                                 right=FALSE)))
  gen_Bal_val_SS$SS <- factor(gen_Bal_val_SS$SS,
                              levels = c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2"))
  
  melt_Bal <- melt(gen_Bal_val_SS, id=c("accepted_name", "SS"), na.rm = TRUE)
  melt_Bal <- na.omit(melt_Bal)
  cast_Bal <- cast(melt_Bal, accepted_name~SS, length)
  
  
  
  #####################################################################################################
  ## all but LAU and BAL
  gen_cont_val <- gen_cont_diff%>% group_by(1:n()) %>% mutate(Ma=ifelse(diff<=56,runif(1,min_ma,max_ma),NA))
  gen_cont_val_SS <- as.data.frame(gen_cont_val %>% group_by(Ma) %>% mutate(SS=cut(Ma,
                                                                                   breaks = base_SS, 
                                                                                   labels = rev(name_SS),
                                                                                   right=FALSE)))
  gen_cont_val_SS$SS <- factor(gen_cont_val_SS$SS,
                               levels = c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2"))
  
  melt_cont <- melt(gen_cont_val_SS, id=c("accepted_name", "SS"), na.rm = TRUE)
  melt_cont <- na.omit(melt_cont)
  cast_cont <- cast(melt_cont, accepted_name~SS, length)
  
  
  overlap <- as.factor(Reduce(intersect, list(cast_Bal$accepted_name, cast_Lau$accepted_name)))
  overlap_Bcont <- as.factor(Reduce(intersect, list(cast_Bal$accepted_name, cast_cont$accepted_name)))
  overlap_Lcont <- as.factor(Reduce(intersect, list(cast_Lau$accepted_name, cast_cont$accepted_name)))
  
  
  
  inp_L <- cast_Lau %>% filter(!accepted_name %in% overlap)
  inp_L <- inp_L  %>% filter(!accepted_name %in% overlap_Lcont)
  inp_Lmat <- data.frame(inp_L[,1],ifelse(inp_L[,2:12]>=1,1,0))
  
  inp_L <- as.matrix(inp_L[,2:12])
  inp_L <- ifelse(inp_L >=1, 1, 0)
  inp_L <- unite(as.data.frame(inp_L), "ch", c(1:11), sep = "")
  inp_L <- cbind(inp_L, geoplate = "NAC", ";")
  
  inp_B <- cast_Bal %>% filter(!accepted_name %in% overlap)
  inp_B <- inp_B  %>% filter(!accepted_name %in% overlap_Bcont)
  inp_Bmat <- data.frame(inp_B[,1],ifelse(inp_B[,2:12]>=1,1,0))
  
  inp_B <- as.matrix(inp_B[,2:12])
  inp_B <- ifelse(inp_B >=1, 1, 0)
  inp_B <- unite(as.data.frame(inp_B), "ch", c(1:11), sep = "")
  inp_B <- cbind(inp_B, geoplate = "BAL", ";")
  
  
  inp <- rbind.data.frame(inp_B,inp_L)
  
  proc.POPAN <- process.data(inp, model= "POPAN", groups = "geoplate",
                             time.intervals = c(11.8,9.65,7.7,5.2,5.8,7.15,6.6,4.6,2.2,2.65)) ## time.intervals from one sampling occasion (midpoint) to the next midpoint
  POP.Time.plate <- mark(proc.POPAN,
                         model.parameters=list(Phi=Phi.Time.plate, p=p.Time.plate, pent=pent.Time.plate, N=N.plate))#, options="SIMANNEAL")
  
  
  reruns.popan[[i]]=POP.Time.plate$results$derived$`N Population Size`
  
}

meanN_results <- rowMeans(do.call(cbind, lapply(reruns.popan, function(x) x$estimate)))

estimates.PC.N <- do.call(cbind, lapply(reruns.popan, function(x) x$estimate))
lcl.PC.N <- do.call(cbind, lapply(reruns.popan, function(x) x$lcl))
ucl.PC.N <- do.call(cbind, lapply(reruns.popan, function(x) x$ucl))


N_B <- estimates.PC.N[1:11,]
N_B_lcI <- lcl.PC.N[1:11,]
N_B_ucI <- ucl.PC.N[1:11,]

N_L <- estimates.PC.N[12:22,]
N_L_lcI <- lcl.PC.N[12:22,]
N_L_ucI <- ucl.PC.N[12:22,]

mean.N.B <- meanN_results[1:11]
mean.N.L <- meanN_results[12:22]


#########################################################################################################
###### Figure 3 #########################################################################################

par(mfrow=c(1,1), mai = c(0.7,0.7,0.2,0.2))

plot(Stagemidpoints, N_B[3:9,2],
     type = "b", 
     ylim = c(-35,400),
     xlim = rev(c(444.18,485.4)),
     axes = F,
     xlab = "",
     ylab = "")


tscales.Ord(400, -5, -35)

lines(Stagemidpoints-0.2, N_B[3:9,2], type = "b", pch=19, lwd=0.8)
lines(Stagemidpoints[1]-1, N_B[3,2]-25, type = "p", pch="B")
arrows(x0=Stagemidpoints-0.2, y0=N_B_lcI[3:9,2], x1=Stagemidpoints-0.2, y1=N_B_ucI[3:9,2],
       length=0.02, lwd = 0.8, angle = 90, code = 3)
lines(Stagemidpoints+0.2, N_L[3:9,2], type = "b", pch=17, lwd=1.2, lty=2)
lines(Stagemidpoints[1]+1, N_L[3,2]+20, type = "p", pch="L")
arrows(x0=Stagemidpoints+0.2, y0=N_L_lcI[3:9,2], x1=Stagemidpoints+0.2, y1=N_L_ucI[3:9,2],
       length=0.02, lwd = 1.2, angle = 90, lty=2, code = 3)


axis(1, col = 'grey75', line = 0.5, at = seq(445,485,10) )
axis(2, col = 'grey75', line = -0.2, at = seq(0, 400, 100))

mtext("Age (Ma)", side = 1, line = 2.5)
mtext("Genus richness", side = 2, line = 2)

######################################################################################
############ Figure S8 ###############################################################

par(mfrow=c(1,1), mai = c(0.7,0.7,0.2,0.2))

plot(Stagemidpoints, mean.N.B[3:9],
     type = "b", 
     ylim = c(-35,400),
     xlim = rev(c(444.18,485.4)),
     axes = F,
     xlab = "",
     ylab = "")
tscales.Ord(400, -5, -35)


for (i in 1:100){
  lines(Stagemidpoints, N_B_lcI[3:9,i], col = "grey55")}
for (i in 1:100){
  lines(Stagemidpoints, N_B_ucI[3:9,i], col = "grey55")}

for (i in 1:100){
  lines(Stagemidpoints, N_L_lcI[3:9,i], col = "grey35")}
for (i in 1:100){
  lines(Stagemidpoints, N_L_ucI[3:9,i], col = "grey35")}


lines(Stagemidpoints,mean.N.B[3:9],type="b", lwd=1.5, pch=19)
lines(Stagemidpoints,mean.N.L[3:9],type="b", lwd=1.5, pch=17, lty=2)


axis(1, col = 'grey75', line = 0.5, at = seq(445,485,10) )
axis(2, col = 'grey75', line = -0.2, at = seq(0, 450, 100))

mtext("Age (Ma)", side = 1, line = 2.5)
mtext("Genus richness", side = 2, line = 2)

######################################################################################
######################################################################################
## Phyla endemic per continent
## Table S7 and S8
###########################################################################
genusB <- filter(gen_Bal, grepl("genus", accepted_rank))
genus_OrdB <- filter(genusB, max_ma < 481.55 & max_ma > 443.8 | min_ma < 481.55 & min_ma > 443.8)
genus_phylaB <- unique(genus_OrdB[,c("accepted_name","phylum")])

###########################################################################
genusL <- filter(gen_Lau, grepl("genus", accepted_rank))
genus_OrdL <- filter(genusL, max_ma < 481.55 & max_ma > 443.8 | min_ma < 481.55 & min_ma > 443.8)
genus_phylaL <- unique(genus_OrdL[,c("accepted_name","phylum")])

###########################################################################
genuscont <- filter(gen_cont, grepl("genus",accepted_rank))
genus_Ord_cont <- filter(genuscont, max_ma < 481.55 & max_ma > 443.8 | min_ma < 481.55 & min_ma > 443.8)
genus_phyla_cont <- unique(genus_Ord_cont[,c("accepted_name","phylum")])


## to exclude any overlapping genera:
overlap_BL <- Reduce(intersect, list(genus_phylaB$accepted_name,genus_phylaL$accepted_name))
overlap_BC <- Reduce(intersect, list(genus_phylaB$accepted_name,genus_phyla_cont$accepted_name))
overlap_LC <- Reduce(intersect, list(genus_phylaL$accepted_name,genus_phyla_cont$accepted_name))

overlap_all <- c(overlap_BL,overlap_BC,overlap_LC)
overlap_all <- unique(overlap_all)

end_Lau <- genus_phylaL %>% filter(!accepted_name %in% overlap_all)
end_Bal <- genus_phylaB %>% filter(!accepted_name %in% overlap_all)


nend_L <- aggregate(list(end_Lau$accepted_name),
                    list(end_Lau$phylum),
                    length)


nend_B <- aggregate(list(end_Bal$accepted_name),
                    list(end_Bal$phylum),
                    length)

Bal_phy_N <- colSums(nend_B[2], na.rm = T) ## total amount of "endemic" genera during the Ordovician on Baltica
Lau_phy_N <- colSums(nend_L[2], na.rm = T) ## total amount of "endemic" genera during the Ordovician on Laurentia

nend_B[,2] <- nend_B[,2]*100/colSums(nend_B[2], na.rm = T)
nend_B <- nend_B %>% filter(!Group.1 %in% c("","Problematica"))

nend_L[,2] <- nend_L[,2]*100/colSums(nend_L[2], na.rm = T)
nend_L <- nend_L %>% filter(!Group.1 %in% c("","Problematica"))

colnames(nend_B)[1] <- "Phylum"
colnames(nend_B)[2] <- "percent"
colnames(nend_L)[1] <- "Phylum"
colnames(nend_L)[2] <- "percent"


end_B <- arrange(nend_B,desc(nend_B[,2]))
end_B$group <- ifelse(end_B[,2]>1.5, 1, 2)
end_B1 <- end_B[end_B$group==1,]

end_B1$Phylum <- factor(end_B1$Phylum,
                        levels = c(levels(end_B1$Phylum),"Others"))


end_B1[nrow(end_B1)+1,1] <- c("Others")
end_B1[nrow(end_B1),2] <- c(100-colSums(end_B1[2], na.rm = T))
endemics_B <- end_B1[,1:2]

end_L <- arrange(nend_L,desc(nend_L[,2]))
end_L$group <- ifelse(end_L[,2]>1.5, 1, 2)
end_L1 <- end_L[end_L$group==1,]

end_L1$Phylum <- factor(end_L1$Phylum,
                        levels = c(levels(end_L1$Phylum),"Others"))

end_L1[nrow(end_L1)+1,1] <- c("Others")
end_L1[nrow(end_L1),2] <- c(100-colSums(end_L1[2], na.rm = T))
endemics_L <- end_L1[,1:2]

df_endemics <- rbind.data.frame(cbind.data.frame(endemics_B,"id"="Baltica"), cbind.data.frame(endemics_L,"id"="Laurentia"))


ggplot(df_endemics, aes(x="", y=percent, fill=Phylum)) + 
  geom_bar(stat="identity", colour = "white") +
  facet_wrap(~id) + 
  theme_minimal() +
  scale_fill_grey(start = 0, end = 0.9) +
  coord_polar("y", start=0)


########################################################################################################
########################################################################################################
## model comparison
## creating input files for Laurentia/Baltica data as above
## Table S2

rm(list=ls())
genus <- read.csv("PBDB_Ord_1.csv", sep = ",", header=T)
genus <- filter(genus, grepl("genus", accepted_rank))

########
## Pradel seniority model specificaitons
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

##############
## POPAN model specifications
Phi.Time.plate <- list(formula=~time+geoplate+time*geoplate)
p.Time.plate <- list(formula=~time+geoplate+time*geoplate)
pent.Time.plate <- list(formula=~time+geoplate+time*geoplate)
N.time <- list(formula=~time)
N.const <- list(formula=~1)
N.plate <- list(formula=~geoplate)

#####################################################################################################
## Laurentia


gen_Lau <- genus %>% filter(geoplate %in% c("101"))
gen_Lau_diff <- cbind(gen_Lau[,c("occurrence_no","accepted_name","max_ma","min_ma","geoplate","ref_pubyr")],
                      "diff"=gen_Lau$max_ma-gen_Lau$min_ma)



##  Baltica
gen_Bal <- genus %>% filter(str_detect(geoplate,'302')|
                              str_detect(cc, 'SE')|
                              str_detect(cc, 'EE')|
                              str_detect(cc, 'FI')|
                              str_detect(cc, 'NO')|
                              str_detect(state, "St. Petersburg"),
                            !str_detect(state,'Spitsbergen|Svalbard'),
                            !str_detect(geoplate, '311'))

gen_Bal_diff <- cbind(gen_Bal[,c("occurrence_no","accepted_name","max_ma","min_ma","geoplate","ref_pubyr")],
                      "diff"=gen_Bal$max_ma-gen_Bal$min_ma)


## extract all data EXCEPT from Baltica or Laurentia
gen_cont <- genus %>% filter(!str_detect(geoplate, "101|302"),
                             !str_detect(cc, "SE|EE|FI|NO"),
                             !str_detect(state, "St. Petersburg")|
                               str_detect(state,'Spitsbergen|Svalbard')|
                               str_detect(geoplate, "311"))



gen_cont_diff <- cbind(gen_cont[,c("occurrence_no","accepted_name","max_ma","min_ma","geoplate","ref_pubyr")],
                       "diff"=gen_cont$max_ma-gen_cont$min_ma)


gen_Lau_val <- gen_Lau_diff%>% group_by(1:n()) %>% mutate(Ma=ifelse(diff<=12,runif(1,min_ma,max_ma),NA))
base_SS <- c(509,497,485.4,477.7,470,467.3,458.4,453,445.2,443.8,440.8,438.5)
name_SS<-c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2")
gen_Lau_val_SS <- as.data.frame(gen_Lau_val %>% group_by(Ma) %>% mutate(SS=cut(Ma,
                                                                               breaks = base_SS,
                                                                               labels = rev(name_SS),
                                                                               right=FALSE)))
gen_Lau_val_SS$SS <- factor(gen_Lau_val_SS$SS,
                            levels = c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2"))
melt_Lau <- melt(gen_Lau_val_SS, id=c("accepted_name", "SS"), na.rm = TRUE)
melt_Lau <- na.omit(melt_Lau)
cast_Lau <- cast(melt_Lau, accepted_name~SS, length)
#####################################################################################################
## Baltic Shield

gen_Bal_val <- gen_Bal_diff%>% group_by(1:n()) %>% mutate(Ma=ifelse(diff<=12,runif(1,min_ma,max_ma),NA))
gen_Bal_val_SS <- as.data.frame(gen_Bal_val %>% group_by(Ma) %>% mutate(SS=cut(Ma,
                                                                               breaks = base_SS,
                                                                               labels = rev(name_SS),
                                                                               right=FALSE)))
gen_Bal_val_SS$SS <- factor(gen_Bal_val_SS$SS,
                            levels = c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2"))

melt_Bal <- melt(gen_Bal_val_SS, id=c("accepted_name", "SS"), na.rm = TRUE)
melt_Bal <- na.omit(melt_Bal)
cast_Bal <- cast(melt_Bal, accepted_name~SS, length)


#####################################################################################################
#####################################################################################################
gen_cont_val <- gen_cont_diff%>% group_by(1:n()) %>% mutate(Ma=ifelse(diff<=56,runif(1,min_ma,max_ma),NA))
gen_cont_val_SS <- as.data.frame(gen_cont_val %>% group_by(Ma) %>% mutate(SS=cut(Ma,
                                                                                 breaks = base_SS,
                                                                                 labels = rev(name_SS),
                                                                                 right=FALSE)))
gen_cont_val_SS$SS <- factor(gen_cont_val_SS$SS,
                             levels = c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2"))

melt_cont <- melt(gen_cont_val_SS, id=c("accepted_name", "SS"), na.rm = TRUE)
melt_cont <- na.omit(melt_cont)
cast_cont <- cast(melt_cont, accepted_name~SS, length)



overlap <- as.factor(Reduce(intersect, list(cast_Bal$accepted_name, cast_Lau$accepted_name)))
overlap_Bcont <- as.factor(Reduce(intersect, list(cast_Bal$accepted_name, cast_cont$accepted_name)))
overlap_Lcont <- as.factor(Reduce(intersect, list(cast_Lau$accepted_name, cast_cont$accepted_name)))


inp_L <- cast_Lau %>% filter(!accepted_name %in% overlap)
inp_L <- inp_L  %>% filter(!accepted_name %in% overlap_Lcont)
inp_Lmat <- data.frame(inp_L[,1],ifelse(inp_L[,2:12]>=1,1,0))

inp_L <- as.matrix(inp_L[,2:12])
inp_L <- ifelse(inp_L >=1, 1, 0)
inp_L <- unite(as.data.frame(inp_L), "ch", c(1:11), sep = "")
inp_L <- cbind(inp_L, geoplate = "NAC", ";")

inp_B <- cast_Bal %>% filter(!accepted_name %in% overlap)
inp_B <- inp_B  %>% filter(!accepted_name %in% overlap_Bcont)
inp_Bmat <- data.frame(inp_B[,1],ifelse(inp_B[,2:12]>=1,1,0))

inp_B <- as.matrix(inp_B[,2:12])
inp_B <- ifelse(inp_B >=1, 1, 0)
inp_B <- unite(as.data.frame(inp_B), "ch", c(1:11), sep = "")
inp_B <- cbind(inp_B, geoplate = "BAL", ";")


inp <- rbind.data.frame(inp_B, inp_L)
inp$geoplate <- factor(inp$geoplate,
                       levels = c("BAL", "NAC"))


##########################################################################################
##########################################################################################
##  model comparison
#
proc.pradsen  <- process.data(inp, model= "Pradsen", groups = "geoplate")
# proc.pradlambda <- process.data(inp, model= "Pradlambda", groups = "geoplate")#Pradel's survival and "growth"

time <- mark(proc.pradsen,
             model.parameters=list(Phi=Phi.time, p=p.time, Gamma=Gamma.time))
time.const <- mark(proc.pradsen,
                   model.parameters = list(Phi=Phi.const, p=p.const, Gamma=Gamma.const))
phi.plate.p.const <- mark(proc.pradsen,
                          model.parameters = list(Phi=Phi.plate, p=p.const, Gamma=Gamma.const))
phi.plate.p.plate <- mark(proc.pradsen,
                          model.parameters = list(Phi=Phi.plate, p=p.plate, Gamma=Gamma.const))
phi.plate.p.time <- mark(proc.pradsen,
                         model.parameters = list(Phi=Phi.plate, p=p.time, Gamma=Gamma.const))
phi.plate.p.time.gamma.time <- mark(proc.pradsen,
                                    model.parameters = list(Phi=Phi.plate, p=p.time, Gamma=Gamma.time))
phi.plate.p.plate.gamma.time <- mark(proc.pradsen,
                                     model.parameters = list(Phi=Phi.plate, p=p.plate, Gamma=Gamma.time))

phi.tp.p.tp.gamma.t <- mark(proc.pradsen,
                            model.parameters = list(Phi=Phi.time.plate, p=p.time.plate, Gamma=Gamma.time))

phi.tp.p.t.gamma.tp <- mark(proc.pradsen,
                            model.parameters = list(Phi=Phi.time.plate, p=p.time, Gamma=Gamma.time.plate))

phi.t.p.tp.g.tp <- mark(proc.pradsen,
                        model.parameters = list(Phi=Phi.time, p=p.time.plate, Gamma=Gamma.time.plate))

time.plate <- mark(proc.pradsen,
                   model.parameters = list(Phi=Phi.time.plate, p=p.time.plate, Gamma=Gamma.time.plate))

phi.Tp.p.tp.gamma.Tp <- mark(proc.pradsen,
                             model.parameters = list(Phi=Phi.Time.plate, p=p.time.plate, Gamma=L.Time.plate))

phi.tp.p.Tp.gamma.tp <- mark(proc.pradsen,
                             model.parameters = list(Phi=Phi.time.plate, p=p.Time.plate, Gamma=L.time.plate))

Time.plate <- mark(proc.pradsen,
                   model.parameters = list(Phi=Phi.Time.plate, p=p.Time.plate, Gamma=L.Time.plate))

plate <- mark(proc.pradsen, model.parameters = list(Phi=Phi.plate, p=p.plate, Gamma= Gamma.plate))

## compare all models in following table
comparison <- collect.models(type = "Pradsen")
comparison

