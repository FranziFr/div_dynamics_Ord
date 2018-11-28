library(reshape2)
library(reshape)
library(dplyr)
library(tidyr)
library(RMark)
library(ggplot2)



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
Phi.time.env <- list(formula=~time+environment)
p.time.env <- list(formula=~time+environment)
Gamma.time.env <- list(formula=~time+environment)
L.time.env <- list(formula=~time+environment)

Phi.Time.env <- list(formula=~time+environment+time*environment)
p.Time.env <- list(formula=~time+environment+time*environment)
Gamma.Time.env <- list(formula=~time+environment+time*environment)
L.Time.env <- list(formula=~time+environment+time*environment)

Phi.const <- list(formula=~1)
p.const <- list(formula=~1)
Gamma.const <- list(formula=~1)
L.const <- list(formula=~1)

Phi.env <- list(formula=~environment)
p.env <- list(formula=~environment)
Gamma.env <- list(formula=~environment)
L.env <- list(formula=~environment)

Phi.Time.env <- list(formula=~time+environment+time*environment)
p.Time.env <- list(formula=~time+environment+time*environment)

#####################################################################################################
## onshore
## offshore

gen_on <- genus %>% filter(environment %in% c("foreshore","coastal indet.",
                                              "marginal marine indet.",
                                              "carbonate indet.","peritidal",
                                              "sand shoal", "shoreface"))



gen_on_diff <- cbind(gen_on[,c("occurrence_no","accepted_name","phylum","max_ma","min_ma","geoplate")],
                     "diff"=gen_on$max_ma-gen_on$min_ma)


gen_off <- genus %>% filter(environment %in% c("offshore indet.","offshore","offshore ramp","offshore shelf"))

gen_off_diff <- cbind(gen_off[,c("occurrence_no","accepted_name","phylum","max_ma","min_ma","geoplate")],
                      "diff"=gen_off$max_ma-gen_off$min_ma)


gen_else <- genus %>% filter(!environment %in% c("offshore indet.","offshore","offshore ramp","offshore shelf",
                                                 "foreshore","coastal indet.",
                                                 "marginal marine indet.",
                                                 "carbonate indet.","peritidal",
                                                 "sand shoal", "shoreface"))


gen_else_diff <- cbind(gen_else[,c("occurrence_no","accepted_name","phylum","max_ma","min_ma","geoplate")],
                       "diff"=gen_else$max_ma-gen_else$min_ma)

##############################################################################
## bases and names of stages for later assignment to Stages 
base_SS <- c(509,497,485.4,477.7,470,467.3,458.4,453,445.2,443.8,440.8,438.5)
name_SS<-c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2")


#########################################################################################################
## Pradel model
## Pradsen parametrization


reruns.oo=list()
reruns.lambda.oo=list()


for (i in 1:100){
  
  tryCatch({
    
    ## creation of presence-absence matrices for the different environments,
    ## as it was done for the Palaeocontinents
    ## commands are as described for Laurentia
    
    #####
    ## onshore
    gen_on_val <- gen_on_diff%>% group_by(1:n()) %>% mutate(Ma=ifelse(diff<=12,runif(1,min_ma,max_ma),NA))
    gen_on_val_SS <- as.data.frame(gen_on_val %>% group_by(Ma) %>% mutate(SS=cut(Ma, 
                                                                                 breaks = base_SS, 
                                                                                 labels = rev(name_SS),
                                                                                 right=FALSE)))
    gen_on_val_SS$SS <- factor(gen_on_val_SS$SS,
                               levels = c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2"))
    melt_on <- melt(gen_on_val_SS, id=c("accepted_name", "SS"), na.rm = TRUE)
    melt_on <- na.omit(melt_on)
    cast_on <- cast(melt_on, accepted_name~SS, length)
    
    
    #####
    ## offshore
    gen_off_val <- gen_off_diff%>% group_by(1:n()) %>% mutate(Ma=ifelse(diff<=12,runif(1,min_ma,max_ma),NA))
    gen_off_val_SS <- as.data.frame(gen_off_val %>% group_by(Ma) %>% mutate(SS=cut(Ma, 
                                                                                   breaks = base_SS, 
                                                                                   labels = rev(name_SS),
                                                                                   right=FALSE)))
    gen_off_val_SS$SS <- factor(gen_off_val_SS$SS,
                                levels = c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2"))
    
    melt_off <- melt(gen_off_val_SS, id=c("accepted_name", "SS"), na.rm = TRUE)
    melt_off <- na.omit(melt_off)
    cast_off <- cast(melt_off, accepted_name~SS, length)
    
    #####
    ## all other environments
    gen_else_val <- gen_else_diff%>% group_by(1:n()) %>% mutate(Ma=ifelse(diff<=12,runif(1,min_ma,max_ma),NA))
    gen_else_val_SS <- as.data.frame(gen_else_val %>% group_by(Ma) %>% mutate(SS=cut(Ma,
                                                                                     breaks = base_SS, 
                                                                                     labels = rev(name_SS),
                                                                                     right=FALSE)))
    gen_else_val_SS$SS <- factor(gen_else_val_SS$SS,
                                 levels = c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2"))
    
    melt_else <- melt(gen_else_val_SS, id=c("accepted_name", "SS"), na.rm = TRUE)
    melt_else <- na.omit(melt_else)
    cast_else <- cast(melt_else, accepted_name~SS, length)
    
    
    
    ##########
    ## finding overlapping genera
    overlap <- as.factor(Reduce(intersect, list(cast_on$accepted_name, cast_off$accepted_name)))
    overlap_onelse <- as.factor(Reduce(intersect,list(cast_on$accepted_name, cast_else$accepted_name)))
    overlap_offelse <-   as.factor(Reduce(intersect, list(cast_off$accepted_name, cast_else$accepted_name)))
    
    ## genera that occur in all three groups
    overlap_ALL <- as.factor(Reduce(intersect, list(cast_on$accepted_name, cast_off$accepted_name, cast_else$accepted_name)))
    
    ###########
    ## filtering and creating an input file for all genera that occur in ONSHORE areas ONLY
    inp_on <- cast_on %>% filter(!accepted_name %in% overlap)
    inp_on <- inp_on %>% filter(!accepted_name %in% overlap_onelse)
    inp_on <- as.matrix(inp_on[,2:12])
    inp_on <- ifelse(inp_on >=1, 1, 0)
    data_on1 <- inp_on                                                    ## for later check of distribution of observations
    inp_on <- unite(as.data.frame(inp_on), "ch", c(1:11), sep = "")
    inp_on <- cbind(inp_on, environment = "onshore", ";") ### FOR FINAL INPUT DATASET
    
    #####
    ## filtering and creating an input file for all genera the occur in OFFSHORE areas ONLY
    inp_off <- cast_off %>% filter(!accepted_name %in% overlap)
    inp_off <- cast_off %>% filter(!accepted_name %in% overlap_offelse)
    inp_off <- as.matrix(inp_off[,2:12])
    inp_off <- ifelse(inp_off >=1, 1, 0)
    data_off1 <- inp_off                                                   ## for later check of distribution of observations
    inp_off <- unite(as.data.frame(inp_off), "ch", c(1:11), sep = "")
    inp_off <- cbind(inp_off, environment = "offshore", ";") ### FOR FINAL INPUT DATASET
    
    #####
    ## filtering ALL genera that occur both onshore and offshore, but considering only data from onshore 
    on_OL <- cast_on %>% filter(!accepted_name %in% overlap_onelse)
    on_OL <- on_OL %>% filter(accepted_name %in% overlap)

    ## filtering ALL genera that occur both onshore and offshore, but considering only data from offshore
    off_OL <- cast_off %>% filter(!accepted_name %in% overlap_offelse)
    off_OL <- off_OL %>% filter(accepted_name %in% overlap)
    
    #####
    ## all genera from overlap, which were just extracted and occur only in onshore areas
    ## 1 for presence and 0 for absence
    inp_on_ol <- as.matrix(on_OL[,2:12])
    inp_on_ol <- ifelse(inp_on_ol >=1, 1, 0)
    
    ## all genera from overlap, which only occur in offshore areas
    ## 2 for presence and 0 for absence
    inp_off_ol <- as.matrix(off_OL[,2:12])
    inp_off_ol <- ifelse(inp_off_ol >=1, 2, 0)
    
    ## following matrix adds up the previous two data.frames
    ## 1 = presence onshore
    ## 2 = presence offshore
    ## 3 = presence onshore AND offshore
    inp_onoff <- as.matrix(inp_off_ol)+as.matrix(inp_on_ol) ### FOR MIG-DATA
    
    #####
    ## filtering ALL genera that occur both onshore and anywhere else, but considering only data from onshore
    onelse_OL <- cast_on %>% filter(!accepted_name %in% overlap)
    onelse_OL <- onelse_OL %>% filter(accepted_name %in% overlap_onelse)
    
    ## filtering ALL genera that occur both onshore and anywhere else, but considering only data from anywhere else
    elseon_OL <- cast_else %>% filter(!accepted_name %in% overlap)
    elseon_OL <- elseon_OL %>% filter(accepted_name %in% overlap_onelse)  
    
    #####
    ## data extracted from onshore gets assigned 1 for prensence (as above)
    inp_onelse_ol <- as.matrix(onelse_OL[,2:12])
    inp_onelse_ol <- ifelse(inp_onelse_ol >=1, 1, 0)
    
    ## data extracted from anywhere else gets assigned 10 for presence
    inp_elseon_ol <- as.matrix(elseon_OL[,2:12])
    inp_elseon_ol <- ifelse(inp_elseon_ol >=1, 10, 0)
    
    ## following matrix adds up on and onelse matrices
    ## 1 = presence onshore
    ## 10 = presence anywhere else but onshore and offshore
    ## 11 = presence onshore and anywhere else
    ## 0 = absence
    inp_onelse <- as.matrix(inp_onelse_ol)+as.matrix(inp_elseon_ol) ### FOR MIG-DATA
    
    
    #####
    ## filtering ALL genera that occur both offshore and anywhere else, but considering only data from offshore
    offelse_OL <- cast_off %>% filter(!accepted_name %in% overlap)
    offelse_OL <- offelse_OL %>% filter(accepted_name %in% overlap_offelse)
    
    ## filterin ALL genera that occur both offshore and anywhere else, but considering only data from anywhere else
    elseoff_OL <- cast_else %>% filter(!accepted_name %in% overlap)
    elseoff_OL <- elseoff_OL %>% filter(accepted_name %in% overlap_offelse)  
    
    ## data extracted from offshore gets assigned 2 for presence
    inp_offelse_ol <- as.matrix(offelse_OL[,2:12])
    inp_offelse_ol <- ifelse(inp_offelse_ol >=1, 2, 0)
    
    ## data extracted from anywhere else gets assigned 10 for presence
    inp_elseoff_ol <- as.matrix(elseoff_OL[,2:12])
    inp_elseoff_ol <- ifelse(inp_elseoff_ol >=1, 10, 0)
    
    ## following matrix adds up off and offelse matrices
    ## 2 = presence offshore
    ## 10 = presence anywhere else but onshore and offshore
    ## 12 = presence offshore and anywhere else 
    ## 0 = absence
    inp_offelse <- as.matrix(inp_elseoff_ol)+as.matrix(inp_offelse_ol) ### FOR MIG-DATA
    
    #####
    ## filtering ALL genera that occur both onshore, offshore AND anywhere else
    ## but considering only data from onshore first
    onoffelse_OL <- cast_on %>% filter(accepted_name %in% overlap_ALL)
    inp_onoffelse_ol <- as.matrix(onoffelse_OL[,2:12])
    inp_onoffelse_ol <- ifelse(inp_onoffelse_ol >=1,1,0)
    
    ## considering only data from offshore
    offonelse_OL <- cast_off %>% filter(accepted_name %in% overlap_ALL)
    inp_offonelse_ol <- as.matrix(offonelse_OL[,2:12])
    inp_offonelse_ol <- ifelse(inp_offonelse_ol >=1,2,0)
    
    ## considering only data from environments outside onshore offshore
    elseonoff_OL <- cast_else %>% filter(accepted_name %in% overlap_ALL)
    inp_elseonoff_ol <- as.matrix(elseonoff_OL[,2:12])
    inp_elseonoff_ol <- ifelse(inp_elseonoff_ol >=1,10,0)
    
    inp_onoffelse <- as.matrix(inp_onoffelse_ol)+as.matrix(inp_offonelse_ol)+as.matrix(inp_elseonoff_ol) ### FOR MIG-DATA
    
    
    ## binding all presence/absence matrices w/codes together
    ## substitution of numbers with places of occurrences
    test_mig <- rbind.data.frame(inp_onoff,inp_onelse,inp_offelse,inp_onoffelse) ### EVERYTHING FOR MIG-DATA GOES IN HERE
    test_mig[test_mig == 3] <- "on+off"
    test_mig[test_mig == 2] <- "off"
    test_mig[test_mig == 1] <- "on"
    test_mig[test_mig == 11] <- "on+else"
    test_mig[test_mig == 12] <- "off+else"
    test_mig[test_mig == 13] <- "on+off+else"
    test_mig[test_mig == 10] <- "else"
    test_mig[test_mig == 0] <- NA
    
    
    test <- as.data.frame(test_mig)
    
    ## looking for first occurence (in which time and which place)
    first <-   setNames(data.frame(t(apply(test[1:11], 1, function(x) {
      ind <- which(!is.na((x)))[1]
      c(ind, x[ind])
    }))),
    c("Stage", "first"))
    
    test_first <- cbind(test, "first"=first$first)
    
    ## extracting ONLY genera that have their first occurrence in EITHER onshore OR offshore areas
    extract_inp_first <- filter(test_first, first=="on" | first=="off")
    
    ## if first =="on" --> genus is treated as an onshore genus
    ## if first =="off" -> genus is treated as an offshore genus
    extract_inp_first$environment <- ifelse(extract_inp_first$first=="on", "onshore", "offshore")
    
    ## coding of all presences as 1 and all absences as 0
    inp_mig <- ifelse(!is.na(extract_inp_first[1:11]),1,0)
    data_onoff <- cbind.data.frame(inp_mig, environment = extract_inp_first$environment)
    inp_mig <- unite(as.data.frame(inp_mig),"ch",c(1:11),sep="")
    inp_mig <- cbind(inp_mig, environment = extract_inp_first$environment, ";") ### FOR FINAL INPUT DATASET
    
    
    inp <- rbind.data.frame(inp_on,inp_off,inp_mig)
    
    
    proc.pradsen <- process.data(inp, model= "Pradsen", groups = "environment")
    Time.env <- mark(proc.pradsen, model.parameters = list(Phi=Phi.Time.env, p=p.Time.env, Gamma=Gamma.Time.env))
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  reruns.oo[[i]]=Time.env$results$real
  
}


## Pradlambda parametrization 
## commands to create the input file are as above for the Pradsen parametrization

for(i in 1:100){
  
  tryCatch({

    ## onshore
    gen_on_val <- gen_on_diff%>% group_by(1:n()) %>% mutate(Ma=ifelse(diff<=12,runif(1,min_ma,max_ma),NA))
    gen_on_val_SS <- as.data.frame(gen_on_val %>% group_by(Ma) %>% mutate(SS=cut(Ma, 
                                                                                 breaks = base_SS, 
                                                                                 labels = rev(name_SS),
                                                                                 right=FALSE)))
    gen_on_val_SS$SS <- factor(gen_on_val_SS$SS,
                               levels = c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2"))
    melt_on <- melt(gen_on_val_SS, id=c("accepted_name", "SS"), na.rm = TRUE)
    melt_on <- na.omit(melt_on)
    cast_on <- cast(melt_on, accepted_name~SS, length)
    
    
    #####
    ## offshore
    gen_off_val <- gen_off_diff%>% group_by(1:n()) %>% mutate(Ma=ifelse(diff<=12,runif(1,min_ma,max_ma),NA))
    gen_off_val_SS <- as.data.frame(gen_off_val %>% group_by(Ma) %>% mutate(SS=cut(Ma, 
                                                                                   breaks = base_SS, 
                                                                                   labels = rev(name_SS),
                                                                                   right=FALSE)))
    gen_off_val_SS$SS <- factor(gen_off_val_SS$SS,
                                levels = c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2"))
    
    melt_off <- melt(gen_off_val_SS, id=c("accepted_name", "SS"), na.rm = TRUE)
    melt_off <- na.omit(melt_off)
    cast_off <- cast(melt_off, accepted_name~SS, length)
    
    #####
    ## all other environments
    gen_else_val <- gen_else_diff%>% group_by(1:n()) %>% mutate(Ma=ifelse(diff<=12,runif(1,min_ma,max_ma),NA))
    gen_else_val_SS <- as.data.frame(gen_else_val %>% group_by(Ma) %>% mutate(SS=cut(Ma,
                                                                                     breaks = base_SS, 
                                                                                     labels = rev(name_SS),
                                                                                     right=FALSE)))
    gen_else_val_SS$SS <- factor(gen_else_val_SS$SS,
                                 levels = c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2"))
    
    melt_else <- melt(gen_else_val_SS, id=c("accepted_name", "SS"), na.rm = TRUE)
    melt_else <- na.omit(melt_else)
    cast_else <- cast(melt_else, accepted_name~SS, length)
    
    
    
    ##########
    ## finding overlapping genera
    overlap <- as.factor(Reduce(intersect, list(cast_on$accepted_name, cast_off$accepted_name)))
    overlap_onelse <- as.factor(Reduce(intersect,list(cast_on$accepted_name, cast_else$accepted_name)))
    overlap_offelse <-   as.factor(Reduce(intersect, list(cast_off$accepted_name, cast_else$accepted_name)))
    
    ## genera that occur in all three groups
    overlap_ALL <- as.factor(Reduce(intersect, list(cast_on$accepted_name, cast_off$accepted_name, cast_else$accepted_name)))
    
    ###########
    ## filtering and creating an input file for all genera that occur in ONSHORE areas ONLY
    inp_on <- cast_on %>% filter(!accepted_name %in% overlap)
    inp_on <- inp_on %>% filter(!accepted_name %in% overlap_onelse)
    inp_on <- as.matrix(inp_on[,2:12])
    inp_on <- ifelse(inp_on >=1, 1, 0)
    data_on1 <- inp_on                                                    ## for later check of distribution of observations
    inp_on <- unite(as.data.frame(inp_on), "ch", c(1:11), sep = "")
    inp_on <- cbind(inp_on, environment = "onshore", ";") ### FOR FINAL INPUT DATASET
    
    #####
    ## filtering and creating an input file for all genera the occur in OFFSHORE areas ONLY
    inp_off <- cast_off %>% filter(!accepted_name %in% overlap)
    inp_off <- cast_off %>% filter(!accepted_name %in% overlap_offelse)
    inp_off <- as.matrix(inp_off[,2:12])
    inp_off <- ifelse(inp_off >=1, 1, 0)
    data_off1 <- inp_off                                                   ## for later check of distribution of observations
    inp_off <- unite(as.data.frame(inp_off), "ch", c(1:11), sep = "")
    inp_off <- cbind(inp_off, environment = "offshore", ";") ### FOR FINAL INPUT DATASET
    
    #####
    ## filtering ALL genera that occur both onshore and offshore, but considering only data from onshore 
    on_OL <- cast_on %>% filter(!accepted_name %in% overlap_onelse)
    on_OL <- on_OL %>% filter(accepted_name %in% overlap)
    
    ## filtering ALL genera that occur both onshore and offshore, but considering only data from offshore
    off_OL <- cast_off %>% filter(!accepted_name %in% overlap_offelse)
    off_OL <- off_OL %>% filter(accepted_name %in% overlap)
    
    #####
    ## all genera from overlap, which were just extracted and occur only in onshore areas
    ## 1 for presence and 0 for absence
    inp_on_ol <- as.matrix(on_OL[,2:12])
    inp_on_ol <- ifelse(inp_on_ol >=1, 1, 0)
    
    ## all genera from overlap, which only occur in offshore areas
    ## 2 for presence and 0 for absence
    inp_off_ol <- as.matrix(off_OL[,2:12])
    inp_off_ol <- ifelse(inp_off_ol >=1, 2, 0)
    
    ## following matrix adds up the previous two data.frames
    ## 1 = presence onshore
    ## 2 = presence offshore
    ## 3 = presence onshore AND offshore
    inp_onoff <- as.matrix(inp_off_ol)+as.matrix(inp_on_ol) ### FOR MIG-DATA
    
    #####
    ## filtering ALL genera that occur both onshore and anywhere else, but considering only data from onshore
    onelse_OL <- cast_on %>% filter(!accepted_name %in% overlap)
    onelse_OL <- onelse_OL %>% filter(accepted_name %in% overlap_onelse)
    
    ## filtering ALL genera that occur both onshore and anywhere else, but considering only data from anywhere else
    elseon_OL <- cast_else %>% filter(!accepted_name %in% overlap)
    elseon_OL <- elseon_OL %>% filter(accepted_name %in% overlap_onelse)  
    
    #####
    ## data extracted from onshore gets assigned 1 for prensence (as above)
    inp_onelse_ol <- as.matrix(onelse_OL[,2:12])
    inp_onelse_ol <- ifelse(inp_onelse_ol >=1, 1, 0)
    
    ## data extracted from anywhere else gets assigned 10 for presence
    inp_elseon_ol <- as.matrix(elseon_OL[,2:12])
    inp_elseon_ol <- ifelse(inp_elseon_ol >=1, 10, 0)
    
    ## following matrix adds up on and onelse matrices
    ## 1 = presence onshore
    ## 10 = presence anywhere else but onshore and offshore
    ## 11 = presence onshore and anywhere else
    ## 0 = absence
    inp_onelse <- as.matrix(inp_onelse_ol)+as.matrix(inp_elseon_ol) ### FOR MIG-DATA
    
    
    #####
    ## filtering ALL genera that occur both offshore and anywhere else, but considering only data from offshore
    offelse_OL <- cast_off %>% filter(!accepted_name %in% overlap)
    offelse_OL <- offelse_OL %>% filter(accepted_name %in% overlap_offelse)
    
    ## filterin ALL genera that occur both offshore and anywhere else, but considering only data from anywhere else
    elseoff_OL <- cast_else %>% filter(!accepted_name %in% overlap)
    elseoff_OL <- elseoff_OL %>% filter(accepted_name %in% overlap_offelse)  
    
    ## data extracted from offshore gets assigned 2 for presence
    inp_offelse_ol <- as.matrix(offelse_OL[,2:12])
    inp_offelse_ol <- ifelse(inp_offelse_ol >=1, 2, 0)
    
    ## data extracted from anywhere else gets assigned 10 for presence
    inp_elseoff_ol <- as.matrix(elseoff_OL[,2:12])
    inp_elseoff_ol <- ifelse(inp_elseoff_ol >=1, 10, 0)
    
    ## following matrix adds up off and offelse matrices
    ## 2 = presence offshore
    ## 10 = presence anywhere else but onshore and offshore
    ## 12 = presence offshore and anywhere else 
    ## 0 = absence
    inp_offelse <- as.matrix(inp_elseoff_ol)+as.matrix(inp_offelse_ol) ### FOR MIG-DATA
    
    #####
    ## filtering ALL genera that occur both onshore, offshore AND anywhere else
    ## but considering only data from onshore first
    onoffelse_OL <- cast_on %>% filter(accepted_name %in% overlap_ALL)
    inp_onoffelse_ol <- as.matrix(onoffelse_OL[,2:12])
    inp_onoffelse_ol <- ifelse(inp_onoffelse_ol >=1,1,0)
    
    ## considering only data from offshore
    offonelse_OL <- cast_off %>% filter(accepted_name %in% overlap_ALL)
    inp_offonelse_ol <- as.matrix(offonelse_OL[,2:12])
    inp_offonelse_ol <- ifelse(inp_offonelse_ol >=1,2,0)
    
    ## considering only data from environments outside onshore offshore
    elseonoff_OL <- cast_else %>% filter(accepted_name %in% overlap_ALL)
    inp_elseonoff_ol <- as.matrix(elseonoff_OL[,2:12])
    inp_elseonoff_ol <- ifelse(inp_elseonoff_ol >=1,10,0)
    
    inp_onoffelse <- as.matrix(inp_onoffelse_ol)+as.matrix(inp_offonelse_ol)+as.matrix(inp_elseonoff_ol) ### FOR MIG-DATA
    
    
    ## binding all presence/absence matrices w/codes together
    ## substitution of numbers with places of occurrences
    test_mig <- rbind.data.frame(inp_onoff,inp_onelse,inp_offelse,inp_onoffelse) ### EVERYTHING FOR MIG-DATA GOES IN HERE
    test_mig[test_mig == 3] <- "on+off"
    test_mig[test_mig == 2] <- "off"
    test_mig[test_mig == 1] <- "on"
    test_mig[test_mig == 11] <- "on+else"
    test_mig[test_mig == 12] <- "off+else"
    test_mig[test_mig == 13] <- "on+off+else"
    test_mig[test_mig == 10] <- "else"
    test_mig[test_mig == 0] <- NA
    
    
    test <- as.data.frame(test_mig)
    
    first <-   setNames(data.frame(t(apply(test[1:11], 1, function(x) {
      ind <- which(!is.na((x)))[1]
      c(ind, x[ind])
    }))),
    c("Stage", "first"))
    
    test_first <- cbind(test, "first"=first$first)

    
    extract_inp_first <- filter(test_first, first=="on" | first=="off")
    
    ## define here if I want to have separate groups for onshore migration, offshore migration or include them in on/offshore
    extract_inp_first$environment <- ifelse(extract_inp_first$first=="on", "onshore", "offshore")
    inp_mig <- ifelse(!is.na(extract_inp_first[1:11]),1,0)
    inp_mig <- unite(as.data.frame(inp_mig),"ch",c(1:11),sep="")
    inp_mig <- cbind(inp_mig, environment = extract_inp_first$environment, ";") ### FOR FINAL INPUT DATASET
    
    inp <- rbind.data.frame(inp_on,inp_off,inp_mig)
    
    proc.pradlambda <- process.data(inp, model= "Pradlambda", groups = "environment")
    TimeD.env <- mark(proc.pradlambda, model.parameters = list(Phi=Phi.Time.env, p=p.Time.env, Lambda=L.Time.env))
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  reruns.lambda.oo[[i]]=TimeD.env$results$real
  
}

###########################################################################################
##############################################################################
###averaging from results and results.lambda##################################

mean_results <- rowMeans(do.call(cbind, lapply(reruns.oo, function(x) x$estimate)))

estimates.oo <- do.call(cbind, lapply(reruns.oo, function(x) x$estimate))
lcl.oo <- do.call(cbind, lapply(reruns.oo, function(x) x$lcl))
ucl.oo <- do.call(cbind, lapply(reruns.oo, function(x) x$ucl))

mean_results.lambda <- rowMeans(do.call(cbind, lapply(reruns.lambda.oo, function(x) x$estimate)))

estimates.oo.lambda <- do.call(cbind, lapply(reruns.lambda.oo, function(x) x$estimate))
lcl.oo.lambda <- do.call(cbind, lapply(reruns.lambda.oo, function(x) x$lcl))
ucl.oo.lambda <- do.call(cbind, lapply(reruns.lambda.oo, function(x) x$ucl))

estimates.oo[estimates.oo==0|estimates.oo==1] <- NA
estimates.oo[ucl.oo-lcl.oo==0|ucl.oo-lcl.oo==1]<-NA
lcl.oo[is.na(estimates.oo)] <- NA
ucl.oo[is.na(estimates.oo)] <- NA

estimates.oo.lambda[ucl.oo.lambda-lcl.oo.lambda>10]<-NA
estimates.oo.lambda[ucl.oo.lambda-lcl.oo.lambda==0|ucl.oo.lambda-lcl.oo.lambda==1]<-NA
lcl.oo.lambda[is.na(estimates.oo.lambda)] <- NA
ucl.oo.lambda[is.na(estimates.oo.lambda)] <- NA


#################################################################################
##############################################################################
##############################################################################
## rate estimation from results and results.lambda.
phi_on <- estimates.oo[2:8,]
phiu_on <- ucl.oo[2:8,]
phil_on <- lcl.oo[2:8,]

phi_off <- estimates.oo[12:18,]
phiu_off <- ucl.oo[12:18,]
phil_off <- lcl.oo[12:18,]

p_on <- estimates.oo[23:29,]
pu_on <- ucl.oo[23:29,]
pl_on <- lcl.oo[23:29,]

p_off <- estimates.oo[34:40,]
pu_off <- ucl.oo[34:40,]
pl_off <- lcl.oo[34:40,]

gam_on <- estimates.oo[44:50,]
gamu_on <- ucl.oo[44:50,]
gaml_on <- lcl.oo[44:50,]

gam_off <- estimates.oo[54:60,]
gamu_off <- ucl.oo[54:60,]
gaml_off <- lcl.oo[54:60,]

div_on <- estimates.oo.lambda[44:50,]
divu_on <- ucl.oo.lambda[44:50,]
divl_on <- lcl.oo.lambda[44:50,]

div_off <- estimates.oo.lambda[54:60,]
divu_off <- ucl.oo.lambda[54:60,]
divl_off <- lcl.oo.lambda[54:60,]


## mean results
mean.phi.on <- mean_results[2:8]
mean.p.on <- mean_results[23:29]
mean.gam.on <- mean_results[44:50]
mean.div.on <- mean_results.lambda[44:50]

mean.phi.off <- mean_results[12:18]
mean.p.off <- mean_results[34:40]
mean.gam.off <- mean_results[54:60]
mean.div.off <- mean_results.lambda[54:60]


# following vector is time between midpoints of two neighbouring stageslices
t <- c(9.65,7.7,5.2,5.8,7.15,6.6,4.6)
# following vector is time per interval/stageslice ## used to calculate sampling rate (p)
tp <- c(7.7,7.7,2.7,8.9,5.4,7.8,1.4)

on_origprob <- 1-gam_on
on_origCIl <- 1-gaml_on
on_origCIu <- 1-gamu_on
on_Orig_rate <- -log(1-on_origprob)/t
on_Orig_rate_CIl <- -log(1-on_origCIl)/t
on_Orig_rate_CIu <- -log(1-on_origCIu)/t

off_origprob <- 1-gam_off
off_origCIl <- 1-gaml_off
off_origCIu <- 1-gamu_off
off_Orig_rate <- -log(1-off_origprob)/t
off_Orig_rate_CIl <- -log(1-off_origCIl)/t
off_Orig_rate_CIu <- -log(1-off_origCIu)/t


on_extinctprob <- 1-phi_on
on_extinctCIl <- 1-phil_on
on_extinctCIu <- 1-phiu_on
on_Ext_rate <- -log(1-on_extinctprob)/t
on_Ext_rate_CIl <- -log(1-on_extinctCIl)/t
on_Ext_rate_CIu <- -log(1-on_extinctCIu)/t

off_extinctprob <- 1-phi_off
off_extinctCIl <- 1-phil_off
off_extinctCIu <- 1-phiu_off
off_Ext_rate <- -log(1-off_extinctprob)/t
off_Ext_rate_CIl <- -log(1-off_extinctCIl)/t
off_Ext_rate_CIu <- -log(1-off_extinctCIu)/t


rate_p_on <- -log(1-p_on)/tp
ratel_p_on <- -log(1-pl_on)/tp
rateu_p_on <- -log(1-pu_on)/tp

rate_p_off <- -log(1-p_off)/tp
ratel_p_off <- -log(1-pl_off)/tp
rateu_p_off <- -log(1-pu_off)/tp


## mean rates
mean.origprob.on <- 1-mean.gam.on
mean.origrate.on <- -log(1-mean.origprob.on)/t
mean.extprob.on <- 1-mean.phi.on
mean.extrate.on <- -log(1-mean.extprob.on)/t
mean.prate.on <- -log(1-mean.p.on)/tp

mean.origprob.off <- 1-mean.gam.off
mean.origrate.off <- -log(1-mean.origprob.off)/t
mean.extprob.off <- 1-mean.phi.off
mean.extrate.off <- -log(1-mean.extprob.off)/t
mean.prate.off <- -log(1-mean.p.off)/tp



##################################################################################
## run env.model #################################################################
# 
# 
# genus_Ord <- filter(genus, max_ma < 481.55 & max_ma > 443.8 | min_ma < 481.55 & min_ma >= 443.8)
# 
# gen_on <- genus_Ord %>% filter(environment %in% c("foreshore","coastal indet.",
#                                                   "marginal marine indet.",
#                                                   "carbonate indet.","peritidal",
#                                                   "sand shoal", "shoreface"))
# gen_on_diff <- cbind(gen_on[,c("occurrence_no","accepted_name","max_ma","min_ma","geoplate","phylum")],
#                      "diff"=gen_on$max_ma-gen_on$min_ma)
# 
# 
# gen_off <- genus_Ord %>% filter(environment %in% c("offshore indet.","offshore","offshore ramp","offshore shelf"))
# gen_off_diff <- cbind(gen_off[,c("occurrence_no","accepted_name","max_ma","min_ma","geoplate","phylum")],
#                       "diff"=gen_off$max_ma-gen_off$min_ma)
# 
# 
# 
# 
# gen_else <- genus_Ord %>% filter(!environment %in% c("offshore indet.","offshore","offshore ramp","offshore shelf",
#                                                      "foreshore","coastal indet.",
#                                                      "marginal marine indet.",
#                                                      "carbonate indet.","peritidal",
#                                                      "sand shoal", "shoreface"))
# gen_else_diff <- cbind(gen_else[,c("occurrence_no","accepted_name","phylum","max_ma","min_ma","geoplate","phylum")],
#                        "diff"=gen_else$max_ma-gen_else$min_ma)
# 
# 
# 
# gen_on_val <- gen_on_diff%>% group_by(1:n()) %>% mutate(Ma=ifelse(diff<=12,runif(1,min_ma,max_ma),NA))
# base_SS <- c(509,497,485.4,477.7,470,467.3,458.4,453,445.2,443.8,440.8,438.5)
# name_SS<-c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2")
# gen_on_val_SS <- as.data.frame(gen_on_val %>% group_by(Ma) %>% mutate(SS=cut(Ma, 
#                                                                              breaks = base_SS, 
#                                                                              labels = rev(name_SS),
#                                                                              right=FALSE)))
# gen_on_val_SS$SS <- factor(gen_on_val_SS$SS,
#                            levels = c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2"))
# melt_on <- melt(gen_on_val_SS, id=c("accepted_name", "SS"), na.rm = TRUE)
# melt_on <- na.omit(melt_on)
# cast_on <- cast(melt_on, accepted_name~SS, length)
# 
# 
# #####
# gen_off_val <- gen_off_diff%>% group_by(1:n()) %>% mutate(Ma=ifelse(diff<=12,runif(1,min_ma,max_ma),NA))
# gen_off_val_SS <- as.data.frame(gen_off_val %>% group_by(Ma) %>% mutate(SS=cut(Ma, 
#                                                                                breaks = base_SS, 
#                                                                                labels = rev(name_SS),
#                                                                                right=FALSE)))
# gen_off_val_SS$SS <- factor(gen_off_val_SS$SS,
#                             levels = c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2"))
# 
# melt_off <- melt(gen_off_val_SS, id=c("accepted_name", "SS"), na.rm = TRUE)
# melt_off <- na.omit(melt_off)
# cast_off <- cast(melt_off, accepted_name~SS, length)
# 
# #####
# gen_else_val <- gen_else_diff%>% group_by(1:n()) %>% mutate(Ma=ifelse(diff<=12,runif(1,min_ma,max_ma),NA))
# gen_else_val_SS <- as.data.frame(gen_else_val %>% group_by(Ma) %>% mutate(SS=cut(Ma,
#                                                                                  breaks = base_SS, 
#                                                                                  labels = rev(name_SS),
#                                                                                  right=FALSE)))
# gen_else_val_SS$SS <- factor(gen_else_val_SS$SS,
#                              levels = c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2"))
# 
# melt_else <- melt(gen_else_val_SS, id=c("accepted_name", "SS"), na.rm = TRUE)
# melt_else <- na.omit(melt_else)
# cast_else <- cast(melt_else, accepted_name~SS, length)
# 
# # #
# ##########
# ## finding overlapping genera
# overlap <- as.factor(Reduce(intersect, list(cast_on$accepted_name, cast_off$accepted_name)))
# overlap_onelse <- as.factor(Reduce(intersect,list(cast_on$accepted_name, cast_else$accepted_name)))
# overlap_offelse <-   as.factor(Reduce(intersect, list(cast_off$accepted_name, cast_else$accepted_name)))
# 
# ## genera that occur in all three groups
# overlap_ALL <- as.factor(Reduce(intersect, list(cast_on$accepted_name, cast_off$accepted_name, cast_else$accepted_name)))
# 
# ###########
# ## filtering and creating an input file for all genera that occur in ONSHORE areas ONLY
# inp_on <- cast_on %>% filter(!accepted_name %in% overlap)
# inp_on <- inp_on %>% filter(!accepted_name %in% overlap_onelse)
# genus_on <- inp_on[1]           ## for phylum distribution check
# inp_on <- as.matrix(inp_on[,2:8])
# inp_on <- ifelse(inp_on >=1, 1, 0)
# data_on1 <- inp_on                                                    ## for later check of distribution of observations
# inp_on <- unite(as.data.frame(inp_on), "ch", c(1:7), sep = "")
# inp_on <- cbind(inp_on, environment = "onshore", ";") ### FOR FINAL INPUT DATASET
# 
# #####
# ## filtering and creating an input file for all genera the occur in OFFSHORE areas ONLY
# inp_off <- cast_off %>% filter(!accepted_name %in% overlap)
# inp_off <- cast_off %>% filter(!accepted_name %in% overlap_offelse)
# genus_off <- inp_off[1]
# inp_off <- as.matrix(inp_off[,2:8])
# inp_off <- ifelse(inp_off >=1, 1, 0)
# data_off1 <- inp_off                                                   ## for later check of distribution of observations
# inp_off <- unite(as.data.frame(inp_off), "ch", c(1:7), sep = "")
# inp_off <- cbind(inp_off, environment = "offshore", ";") ### FOR FINAL INPUT DATASET
# 
# #####
# ## filtering ALL genera that occur both onshore and offshore, but considering only data from onshore 
# on_OL <- cast_on %>% filter(!accepted_name %in% overlap_onelse)
# on_OL <- on_OL %>% filter(accepted_name %in% overlap)
# 
# ## filtering ALL genera that occur both onshore and offshore, but considering only data from offshore
# off_OL <- cast_off %>% filter(!accepted_name %in% overlap_offelse)
# off_OL <- off_OL %>% filter(accepted_name %in% overlap)
# 
# #####
# ## all genera from overlap, which were just extracted and occur only in onshore areas
# ## 1 for presence and 0 for absence
# inp_on_ol <- as.matrix(on_OL[,2:8])
# inp_on_ol <- ifelse(inp_on_ol >=1, 1, 0)
# 
# ## all genera from overlap, which only occur in offshore areas
# ## 2 for presence and 0 for absence
# inp_off_ol <- as.matrix(off_OL[,2:8])
# inp_off_ol <- ifelse(inp_off_ol >=1, 2, 0)
# 
# ## following matrix adds up the previous two data.frames
# ## 1 = presence onshore
# ## 2 = presence offshore
# ## 3 = presence onshore AND offshore
# inp_onoff <- as.matrix(inp_off_ol)+as.matrix(inp_on_ol) ### FOR MIG-DATA
# phyla_onoff <- cbind.data.frame(on_OL[1],inp_onoff) ### FOR LATER PHYLA IN ORD/ON AND OFF
# 
# #####
# ## filtering ALL genera that occur both onshore and anywhere else, but considering only data from onshore
# onelse_OL <- cast_on %>% filter(!accepted_name %in% overlap)
# onelse_OL <- onelse_OL %>% filter(accepted_name %in% overlap_onelse)
# 
# ## filtering ALL genera that occur both onshore and anywhere else, but considering only data from anywhere else
# elseon_OL <- cast_else %>% filter(!accepted_name %in% overlap)
# elseon_OL <- elseon_OL %>% filter(accepted_name %in% overlap_onelse)  
# 
# #####
# ## data extracted from onshore gets assigned 1 for prensence (as above)
# inp_onelse_ol <- as.matrix(onelse_OL[,2:8])
# inp_onelse_ol <- ifelse(inp_onelse_ol >=1, 1, 0)
# 
# ## data extracted from anywhere else gets assigned 10 for presence
# inp_elseon_ol <- as.matrix(elseon_OL[,2:8])
# inp_elseon_ol <- ifelse(inp_elseon_ol >=1, 10, 0)
# 
# ## following matrix adds up on and onelse matrices
# ## 1 = presence onshore
# ## 10 = presence anywhere else but onshore and offshore
# ## 11 = presence onshore and anywhere else
# ## 0 = absence
# inp_onelse <- as.matrix(inp_onelse_ol)+as.matrix(inp_elseon_ol) ### FOR MIG-DATA
# phyla_onelse <- cbind.data.frame(elseon_OL[1],inp_onelse) ### FOR LATER PHYLA IN ORD/ON AND OFF
# 
# 
# 
# #####
# ## filtering ALL genera that occur both offshore and anywhere else, but considering only data from offshore
# offelse_OL <- cast_off %>% filter(!accepted_name %in% overlap)
# offelse_OL <- offelse_OL %>% filter(accepted_name %in% overlap_offelse)
# 
# ## filterin ALL genera that occur both offshore and anywhere else, but considering only data from anywhere else
# elseoff_OL <- cast_else %>% filter(!accepted_name %in% overlap)
# elseoff_OL <- elseoff_OL %>% filter(accepted_name %in% overlap_offelse)  
# 
# ## data extracted from offshore gets assigned 2 for presence
# inp_offelse_ol <- as.matrix(offelse_OL[,2:8])
# inp_offelse_ol <- ifelse(inp_offelse_ol >=1, 2, 0)
# 
# ## data extracted from anywhere else gets assigned 10 for presence
# inp_elseoff_ol <- as.matrix(elseoff_OL[,2:8])
# inp_elseoff_ol <- ifelse(inp_elseoff_ol >=1, 10, 0)
# 
# ## following matrix adds up off and offelse matrices
# ## 2 = presence offshore
# ## 10 = presence anywhere else but onshore and offshore
# ## 12 = presence offshore and anywhere else 
# ## 0 = absence
# inp_offelse <- as.matrix(inp_elseoff_ol)+as.matrix(inp_offelse_ol) ### FOR MIG-DATA
# phyla_offelse <- cbind.data.frame(offelse_OL[1],inp_offelse) ### FOR LATER PHYLA IN ORD/ON AND OFF
# 
# 
# #####
# ## filtering ALL genera that occur both onshore, offshore AND anywhere else
# ## but considering only data from onshore first
# onoffelse_OL <- cast_on %>% filter(accepted_name %in% overlap_ALL)
# inp_onoffelse_ol <- as.matrix(onoffelse_OL[,2:8])
# inp_onoffelse_ol <- ifelse(inp_onoffelse_ol >=1,1,0)
# 
# ## considering only data from offshore
# offonelse_OL <- cast_off %>% filter(accepted_name %in% overlap_ALL)
# inp_offonelse_ol <- as.matrix(offonelse_OL[,2:8])
# inp_offonelse_ol <- ifelse(inp_offonelse_ol >=1,2,0)
# 
# ## considering only data from environments outside onshore offshore
# elseonoff_OL <- cast_else %>% filter(accepted_name %in% overlap_ALL)
# inp_elseonoff_ol <- as.matrix(elseonoff_OL[,2:8])
# inp_elseonoff_ol <- ifelse(inp_elseonoff_ol >=1,10,0)
# 
# inp_onoffelse <- as.matrix(inp_onoffelse_ol)+as.matrix(inp_offonelse_ol)+as.matrix(inp_elseonoff_ol) ### FOR MIG-DATA
# phyla_onoffelse <- cbind.data.frame(onoffelse_OL[1],inp_onoffelse) ### FOR LATER PHYLA IN ORD/ON AND OFF
# 
# 
# ## binding all presence/absence matrices w/codes together
# ## substitution of numbers with places of occurrences
# test_mig <- rbind.data.frame(inp_onoff,inp_onelse,inp_offelse,inp_onoffelse) ### EVERYTHING FOR MIG-DATA GOES IN HERE
# test_mig[test_mig == 3] <- "on+off"
# test_mig[test_mig == 2] <- "off"
# test_mig[test_mig == 1] <- "on"
# test_mig[test_mig == 11] <- "on+else"
# test_mig[test_mig == 12] <- "off+else"
# test_mig[test_mig == 13] <- "on+off+else"
# test_mig[test_mig == 10] <- "else"
# test_mig[test_mig == 0] <- NA
# 
# 
# test <- as.data.frame(test_mig)
# 
# first <-   setNames(data.frame(t(apply(test[1:7], 1, function(x) {
#   ind <- which(!is.na((x)))[1]
#   c(ind, x[ind])
# }))),
# c("Stage", "first"))
# 
# test_first <- cbind(test, "first"=first$first)
# 
# 
# extract_inp_first <- filter(test_first, first=="on" | first=="off")
# 
# ## depending on how I assign groups I could define onshore-migration and offshore-migration here! first=="on", "onmig", "offmig") 
# extract_inp_first$environment <- ifelse(extract_inp_first$first=="on", "onshore", "offshore")
# inp_mig <- ifelse(!is.na(extract_inp_first[1:7]),1,0)
# inp_mig <- unite(as.data.frame(inp_mig),"ch",c(1:7),sep="")
# inp_mig <- cbind(inp_mig, environment = extract_inp_first$environment, ";") ### FOR FINAL INPUT DATASET
# 
# 
# ## for all first and seconds assigned to on/off
# inp <- rbind.data.frame(inp_on,inp_off,inp_mig)
# 
# 
# inp <- inp[inp$ch != "0000000",]
# 
# ###############
# 
# 
# proc.pradsen <- process.data(inp, model= "Pradsen", groups = "environment",
#                              time.intervals = c(7.7,5.2,5.8,7.15,6.6,4.6))#Pradel's survival and "growth"
# 
# ## time.intervals = c(12,12,11.6,7.7,7.7,2.7,8.9,5.4,7.8,1.4,3,2.3))#Pradel's survival and "growth"
# 
# envirn <- mark(proc.pradsen, model.parameters = list(Phi=Phi.env, p=p.env, Gamma= Gamma.env))
# 
# 
# proc.pradlambda <- process.data(inp, model= "Pradlambda", groups = "environment",
#                                 time.intervals = c(7.7,5.2,5.8,7.15,6.6,4.6))#Pradel's survival and "growth"
# 
# envirnD <-mark(proc.pradlambda, model.parameters = list(Phi=Phi.env, p=p.env, Lambda= L.env))
# 
# 
# 
# ##############################################################################
# ### results from environment model ###########################################
# phi_onA <- envirn$results$real$estimate[1]
# phiu_onA <- envirn$results$real$ucl[1]
# phil_onA <- envirn$results$real$lcl[1]
# 
# phi_offA <- envirn$results$real$estimate[2]
# phiu_offA <- envirn$results$real$ucl[2]
# phil_offA <- envirn$results$real$lcl[2]
# 
# p_onA <- envirn$results$real$estimate[3]
# pu_onA <- envirn$results$real$ucl[3]
# pl_onA <- envirn$results$real$lcl[3]
# 
# p_offA <- envirn$results$real$estimate[4]
# pu_offA <- envirn$results$real$ucl[4]
# pl_offA <- envirn$results$real$lcl[4]
# 
# gam_onA <- envirn$results$real$estimate[5]
# gamu_onA <- envirn$results$real$ucl[5]
# gaml_onA <- envirn$results$real$lcl[5]
# 
# gam_offA <- envirn$results$real$estimate[6]
# gamu_offA <- envirn$results$real$ucl[6]
# gaml_offA <- envirn$results$real$lcl[6]
# 
# div_onA <- envirnD$results$real$estimate[5]
# divu_onA <- envirnD$results$real$ucl[5]
# divl_onA <- envirnD$results$real$lcl[5]
# 
# div_offA <- envirnD$results$real$estimate[6]
# divu_offA <- envirnD$results$real$ucl[6]
# divl_offA <- envirnD$results$real$lcl[6]
# 
# 
# on_origprobA <- 1-gam_onA
# on_origCIlA <- 1-gaml_onA
# on_origCIuA <- 1-gamu_onA
# 
# 
# off_origprobA <- 1-gam_offA
# off_origCIlA <- 1-gaml_offA
# off_origCIuA <- 1-gamu_offA
# 
# 
# on_extinctprobA <- 1-phi_onA
# on_extinctCIlA <- 1-phil_onA
# on_extinctCIuA <- 1-phiu_onA
# 
# off_extinctprobA <- 1-phi_offA
# off_extinctCIlA <- 1-phil_offA
# off_extinctCIuA <- 1-phiu_offA
# # 
##############################################################################
####### plot # FIGURE 4 ## (origination rates) ##############################
Stagebase <-c(485.4,477.7,470,467.3,458.4,453,445.2)
Stagemidpoints <- c(481.55,473.85,468.65,462.85,455.7,449.1,444.5)


tpx <- c("Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Sil")
per <- data.frame(c(485.4,477.7,470,467.3,458.4,453,445.2,443.8))
Ord.Stage <- cbind.data.frame(tpx,"Stage"=per)

##
## timescale function from
## http://simpson-carl.github.io/articles/15/timescales.to.base
## continue here on monday <- data can be downloaded.. :)

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
####### plot # FIGURE 4 ## (origination rates) ##############################

par(mfrow=c(2,2), mar = c(3,1,0.1,1), oma = c(0,2,0,0))
plot(Stagebase-0.5, on_Orig_rate[,27], type = "b",
     pch=19,
     ylim = c(-0.05,0.35),
     xlim = rev(c(444.18,485.4)),
     axes = F,
     xlab = "",
     ylab = "")

tscales.Ord(0.35, 0, -0.05)

lines(Stagebase-0.5, on_Orig_rate[,27], type = "b", pch=19, lwd=0.8)
text(Stagebase[1]-2, on_Orig_rate[1,27]+0.02, labels="on", lwd=0.8)
arrows(x0=Stagebase-0.5, y0=on_Orig_rate_CIl[,27], x1=Stagebase-0.5, y1=on_Orig_rate_CIu[,27], length=0.02, lwd = 0.8, angle = 90, code = 3)
lines(Stagebase, off_Orig_rate[,27], type = "b",  pch=17, lwd=1.2, lty=2)
text(Stagebase[1]-2, off_Orig_rate[1,27]-0.02,labels="off",lwd=1.2)
arrows(x0=Stagebase, y0=off_Orig_rate_CIl[,27], x1=Stagebase, y1=off_Orig_rate_CIu[,27], length=0.02, lwd = 1.2, angle = 90, code = 3, lty=2)

legend("topleft", legend="A", bty="n", cex = 1.25)


axis(2, col = 'grey75', line = -0.2, at = seq(0, 0.35, 0.05))

mtext("Origination events per myr", side = 2, line = 2, cex = 0.75)

############################################################################
plot(Stagebase-0.5, on_Ext_rate[,27], type = "b",
     pch=17,
     ylim = c(-0.02,0.25),
     xlim = rev(c(444.18,485.4)),
     axes = F,
     xlab = "",
     ylab = "")

tscales.Ord(0.25, 0, -0.02)

lines(Stagebase-0.5, on_Ext_rate[,27], type = "b",lwd=0.8, pch=19)
arrows(x0=Stagebase-0.5, y0=on_Ext_rate_CIl[,27], x1=Stagebase-0.5, y1=on_Ext_rate_CIu[,27], length=0.02, lwd = 0.8, angle = 90, code = 3)
lines(Stagebase, off_Ext_rate[,27], type = "b",  pch=17, lwd=1.2, lty=2)
arrows(x0=Stagebase, y0=off_Ext_rate_CIl[,27], x1=Stagebase, y1=off_Ext_rate_CIu[,27], length=0.02, lwd = 1.2, angle = 90, code = 3, lty=2)

legend("topleft", legend="B", bty="n", cex = 1.25)



axis(2, col = 'grey75', line = -0.2, at = seq(0, 0.25, 0.05))

mtext("Extinction events per myr", side = 2, line = 2, cex = 0.75)

##############################################################################
plot(Stagebase-0.5, div_on[,27]-1, type = "b",
     pch=17,
     ylim = c(-1.5,10),
     xlim = rev(c(444.18,485.4)),
     axes = F,
     xlab = "",
     ylab = "")

tscales.Ord(10, -1, -1.5)
abline(h = 0, col="darkgrey")


lines(Stagebase-0.5, div_on[,27]-1, type = "b", pch=19, lwd=0.8)
arrows(x0=Stagebase-0.5, y0=divl_on[,27]-1, x1=Stagebase-0.5, y1=divu_on[,27]-1, length=0.02, lwd = 0.8, angle = 90, code = 3)
lines(Stagebase,div_off[,27]-1, type = "b",  pch=17, lwd=1.2, lty=2)
arrows(x0=Stagebase, y0=divl_off[,27]-1, x1=Stagebase, y1=divu_off[,27]-1, length=0.02, lwd = 1.2, angle = 90, code = 3, lty=2)

legend("topleft", legend="C", bty="n", cex = 1.25)


axis(1, col = 'grey75', line = -0.1, at = seq(445,485,10))
axis(2, col = 'grey75', line = -0.2, at = seq(-1,10, 2))

mtext("Age (Ma)", side = 1, line = 2, cex = 0.75)
mtext("Net diversification rate per myr", side = 2, line = 2, cex = 0.75)

##############################################################################
plot(Stagemidpoints-0.2, rate_p_on[,27], type = 'b',
     pch=17,
     xlim = rev(c(444.18,485.4)),
     ylim = c(-0.1,1),
     axes = F,
     xlab = "", ylab = "")

tscales.Ord(1, 0, -0.1)

lines(Stagemidpoints-0.2, rate_p_on[,27], type = "b",pch=19, lwd=0.8)
arrows(x0=Stagemidpoints-0.2, y0=ratel_p_on[,27], x1=Stagemidpoints-0.2, y1=rateu_p_on[,27], length=0.02, lwd = 0.8, angle = 90, code = 3)
lines(Stagemidpoints+0.2, rate_p_off[,27], type = "b" ,  pch=17, lwd=1.2, lty=2)
arrows(x0=Stagemidpoints+0.2, y0=ratel_p_off[,27], x1=Stagemidpoints+0.2, y1=rateu_p_off[,27], length=0.02, lwd = 1.2, angle = 90, code = 3,lty=2)

legend("topleft", legend="D", bty="n", cex = 1.25)

axis(1, col = 'grey75', line = 0.15, at=seq(445,485,10))
axis(2, col = 'grey75', line = -0.15, at = seq(0, 1, 0.2))

mtext("Age (Ma)", side = 1, line = 2, cex = 0.75)
mtext("Sampling events per myr", side = 2, line = 2, cex = 0.75)


##############################################################################
####### Figure S9 ############################################################
####### replicate plot from 100 runs  ########################################
par(mfrow=c(2,2), mar = c(3,1,0.1,1), oma = c(0,2,0,0))
plot(Stagebase-0.5, on_Orig_rate[,27], type = "b",
     pch=19,
     # main = "Origination Rate",
     ylim = c(-0.05,1),
     xlim = rev(c(444.18,485.4)),
     axes = F,
     xlab = "",
     ylab = "")

tscales.Ord(1, 0, -0.05)


for (i in 1:100){
  lines(Stagebase-0.2, on_Orig_rate_CIl[,i], col = "grey70")}
for (i in 1:100){
  lines(Stagebase-0.2, on_Orig_rate_CIu[,i], col = "grey70")
}  


for (i in 1:100){
  lines(Stagebase, off_Orig_rate_CIl[,i], col = "grey55")}
for (i in 1:100){
  lines(Stagebase, off_Orig_rate_CIu[,i], col = "grey55")}


lines(Stagebase-0.5, mean.origrate.on, type = "b", pch=19, lwd=0.8)
text(Stagebase[1]-2, mean.origrate.on[1]+0.02, labels="on", lwd=0.8)
lines(Stagebase, mean.origrate.off, type = "b",  pch=17, lwd=1.2, lty=2)
text(Stagebase[1]-2, mean.origrate.off[1]-0.02,labels="off",lwd=1.2)

legend("topleft", legend="A", bty="n", cex = 1.25)

axis(2, col = 'grey75', line = -0.2, at = seq(0, 1, 0.2))

mtext("Origination events per myr", side = 2, line = 2, cex = 0.75)

############################################################################
plot(Stagebase-0.5, mean.extrate.on, type = "b",
     pch=17,
     ylim = c(-0.02,0.5),
     xlim = rev(c(444.18,485.4)),
     axes = F,
     xlab = "",
     ylab = "")

tscales.Ord(0.5, 0, -0.02)

for (i in 1:100){
  lines(Stagebase-0.2, on_Ext_rate_CIl[,i], col = "grey70")}
for (i in 1:100){
  lines(Stagebase-0.2, on_Ext_rate_CIu[,i], col = "grey70")
}  


for (i in 1:100){
  lines(Stagebase, off_Ext_rate_CIl[,i], col = "grey55")}
for (i in 1:100){
  lines(Stagebase, off_Ext_rate_CIu[,i], col = "grey55")}


lines(Stagebase-0.5, mean.extrate.on, type = "b",lwd=0.8, pch=19)
lines(Stagebase, mean.extrate.off, type = "b",  pch=17, lwd=1.2, lty=2)

legend("topleft", legend="B", bty="n", cex = 1.25)

axis(2, col = 'grey75', line = -0.2, at = seq(0, 0.5, 0.1))

mtext("Extinction events per myr", side = 2, line = 2, cex = 0.75)

##############################################################################
plot(Stagebase-0.5, mean.div.on-1, type = "b",
     pch=17,
     ylim = c(-1.5,10),
     xlim = rev(c(444.18,485.4)),
     axes = F,
     xlab = "",
     ylab = "")

tscales.Ord(10, -1, -1.5)
abline(h = 0, col="darkgrey")

for (i in 1:100){
  lines(Stagebase-0.2, divl_on[,i]-1, col = "grey70")}
for (i in 1:100){
  lines(Stagebase-0.2, divu_on[,i]-1, col = "grey70")
}  


for (i in 1:100){
  lines(Stagebase, divl_off[,i]-1, col = "grey55")}
for (i in 1:100){
  lines(Stagebase, divu_off[,i]-1, col = "grey55")}




lines(Stagebase-0.5, mean.div.on-1, type = "b", pch=19, lwd=0.8)
lines(Stagebase, mean.div.off-1, type = "b",  pch=17, lwd=1.2, lty=2)

legend("topleft", legend="C", bty="n", cex = 1.25)

axis(1, col = 'grey75', line = -0.1, at = seq(445,485,10))
axis(2, col = 'grey75', line = -0.2, at = seq(-1,10, 2))

mtext("Age (Ma)", side = 1, line = 2, cex = 0.75)
mtext("Net diversification rate per myr", side = 2, line = 2, cex = 0.75)

##############################################################################
plot(Stagemidpoints-0.2, mean.prate.on, type = 'b',
     pch=17,
     xlim = rev(c(444.18,485.4)),
     ylim = c(-0.1,1),
     axes = F,
     xlab = "", ylab = "")

tscales.Ord(1, 0, -0.1)


for (i in 1:100){
  lines(Stagemidpoints-0.2, ratel_p_on[,i], col = "grey70")}
for (i in 1:100){
  lines(Stagemidpoints-0.2, rateu_p_on[,i], col = "grey70")
}  


for (i in 1:100){
  lines(Stagemidpoints, ratel_p_off[,i], col = "grey55")}
for (i in 1:100){
  lines(Stagemidpoints, rateu_p_off[,i], col = "grey55")}



lines(Stagemidpoints-0.2, mean.prate.on, type = "b",pch=19, lwd=0.8)
lines(Stagemidpoints+0.2, mean.prate.off, type = "b" ,  pch=17, lwd=1.2, lty=2)

legend("topleft", legend="D", bty="n", cex = 1.25)

axis(1, col = 'grey75', line = 0.15, at=seq(445,485,10))
axis(2, col = 'grey75', line = -0.15, at = seq(0, 1, 0.2))

mtext("Age (Ma)", side = 1, line = 2, cex = 0.75)
mtext("Sampling events per myr", side = 2, line = 2, cex = 0.75)




##############################################################################
#############################################################################
## phyla represented by genera onshore/offshore
## 
## run data-transformation for env.model before plotting phyla !!!
## otherwise input data for this plot is not generated
##
## finding unique genus-phylum combinations per onshore/offshore

genusg <- filter(genus, grepl("genus", accepted_rank))
genus_Ord <- filter(genusg, max_ma < 481.55 & max_ma > 443.8 | min_ma < 481.55 & min_ma > 443.8)
genus_phyla <- unique(genus_Ord[,c("accepted_name","phylum")])


###########################################################################
# genus_on <- filter(gen_on_diff, grepl("genus", accepted_rank)) # not needed here because we filtered only "genus" rank already before
genus_phylon <- filter(gen_on_diff, max_ma < 481.55 & max_ma > 443.8 | min_ma < 481.55 & min_ma > 443.8)
genus_phylon <- unique(genus_phylon[,c("accepted_name","phylum")])

###########################################################################
genus_phyloff <- filter(gen_off_diff, max_ma < 481.55 & max_ma > 443.8 | min_ma < 481.55 & min_ma > 443.8)
genus_phyloff <- unique(genus_phyloff[,c("accepted_name","phylum")])

genus_phylelse <- filter(gen_else_diff, max_ma < 481.55 & max_ma > 443.8 | min_ma < 481.55 & min_ma > 443.8)
genus_phylelse <- unique(genus_phylelse[,c("accepted_name","phylum")])


##########
overlap_phyla_onoff <- Reduce(intersect, list(genus_phylon$accepted_name,genus_phyloff$accepted_name))
overlap_phyla_onelse <- Reduce(intersect, list(genus_phylon$accepted_name,genus_phylelse$accepted_name))
overlap_phyla_offelse <- Reduce(intersect, list(genus_phyloff$accepted_name,genus_phylelse$accepted_name))

## for getting total amount of overlapping genera onshore and offshore
overlap_onoffall <- c(overlap_phyla_onoff,overlap_phyla_onelse,overlap_phyla_offelse)
overlap_onoffall <- unique(overlap_onoffall)


## pure offs and ons
phyl_off <- genus_phyloff %>% filter(!accepted_name %in% overlap_phyla_onoff)
phyl_off <- phyl_off %>% filter(!accepted_name %in% overlap_phyla_offelse)

phyl_on <- genus_phylon %>% filter(!accepted_name %in% overlap_phyla_onoff)
phyl_on <- phyl_on %>% filter(!accepted_name %in% overlap_phyla_onelse)


## overlapping data
phyla_mig <- rbind.data.frame(phyla_offelse,phyla_onelse,phyla_onoff,phyla_onoffelse) ### phyla-datasets from above (env.model)
phyla_mig[phyla_mig == 3] <- "on+off"
phyla_mig[phyla_mig == 2] <- "off"
phyla_mig[phyla_mig == 1] <- "on"
phyla_mig[phyla_mig == 11] <- "on+else"
phyla_mig[phyla_mig == 12] <- "off+else"
phyla_mig[phyla_mig == 13] <- "on+off+else"
phyla_mig[phyla_mig == 10] <- "else"
phyla_mig[phyla_mig == 0] <- NA


phyla_mig <- as.data.frame(phyla_mig)

first <-   setNames(data.frame(t(apply(phyla_mig[2:8], 1, function(x) {
  ind <- which(!is.na((x)))[1]
  c(ind, x[ind])
}))),
c("Stage", "first"))

phyla_first <- cbind(phyla_mig, "first"=first$first)


phyla_first <- filter(phyla_first, first=="on" | first=="off")
phyla_first$environment <- ifelse(phyla_first$first=="on", "onshore", "offshore")


extract_inp_first <- cbind.data.frame("accepted_name"=phyla_first$accepted_name,"environment"=phyla_first$environment)

phyla_migon <- filter(extract_inp_first, environment=="onshore")
phyla_migon <- genus_phyla %>% filter(accepted_name %in% phyla_migon$accepted_name)

phyla_migoff <- filter(extract_inp_first, environment=="offshore")
phyla_migoff <- genus_phyla %>% filter(accepted_name %in% phyla_migoff$accepted_name)


## combine pure and migrations
phyla_onshore <- rbind.data.frame(phyl_on, phyla_migon)
phyla_offshore <- rbind.data.frame(phyl_off, phyla_migoff)

##########


nend_on <- aggregate(list(phyla_onshore$accepted_name),
                     list(phyla_onshore$phylum),
                     length)


nend_off <- aggregate(list(phyla_offshore$accepted_name),
                      list(phyla_offshore$phylum),
                      length)


on_phy_N <- colSums(nend_on[2], na.rm = T) ## total amount of "endemic" genera during the Ordovician onshore
off_phy_N <- colSums(nend_off[2], na.rm = T) ## total amount of "endemic" genera during the Ordovician offshore

nend_on[,2] <- nend_on[,2]*100/colSums(nend_on[2], na.rm = T)
nend_on <- nend_on %>% filter(!Group.1 %in% c("","Problematica"))

nend_off[,2] <- nend_off[,2]*100/colSums(nend_off[2], na.rm = T)
nend_off <- nend_off %>% filter(!Group.1 %in% c("","Problematica"))

colnames(nend_on)[1] <- "Phylum"
colnames(nend_on)[2] <- "percent"
colnames(nend_off)[1] <- "Phylum"
colnames(nend_off)[2] <- "percent"


end_on <- arrange(nend_on,desc(nend_on[,2]))
end_on$group <- ifelse(end_on[,2]>1.5, 1, 2)
end_on1 <- end_on[end_on$group==1,]

end_on1$Phylum <- factor(end_on1$Phylum,
                        levels = c(levels(end_on1$Phylum),"Others"))


end_on1[nrow(end_on1)+1,1] <- c("Others")
end_on1[nrow(end_on1),2] <- c(100-colSums(end_on1[2], na.rm = T))
phyla_on <- end_on1[,1:2]

end_off <- arrange(nend_off,desc(nend_off[,2]))
end_off$group <- ifelse(end_off[,2]>1.5, 1, 2)
end_off1 <- end_off[end_off$group==1,]

end_off1$Phylum <- factor(end_off1$Phylum,
                        levels = c(levels(end_off1$Phylum),"Others"))

end_off1[nrow(end_off1)+1,1] <- c("Others")
end_off1[nrow(end_off1),2] <- c(100-colSums(end_off1[2], na.rm = T))
phyla_off <- end_off1[,1:2]

df_phyla_onoff <- rbind.data.frame(cbind.data.frame(phyla_on,"id"="onshore"), cbind.data.frame(phyla_off,"id"="offshore"))


ggplot(df_phyla_onoff, aes(x="", y=percent, fill=Phylum)) + 
  geom_bar(stat="identity", colour = "white") +
  facet_wrap(~id)+ 
  theme_minimal() +
  scale_fill_grey(start = 0, end = 0.9) +
  coord_polar("y", start=0) +
  ggtitle("Phyla represented by Genera offshore/onshore")


#########################################################################################
#########################################################################################
#########################################################################################
## model comparison
## Table S3

rm(list=ls())

genus <- read.csv("PBDB_Ord_1.csv", sep = ",", header=T)
genus <- filter(genus, grepl("genus", accepted_rank))


########MARK MODELS
Phi.time <- list(formula=~time)
Gamma.time <- list(formula=~time)
p.time <- list(formula=~time)
L.time <- list(formula=~time)
Phi.time.env <- list(formula=~time+environment)
p.time.env <- list(formula=~time+environment)
Gamma.time.env <- list(formula=~time+environment)
L.time.env <- list(formula=~time+environment)

Phi.Time.env <- list(formula=~time+environment+time*environment)
p.Time.env <- list(formula=~time+environment+time*environment)
Gamma.Time.env <- list(formula=~time+environment+time*environment)
L.Time.env <- list(formula=~time+environment+time*environment)

Phi.const <- list(formula=~1)
p.const <- list(formula=~1)
Gamma.const <- list(formula=~1)
L.const <- list(formula=~1)

Phi.env <- list(formula=~environment)
p.env <- list(formula=~environment)
Gamma.env <- list(formula=~environment)
L.env <- list(formula=~environment)

Phi.Time.env <- list(formula=~time+environment+time*environment)
p.Time.env <- list(formula=~time+environment+time*environment)
pent.Time.env <- list(formula=~time+environment+time*environment)


##############
#####################################################################################################
## onshore
## offshore

gen_on <- genus %>% filter(environment %in% c("foreshore","coastal indet.",
                                              "marginal marine indet.",
                                              "carbonate indet.","peritidal",
                                              "sand shoal", "shoreface"))



gen_on_diff <- cbind(gen_on[,c("occurrence_no","accepted_name","phylum","max_ma","min_ma","geoplate")],
                     "diff"=gen_on$max_ma-gen_on$min_ma)


gen_off <- genus %>% filter(environment %in% c("offshore indet.","offshore","offshore ramp","offshore shelf"))


gen_off_diff <- cbind(gen_off[,c("occurrence_no","accepted_name","phylum","max_ma","min_ma","geoplate")],
                      "diff"=gen_off$max_ma-gen_off$min_ma)


gen_else <- genus %>% filter(!environment %in% c("offshore indet.","offshore","offshore ramp","offshore shelf",
                                                 "foreshore","coastal indet.",
                                                 "marginal marine indet.",
                                                 "carbonate indet.","peritidal",
                                                 "sand shoal", "shoreface"))


gen_else_diff <- cbind(gen_else[,c("occurrence_no","accepted_name","phylum","max_ma","min_ma","geoplate")],
                       "diff"=gen_else$max_ma-gen_else$min_ma)




## two vectors describing bases and names of stages for later assignment to Stages (command cut)
base_SS <- c(509,497,485.4,477.7,470,467.3,458.4,453,445.2,443.8,440.8,438.5)
name_SS<-c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2")

gen_on_val <- gen_on_diff%>% group_by(1:n()) %>% mutate(Ma=ifelse(diff<=12,runif(1,min_ma,max_ma),NA))
gen_on_val_SS <- as.data.frame(gen_on_val %>% group_by(Ma) %>% mutate(SS=cut(Ma, 
                                                                             breaks = base_SS, 
                                                                             labels = rev(name_SS),
                                                                             right=FALSE)))
gen_on_val_SS$SS <- factor(gen_on_val_SS$SS,
                           levels = c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2"))
melt_on <- melt(gen_on_val_SS, id=c("accepted_name", "SS"), na.rm = TRUE)
melt_on <- na.omit(melt_on)
cast_on <- cast(melt_on, accepted_name~SS, length)


#####
gen_off_val <- gen_off_diff%>% group_by(1:n()) %>% mutate(Ma=ifelse(diff<=12,runif(1,min_ma,max_ma),NA))
gen_off_val_SS <- as.data.frame(gen_off_val %>% group_by(Ma) %>% mutate(SS=cut(Ma, 
                                                                               breaks = base_SS, 
                                                                               labels = rev(name_SS),
                                                                               right=FALSE)))
gen_off_val_SS$SS <- factor(gen_off_val_SS$SS,
                            levels = c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2"))

melt_off <- melt(gen_off_val_SS, id=c("accepted_name", "SS"), na.rm = TRUE)
melt_off <- na.omit(melt_off)
cast_off <- cast(melt_off, accepted_name~SS, length)

#####
gen_else_val <- gen_else_diff%>% group_by(1:n()) %>% mutate(Ma=ifelse(diff<=12,runif(1,min_ma,max_ma),NA))
gen_else_val_SS <- as.data.frame(gen_else_val %>% group_by(Ma) %>% mutate(SS=cut(Ma,
                                                                                 breaks = base_SS, 
                                                                                 labels = rev(name_SS),
                                                                                 right=FALSE)))
gen_else_val_SS$SS <- factor(gen_else_val_SS$SS,
                             levels = c("C3","C4","Tr", "Fl","Dp","Dw","Sa","Ka","Hi","Si1","Si2"))

melt_else <- melt(gen_else_val_SS, id=c("accepted_name", "SS"), na.rm = TRUE)
melt_else <- na.omit(melt_else)
cast_else <- cast(melt_else, accepted_name~SS, length)


# #
##########
## finding overlapping genera
overlap <- as.factor(Reduce(intersect, list(cast_on$accepted_name, cast_off$accepted_name)))
overlap_onelse <- as.factor(Reduce(intersect,list(cast_on$accepted_name, cast_else$accepted_name)))
overlap_offelse <-   as.factor(Reduce(intersect, list(cast_off$accepted_name, cast_else$accepted_name)))

## genera that occur in all three groups
overlap_ALL <- as.factor(Reduce(intersect, list(cast_on$accepted_name, cast_off$accepted_name, cast_else$accepted_name)))

###########
## filtering and creating an input file for all genera that occur in ONSHORE areas ONLY
inp_on <- cast_on %>% filter(!accepted_name %in% overlap)
inp_on <- inp_on %>% filter(!accepted_name %in% overlap_onelse)
inp_on <- as.matrix(inp_on[,2:12])
inp_on <- ifelse(inp_on >=1, 1, 0)
data_on1 <- inp_on                                                    ## for later check of distribution of observations
inp_on <- unite(as.data.frame(inp_on), "ch", c(1:11), sep = "")
inp_on <- cbind(inp_on, environment = "onshore", ";") ### FOR FINAL INPUT DATASET

#####
## filtering and creating an input file for all genera the occur in OFFSHORE areas ONLY
inp_off <- cast_off %>% filter(!accepted_name %in% overlap)
inp_off <- cast_off %>% filter(!accepted_name %in% overlap_offelse)
inp_off <- as.matrix(inp_off[,2:12])
inp_off <- ifelse(inp_off >=1, 1, 0)
data_off1 <- inp_off                                                   ## for later check of distribution of observations
inp_off <- unite(as.data.frame(inp_off), "ch", c(1:11), sep = "")
inp_off <- cbind(inp_off, environment = "offshore", ";") ### FOR FINAL INPUT DATASET

#####
## filtering ALL genera that occur both onshore and offshore, but considering only data from onshore 
on_OL <- cast_on %>% filter(!accepted_name %in% overlap_onelse)
on_OL <- on_OL %>% filter(accepted_name %in% overlap)

## filtering ALL genera that occur both onshore and offshore, but considering only data from offshore
off_OL <- cast_off %>% filter(!accepted_name %in% overlap_offelse)
off_OL <- off_OL %>% filter(accepted_name %in% overlap)

#####
## all genera from overlap, which were just extracted and occur only in onshore areas
## 1 for presence and 0 for absence
inp_on_ol <- as.matrix(on_OL[,2:12])
inp_on_ol <- ifelse(inp_on_ol >=1, 1, 0)

## all genera from overlap, which only occur in offshore areas
## 2 for presence and 0 for absence
inp_off_ol <- as.matrix(off_OL[,2:12])
inp_off_ol <- ifelse(inp_off_ol >=1, 2, 0)

## following matrix adds up the previous two data.frames
## 1 = presence onshore
## 2 = presence offshore
## 3 = presence onshore AND offshore
inp_onoff <- as.matrix(inp_off_ol)+as.matrix(inp_on_ol) ### FOR MIG-DATA

#####
## filtering ALL genera that occur both onshore and anywhere else, but considering only data from onshore
onelse_OL <- cast_on %>% filter(!accepted_name %in% overlap)
onelse_OL <- onelse_OL %>% filter(accepted_name %in% overlap_onelse)

## filtering ALL genera that occur both onshore and anywhere else, but considering only data from anywhere else
elseon_OL <- cast_else %>% filter(!accepted_name %in% overlap)
elseon_OL <- elseon_OL %>% filter(accepted_name %in% overlap_onelse)  

#####
## data extracted from onshore gets assigned 1 for prensence (as above)
inp_onelse_ol <- as.matrix(onelse_OL[,2:12])
inp_onelse_ol <- ifelse(inp_onelse_ol >=1, 1, 0)

## data extracted from anywhere else gets assigned 10 for presence
inp_elseon_ol <- as.matrix(elseon_OL[,2:12])
inp_elseon_ol <- ifelse(inp_elseon_ol >=1, 10, 0)

## following matrix adds up on and onelse matrices
## 1 = presence onshore
## 10 = presence anywhere else but onshore and offshore
## 11 = presence onshore and anywhere else
## 0 = absence
inp_onelse <- as.matrix(inp_onelse_ol)+as.matrix(inp_elseon_ol) ### FOR MIG-DATA


#####
## filtering ALL genera that occur both offshore and anywhere else, but considering only data from offshore
offelse_OL <- cast_off %>% filter(!accepted_name %in% overlap)
offelse_OL <- offelse_OL %>% filter(accepted_name %in% overlap_offelse)

## filterin ALL genera that occur both offshore and anywhere else, but considering only data from anywhere else
elseoff_OL <- cast_else %>% filter(!accepted_name %in% overlap)
elseoff_OL <- elseoff_OL %>% filter(accepted_name %in% overlap_offelse)  

## data extracted from offshore gets assigned 2 for presence
inp_offelse_ol <- as.matrix(offelse_OL[,2:12])
inp_offelse_ol <- ifelse(inp_offelse_ol >=1, 2, 0)

## data extracted from anywhere else gets assigned 10 for presence
inp_elseoff_ol <- as.matrix(elseoff_OL[,2:12])
inp_elseoff_ol <- ifelse(inp_elseoff_ol >=1, 10, 0)

## following matrix adds up off and offelse matrices
## 2 = presence offshore
## 10 = presence anywhere else but onshore and offshore
## 12 = presence offshore and anywhere else 
## 0 = absence
inp_offelse <- as.matrix(inp_elseoff_ol)+as.matrix(inp_offelse_ol) ### FOR MIG-DATA

#####
## filtering ALL genera that occur both onshore, offshore AND anywhere else
## but considering only data from onshore first
onoffelse_OL <- cast_on %>% filter(accepted_name %in% overlap_ALL)
inp_onoffelse_ol <- as.matrix(onoffelse_OL[,2:12])
inp_onoffelse_ol <- ifelse(inp_onoffelse_ol >=1,1,0)

## considering only data from offshore
offonelse_OL <- cast_off %>% filter(accepted_name %in% overlap_ALL)
inp_offonelse_ol <- as.matrix(offonelse_OL[,2:12])
inp_offonelse_ol <- ifelse(inp_offonelse_ol >=1,2,0)

## considering only data from environments outside onshore offshore
elseonoff_OL <- cast_else %>% filter(accepted_name %in% overlap_ALL)
inp_elseonoff_ol <- as.matrix(elseonoff_OL[,2:12])
inp_elseonoff_ol <- ifelse(inp_elseonoff_ol >=1,10,0)

inp_onoffelse <- as.matrix(inp_onoffelse_ol)+as.matrix(inp_offonelse_ol)+as.matrix(inp_elseonoff_ol) ### FOR MIG-DATA


## binding all presence/absence matrices w/codes together
## substitution of numbers with places of occurrences
test_mig <- rbind.data.frame(inp_onoff,inp_onelse,inp_offelse,inp_onoffelse) ### EVERYTHING FOR MIG-DATA GOES IN HERE
test_mig[test_mig == 3] <- "on+off"
test_mig[test_mig == 2] <- "off"
test_mig[test_mig == 1] <- "on"
test_mig[test_mig == 11] <- "on+else"
test_mig[test_mig == 12] <- "off+else"
test_mig[test_mig == 13] <- "on+off+else"
test_mig[test_mig == 10] <- "else"
test_mig[test_mig == 0] <- NA


test <- as.data.frame(test_mig)

first <-   setNames(data.frame(t(apply(test[1:11], 1, function(x) {
  ind <- which(!is.na((x)))[1]
  c(ind, x[ind])
}))),
c("Stage", "first"))

test_first <- cbind(test, "first"=first$first)


extract_inp_first <- filter(test_first, first=="on" | first=="off")

extract_inp_first$environment <- ifelse(extract_inp_first$first=="on", "onshore", "offshore")
inp_mig <- ifelse(!is.na(extract_inp_first[1:11]),1,0)
data_onoff <- cbind.data.frame(inp_mig, environment = extract_inp_first$environment)
data_on2 <- filter(data_onoff, environment=="onshore")          ## for later check of distribution of observations
data_off2 <- filter(data_onoff, environment=="offshore")        ## for later check of distribution of observations
inp_mig <- unite(as.data.frame(inp_mig),"ch",c(1:11),sep="")
inp_mig <- cbind(inp_mig, environment = extract_inp_first$environment, ";") ### FOR FINAL INPUT DATASET


inp <- rbind.data.frame(inp_on,inp_off,inp_mig)



##############################################################################
##############################################################################
##############################################################################
##  model comparison for Pradsen
proc.pradsen  <- process.data(inp, model= "Pradsen", groups = "environment")

time <- mark(proc.pradsen,
             model.parameters=list(Phi=Phi.time, p=p.time, Gamma=Gamma.time))

time.const <- mark(proc.pradsen,
                   model.parameters = list(Phi=Phi.const, p=p.const, Gamma=Gamma.const))

phi.env.p.const <- mark(proc.pradsen,
                        model.parameters = list(Phi=Phi.env, p=p.const, Gamma=Gamma.const))

phi.env.p.plate <- mark(proc.pradsen,
                        model.parameters = list(Phi=Phi.env, p=p.env, Gamma=Gamma.const))

phi.env.p.time <- mark(proc.pradsen,
                       model.parameters = list(Phi=Phi.env, p=p.time, Gamma=Gamma.const))

phi.env.p.time.gamma.time <- mark(proc.pradsen,
                                  model.parameters = list(Phi=Phi.env, p=p.time, Gamma=Gamma.time))

phi.env.p.env.gamma.time <- mark(proc.pradsen,
                                 model.parameters = list(Phi=Phi.env, p=p.env, Gamma=Gamma.time))

phi.te.p.te.gamma.t <- mark(proc.pradsen,
                            model.parameters = list(Phi=Phi.time.env, p=p.time.env, Gamma=Gamma.time))

phi.te.p.t.gamma.te <- mark(proc.pradsen,
                            model.parameters = list(Phi=Phi.time.env, p=p.time, Gamma=Gamma.time.env))

phi.t.p.te.g.te <- mark(proc.pradsen,
                        model.parameters = list(Phi=Phi.time, p=p.time.env, Gamma=Gamma.time.env))

time.env <- mark(proc.pradsen,
                 model.parameters = list(Phi=Phi.time.env, p=p.time.env, Gamma=Gamma.time.env))

Time.env <- mark(proc.pradsen,
                 model.parameters = list(Phi=Phi.Time.env, p=p.Time.env, Gamma=Gamma.Time.env))

phi.Te.p.te.gamma.Te <- mark(proc.pradsen,
                             model.parameters = list(Phi=Phi.Time.env, p=p.time.env, Gamma=Gamma.Time.env))

phi.te.p.Te.gamma.Te <- mark(proc.pradsen,
                             model.parameters = list(Phi=Phi.time.env, p=p.Time.env, Gamma=Gamma.Time.env))

phi.Te.p.Te.gamma.te <- mark(proc.pradsen,
                             model.parameters = list(Phi=Phi.Time.env, p=p.Time.env, Gamma=Gamma.time.env))

env <- mark(proc.pradsen, model.parameters = list(Phi=Phi.env, p=p.env, Gamma= Gamma.env))

## compare all models in following table
comparison <- collect.models(type = "Pradsen")
comparison


