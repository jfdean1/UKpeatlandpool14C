######################################################################
##### Data analysis and presentation for the manuscript ##############
##### JF Dean et al., Peatland pools are tightly coupled to the ######
##### contemporary carbon cycle. #####################################
######################################################################

##### All data files used below are included as individual .csv files in the associated Github repository #####

# packages #
library(ggplot2)
library(rootSolve)
library(dplyr)
library(RColorBrewer)
library(forcats)
library(ggpubr)
library(cowplot)
library(gridExtra)
library(conover.test)


##### 14C exponential age model #####
# adapted from Dean et al. (2018) https://doi.org/10.1088/1748-9326/aaa1fe

## Reading data
atmdata<-read.table("Atmospheric14CO2_2014-10ky.csv",header=T,sep=",")

## make a continuous function of data (function calls for any year)
datfun<-splinefun(atmdata[,1],atmdata[,3])

## function to calculate F14C with exponential distribution (fraction of given year in a given sample)
cdiste<-function(Y,k,MY){
  datfun(Y)*dexp(MY-Y,rate=k)
}

## function to check if the distribution integrates to 1 for chosen time interval
tdiste<-function(Y,k,MY){
  dexp(MY-Y,rate=k)
}

## sequence of K values of exponential distribution
kse<- exp(seq(-10,log(2),0.01))
coute<-vector()
toute<-vector()

## loop over all K values
for(i in 1:length(kse)){
  ##try to integrate over time interval of input data #NOTE I ADDED THE "stop.on.error = FALSE" IN THE INTEGRATE FUNCTION OTHERWISE LOTS OF ROUNDING ERRORS
  coutt<-try(integrate(cdiste,min(atmdata[,1]),max(atmdata[,1]),k=kse[i],MY=max(atmdata[,1]),subdivisions=10*length(atmdata[,1]),
                       stop.on.error = FALSE)[[1]])
  ##if integration succeeded save data
  if(is.numeric(coutt))coute[i]<-coutt
  toutt<-try(integrate(tdiste,min(atmdata[,1]),max(atmdata[,1]),k=kse[i],MY=max(atmdata[,1]),subdivisions=10*length(atmdata[,1]),
                       stop.on.error = FALSE)[[1]])
  if(is.numeric(toutt))toute[i]<-toutt
}

## plot
plot( 1/kse[toute>0.99]/log(2),coute[toute>0.99],typ="l",lwd=2,xlab="Median age",ylab="F14C",log="x")
plot( 1/kse[toute>0.99]/log(2),coute[toute>0.99],typ="l",lwd=2,xlab="Median age",ylab="F14C")

## read data
data<-read.table("Pools14Cdata.csv",sep=",", header=T)
datums<-strptime(data[,2],"%d/%m/%Y",tz="UTC")

# #kse as function of concentration
agefune<-splinefun(kse[toute>0.99],coute[toute>0.99])

## kse as function of concentration (log scale)
agefunelog<-splinefun(log(kse[toute>0.99]),coute[toute>0.99])

## function that calculates difference between a concentration belonging to a certain age (1/ks) and the sample
agefun2e<-function(kv,cm){
  return(agefune(kv)-cm)
}

agefun2elog<-function(kv,cm){
  return(agefunelog(kv)-cm)
}

## median age (average) ##log scale
## for each sample (data[,5]) find the ages that have a concentration equal to the measured concentration
outexpl<-sapply(data[,5],FUN=function(x)sort((1/exp(uniroot.all(agefun2elog,log(c(exp(-10),max(kse[toute>0.99]))),cm=x)))*log(2)))

## median age (average)
outexp<-sapply(data[,5],FUN=function(x)sort((1/uniroot.all(agefun2e,c(0,max(kse[toute>0.99])),cm=x))*log(2)))
## median age + 1 SD
outexpp<-sapply(data[,5]+data[,6],FUN=function(x)sort((1/uniroot.all(agefun2e,c(0,max(kse[toute>0.99])),cm=x))*log(2)))
## median age - 1 SD
outexpm<-sapply(data[,5]-data[,6],FUN=function(x)sort((1/uniroot.all(agefun2e,c(0,max(kse[toute>0.99])),cm=x))*log(2)))

## fraction < 10 Y
outexp2<-sapply(data[,5],FUN=function(x){
  ks<-uniroot.all(agefun2e,c(0,max(kse[toute>0.99])),cm=x)
  return(sort((1-exp(-ks*10))*100))
})

## fraction < 50 Y
outexp3<-sapply(data[,5],FUN=function(x){
  ks<-uniroot.all(agefun2e,c(0,max(kse[toute>0.99])),cm=x)
  return(sort((1-exp(-ks*50))*100))
})

## fraction < 100 Y
outexp4<-sapply(data[,5],FUN=function(x){
  ks<-uniroot.all(agefun2e,c(0,max(kse[toute>0.99])),cm=x)
  return(sort((1-exp(-ks*100))*100))
})

## fraction < 300 Y
outexp5<-sapply(data[,5],FUN=function(x){
  ks<-uniroot.all(agefun2e,c(0,max(kse[toute>0.99])),cm=x)
  return(sort((1-exp(-ks*300))*100))
})

## fraction < 500 Y
outexp6<-sapply(data[,5],FUN=function(x){
  ks<-uniroot.all(agefun2e,c(0,max(kse[toute>0.99])),cm=x)
  return(sort((1-exp(-ks*500))*100))
})

## fraction < 1000 Y
outexp7<-sapply(data[,5],FUN=function(x){
  ks<-uniroot.all(agefun2e,c(0,max(kse[toute>0.99])),cm=x)
  return(sort((1-exp(-ks*1000))*100))
})

## write output
zz <- file("ranges2014-10ky.dat", "w")
for( i in 1:length(outexp)){
  if(length(outexp[[i]])==1)writeLines(paste(i,datums[i],data[i,3],data[i,5],data[i,6],NA,NA,NA,round(outexp[[i]][1],2),round(outexpp[[i]][1],2),round(outexpm[[i]][1],2),round(outexp2[[i]][1],2),round(outexp3[[i]][1],2),round(outexp4[[i]][1],4),round(outexp5[[i]][1],4),round(outexp6[[i]][1],4),round(outexp7[[i]][1],4),sep=";"), zz)
  if(length(outexp[[i]])==2)writeLines(paste(i,datums[i],data[i,3],data[i,5],data[i,6],round(outexp[[i]][1],2),round(outexpp[[i]][1],2),round(outexpm[[i]][1],2),round(outexp[[i]][2],2),round(outexpp[[i]][2],2),round(outexpm[[i]][2],2),round(outexp2[[i]][1],2),round(outexp3[[i]][1],2),round(outexp4[[i]][1],4),round(outexp5[[i]][1],4),round(outexp6[[i]][1],4),round(outexp7[[i]][1],4),sep=";"), zz)
}
close(zz)

##### model data is processed into a csv then re-uploaded below #####

## Data from model output
F14C_data <- read.csv("14C_modeloutput_10ky.csv")
names(F14C_data)

## Data for correlations
All_data = read.csv("All_corr_data.csv")
names(All_data)








#### Figure 1 ####

## Figure 1A ##
# F14C by C species across all sites

# Kruskal-Wallis test, Conover-Iman posthoc
kruskal.test(F14C ~ C_species_n, data = F14C_data) # p = 2.655e-07
conover.test(F14C_data$F14C , F14C_data$C_species_n, altp=T) #excluded where n < 2

# note:
# DOC = C_species_n "a"
# POC = C_species_n "b"
# CO2 = C_species_n "c"
# CO2eb = C_species_n "d"
# CH4eb = C_species_n "e"
# CH4 = C_species_n "f"
# Sed = C_species_n "g"

# subset by C-species
F14C_data_DOC <- F14C_data[which(F14C_data$C_species_n=="a"),]
F14C_data_POC <- F14C_data[which(F14C_data$C_species_n=="b"),]
F14C_data_CO2 <- F14C_data[which(F14C_data$C_species_n=="c"),]
F14C_data_CO2eb <- F14C_data[which(F14C_data$C_species_n=="d"),]
F14C_data_CH4eb <- F14C_data[which(F14C_data$C_species_n=="e"),]
F14C_data_CH4 <- F14C_data[which(F14C_data$C_species_n=="f"),]
F14C_data_Sed <- F14C_data[which(F14C_data$C_species_n=="g"),]

# descriptive stats
median(F14C_data_DOC$F14C) # 1.0663
min(F14C_data_DOC$F14C) # 0.9690
max(F14C_data_DOC$F14C) # 1.1019

median(F14C_data_POC$F14C) # 1.0050
min(F14C_data_POC$F14C) # 0.9370
max(F14C_data_POC$F14C) # 1.0426

median(F14C_data_CO2$F14C) # 1.0294
min(F14C_data_CO2$F14C) # 0.9420
max(F14C_data_CO2$F14C) # 1.0550

median(F14C_data_CO2eb$F14C) # 1.0169

median(F14C_data_CH4eb$F14C) # 1.0279 
min(F14C_data_CH4eb$F14C) # 0.8732 
max(F14C_data_CH4eb$F14C) # 1.0758

median(F14C_data_CH4$F14C) # 1.0523

median(F14C_data_Sed$F14C) # 0.7608
min(F14C_data_Sed$F14C) # 0.7371
max(F14C_data_Sed$F14C) # 0.9274


# get max F14C values for each C species
F14CbySpecies_All.MaxValue = as.vector(c(max(F14C_data_DOC$F14C),
                                         max(F14C_data_POC$F14C),
                                         max(F14C_data_CO2$F14C),
                                         max(F14C_data_CO2eb$F14C),
                                         max(F14C_data_CH4eb$F14C),
                                         max(F14C_data_CH4$F14C),
                                         max(F14C_data_Sed$F14C)))

F14CbySpecies_All.SigLett_groups = c("a","c","b"," ","bc"," ","c")

# build data frame of values and sig-letts for plot
F14CbySpecies_All.MaxValue.df = data.frame(F14CbySpecies_All.SigLett_groups,
                                           F14CbySpecies_All.MaxValue,
                                           c("a","b","c","d","e","f","g"))
F14CbySpecies_All.MaxValue.df_names <- c("SigLett", "MaxValue", "C_species_n")
names(F14CbySpecies_All.MaxValue.df) <- F14CbySpecies_All.MaxValue.df_names

Fig_F14CbySpecies_All = ggplot(F14C_data, aes(x = C_species_n, y = F14C, fill = C_species_n)) +
  geom_boxplot(outlier.shape=NA)  +
  geom_jitter(width = 0.2, height = 0, shape=16, color="black", size=1.2, alpha=0.7) +
  scale_fill_manual(values=c("#FFFFFF", "#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF")) +
  scale_y_continuous(limits = c(0.7,1.12), breaks = c(0.7,0.8,0.9,1,1.1)) +
  scale_x_discrete(labels=c('DOC','POC',expression('CO'[2]),expression('CO'[2]*' eb.'),expression('CH'[4]*' eb.'),expression('CH'[4]),'Sed.')) +
  xlab('C forms') + ylab(expression('F'^14*'C (fraction modern)')) +
  annotate("text", x=1, y=0.95, size=2.5, label= "253 yBP") +
  annotate("text", x=2, y=0.92, size=2.5, label= "525 yBP") +
  annotate("text", x=3, y=0.9275, size=2.5, label= "480 yBP") +
  annotate("text", x=5, y=0.855, size=2.5, label= "1090 yBP") +
  annotate("text", x=7, y=0.72, size=2.5, label= "2451 yBP") +
  annotate("text", x=6.9, y=1.01, size=2.25, label=expression(italic('F'^14*'C modern')), color = "snow3") +
  geom_text(data = F14CbySpecies_All.MaxValue.df,
            aes(x = c(1:7), y = 0.01 + MaxValue,  label = SigLett), vjust=0 , size = 2.5) +
  geom_hline(yintercept=1, linetype="dashed", color = "snow3") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position = "none")
Fig_F14CbySpecies_All


## Figure 1B ##
# Age contributions to mixtures by C species

## plotting age contributions from distributions

# first calculate values...
# 0-9 years
p0_9y <- F14C_data[,23]

# 10-49 years
p10_49y <- (F14C_data[,24] - F14C_data[,23])

# 50-99 years
p50_99y <- (F14C_data[,25] - F14C_data[,24])

# 100-299 years
p100_299y <- (F14C_data[,26] - F14C_data[,25])

# 300-499 years
p300_499y <- (F14C_data[,27] - F14C_data[,26])

# 500-999 years
p500_999y <- (F14C_data[,28] - F14C_data[,27])

# 1000+ years
# create function which is 100-value, value being % > p_1000y
minus <- function(value) {
  x <- 100-(value)
  return(x)
}
# calculate
p1000y_plus <- sapply(F14C_data[,28], minus)

# add to dataframe
F14C_data$p0_9y <- c(p0_9y)
F14C_data$p10_49y <- c(p10_49y)
F14C_data$p50_99y <- c(p50_99y)
F14C_data$p100_299y <- c(p100_299y)
F14C_data$p300_499y <- c(p300_499y)
F14C_data$p500_999y <- c(p500_999y)
F14C_data$p1000y_plus <- c(p1000y_plus)

#... then resample new dataframe

## Subset for C Species

F14C_data_DOC <- F14C_data[which(F14C_data$C_species_n=="a"),]
F14C_data_POC <- F14C_data[which(F14C_data$C_species_n=="b"),]
F14C_data_CO2 <- F14C_data[which(F14C_data$C_species_n=="c"),]
F14C_data_CO2eb <- F14C_data[which(F14C_data$C_species_n=="d"),]
F14C_data_CH4eb <- F14C_data[which(F14C_data$C_species_n=="e"),]
F14C_data_CH4 <- F14C_data[which(F14C_data$C_species_n=="f"),]
F14C_data_Sed <- F14C_data[which(F14C_data$C_species_n=="g"),]

DOC_ageclasses_calc <- F14C_data_DOC[,c(29:35)]
POC_ageclasses_calc <- F14C_data_POC[,c(29:35)]
CO2_ageclasses_calc <- F14C_data_CO2[,c(29:35)]
CO2eb_ageclasses_calc <- F14C_data_CO2eb[,c(29:35)]
CH4eb_ageclasses_calc <- F14C_data_CH4eb[,c(29:35)]
CH4_ageclasses_calc <- F14C_data_CH4[,c(29:35)]
Sed_ageclasses_calc <- F14C_data_Sed[,c(29:35)]

# Pull apart this split for a table (i.e. potential contributions for different aged soil layers)
A1 <- as.vector(colMeans(DOC_ageclasses_calc[sapply(DOC_ageclasses_calc, is.numeric)]))
B1 <- as.vector(colMeans(POC_ageclasses_calc[sapply(POC_ageclasses_calc, is.numeric)]))
C1 <- as.vector(colMeans(CO2_ageclasses_calc[sapply(CO2_ageclasses_calc, is.numeric)]))
D1 <- as.vector(colMeans(CO2eb_ageclasses_calc[sapply(CO2eb_ageclasses_calc, is.numeric)]))
E1 <- as.vector(colMeans(CH4eb_ageclasses_calc[sapply(CH4eb_ageclasses_calc, is.numeric)]))
F1 <- as.vector(colMeans(CH4_ageclasses_calc[sapply(CH4_ageclasses_calc, is.numeric)]))
G1 <- as.vector(colMeans(Sed_ageclasses_calc[sapply(Sed_ageclasses_calc, is.numeric)]))

Lab1 <- as.vector(c("0-9 yBSD","10-49 yBSD","50-99 yBSD","100-299 yBSD","300-499 yBSD","500-999 yBSD",">1000 yBSD",
                    "0-9 yBSD","10-49 yBSD","50-99 yBSD","100-299 yBSD","300-499 yBSD","500-999 yBSD",">1000 yBSD",
                    "0-9 yBSD","10-49 yBSD","50-99 yBSD","100-299 yBSD","300-499 yBSD","500-999 yBSD",">1000 yBSD",
                    "0-9 yBSD","10-49 yBSD","50-99 yBSD","100-299 yBSD","300-499 yBSD","500-999 yBSD",">1000 yBSD",
                    "0-9 yBSD","10-49 yBSD","50-99 yBSD","100-299 yBSD","300-499 yBSD","500-999 yBSD",">1000 yBSD",
                    "0-9 yBSD","10-49 yBSD","50-99 yBSD","100-299 yBSD","300-499 yBSD","500-999 yBSD",">1000 yBSD",
                    "0-9 yBSD","10-49 yBSD","50-99 yBSD","100-299 yBSD","300-499 yBSD","500-999 yBSD",">1000 yBSD"))

C_species_n <- as.vector(c("a","a","a","a","a","a","a",
                         "b","b","b","b","b","b","b",
                         "c","c","c","c","c","c","c",
                         "d","d","d","d","d","d","d",
                         "e","e","e","e","e","e","e",
                         "f","f","f","f","f","f","f",
                         "g","g","g","g","g","g","g"))
age_sum <- as.vector(c(A1, B1, C1, D1, E1, F1, G1))
age_class_sum <- data.frame(Lab1, C_species_n, age_sum)
age_class_sum_names <- c("Age_class", "C_species_n", "age_sum")

age_class_sum <- data.frame(Lab1, A1, B1, C1, D1, E1, F1, G1)
age_class_sum_names <- c("Age_class", "DOC", "POC", "CO2", "CO2_eb", "CH4_eb", "CH4", "Sed")
names(age_class_sum) <- age_class_sum_names

Fig_age_class_barplot = ggplot(age_class_sum, aes(fill=factor(Lab1, unique(Lab1)), x=C_species_n, y=age_sum)) + 
  scale_fill_manual(values=c("#C7EAE5","#F5F5F5","#F6E8C3","#DFC27D","#BF812D","#8C510A","#543005"), name="Age class") +
  scale_x_discrete(labels=c('DOC','POC',expression('CO'[2]),expression('CO'[2]*' eb.'),expression('CH'[4]*' eb.'),expression('CH'[4]),'Sed.')) +
  xlab('C forms') + ylab('Predicted contribution to age mixture (%)') +
  geom_bar(position="stack", stat="identity", color="black") + labs(color=NULL) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Fig_age_class_barplot


## combine plots ##
ggdraw() +
  draw_plot(Fig_F14CbySpecies_All, x = 0, y = 0, width = .42, height = 1) +
  draw_plot(Fig_age_class_barplot, x = .45, y = 0, width = .55, height = 1) +
  draw_plot_label(label = c("a", "b"), size = 15,
                  x = c(0, 0.445), y = c(1, 1))
ggsave("Fig_1.png", width = 25, height = 10, units = c("cm"), dpi = 600)







##### Figure 2 #####
#  F14C by site (excl. sed)

## remove sediment values
F14C_data_no_sed <- F14C_data[-c(50,51,52),]

plot_F14CbyLocation <- ggplot(F14C_data_no_sed, aes(x=Location_n, y=F14C, fill = Location_n)) + 
  geom_boxplot(outlier.shape=NA)

# Kruskal-Wallis test, Conover-Iman posthoc
kruskal.test(F14C ~ Location_n, data = F14C_data_no_sed) # p = 0.009112
conover.test(F14C_data_no_sed$F14C, F14C_data_no_sed$Location_n, altp=T) 

# note:
# Cross Lochs = Location_n "a"
# Loch Leir = Location_n "b"
# Munsary = Location_n "c"
# Silver Flowe = Location_n "d"
# Garron Plateau = Location_n "e"
# Slieveanorra = Location_n "f"

#subset by location
F14C_data_CL <- F14C_data_no_sed[which(F14C_data_no_sed$Location_n=="a"),]
F14C_data_LL <- F14C_data_no_sed[which(F14C_data_no_sed$Location_n=="b"),]
F14C_data_MUN <- F14C_data_no_sed[which(F14C_data_no_sed$Location_n=="c"),]
F14C_data_SF <- F14C_data_no_sed[which(F14C_data_no_sed$Location_n=="d"),]
F14C_data_GP <- F14C_data_no_sed[which(F14C_data_no_sed$Location_n=="e"),]
F14C_data_SL <- F14C_data_no_sed[which(F14C_data_no_sed$Location_n=="f"),]

# descriptive stats
median(F14C_data_CL$F14C) # 1.0218
min(F14C_data_CL$F14C) # 0.8732
max(F14C_data_CL$F14C) # 1.0962

median(F14C_data_LL$F14C) # 1.0187
min(F14C_data_LL$F14C) # 0.9370
max(F14C_data_LL$F14C) # 1.0807

median(F14C_data_MUN$F14C) # 1.0696
min(F14C_data_MUN$F14C) # 1.0110
max(F14C_data_MUN$F14C) # 1.0799

median(F14C_data_SF$F14C) # 1.0448
min(F14C_data_SF$F14C) #  1.0231
max(F14C_data_SF$F14C) # 1.0660

median(F14C_data_GP$F14C) # 1.0671
min(F14C_data_GP$F14C) # 1.0326
max(F14C_data_GP$F14C) # 1.0939

median(F14C_data_SL$F14C) # 1.0378
min(F14C_data_SL$F14C) # 1.0305
max(F14C_data_SL$F14C) # 1.1019

# get max F14C value of each location
F14CbyLocation.MaxValue = as.vector(c(max(F14C_data_CL$F14C),
                                      max(F14C_data_LL$F14C),
                                      max(F14C_data_MUN$F14C),
                                      max(F14C_data_SF$F14C),
                                      max(F14C_data_GP$F14C),
                                      max(F14C_data_SL$F14C)))

F14CbySpecies_All.SigLett_groups = c("bc","c","abc","abc","a","ac")

# build data frame of values and sigletts for plot
F14CbyLocation_All.MaxValue.df = data.frame(F14CbySpecies_All.SigLett_groups,
                                            F14CbyLocation.MaxValue,
                                           c("a","b","c","d","e","f"))
F14CbyLocation_All.MaxValue.df_names <- c("SigLett", "MaxValue", "Location_n")
names(F14CbyLocation_All.MaxValue.df) <- F14CbyLocation_All.MaxValue.df_names

Fig_F14CbyLocation <- plot_F14CbyLocation + geom_jitter(width = 0.2, height = 0, shape=16, color="black", size=1.2, alpha=0.7) +
  scale_fill_manual(values=c("#FFFFFF", "#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF")) +
  scale_y_continuous(limits = c(0.85,1.12), breaks = c(0.85,0.9,0.95,1,1.05,1.1)) +
  scale_x_discrete(labels=c("Cross Lochs", "Loch Leir", "Munsary", "Silver Flowe", "Garron Plateau", "Slieveanorra")) +
  xlab('Location') + ylab(expression('F'^14*'C (fraction modern)')) +
  annotate("text", x=1, y=0.862, size=2.5, label= "1090 yBP") +
  annotate("text", x=2, y=0.926, size=2.5, label= "525 yBP") +
  annotate("text", x=6, y=1.006, size=2.25, label=expression(italic('F'^14*'C modern')), color = "snow3") +
  geom_hline(yintercept=1, linetype="dashed", color = "snow3") +
  geom_text(data = F14CbyLocation_All.MaxValue.df,
            aes(x = c(1:6), y = 0.01 + MaxValue,  label = SigLett), vjust=0 , size = 2.5) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(legend.position = "none")
Fig_F14CbyLocation

ggsave("Fig_2.png", width = 10, height = 12, units = c("cm"), dpi = 600)








##### Figure 3 #####

## Figure 3A ##
# Natural & Restoration pool type comparison F14C

plot_F14CbyPoolType <- ggplot(All_data, aes(x=C_species_n, y=F14C, fill = Pool_type)) + 
  geom_boxplot(outlier.shape=NA) + facet_wrap(~Pool_type, scale='free') +
  geom_point(position=position_jitterdodge(0.5), shape=16, color="black", size=1.2, alpha=0.7)  

## posthoc stats: unpaired two-sample Wilcoxon test (by pool type)
wilcox.test(All_data$F14C ~ All_data$Pool_type, paired = FALSE) # p = 0.1642

# selections
All_data_nat = (All_data[which(All_data$Pool_type=="Natural"),])
All_data_res = (All_data[which(All_data$Pool_type=="Restoration"),])
All_data_nat.F14Cmean = mean(All_data_nat$F14C)
All_data_res.F14Cmean = mean(All_data_res$F14C)
plot_F14CbyPoolType.means = data.frame(Pool_type = c("Natural", "Restoration"), F14C = c(All_data_nat.F14Cmean, All_data_res.F14Cmean))

DOC = All_data[which(All_data$C_species=="DOC"),]
POC = All_data[which(All_data$C_species=="POC"),]
CO2 = All_data[which(All_data$C_species=="CO2"),]

## posthoc stats: unpaired two-samples Wilcoxon test (by pool type and C species)
wilcox.test(DOC$F14C ~ DOC$Pool_type, paired = FALSE) # p = 0.03053
wilcox.test(POC$F14C ~ POC$Pool_type, paired = FALSE) # p = 1
wilcox.test(CO2$F14C ~ CO2$Pool_type, paired = FALSE) # p = 0.06355

F14CbyPoolType.SigLett <- c('b','d','c','a','d','c')
F14CbyPoolType.MaxValue.df = data.frame(F14CbyPoolType.SigLett, c(1.10,1.04,1.06,1.10,1.04,1.06),
                                        c("Natural","Natural","Natural","Restoration","Restoration","Restoration"))
F14CbyPoolType.MaxValue.df_names <- c("SigLett", "MaxValue", "Pool_type")
names(F14CbyPoolType.MaxValue.df) <- F14CbyPoolType.MaxValue.df_names

F14CbyPoolType.yBP.df = data.frame(c("253 yBP", "295 yBP", "480 yBP", "525 yBP", "175 yBP"),
                                   c(0.96, 0.955, 0.935, 0.93, 0.97),
                                   c("Natural", "Natural", "Natural", "Restoration", "Restoration"))
F14CbyPoolType.yBP.df_names = c("Age", "Y", "Pool_type")
names(F14CbyPoolType.yBP.df) = F14CbyPoolType.yBP.df_names

Fig_F14CbyPoolType <- plot_F14CbyPoolType +
  scale_fill_manual(values=c("#FFFFFF", "#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF")) +
  xlab('C forms') + ylab(expression('F'^14*'C (fraction modern)')) +
  scale_y_continuous(limits = c(0.93,1.11), breaks = c(0.9,0.95,1,1.05,1.1)) +
  scale_x_discrete(labels=c("DOC", "POC", expression('CO'[2]))) +
  geom_hline(aes(yintercept=F14C), plot_F14CbyPoolType.means, linetype="dashed", color = "black") +
  geom_hline(yintercept=1, linetype="dashed", color = "snow3") +
  geom_text(data = F14CbyPoolType.yBP.df,
            aes(x = c(1,2,3,2,3), y = Y,  label = Age), vjust=0 , size = 2.5) +
  geom_text(data = F14CbyPoolType.MaxValue.df,
            aes(x = c(1,2,3,1,2,3), y = 0.01 + MaxValue,  label = SigLett), vjust=0 , size = 2.5) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none")
Fig_F14CbyPoolType


## Figure 3B ##
# Age contributions to mixtures by pool type

## Subset for pool type
F14C_data_nat <- F14C_data[which(F14C_data$Pool_type=="Natural"),]
F14C_data_nat_DOC <- F14C_data_nat[which(F14C_data_nat$C_species_n=="a"),]
F14C_data_nat_POC <- F14C_data_nat[which(F14C_data_nat$C_species_n=="b"),]
F14C_data_nat_CO2 <- F14C_data_nat[which(F14C_data_nat$C_species_n=="c"),]
nat_DOC_ageclasses_calc <- F14C_data_nat_DOC[,c(29:35)]
nat_POC_ageclasses_calc <- F14C_data_nat_POC[,c(29:35)]
nat_CO2_ageclasses_calc <- F14C_data_nat_CO2[,c(29:35)]

F14C_data_res <- F14C_data[which(F14C_data$Pool_type=="Restoration"),]
F14C_data_res_DOC <- F14C_data_res[which(F14C_data_res$C_species_n=="a"),]
F14C_data_res_POC <- F14C_data_res[which(F14C_data_res$C_species_n=="b"),]
F14C_data_res_CO2 <- F14C_data_res[which(F14C_data_res$C_species_n=="c"),]
res_DOC_ageclasses_calc <- F14C_data_res_DOC[,c(29:35)]
res_POC_ageclasses_calc <- F14C_data_res_POC[,c(29:35)]
res_CO2_ageclasses_calc <- F14C_data_res_CO2[,c(29:35)]

# descriptive stats
median(F14C_data_nat_DOC$F14C) # 1.0554
min(F14C_data_nat_DOC$F14C) # 0.9690
max(F14C_data_nat_DOC$F14C) # 1.1019

median(F14C_data_nat_POC$F14C) # 1.0052
min(F14C_data_nat_POC$F14C) # 0.9640
max(F14C_data_nat_POC$F14C) # 1.0427

median(F14C_data_nat_CO2$F14C) # 1.0231
min(F14C_data_nat_CO2$F14C) # 0.9420
max(F14C_data_nat_CO2$F14C) # 1.0550

median(F14C_data_res_DOC$F14C) # 1.073
min(F14C_data_res_DOC$F14C) # 1.0466
max(F14C_data_res_DOC$F14C) # 1.0962

median(F14C_data_res_POC$F14C) # 1.0040
min(F14C_data_res_POC$F14C) # 0.9370
max(F14C_data_res_POC$F14C) # 1.0335

median(F14C_data_res_CO2$F14C) # 1.0342
min(F14C_data_res_CO2$F14C) # 0.9785
max(F14C_data_res_CO2$F14C) # 1.0474

# Pull subsets apart for a table (i.e. potential contributions for different aged soil layers)
A2 <- as.vector(colMeans(nat_DOC_ageclasses_calc[sapply(nat_DOC_ageclasses_calc, is.numeric)]))
B2 <- as.vector(colMeans(res_DOC_ageclasses_calc[sapply(res_DOC_ageclasses_calc, is.numeric)]))
C2 <- as.vector(colMeans(nat_POC_ageclasses_calc[sapply(nat_POC_ageclasses_calc, is.numeric)]))
D2 <- as.vector(colMeans(res_POC_ageclasses_calc[sapply(res_POC_ageclasses_calc, is.numeric)]))
E2 <- as.vector(colMeans(nat_CO2_ageclasses_calc[sapply(nat_CO2_ageclasses_calc, is.numeric)]))
F2 <- as.vector(colMeans(res_CO2_ageclasses_calc[sapply(res_CO2_ageclasses_calc, is.numeric)]))

Lab2 <- as.vector(c("0-9 yBSD","10-49 yBSD","50-99 yBSD","100-299 yBSD","300-499 yBSD","500-999 yBSD",">1000 yBSD",
                    "0-9 yBSD","10-49 yBSD","50-99 yBSD","100-299 yBSD","300-499 yBSD","500-999 yBSD",">1000 yBSD",
                    "0-9 yBSD","10-49 yBSD","50-99 yBSD","100-299 yBSD","300-499 yBSD","500-999 yBSD",">1000 yBSD",
                    "0-9 yBSD","10-49 yBSD","50-99 yBSD","100-299 yBSD","300-499 yBSD","500-999 yBSD",">1000 yBSD",
                    "0-9 yBSD","10-49 yBSD","50-99 yBSD","100-299 yBSD","300-499 yBSD","500-999 yBSD",">1000 yBSD",
                    "0-9 yBSD","10-49 yBSD","50-99 yBSD","100-299 yBSD","300-499 yBSD","500-999 yBSD",">1000 yBSD"))

C_species_n2 <- as.vector(c("a","a","a","a","a","a","a",
                            "b","b","b","b","b","b","b",
                            "c","c","c","c","c","c","c",
                            "d","d","d","d","d","d","d",
                            "e","e","e","e","e","e","e",
                            "f","f","f","f","f","f","f"))
age_sum2 <- as.vector(c(A2, B2, C2, D2, E2, F2))
age_class_sum2 <- data.frame(Lab2, C_species_n2, age_sum2)
age_class_sum_names2 <- c("Age_class", "C_species_n", "age_sum")

age_class_sum2 <- data.frame(Lab2, A2, B2, C2, D2, E2, F2)
age_class_sum_names2 <- c("Age_class", "DOC Nat", "DOC Res", "POC Nat", "POC Res", "CO2 Nat", "CO2 Res")
names(age_class_sum2) <- age_class_sum_names2

Fig_age_class_barplot2 = ggplot(age_class_sum2, aes(fill=factor(Lab2, unique(Lab2)), x=C_species_n2, y=age_sum2)) + 
  scale_fill_manual(values=c("#C7EAE5","#F5F5F5","#F6E8C3","#DFC27D","#BF812D","#8C510A","#543005"), name="Age class") +
  scale_x_discrete(labels=c('DOC (Nat.)','DOC (Res.)','POC (Nat.)','POC (Res.)',
                            expression('CO'[2]*' (Nat.)'),expression('CO'[2]*' (Res.)'))) +
  xlab('C forms (pool type)') + ylab('Predicted contribution to age mixture (%)') +
  geom_bar(position="stack", stat="identity", color="black") + labs(color=NULL) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                     axis.title.x = element_text(vjust=-0.5))
Fig_age_class_barplot2


## combine plots ##
ggdraw() +
  draw_plot(Fig_F14CbyPoolType, x = 0, y = 0, width = .43, height = 1) +
  draw_plot(Fig_age_class_barplot2, x = .46, y = 0, width = .54, height = 1) +
  draw_plot_label(label = c("a", "b"), size = 15,
                  x = c(0, 0.455), y = c(1, 1))

ggsave("Fig_3.png", width = 22, height = 10, units = c("cm"), dpi = 600)








##### Figure 4 #####
# Natural & Restoration pool type comparison - d13C

All_data = read.csv("All_corr_data.csv")

All_data$C_species_n = as.factor(as.character(All_data$C_species_n))
levels(All_data$C_species_n) <- c("DOC", "POC",  expression("CO"[2]))

fct_relevel(All_data$C_species_n)
fct_relevel(All_data$C_species_n, "DOC", "POC")

All_data$Pool_type = as.factor(as.character(All_data$Pool_type), levels=All_data$Pool_type)
levels(All_data$Pool_type) <- c("Natural", "Restoration")

plot_d13CbyPoolType2 <- ggplot(All_data, aes(x=Pool_type, y=d13C, fill = C_species_n)) + 
  geom_boxplot(outlier.shape=NA) + facet_wrap(~C_species_n, scale='free', labeller = label_parsed) +
  geom_point(position=position_jitterdodge(0.5), shape=16, color="black", size=1.2, alpha=0.7)

# Kruskal-Wallis test, Conover-Iman posthoc
kruskal.test(d13C ~ C_species, data = All_data) # p = 4.066e-13
conover.test(All_data$d13C, All_data$C_species, altp=T) 

# selections
All_data_nat = (All_data[which(All_data$Pool_type=="Natural"),])
All_data_res = (All_data[which(All_data$Pool_type=="Restoration"),])

DOC = All_data[which(All_data$C_species=="DOC"),]
POC = All_data[which(All_data$C_species=="POC"),]
CO2 = All_data[which(All_data$C_species=="CO2"),]

# descriptive stats
All_data_nat.d13Cmean = mean(All_data_nat$d13C, na.rm=T)
All_data_res.d13Cmean = mean(All_data_res$d13C, na.rm=T)
plot_d13CbyPoolType.means = data.frame(Pool_type = c("Natural", "Restoration"), d13C = c(All_data_nat.d13Cmean, All_data_res.d13Cmean))

All_data_DOC.d13Cmean = mean(DOC$d13C, na.rm=T)
min(DOC$d13C, na.rm=T) # -29.1
max(DOC$d13C, na.rm=T) # -23.4

All_data_POC.d13Cmean = mean(POC$d13C, na.rm=T)
min(POC$d13C, na.rm=T) # -29.9
max(POC$d13C, na.rm=T) # -24.3

All_data_CO2.d13Cmean = mean(CO2$d13C, na.rm=T)
min(CO2$d13C, na.rm=T) # -20.5
max(CO2$d13C, na.rm=T) # -9.4

plot_d13CbyPoolType.means2 = data.frame(C_species_n = c("a", "b", "c"), d13C = c(All_data_DOC.d13Cmean, All_data_POC.d13Cmean, All_data_CO2.d13Cmean))
plot_d13CbyPoolType.means2$C_species_n = as.factor(as.character(plot_d13CbyPoolType.means2$C_species_n))
levels(plot_d13CbyPoolType.means2$C_species_n) <- c("DOC", "POC",  expression("CO"[2]))

## posthoc stats: posthoc stats: unpaired two-samples Wilcoxon test (by pool type and C species)
wilcox.test(DOC$d13C ~ DOC$Pool_type, paired = FALSE) # p = 3.647e-09
wilcox.test(POC$d13C ~ POC$Pool_type, paired = FALSE) # p = 0.00055
wilcox.test(CO2$d13C ~ CO2$Pool_type, paired = FALSE) # p = 0.0001459


Fig_d13CbyPoolType2 <- plot_d13CbyPoolType2 + 
  scale_fill_manual(values=c("#FFFFFF", "#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF")) +
  xlab('Pool type') + ylab(expression(paste(delta^13*'C (\u2030)'))) +
  scale_y_continuous(limits = c(-30,-8), breaks = c(-30,-25,-20,-15,-10)) +
  scale_x_discrete(labels=c("Natural", "Restoration")) +
  geom_hline(aes(yintercept=d13C), plot_d13CbyPoolType.means2, linetype="dashed", color = "black") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none")
Fig_d13CbyPoolType2

d13CbyPoolType.MaxValue.df2 = data.frame(d13C = c(-22,-25,-8.5,-27,-23,-10),
                                         Pool_type = as.character(c("Natural","Restoration","Natural","Restoration","Natural","Restoration")),
                                         lab = c('cd','cf','c','ce','a','b'),
                                         C_species_n = as.character(c("DOC","POC",(expression("CO"[2])))))
d13CbyPoolType.MaxValue.df2$C_species_n = as.factor(as.character(d13CbyPoolType.MaxValue.df2$C_species))
d13CbyPoolType.MaxValue.df2
# geom_text (below) was changing the order of facets in the plot
# so "fct_relevel" was used to reorder the names in the dataframe used for annotating the plot
# but this changed the annotation ordering, so the lettering and positioning values had to be ordered manually below
fct_relevel(d13CbyPoolType.MaxValue.df2$C_species_n)
fct_relevel(d13CbyPoolType.MaxValue.df2$C_species_n, "DOC", "POC")
Fig_d13CbyPoolType2 + geom_text(data = d13CbyPoolType.MaxValue.df2, label = c('cd','ce','a','cf','c','b'), vjust=0 , size = 2.5)

ggsave("Fig_4.png", width = 14, height = 10, units = c("cm"), dpi = 600)







##### Figure 5 #####
# F14C v d13C by C-species

## Data
All_data = read.csv("All_corr_data.csv")

# selections
DOC = All_data[which(All_data$C_species=="DOC"),]
POC = All_data[which(All_data$C_species=="POC"),]
CO2 = All_data[which(All_data$C_species=="CO2"),]


# DOC #
nat.DOC = DOC[which(DOC$Pool_type=="Natural"),]
res.DOC = DOC[which(DOC$Pool_type=="Restoration"),]

DOC.lm = lm(DOC$F14C ~ DOC$d13C)
summary(DOC.lm)$r.squared # R^2 = 0.203898
summary(DOC.lm)$coefficients[,4] # p-value = 0.007347662***
summary(DOC.lm)$coefficients[] # y = -0.01071184x + 0.76129323 

nat.DOC.lm = lm(nat.DOC$F14C ~ nat.DOC$d13C)
summary(nat.DOC.lm)$r.squared # R^2 = 0.07468534
summary(nat.DOC.lm)$coefficients[,4] # p-value = 0.21847419
summary(nat.DOC.lm)$coefficients[] # y = -0.008466108x + 0.819663097

res.DOC.lm = lm(res.DOC$F14C ~ res.DOC$d13C)
summary(res.DOC.lm)$r.squared # R^2 = 0.003179953
summary(res.DOC.lm)$coefficients[,4] # p-value = 8.618118e-01
summary(res.DOC.lm)$coefficients[] # y = 0.004658142x + 1.207576355 

Fig_DOC_d13CvF14C = ggplot(data = DOC,aes(x = d13C, y = F14C, group = Pool_type, color = Pool_type)) + 
  geom_point() + 
  xlab(expression(paste(delta^13*'C-DOC (\u2030)'))) + ylab(expression('F'^14*'C-DOC (fraction modern)')) +
  scale_x_continuous(limits = c(-30,-22.5), breaks = c(-30,-27.5,-25,-22.5)) +
  scale_y_continuous(limits = c(0.93,1.11), breaks = c(0.9,0.95,1,1.05,1.1)) +
  scale_color_manual(values = c("Natural" = "black", "Restoration" = "grey65")) +
  annotate("text", x=-24, y=0.94, size=2.5, label=expression('R'^2*' = 0.08, p = 0.219 '), color = "black") +
  annotate("text", x=-24, y=0.93, size=2.5, label=expression('R'^2*' = 0.00, p = 0.862 '), color = "grey65") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title=element_blank()) + theme(legend.position = c(0.21, 0.12)) +
  theme(legend.background = element_rect(fill='transparent'))
Fig_DOC_d13CvF14C


# POC #
nat.POC = POC[which(POC$Pool_type=="Natural"),]
res.POC = POC[which(POC$Pool_type=="Restoration"),]

POC.lm = lm(POC$F14C ~ POC$d13C)
summary(POC.lm)$r.squared # R^2 = 0.01152022
summary(POC.lm)$coefficients[,4] # p-value = 0.64330555
summary(POC.lm)$coefficients[] # y = 0.001700538x + 1.047021799

nat.POC.lm = lm(nat.POC$F14C ~ nat.POC$d13C)
summary(nat.POC.lm)$r.squared # R^2 = 0.1974546
summary(nat.POC.lm)$coefficients[,4] # p-value = 0.17091164
summary(nat.POC.lm)$coefficients[] # y = 0.007590842x + 1.197036131

res.POC.lm = lm(res.POC$F14C ~ res.POC$d13C)
summary(res.POC.lm)$r.squared # R^2 = 0.009565583
summary(res.POC.lm)$coefficients[,4] # p-value = 0.7880890
summary(res.POC.lm)$coefficients[] # y = -0.003096849x + 0.914765301

Fig_POC_d13CvF14C = ggplot(data = POC,aes(x = d13C, y = F14C, group = Pool_type, color = Pool_type)) + 
  geom_point() +
  xlab(expression(paste(delta^13*'C-POC (\u2030)'))) + ylab(expression('F'^14*'C-POC (fraction modern)')) +
  scale_x_continuous(limits = c(-30,-22.5), breaks = c(-30,-27.5,-25,-22.5)) +
  scale_y_continuous(limits = c(0.93,1.11), breaks = c(0.9,0.95,1,1.05,1.1)) +
  scale_color_manual(values = c("Natural" = "black", "Restoration" = "grey65")) +
  annotate("text", x=-24, y=0.94, size=2.5, label=expression('R'^2*' = 0.20, p = 0.171 '), color = "black") +
  annotate("text", x=-24, y=0.93, size=2.5, label=expression('R'^2*' = 0.01, p = 0.788 '), color = "grey65") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position = "none")
Fig_POC_d13CvF14C


# CO2 #
nat.CO2 = CO2[which(CO2$Pool_type=="Natural"),]
res.CO2 = CO2[which(CO2$Pool_type=="Restoration"),]

CO2.lm = lm(CO2$F14C ~ CO2$d13C)
summary(CO2.lm)$r.squared # R^2 = 0.4748007
summary(CO2.lm)$coefficients[,4] # p-value = 3.573178e-05***
summary(CO2.lm)$coefficients[] # y = -0.00393726x + 0.97045997

nat.CO2.lm = lm(nat.CO2$F14C ~ nat.CO2$d13C) 
summary(nat.CO2.lm)$r.squared # R^2 = 0.7564363
summary(nat.CO2.lm)$coefficients[,4] # p-value = 2.748718e-06***
summary(nat.CO2.lm)$coefficients[] # y = -0.007304258x + 0.931297047

res.CO2.lm = lm(res.CO2$F14C ~ res.CO2$d13C)
summary(res.CO2.lm)$r.squared # R^2 = 0.4670343
summary(res.CO2.lm)$coefficients[,4] # p-value = 0.02043317***
summary(res.CO2.lm)$coefficients[] # y = -0.004815698x + 0.949419334

Fig_CO2_d13CvF14C = ggplot(data = CO2,aes(x = d13C, y = F14C, group = Pool_type, color = Pool_type)) + 
  geom_point() +
  xlab(expression(paste(delta^13*'C-CO'[2]*' (\u2030)'))) + ylab(expression('F'^14*'C-CO'[2]*' (fraction modern)')) +
  geom_smooth(method=lm, se = FALSE, linewidth = 0.5) +
  scale_x_continuous(limits = c(-22,-8), breaks = c(-20,-15,-10)) +
  scale_y_continuous(limits = c(0.93,1.11), breaks = c(0.9,0.95,1,1.05,1.1)) +
  scale_color_manual(values = c("Natural" = "black", "Restoration" = "grey65")) +
  annotate("text", x=-16, y=1.07, size=2.5, label=expression('y = -0.0073x + 0.9313'), color = "black") +
  annotate("text", x=-18, y=1.01, size=2.5, label=expression('y = -0.0048x + 0.9494'), color = "grey65") +
  annotate("text", x=-11, y=0.94, size=2.5, label=expression('R'^2*' = 0.76, p = 0.000*'), color = "black") +
  annotate("text", x=-11, y=0.93, size=2.5, label=expression('R'^2*' = 0.47, p = 0.020*'), color = "grey65") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position = "none")
Fig_CO2_d13CvF14C

# Combine panels #
ggdraw() +
  draw_plot(Fig_DOC_d13CvF14C, x = 0, y = 0.2, width = .33, height = 0.8) +
  draw_plot(Fig_POC_d13CvF14C, x = .33, y = 0.2, width = .33, height = 0.8) +
  draw_plot(Fig_CO2_d13CvF14C, x = .66, y = 0.2, width = .33, height = 0.8) +
  draw_plot_label(label = c("a", "b", "c"), size = 15,
                  x = c(0, 0.33, 0.66), y = c(1, 1, 1))

ggsave("Fig_5.png", width = 25, height = 10, units = c("cm"), dpi = 600)







##### Figure 6 #####
# C-form vs C-form by isotope


## CO2 v DOC ##

# d13C #
CO2vDOC_d13C = read.csv("CO2vDOC_d13C.csv")
nat.CO2vDOC_d13C = CO2vDOC_d13C[which(CO2vDOC_d13C$Pool_type=="Natural"),]
res.CO2vDOC_d13C = CO2vDOC_d13C[which(CO2vDOC_d13C$Pool_type=="Restoration"),]

CO2vDOC_d13C.lm = lm(CO2vDOC_d13C$CO2_d13C ~ CO2vDOC_d13C$DOC_d13C)
summary(CO2vDOC_d13C.lm)$r.squared # R^2 = 0.2985788
summary(CO2vDOC_d13C.lm)$coefficients[,4] # p-value = 0.0021639***
summary(CO2vDOC_d13C.lm)$coefficients[] # y = 1.157514x + 17.421113

nat.CO2vDOC_d13C.lm = lm(nat.CO2vDOC_d13C$CO2_d13C ~ nat.CO2vDOC_d13C$DOC_d13C)
summary(nat.CO2vDOC_d13C.lm)$r.squared # R^2 = 0.002364946
summary(nat.CO2vDOC_d13C.lm)$coefficients[,4] # p-value = 0.8480375
summary(nat.CO2vDOC_d13C.lm)$coefficients[] # y = 0.07040938x - 10.85266281

res.CO2vDOC_d13C.lm = lm(res.CO2vDOC_d13C$CO2_d13C ~ res.CO2vDOC_d13C$DOC_d13C)
summary(res.CO2vDOC_d13C.lm)$r.squared # R^2 = 0.5194643
summary(res.CO2vDOC_d13C.lm)$coefficients[,4] # p-value = 0.01233551***
summary(res.CO2vDOC_d13C.lm)$coefficients[] # y = 4.382868x + 108.692992

Fig_CO2vDOC_d13C = ggplot(data = CO2vDOC_d13C,aes(x = DOC_d13C, y = CO2_d13C, group = Pool_type, color = Pool_type)) + 
  geom_point() + #geom_smooth(method=lm, se = FALSE, linewidth = 0.5) +
  geom_smooth(method=lm, se = FALSE, linewidth = 0.5, col = "mediumblue", aes(group=FALSE)) +
  geom_smooth(data = subset(CO2vDOC_d13C, Pool_type=="Restoration"), method=lm, se = FALSE, linewidth = 0.5) +
  xlab(expression(paste(delta^13*'C-DOC (\u2030)'))) + ylab(expression(paste(delta^13*'C-CO'[2]*' (\u2030)'))) +
  scale_x_continuous(limits = c(-32,-18), breaks = c(-30,-25,-20)) +
  scale_y_continuous(limits = c(-22,-8), breaks = c(-20,-15,-10)) +
  scale_color_manual(values = c("Natural" = "black", "Restoration" = "grey65")) +
  annotate("text", x=-21, y=-11, size=2.5, label=expression('y = 1.2x + 17.4'), color = "mediumblue") +
  annotate("text", x=-25, y=-18, size=2.5, label=expression('y = 4.4x + 108.7'), color = "grey65") +
  annotate("text", x=-21, y=-19.5, size=2.5, label=expression('R'^2*' = 0.30, p = 0.002*'), color = "mediumblue") +
  annotate("text", x=-21, y=-20.5, size=2.5, label=expression('R'^2*' = 0.00, p = 0.848 '), color = "black") +
  annotate("text", x=-21, y=-21.5, size=2.5, label=expression('R'^2*' = 0.52, p = 0.012*'), color = "grey65") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position = "none")
Fig_CO2vDOC_d13C


# F14C #
CO2vDOC_F14C = read.csv("CO2vDOC_F14C.csv")
nat.CO2vDOC_F14C = CO2vDOC_F14C[which(CO2vDOC_F14C$Pool_type=="Natural"),]
res.CO2vDOC_F14C = CO2vDOC_F14C[which(CO2vDOC_F14C$Pool_type=="Restoration"),]

CO2vDOC_F14C.lm = lm(CO2vDOC_F14C$CO2_F14C ~ CO2vDOC_F14C$DOC_F14C)
summary(CO2vDOC_F14C.lm)$r.squared # R^2 = 0.1535156
summary(CO2vDOC_F14C.lm)$coefficients[,4] # p-value = 0.03224995***
summary(CO2vDOC_F14C.lm)$coefficients[] # y = 0.2678998x + 0.7414279

nat.CO2vDOC_F14C.lm = lm(nat.CO2vDOC_F14C$CO2_F14C ~ nat.CO2vDOC_F14C$DOC_F14C)
summary(nat.CO2vDOC_F14C.lm)$r.squared # R^2 = 0.1438807
summary(nat.CO2vDOC_F14C.lm)$coefficients[,4] # p-value = 0.1092264 
summary(nat.CO2vDOC_F14C.lm)$coefficients[] # y = 0.2454520x + 0.7633348

res.CO2vDOC_F14C.lm = lm(res.CO2vDOC_F14C$CO2_F14C ~ res.CO2vDOC_F14C$DOC_F14C)
summary(res.CO2vDOC_F14C.lm)$r.squared # R^2 = 0.01834884
summary(res.CO2vDOC_F14C.lm)$coefficients[,4] # p-value = 0.691284541 
summary(res.CO2vDOC_F14C.lm)$coefficients[] # y = 0.1751647x + 0.8435962

Fig_CO2vDOC_F14C = ggplot(data = CO2vDOC_F14C,aes(x = DOC_F14C, y = CO2_F14C, group = Pool_type, color = Pool_type)) + 
  geom_point() + 
  geom_smooth(method=lm, se = FALSE, linewidth = 0.5, col = "mediumblue", aes(group=FALSE)) + 
  xlab(expression('F'^14*'C-DOC')) + ylab(expression('F'^14*'C-CO'[2])) +
  scale_x_continuous(limits = c(0.9, 1.2), breaks = c(1.2, 1.1, 1.0, 0.9)) +
  scale_y_continuous(limits = c(0.9, 1.2), breaks = c(1.2, 1.1, 1.0, 0.9)) +
  scale_color_manual(values = c("Natural" = "black", "Restoration" = "grey65")) +
  annotate("text", x=0.97, y=0.96, size=2.5, label=expression('y = 0.27x + 0.74'), color = "mediumblue") +
  annotate("text", x=1.14, y=1.19, size=2.5, label=expression('R'^2*' = 0.15, p = 0.032*'), color = "mediumblue") +
  annotate("text", x=1.14, y=1.17, size=2.5, label=expression('R'^2*' = 0.14, p = 0.109 '), color = "black") +
  annotate("text", x=1.14, y=1.15, size=2.5, label=expression('R'^2*' = 0.02, p = 0.691 '), color = "grey65") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title=element_blank()) + theme(legend.position = c(0.23, 0.89)) +
  theme(legend.background = element_rect(fill='transparent'))
Fig_CO2vDOC_F14C


## CO2 v POC ##

# d13C #
CO2vPOC_d13C = read.csv("CO2vPOC_d13C.csv")
nat.CO2vPOC_d13C = CO2vPOC_d13C[which(CO2vPOC_d13C$Pool_type=="Natural"),]
res.CO2vPOC_d13C = CO2vPOC_d13C[which(CO2vPOC_d13C$Pool_type=="Restoration"),]

CO2vPOC_d13C.lm = lm(CO2vPOC_d13C$CO2_d13C ~ CO2vPOC_d13C$POC_d13C)
summary(CO2vPOC_d13C.lm)$r.squared # R^2 = 0.430153
summary(CO2vPOC_d13C.lm)$coefficients[,4] # p-value = 0.00229562***
summary(CO2vPOC_d13C.lm)$coefficients[] # y = 1.537712x + 27.546912

nat.CO2vPOC_d13C.lm = lm(nat.CO2vPOC_d13C$CO2_d13C ~ nat.CO2vPOC_d13C$POC_d13C)
summary(nat.CO2vPOC_d13C.lm)$r.squared # R^2 = 0.006259433
summary(nat.CO2vPOC_d13C.lm)$coefficients[,4] # p-value = 0.8396633847 
summary(nat.CO2vPOC_d13C.lm)$coefficients[] # y = 0.09580662x - 9.06170756

res.CO2vPOC_d13C.lm = lm(res.CO2vPOC_d13C$CO2_d13C ~ res.CO2vPOC_d13C$POC_d13C)
summary(res.CO2vPOC_d13C.lm)$r.squared # R^2 = 0.1683534
summary(res.CO2vPOC_d13C.lm)$coefficients[,4] # p-value = 0.2389032
summary(res.CO2vPOC_d13C.lm)$coefficients[] # y = 1.161664x + 15.960551

Fig_CO2vPOC_d13C = ggplot(data = CO2vPOC_d13C,aes(x = POC_d13C, y = CO2_d13C, group = Pool_type, color = Pool_type)) + 
  geom_point() + 
  xlab(expression(paste(delta^13*'C-POC (\u2030)'))) + ylab(expression(paste(delta^13*'C-CO'[2]*' (\u2030)'))) +
  geom_smooth(method=lm, se = FALSE, linewidth = 0.5, col = "mediumblue", aes(group=FALSE)) + 
  scale_x_continuous(limits = c(-32,-18), breaks = c(-30,-25,-20)) +
  scale_y_continuous(limits = c(-22,-8), breaks = c(-20,-15,-10)) +
  scale_color_manual(values = c("Natural" = "black", "Restoration" = "grey65")) +
  annotate("text", x=-21.5, y=-11, size=2.5, label=expression('y = 1.5x + 27.6'), color = "mediumblue") +
  annotate("text", x=-21, y=-19.5, size=2.5, label=expression('R'^2*' = 0.43, p = 0.002*'), color = "mediumblue") +
  annotate("text", x=-21, y=-20.5, size=2.5, label=expression('R'^2*' = 0.01, p = 0.840 '), color = "black") +
  annotate("text", x=-21, y=-21.5, size=2.5, label=expression('R'^2*' = 0.17, p = 0.239 '), color = "grey65") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position = "none")
Fig_CO2vPOC_d13C

# F14C #
CO2vPOC_F14C = read.csv("CO2vPOC_F14C.csv")
nat.CO2vPOC_F14C = CO2vPOC_F14C[which(CO2vPOC_F14C$Pool_type=="Natural"),]
res.CO2vPOC_F14C = CO2vPOC_F14C[which(CO2vPOC_F14C$Pool_type=="Restoration"),]

CO2vPOC_F14C.lm = lm(CO2vPOC_F14C$CO2_F14C ~ CO2vPOC_F14C$POC_F14C)
summary(CO2vPOC_F14C.lm)$r.squared # R^2 = 0.01553693
summary(CO2vPOC_F14C.lm)$coefficients[,4] # p-value = 0.6005597050 
summary(CO2vPOC_F14C.lm)$coefficients[] # y = -0.1120056x + 1.1318664

nat.CO2vPOC_F14C.lm = lm(nat.CO2vPOC_F14C$CO2_F14C ~ nat.CO2vPOC_F14C$POC_F14C)
summary(nat.CO2vPOC_F14C.lm)$r.squared # R^2 = 0.4380042
summary(nat.CO2vPOC_F14C.lm)$coefficients[,4] # p-value = 0.0521723213 
summary(nat.CO2vPOC_F14C.lm)$coefficients[] # y = -0.8452345x + 1.8541131

res.CO2vPOC_F14C.lm = lm(res.CO2vPOC_F14C$CO2_F14C ~ res.CO2vPOC_F14C$POC_F14C)
summary(res.CO2vPOC_F14C.lm)$r.squared # R^2 = 0.0882664
summary(res.CO2vPOC_F14C.lm)$coefficients[,4] # p-value = 0.3749587 
summary(res.CO2vPOC_F14C.lm)$coefficients[] # y = 0.1692520x + 0.8629164

Fig_CO2vPOC_F14C = ggplot(data = CO2vPOC_F14C,aes(x = POC_F14C, y = CO2_F14C, group = Pool_type, color = Pool_type)) + 
  geom_point() + 
  xlab(expression('F'^14*'C-POC')) + ylab(expression('F'^14*'C-CO'[2])) +
  scale_x_continuous(limits = c(0.9, 1.2), breaks = c(1.2, 1.1, 1.0, 0.9)) +
  scale_y_continuous(limits = c(0.9, 1.2), breaks = c(1.2, 1.1, 1.0, 0.9)) +
  scale_color_manual(values = c("Natural" = "black", "Restoration" = "grey65")) +
  annotate("text", x=1.14, y=1.19, size=2.5, label=expression('R'^2*' = 0.02, p = 0.601 '), color = "mediumblue") +
  annotate("text", x=1.14, y=1.17, size=2.5, label=expression('R'^2*' = 0.44, p = 0.052 '), color = "black") +
  annotate("text", x=1.14, y=1.15, size=2.5, label=expression('R'^2*' = 0.09, p = 0.375 '), color = "grey65") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position = "none")
Fig_CO2vPOC_F14C


## DOC vs POC ##

# d13C #
DOCvPOC_d13C = read.csv("DOCvPOC_d13C.csv")
nat.DOCvPOC_d13C = DOCvPOC_d13C[which(DOCvPOC_d13C$Pool_type=="Natural"),]
res.DOCvPOC_d13C = DOCvPOC_d13C[which(DOCvPOC_d13C$Pool_type=="Restoration"),]

DOCvPOC_d13C.lm = lm(DOCvPOC_d13C$DOC_d13C ~ DOCvPOC_d13C$POC_d13C)
summary(DOCvPOC_d13C.lm)$r.squared # R^2 = 0.4861462
summary(DOCvPOC_d13C.lm)$coefficients[,4] # p-value = 0.0003104893***
summary(DOCvPOC_d13C.lm)$coefficients[] # y = 0.5611323x - 12.5347830

nat.DOCvPOC_d13C.lm = lm(nat.DOCvPOC_d13C$DOC_d13C ~ nat.DOCvPOC_d13C$POC_d13C)
summary(nat.DOCvPOC_d13C.lm)$r.squared # R^2 = 0.1652272
summary(nat.DOCvPOC_d13C.lm)$coefficients[,4] # p-value = 0.18977391 
summary(nat.DOCvPOC_d13C.lm)$coefficients[] # y = 0.3475535x - 17.7934844

res.DOCvPOC_d13C.lm = lm(res.DOCvPOC_d13C$DOC_d13C ~ res.DOCvPOC_d13C$POC_d13C)
summary(res.DOCvPOC_d13C.lm)$r.squared # R^2 = 0.03635343
summary(res.DOCvPOC_d13C.lm)$coefficients[,4] # p-value = 0.597753124
summary(res.DOCvPOC_d13C.lm)$coefficients[] # y = 0.09295355x - 26.03218233

Fig_DOCvPOC_d13C = ggplot(data = DOCvPOC_d13C,aes(x = POC_d13C, y = DOC_d13C, group = Pool_type, color = Pool_type)) + 
  geom_point() +
  geom_smooth(method=lm, se = FALSE, linewidth = 0.5, col = "mediumblue", aes(group=FALSE)) +
  xlab(expression(paste(delta^13*'C-POC (\u2030)'))) + ylab(expression(paste(delta^13*'C-DOC (\u2030)'))) +
  scale_x_continuous(limits = c(-32,-18), breaks = c(-30,-25,-20)) +
  scale_y_continuous(limits = c(-32,-18), breaks = c(-30,-25,-20)) +
  scale_color_manual(values = c("Natural" = "black", "Restoration" = "grey65")) +
  annotate("text", x=-20.75, y=-26, size=2.5, label=expression('y = 0.6x - 12.5'), color = "mediumblue") +
  annotate("text", x=-21, y=-29.5, size=2.5, label=expression('R'^2*' = 0.49, p < 0.001*'), color = "mediumblue") +
  annotate("text", x=-21, y=-30.5, size=2.5, label=expression('R'^2*' = 0.16, p = 0.190 '), color = "black") +
  annotate("text", x=-21, y=-31.5, size=2.5, label=expression('R'^2*' = 0.04, p = 0.598 '), color = "grey65") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position = "none")
Fig_DOCvPOC_d13C

# F14C #
DOCvPOC_F14C = read.csv("DOCvPOC_F14C.csv")
nat.DOCvPOC_F14C = DOCvPOC_F14C[which(DOCvPOC_F14C$Pool_type=="Natural"),]
res.DOCvPOC_F14C = DOCvPOC_F14C[which(DOCvPOC_F14C$Pool_type=="Restoration"),]

DOCvPOC_F14C.lm = lm(DOCvPOC_F14C$DOC_F14C ~ DOCvPOC_F14C$POC_F14C)
summary(DOCvPOC_F14C.lm)$r.squared # R^2 = 0.05556653
summary(DOCvPOC_F14C.lm)$coefficients[,4] # p-value = 0.2909289464 
summary(DOCvPOC_F14C.lm)$coefficients[] # y = 0.3021755x + 0.7480413

nat.DOCvPOC_F14C.lm = lm(nat.DOCvPOC_F14C$DOC_F14C ~ nat.DOCvPOC_F14C$POC_F14C)
summary(nat.DOCvPOC_F14C.lm)$r.squared # R^2 = 0.4296445
summary(nat.DOCvPOC_F14C.lm)$coefficients[,4] # p-value = 0.028560965***
summary(nat.DOCvPOC_F14C.lm)$coefficients[] # y = 1.10488337x - 0.07883028

res.DOCvPOC_F14C.lm = lm(res.DOCvPOC_F14C$DOC_F14C ~ res.DOCvPOC_F14C$POC_F14C)
summary(res.DOCvPOC_F14C.lm)$r.squared # R^2 = 0.03592208
summary(res.DOCvPOC_F14C.lm)$coefficients[,4] # p-value = 0.5767330 
summary(res.DOCvPOC_F14C.lm)$coefficients[] # y = 0.08349773x + 0.98943662

Fig_DOCvPOC_F14C = ggplot(data = DOCvPOC_F14C,aes(x = POC_F14C, y = DOC_F14C, group = Pool_type, color = Pool_type)) + 
  geom_point() + 
  geom_smooth(data = subset(DOCvPOC_F14C, Pool_type=="Natural"), method=lm, se = FALSE, linewidth = 0.5) +
  xlab(expression('F'^14*'C-POC')) + ylab(expression('F'^14*'C-DOC')) +
  scale_x_continuous(limits = c(0.9, 1.2), breaks = c(1.2, 1.1, 1.0, 0.9)) +
  scale_y_continuous(limits = c(0.9, 1.2), breaks = c(1.2, 1.1, 1.0, 0.9)) +
  scale_color_manual(values = c("Natural" = "black", "Restoration" = "grey65")) +
  annotate("text", x=1.1, y=1.01, size=2.5, label=expression('y = 1.1049x - 0.0788'), color = "black") +
  annotate("text", x=1.14, y=1.19, size=2.5, label=expression('R'^2*' = 0.06, p = 0.291 '), color = "mediumblue") +
  annotate("text", x=1.14, y=1.17, size=2.5, label=expression('R'^2*' = 0.43, p = 0.029*'), color = "black") +
  annotate("text", x=1.14, y=1.15, size=2.5, label=expression('R'^2*' = 0.04, p = 0.577 '), color = "grey65") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position = "none")
Fig_DOCvPOC_F14C


# Combine panels #
plot_grid(Fig_CO2vDOC_F14C, Fig_CO2vDOC_d13C,
          Fig_CO2vPOC_F14C, Fig_CO2vPOC_d13C,
          Fig_DOCvPOC_F14C, Fig_DOCvPOC_d13C,
          labels = c("a", "b", "c", "d", "e", "f"),
          ncol = 2, nrow=3)

ggsave("Fig_6.png", width = 15, height = 20, units = c("cm"), dpi = 600)





#### Figure 7 ####
# 14C - d13C space

Fig_F14Cd13Cspace = ggplot(F14C_data, aes(x = d13C, y = F14C, shape = factor(Pool_type), color = factor(C_species_n))) +
  geom_point(alpha=0.65) + 
  labs(shape = "Pool type", colour = "C forms") +
  scale_color_manual(values = c('black', 'gray65', 'dodgerblue', 'dodgerblue4', 'goldenrod1', 'goldenrod3', 'saddlebrown'),
                     labels = c(expression('DOC ('*italic('n')*' = 34)'), expression('POC ('*italic('n')*' = 21)'),
                                expression('CO'[2]*' ('*italic('n')*' = 30)'), 
                                expression('CO'[2]*' eb. ('*italic('n')*' = 1)'),
                                expression('CH'[4]*' eb. ('*italic('n')*' = 3)'), 
                                expression('CH'[4]*' ('*italic('n')*' = 0)'),
                                expression('Sed. ('*italic('n')*' = 3)'))) +
  scale_shape_manual(values = c(16, 17),
                     labels = c(expression('Natural ('*italic('n')*' = 57)'),
                                expression('Restoration ('*italic('n')*' = 34)'))) +
  coord_cartesian(xlim=c(-70, 0), ylim=c(0.7, 1.12)) +
  xlab(expression(paste(delta^13*'C (\u2030)'))) + ylab(expression('F'^14*'C (fraction modern)')) +
  annotate("text", x=-68, y=1.01, size=2.5, label=expression(italic('F'^14*'C modern')), color = "snow3") +
  geom_hline(yintercept=1, linetype="dashed", linewidth = 0.25, color = "snow3") +
  annotate("text", x=-3.7, y=1.0575, size=3, label='atmospheric', color = "black") +
  annotate("text", x=-3.5, y=1.0385, size=3, label=expression('CO'[2]*' 2013-15'), color = "black") +
  geom_rect(aes(xmin = -9.5, xmax = -8.5, ymin = 1.0234, ymax = 1.0301), linewidth = 0.15, fill = "white", alpha = 0, color = "black") +
  annotate("text", x=-40, y=0.983, size=3, label='modern plant material', color = "black") +
  annotate("text", x=-35.25, y=0.968, size=3, label='1950-2012', color = "black") +
  geom_rect(aes(xmin = -30, xmax = -23, ymin = 0.9740, ymax = 1.8495), linewidth = 0.15, fill = "white", alpha = 0, color = "black") +
  annotate("text", x=-49.75, y=1.115, size=3, label=expression('peak atmospheric '^14*'CO'[2]*' in 1964 = 1.8495'), color = "black") +
  geom_segment(aes(x = -31, y = 1.1, xend = -31, yend = 1.13), arrow = arrow(length=unit(.125, 'cm')), linewidth = 0.3, color = "black") +
  annotate("text", x=-37.5, y=0.85, size=3, label='0.5-1 m deep peat', color = "black") +
  annotate("text", x=-39.5, y=0.835, size=3, label='~1,000-3,000 years old', color = "black") +
  geom_rect(aes(xmin = -29, xmax = -26, ymin = 0.7035, ymax = 0.8788), linewidth = 0.15, fill = "white", alpha = 0, color = "black") + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.text.align = 0, legend.text=element_text(size = 9), legend.title=element_text(size = 10)) +
  theme(legend.spacing.y = unit(0.05, 'cm'))  +
  theme(legend.background = element_rect(fill='transparent'))
Fig_F14Cd13Cspace

ggsave("Fig_7.png", width = 18, height = 10, units = c("cm"), dpi = 600)







######## Supplementary figures ########


##### Table S3 #####
# Correlation analysis for Cross Lochs sites only, because other locations don't have supporting data (e.g. C concentrations)

# load data
CL_data = All_data[which(All_data$Location=="Cross Lochs"),]

# selections
DOC.CL = CL_data[which(CL_data$C_species=="DOC"),]
POC.CL = CL_data[which(CL_data$C_species=="POC"),]
CO2.CL = CL_data[which(CL_data$C_species=="CO2"),]

# correlation matrices
DOC.CLmatrix = matrix(data = c(DOC.CL$F14C, DOC.CL$d13C, DOC.CL$DOC_conc, DOC.CL$POC_conc,
                               DOC.CL$CO2_conc, DOC.CL$CH4_conc, DOC.CL$pH, DOC.CL$EC,
                               DOC.CL$DO, DOC.CL$Temp, DOC.CL$WTD),
                      nrow = 18, ncol = 11, byrow = FALSE)
colnames(DOC.CLmatrix) = c("F14C", "d13C", "DOC_conc", "POC_conc", "CO2_conc", "CH4_conc", "pH", "EC", "DO", "Temp", "WTD")

POC.CLmatrix = matrix(data = c(POC.CL$F14C, POC.CL$d13C, POC.CL$DOC_conc, POC.CL$POC_conc,
                               POC.CL$CO2_conc, POC.CL$CH4_conc, POC.CL$pH, POC.CL$EC,
                               POC.CL$DO, POC.CL$Temp, POC.CL$WTD),
                      nrow = 17, ncol = 11, byrow = FALSE)
colnames(POC.CLmatrix) = c("F14C", "d13C", "DOC_conc", "POC_conc", "CO2_conc", "CH4_conc", "pH", "EC", "DO", "Temp", "WTD")

CO2.CLmatrix = matrix(data = c(CO2.CL$F14C, CO2.CL$d13C, CO2.CL$DOC_conc, CO2.CL$POC_conc,
                               CO2.CL$CO2_conc, CO2.CL$CH4_conc, CO2.CL$pH, CO2.CL$EC,
                               CO2.CL$DO, CO2.CL$Temp, CO2.CL$WTD),
                      nrow = 17, ncol = 11, byrow = FALSE)
colnames(CO2.CLmatrix) = c("F14C", "d13C", "DOC_conc", "POC_conc", "CO2_conc", "CH4_conc", "pH", "EC", "DO", "Temp", "WTD")

# correlations
DOC.CLcor <- rcorr(DOC.CLmatrix, type="spearman")
write.csv(DOC.CLcor$P, file = "DOC.CLcor_P.csv")
write.csv(DOC.CLcor$r, file = "DOC.CLcor_r.csv")
DOC.CLcor$n

POC.CLcor <- rcorr(POC.CLmatrix, type="spearman")
write.csv(POC.CLcor$P, file = "POC.CLcor_P.csv")
write.csv(POC.CLcor$r, file = "POC.CLcor_r.csv")
POC.CLcor$n

CO2.CLcor <- rcorr(CO2.CLmatrix, type="spearman")
write.csv(CO2.CLcor$P, file = "CO2.CLcor_P.csv")
write.csv(CO2.CLcor$r, file = "CO2.CLcor_r.csv")
CO2.CLcor$n

# values manually extracted to Table S3 in the Supplementary Information file.






##### Figure S1 #####
# Physical parameters of the two different pool types (where available)

Pool_phys = read.csv('Pool_physical_parameters.csv')

# Depth (cm) 
Fig_pooldepth <- ggplot(Pool_phys, aes(x=Pool_type, y=Depth, fill= Pool_type)) + 
  geom_boxplot(outlier.shape=NA, color = "snow3") +
  geom_jitter(width = 0.1, height = 0, shape=16, color="black", size=2, alpha=0.7) +
  scale_fill_manual(values=c("#FFFFFF", "#FFFFFF")) +
  scale_x_discrete(labels = c(expression('Natural ('*italic(n)*' = 19)'), expression('Restoration ('*italic(n)*' = 7)'))) +
  scale_y_continuous(limits = c(5, 60), breaks = c(10, 20, 30, 40, 50, 60)) +
  xlab('Pool type') + ylab('Pool depth (cm)') +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position = "none") +
  theme (axis.title=element_text(size=8), axis.text = element_text(size = 7))
Fig_pooldepth

# Pool surface area (m^2)
Fig_poolarea <- ggplot(Pool_phys, aes(x=Pool_type, y=Area, fill= Pool_type)) + 
  geom_boxplot(outlier.shape=NA, color = "snow3") +
  geom_jitter(width = 0.1, height = 0, shape=16, color="black", size=2, alpha=0.7) +
  scale_fill_manual(values=c("#FFFFFF", "#FFFFFF")) +
  scale_x_discrete(labels = c(expression('Natural ('*italic(n)*' = 19)'), expression('Restoration ('*italic(n)*' = 9)'))) +
  scale_y_continuous(trans='log2', limits = c(1, 3000), breaks = c(1, 10, 100, 1000)) +
  xlab('Pool type') + ylab(expression('Pool surface area (m'^2*')')) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position = "none") +
  theme (axis.title=element_text(size=8), axis.text = element_text(size = 7))
Fig_poolarea

# r_time
Fig_restime <- ggplot(Pool_phys, aes(x=Pool_type, y=r_time, fill= Pool_type)) + 
  geom_boxplot(outlier.shape=NA, color = "snow3") +
  geom_jitter(width = 0.1, height = 0, shape=16, color="black", size=2, alpha=0.7) +
  scale_fill_manual(values=c("#FFFFFF", "#FFFFFF")) +
  scale_x_discrete(labels = c(expression('Natural ('*italic(n)*' = 6)'), expression('Restoration ('*italic(n)*' = 6)'))) +
  xlab('Pool type') + ylab('Mean residence time (days)') +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position = "none") +
  theme (axis.title=element_text(size=8), axis.text = element_text(size = 7))
Fig_restime

# d_area
Fig_d_area <- ggplot(Pool_phys, aes(x=Pool_type, y=d_area, fill= Pool_type)) + 
  geom_boxplot(outlier.shape=NA, color = "snow3") +
  geom_jitter(width = 0.1, height = 0, shape=16, color="black", size=2, alpha=0.7) +
  scale_fill_manual(values=c("#FFFFFF", "#FFFFFF")) +
  scale_x_discrete(labels = c(expression('Natural ('*italic(n)*' = 6)'), expression('Restoration ('*italic(n)*' = 6)'))) +
  scale_y_continuous(limits = c(0, 1300), breaks = c(0,250,500,750,1000,1250)) +
  xlab('Pool type') + ylab(expression('Upslope drainage area (m'^2*')')) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position = "none") +
  theme (axis.title=element_text(size=8), axis.text = element_text(size = 7))
Fig_d_area

# Combine panels #
plot_grid(Fig_pooldepth, Fig_poolarea, Fig_restime, Fig_d_area,
          labels = c("a", "b", "c", "d"),
          ncol = 2, nrow = 2)

ggsave("Fig_S1.jpeg", width = 12.5, height = 12.5, units = c("cm"), dpi = 600)






#### Figure S2 ####
# Importance of sampling campaign timing on isotope observations - "by date" 

## F14C ##

# All C forms by date

## reload data
F14C_data_no_sed <- F14C_data[-c(50,51,52),]

## label factor with names rather than letters but keep order
F14C_data_no_sed$C_species_n = as.factor(as.character(F14C_data_no_sed$C_species_n))
levels(F14C_data_no_sed$C_species_n) <- c("DOC", "POC",  expression("CO"[2]), expression("CO"[2]*"-eb."),
                                          expression("CH"[4]*"-eb."), expression("CH"[4]))

plot_F14CbyDatebyCspecies <- ggplot(F14C_data_no_sed, aes(x=C_species_n, y=F14C, fill = C_species_n)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width = 0.1, height = 0, shape=16, color="black", size=1.2, alpha=0.7) +
  facet_wrap(~Date_n, scale='free', labeller = label_parsed) +
  scale_fill_manual(values=c("#FFFFFF", "#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF")) +
  scale_y_continuous(limits = c(0.85,1.11), breaks = c(0.85,0.9,0.95,1,1.05,1.1)) +
  scale_x_discrete(labels=c("DOC", "POC",  expression("CO"[2]), expression("CO"[2]*"-eb."),
                            expression("CH"[4]*"-eb."), expression("CH"[4]))) +
  xlab('C form') + ylab(expression('F'^14*'C (fraction modern)')) +
  geom_hline(yintercept=1, linetype="dashed", color = "snow3") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(legend.position = "none")
plot_F14CbyDatebyCspecies


# DOC, POC and CO2 date comparison

## remove CO2-eb and CH4 values
F14C_data_no_eb = F14C_data[-c(27,28,29,39,50,51,52,94,95,96,97),]

## label factor with names rather than letters but keep order
F14C_data_no_eb$C_species_n = as.factor(as.character(F14C_data_no_eb$C_species_n))
levels(F14C_data_no_eb$C_species_n) <- c("DOC", "POC",  expression("CO"[2]))

plot_F14CbyCspecies_datecomparison = ggplot(F14C_data_no_eb, aes(x=Date_n, y=F14C, fill = Date_n)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width = 0.1, height = 0, shape=16, color="black", size=1.2, alpha=0.7) +
  facet_wrap(~C_species_n, scale='free', labeller = label_parsed) +
  scale_fill_manual(values=c("#FFFFFF", "#FFFFFF","#FFFFFF")) +
  scale_y_continuous(limits = c(0.85,1.11), breaks = c(0.85,0.9,0.95,1,1.05,1.1)) +
  xlab('Sampling campaign') + ylab(expression('F'^14*'C (fraction modern)')) +
  geom_hline(yintercept=1, linetype="dashed", color = "snow3") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(legend.position = "none")
plot_F14CbyCspecies_datecomparison


# Location by date

## reload data
F14C_data_no_sed <- F14C_data[-c(50,51,52),]

## label factor with names rather than letters but keep order
F14C_data_no_sed$Location_n = as.factor(as.character(F14C_data_no_sed$Location_n))
levels(F14C_data_no_sed$Location_n) <- c("Cross Lochs", "Loch Leir", "Munsary", "Silver Flowe", "Garron Plateau", "Slieveanorra")

plot_F14CbyDatebyLocation <- ggplot(F14C_data_no_sed, aes(x=Location_n, y=F14C, fill = Location_n)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width = 0.1, height = 0, shape=16, color="black", size=1.2, alpha=0.7) +
  facet_wrap(~Date_n, scale='free', labeller = label_parsed) +
  scale_y_continuous(limits = c(0.85,1.11), breaks = c(0.85,0.9,0.95,1,1.05,1.1)) +
  scale_fill_manual(values=c("#FFFFFF", "#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF")) +
  xlab('Location') + ylab(expression('F'^14*'C (fraction modern)')) +
  geom_hline(yintercept=1, linetype="dashed", color = "snow3") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(legend.position = "none")
plot_F14CbyDatebyLocation


# Combine panels #
plot_grid(plot_F14CbyDatebyCspecies,
          plot_F14CbyCspecies_datecomparison,
          plot_F14CbyDatebyLocation,
          labels = c("a", "b", "c"),
          ncol = 1, nrow=3)

ggsave("Fig_S2.jpeg", width = 14, height = 25, units = c("cm"), dpi = 600)






#### Figure S3 ####
# Importance of sampling campaign timing on isotope observations - "by date" 

## d13C ##

# DOC, POC and CO2 date comparison

## reload data
All_data = read.csv("All_corr_data.csv")

## label factor with names rather than letters but keep order
All_data$C_species_n = as.factor(as.character(All_data$C_species_n))
levels(All_data$C_species_n) <- c("DOC", "POC",  expression("CO"[2]))

plot_d13CCbyCspecies_datecomparison = ggplot(All_data, aes(x=Date_n, y=d13C, fill = Date_n)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width = 0.1, height = 0, shape=16, color="black", size=1.2, alpha=0.7) +
  facet_wrap(~C_species_n, scale='free', labeller = label_parsed) +
  scale_fill_manual(values=c("#FFFFFF", "#FFFFFF","#FFFFFF")) +
  scale_y_continuous(limits = c(-31,-9), breaks = c(-30,-25,-20,-15,-10)) +
  xlab('Sampling campaign') + ylab(expression(paste(delta^13*'C (\u2030)'))) +
  geom_hline(yintercept=1, linetype="dashed", color = "snow3") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(legend.position = "none")
plot_d13CCbyCspecies_datecomparison

ggsave("Fig_S3.jpeg", width = 14, height = 8, units = c("cm"), dpi = 600)





#### Figure S4 ####

# Atmospheric F14C-CO2 and age distribution solutions

FigF14CvsyBSD = ggplot(data = F14C_data_no_sed, aes(x = old_soln, y = F14C, group = C_species_n, color = Pool_type)) +
  geom_point() +
  geom_point(data = F14C_data_no_sed, aes(x = young_soln, y = F14C, group = C_species_n, color = Pool_type)) +
  geom_line(data = atmdata, aes(x = yBSD, y = F14C, group = NA), color = 'mediumblue') +
  scale_x_continuous(limits = c(1,10000), breaks = c(1,10,100,1000,10000), trans='log10') +
  scale_y_continuous(limits = c(0.35,1.85), breaks = c(0.4,0.6,0.8,1,1.2,1.4,1.6,1.8)) +
  xlab("Mean modelled age (yBSD)") + ylab(expression('F'^14*'C (fraction modern)')) + 
  scale_color_manual(values = c("Natural" = "black", "Restoration" = "grey65")) +
  annotate("text", x=3, y=1.25, size=2.5, label=expression('young solutions, '*italic(n)*' = 49'), color = "black") +
  annotate("text", x=3.2, y=1.175, size=2.5, label=expression('natural pools '*italic(n)*' = 27 (55%)'), color = "black") +
  annotate("text", x=200, y=1.15, size=2.5, label=expression('old solutions, '*italic(n)*' = 94'), color = "black") +
  annotate("text", x=150, y=1.6, size=2.5, label=expression('Atmospheric F'^14*'C-CO'[2]), color = "mediumblue") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title=element_blank()) + theme(legend.position = c(0.1, 0.12)) +
  theme(legend.background = element_rect(fill='transparent'))
FigF14CvsyBSD

ggsave("Fig_S4.jpeg", width = 16, height = 8, units = c("cm"), dpi = 600)






##### Figure S5 #####
# DOM correlation analysis plus figure panels

# load data
DOM = read.csv("DOMvDOC.csv")

# correlations
d13CvE4E6.lm = lm(DOM$d13C ~ DOM$E4E6)
summary(d13CvE4E6.lm)$r.squared # R2 = 0.002253518
summary(d13CvE4E6.lm)$coefficients[,4] # p = 0.8963904
summary(d13CvE4E6.lm)$coefficients[] # y = 0.01602381x - 28.05859446

F14CvE4E6.lm = lm(DOM$F14C ~ DOM$E4E6)
summary(F14CvE4E6.lm)$r.squared # R2 = 0.01535382
summary(F14CvE4E6.lm)$coefficients[,4] # p = 0.7330692 
summary(F14CvE4E6.lm)$coefficients[] # y = 0.001286215x + 1.041363234

d13CvSUVA.lm = lm(DOM$d13C ~ DOM$SUVA)
summary(d13CvSUVA.lm)$r.squared # R2 = 0.8803434
summary(d13CvSUVA.lm)$coefficients[,4] # p = 6.362815e-06***
summary(d13CvSUVA.lm)$coefficients[] # y = -1.28446x - 23.82263

F14CvSUVA.lm = lm(DOM$F14C ~ DOM$SUVA)
summary(F14CvSUVA.lm)$r.squared # R2 = 0.7395205
summary(F14CvSUVA.lm)$coefficients[,4] # p = 0.0003337715***
summary(F14CvSUVA.lm)$coefficients[] # y = 0.03369893x + 0.94222292

Fig_DOCd13CvSUVA = ggplot(data = DOM, aes(x = SUVA, y = d13C, group = Pool_type, color = Pool_type)) + 
  geom_point() + 
  geom_smooth(method=lm, se = FALSE, linewidth = 0.5, col = "mediumblue", aes(group=FALSE)) +
  xlab('SUVA') + ylab(expression(paste(delta^13*'C-DOC (\u2030)'))) +
  scale_y_continuous(limits = c(-30,-25), breaks = c(-30, -29, -28, -27, -26, -25, -25)) +
  scale_color_manual(values = c("Natural" = "black", "Restoration" = "grey65")) +
  annotate("text", x=2.5, y=-28.3, size=2.5, label=expression('y = -1.3x - 23.8'), color = "mediumblue") +
  annotate("text", x=2.3, y=-29.9, size=2.5, label=expression('R'^2*' = 0.88, p < 0.001*'), color = "mediumblue") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position = "none") 
Fig_DOCd13CvSUVA

Fig_DOCF14CvSUVA = ggplot(data = DOM, aes(x = SUVA, y = F14C, group = Pool_type, color = Pool_type)) + 
  geom_point() + 
  geom_smooth(method=lm, se = FALSE, linewidth = 0.5, col = "mediumblue", aes(group=FALSE)) +
  xlab('SUVA') + ylab(expression('F'^14*'C-DOC')) +
  scale_y_continuous(limits = c(0.9, 1.1), breaks = c(0.9, 0.95, 1, 1.05, 1.1)) +
  scale_color_manual(values = c("Natural" = "black", "Restoration" = "grey65")) +
  annotate("text", x=3, y=0.975, size=2.5, label=expression('y = 0.0337x + 0.9422'), color = "mediumblue") +
  annotate("text", x=3.55, y=0.905, size=2.5, label=expression('R'^2*' = 0.74, p < 0.001*'), color = "mediumblue") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title=element_blank()) + theme(legend.position = c(0.275, 0.89)) +
  theme(legend.background = element_rect(fill='transparent'))
Fig_DOCF14CvSUVA


# Combine panels #
plot_grid(Fig_DOCF14CvSUVA, Fig_DOCd13CvSUVA,
          labels = c("a", "b"),
          ncol = 2, nrow = 1)

ggsave("Fig_S5.jpeg", width = 12.5, height = 6.25, units = c("cm"), dpi = 600)





#### Figure S6 ####
# Keeling and Miller-Tans plots

# Isotope source estimates: d13C

# Keeling

CO2 = All_data[which(All_data$C_species=="CO2"),]
Keeling_x = as.numeric(1/CO2[,18])
CO2$Keeling_x = c(Keeling_x)

nat.CO2 = CO2[which(CO2$Pool_type=="Natural"),]
res.CO2 = CO2[which(CO2$Pool_type=="Restoration"),]

Keeling_CO2.lm = lm(CO2$d13C ~ CO2$Keeling_x)
summary(Keeling_CO2.lm)$r.squared # R^2 = 0.6957015
summary(Keeling_CO2.lm)$coefficients[,4] # p-value = 0.00005912619***
summary(Keeling_CO2.lm)$coefficients[] # y = 2.64526x - 18.54693

Keeling_CO2_nat.lm = lm(nat.CO2$d13C ~ nat.CO2$Keeling_x)
summary(Keeling_CO2_nat.lm)$r.squared # R^2 = 0.1674438
summary(Keeling_CO2_nat.lm)$coefficients[,4] # p-value = 0.361990532
summary(Keeling_CO2_nat.lm)$coefficients[] # y = 1.182897x - 14.629201

Keeling_CO2_res.lm = lm(res.CO2$d13C ~ res.CO2$Keeling_x)
summary(Keeling_CO2_res.lm)$r.squared # R^2 = 0.2964396
summary(Keeling_CO2_res.lm)$coefficients[,4] # p-value = 0.1296103
summary(Keeling_CO2_res.lm)$coefficients[] # y = 2.607486x - 18.697676

Fig_Keeling_CO2 = ggplot(CO2, aes(x = Keeling_x, y = d13C, group = Pool_type, color = Pool_type)) +
  geom_point() +
  xlab(expression('1/CO'[2]*' (mg C/L)')) + ylab(expression(paste(delta^13*'C-CO'[2]*' (\u2030)'))) +
  geom_smooth(method=lm, se = FALSE, linewidth = 0.5, col = "mediumblue", aes(group=FALSE)) +
  scale_color_manual(values = c("Natural" = "black", "Restoration" = "grey65")) +
  annotate("text", x=2.5, y=-16, size=2.5, label=expression('intercept = -18.6 \u2030'), color = "mediumblue") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title=element_blank()) + theme(legend.position = "none") +
  theme(legend.background = element_rect(fill='transparent')) +
  theme (axis.title=element_text(size=9), axis.text = element_text(size = 8))
Fig_Keeling_CO2

# Miller-Tans

MillerTans_y = as.numeric(CO2[,10] * CO2[,18])
CO2$MillerTans_y = c(MillerTans_y)

MillerTans_CO2.lm = lm(CO2$MillerTans_y ~ CO2$CO2_conc)
summary(MillerTans_CO2.lm)$r.squared # R^2 = 0.9942459
summary(MillerTans_CO2.lm)$coefficients[,4] # p-value = 4.386138e-17***
summary(MillerTans_CO2.lm)$coefficients[] # y = -18.569474x -18.569474 + 2.519002

MillerTans_CO2_nat.lm = lm(nat.CO2$MillerTans_y ~ nat.CO2$CO2_conc)
summary(MillerTans_CO2_nat.lm)$r.squared # R^2 = 0.8501815
summary(MillerTans_CO2_nat.lm)$coefficients[,4] # p-value = 0.003123012***
summary(MillerTans_CO2_nat.lm)$coefficients[] # y = -15.816170x + 1.665006

MillerTans_CO2_res.lm = lm(res.CO2$MillerTans_y ~ res.CO2$CO2_conc)
summary(MillerTans_CO2_res.lm)$r.squared # R^2 = 0.9900906
summary(MillerTans_CO2_res.lm)$coefficients[,4] # p-value = 2.829911e-08***
summary(MillerTans_CO2_res.lm)$coefficients[] # y = -18.405577x + 1.828463

Fig_MillerTans_CO2 = ggplot(CO2, aes(x = CO2_conc, y = MillerTans_y, group = Pool_type, color = Pool_type)) +
  geom_point() +
  xlab(expression('CO'[2]*' (mg C/L)')) + ylab(expression(paste(delta^13*'C-CO'[2]*' (\u2030) x CO'[2]*' (mg C/L)'))) +
  #geom_smooth(method=lm, se = FALSE, linewidth = 0.5, col = "mediumblue", aes(group=FALSE)) +
  geom_smooth(method=lm, se = FALSE, linewidth = 0.5) +
  scale_color_manual(values = c("Natural" = "black", "Restoration" = "grey65")) +
  annotate("text", x=1.5, y=-105, size=2.5, label=expression('slope = -18.6 \u2030'), color = "mediumblue") +
  annotate("text", x=2.2, y=-5, size=2.5, label=expression('slope = -15.8 \u2030'), color = "black") +
  annotate("text", x=4, y=-40, size=2.5, label=expression('slope = -18.4 \u2030'), color = "grey65") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title=element_blank()) + theme(legend.position = "none") +
  theme(legend.background = element_rect(fill='transparent')) +
  theme (axis.title=element_text(size=9), axis.text = element_text(size = 8))
Fig_MillerTans_CO2


## Isotope source estimates: F14C

# Keeling

Keeling_F14CO2.lm = lm(CO2$F14C ~ CO2$Keeling_x)
summary(Keeling_F14CO2.lm)$r.squared # R^2 = 0.7517929
summary(Keeling_F14CO2.lm)$coefficients[,4] # p-value = 0.000006631583***
summary(Keeling_F14CO2.lm)$coefficients[] # y = -0.01752729x + 1.04858469

Keeling_F14CO2_nat.lm = lm(nat.CO2$F14C ~ nat.CO2$Keeling_x)
summary(Keeling_F14CO2_nat.lm)$r.squared # R^2 = 0.6797816
summary(Keeling_F14CO2_nat.lm)$coefficients[,4] # p-value = 0.01179949***
summary(Keeling_F14CO2_nat.lm)$coefficients[] # y = -0.02719809x + 1.07674165

Keeling_F14CO2_res.lm = lm(res.CO2$F14C ~ res.CO2$Keeling_x)
summary(Keeling_F14CO2_res.lm)$r.squared # R^2 = 0.3769121
summary(Keeling_F14CO2_res.lm)$coefficients[,4] # p-value = 0.07862482
summary(Keeling_F14CO2_res.lm)$coefficients[] # y = -0.010259x + 1.042465

Fig_Keeling_F14CO2 = ggplot(CO2, aes(x = Keeling_x, y = F14C, group = Pool_type, color = Pool_type)) +
  geom_point() +
  xlab(expression('1/CO'[2]*' (mg C/L)')) + ylab(expression('F'^14*'C-CO'[2])) +
  geom_smooth(method=lm, se = FALSE, linewidth = 0.5, col = "mediumblue", aes(group=FALSE)) +
  geom_smooth(data = subset(CO2, Pool_type=="Natural"), method=lm, se = FALSE, linewidth = 0.5) +
  scale_color_manual(values = c("Natural" = "black", "Restoration" = "grey65")) +
  annotate("text", x=2, y=1.045, size=2.5, label=expression('intercept F'^14*'C = 1.0486'), color = "mediumblue") +
  annotate("text", x=3.1, y=1.032, size=2.5, label=expression('intercept F'^14*'C = 1.0767'), color = "black") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title=element_blank()) + theme(legend.position = c(0.29, 0.15)) +
  theme(legend.background = element_rect(fill='transparent')) +
  theme (axis.title=element_text(size=9), axis.text = element_text(size = 8))
Fig_Keeling_F14CO2

# Miller-Tans

MillerTansF14C_y = as.numeric(CO2[,8] * CO2[,18])
CO2$MillerTansF14C_y = c(MillerTansF14C_y)

MillerTansF14C_CO2.lm = lm(CO2$MillerTansF14C_y ~ CO2$CO2_conc)
summary(MillerTansF14C_CO2.lm)$r.squared # R^2 = 0.9999308
summary(MillerTansF14C_CO2.lm)$coefficients[,4] # p-value = 1.285960e-32***
summary(MillerTansF14C_CO2.lm)$coefficients[] # y = 1.04274423x -18.569474 - 0.01233804

MillerTans_F14CO2_nat.lm = lm(nat.CO2$MillerTansF14C_y ~ nat.CO2$CO2_conc)
summary(MillerTans_F14CO2_nat.lm)$r.squared # R^2 = 0.9982607
summary(MillerTans_F14CO2_nat.lm)$coefficients[,4] # p-value = 1.64525e-09***
summary(MillerTans_F14CO2_nat.lm)$coefficients[] # y = 1.05924974x - 0.02060613

MillerTans_F14CO2_res.lm = lm(res.CO2$MillerTansF14C_y ~ res.CO2$CO2_conc)
summary(MillerTans_F14CO2_res.lm)$r.squared # R^2 = 0.9998801
summary(MillerTans_F14CO2_res.lm)$coefficients[,4] # p-value = 5.487818e-15***
summary(MillerTans_F14CO2_res.lm)$coefficients[] # y = 1.041500029x - 0.007249797

Fig_MillerTansF14C_CO2 = ggplot(CO2, aes(x = CO2_conc, y = MillerTansF14C_y, group = Pool_type, color = Pool_type)) +
  geom_point() +
  xlab(expression('CO'[2]*' (mg C/L)')) + ylab(expression(paste('F'^14*'C-CO'[2]*' x CO'[2]*' (mg C/L)'))) +
  #geom_smooth(method=lm, se = FALSE, linewidth = 0.5, col = "mediumblue", aes(group=FALSE)) +
  geom_smooth(method=lm, se = FALSE, linewidth = 0.5) +
  scale_color_manual(values = c("Natural" = "black", "Restoration" = "grey65")) +
  annotate("text", x=1.75, y=6, size=2.5, label=expression('slope F'^14*'C = 1.0427'), color = "mediumblue") +
  annotate("text", x=2.2, y=0.3, size=2.5, label=expression('slope = F'^14*'C = 1.0593'), color = "black") +
  annotate("text", x=4, y=2.2, size=2.5, label=expression('slope = F'^14*'C = 1.0415'), color = "grey65") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title=element_blank()) + theme(legend.position = "none") +
  theme(legend.background = element_rect(fill='transparent')) +
  theme (axis.title=element_text(size=9), axis.text = element_text(size = 8))
Fig_MillerTansF14C_CO2


# Combine panels #
plot_grid(Fig_Keeling_F14CO2, Fig_Keeling_CO2, 
          Fig_MillerTansF14C_CO2, Fig_MillerTans_CO2,
          labels = c("a", "b", "c", "d"),
          ncol = 2, nrow = 2)

ggsave("Fig_S6.jpeg", width = 12.5, height = 12.5, units = c("cm"), dpi = 600)





#### Figure S7 ####
# DOC v CO2 potential outlier removal 

# d13C

# identify potential outliers
# https://stats.stackexchange.com/questions/164099/removing-outliers-based-on-cooks-distance-in-r-language 

DOCvCO2_d13C_cooksd <- cooks.distance(CO2vDOC_d13C.lm)
sample_size1 = nrow(CO2vDOC_d13C)
plot(DOCvCO2_d13C_cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")
abline(h = 4/sample_size1, col="red")  # add cutoff line
text(x=1:length(DOCvCO2_d13C_cooksd)+1, y=DOCvCO2_d13C_cooksd,
     labels=ifelse(DOCvCO2_d13C_cooksd>4/sample_size1, names(DOCvCO2_d13C_cooksd),""), col="red")  # add labels
# identifies:
# 09. P01 15/09/2015
# 26. SF06 10/09/2015

nat.CO2vDOC_d13C_cooksd <- cooks.distance(nat.CO2vDOC_d13C.lm)
sample_size2 = nrow(nat.CO2vDOC_d13C)
plot(nat.CO2vDOC_d13C_cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")
abline(h = 4/sample_size2, col="red")  # add cutoff line
text(x=1:length(nat.CO2vDOC_d13C_cooksd)+1, y=nat.CO2vDOC_d13C_cooksd,
     labels=ifelse(nat.CO2vDOC_d13C_cooksd>4/sample_size2, names(nat.CO2vDOC_d13C_cooksd),""), col="red")  # add labels
# identifies:
# 02. GP05 08/09/2015

res.CO2vDOC_d13C_cooksd <- cooks.distance(res.CO2vDOC_d13C.lm)
sample_size3 = nrow(res.CO2vDOC_d13C)
plot(res.CO2vDOC_d13C_cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")
abline(h = 4/sample_size3, col="red")  # add cutoff line
text(x=1:length(res.CO2vDOC_d13C_cooksd)+1, y=res.CO2vDOC_d13C_cooksd,
     labels=ifelse(res.CO2vDOC_d13C_cooksd>4/sample_size3, names(res.CO2vDOC_d13C_cooksd),""), col="red")  # add labels
# identifies: NONE

# remove outliers
CO2vDOC_d13C_out = CO2vDOC_d13C[-c(2,9,26),]
nat.CO2vDOC_d13C_out = nat.CO2vDOC_d13C[-c(2),]

CO2vDOC_d13C_out.lm = lm(CO2vDOC_d13C_out$CO2_d13C ~ CO2vDOC_d13C_out$DOC_d13C)
summary(CO2vDOC_d13C_out.lm)$r.squared # R^2 = 0.4441137
summary(CO2vDOC_d13C_out.lm)$coefficients[,4] # p-value = 0.0002015675***
summary(CO2vDOC_d13C_out.lm)$coefficients[] # y = 1.630147x + 30.449111 

nat.CO2vDOC_d13C_out.lm = lm(nat.CO2vDOC_d13C_out$CO2_d13C ~ nat.CO2vDOC_d13C_out$DOC_d13C)
summary(nat.CO2vDOC_d13C_out.lm)$r.squared # R^2 = 0.004816088
summary(nat.CO2vDOC_d13C_out.lm)$coefficients[,4] # p-value = 0.7912721 
summary(nat.CO2vDOC_d13C_out.lm)$coefficients[] # y = -0.08496646x - 14.74058234  

res.CO2vDOC_d13C.lm = lm(res.CO2vDOC_d13C$CO2_d13C ~ res.CO2vDOC_d13C$DOC_d13C)
summary(res.CO2vDOC_d13C_out.lm)$r.squared # R^2 = 0.5194643
summary(res.CO2vDOC_d13C.lm)$coefficients[,4] # p-value = 0.01233551***
summary(res.CO2vDOC_d13C.lm)$coefficients[] # y = 4.382868x + 108.692992

Fig_CO2vDOC_d13C_out = ggplot(data = CO2vDOC_d13C_out, aes(x = DOC_d13C, y = CO2_d13C, group = Pool_type, color = Pool_type)) + 
  geom_point() +
  geom_smooth(method=lm, se = FALSE, linewidth = 0.5, col = "mediumblue", aes(group=FALSE)) +
  geom_smooth(data = subset(CO2vDOC_d13C_out, Pool_type=="Restoration"), method=lm, se = FALSE, linewidth = 0.5) +
  xlab(expression(paste(delta^13*'C-DOC (\u2030)'))) + ylab(expression(paste(delta^13*'C-CO'[2]*' (\u2030)'))) +
  scale_x_continuous(limits = c(-32,-18), breaks = c(-30,-25,-20)) +
  scale_y_continuous(limits = c(-22,-8), breaks = c(-20,-15,-10)) +
  scale_color_manual(values = c("Natural" = "black", "Restoration" = "grey65")) +
  annotate("text", x=-22.5, y=-11, size=2.5, label=expression('y = 1.6x + 30.5'), color = "mediumblue") +
  annotate("text", x=-25.5, y=-17, size=2.5, label=expression('y = 4.4x + 108.7'), color = "grey65") +
  annotate("text", x=-21, y=-19.5, size=2.5, label=expression('R'^2*' = 0.30, p = 0.000*'), color = "mediumblue") +
  annotate("text", x=-21, y=-20.5, size=2.5, label=expression('R'^2*' = 0.01, p = 0.791 '), color = "black") +
  annotate("text", x=-21, y=-21.5, size=2.5, label=expression('R'^2*' = 0.52, p = 0.012*'), color = "grey65") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position = "none")
Fig_CO2vDOC_d13C_out


# F14C

# identify potential outliers

CO2vDOC_F14C_cooksd <- cooks.distance(CO2vDOC_F14C.lm)
sample_size4 = nrow(CO2vDOC_F14C)
plot(CO2vDOC_F14C_cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")
abline(h = 4/sample_size4, col="red")  # add cutoff line
text(x=1:length(CO2vDOC_F14C_cooksd)+1, y=CO2vDOC_F14C_cooksd,
     labels=ifelse(CO2vDOC_F14C_cooksd>4/sample_size4, names(CO2vDOC_F14C_cooksd),""), col="red")  # add labels
# identifies:
# 07. L11 04/11/2014
# 13. P06 20/05/2014

nat.CO2vDOC_F14C_cooksd <- cooks.distance(nat.CO2vDOC_F14C.lm)
sample_size5 = nrow(nat.CO2vDOC_F14C)
plot(nat.CO2vDOC_F14C_cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")
abline(h = 4/sample_size5, col="red")  # add cutoff line
text(x=1:length(nat.CO2vDOC_F14C_cooksd)+1, y=nat.CO2vDOC_F14C_cooksd,
     labels=ifelse(nat.CO2vDOC_F14C_cooksd>4/sample_size5, names(nat.CO2vDOC_F14C_cooksd),""), col="red")  # add labels
# identifies:
# 11. P04 05/11/2014

res.CO2vDOC_F14C_cooksd <- cooks.distance(res.CO2vDOC_F14C.lm)
sample_size6 = nrow(res.CO2vDOC_F14C)
plot(res.CO2vDOC_F14C_cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")
abline(h = 4/sample_size6, col="red")  # add cutoff line
text(x=1:length(res.CO2vDOC_F14C_cooksd)+1, y=res.CO2vDOC_F14C_cooksd,
     labels=ifelse(res.CO2vDOC_F14C_cooksd>4/sample_size6, names(res.CO2vDOC_F14C_cooksd),""), col="red")  # add labels
# identifies:
# 02. L11 04/11/2014

# remove outliers
CO2vDOC_F14C_out = CO2vDOC_F14C[-c(7,13),]
nat.CO2vDOC_F14C_out = nat.CO2vDOC_F14C[-c(11),]
res.CO2vDOC_F14C_out = res.CO2vDOC_F14C[-c(2),]

CO2vDOC_F14C_out.lm = lm(CO2vDOC_F14C_out$CO2_F14C ~ CO2vDOC_F14C_out$DOC_F14C)
summary(CO2vDOC_F14C_out.lm)$r.squared # R^2 = 0.6725614
summary(CO2vDOC_F14C_out.lm)$coefficients[,4] # p-value = 9.243382e-08***
summary(CO2vDOC_F14C_out.lm)$coefficients[] # y = 0.3606807x + 0.6486170 

nat.CO2vDOC_F14C_out.lm = lm(nat.CO2vDOC_F14C_out$CO2_F14C ~ nat.CO2vDOC_F14C_out$DOC_F14C)
summary(nat.CO2vDOC_F14C_out.lm)$r.squared # R^2 = 0.6432995
summary(nat.CO2vDOC_F14C_out.lm)$coefficients[,4] # p-value = 6.233836e-05***
summary(nat.CO2vDOC_F14C_out.lm)$coefficients[] # y = 0.3406241x + 0.6686759 

res.CO2vDOC_F14C_out.lm = lm(res.CO2vDOC_F14C_out$CO2_F14C ~ res.CO2vDOC_F14C_out$DOC_F14C)
summary(res.CO2vDOC_F14C_out.lm)$r.squared # R^2 = 0.4768039
summary(res.CO2vDOC_F14C_out.lm)$coefficients[,4] # p-value = 0.027069122*** 
summary(res.CO2vDOC_F14C_out.lm)$coefficients[] # y = 0.3938443x + 0.6146593  

Fig_CO2vDOC_F14C_out = ggplot(data = CO2vDOC_F14C_out, aes(x = DOC_F14C, y = CO2_F14C, group = Pool_type, color = Pool_type)) + 
  geom_point() + 
  geom_smooth(method=lm, se = FALSE, linewidth = 0.5, col = "mediumblue", aes(group=FALSE)) +
  geom_smooth(method=lm, se = FALSE, linewidth = 0.5) +
  xlab(expression('F'^14*'C-DOC')) + ylab(expression('F'^14*'C-CO'[2])) +
  scale_x_continuous(limits = c(0.9, 1.2), breaks = c(1.2, 1.1, 1.0, 0.9)) +
  scale_y_continuous(limits = c(0.9, 1.2), breaks = c(1.2, 1.1, 1.0, 0.9)) +
  scale_color_manual(values = c("Natural" = "black", "Restoration" = "grey65")) +
  annotate("text", x=0.97, y=0.96, size=2.5, label=expression('y = 0.3606x + 0.6486'), color = "mediumblue") +
  annotate("text", x=1.14, y=1.19, size=2.5, label=expression('R'^2*' = 0.67, p = 0.000*'), color = "mediumblue") +
  annotate("text", x=0.97, y=0.94, size=2.5, label=expression('y = 0.3406x + 0.6687'), color = "black") +
  annotate("text", x=1.14, y=1.17, size=2.5, label=expression('R'^2*' = 0.64, p = 0.000*'), color = "black") +
  annotate("text", x=0.97, y=0.92, size=2.5, label=expression('y = 0.3938x + 0.6147'), color = "grey65") +
  annotate("text", x=1.14, y=1.15, size=2.5, label=expression('R'^2*' = 0.48, p = 0.027*'), color = "grey65") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title=element_blank()) + theme(legend.position = c(0.23, 0.89)) +
  theme(legend.spacing.y = unit(0.05, 'cm')) +
  theme(legend.background = element_rect(fill='transparent'))
Fig_CO2vDOC_F14C_out


# combine panels

plot_grid(Fig_CO2vDOC_F14C_out, Fig_CO2vDOC_d13C_out,
          labels = c("a", "b"),
          ncol = 2, nrow=1)

ggsave("Fig_S7.jpeg", width = 15, height = 7, units = c("cm"), dpi = 600)


####################### END #######################

