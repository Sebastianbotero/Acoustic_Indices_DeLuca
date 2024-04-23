library(ggplot2)
library(ggpubr)
library (iNEXT)
library(reshape2)
library(vegan)
library(caret)
library(corrplot)
library(RColorBrewer)


################ EFunction to estimate biodiversity at a determined sampling converge


biod.calc<-function(x,dates,SC){
  results<-data.frame(Site=NULL,Season=NULL,N=NULL,SampleCoverage=NULL,ObservedS=NULL,Hill_0=NULL,Hill_1=NULL,Hill_2=NULL,sampleCoverage_at_Estimation=NULL)
  community<-data.frame(Season=NA,Site=NA, NMDS1=NA,NMDS2=NA)
  
  for (i in 1:length(dates[,1])){
    
    x1<-x[which(x$Date>=dates[i,2] & x$Date<dates[i,3]),]
    
    
    
    
    for (s in 1:8){
      
      xs<-x1[x1$Site==s,]
      if(length(xs[,1])==0){
        results1<-data.frame(Site=NA,Season=NA,N=NA,SampleCoverage=NA,ObservedS=NA,Hill_0=NA,Hill_1=NA,Hill_2=NA,sampleCoverage_at_Estimation=NA)
        results<-rbind(results, results1)} else {
          
          com<-aggregate(N_individuals ~ Species, data =xs, sum)
          abundance<-com$N_individuals
          
          out1 <- iNEXT(abundance, q=c(0,1,2), datatype="abundance", endpoint = 3*sum(com$N_individuals))
          
          closer<-min(abs(out1$iNextEst$coverage_based$SC-SC))
          Res_SC<-out1$iNextEst$coverage_based[which(abs(out1$iNextEst$coverage_based$SC-SC)==closer),]
          results1<-data.frame(Site=s,Season=dates[i,1],N=out1$DataInfo$n,SampleCoverage=out1$DataInfo$SC,ObservedS=out1$DataInfo$S.obs,
                               Hill_0=Res_SC[Res_SC$Order.q==0,"qD"][1],Hill_1=Res_SC[Res_SC$Order.q==1,"qD"][1],Hill_2=Res_SC[Res_SC$Order.q==2,"qD"][1],
                               sampleCoverage_at_Estimation=Res_SC$SC[1])
          
          
          results<-rbind(results, results1)}
        }
      
      
      comm1<-dcast(data=x1, Site ~ Species, value.var= "N_individuals", fun.aggregate=sum)
      comm<-as.matrix(comm1[,-1])
      row.names(comm)<-comm1$Site
      names(comm)<-names(comm1[-1])
      
      Commdist<-vegdist(wisconsin(comm),"horn")
      Commclus <- hclust(Commdist, "average")
      
      NMDS=metaMDS(wisconsin(comm),trymax=100)
      nmdsAxis<-scores(NMDS)
      community1<-data.frame(Season=dates[i,1],Site=comm1$Site, NMDS1=nmdsAxis$sites[,1],NMDS2=nmdsAxis$sites[,2])
      community<-rbind(community,community1)
      
      

      
      
     
  }
  
  results<-merge(results,community, by= c("Season","Site"))
  return(results)
}

################### Dates for each season


dates<-data.frame(season=c("ALL","Winter","Spring","Summer"),ini=as.Date(c("Jan_1_2022","Jan_1_2022","Mar_20_2022","June_21_2022"),
                                                                         format='%b_%d_%Y'),end=as.Date(c("Aug_30_2022","Mar_20_2022","June_21_2022","Aug_30_2022"),
                                                                                                        format='%b_%d_%Y'))
################################################################
################## Bird diversity by site x season ##############
##################################################################
birds<-read.csv("C://Users//seboc//Box//DeLuca//All_Bird_Data_SciNames.csv")
head(birds)
birds$Date<-paste(birds$Date,"2022", sep="_")
birds$Date<-as.Date(birds$Date,format='%b_%d_%Y')
birds<-birds[birds$Estimate_Method=="Abundance",] # only use abundance method
names(birds)[9]<-"N_individuals"

mean(c(0.9745,0.9745,0.9899,0.9812,0.9706,0.9872,0.9742,0.9590)) ## We will compare biodiversity estimates for each site at the mean sampling coverage (SC)
Birds<-biod.calc(x=birds,dates=dates,SC=0.976)

write.csv(Birds,"Bird_Diversity.csv")




################################################################
################## Herps diversity by site x season ##############
##################################################################
HerpsD<-read.csv("C:/Users/seboc/Box/DeLuca/Final_Databases/Drift_Fence_Herps.csv")
HerpsD$N_individuals<-1
HerpsV<-read.csv("C:/Users/seboc/Box/DeLuca/Final_Databases/Visual_Herps.csv")
HerpsV$Species<-HerpsV$SPECIES
head(HerpsD)
head(HerpsV)

Herps<-rbind(HerpsD[,c(5,6,8,9,20)],HerpsV[,c(2,6,13,8,11)])
head(Herps)
Herps$Date<-paste(Herps$Date,"2022", sep="_")
Herps$Date<-as.Date(Herps$Date,format='%b_%d_%Y')
Herps<-Herps[Herps$Site %in% c("1","2","3","4","5","6","7","8"),]
Herps$Site<-as.numeric(Herps$Site)

Amphibians<-Herps[Herps$Class=="Amphibia",]
unique(Amphibians$Species)


Reptiles<-Herps[Herps$Class=="Reptiles",]
unique(Reptiles$Species)


mean(c(1,0.85,1,1,0.95,0.96,0.96))
Amphibia<-biod.calc(x=Amphibians,dates=dates,SC=0.96)
write.csv(Amphibia,"Amphibians_Diversity.csv")


mean(c(0.902,0.96,0.94,0.95,0.86,0.93,0.86,0.89))
Rep<-biod.calc(x=Reptiles,dates=dates,SC=0.91)
write.csv(Rep,"Reptiles_Diversity.csv")

#########################################################################################################################
############################  Mammals ###################################################################################
################################################################################################################

CTs<-read.csv("C:/Users/seboc/Box/DeLuca/Final_Databases/Camera_trapping_NoHumanact/images.csv")
head(CTs)
CTs<-CTs[CTs$is_blank==0,]
CTs<-CTs[CTs$class=="Mammalia",]
CTs$Species<-paste(CTs$genus, CTs$species, sep=" ")
CTs<-CTs[!CTs$species=="",]
CTs<-CTs[!CTs$Species %in% c(  "Bos taurus", "Equus caballus" ,"Canis familiaris" ,"Sigmodon hispidus"),]
CTs<-CTs[CTs$Site %in% 1:8,]
CTs$Date<-as.Date(CTs$timestamp,format='%m/%d/%Y')

##### asume records to be independent per day
#CTind<-aggregate(number_of_objects ~ Species, data =CTs, sum)
CTind<-aggregate(number_of_objects ~ Species + Site + Date, data =CTs, max)
names(CTind)[4]<-"N_individuals"

mean(c(1,1,0.99,1,1,0.98,1,0.99))
Mammals<-biod.calc(x=CTind,dates=dates,SC=0.99)
write.csv(Mammals,"MammalsCT_Diversity.csv")




#########################################################################################################################
############################  Small terrestrial mammals ###################################################################################
################################################################################################################


STM<-read.csv("C:/Users/seboc/Box/DeLuca/Final_Databases/Mammal_trapping.csv")
head(STM)
STM<-STM[!STM$Species %in% c(NA,"Spilogale putorius"),]
STM$Date<-as.Date(STM$Date,format='%b_%d_%Y')
STM$N_individuals<-1
#STM<-aggregate(N_individuals ~ Species, data =STM, sum)
S.Mammals<-biod.calc(x=STM,dates=dates,SC=0.99)
write.csv(S.Mammals,"Small_Mammals_Diversity.csv")


#########################################################################################################################
############################  All groups ###################################################################################
################################################################################################################


All.verts<-rbind(birds[,c("Date","Site", "Species", "N_individuals")],Herps[,c("Date","Site", "Species", "N_individuals")],
                 CTind[,c("Date","Site", "Species", "N_individuals")])
head(All.verts)
dim(All.verts)

mean(c(0.99,0.98,0.99,0.99,0.97,0.98,0.98,0.97))
all.biod<-biod.calc(x=All.verts,dates=dates,SC=0.98)
write.csv(all.biod,"All_verts_Diversity.csv")


########################  Explore diversity among sites


res<-list.files()
NoAco_div<-list(ALL=read.csv("All_verts_Diversity.csv" ),Amphibians=read.csv("Amphibians_Diversity.csv"),
                Reptiles=read.csv("Reptiles_Diversity.csv" ),MammalsCT=read.csv( "MammalsCT_Diversity.csv" ),
                Birds=read.csv("Bird_Diversity.csv" ))

#### Coeficient of variation
for (i in 1:5){
  
  test<-NoAco_div[[i]]
  
  print(names(NoAco_div)[i])
  print(sapply(test[1:8,c(9,8,7)],function(x){sd(x) / mean(x) * 100}))
  
}


###### Plot
layout(matrix(c(1:15), 5, 3, byrow = T))

par(mar=c(5,2,1,2) + 0.1)

for (i in 1:5){
  test<-NoAco_div[[i]]
  test$Site[1:8]<-paste(test$Site[1:8],c("Forest","Forest","Orchard","Orchard","Scrub","Scrub","Wetland","Wetland"),sep="_")
  barplot(height=test[1:8,7], names=test[1:8,3], main=paste(names(NoAco_div)[i],names(test)[7],sep="_"))
  barplot(height=test[1:8,8], names=test[1:8,3], main=paste(names(NoAco_div)[i],names(test)[8],sep="_"))
  barplot(height=test[1:8,9], names=test[1:8,3], main=paste(names(NoAco_div)[i],names(test)[9],sep="_"))
  
}


########## Correlation among groups
tiff(file="Correlation_Spearman_Apr21_50mBirds.tif",
     width=11.326, height=8.8, units="in", res=1000)

layout(matrix(c(1,2,3), nrow = 1, ncol = 3, byrow = TRUE))
par(oma=c(0,0,2,0),mai=c(0.5, 0, 0.1, 0),mgp=c(1,1,0))
## Hill 0

H0<-NoAco_div[[1]][1:8,c(7,3)]
names(H0)[1]<-names(NoAco_div)[1]

for (i in 2:5) {
  H0a<-NoAco_div[[i]][1:8,c(7,3)]
  names(H0a)[1]<-names(NoAco_div)[i]
  H0<-merge(H0,H0a, by="Site")
  
  
}






corrtax<-function(H) {
  
  M<-matrix(nrow=5,ncol=5)
  row.names(M)<-names(H0[,-1])
  colnames(M)<-names(H0[,-1])
  M.p<-M
  for (i in 1:5){
    for (e in 1:5){
      corr<-corci(H[,names(H)[i]],H[,names(H)[e]],method = "spearman",nboot=5000)
      M[i,e] <-corr$estimate
      M.p[i,e] <-corr$p.value
    }
    
  }
  return(list(correlation=M,pval=M.p))
  }


Hill0<-corrtax(H0[,-1])

#testRes = cor.mtest(H0[,-1], conf.level = 0.95,method="spearman")
#M<-cor(H0[,-1], method="spearman")
corrplot(Hill0$correlation, type="lower", order='original',tl.col = "black",tl.cex=1.3,cl.cex=1.1,
         col=rev(brewer.pal(n=8, name="RdYlBu"))) #, p.mat = Hill0$pval, sig.level = 0.05,insig='label_sig')
text(2,8,"Species richness", adj=c(0.5,1),cex=1.9, font=2)

H1<-NoAco_div[[1]][1:8,c(8,3)]
names(H1)[1]<-names(NoAco_div)[1]

for (i in 2:5) {
  H1a<-NoAco_div[[i]][1:8,c(8,3)]
  names(H1a)[1]<-names(NoAco_div)[i]
  H1<-merge(H1,H1a, by="Site")
  
  
}

#M<-cor(H1[,-1],method="spearman")
#testRes = cor.mtest(H1[,-1], conf.level = 0.95,method="spearman")
Hill1<-corrtax(H1[,-1])

corrplot(Hill1$correlation, type="lower", order='original',tl.col = "black",tl.cex=1.3,cl.cex=1.1,
         col=rev(brewer.pal(n=8, name="RdYlBu"))) #, p.mat = Hill1$pval, sig.level = 0.05,insig='label_sig')
text(2,8,"Shannon's div.", adj=c(0.5,1),cex=1.9, font=2)



H2<-NoAco_div[[1]][1:8,c(11,3)]

names(H2)[1]<-names(NoAco_div)[1]

for (i in 2:5) {
  H2a<-NoAco_div[[i]][1:8,c(11,3)]
  names(H2a)[1]<-names(NoAco_div)[i]
  H2<-merge(H2,H2a, by="Site")
  
  
}

#M<-cor(H2[,-1],method="spearman")
#testRes = cor.mtest(H1[,-1], conf.level = 0.95,method="spearman")

nmds<-corrtax(H2[,-1])
corrplot(nmds$correlation, type="lower", order='original',tl.col = "black",tl.cex=1.3,cl.cex=1.1,
         col=rev(brewer.pal(n=8, name="RdYlBu")))#, p.mat = nmds$pval, sig.level = 0.005,insig='label_sig')
text(2,8,"NMDS I", adj=c(0.5,1),cex=1.9, font=2)
dev.off()

##################################################################
################### Ecoacoustic Indexes ######################
##########################################################
# Aggregate each index by average of 5 min samples from each recording for temporal bins

# Import estimates for each recording

Winter<-read.csv("C:/Users/seboc/Box/Ecoacustics_paper/AI_5minWinter_Final_Jan17.csv")
head(Winter)
Winter$Season<-"Winter"
Summer<-read.csv("C:/Users/seboc/Box/Ecoacustics_paper/AI_5minSummerAR_Final_Jan_16.csv")
head(Summer)
Summer$Season<-"Summer"
Spring<-read.csv("C:/Users/seboc/Box/Ecoacustics_paper/AI_5minSpring_Final_Jan14.csv")
Spring$Season<-"Spring" 
head(Spring)



All_AI<-rbind(Winter, Summer, Spring)
head(All_AI)
All_AI$SiteGeneral<-NA
All_AI[All_AI$Site %in% c("AR-1A","AR-1B" ),"SiteGeneral"]<-1
All_AI[All_AI$Site %in% c("AR-2A","AR-2B" ),"SiteGeneral"]<-2
All_AI[All_AI$Site %in% c("AR-3A","AR-3B" ),"SiteGeneral"]<-3
All_AI[All_AI$Site %in% c("AR-4A","AR-4B" ),"SiteGeneral"]<-4
All_AI[All_AI$Site %in% c("AR-5A","AR-5B" ),"SiteGeneral"]<-5
All_AI[All_AI$Site %in% c("AR-6A","AR-6B" ),"SiteGeneral"]<-6
All_AI[All_AI$Site %in% c("AR-7A","AR-7B" ),"SiteGeneral"]<-7
All_AI[All_AI$Site %in% c("AR-8A","AR-8B", "AR-8A-Jan19-Jan20"),"SiteGeneral"]<-8


### Add the first pca axis as a EA variable
pc <- prcomp(x = All_AI[ , -c(1:7,15,16)],
             center = TRUE, 
             scale. = TRUE)

print(pc)
summary(pc)
All_AI$PCA_I<-pc$x[, 1]

#### Time

xx<-strsplit(All_AI$Date,"_")
All_AI$time<-sapply(xx, function(x){x[length(x)]})
All_AI$time<-gsub(".wav","", All_AI$time)
All_AI$time<-as.numeric(All_AI$time)

hist(All_AI$time,breaks=100)

All_AI_cur<-All_AI[-which(All_AI$time<30000 | All_AI$time>230000),]
All_AI_cur$TimeofDay<-"Morning"
All_AI_cur$TimeofDay[All_AI_cur$time>120000]<-"Evening"

##################################################################################
################# Aggregate for each time bin using the average index value
######################################################################

All_AI_cur_Mstandart<- All_AI_cur
All_AI_cur_Mstandart$M <-  All_AI_cur_Mstandart$M/ max( All_AI_cur_Mstandart$M)
AI_season<-aggregate( .~  SiteGeneral + Season, data =All_AI_cur_Mstandart[,8:17], mean)

AI_winter<-AI_season[AI_season$Season=="Winter",]
names(AI_winter)[-1]<-paste(names(AI_winter)[-1],"_winter", sep="")
AI_Spring<-AI_season[AI_season$Season=="Spring",]
names(AI_Spring)[-1]<-paste(names(AI_Spring)[-1],"_spring", sep="")
AI_Summer<-AI_season[AI_season$Season=="Summer",]
names(AI_Summer)[-1]<-paste(names(AI_Summer)[-1],"_summer", sep="")


AI_TofDay<-aggregate( .~  SiteGeneral + TimeofDay, data =All_AI_cur_Mstandart[,c(8:14,16,17,19)], mean)

AI_morning<-AI_TofDay[AI_TofDay$TimeofDay=="Morning",]
names(AI_morning)[-1]<-paste(names(AI_morning)[-1],"_morning", sep="")
AI_evening<-AI_TofDay[AI_TofDay$TimeofDay=="Evening",]
names(AI_evening)[-1]<-paste(names(AI_evening)[-1],"_evening", sep="")

AI_All<-aggregate( .~  SiteGeneral , data =All_AI_cur_Mstandart[,c(8:14,16,17)], mean)

all.data<-list(AI_winter,AI_Spring,AI_Summer,AI_morning, AI_evening,AI_All )

finalDB<-all.data[[1]]

for (i in 2:6){
  
  finalDB<-merge(finalDB,all.data[[i]],by="SiteGeneral")
  
  
}

names(finalDB)

finalDB<-finalDB[,-c(2,11,20,29,38)]


############ Merge with non acoustic estimates
ALL<-merge(finalDB ,NoAco_div[[1]][1:8,c(3,7,8,11)], by.x="SiteGeneral", by.y="Site")[,-1]

ornames<-names(ALL)[1:48]
replace<-c("BIO_winter",       "ACI_winter",       "NDSI_winter",      "ADI_winter",      
           "AEI_winter",       "H_winter",         "M_winter",         "PCAI_winter",    
           "BIO_spring",       "ACI_spring",       "NDSI_spring",      "ADI_spring",      
           "AEI_spring",       "H_spring",         "M_spring",         "PCAI_spring",    
           "BIO_summer",       "ACI_summer",       "NDSI_summer",      "ADI_summer",      
           "AEI_summer",       "H_summer",         "M_summer",         "PCAI_summer",    
           "BIO_morning",      "ACI_morning",      "NDSI_morning",     "ADI_morning",     
           "AEI_morning",      "H_morning",        "M_morning",        "PCAI_morning",   
           "BIO_evening",      "ACI_evening",      "NDSI_evening",     "ADI_evening",     
           "AEI_evening",      "H_evening",        "M_evening",        "PCAI_evening",   
           "BIO_all", "ACI_all",  "NDSI_all",             "ADI_all",             
           "AEI_all", "H_all","M_all",                "PCAI_all")




names(ALL)[1:48]<-replace
EAindeces<-replace

div<-c("Hill_0" , "Hill_1", "NMDS1")

### Fuction for assessing correlation

corr<-function(dat,EAindeces,div){
  Corr_results<-data.frame(GroundDiv=NULL, IndCom=NULL, Index=NULL ,Time=NULL, rho=NULL,p=NULL,up=NULL,low=NULL, lm_beta=NULL,lm_pval=NULL,lm_rsq=NULL)
  ctrl <- trainControl(method = "LOOCV", number = 8)
  fm<-as.formula( "x~.")
  for (i in 1:length(div)){
    
    x<-dat[,div[i]]
    
    for(e in 1:length(EAindeces)){
      
      y=dat[,EAindeces[e]]
      
      z<-corci(x,y,method = "spearman",nboot=5000)
      
      
      
      mod<-lm(fm, data=data.frame(x,y))
      model <- train(fm, data=data.frame(x,y), method = "lm", trControl = ctrl)
      
      
      Corr_results1<-data.frame(GroundDiv=div[i],IndCom=EAindeces[e], Index=strsplit(EAindeces[e],"_")[[1]][1] ,
                                Time=strsplit(EAindeces[e],"_")[[1]][2],
                                rho=z$estimate,p=z$p.value,up=z$conf.int[2],low=z$conf.int[1],lm_beta=summary(mod)$coefficients[2,1],lm_pval=summary(mod)$coefficients[2,4],lm_rsq=model$results$Rsquared)
      
      Corr_results<-rbind(Corr_results,Corr_results1)
    }
    
    
  }
  
  return(Corr_results)
} 







#### All vertebrates
all.corrs<-corr(ALL,EAindeces,div)

all.corrs[all.corrs$p<=0.05,]

write.csv(all.corrs,"Spearman_All_verts.csv")



################################ Birds #######################


birds<-merge(finalDB,NoAco_div[[5]][1:8,c(3,7,8,11)],by.x="SiteGeneral",by.y="Site")[,-1]
names(birds)[1:48]<-replace
EAindeces<-names(birds)[1:48]
birds.corrs<-corr(birds,EAindeces,div)

birds.corrs[birds.corrs$p<0.06,]

write.csv(birds.corrs,"Spearman_Birds.csv")



################################ Mammals #######################
Mammals<-merge(finalDB,NoAco_div[[4]][1:8,c(3,7,8,11)],by.x="SiteGeneral",by.y="Site")[,-1]
names(Mammals)[1:48]<-replace
mammals.corrs<-corr(Mammals,EAindeces,div)

mammals.corrs[mammals.corrs$p<0.05,]

write.csv(mammals.corrs,"Spearman_Mammals.csv")


###################### Reptiles 
Reptiles<-merge(finalDB,NoAco_div[[3]][1:8,c(3,7,8,11)],by.x="SiteGeneral",by.y="Site")[,-1]
names(Reptiles)[1:48]<-replace

Reptiles.corrs<-corr(Reptiles,EAindeces,div)

Reptiles.corrs[Reptiles.corrs$p<0.05,]

write.csv(Reptiles.corrs,"Spearman_Reptiles.csv")


########################## Amphibians
Amphibians<-merge(finalDB,NoAco_div[[2]][1:8,c(3,7,8,11)],by.x="SiteGeneral",by.y="Site")[,-1]
names(Amphibians)[1:48]<-replace

Amphibians.corrs<-corr( Amphibians,EAindeces,div)

Amphibians.corrs[ Amphibians.corrs$p<0.05,]

write.csv( Amphibians.corrs,"Spearman_Amphibians.csv")

#####################################################
######################Plot results  -forest plot
#################################

#Function to plot confidence intervals
cints<-function(x,ind,col){for (i in 1:length(x[,1])){
  
  lines(x=x[i,c("low","up")],y=(c(x[i,"ind"],x[i,"ind"])*5)+ind,col= yarrr::transparent(col, trans.val = x[i,"col"]),lwd=1.8)
  
  
}}

#Function to create plots
forestplot<-function(corr.results,file){
  
  tiff(file=file,
       width=8.8, height=11.326, units="in", res=1000)
  
  
  
  corr.results$col<-0.65
  corr.results$col[corr.results$p<=0.05]<-0.1
  corr.results$font<-1
  corr.results$font[corr.results$p<=0.05]<-2
  corr.results<-corr.results[order(corr.results$IndCom),]
  
  indices<-data.frame(ind=unique(corr.results$Index),font=1)
  indices$font[indices$ind %in% corr.results$Index[corr.results$p <= 0.05]]  <-2
  
  
  
  dat<-corr.results[corr.results$GroundDiv=="Hill_0",]
  dat<-dat[order(dat$IndCom),]
  dat$ind<-48:1
  
  
  dat1<-corr.results[corr.results$GroundDiv=="Hill_1",]
  dat1<-dat1[order(dat1$IndCom),]
  dat1$ind<-48:1
  
  dat2<-corr.results[corr.results$GroundDiv=="NMDS1",]
  dat2<-dat2[order(dat2$IndCom),]
  dat2$ind<-48:1
  
  dat$font<-1
  dat$font[dat$p <=0.05 | dat1$p <=0.05 | dat2$p <=0.05 ] <- 2
  
  
  
  layout(matrix(c(1,2,2,3,3), nrow = 1, ncol = 5, byrow = TRUE))
  par(oma=c(0,0,2,0),mai=c(0.5, 0, 0.4, 0))
  plot.new()
  plot.window(ylim=c(1,240), xlim=c(0,1),oma=c(0,0,0,0),mar=c(0,0,0,0))
  text(y=c(3,9,15,21,27,33,39,45)*5, x=0.25,labels=rev(indices$ind),font=rev(indices$font), cex=3)
  text(y=rev(seq(5,240,by=5)),  x=0.9,labels=dat$Time,font=dat$font,cex=1.2,srt = 35, pos = 2, xpd = TRUE)
  
  par(mai=c(0.6, 0, 0.4, 0))
  plot.window(xlim=c(-1,1), ylim=c(1,144))
  plot(dat$rho,dat$ind*5,pch=19, cex=1.8,cex.main=2, main= "Diversity",col="purple",axes=F,xlab="rho",cex.lab=2,ylab=NA,xlim=c(-1,1))
  axis(1, at=c(-1,-0.5,0,0.5,1),cex.axis=1.5)
  cints(dat,0,"purple")
  abline(v=0,lwd=2,lty=2)
  axis(2,seq(5,240,by=5), labels=NA, tick = TRUE)
  
  
  
  points(dat1$rho,(dat1$ind*5)+1,pch=19, cex=1.8,col="darkgreen")
  cints(dat1,1,"darkgreen")
  
  par(mai=c(0.6, 0.2, 0.4, 0))
  plot.window(xlim=c(-1,1), ylim=c(1,144))
  plot(dat2$rho,dat2$ind*5,pch=19, cex=1.8, cex.main=2,main= "Composition",col="black",axes=F,xlab="rho",cex.lab=2,ylab=NA,xlim=c(-1,1))
  axis(1, at=c(-1,-0.5,0,0.5,1),cex.axis=1.5)
  cints(dat2,0,"black")
  abline(v=0,lwd=2,lty=2)
  
  
  
  dev.off()
}



######### Plot correlations for all vertebrates

all.verts<-read.csv("Spearman_All_verts.csv")
forestplot(corr.results=all.verts, file="Vertebrate_Community.tif")



######### Plot correlations for Birds

Birds<-read.csv("Spearman_Birds.csv")
forestplot(corr.results=Birds, file="Birds.tif")


######### Plot correlations for Amphibians

Amphibians<-read.csv("Spearman_Amphibians.csv")
forestplot(corr.results=Amphibians, file="Amphibians.tif")


######### Plot correlations for Reptiles

Reptiles<-read.csv("Spearman_Reptiles.csv")
forestplot(corr.results=Reptiles, file="Reptiles.tif")

######### Plot correlations for Mammals

Mammals<-read.csv("Spearman_Mammals.csv")
forestplot(corr.results=Mammals, file="Mammals.tif")








#################### Plot Best Models ########################

results[results$Delta==0,]
par( family="serif",mar=c(7,7,1,1))

### ALL Hill 0 

plot(ALL$M_evening_avg, ALL$Hill_0, xlab= "", pch=19, cex.axis=2,ylim=c(50,90),cex.lab=2.5, ylab=NA)
title(xlab = expression("M evening"),ylab = expression("Richness"), line = 3.5,cex.lab=2) 
mod<-lm(Hill_0 ~ M_evening_avg, data=ALL)
var<-ALL$M_evening
newx <- seq(min(var), max(var), length.out=100)
preds <- as.data.frame(predict(mod, newdata = data.frame(M_evening_avg=newx), interval = 'confidence'))

#add dashed lines for confidence bands
lines(newx, preds[ ,1], lty = 1,lwd=5, col = 'blue')
lines(newx, preds[ ,3], lty = 'dashed', col = 'blue')
lines(newx, preds[ ,2], lty = 'dashed', col = 'blue')
points(ALL$M_evening, ALL$Hill_0,  pch=19, cex=1.8)



################################################### ALL Hill 1


plot(ALL$BIO_winter, ALL$Hill_1, xlab= "", pch=19, cex.axis=2,ylim=c(20,50),cex.lab=2.5, ylab=NA)
title(xlab = expression("Bio winter"),ylab = expression("Shannon's div."), line = 3.5,cex.lab=2)

mod<-lm(Hill_1 ~ BIO_winter_avg, data=ALL)
var<-ALL$BIO_winter
newx <- seq(min(var), max(var), length.out=100)
preds <- as.data.frame(predict(mod, newdata = data.frame(BIO_winter_avg=newx), interval = 'confidence'))

#add dashed lines for confidence bands
lines(newx, preds[ ,1], lty = 1,lwd=5, col = 'blue')
lines(newx, preds[ ,3], lty = 'dashed', col = 'blue')
lines(newx, preds[ ,2], lty = 'dashed', col = 'blue')
points(ALL$BIO_winter, ALL$Hill_1,  pch=19, cex=1.8)

################################################### ALl Hill 2


plot(ALL$BIO_spring_90p, ALL$Hill_2, xlab= "", pch=19, cex.axis=2,ylim=c(9,35),cex.lab=2.5, ylab=NA)
title(xlab = expression("Bio Spring p"["90th"]), line = 3.8,cex.lab=2.5)

mod<-lm(Hill_2 ~ BIO_spring_90p, data=ALL)
var<-ALL$BIO_spring_90p
newx <- seq(min(var), max(var), length.out=100)
preds <- as.data.frame(predict(mod, newdata = data.frame(BIO_spring_90p=newx), interval = 'confidence'))

#add dashed lines for confidence bands
lines(newx, preds[ ,1], lty = 1,lwd=3, col = 'blue')
lines(newx, preds[ ,3], lty = 'dashed', col = 'blue')
lines(newx, preds[ ,2], lty = 'dashed', col = 'blue')
points(ALL$BIO_spring_90p, ALL$Hill_2,  pch=19, cex=1.8)

mtext(side=3,"Hill 0",2,adj=0.15, outer = TRUE, cex = 1.8)
mtext(side=3,"Hill 1",2,adj=0.52, outer = TRUE, cex = 1.8)
mtext(side=3,"Hill 2",2,adj=0.88, outer = TRUE, cex = 1.8)
mtext(side=2,"Effective number of species",2, outer = TRUE, cex = 2)
################################################### ALL NMDS1


plot(ALL$BIO_spring_avg, ALL$NMDS1, pch=19, cex.axis=2,cex.lab=2.5,xlab=NA, ylab=NA)
title(xlab = expression("Bio Spring"),ylab = expression("NMDS I"), line = 3.5,cex.lab=2)

mod<-lm(NMDS1 ~ BIO_spring_avg, data=ALL)
var<-ALL$BIO_spring_avg
newx <- seq(min(var), max(var), length.out=100)
preds <- as.data.frame(predict(mod, newdata = data.frame(BIO_spring_avg=newx), interval = 'confidence'))

#add dashed lines for confidence bands
lines(newx, preds[ ,1], lty = 1,lwd=5, col = 'blue')
lines(newx, preds[ ,3], lty = 'dashed', col = 'blue')
lines(newx, preds[ ,2], lty = 'dashed', col = 'blue')
points(ALL$BIO_spring_avg, ALL$NMDS1,  pch=19, cex=1.8)

############################## ALL Other groups#####
############################################################

#par(mfrow=c(4,3),mgp=c(2,1,0), tcl=-0.5, family="serif", omi=c(0.5,1.2,1,0.5),mar=c(5,3.5,0,0))
par(mar=c(5,5,1,1))
#### Birds
#H 0
plot(birds$M_winter, birds$Hill_0, xlab= "", pch=19, cex.axis=1.8,ylim=c(25,55),cex.lab=2, ylab=NA)
title(xlab = expression("M winter Avg"),ylab = expression("Richness"), line = 3.5,cex.lab=2) 
mod<-lm(Hill_0~M_winter, data=birds)
var<-birds$M_winter
newx <- seq(min(var), max(var), length.out=100)
preds <- as.data.frame(predict(mod, newdata = data.frame(M_winter=newx), interval = 'confidence'))

#add dashed lines for confidence bands
lines(newx, preds[ ,1], lty = 1,lwd=5, col = 'blue')
lines(newx, preds[ ,3], lty = 'dashed', col = 'blue')
lines(newx, preds[ ,2], lty = 'dashed', col = 'blue')
points(birds$M_winter, birds$Hill_0,  pch=19, cex=1.8)
par(mar=c(5,5,1,1))

#H 1
plot(birds$BIO, birds$Hill_1, xlab= "", pch=19, cex.axis=1.8,ylim=c(15,30),cex.lab=2, ylab=NA)
title(xlab = expression("BIO"),ylab = expression("Shannon's div."), line = 3.5,cex.lab=2) 
mod<-lm(Hill_1~BIO, data=birds)
var<-birds$BIO
newx <- seq(min(var), max(var), length.out=100)
preds <- as.data.frame(predict(mod, newdata = data.frame(BIO=newx), interval = 'confidence'))

#add dashed lines for confidence bands
lines(newx, preds[ ,1], lty = 1,lwd=5, col = 'blue')
lines(newx, preds[ ,3], lty = 'dashed', col = 'blue')
lines(newx, preds[ ,2], lty = 'dashed', col = 'blue')
points(birds$BIO, birds$Hill_1,  pch=19, cex=1.8)

#H 2
plot(birds$BIO_summer_90p, birds$Hill_2, xlab= "", pch=19, cex.axis=1.8,ylim=c(5,27),cex.lab=2, ylab=NA)
title(xlab = expression("BIO Summer p"["90th"]), line = 3.5,cex.lab=2) 
mod<-lm(Hill_2~BIO_summer_90p, data=birds)
var<-birds$BIO_summer_90p
newx <- seq(min(var), max(var), length.out=100)
preds <- as.data.frame(predict(mod, newdata = data.frame(BIO_summer_90p=newx), interval = 'confidence'))

#add dashed lines for confidence bands
lines(newx, preds[ ,1], lty = 1,lwd=3, col = 'blue')
lines(newx, preds[ ,3], lty = 'dashed', col = 'blue')
lines(newx, preds[ ,2], lty = 'dashed', col = 'blue')
points(birds$BIO_summer_90p, birds$Hill_2,  pch=19, cex=1.8)


############# Amphibians

#H 0
plot(Amphibians$BIO_evening_avg, Amphibians$Hill_0, xlab= "", pch=19, cex.axis=1.8,ylim=c(1,15),cex.lab=2, ylab=NA)
title(xlab = expression("BIO evening avg"),ylab = expression("Richness"), line = 3.5,cex.lab=2) 
mod<-lm(Hill_0 ~ BIO_evening_avg, data=Amphibians)
var<-Amphibians$BIO_evening_avg
newx <- seq(min(var), max(var), length.out=100)
preds <- as.data.frame(predict(mod, newdata = data.frame(BIO_evening_avg=newx), interval = 'confidence'))

#add dashed lines for confidence bands
lines(newx, preds[ ,1], lty = 1,lwd=5, col = 'blue')
lines(newx, preds[ ,3], lty = 'dashed', col = 'blue')
lines(newx, preds[ ,2], lty = 'dashed', col = 'blue')
points(Amphibians$BIO_evening_avg, Amphibians$Hill_0,  pch=19, cex=1.8)

#H 1 

plot(Amphibians$BIO_evening_avg, Amphibians$Hill_1, xlab= "", pch=19, cex.axis=1.8,ylim=c(1,10),cex.lab=2, ylab=NA)
title(xlab = expression("BIO evening Avg."),ylab = expression("Shannon's div."), line = 3.5,cex.lab=2) 
mod<-lm(Hill_1~BIO_evening_avg, data=Amphibians)
var<-Amphibians$BIO_evening_avg
newx <- seq(min(var), max(var), length.out=100)
preds <- as.data.frame(predict(mod, newdata = data.frame(BIO_evening_avg=newx), interval = 'confidence'))

#add dashed lines for confidence bands
lines(newx, preds[ ,1], lty = 1,lwd=5, col = 'blue')
lines(newx, preds[ ,3], lty = 'dashed', col = 'blue')
lines(newx, preds[ ,2], lty = 'dashed', col = 'blue')
points(Amphibians$BIO_evening_avg, Amphibians$Hill_1,  pch=19, cex=1.8)

#H 2

plot(Amphibians$BIO_evening, Amphibians$Hill_2, xlab= "", pch=19, cex.axis=1.8,ylim=c(1,7),cex.lab=2, ylab=NA)
title(xlab = expression("BIO evening Avg."), line = 3.5,cex.lab=2) 
mod<-lm(Hill_2~BIO_evening, data=Amphibians)
var<-Amphibians$BIO_evening
newx <- seq(min(var), max(var), length.out=100)
preds <- as.data.frame(predict(mod, newdata = data.frame(BIO_evening=newx), interval = 'confidence'))

#add dashed lines for confidence bands
lines(newx, preds[ ,1], lty = 1,lwd=3, col = 'blue')
lines(newx, preds[ ,3], lty = 'dashed', col = 'blue')
lines(newx, preds[ ,2], lty = 'dashed', col = 'blue')
points(Amphibians$BIO_evening, Amphibians$Hill_2,  pch=19, cex=1.8)


################################ Reptiles
#H 0
plot(Reptiles$ACI_summer, Reptiles$Hill_0, xlab= "", pch=19, cex.axis=1.8,ylim=c(1,15),cex.lab=2, ylab=NA)
title(xlab = expression("ACI_summer Avg."),ylab = expression("Richness"), line = 3.5,cex.lab=2) 
mod<-lm(Hill_0~ACI_summer, data=Reptiles)
var<-Reptiles$ACI_summer
newx <- seq(min(var), max(var), length.out=100)
preds <- as.data.frame(predict(mod, newdata = data.frame(ACI_summer=newx), interval = 'confidence'))

#add dashed lines for confidence bands
lines(newx, preds[ ,1], lty = 1,lwd=5, col = 'blue')
lines(newx, preds[ ,3], lty = 'dashed', col = 'blue')
lines(newx, preds[ ,2], lty = 'dashed', col = 'blue')
points(Reptiles$ACI_summer, Reptiles$Hill_0,  pch=19, cex=1.8)

#H 1 

plot(Reptiles$ACI_summer, Reptiles$Hill_1, xlab= "", pch=19, cex.axis=1.8,ylim=c(1,9),cex.lab=2, ylab=NA)
title(xlab = expression("ACI_summer Avg."),ylab = expression("Shannon's div."), line = 3.5,cex.lab=2) 
mod<-lm(Hill_1~ACI_summer, data=Reptiles)
var<-Reptiles$ACI_summer
newx <- seq(min(var), max(var), length.out=100)
preds <- as.data.frame(predict(mod, newdata = data.frame(ACI_summer=newx), interval = 'confidence'))

#add dashed lines for confidence bands
lines(newx, preds[ ,1], lty = 1,lwd=5, col = 'blue')
lines(newx, preds[ ,3], lty = 'dashed', col = 'blue')
lines(newx, preds[ ,2], lty = 'dashed', col = 'blue')
points(Reptiles$ACI_summer, Reptiles$Hill_1,  pch=19, cex=1.8)

#H 2

plot(Reptiles$M_evening, Reptiles$Hill_2, xlab= "", pch=19, cex.axis=1.8,ylim=c(1,7),cex.lab=2, ylab=NA)
title(xlab = expression("M evening Avg."), line = 3.5,cex.lab=2) 
mod<-lm(Hill_2~M_evening, data=Reptiles)
var<-Reptiles$M_evening
newx <- seq(min(var), max(var), length.out=100)
preds <- as.data.frame(predict(mod, newdata = data.frame(M_evening=newx), interval = 'confidence'))

#add dashed lines for confidence bands
lines(newx, preds[ ,1], lty = 1,lwd=3, col = 'blue')
lines(newx, preds[ ,3], lty = 'dashed', col = 'blue')
lines(newx, preds[ ,2], lty = 'dashed', col = 'blue')
points(Reptiles$M_evening, Reptiles$Hill_2,  pch=19, cex=1.8)

###################### Mammals
#H 0
plot(Mammals$BIO_spring , Mammals$Hill_0, xlab= "", pch=19, cex.axis=1.8,ylim=c(6,12),cex.lab=2, ylab=NA)
title(xlab = expression("BIO Spring Avg."),ylab = expression("Richness"), line = 3.5,cex.lab=2) 
mod<-lm(Hill_0~BIO_spring , data=Mammals)
var<-Mammals$BIO_spring 
newx <- seq(min(var), max(var), length.out=100)
preds <- as.data.frame(predict(mod, newdata = data.frame(BIO_spring =newx), interval = 'confidence'))

#add dashed lines for confidence bands
lines(newx, preds[ ,1], lty = 1,lwd=5, col = 'blue')
lines(newx, preds[ ,3], lty = 'dashed', col = 'blue')
lines(newx, preds[ ,2], lty = 'dashed', col = 'blue')
points(Mammals$BIO_spring , Mammals$Hill_0,  pch=19, cex=1.8)

#H 1 

plot(Mammals$BIO_spring, Mammals$Hill_1, xlab= "", pch=19, cex.axis=1.8,ylim=c(1,9),cex.lab=2, ylab=NA)
title(xlab = expression("BIO Spring Avg."),ylab = expression("Shannon's div."), line = 3.5,cex.lab=2) 
mod<-lm(Hill_1~BIO_spring, data=Mammals)
var<-Mammals$BIO_spring
newx <- seq(min(var), max(var), length.out=100)
preds <- as.data.frame(predict(mod, newdata = data.frame(BIO_spring=newx), interval = 'confidence'))

#add dashed lines for confidence bands
lines(newx, preds[ ,1], lty = 1,lwd=5, col = 'blue')
lines(newx, preds[ ,3], lty = 'dashed', col = 'blue')
lines(newx, preds[ ,2], lty = 'dashed', col = 'blue')
points(Mammals$BIO_spring, Mammals$Hill_1,  pch=19, cex=1.8)

#H 2

plot(Mammals$H_morning_90p, Mammals$Hill_2, xlab= "", pch=19, cex.axis=1.8,ylim=c(1,7),cex.lab=2, ylab=NA)
title(xlab = expression("H morning p"["90th"]), line = 3.5,cex.lab=2) 
mod<-lm(Hill_2~H_morning_90p, data=Mammals)
var<-Mammals$H_morning_90p
newx <- seq(min(var), max(var), length.out=100)
preds <- as.data.frame(predict(mod, newdata = data.frame(H_morning_90p=newx), interval = 'confidence'))

#add dashed lines for confidence bands
lines(newx, preds[ ,1], lty = 1,lwd=3, col = 'blue')
lines(newx, preds[ ,3], lty = 'dashed', col = 'blue')
lines(newx, preds[ ,2], lty = 'dashed', col = 'blue')
points(Mammals$H_morning_90p, Mammals$Hill_2,  pch=19, cex=1.8)



mtext(side=3,"Hill 0",2,adj=0.15, outer = TRUE, cex = 1.8)
mtext(side=3,"Hill 1",2,adj=0.52, outer = TRUE, cex = 1.8)
mtext(side=3,"Hill 2",2,adj=0.88, outer = TRUE, cex = 1.8)
mtext(side=2,"Birds",1, adj=1,outer = TRUE, cex = 1.8)
mtext(side=2,"Amphibians",1, adj=0.7,outer = TRUE, cex = 1.8)
mtext(side=2,"Reptiles",1, adj=0.4,outer = TRUE, cex = 1.8)
mtext(side=2,"Mammals",1, adj=0.1,outer = TRUE, cex = 1.8)
mtext(side=2,"Effective number of species",5, adj=0.8,outer = TRUE, cex = 2)


