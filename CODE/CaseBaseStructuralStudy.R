##################################################################################################
# Script for Statistical Analysis and Visualization of Neuroimaging Data
# Author: Claudio Aleman Morillo    email: caleman@us.es
# Date: July 2024
# Version: 1.0

# Introduction:
#This script performs a linear regression looking for the relationship between
#the baseline centiles of the cohort and covariates such as age, sex or medication.


#Script inputs:
#  1- File with the cohort data (age, sex, medication, regional centiles)
#Script outputs:
#  1- Dataframe that contains all the t-values of each covariate for each cortical region and if this value is significant (p-value<0.05). 
#  2- .png image files that represent the t values previously found on an image of the cerebral cortex.
###################################################################################################
#functions:

# 1)function regional_brainmap_representation_borders
#Function to represent the data obtained.
#Input Parameters:
# 1- vector of the 34 cortical t-values.
# 2- image title
# 3-maximum value of the representable t-value scale
# 4-minimum value of the representable t-value scale
# 5-mean value of the representable t-value scale

#Output:
#  1- Image


# 2) Function that performs linear interpolation using the data from the indicated column as the 
#    dependent variable (Y).
# Input Parameters:
#  1- dataframe with all the data


#output:
# Dataframe:
# 1-dataframe containing the p-values of each covariate for each regional cortical interpolation.
###################################################################################################


# Clear the workspace
rm(list=ls())


# Load necessary libraries
library(dplyr)
library(lmerTest)        
library(ggplot2)
library(ggseg)

# Function to perform generalized linear model analysis
LocalClinicalAnalisis_GeneralizedLinearModelGauss_V30<-function(D){
  
  names_Objetives <- c( 'bankssts',
                        'caudalanteriorcingulate',
                        'caudalmiddlefrontal',
                        'cuneus',
                        'entorhinal',
                        'fusiform',
                        'inferiorparietal',
                        'inferiortemporal',
                        'lateraloccipital',
                        'lateralorbitofrontal',
                        'lingual',
                        'medialorbitofrontal',
                        'middletemporal',
                        'paracentral',
                        'parsorbitalis',
                        'parstriangularis',
                        'pericalcarine',
                        'postcentral',
                        'posteriorcingulate',
                        'precentral',
                        'precuneus',
                        'rostralanteriorcingulate',
                        'rostralmiddlefrontal',
                        'superiorfrontal',
                        'superiorparietal',
                        'superiortemporal',
                        'supramarginal',
                        'temporalpole',
                        'insula',
                        'isthmuscingulate',
                        'parahippocampal',
                        'parsopercularis',
                        'frontalpole',
                        'transversetemporal')
  
  
  v<-c('dx','Age_inclusion','Sex','residual_etiv')
  
  res <- data.frame(matrix(0,ncol=length(names_Objetives),nrow=length(v)))
  names(res)<-names_Objetives
  row.names(res)<-v
  
  resPValues <- data.frame(matrix(0,ncol=length(names_Objetives),nrow=length(v)))
  names(resPValues)<-names_Objetives
  row.names(resPValues)<-v
  
  n<-length(v)+1
  print('nnn')
  print(n)
  print('----')
  
  
  
  
  for(i in 1:length(names_Objetives)){
    
    formulaText<-paste0(names_Objetives[[i]],'~1+dx+Age_inclusion+Sex+residual_etiv',sep='')
    
    print(formulaText)
    linealModel <- lm(formulaText,data=D)
    print(summary(linealModel))
    
    
    res[,i]<-(as.numeric((coef(summary(linealModel)))[2:n,3]))
    resPValues[,i]<-(as.numeric((coef(summary(linealModel)))[2:n,4]))
    
  }
  
  
  
  
  
  

  res<-t(res)
  resPValues<-t(resPValues)
  print(resPValues)

  for(z in 1:ncol(res)){
    resPValues[,z]<- p.adjust(resPValues[,z], method="fdr")
  }


  
  
  for(i in 1:nrow(res)){
    
    for(j in 1:ncol(res)){
      if(resPValues[i,j]>0.05){
        #res[i,j]<-NaN
        resPValues[i,j]=0
        
      }else{
        resPValues[i,j]=1
      }
      
    }
    
  }
  
  resD<-data.frame(res)
  respD<-data.frame(resPValues)
  colnames(respD)<-paste0(colnames(respD),'_pValues')
 
  resD$namesArea<-rownames(resD)
  respD$namesArea<-rownames(respD)
  
  print(respD)
  print(resD)
  return(merge(resD,respD,by='namesArea'))
  
  
  
  
}

# Function to create a brain map visualization
regional_brainmap_representation_borders <- function(data,title,sup_lim,inf_lim,midd_p){
  
  
  
  zones <- c( 'bankssts',
              'caudal anterior cingulate',
              'caudal middle frontal',
              'cuneus',
              'entorhinal',
              'frontal pole',
              'fusiform',
              'inferior parietal',
              'inferior temporal',
              'insula',
              'isthmus cingulate',
              'lateral occipital',
              'lateral orbitofrontal',
              'lingual',
              'medial orbitofrontal',
              'middle temporal',
              'paracentral',
              'parahippocampal',
              'pars opercularis',
              'pars orbitalis',
              'pars triangularis',
              'pericalcarine',
              'postcentral',
              'posterior cingulate',
              'precentral',
              'precuneus',
              'rostral anterior cingulate',
              'rostral middle frontal',
              'superior frontal',
              'superior parietal',
              'superior temporal',
              'supramarginal',
              'temporal pole',
              'transverse temporal')
  
  if(nrow(data) == 34){
    someData <- tibble(
      region = zones, 
      values = data$values,
      significant = data$sig,
      groups = c(rep(title, nrow(data)))
    )
    
  }else if (nrow(data) == 68){
    someData <- tibble(
      label = c(paste0('lh_',gsub(' ','',zones)),paste0('rh_',gsub(' ','',zones))),
      values = data$values,
      significant = data$sig,
      groups = c(rep(title, nrow(data)))
    )
  }
  
  if (any(someData$significant) && any(!someData$significant)){
    borders <- c(rgb(0, 0, 0, alpha = 0.1), "black")
  }else if (all(someData$significant)){
    borders <- rgb(0, 0, 0, alpha = 1)
  }else{
    borders <- rgb(0, 0, 0, alpha = 0.1)
  }
  
  
  # Brain map configuration
  someData %>%
    group_by(groups) %>%
    ggplot() + scale_fill_gradient2(
      
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = midd_p,
      limits=c(inf_lim,sup_lim),
      space = "Lab",
      na.value = "grey50",
     
      guide = guide_colorbar(title = "values",
                             
                             direction = 'vertical',
                             barheight = 8,
                             label.theme = element_text(color = "black"),
                             title.theme = element_text(color = "black")),
      
      aesthetics = "fill"
    ) +
    geom_brain(
      mapping = aes(fill = values, col = significant,alpha = significant, size = significant),
      atlas = dk, 
      position = position_brain(hemi ~ side),
      show.legend = TRUE) +
    scale_colour_manual(values=borders) +
   
    scale_alpha_manual(values= c(0.6,1)) +                   
    scale_size_manual(values= c(0.9,1)) + 
    facet_wrap(~groups)+theme_brain2(plot.background = "white")
  
}





# Load and preprocess data
# replace  the  'path_to_datafile' path with the path to centiles data file
DataBaseGlobal<-'path_to_datafile'
centils<-read.csv(DataBaseGlobal)
centils$CPZ_equivalent<-ifelse(centils$Assessment==1,0,centils$CPZ_equivalent)
centils$dx<-ifelse(centils$Subject>=1000,0,1)
centils$CPZ_equivalent<-ifelse(centils$dx==0,0.0,centils$CPZ_equivalent)
centils$dx<-as.factor(centils$dx)
centils$Age_inclusion<-scale(centils$Age_inclusion)
centils$Treatment_Time<-scale(centils$Treatment_Time)
centils$CPZ_equivalent<-scale(centils$CPZ_equivalent)
centils$residual_etiv<-scale(centils$residual_etiv)
centils$Sex<-as.factor(centils$Sex)
centilsBaseCase<-filter(centils,Assessment%in%c(1))

# Perform analysis
clinicalDataSetDataInterpolationResults<-LocalClinicalAnalisis_GeneralizedLinearModelGauss_V30(centilsBaseCase)

print(clinicalDataSetDataInterpolationResults)


# Generate and save brain maps
var_name<-' Base Case Analisys '
for (i in 2:5){
  varNameCol<-colnames(clinicalDataSetDataInterpolationResults)[i]
  varNameColpValue<-paste0(colnames(clinicalDataSetDataInterpolationResults)[i],'_pValues')
  dataRepresent<-data.frame(clinicalDataSetDataInterpolationResults[,c(varNameCol,varNameColpValue)])
  
  colnames(dataRepresent)<-c('values','sig')
  dataRepresent[,'sig']<-as.logical(dataRepresent[,'sig'])
  regional_brainmap_representation_borders(dataRepresent,paste0(var_name,'-> Fixed effect: ',varNameCol,sep=''),10,-10,0)
  #replace  the  'path_to_file' path with the path where you want to save the images
  ggsave(paste0('path_to_file/BaseCaseStructural_TValues_',varNameCol,'_.png',sep=''),dpi=400,height=5,width=8)
  
}












