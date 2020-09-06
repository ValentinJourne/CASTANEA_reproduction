#Sys.setenv(JAVA_HOME = '/usr/lib/jvm/jdk1.8.0_251')

###########################
#load packages 
library(chron)
library(ggplot2)
library(cowplot)
library(hydroGOF)
library(stringr)
library(tidyr)
library(dplyr)
library(directlabels)
library(chron)
library(rJava)
library(readr)
library(reshape2)
library(tibble)
library(MASS)

#for parallel loop
library(foreach)
library(doParallel)
library(lubridate)
library(purrr)
###########################


############################
# define user path
user<- "valentin_unix_pers"
###########################

###########################
#specify path of capsis and inventory data 
###########################

if(user == "valentin_unix_pers"){
  localCapsis="/media/journe/DDlabo1/Applications/capsis4"
  pathToParametersFilesDirectory <- paste0(localCapsis,"/data/castaneaonly") #path for data init castanea
  capsisVarPath= paste0(localCapsis,"/var/") #var path where we found castanaea outputs simulations
  capsisRPath= paste0(localCapsis,"/src/castaneaonly/R/") #path for code translation java to R
  pathInventories= "/media/journe/DDlabo1/Thesis/Data/BDD/inventories/" #bdd for data 
  #specify path output 
  outResults <- paste("results_outputs_simAll")
  if(!file.exists(outResults))dir.create(outResults)
  pathSim = outResults
}

############################
#path for inventories for init.
###########################
DryadeInv2013 = read.table(paste0(pathInventories,"DryadeInv2013.csv"),header=TRUE,sep=";",dec=".") 
DryadeInv2007_TCN = read.table(paste0(pathInventories,"DryadeInv2007_TCN.csv"),header=TRUE,sep=";",dec=".") 

###########################
#create inventories for Castanea run 
###########################

inventories = subset(DryadeInv2007_TCN,!is.na(DryadeInv2007_TCN$c130_mm))
inventories$dbh = inventories$c130_mm/pi/10
aGF = 1.27
bGF = 0.745 
woodDensity = 450
inventories$height = aGF*inventories$dbh^bGF
inventories$treeVolume = 0.486*(inventories$dbh/100/2)^2*pi*inventories$height
sumStand = aggregate(inventories$treeVolume,by=list(inventories$placette), FUN="sum")
meanStand = aggregate(inventories,by=list(inventories$placette), FUN="mean")
lengthStand = aggregate(inventories$placette,by=list(inventories$placette), FUN="length")
Vha = sumStand$x/400*10000
Nha = lengthStand$x/400*10000
dbh = meanStand$dbh
mean(Nha)

#no sylvicularal scenario 
scenarioFileName = "no"

#specify if outputs variables neeeded
pathToCapsisDirectory <- localCapsis
dailyVariablesNames <-  NULL
yearlyVariablesNames <-  NULL


#use template old inventory
inventoryFileNameIn = paste0(pathToParametersFilesDirectory,"/inventories/reproduction/Old_inventory/Inventory_Repro_Abies_TC.txt")


######################################
#test submodel of reproduction based on sigmoid
######################################

inventoryFileNameOut = paste0(pathToParametersFilesDirectory,"/inventories/reproduction/","Inventory_Repro_Abies_TC_SIGMOID.txt")
inventory = readLines(inventoryFileNameIn)
base = inventory[60]
firstPart = str_sub(base,1,5)
secondPart = str_sub(base,9,106)
thirdPart = str_sub(base,123,156)
lambda = 0.02
tronviv = 0.1
rateOfSeedJourne = 0.5
reservesToReproduce = 200
seedMass_mu = 2.11
seedMass_spe = 50.5
seedMortality =0.0032

for( k in c(1:16)){
  inventory[59+k]= paste0(firstPart,as.character(k),"\t",as.character(1),secondPart,as.character(dbh[k]),"\t",as.character(Nha[k]),"\t",as.character(Vha[k]), thirdPart, as.character(rateOfSeedJourne)," ",as.character(reservesToReproduce)," ",as.character(seedMass_mu)," ",as.character(seedMass_spe)," ",as.character(seedMortality)," ", as.character(lambda)," ",as.character(tronviv))
}

optionalVariables = "optionalVariables= LMA nitrogen aGF bGF rateOfSeedJourne reservesToReproduce seedMass_mu seedMass_spe seedMortality lambdaSeed tronviv"

inventory[5] = "year = 1959"
inventory[21] = "typeOfVegetation= TYPE_VEG_STAND"
inventory[24] = "initWoodByVolume = TRUE"
#inventory[27]="fit2018FileName = fit2018/UniChillRenecofor.fit2018"
inventory[27] = "fit2018FileName = fit2018/UniChillThreshold_repro.fit2018" #update renefocoer analysis on Mt Ventoux
inventory[40] = "initRepro= REPRO_INIT_STORAGE_SIGMOID"
inventory[47] = "parameterPot = -10.67"
inventory[51] = "logPrefix= TC_sigmoid_"       
inventory[55] = optionalVariables

#compile a text file 
writeLines(inventory,inventoryFileNameOut)

# Loads rJava package

source(paste0(capsisRPath,"castaneaOne.R"))

#specify initiation and period for simulation 
fromYear <- 1959 #2007
enYear <- 2015
numberOfYears <- enYear-fromYear+1

flistToRemove <- list.files(capsisVarPath, pattern="TC_")
file.remove(paste0(capsisVarPath,flistToRemove))

#specify climate file and species file 
speciesFileName <- paste0(pathToParametersFilesDirectory,"/species/CastaneaSpecies4.txt")
climateFileName <- paste0(pathToParametersFilesDirectory,"/climate/Dvx3_59to2015.txt")
pathToParametersFilesDirectory <- paste0(localCapsis,"/data/castaneaonly") 
inventoryFileName <- paste0(pathToParametersFilesDirectory,"/inventories/reproduction/Inventory_Repro_Abies_TC_SIGMOID.txt")

#test run part 
#run simulation 
simTC = castaneaOne(pathToCapsisDirectory,inventoryFileNameOut,speciesFileName,climateFileName,scenarioFileName,yearlyVariablesNames,dailyVariablesNames,fromYear,numberOfYears)

flistToRemove <- list.files(capsisVarPath, pattern="lck")
file.remove(paste0(capsisVarPath,flistToRemove))
flist <- list.files(capsisVarPath, pattern="TC_")
file.copy(paste0(capsisVarPath,flist), pathSim, overwrite=T)

yealySimoid = read.table(paste0(pathSim,"TC_sigmoid_yearlyResults.log"), sep=";", head=T,dec=".")


######################################
#test submodel of reproduction based on threshold 
######################################

inventoryFileNameOut = paste0(pathToParametersFilesDirectory,"/inventories/reproduction/","Inventory_Repro_Abies_TC_THRESHOLD.txt")
inventory = readLines(inventoryFileNameIn)
base = inventory[60]
firstPart = str_sub(base,1,5)
secondPart = str_sub(base,9,106)
thirdPart = str_sub(base,123,nchar(base))

lambda = 0.02
tronviv = 0.1

for( k in c(1:16)){
  inventory[59+k] = paste0(firstPart,as.character(k),"\t",as.character(1),secondPart,as.character(dbh[k]),"\t",as.character(Nha[k]),"\t",as.character(Vha[k]),thirdPart," ",as.character(lambda)," ",as.character(tronviv))
  
}

optionalVariables = "optionalVariables= LMA nitrogen aGF bGF rateOfSeedJourne reservesToReproduce seedMass_mu seedMass_spe seedMortality lambdaSeed tronviv"

inventory[5] = "year = 2007"
inventory[21] = "typeOfVegetation= TYPE_VEG_STAND"
inventory[24] = "initWoodByVolume = TRUE"
inventory[27] = "fit2018FileName = fit2018/UniChillRenecofor.fit2018"
inventory[40] = "initRepro= REPRO_INIT_STORAGE_THRESHOLD"
inventory[47] = "parameterPot = -10.67"
inventory[51] = "logPrefix= TC_threshold_"       
inventory[55] = optionalVariables

writeLines(inventory,inventoryFileNameOut)

simTC = castaneaOne(pathToCapsisDirectory,inventoryFileNameOut,speciesFileName,climateFileName,scenarioFileName,yearlyVariablesNames,dailyVariablesNames,fromYear,numberOfYears)

flistToRemove <- list.files(capsisVarPath, pattern="lck")
file.remove(paste0(capsisVarPath,flistToRemove))
flist <- list.files(capsisVarPath, pattern="TC_")


file.copy(paste0(capsisVarPath,flist), pathSim, overwrite=T)

yealyThreshold = read.table(paste0(pathSim,"TC_threshold_yearlyResults.log"), sep=";", head=T,dec=".")
