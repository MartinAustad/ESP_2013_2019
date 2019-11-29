####### DATA PREPARATION FROM ACCESS DATABASE ################################################################

####### for CMR analysis of Filfla ESP 2013 & 2019 data using Openpopscr package #############################





library(RODBC)

library(data.table)

library(tidyverse)

library(sf)



#Access database with data from both years 



setwd("C:\\Users\\Martin\\Documents\\working_folder\\CaptureMarkRecapture\\stormies\\ESP_2019_analysis")

setwd("C:\\STEFFEN\\RSPB\\Malta\\Raw_data")



SP<-odbcConnectAccess2007('ESP_Malta_2013_2019.accdb')

ESP_comb <- sqlQuery(SP, "SELECT * FROM SECR_captures")  

trap_comb <- sqlQuery(SP, "SELECT * FROM SECR_input_detectors")

effort <- sqlQuery(SP, "SELECT * FROM SECR_input_effort")

odbcClose(SP)





#########################################################################

# ENSURE THAT RING REPLACEMENTS ARE INCORPORATED INTO RECORDS

# adapted from Steffen's YESH script

#########################################################################

setwd("C:\\STEFFEN\\RSPB\\Malta\\Analysis\\ESP_abundance_survival\\ESP_2013_2019")

rings <- fread("ESP_ring_replacment.csv")



## CREATE REPLACEMENT LIST

replist<-rings %>% #dplyr::filter(Replacement_ring!="") %>%
  
  rename(orig=Origninal_ring, repl=Replacement_ring) %>%
  
  select(orig,repl)



## UPDATE RECORDS

ESP<- ESP_comb %>%
  
  mutate(Ring_Nr=ifelse(Ring_Nr %in% replist$repl,as.character(replist$orig[match(Ring_Nr,replist$repl)]),as.character(Ring_Nr)))



#remove duplicate values manually

ESP <- ESP[c(-7678),]

####filter out chicks in 2013 and then remove age vector

ESP <- ESP %>%
  filter(Age!=1) #9 records of chicks removed (euring code age 1) #7 ind - two of which retrapped as adults in 2019

ESP <- ESP[,-5]


#save for secr capthis

write.table(ESP, "ESP_secr_input_robust.txt", row.names=F, col.names=F, sep=" ")





########preparation of trap#################



trap_sf<-st_sfc(st_multipoint(as.matrix(trap_comb[,3:4]), dim="XY"))

st_crs(trap_sf)=4326

trap_sf<-st_transform(trap_sf, crs=32633) 

trap_comb[,3:4]<-st_coordinates(trap_sf)[,1:2]

trap<-trap_comb[order(trap_comb$Net_ID, decreasing=F),c(1,3,4,5:31)]

trap[is.na(trap)]<-0

names(trap)[2:3]<-c('x','y')



#save for secr

write.table(trap, "ESP_trap_locations_secr_robust.txt", row.names=F, col.names=F, sep=" ", quote=F)

##########preperation of effort - for covariate#####



effortX <- effort[,-(2:3)]

effortX <- effortX %>%
  group_by(EncOcc) %>%
  summarise_all(funs(sum(., na.rm = TRUE)))

effortX <- effortX[,-1]

write.table(effortX, "ESP_effort_robust.txt", row.names=F, col.names=F, sep= " ", quote=F)

#################################################################

############CMR analysis of Filfla ESP 2013 & 2019 data using Openpopscr package###############################



library(R6)

library(openpopscr)

library(secr)

library(sp)

library(data.table)

library(maptools)

library(rgdal)

library(dplyr)

library(foreign)



#source of code for openpopscr package: 

#https://github.com/r-glennie/openpopscr/blob/master/vignettes/ScrData.Rmd

#https://github.com/r-glennie/openpopscr/blob/master/vignettes/ScrModel.Rmd

#https://github.com/r-glennie/openpopscr/blob/master/inst/examples/10_js_robust.R



#script from Steffen Oppel's SECR_density_estimation integrated as necessary

#compiled by Martin Austad Nov 2019

# modified by Steffen Oppel to run only robust design model





## MAIN QUESTIONS: how do you specify the detection function (or is a hazard rate fixed)?

## How to extract output? Why are density estimates reported on link scale and with get_par different?

## How to allow density to vary between two primary occasions?







################################################################################################################

###################### LOAD SECR INPUT  ######################################################

################################################################################################################



## call RODBC script in 32-bit R to extract data and write them into text file

## only need to re-run after changes to database

system(paste0(Sys.getenv("R_HOME"), "/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\Malta\\Analysis\\ESP_abundance_survival\\ESP_2013_2019\\RODBC_ESP.r")), wait = FALSE, invisible = FALSE)





## LOAD SPATIAL EXTENT OF ISLAND

setwd("C:\\Users\\Martin\\Documents\\stormiesCMR")

setwd("C:\\STEFFEN\\RSPB\\Malta\\Analysis\\ESP_abundance_survival\\ESP_2013_2019")



#Filflapoly <- readShapePoly("Filfla.shp")

#plot(Filflapoly) #this includes the plateau



Boulderscreepoly <- readShapePoly("boulderscree.shp")

plot(Boulderscreepoly)





### CREATE SCR DATA OBJECT

ESP_CH_comb <- read.capthist("ESP_secr_input_robust.txt",
                             
                             "ESP_trap_locations_secr_robust.txt",
                             
                             detector="count", verify=T, binary.usage = FALSE)



trap <- traps(ESP_CH_comb)

detectors <- read.traps(data = trap, detector = "count") 

mesh<-make.mask(detectors, type='trapbuffer', spacing=10, buffer=80, poly=Boulderscreepoly)



plot(ESP_CH_comb, tracks=T, varycol=FALSE)   # tracks=T specifies that you want lines to connect subsequent sightings of the same individual at different detectors

plot(mesh)

plot(ESP_CH_comb, tracks=T, varycol=FALSE, add=T)

plot(trap, add=T)

################################################################################################################

###################### SET UP TEMPORAL STRUCTURE OF MODEL  #####################################################

################################################################################################################



### PRIMARY OCCASIONS: 2 (2013 and 2019)

## vector of primary occasion number for each secondary occasion

primary <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2)





### Occasions may occur at irregular intervals. This information can be included in the ScrData object when it is created.

## FOR ROBUST DESIGN ONLY TIME FOR PRIMARY OCCASIONS

#time<-c(0,6) ## specify time in units over which we want to know survival (annual)

time <-c(2013, 2019)

### CREATE SCR DATA OBJECT

scrdat <- ScrData$new(ESP_CH_comb, mesh, primary = primary, time = time) 

scrdat # plot & summary



scrdat$n() #unique individuals seen 



scrdat$n_occasions() # get number of occasions in the survey



scrdat$area() # get total area of the mesh in square kilometres 



scrdat$time()







### ADD EFFORT COVARIATE ###

effort<-fread("ESP_effort_robust.txt")

#effort<-fread("ESP_trap_locations_secr_robust.txt")

#effort <- lapply(effort, as.numeric)

#effort <- as.data.frame(effort)

#scrdat$add_covariate("effort", t(as.matrix(effortY[,4:30])), "kj") ## transose matrix as occasions need to be in rows and detectors in columns

scrdat$add_covariate("effort", as.matrix(effort[,1:9]), "kj")


################################################################################################################

###################### FIT SIMPLE ROBUST DESIGN MODEL  #####################################################

################################################################################################################



####fit simple model withg intercept only

## specify formulas for each parameter

par <- list(lambda0 ~ effort,   ## this should have an effort correction (i.e. net metre hours)
            
            sigma ~ 1, 
            
            phi ~ 1, 
            
            beta ~ 1,
            
            D ~ 1)



start <- get_start_values(scrdat, model = "JsModel")

start



# create model object 

ESP_JS_RD_model <- JsModel$new(par, scrdat, start)

#with both effort versions:
#Warning message:
  #In (function (..., row.names = NULL, check.rows = FALSE, check.names = TRUE,  :
                 # row names were found from a short variable and have been discarded


# fit model 

ESP_JS_RD_model$fit()

#with effort in netlenght*hours #
#Warning message:
  #In ESP_JS_RD_model$fit() : model failed to converge with nlm code 3

#D estimate values smaller (still make sense) and #lamda0.effort CI does not pass zero
# see results 

ESP_JS_RD_model



ESP_JS_RD_model$get_par("lambda0", k = 1) # capture probability at each trap

ESP_JS_RD_model$get_par("sigma", k = 1) # home range radius/scale parameter around each trap

ESP_JS_RD_model$get_par("phi", k = 1) # annual survival from o

ESP_JS_RD_model$get_par("beta", k = 1)

ESP_JS_RD_model$get_par("D", k=1)



AIC(ESP_JS_RD_model) 



########### EXTRACT OUTPUT FROM MODEL ####

as.list(ESP_JS_RD_model)



## for occasion 1 = 2013

(1198000/100)*3.92

(1133000/100)*3.92

(1268000/100)*3.92



## for occasion 2 = 2019

(1509000/100)*3.92

(1420000/100)*3.92

(1604000/100)*3.92

## for occasion 1 = 2013
(1120000/100)*2.75

(1059000/100)*2.75

(1185000/100)*2.75

## for occasion 2 = 2019
(1406000/100)*2.75

(1322000/100)*2.75

(1494000/100)*2.75









################################################################################################################

###################### FIT TRANSIENT MODEL  #####################################################

################################################################################################################



form <- list(lambda0 ~ effort, 
             
             sigma ~ 1, 
             
             beta ~ 1, 
             
             phi ~ 1, 
             
             sd ~ 1, 
             
             D ~ 1)





start <- get_start_values(scrdat, model = "JsTransientModel")

start



ESP_JS_transient_model  <- JsTransientModel$new(form, scrdat, start)



ESP_JS_transient_model$fit()



ESP_JS_transient_model



# look at parameters on response scale 

ESP_JS_transient_model$get_par("lambda0", k = 1)

ESP_JS_transient_model$get_par("phi", k = 1)

ESP_JS_transient_model$get_par("D")



AIC(ESP_JS_transient_model)

