
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

# Filflapoly <- readShapePoly("Filfla.shp")
# plot(Filflapoly) #this includes the plateau

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
time<-c(0,6) ## specify time in units over which we want to know survival (annual)
#time19 <- c(1, 20, 26, 28, 30, 33, 41, 42, 46, 47, 54, 64)
#time13 <- time <- c(1, 4, 5, 8, 18, 19, 20, 21, 26, 28, 30, 43, 51, 61, 71) 
#time19 <- c(1, 20, 26, 28, 30, 33, 41, 42, 46, 47, 54, 64)
#time <- c(0, 3, 4, 7, 17, 18, 19, 20, 25, 27, 29, 42, 50, 60, 70, 2114, 2134, 2140, 2142, 2144, 2147, 2155, 2156, 2160, 2161, 2168, 2178)



### CREATE SCR DATA OBJECT
scrdat <- ScrData$new(ESP_CH_comb, mesh, primary = primary, time = time) 
scrdat # plot & summary

scrdat$n() #unique individuals seen 

scrdat$n_occasions() # get number of occasions in the survey

scrdat$area() # get total area of the mesh in square kilometres 

scrdat$time()



### ADD EFFORT COVARIATE ###
effort<-fread("ESP_trap_locations_secr_robust.txt")
scrdat$add_covariate("effort", t(as.matrix(effort[,4:30])), "kj") ## transose matrix as occasions need to be in rows and detectors in columns



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

# fit model 
ESP_JS_RD_model$fit()

# see results 
ESP_JS_RD_model

ESP_JS_RD_model$get_par("lambda0", k = 1) # capture probability at each trap
ESP_JS_RD_model$get_par("sigma", k = 1) # home range radius/scale parameter around each trap
ESP_JS_RD_model$get_par("phi", k = 1) # annual survival from o
ESP_JS_RD_model$get_par("beta", k = 1)
ESP_JS_RD_model$get_par("D")

AIC(ESP_JS_RD_model) #131844.4

########### EXTRACT OUTPUT FROM MODEL ####
as.list(ESP_JS_RD_model)

## for occasion 1 = 2013
(1198000/100)*3.92
(1133000/100)*3.92
(1268000/100)*3.92

## for occasion 1 = 2019
(1509000/100)*3.92
(1420000/100)*3.92
(1604000/100)*3.92






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

