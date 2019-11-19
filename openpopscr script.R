
############CMR analysis of Filfla ESP 2013 & 2019 data using Openpopscr package###############################


library(openpopscr)
library(secr)
library(RODBC)
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

################################################################################################################
###################### LOAD SECR INPUT  ######################################################
################################################################################################################
setwd("C:\\Users\\Martin\\Documents\\stormiesCMR")

Filflapoly <- readShapePoly("Filfla.shp")
plot(Filflapoly) #this includes the plateau

ESP <- fread("ESP_secr_input.txt") #generated previously from Steffen's SECR_density_estimation script


#########################################################################
####ENSURE THAT RING REPLACEMENTS ARE INCORPORATED INTO RECORDS
# adapted from Steffen's YESH script
#########################################################################

length(unique(ESP$V2))
rings <- fread("ESP_ring_replacment.csv")


## CREATE REPLACEMENT LIST

replist<-rings %>% #dplyr::filter(Replacement_ring!="") %>%
  
  rename(orig=Origninal_ring, repl=Replacement_ring) %>%
  
  select(orig,repl)

## UPDATE RECORDS

ESP<- ESP %>%
  
  mutate(V2=ifelse(V2 %in% replist$repl,replist$orig[match(V2,replist$repl)],V2))

length(unique(ESP$V2)) 

View(ESP)
#remove duplicate values manually #found by filter
ESP <- ESP[c(-3534),]

#save for secr capthis
write.table(ESP, "ESP_secr_input_rr.txt", row.names=F, col.names=F, sep=" ")


###############Preparing data#########################################
#count detector: means that each detector records how many times each individual was seen by that detector
#proximity these detectors only record whether an individual was seen at least once or not. 
#Multi-catch detectors can also be used where individuals may only be detected on a single detector for each occasion, but a detector can detect multiple individuals.

ESP_CH <- read.capthist("ESP_secr_input_rr.txt", "ESP_trap_locations_secr.txt", detector="count", verify=T, binary.usage = FALSE)

trap <- traps(ESP_CH)
detectors <- read.traps(data = trap, detector = "count") 
#mesh <- make.mask(detectors)
mesh<-make.mask(detectors, type='trapbuffer', spacing=10, buffer=80, poly=Filflapoly)


plot(mesh)
plot(ESP_CH, tracks=T, varycol=FALSE, add=T)   # tracks=T specifies that you want lines to connect subsequent sightings of the same individual at different detectors
plot(trap, add=T)

#Occasions may occur at irregular intervals. This information can be included in the ScrData object when it is created.
time <- c(1, 20, 26, 28, 30, 33, 41, 42, 46, 47, 54, 64)

  
scrdat <- ScrData$new(ESP_CH, mesh, time = time)
scrdat # plot & summary

scrdat$n() #unique individuals seen 

scrdat$n_occasions() # get number of occasions in the survey

scrdat$area() # get total area of the mesh in square kilometres 

scrdat$time()

################# Fitting a model #####################################
#####ScrModel: Spatial Capture-Recapture Model Object##################


# set each parameter to be a constant 
form <- list(lambda0 ~ 1, 
             sigma  ~ 1, 
             D ~ 1)
# get some starting values based on data 
start <- get_start_values(scrdat)
start 

mod <- ScrModel$new(form, scrdat, start) #create model object

mod$fit()
mod$par()

# Reading model output
#mod
#The package <code>openpopscr</code>
#uses the log-link function for all of these parameters. So, you would take
#the exponential of these reported numbers to obtain the parameter values on the
#response scale. 
mod$get_par("lambda0", k = 1, j = 1)
mod$get_par("sigma", k = 1, j = 1)
mod$get_par("D") #1799109/km2 

mod$estimates()
AIC(mod) #44689.12

####for 2013
setwd("C:\\Users\\Martin\\Documents\\stormiesCMR\\ESP_2019_analysis\\2013")

ESP_CH13 <- read.capthist("ESP_secr_input.txt", "ESP_trap_locations_secr.txt", detector="count", verify=T, binary.usage = FALSE)

trap <- traps(ESP_CH13)
detectors13 <- read.traps(data = trap, detector = "count")

#Occasions may occur at irregular intervals. This information can be included in the ScrData object when it is created.
time13 <- c(1, 4, 5, 8, 18, 19, 20, 21, 26, 28, 30, 43, 51, 61, 71) 

scrdat13 <- ScrData$new(ESP_CH13, mesh, time = time13)
scrdat13 # plot & summary

scrdat13$n() #unique individuals seen 

scrdat13$n_occasions() # get number of occasions in the survey

scrdat13$area() # get total area of the mesh in square kilometres 

scrdat13$time()

#2013####ScrModel: Spatial Capture-Recapture Model Object##################
# set number of occasions to simulate
n_occasions <- 15

# set each parameter to be a constant 
form <- list(lambda0 ~ 1, 
             sigma  ~ 1, 
             D ~ 1)
# get some starting values based on data 
start13 <- get_start_values(scrdat13)
start13

mod13 <- ScrModel$new(form, scrdat13, start13) #create model object

mod13$fit()
mod13$par()
mod13$get_par("lambda0", k = 1, j = 1)
mod13$get_par("sigma", k = 1, j = 1)
mod13$get_par("D") #748990.9/km2 ie 7489.909/ha

mod13$estimates()


###########################################################################
#####################"JsModel: Jolly-SeberSCR Model Object"################
###########################################################################

#tried only for 2019

# set true parameters #not sure what this means 
#true_par <- list(lambda0 = 0.2, sigma = 20, phi = 0.8, beta = 0.2, D = 1000)

n_occasions <- 12

#simulatescr
#scrdat <- simulate_js_openscr(scrdat, n_occasions, detectors, mask) #seed = 52381

form <- list(lambda0 ~ 1, 
             sigma  ~ 1,
             phi ~ 1, 
             beta ~ 1,
             D ~ 1) # D inserted by MA cause of this error Error in names(private$form_) <- c(private$detfn_$pars(), "phi", "beta",  :'names' attribute [5] must be the same length as the vector [4]

start <- get_start_values(scrdat, model = "JsModel")
start
mod <- JsModel$new(form, scrdat, start)

mod$fit()
#Completed model fitting in 15.558 mins 
#Checking convergence.......converged 
#Computing variances.......done
#Inferring density..........done
#Computing variances..........done
#Computing confidence intervals..........done
#Warning messages:
 # 1: In self$set_mle(mle, V, llk) :
  #Variance estimates not reliable, do a bootstrap.
#2: In sqrt(diag(private$V_)) : NaNs produced
#3: In sqrt(diag(V)) : NaNs produced
#4: In sqrt(diag(V)) : NaNs produced

mod


mod$get_par("lambda0", k = 1, j = 1)
mod$get_par("sigma", k = 1, j = 1)
mod$get_par("phi", k = 1)
mod$get_par("beta", k = 1)
mod$get_par("D") #2477978/km2 ie 24779.78/ha

mod$estimates()
AIC(mod)

###########################2013 & 2019 combined: Robust design #########################


############################Preperation of data in SECR format ##########

#Access database with data from both years 

setwd("C:\\Users\\Martin\\Documents\\working_folder\\CaptureMarkRecapture\\stormies")

SP<-odbcConnectAccess2007('ESP_Malta_2013_2019.accdb')
ESP_comb <- sqlQuery(SP, "SELECT * FROM SECR_captures")  
trap_comb <- sqlQuery(SP, "SELECT * FROM SECR_input_detectors")
odbcClose(SP)

fwrite(ESP_comb, "ESP_Secr_20132019.csv") #can skip but no Microsoft tools on laptop i ran openpopscr
fwrite(trap_comb, "Trap_20132019.csv")  #an skip but no Microsoft tools on laptop i ran openpopscr 

#no Access on this laptop 
setwd("C:\\Users\\Martin\\Documents\\stormiesCMR")

ESP1319 <- fread("ESP_Secr_20132019.csv")

ESP <- ESP1319
length(unique(ESP$Ring_Nr))

trap1319 <-fread("Trap_20132019.csv")
#trap <- trap1319

#########################################################################

# ENSURE THAT RING REPLACEMENTS ARE INCORPORATED INTO RECORDS
# adapted from Steffen's YESH script
#########################################################################

rings <- fread("ESP_ring_replacment.csv")


## CREATE REPLACEMENT LIST

replist<-rings %>% #dplyr::filter(Replacement_ring!="") %>%
  
  rename(orig=Origninal_ring, repl=Replacement_ring) %>%
  
  select(orig,repl)

## UPDATE RECORDS

ESP<- ESP %>%
  
  mutate(Ring_Nr=ifelse(Ring_Nr %in% replist$repl,replist$orig[match(Ring_Nr,replist$repl)],Ring_Nr))

length(unique(ESP$Ring_Nr)) 

View(ESP)
#remove duplicate values manually
ESP <- ESP[c(-7678),]

#save for secr capthis
write.table(ESP, "ESP_secr_input_robust.txt", row.names=F, col.names=F, sep=" ")

########preperation of trap#################

trap[is.na(trap)]<-0

### CONVERT COORDINATES INTO UTM COORDINATES
LongLatToUTM<-function(x,y,zone){
  xy <- data.frame(ID = 1:length(x), X = x, Y = y)
  coordinates(xy) <- c("X", "Y")
  proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
  res <- spTransform(xy, CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep='')))
  return(as.data.frame(res))
}

trap <- as.data.frame(trap)
x<-c( trap$LONG)
y<-c( trap$LAT)
trap[,32:34] <- LongLatToUTM(x,y,33)

trap<-trap[order(trap$Net_ID, decreasing=F),c(1,33,34,5:31)]

#save for secr
write.table(trap, "ESP_trap_locations_secr_robust.txt", row.names=F, col.names=F, sep=" ", quote=F)
#################################################################


ESP_CH_comb <- read.capthist("ESP_secr_input_robust.txt", "ESP_trap_locations_secr_robust.txt", detector="count", verify=T, binary.usage = FALSE)

trap <- traps(ESP_CH_comb)
detectors <- read.traps(data = trap, detector = "count") 
#mesh <- make.mask(detectors)
mesh<-make.mask(detectors, type='trapbuffer', spacing=10, buffer=80, poly=Filflapoly)

plot(ESP_CH_comb, tracks=T, varycol=FALSE)   # tracks=T specifies that you want lines to connect subsequent sightings of the same individual at different detectors
plot(mesh)
plot(ESP_CH_comb, tracks=T, varycol=FALSE, add=T)
plot(trap, add=T)

#Occasions may occur at irregular intervals. This information can be included in the ScrData object when it is created.
#time19 <- c(1, 20, 26, 28, 30, 33, 41, 42, 46, 47, 54, 64)
#time13 <- time <- c(1, 4, 5, 8, 18, 19, 20, 21, 26, 28, 30, 43, 51, 61, 71) 

#merge time: 
#time19 <- c(1, 20, 26, 28, 30, 33, 41, 42, 46, 47, 54, 64)
#time <- c(0, 3, 4, 7, 17, 18, 19, 20, 25, 27, 29, 42, 50, 60, 70, 2114, 2134, 2140, 2142, 2144, 2147, 2155, 2156, 2160, 2161, 2168, 2178)

time <- c(2013, 2019)
primary <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2)
scrdat <- ScrData$new(ESP_CH_comb, mesh, primary = primary, time = time)
scrdat # plot & summary

scrdat$n() #unique individuals seen 

scrdat$n_occasions() # get number of occasions in the survey

scrdat$area() # get total area of the mesh in square kilometres 

scrdat$time()

####fit robust model

par <- list(lambda0 ~ 1, 
sigma ~ 1, 
phi ~ 1, 
beta ~ 1,
D ~ 1)

start <- get_start_values(scrdat, model = "JsModel")
start
# create model object 
oo <- JsModel$new(par, scrdat, start)

# compute initial likelihood 
oo$calc_llk()

# fit model 
oo$fit()

#Completed model fitting in 8.253749 mins 
#Checking convergence.......converged 
#Computing variances.......done
#Inferring density..........done
#Computing variances..........done
#Computing confidence intervals..........done

# see results 
oo

oo$get_par("lambda0", k = 1)
oo$get_par("sigma", k = 1)
oo$get_par("phi", k = 1)
oo$get_par("beta", k = 1)
oo$get_par("D")

oo$get_par("phi")

AIC(oo) #131844.4

#2019:
  #100ha - 992600  birds
  #1ha - 9926 birds
  #3.92ha - 9926*3.92 = 38909.92 birds
  
  #2013: 
  #100ha - 783700 birds
  #1ha - 7837 birds
  #3.92 - 7837*3.92 = 30721.04 birds

##################a transient model on combined data######


form <- list(lambda0 ~ 1, 
     sigma ~ 1, 
     beta ~ 1, 
     phi ~ 1, 
     sd ~ 1, 
     D ~ 1)


start <- get_start_values(scrdat, model = "JsTransientModel")
start

trans <- JsTransientModel$new(form, scrdat, start)

trans$fit()

trans
# look at parameters on response scale 
trans$get_par("lambda0", k = 1)
trans$get_par("sigma", k = 1)
trans$get_par("sd", k = 1)
trans$get_par("phi", k = 1)
trans$get_par("beta", k = 1)
trans$get_par("D")

AIC(trans)

#Completed model fitting in 2.227161 hours 
#Checking convergence.......converged 
#Computing variances.......done
#Inferring density..........done
#Computing variances..........done
#Computing confidence intervals..........done
#Warning messages:
#  1: In self$set_mle(mle, V, llk) :
 # Variance estimates not reliable, do a bootstrap.
#2: In sqrt(diag(private$V_)) : NaNs produced
#3: In sqrt(diag(V)) : NaNs produced
#4: In sqrt(diag(V)) : NaNs produced

###To Do 15/11/2019
######repeat all models with boulderscree mask
######check if parameters should be constant

#####################################################################
################repeat robust model on boulder scree mask############
#####################################################################

#reasoning: previously the density of ESP on filfla was calculated for the whole islet
#but in reality ESP only breed in the boulder scree and not on the plateau (circa 1ha)

Boulderscreepoly <- readShapePoly("boulderscree.shp")
plot(Boulderscreepoly)

mesh<-make.mask(detectors, type='trapbuffer', spacing=10, buffer=80, poly=Boulderscreepoly)

#plot(ESP_CH_comb, tracks=T, varycol=FALSE)   # tracks=T specifies that you want lines to connect subsequent sightings of the same individual at different detectors

plot(mesh)
plot(ESP_CH_comb, tracks=T, varycol=FALSE, add=T)
plot(trap, add=T)

time <- c(2013, 2019)
primary <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2)
scrdat <- ScrData$new(ESP_CH_comb, mesh, primary = primary, time = time)
scrdat # plot & summary

scrdat$n() #unique individuals seen 

scrdat$n_occasions() # get number of occasions in the survey

scrdat$area() # get total area of the mesh in square kilometres 
#0.0275

scrdat$time()

####fit robust model

par <- list(lambda0 ~ 1, 
            sigma ~ 1, 
            phi ~ 1, 
            beta ~ 1,
            D ~ 1)

start <- get_start_values(scrdat, model = "JsModel")
start
# create model object 
oo <- JsModel$new(par, scrdat, start)

# compute initial likelihood 
oo$calc_llk()

# fit model 
oo$fit()

#Completed model fitting in 2.197315 mins 
#Checking convergence.......converged 
#Computing variances.......done
#Inferring density..........done
#Computing variances..........done
#Computing confidence intervals..........done

oo

oo$get_par("lambda0", k = 1)
oo$get_par("sigma", k = 1)
oo$get_par("phi", k = 1)
oo$get_par("beta", k = 1)
oo$get_par("D")

oo$get_par("phi")

AIC(oo) # 131695.6

#boulder scree model:
#2019:
#100ha - 1507000 birds
#1ha - 15070 birds
#2.75ha - 15070*2.75 = 41442.5 birds

#2013: 
#100ha - 1202000 birds
#1ha - 12020 birds
#2.75 - 12020*2.75 = 33055 birds

###########################
