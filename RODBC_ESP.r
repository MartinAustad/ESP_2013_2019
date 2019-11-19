####### DATA PREPARATION FROM ACCESS DATABASE ################################################################
####### for CMR analysis of Filfla ESP 2013 & 2019 data using Openpopscr package #############################


library(RODBC)
library(data.table)
library(tidyverse)
library(sf)

#Access database with data from both years 

#setwd("C:\\Users\\Martin\\Documents\\working_folder\\CaptureMarkRecapture\\stormies")
setwd("C:\\STEFFEN\\RSPB\\Malta\\Raw_data")

SP<-odbcConnectAccess2007('ESP_Malta_2013_2019.accdb')
ESP_comb <- sqlQuery(SP, "SELECT * FROM SECR_captures")  
trap_comb <- sqlQuery(SP, "SELECT * FROM SECR_input_detectors")
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
#################################################################

