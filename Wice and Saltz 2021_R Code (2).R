# SELECTION ON HERITABLE SOCIAL NETWORK POSITIONS IS CONTEXT-DEPENDENT IN DROSOPHILA MELANOGASTER#
# CREATED IN R VERSION 3.6.2 "DARK AND STORMY NIGHT", MAC OSX 10.14.3 #
# ERIC WESLEY WICE #
# FEBRUARY 6, 2021 #

# LOAD LIBRARIES
library(readr)
library(readxl)
library(reshape2)
library(data.table)
library(stringr)
library(igraph)
library(lme4)
library(car)
library(scales)
library(plyr)
library(hms)
library(R2admb)
library(glmmADMB)
library(lmerTest)
library(tidyverse)
library(DHARMa)
library(glmmTMB)
library(MuMIn)
library(GGally)
library(survival)
library(coxme)
library(survminer)

# SET WORKING DIRECTORY TO FOLDER CONTAINING ALL DATA #
# LOAD IN VIDEO INFORMATION #
videoinformation <- read_csv("videoinformation.csv", 
                             col_types = cols(video = col_character(), 
                                              videodate = col_date(format = "%m/%d/%y"), 
                                              trial = col_character(), 
                                              day = col_number(), 
                                              time = col_time(format = "%H:%M"), 
                                              lightslefton = col_character(), 
                                              condensation = col_character(), 
                                              flystuck = col_character(), 
                                              frameshift = col_character(), 
                                              flicker = col_character(), 
                                              trackable = col_character(), 
                                              doubleexposed = col_character(), 
                                              potentialswaps = col_number(), 
                                              actualswaps = col_number(), 
                                              glitch = col_character(), 
                                              useable = col_character(), 
                                              frames = col_number(), 
                                              pxmm = col_number(), 
                                              femaleNR = col_number(), 
                                              femaleCO = col_number(), 
                                              femaleCY = col_number(), 
                                              femaleLG = col_number(), 
                                              femaleCG = col_number(), 
                                              femaleUB = col_number(), 
                                              femaleBB = col_number(), 
                                              femaleBP = col_number(), 
                                              femalePP = col_number(), 
                                              femaleTW = col_number(), 
                                              maleNR = col_number(), 
                                              maleCO = col_number(), 
                                              maleCY = col_number(), 
                                              maleLG = col_number(), 
                                              maleCG = col_number(), 
                                              maleUB = col_number(), 
                                              maleBB = col_number(), 
                                              maleBP = col_number(), 
                                              malePP = col_number(), 
                                              maleTW  = col_number()))
videoinformation <- videoinformation[c(1:5,16:38)]
videoinformation <- videoinformation[videoinformation$useable == "Yes",]
videoinformation <- videoinformation[complete.cases(videoinformation$useable),]

# CREATE LISTS AND DATAFRAME FOR MEAN FLY SIZE CALCULATION TO REFER TO #
videos <- unique(videoinformation$video)
flies <- unlist(lapply(1:20, function(x) {return(paste("fly", x, sep=""))}))
sizes <- data.frame(matrix(ncol = 1, nrow = 0))
colnames(sizes) <- c("major_axis_len")

# LOAD IN SIZE DATA AND CALCULATE MEAN FLY SIZE #
for (currentvideo in videos) {
  currentvideodata = paste(currentvideo,"-trackfeat.xls", sep="")
  for (currentfly in flies) {
    currentflydata <- read_excel(currentvideodata, sheet = currentfly, range = cell_cols(4))
    colnames(currentflydata) <- c("major_axis_len")
    currentflydata$video <- currentvideo
    sizes <- rbind(sizes, currentflydata)
  }}
rm(currentflydata,currentfly,currentvideo,currentvideodata)
pxmm <- data.frame(videoinformation$video, videoinformation$pxmm)
colnames(pxmm) <- c("video", "pxmm")
sizes <- merge(sizes, pxmm, by = c("video"))
rm(pxmm)
sizes$major_axis_len <- (sizes$major_axis_len/sizes$pxmm)
sizes$pxmm <- NULL
meansize <- mean(sizes$major_axis_len, na.rm=TRUE)
rm(sizes)

# LOAD IN TRACKING DATA FROM EACH VIDEO, AND CALCULATE ADJACENCY MATRICES, SIZE, AND ACTIVITY INFO #
headers <- c("pos_x", "pos_y", "ori")
activity <- data.frame(matrix(ncol = 1, nrow = 0))
colnames(activity) <- c("activity")
for (currentvideo in videos) {
  currentvideodata = paste(currentvideo,"-trackfeat.xls", sep="")
  # LOAD IN DATA FOR ALL 20 FLIES # -----------------------------------------
  fly1 <-  read_excel(currentvideodata, sheet = "fly1",  range = cell_cols(1:3))
  fly2 <-  read_excel(currentvideodata, sheet = "fly2",  range = cell_cols(1:3))
  fly3 <-  read_excel(currentvideodata, sheet = "fly3",  range = cell_cols(1:3))
  fly4 <-  read_excel(currentvideodata, sheet = "fly4",  range = cell_cols(1:3))
  fly5 <-  read_excel(currentvideodata, sheet = "fly5",  range = cell_cols(1:3))
  fly6 <-  read_excel(currentvideodata, sheet = "fly6",  range = cell_cols(1:3))
  fly7 <-  read_excel(currentvideodata, sheet = "fly7",  range = cell_cols(1:3))
  fly8 <-  read_excel(currentvideodata, sheet = "fly8",  range = cell_cols(1:3))
  fly9 <-  read_excel(currentvideodata, sheet = "fly9",  range = cell_cols(1:3))
  fly10 <- read_excel(currentvideodata, sheet = "fly10", range = cell_cols(1:3))
  fly11 <- read_excel(currentvideodata, sheet = "fly11", range = cell_cols(1:3))
  fly12 <- read_excel(currentvideodata, sheet = "fly12", range = cell_cols(1:3))
  fly13 <- read_excel(currentvideodata, sheet = "fly13", range = cell_cols(1:3))
  fly14 <- read_excel(currentvideodata, sheet = "fly14", range = cell_cols(1:3))
  fly15 <- read_excel(currentvideodata, sheet = "fly15", range = cell_cols(1:3))
  fly16 <- read_excel(currentvideodata, sheet = "fly16", range = cell_cols(1:3))
  fly17 <- read_excel(currentvideodata, sheet = "fly17", range = cell_cols(1:3))
  fly18 <- read_excel(currentvideodata, sheet = "fly18", range = cell_cols(1:3))
  fly19 <- read_excel(currentvideodata, sheet = "fly19", range = cell_cols(1:3))
  fly20 <- read_excel(currentvideodata, sheet = "fly20", range = cell_cols(1:3))
  # RENAME VARIABLES WITHOUT SPACES # ---------------------------------------
  colnames(fly1)  <- headers
  colnames(fly2)  <- headers
  colnames(fly3)  <- headers
  colnames(fly4)  <- headers
  colnames(fly5)  <- headers
  colnames(fly6)  <- headers
  colnames(fly7)  <- headers
  colnames(fly8)  <- headers
  colnames(fly9)  <- headers
  colnames(fly10) <- headers
  colnames(fly11) <- headers
  colnames(fly12) <- headers
  colnames(fly13) <- headers
  colnames(fly14) <- headers
  colnames(fly15) <- headers
  colnames(fly16) <- headers
  colnames(fly17) <- headers
  colnames(fly18) <- headers
  colnames(fly19) <- headers
  colnames(fly20) <- headers
  # CONVERT ALL VARIABLES TO NUMERIC CLASS # ----------------------------------
  fly1 <- as.data.frame(sapply (fly1, as.numeric))
  fly2 <- as.data.frame(sapply (fly2, as.numeric))
  fly3 <- as.data.frame(sapply (fly3, as.numeric))
  fly4 <- as.data.frame(sapply (fly4, as.numeric))
  fly5 <- as.data.frame(sapply (fly5, as.numeric))
  fly6 <- as.data.frame(sapply (fly6, as.numeric))
  fly7 <- as.data.frame(sapply (fly7, as.numeric))
  fly8 <- as.data.frame(sapply (fly8, as.numeric))
  fly9 <- as.data.frame(sapply (fly9, as.numeric))
  fly10 <- as.data.frame(sapply (fly10, as.numeric))
  fly11 <- as.data.frame(sapply (fly11, as.numeric))
  fly12 <- as.data.frame(sapply (fly12, as.numeric))
  fly13 <- as.data.frame(sapply (fly13, as.numeric))
  fly14 <- as.data.frame(sapply (fly14, as.numeric))
  fly15 <- as.data.frame(sapply (fly15, as.numeric))
  fly16 <- as.data.frame(sapply (fly16, as.numeric))
  fly17 <- as.data.frame(sapply (fly17, as.numeric))
  fly18 <- as.data.frame(sapply (fly18, as.numeric))
  fly19 <- as.data.frame(sapply (fly19, as.numeric))
  fly20 <- as.data.frame(sapply (fly20, as.numeric))
  
  # CONVERT pos_x AND pos_y TO MM INSTEAD OF PIXELS #############
  currentvideopxmm <- videoinformation$pxmm[videoinformation$video == currentvideo]
  fly1$pos_x <- (fly1$pos_x / currentvideopxmm)
  fly2$pos_x <- (fly2$pos_x / currentvideopxmm)
  fly3$pos_x <- (fly3$pos_x / currentvideopxmm)
  fly4$pos_x <- (fly4$pos_x / currentvideopxmm)
  fly5$pos_x <- (fly5$pos_x / currentvideopxmm)
  fly6$pos_x <- (fly6$pos_x / currentvideopxmm)
  fly7$pos_x <- (fly7$pos_x / currentvideopxmm)
  fly8$pos_x <- (fly8$pos_x / currentvideopxmm)
  fly9$pos_x <- (fly9$pos_x / currentvideopxmm)
  fly10$pos_x <- (fly10$pos_x / currentvideopxmm)
  fly11$pos_x <- (fly11$pos_x / currentvideopxmm)
  fly12$pos_x <- (fly12$pos_x / currentvideopxmm)
  fly13$pos_x <- (fly13$pos_x / currentvideopxmm)
  fly14$pos_x <- (fly14$pos_x / currentvideopxmm)
  fly15$pos_x <- (fly15$pos_x / currentvideopxmm)
  fly16$pos_x <- (fly16$pos_x / currentvideopxmm)
  fly17$pos_x <- (fly17$pos_x / currentvideopxmm)
  fly18$pos_x <- (fly18$pos_x / currentvideopxmm)
  fly19$pos_x <- (fly19$pos_x / currentvideopxmm)
  fly20$pos_x <- (fly20$pos_x / currentvideopxmm)
  fly1$pos_y <- (fly1$pos_y / currentvideopxmm)
  fly2$pos_y <- (fly2$pos_y / currentvideopxmm)
  fly3$pos_y <- (fly3$pos_y / currentvideopxmm)
  fly4$pos_y <- (fly4$pos_y / currentvideopxmm)
  fly5$pos_y <- (fly5$pos_y / currentvideopxmm)
  fly6$pos_y <- (fly6$pos_y / currentvideopxmm)
  fly7$pos_y <- (fly7$pos_y / currentvideopxmm)
  fly8$pos_y <- (fly8$pos_y / currentvideopxmm)
  fly9$pos_y <- (fly9$pos_y / currentvideopxmm)
  fly10$pos_y <- (fly10$pos_y / currentvideopxmm)
  fly11$pos_y <- (fly11$pos_y / currentvideopxmm)
  fly12$pos_y <- (fly12$pos_y / currentvideopxmm)
  fly13$pos_y <- (fly13$pos_y / currentvideopxmm)
  fly14$pos_y <- (fly14$pos_y / currentvideopxmm)
  fly15$pos_y <- (fly15$pos_y / currentvideopxmm)
  fly16$pos_y <- (fly16$pos_y / currentvideopxmm)
  fly17$pos_y <- (fly17$pos_y / currentvideopxmm)
  fly18$pos_y <- (fly18$pos_y / currentvideopxmm)
  fly19$pos_y <- (fly19$pos_y / currentvideopxmm)
  fly20$pos_y <- (fly20$pos_y / currentvideopxmm)
  
  # CREATE A DATAFRAME THAT CALCULATES THE DISTANCE BETWEEN EVERY PAIRWISE COMBOINATION OF FLIES # --------
  distances = data.frame(sqrt(((fly1$pos_x - fly2$pos_x)^2) + ((fly1$pos_y - fly2$pos_y)^2)))
  colnames(distances) <- c("dist1.2")
  distances$dist1.3 <- sqrt(((fly1$pos_x - fly3$pos_x)^2) + ((fly1$pos_y - fly3$pos_y)^2))
  distances$dist1.4 <- sqrt(((fly1$pos_x - fly4$pos_x)^2) + ((fly1$pos_y - fly4$pos_y)^2))
  distances$dist1.5 <- sqrt(((fly1$pos_x - fly5$pos_x)^2) + ((fly1$pos_y - fly5$pos_y)^2))
  distances$dist1.6 <- sqrt(((fly1$pos_x - fly6$pos_x)^2) + ((fly1$pos_y - fly6$pos_y)^2))
  distances$dist1.7 <- sqrt(((fly1$pos_x - fly7$pos_x)^2) + ((fly1$pos_y - fly7$pos_y)^2))
  distances$dist1.8 <- sqrt(((fly1$pos_x - fly8$pos_x)^2) + ((fly1$pos_y - fly8$pos_y)^2))
  distances$dist1.9 <- sqrt(((fly1$pos_x - fly9$pos_x)^2) + ((fly1$pos_y - fly9$pos_y)^2))
  distances$dist1.10 <- sqrt(((fly1$pos_x - fly10$pos_x)^2) + ((fly1$pos_y - fly10$pos_y)^2))
  distances$dist1.11 <- sqrt(((fly1$pos_x - fly11$pos_x)^2) + ((fly1$pos_y - fly11$pos_y)^2))
  distances$dist1.12 <- sqrt(((fly1$pos_x - fly12$pos_x)^2) + ((fly1$pos_y - fly12$pos_y)^2))
  distances$dist1.13 <- sqrt(((fly1$pos_x - fly13$pos_x)^2) + ((fly1$pos_y - fly13$pos_y)^2))
  distances$dist1.14 <- sqrt(((fly1$pos_x - fly14$pos_x)^2) + ((fly1$pos_y - fly14$pos_y)^2))
  distances$dist1.15 <- sqrt(((fly1$pos_x - fly15$pos_x)^2) + ((fly1$pos_y - fly15$pos_y)^2))
  distances$dist1.16 <- sqrt(((fly1$pos_x - fly16$pos_x)^2) + ((fly1$pos_y - fly16$pos_y)^2))
  distances$dist1.17 <- sqrt(((fly1$pos_x - fly17$pos_x)^2) + ((fly1$pos_y - fly17$pos_y)^2))
  distances$dist1.18 <- sqrt(((fly1$pos_x - fly18$pos_x)^2) + ((fly1$pos_y - fly18$pos_y)^2))
  distances$dist1.19 <- sqrt(((fly1$pos_x - fly19$pos_x)^2) + ((fly1$pos_y - fly19$pos_y)^2))
  distances$dist1.20 <- sqrt(((fly1$pos_x - fly20$pos_x)^2) + ((fly1$pos_y - fly20$pos_y)^2))
  distances$dist2.3 <- sqrt(((fly2$pos_x - fly3$pos_x)^2) + ((fly2$pos_y - fly3$pos_y)^2))
  distances$dist2.4 <- sqrt(((fly2$pos_x - fly4$pos_x)^2) + ((fly2$pos_y - fly4$pos_y)^2))
  distances$dist2.5 <- sqrt(((fly2$pos_x - fly5$pos_x)^2) + ((fly2$pos_y - fly5$pos_y)^2))
  distances$dist2.6 <- sqrt(((fly2$pos_x - fly6$pos_x)^2) + ((fly2$pos_y - fly6$pos_y)^2))
  distances$dist2.7 <- sqrt(((fly2$pos_x - fly7$pos_x)^2) + ((fly2$pos_y - fly7$pos_y)^2))
  distances$dist2.8 <- sqrt(((fly2$pos_x - fly8$pos_x)^2) + ((fly2$pos_y - fly8$pos_y)^2))
  distances$dist2.9 <- sqrt(((fly2$pos_x - fly10$pos_x)^2) + ((fly2$pos_y - fly9$pos_y)^2))
  distances$dist2.10 <- sqrt(((fly2$pos_x - fly10$pos_x)^2) + ((fly2$pos_y - fly10$pos_y)^2))
  distances$dist2.11 <- sqrt(((fly2$pos_x - fly11$pos_x)^2) + ((fly2$pos_y - fly11$pos_y)^2))
  distances$dist2.12 <- sqrt(((fly2$pos_x - fly12$pos_x)^2) + ((fly2$pos_y - fly12$pos_y)^2))
  distances$dist2.13 <- sqrt(((fly2$pos_x - fly13$pos_x)^2) + ((fly2$pos_y - fly13$pos_y)^2))
  distances$dist2.14 <- sqrt(((fly2$pos_x - fly14$pos_x)^2) + ((fly2$pos_y - fly14$pos_y)^2))
  distances$dist2.15 <- sqrt(((fly2$pos_x - fly15$pos_x)^2) + ((fly2$pos_y - fly15$pos_y)^2))
  distances$dist2.16 <- sqrt(((fly2$pos_x - fly16$pos_x)^2) + ((fly2$pos_y - fly16$pos_y)^2))
  distances$dist2.17 <- sqrt(((fly2$pos_x - fly17$pos_x)^2) + ((fly2$pos_y - fly17$pos_y)^2))
  distances$dist2.18 <- sqrt(((fly2$pos_x - fly18$pos_x)^2) + ((fly2$pos_y - fly18$pos_y)^2))
  distances$dist2.19 <- sqrt(((fly2$pos_x - fly19$pos_x)^2) + ((fly2$pos_y - fly19$pos_y)^2))
  distances$dist2.20 <- sqrt(((fly2$pos_x - fly20$pos_x)^2) + ((fly2$pos_y - fly20$pos_y)^2))
  distances$dist3.4 <- sqrt(((fly3$pos_x - fly4$pos_x)^2) + ((fly3$pos_y - fly4$pos_y)^2))
  distances$dist3.5 <- sqrt(((fly3$pos_x - fly5$pos_x)^2) + ((fly3$pos_y - fly5$pos_y)^2))
  distances$dist3.6 <- sqrt(((fly3$pos_x - fly6$pos_x)^2) + ((fly3$pos_y - fly6$pos_y)^2))
  distances$dist3.7 <- sqrt(((fly3$pos_x - fly7$pos_x)^2) + ((fly3$pos_y - fly7$pos_y)^2))
  distances$dist3.8 <- sqrt(((fly3$pos_x - fly8$pos_x)^2) + ((fly3$pos_y - fly8$pos_y)^2))
  distances$dist3.9 <- sqrt(((fly3$pos_x - fly9$pos_x)^2) + ((fly3$pos_y - fly9$pos_y)^2))
  distances$dist3.10 <- sqrt(((fly3$pos_x - fly10$pos_x)^2) + ((fly3$pos_y - fly10$pos_y)^2))
  distances$dist3.11 <- sqrt(((fly3$pos_x - fly11$pos_x)^2) + ((fly3$pos_y - fly11$pos_y)^2))
  distances$dist3.12 <- sqrt(((fly3$pos_x - fly12$pos_x)^2) + ((fly3$pos_y - fly12$pos_y)^2))
  distances$dist3.13 <- sqrt(((fly3$pos_x - fly13$pos_x)^2) + ((fly3$pos_y - fly13$pos_y)^2))
  distances$dist3.14 <- sqrt(((fly3$pos_x - fly14$pos_x)^2) + ((fly3$pos_y - fly14$pos_y)^2))
  distances$dist3.15 <- sqrt(((fly3$pos_x - fly15$pos_x)^2) + ((fly3$pos_y - fly15$pos_y)^2))
  distances$dist3.16 <- sqrt(((fly3$pos_x - fly16$pos_x)^2) + ((fly3$pos_y - fly16$pos_y)^2))
  distances$dist3.17 <- sqrt(((fly3$pos_x - fly17$pos_x)^2) + ((fly3$pos_y - fly17$pos_y)^2))
  distances$dist3.18 <- sqrt(((fly3$pos_x - fly18$pos_x)^2) + ((fly3$pos_y - fly18$pos_y)^2))
  distances$dist3.19 <- sqrt(((fly3$pos_x - fly19$pos_x)^2) + ((fly3$pos_y - fly19$pos_y)^2))
  distances$dist3.20 <- sqrt(((fly3$pos_x - fly20$pos_x)^2) + ((fly3$pos_y - fly20$pos_y)^2))
  distances$dist4.5 <- sqrt(((fly4$pos_x - fly5$pos_x)^2) + ((fly4$pos_y - fly5$pos_y)^2))
  distances$dist4.6 <- sqrt(((fly4$pos_x - fly6$pos_x)^2) + ((fly4$pos_y - fly6$pos_y)^2))
  distances$dist4.7 <- sqrt(((fly4$pos_x - fly7$pos_x)^2) + ((fly4$pos_y - fly7$pos_y)^2))
  distances$dist4.8 <- sqrt(((fly4$pos_x - fly8$pos_x)^2) + ((fly4$pos_y - fly8$pos_y)^2))
  distances$dist4.9 <- sqrt(((fly4$pos_x - fly9$pos_x)^2) + ((fly4$pos_y - fly9$pos_y)^2))
  distances$dist4.10 <- sqrt(((fly4$pos_x - fly10$pos_x)^2) + ((fly4$pos_y - fly10$pos_y)^2))
  distances$dist4.11 <- sqrt(((fly4$pos_x - fly11$pos_x)^2) + ((fly4$pos_y - fly11$pos_y)^2))
  distances$dist4.12 <- sqrt(((fly4$pos_x - fly12$pos_x)^2) + ((fly4$pos_y - fly12$pos_y)^2))
  distances$dist4.13 <- sqrt(((fly4$pos_x - fly13$pos_x)^2) + ((fly4$pos_y - fly13$pos_y)^2))
  distances$dist4.14 <- sqrt(((fly4$pos_x - fly14$pos_x)^2) + ((fly4$pos_y - fly14$pos_y)^2))
  distances$dist4.15 <- sqrt(((fly4$pos_x - fly15$pos_x)^2) + ((fly4$pos_y - fly15$pos_y)^2))
  distances$dist4.16 <- sqrt(((fly4$pos_x - fly16$pos_x)^2) + ((fly4$pos_y - fly16$pos_y)^2))
  distances$dist4.17 <- sqrt(((fly4$pos_x - fly17$pos_x)^2) + ((fly4$pos_y - fly17$pos_y)^2))
  distances$dist4.18 <- sqrt(((fly4$pos_x - fly18$pos_x)^2) + ((fly4$pos_y - fly18$pos_y)^2))
  distances$dist4.19 <- sqrt(((fly4$pos_x - fly19$pos_x)^2) + ((fly4$pos_y - fly19$pos_y)^2))
  distances$dist4.20 <- sqrt(((fly4$pos_x - fly20$pos_x)^2) + ((fly4$pos_y - fly20$pos_y)^2))
  distances$dist5.6 <- sqrt(((fly5$pos_x - fly6$pos_x)^2) + ((fly5$pos_y - fly6$pos_y)^2))
  distances$dist5.7 <- sqrt(((fly5$pos_x - fly7$pos_x)^2) + ((fly5$pos_y - fly7$pos_y)^2))
  distances$dist5.8 <- sqrt(((fly5$pos_x - fly8$pos_x)^2) + ((fly5$pos_y - fly8$pos_y)^2))
  distances$dist5.9 <- sqrt(((fly5$pos_x - fly9$pos_x)^2) + ((fly5$pos_y - fly9$pos_y)^2))
  distances$dist5.10 <- sqrt(((fly5$pos_x - fly10$pos_x)^2) + ((fly5$pos_y - fly10$pos_y)^2))
  distances$dist5.11 <- sqrt(((fly5$pos_x - fly11$pos_x)^2) + ((fly5$pos_y - fly11$pos_y)^2))
  distances$dist5.12 <- sqrt(((fly5$pos_x - fly12$pos_x)^2) + ((fly5$pos_y - fly12$pos_y)^2))
  distances$dist5.13 <- sqrt(((fly5$pos_x - fly13$pos_x)^2) + ((fly5$pos_y - fly13$pos_y)^2))
  distances$dist5.14 <- sqrt(((fly5$pos_x - fly14$pos_x)^2) + ((fly5$pos_y - fly14$pos_y)^2))
  distances$dist5.15 <- sqrt(((fly5$pos_x - fly15$pos_x)^2) + ((fly5$pos_y - fly15$pos_y)^2))
  distances$dist5.16 <- sqrt(((fly5$pos_x - fly16$pos_x)^2) + ((fly5$pos_y - fly16$pos_y)^2))
  distances$dist5.17 <- sqrt(((fly5$pos_x - fly17$pos_x)^2) + ((fly5$pos_y - fly17$pos_y)^2))
  distances$dist5.18 <- sqrt(((fly5$pos_x - fly18$pos_x)^2) + ((fly5$pos_y - fly18$pos_y)^2))
  distances$dist5.19 <- sqrt(((fly5$pos_x - fly19$pos_x)^2) + ((fly5$pos_y - fly19$pos_y)^2))
  distances$dist5.20 <- sqrt(((fly5$pos_x - fly20$pos_x)^2) + ((fly5$pos_y - fly20$pos_y)^2))
  distances$dist6.7 <- sqrt(((fly6$pos_x - fly7$pos_x)^2) + ((fly6$pos_y - fly7$pos_y)^2))
  distances$dist6.8 <- sqrt(((fly6$pos_x - fly8$pos_x)^2) + ((fly6$pos_y - fly8$pos_y)^2))
  distances$dist6.9 <- sqrt(((fly6$pos_x - fly9$pos_x)^2) + ((fly6$pos_y - fly9$pos_y)^2))
  distances$dist6.10 <- sqrt(((fly6$pos_x - fly10$pos_x)^2) + ((fly6$pos_y - fly10$pos_y)^2))
  distances$dist6.11 <- sqrt(((fly6$pos_x - fly11$pos_x)^2) + ((fly6$pos_y - fly11$pos_y)^2))
  distances$dist6.12 <- sqrt(((fly6$pos_x - fly12$pos_x)^2) + ((fly6$pos_y - fly12$pos_y)^2))
  distances$dist6.13 <- sqrt(((fly6$pos_x - fly13$pos_x)^2) + ((fly6$pos_y - fly13$pos_y)^2))
  distances$dist6.14 <- sqrt(((fly6$pos_x - fly14$pos_x)^2) + ((fly6$pos_y - fly14$pos_y)^2))
  distances$dist6.15 <- sqrt(((fly6$pos_x - fly15$pos_x)^2) + ((fly6$pos_y - fly15$pos_y)^2))
  distances$dist6.16 <- sqrt(((fly6$pos_x - fly16$pos_x)^2) + ((fly6$pos_y - fly16$pos_y)^2))
  distances$dist6.17 <- sqrt(((fly6$pos_x - fly17$pos_x)^2) + ((fly6$pos_y - fly17$pos_y)^2))
  distances$dist6.18 <- sqrt(((fly6$pos_x - fly18$pos_x)^2) + ((fly6$pos_y - fly18$pos_y)^2))
  distances$dist6.19 <- sqrt(((fly6$pos_x - fly19$pos_x)^2) + ((fly6$pos_y - fly19$pos_y)^2))
  distances$dist6.20 <- sqrt(((fly6$pos_x - fly20$pos_x)^2) + ((fly6$pos_y - fly20$pos_y)^2))
  distances$dist7.8 <- sqrt(((fly7$pos_x - fly8$pos_x)^2) + ((fly7$pos_y - fly8$pos_y)^2))
  distances$dist7.9 <- sqrt(((fly7$pos_x - fly9$pos_x)^2) + ((fly7$pos_y - fly9$pos_y)^2))
  distances$dist7.10 <- sqrt(((fly7$pos_x - fly10$pos_x)^2) + ((fly7$pos_y - fly10$pos_y)^2))
  distances$dist7.11 <- sqrt(((fly7$pos_x - fly11$pos_x)^2) + ((fly7$pos_y - fly11$pos_y)^2))
  distances$dist7.12 <- sqrt(((fly7$pos_x - fly12$pos_x)^2) + ((fly7$pos_y - fly12$pos_y)^2))
  distances$dist7.13 <- sqrt(((fly7$pos_x - fly13$pos_x)^2) + ((fly7$pos_y - fly13$pos_y)^2))
  distances$dist7.14 <- sqrt(((fly7$pos_x - fly14$pos_x)^2) + ((fly7$pos_y - fly14$pos_y)^2))
  distances$dist7.15 <- sqrt(((fly7$pos_x - fly15$pos_x)^2) + ((fly7$pos_y - fly15$pos_y)^2))
  distances$dist7.16 <- sqrt(((fly7$pos_x - fly16$pos_x)^2) + ((fly7$pos_y - fly16$pos_y)^2))
  distances$dist7.17 <- sqrt(((fly7$pos_x - fly17$pos_x)^2) + ((fly7$pos_y - fly17$pos_y)^2))
  distances$dist7.18 <- sqrt(((fly7$pos_x - fly18$pos_x)^2) + ((fly7$pos_y - fly18$pos_y)^2))
  distances$dist7.19 <- sqrt(((fly7$pos_x - fly19$pos_x)^2) + ((fly7$pos_y - fly19$pos_y)^2))
  distances$dist7.20 <- sqrt(((fly7$pos_x - fly20$pos_x)^2) + ((fly7$pos_y - fly20$pos_y)^2))
  distances$dist8.9 <- sqrt(((fly8$pos_x - fly9$pos_x)^2) + ((fly8$pos_y - fly9$pos_y)^2))
  distances$dist8.10 <- sqrt(((fly8$pos_x - fly10$pos_x)^2) + ((fly8$pos_y - fly10$pos_y)^2))
  distances$dist8.11 <- sqrt(((fly8$pos_x - fly11$pos_x)^2) + ((fly8$pos_y - fly11$pos_y)^2))
  distances$dist8.12 <- sqrt(((fly8$pos_x - fly12$pos_x)^2) + ((fly8$pos_y - fly12$pos_y)^2))
  distances$dist8.13 <- sqrt(((fly8$pos_x - fly13$pos_x)^2) + ((fly8$pos_y - fly13$pos_y)^2))
  distances$dist8.14 <- sqrt(((fly8$pos_x - fly14$pos_x)^2) + ((fly8$pos_y - fly14$pos_y)^2))
  distances$dist8.15 <- sqrt(((fly8$pos_x - fly15$pos_x)^2) + ((fly8$pos_y - fly15$pos_y)^2))
  distances$dist8.16 <- sqrt(((fly8$pos_x - fly16$pos_x)^2) + ((fly8$pos_y - fly16$pos_y)^2))
  distances$dist8.17 <- sqrt(((fly8$pos_x - fly17$pos_x)^2) + ((fly8$pos_y - fly17$pos_y)^2))
  distances$dist8.18 <- sqrt(((fly8$pos_x - fly18$pos_x)^2) + ((fly8$pos_y - fly18$pos_y)^2))
  distances$dist8.19 <- sqrt(((fly8$pos_x - fly19$pos_x)^2) + ((fly8$pos_y - fly19$pos_y)^2))
  distances$dist8.20 <- sqrt(((fly8$pos_x - fly20$pos_x)^2) + ((fly8$pos_y - fly20$pos_y)^2))
  distances$dist9.10 <- sqrt(((fly9$pos_x - fly10$pos_x)^2) + ((fly9$pos_y - fly10$pos_y)^2))
  distances$dist9.11 <- sqrt(((fly9$pos_x - fly11$pos_x)^2) + ((fly9$pos_y - fly11$pos_y)^2))
  distances$dist9.12 <- sqrt(((fly9$pos_x - fly12$pos_x)^2) + ((fly9$pos_y - fly12$pos_y)^2))
  distances$dist9.13 <- sqrt(((fly9$pos_x - fly13$pos_x)^2) + ((fly9$pos_y - fly13$pos_y)^2))
  distances$dist9.14 <- sqrt(((fly9$pos_x - fly14$pos_x)^2) + ((fly9$pos_y - fly14$pos_y)^2))
  distances$dist9.15 <- sqrt(((fly9$pos_x - fly15$pos_x)^2) + ((fly9$pos_y - fly15$pos_y)^2))
  distances$dist9.16 <- sqrt(((fly9$pos_x - fly16$pos_x)^2) + ((fly9$pos_y - fly16$pos_y)^2))
  distances$dist9.17 <- sqrt(((fly9$pos_x - fly17$pos_x)^2) + ((fly9$pos_y - fly17$pos_y)^2))
  distances$dist9.18 <- sqrt(((fly9$pos_x - fly18$pos_x)^2) + ((fly9$pos_y - fly18$pos_y)^2))
  distances$dist9.19 <- sqrt(((fly9$pos_x - fly19$pos_x)^2) + ((fly9$pos_y - fly19$pos_y)^2))
  distances$dist9.20 <- sqrt(((fly9$pos_x - fly20$pos_x)^2) + ((fly9$pos_y - fly20$pos_y)^2))
  distances$dist10.11 <- sqrt(((fly10$pos_x - fly11$pos_x)^2) + ((fly10$pos_y - fly11$pos_y)^2))
  distances$dist10.12 <- sqrt(((fly10$pos_x - fly12$pos_x)^2) + ((fly10$pos_y - fly12$pos_y)^2))
  distances$dist10.13 <- sqrt(((fly10$pos_x - fly13$pos_x)^2) + ((fly10$pos_y - fly13$pos_y)^2))
  distances$dist10.14 <- sqrt(((fly10$pos_x - fly14$pos_x)^2) + ((fly10$pos_y - fly14$pos_y)^2))
  distances$dist10.15 <- sqrt(((fly10$pos_x - fly15$pos_x)^2) + ((fly10$pos_y - fly15$pos_y)^2))
  distances$dist10.16 <- sqrt(((fly10$pos_x - fly16$pos_x)^2) + ((fly10$pos_y - fly16$pos_y)^2))
  distances$dist10.17 <- sqrt(((fly10$pos_x - fly17$pos_x)^2) + ((fly10$pos_y - fly17$pos_y)^2))
  distances$dist10.18 <- sqrt(((fly10$pos_x - fly18$pos_x)^2) + ((fly10$pos_y - fly18$pos_y)^2))
  distances$dist10.19 <- sqrt(((fly10$pos_x - fly19$pos_x)^2) + ((fly10$pos_y - fly19$pos_y)^2))
  distances$dist10.20 <- sqrt(((fly10$pos_x - fly20$pos_x)^2) + ((fly10$pos_y - fly20$pos_y)^2))
  distances$dist11.12 <- sqrt(((fly11$pos_x - fly12$pos_x)^2) + ((fly11$pos_y - fly12$pos_y)^2))
  distances$dist11.13 <- sqrt(((fly11$pos_x - fly13$pos_x)^2) + ((fly11$pos_y - fly13$pos_y)^2))
  distances$dist11.14 <- sqrt(((fly11$pos_x - fly14$pos_x)^2) + ((fly11$pos_y - fly14$pos_y)^2))
  distances$dist11.15 <- sqrt(((fly11$pos_x - fly15$pos_x)^2) + ((fly11$pos_y - fly15$pos_y)^2))
  distances$dist11.16 <- sqrt(((fly11$pos_x - fly16$pos_x)^2) + ((fly11$pos_y - fly16$pos_y)^2))
  distances$dist11.17 <- sqrt(((fly11$pos_x - fly17$pos_x)^2) + ((fly11$pos_y - fly17$pos_y)^2))
  distances$dist11.18 <- sqrt(((fly11$pos_x - fly18$pos_x)^2) + ((fly11$pos_y - fly18$pos_y)^2))
  distances$dist11.19 <- sqrt(((fly11$pos_x - fly19$pos_x)^2) + ((fly11$pos_y - fly19$pos_y)^2))
  distances$dist11.20 <- sqrt(((fly11$pos_x - fly20$pos_x)^2) + ((fly11$pos_y - fly20$pos_y)^2))
  distances$dist12.13 <- sqrt(((fly12$pos_x - fly13$pos_x)^2) + ((fly12$pos_y - fly13$pos_y)^2))
  distances$dist12.14 <- sqrt(((fly12$pos_x - fly14$pos_x)^2) + ((fly12$pos_y - fly14$pos_y)^2))
  distances$dist12.15 <- sqrt(((fly12$pos_x - fly15$pos_x)^2) + ((fly12$pos_y - fly15$pos_y)^2))
  distances$dist12.16 <- sqrt(((fly12$pos_x - fly16$pos_x)^2) + ((fly12$pos_y - fly16$pos_y)^2))
  distances$dist12.17 <- sqrt(((fly12$pos_x - fly17$pos_x)^2) + ((fly12$pos_y - fly17$pos_y)^2))
  distances$dist12.18 <- sqrt(((fly12$pos_x - fly18$pos_x)^2) + ((fly12$pos_y - fly18$pos_y)^2))
  distances$dist12.19 <- sqrt(((fly12$pos_x - fly19$pos_x)^2) + ((fly12$pos_y - fly19$pos_y)^2))
  distances$dist12.20 <- sqrt(((fly12$pos_x - fly20$pos_x)^2) + ((fly12$pos_y - fly20$pos_y)^2))
  distances$dist13.14 <- sqrt(((fly13$pos_x - fly14$pos_x)^2) + ((fly13$pos_y - fly14$pos_y)^2))
  distances$dist13.15 <- sqrt(((fly13$pos_x - fly15$pos_x)^2) + ((fly13$pos_y - fly15$pos_y)^2))
  distances$dist13.16 <- sqrt(((fly13$pos_x - fly16$pos_x)^2) + ((fly13$pos_y - fly16$pos_y)^2))
  distances$dist13.17 <- sqrt(((fly13$pos_x - fly17$pos_x)^2) + ((fly13$pos_y - fly17$pos_y)^2))
  distances$dist13.18 <- sqrt(((fly13$pos_x - fly18$pos_x)^2) + ((fly13$pos_y - fly18$pos_y)^2))
  distances$dist13.19 <- sqrt(((fly13$pos_x - fly19$pos_x)^2) + ((fly13$pos_y - fly19$pos_y)^2))
  distances$dist13.20 <- sqrt(((fly13$pos_x - fly20$pos_x)^2) + ((fly13$pos_y - fly20$pos_y)^2))
  distances$dist14.15 <- sqrt(((fly14$pos_x - fly15$pos_x)^2) + ((fly14$pos_y - fly15$pos_y)^2))
  distances$dist14.16 <- sqrt(((fly14$pos_x - fly16$pos_x)^2) + ((fly14$pos_y - fly16$pos_y)^2))
  distances$dist14.17 <- sqrt(((fly14$pos_x - fly17$pos_x)^2) + ((fly14$pos_y - fly17$pos_y)^2))
  distances$dist14.18 <- sqrt(((fly14$pos_x - fly18$pos_x)^2) + ((fly14$pos_y - fly18$pos_y)^2))
  distances$dist14.19 <- sqrt(((fly14$pos_x - fly19$pos_x)^2) + ((fly14$pos_y - fly19$pos_y)^2))
  distances$dist14.20 <- sqrt(((fly14$pos_x - fly20$pos_x)^2) + ((fly14$pos_y - fly20$pos_y)^2))
  distances$dist15.16 <- sqrt(((fly15$pos_x - fly16$pos_x)^2) + ((fly15$pos_y - fly16$pos_y)^2))
  distances$dist15.17 <- sqrt(((fly15$pos_x - fly17$pos_x)^2) + ((fly15$pos_y - fly17$pos_y)^2))
  distances$dist15.18 <- sqrt(((fly15$pos_x - fly18$pos_x)^2) + ((fly15$pos_y - fly18$pos_y)^2))
  distances$dist15.19 <- sqrt(((fly15$pos_x - fly19$pos_x)^2) + ((fly15$pos_y - fly19$pos_y)^2))
  distances$dist15.20 <- sqrt(((fly15$pos_x - fly20$pos_x)^2) + ((fly15$pos_y - fly20$pos_y)^2))
  distances$dist16.17 <- sqrt(((fly16$pos_x - fly17$pos_x)^2) + ((fly16$pos_y - fly17$pos_y)^2))
  distances$dist16.18 <- sqrt(((fly16$pos_x - fly18$pos_x)^2) + ((fly16$pos_y - fly18$pos_y)^2))
  distances$dist16.19 <- sqrt(((fly16$pos_x - fly19$pos_x)^2) + ((fly16$pos_y - fly19$pos_y)^2))
  distances$dist16.20 <- sqrt(((fly16$pos_x - fly20$pos_x)^2) + ((fly16$pos_y - fly20$pos_y)^2))
  distances$dist17.18 <- sqrt(((fly17$pos_x - fly18$pos_x)^2) + ((fly17$pos_y - fly18$pos_y)^2))
  distances$dist17.19 <- sqrt(((fly17$pos_x - fly19$pos_x)^2) + ((fly17$pos_y - fly19$pos_y)^2))
  distances$dist17.20 <- sqrt(((fly17$pos_x - fly20$pos_x)^2) + ((fly17$pos_y - fly20$pos_y)^2))
  distances$dist18.19 <- sqrt(((fly18$pos_x - fly19$pos_x)^2) + ((fly18$pos_y - fly19$pos_y)^2))
  distances$dist18.20 <- sqrt(((fly18$pos_x - fly20$pos_x)^2) + ((fly18$pos_y - fly20$pos_y)^2))
  distances$dist19.20 <- sqrt(((fly19$pos_x - fly20$pos_x)^2) + ((fly19$pos_y - fly20$pos_y)^2))
  
  # CREATE A DATAFRAME WHERE IF DISTANCE BETWEEN TWO FLIES < (2.5 X MEANSIZE) THEN 1, IF NO = 0 # --------
  distanceinteraction = data.frame(ifelse(distances$dist1.2 < (2.5*meansize), c(1), c(0)))
  colnames(distanceinteraction) <- c("int1.2")
  distanceinteraction$int1.3 <- ifelse(distances$dist1.3 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int1.4 <- ifelse(distances$dist1.4 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int1.5 <- ifelse(distances$dist1.5 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int1.6 <- ifelse(distances$dist1.6 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int1.7 <- ifelse(distances$dist1.7 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int1.8 <- ifelse(distances$dist1.8 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int1.9 <- ifelse(distances$dist1.9 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int1.10 <- ifelse(distances$dist1.10 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int1.11 <- ifelse(distances$dist1.11 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int1.12 <- ifelse(distances$dist1.12 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int1.13 <- ifelse(distances$dist1.13 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int1.14 <- ifelse(distances$dist1.14 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int1.15 <- ifelse(distances$dist1.15 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int1.16 <- ifelse(distances$dist1.16 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int1.17 <- ifelse(distances$dist1.17 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int1.18 <- ifelse(distances$dist1.18 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int1.19 <- ifelse(distances$dist1.19 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int1.20 <- ifelse(distances$dist1.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int2.1 <- ifelse(distances$dist1.2 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int2.3 <- ifelse(distances$dist2.3 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int2.4 <- ifelse(distances$dist2.4 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int2.5 <- ifelse(distances$dist2.5 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int2.6 <- ifelse(distances$dist2.6 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int2.7 <- ifelse(distances$dist2.7 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int2.8 <- ifelse(distances$dist2.8 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int2.9 <- ifelse(distances$dist2.9 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int2.10 <- ifelse(distances$dist2.10 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int2.11 <- ifelse(distances$dist2.11 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int2.12 <- ifelse(distances$dist2.12 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int2.13 <- ifelse(distances$dist2.13 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int2.14 <- ifelse(distances$dist2.14 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int2.15 <- ifelse(distances$dist2.15 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int2.16 <- ifelse(distances$dist2.16 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int2.17 <- ifelse(distances$dist2.17 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int2.18 <- ifelse(distances$dist2.18 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int2.19 <- ifelse(distances$dist2.19 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int2.20 <- ifelse(distances$dist2.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int3.1 <- ifelse(distances$dist1.3 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int3.2 <- ifelse(distances$dist2.3 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int3.4 <- ifelse(distances$dist3.4 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int3.5 <- ifelse(distances$dist3.5 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int3.6 <- ifelse(distances$dist3.6 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int3.7 <- ifelse(distances$dist3.7 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int3.8 <- ifelse(distances$dist3.8 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int3.9 <- ifelse(distances$dist3.9 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int3.10 <- ifelse(distances$dist3.10 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int3.11 <- ifelse(distances$dist3.11 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int3.12 <- ifelse(distances$dist3.12 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int3.13 <- ifelse(distances$dist3.13 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int3.14 <- ifelse(distances$dist3.14 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int3.15 <- ifelse(distances$dist3.15 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int3.16 <- ifelse(distances$dist3.16 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int3.17 <- ifelse(distances$dist3.17 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int3.18 <- ifelse(distances$dist3.18 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int3.19 <- ifelse(distances$dist3.19 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int3.20 <- ifelse(distances$dist3.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int4.1 <- ifelse(distances$dist1.4 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int4.2 <- ifelse(distances$dist2.4 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int4.3 <- ifelse(distances$dist3.4 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int4.5 <- ifelse(distances$dist4.5 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int4.6 <- ifelse(distances$dist4.6 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int4.7 <- ifelse(distances$dist4.7 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int4.8 <- ifelse(distances$dist4.8 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int4.9 <- ifelse(distances$dist4.9 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int4.10 <- ifelse(distances$dist4.10 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int4.11 <- ifelse(distances$dist4.11 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int4.12 <- ifelse(distances$dist4.12 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int4.13 <- ifelse(distances$dist4.13 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int4.14 <- ifelse(distances$dist4.14 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int4.15 <- ifelse(distances$dist4.15 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int4.16 <- ifelse(distances$dist4.16 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int4.17 <- ifelse(distances$dist4.17 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int4.18 <- ifelse(distances$dist4.18 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int4.19 <- ifelse(distances$dist4.19 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int4.20 <- ifelse(distances$dist4.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int5.1 <- ifelse(distances$dist1.5 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int5.2 <- ifelse(distances$dist2.5 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int5.3 <- ifelse(distances$dist3.5 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int5.4 <- ifelse(distances$dist4.5 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int5.6 <- ifelse(distances$dist5.6 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int5.7 <- ifelse(distances$dist5.7 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int5.8 <- ifelse(distances$dist5.8 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int5.9 <- ifelse(distances$dist5.9 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int5.10 <- ifelse(distances$dist5.10 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int5.11 <- ifelse(distances$dist5.11 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int5.12 <- ifelse(distances$dist5.12 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int5.13 <- ifelse(distances$dist5.13 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int5.14 <- ifelse(distances$dist5.14 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int5.15 <- ifelse(distances$dist5.15 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int5.16 <- ifelse(distances$dist5.16 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int5.17 <- ifelse(distances$dist5.17 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int5.18 <- ifelse(distances$dist5.18 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int5.19 <- ifelse(distances$dist5.19 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int5.20 <- ifelse(distances$dist5.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int6.1 <- ifelse(distances$dist1.6 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int6.2 <- ifelse(distances$dist2.6 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int6.3 <- ifelse(distances$dist3.6 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int6.4 <- ifelse(distances$dist4.6 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int6.5 <- ifelse(distances$dist5.6 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int6.7 <- ifelse(distances$dist6.7 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int6.8 <- ifelse(distances$dist6.8 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int6.9 <- ifelse(distances$dist6.9 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int6.10 <- ifelse(distances$dist6.10 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int6.11 <- ifelse(distances$dist6.11 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int6.12 <- ifelse(distances$dist6.12 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int6.13 <- ifelse(distances$dist6.13 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int6.14 <- ifelse(distances$dist6.14 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int6.15 <- ifelse(distances$dist6.15 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int6.16 <- ifelse(distances$dist6.16 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int6.17 <- ifelse(distances$dist6.17 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int6.18 <- ifelse(distances$dist6.18 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int6.19 <- ifelse(distances$dist6.19 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int6.20 <- ifelse(distances$dist6.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int7.1 <- ifelse(distances$dist1.7 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int7.2 <- ifelse(distances$dist2.7 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int7.3 <- ifelse(distances$dist3.7 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int7.4 <- ifelse(distances$dist4.7 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int7.5 <- ifelse(distances$dist5.7 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int7.6 <- ifelse(distances$dist6.7 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int7.8 <- ifelse(distances$dist7.8 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int7.9 <- ifelse(distances$dist7.9 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int7.10 <- ifelse(distances$dist7.10 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int7.11 <- ifelse(distances$dist7.11 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int7.12 <- ifelse(distances$dist7.12 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int7.13 <- ifelse(distances$dist7.13 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int7.14 <- ifelse(distances$dist7.14 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int7.15 <- ifelse(distances$dist7.15 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int7.16 <- ifelse(distances$dist7.16 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int7.17 <- ifelse(distances$dist7.17 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int7.18 <- ifelse(distances$dist7.18 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int7.19 <- ifelse(distances$dist7.19 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int7.20 <- ifelse(distances$dist7.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int8.1 <- ifelse(distances$dist1.8 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int8.2 <- ifelse(distances$dist2.8 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int8.3 <- ifelse(distances$dist3.8 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int8.4 <- ifelse(distances$dist4.8 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int8.5 <- ifelse(distances$dist5.8 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int8.6 <- ifelse(distances$dist6.8 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int8.7 <- ifelse(distances$dist7.8 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int8.9 <- ifelse(distances$dist8.9 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int8.10 <- ifelse(distances$dist8.10 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int8.11 <- ifelse(distances$dist8.11 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int8.12 <- ifelse(distances$dist8.12 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int8.13 <- ifelse(distances$dist8.13 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int8.14 <- ifelse(distances$dist8.14 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int8.15 <- ifelse(distances$dist8.15 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int8.16 <- ifelse(distances$dist8.16 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int8.17 <- ifelse(distances$dist8.17 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int8.18 <- ifelse(distances$dist8.18 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int8.19 <- ifelse(distances$dist8.19 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int8.20 <- ifelse(distances$dist8.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int9.1 <- ifelse(distances$dist1.9 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int9.2 <- ifelse(distances$dist2.9 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int9.3 <- ifelse(distances$dist3.9 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int9.4 <- ifelse(distances$dist4.9 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int9.5 <- ifelse(distances$dist5.9 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int9.6 <- ifelse(distances$dist6.9 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int9.7 <- ifelse(distances$dist7.9 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int9.8 <- ifelse(distances$dist8.9 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int9.10 <- ifelse(distances$dist9.10 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int9.11 <- ifelse(distances$dist9.11 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int9.12 <- ifelse(distances$dist9.12 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int9.13 <- ifelse(distances$dist9.13 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int9.14 <- ifelse(distances$dist9.14 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int9.15 <- ifelse(distances$dist9.15 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int9.16 <- ifelse(distances$dist9.16 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int9.17 <- ifelse(distances$dist9.17 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int9.18 <- ifelse(distances$dist9.18 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int9.19 <- ifelse(distances$dist9.19 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int9.20 <- ifelse(distances$dist9.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int10.1 <- ifelse(distances$dist1.10 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int10.2 <- ifelse(distances$dist2.10 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int10.3 <- ifelse(distances$dist3.10 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int10.4 <- ifelse(distances$dist4.10 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int10.5 <- ifelse(distances$dist5.10 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int10.6 <- ifelse(distances$dist6.10 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int10.7 <- ifelse(distances$dist7.10 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int10.8 <- ifelse(distances$dist8.10 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int10.9 <- ifelse(distances$dist9.10 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int10.11 <- ifelse(distances$dist10.11 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int10.12 <- ifelse(distances$dist10.12 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int10.13 <- ifelse(distances$dist10.13 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int10.14 <- ifelse(distances$dist10.14 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int10.15 <- ifelse(distances$dist10.15 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int10.16 <- ifelse(distances$dist10.16 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int10.17 <- ifelse(distances$dist10.17 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int10.18 <- ifelse(distances$dist10.18 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int10.19 <- ifelse(distances$dist10.19 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int10.20 <- ifelse(distances$dist10.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int11.1 <- ifelse(distances$dist1.11 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int11.2 <- ifelse(distances$dist2.11 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int11.3 <- ifelse(distances$dist3.11 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int11.4 <- ifelse(distances$dist4.11 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int11.5 <- ifelse(distances$dist5.11 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int11.6 <- ifelse(distances$dist6.11 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int11.7 <- ifelse(distances$dist7.11 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int11.8 <- ifelse(distances$dist8.11 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int11.9 <- ifelse(distances$dist9.11 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int11.10 <- ifelse(distances$dist10.11 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int11.12 <- ifelse(distances$dist11.12 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int11.13 <- ifelse(distances$dist11.13 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int11.14 <- ifelse(distances$dist11.14 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int11.15 <- ifelse(distances$dist11.15 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int11.16 <- ifelse(distances$dist11.16 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int11.17 <- ifelse(distances$dist11.17 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int11.18 <- ifelse(distances$dist11.18 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int11.19 <- ifelse(distances$dist11.19 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int11.20 <- ifelse(distances$dist11.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int12.1 <- ifelse(distances$dist1.12 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int12.2 <- ifelse(distances$dist2.12 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int12.3 <- ifelse(distances$dist3.12 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int12.4 <- ifelse(distances$dist4.12 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int12.5 <- ifelse(distances$dist5.12 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int12.6 <- ifelse(distances$dist6.12 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int12.7 <- ifelse(distances$dist7.12 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int12.8 <- ifelse(distances$dist8.12 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int12.9 <- ifelse(distances$dist9.12 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int12.10 <- ifelse(distances$dist10.12 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int12.11 <- ifelse(distances$dist11.12 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int12.13 <- ifelse(distances$dist12.13 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int12.14 <- ifelse(distances$dist12.14 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int12.15 <- ifelse(distances$dist12.15 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int12.16 <- ifelse(distances$dist12.16 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int12.17 <- ifelse(distances$dist12.17 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int12.18 <- ifelse(distances$dist12.18 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int12.19 <- ifelse(distances$dist12.19 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int12.20 <- ifelse(distances$dist12.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int13.1 <- ifelse(distances$dist1.13 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int13.2 <- ifelse(distances$dist2.13 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int13.3 <- ifelse(distances$dist3.13 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int13.4 <- ifelse(distances$dist4.13 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int13.5 <- ifelse(distances$dist5.13 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int13.6 <- ifelse(distances$dist6.13 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int13.7 <- ifelse(distances$dist7.13 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int13.8 <- ifelse(distances$dist8.13 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int13.9 <- ifelse(distances$dist9.13 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int13.10 <- ifelse(distances$dist10.13 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int13.11 <- ifelse(distances$dist11.13 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int13.12 <- ifelse(distances$dist12.13 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int13.14 <- ifelse(distances$dist13.14 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int13.15 <- ifelse(distances$dist13.15 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int13.16 <- ifelse(distances$dist13.16 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int13.17 <- ifelse(distances$dist13.17 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int13.18 <- ifelse(distances$dist13.18 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int13.19 <- ifelse(distances$dist13.19 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int13.20 <- ifelse(distances$dist13.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int14.1 <- ifelse(distances$dist1.14 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int14.2 <- ifelse(distances$dist2.14 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int14.3 <- ifelse(distances$dist3.14 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int14.4 <- ifelse(distances$dist4.14 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int14.5 <- ifelse(distances$dist5.14 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int14.6 <- ifelse(distances$dist6.14 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int14.7 <- ifelse(distances$dist7.14 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int14.8 <- ifelse(distances$dist8.14 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int14.9 <- ifelse(distances$dist9.14 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int14.10 <- ifelse(distances$dist10.14 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int14.11 <- ifelse(distances$dist11.14 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int14.12 <- ifelse(distances$dist12.14 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int14.13 <- ifelse(distances$dist13.14 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int14.15 <- ifelse(distances$dist14.15 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int14.16 <- ifelse(distances$dist14.16 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int14.17 <- ifelse(distances$dist14.17 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int14.18 <- ifelse(distances$dist14.18 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int14.19 <- ifelse(distances$dist14.19 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int14.20 <- ifelse(distances$dist14.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int15.1 <- ifelse(distances$dist1.15 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int15.2 <- ifelse(distances$dist2.15 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int15.3 <- ifelse(distances$dist3.15 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int15.4 <- ifelse(distances$dist4.15 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int15.5 <- ifelse(distances$dist5.15 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int15.6 <- ifelse(distances$dist6.15 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int15.7 <- ifelse(distances$dist7.15 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int15.8 <- ifelse(distances$dist8.15 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int15.9 <- ifelse(distances$dist9.15 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int15.10 <- ifelse(distances$dist10.15 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int15.11 <- ifelse(distances$dist11.15 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int15.12 <- ifelse(distances$dist12.15 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int15.13 <- ifelse(distances$dist13.15 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int15.14 <- ifelse(distances$dist14.15 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int15.16 <- ifelse(distances$dist15.16 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int15.17 <- ifelse(distances$dist15.17 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int15.18 <- ifelse(distances$dist15.18 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int15.19 <- ifelse(distances$dist15.19 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int15.20 <- ifelse(distances$dist15.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int16.1 <- ifelse(distances$dist1.16 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int16.2 <- ifelse(distances$dist2.16 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int16.3 <- ifelse(distances$dist3.16 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int16.4 <- ifelse(distances$dist4.16 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int16.5 <- ifelse(distances$dist5.16 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int16.6 <- ifelse(distances$dist6.16 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int16.7 <- ifelse(distances$dist7.16 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int16.8 <- ifelse(distances$dist8.16 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int16.9 <- ifelse(distances$dist9.16 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int16.10 <- ifelse(distances$dist10.16 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int16.11 <- ifelse(distances$dist11.16 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int16.12 <- ifelse(distances$dist12.16 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int16.13 <- ifelse(distances$dist13.16 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int16.14 <- ifelse(distances$dist14.16 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int16.15 <- ifelse(distances$dist15.16 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int16.17 <- ifelse(distances$dist16.17 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int16.18 <- ifelse(distances$dist16.18 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int16.19 <- ifelse(distances$dist16.19 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int16.20 <- ifelse(distances$dist16.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int17.1 <- ifelse(distances$dist1.17 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int17.2 <- ifelse(distances$dist2.17 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int17.3 <- ifelse(distances$dist3.17 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int17.4 <- ifelse(distances$dist4.17 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int17.5 <- ifelse(distances$dist5.17 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int17.6 <- ifelse(distances$dist6.17 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int17.7 <- ifelse(distances$dist7.17 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int17.8 <- ifelse(distances$dist8.17 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int17.9 <- ifelse(distances$dist9.17 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int17.10 <- ifelse(distances$dist10.17 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int17.11 <- ifelse(distances$dist11.17 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int17.12 <- ifelse(distances$dist12.17 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int17.13 <- ifelse(distances$dist13.17 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int17.14 <- ifelse(distances$dist14.17 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int17.15 <- ifelse(distances$dist15.17 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int17.16 <- ifelse(distances$dist16.17 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int17.18 <- ifelse(distances$dist17.18 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int17.19 <- ifelse(distances$dist17.19 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int17.20 <- ifelse(distances$dist17.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int18.1 <- ifelse(distances$dist1.18 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int18.2 <- ifelse(distances$dist2.18 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int18.3 <- ifelse(distances$dist3.18 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int18.4 <- ifelse(distances$dist4.18 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int18.5 <- ifelse(distances$dist5.18 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int18.6 <- ifelse(distances$dist6.18 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int18.7 <- ifelse(distances$dist7.18 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int18.8 <- ifelse(distances$dist8.18 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int18.9 <- ifelse(distances$dist9.18 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int18.10 <- ifelse(distances$dist10.18 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int18.11 <- ifelse(distances$dist11.18 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int18.12 <- ifelse(distances$dist12.18 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int18.13 <- ifelse(distances$dist13.18 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int18.14 <- ifelse(distances$dist14.18 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int18.15 <- ifelse(distances$dist15.18 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int18.16 <- ifelse(distances$dist16.18 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int18.17 <- ifelse(distances$dist17.18 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int18.19 <- ifelse(distances$dist18.19 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int18.20 <- ifelse(distances$dist18.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int19.1 <- ifelse(distances$dist1.19 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int19.2 <- ifelse(distances$dist2.19 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int19.3 <- ifelse(distances$dist3.19 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int19.4 <- ifelse(distances$dist4.19 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int19.5 <- ifelse(distances$dist5.19 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int19.6 <- ifelse(distances$dist6.19 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int19.7 <- ifelse(distances$dist7.19 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int19.8 <- ifelse(distances$dist8.19 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int19.9 <- ifelse(distances$dist9.19 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int19.10 <- ifelse(distances$dist10.19 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int19.11 <- ifelse(distances$dist11.19 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int19.12 <- ifelse(distances$dist12.19 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int19.13 <- ifelse(distances$dist13.19 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int19.14 <- ifelse(distances$dist14.19 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int19.15 <- ifelse(distances$dist15.19 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int19.16 <- ifelse(distances$dist16.19 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int19.17 <- ifelse(distances$dist17.19 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int19.18 <- ifelse(distances$dist18.19 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int19.20 <- ifelse(distances$dist19.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int20.1 <- ifelse(distances$dist1.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int20.2 <- ifelse(distances$dist2.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int20.3 <- ifelse(distances$dist3.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int20.4 <- ifelse(distances$dist4.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int20.5 <- ifelse(distances$dist5.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int20.6 <- ifelse(distances$dist6.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int20.7 <- ifelse(distances$dist7.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int20.8 <- ifelse(distances$dist8.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int20.9 <- ifelse(distances$dist9.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int20.10 <- ifelse(distances$dist10.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int20.11 <- ifelse(distances$dist11.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int20.12 <- ifelse(distances$dist12.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int20.13 <- ifelse(distances$dist13.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int20.14 <- ifelse(distances$dist14.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int20.15 <- ifelse(distances$dist15.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int20.16 <- ifelse(distances$dist16.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int20.17 <- ifelse(distances$dist17.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int20.18 <- ifelse(distances$dist18.20 < (2.5*meansize), c(1), c(0))
  distanceinteraction$int20.19 <- ifelse(distances$dist19.20 < (2.5*meansize), c(1), c(0))
  
  # CONVERT FLY ORIENTATIONS TO RADIANS # -------------------------------------
  fly1$oriradians <- ifelse(fly1$ori >= 0, c((fly1$ori*-1)+(2*pi)), c(fly1$ori*-1))
  fly2$oriradians <- ifelse(fly2$ori >= 0, c((fly2$ori*-1)+(2*pi)), c(fly2$ori*-1))
  fly3$oriradians <- ifelse(fly3$ori >= 0, c((fly3$ori*-1)+(2*pi)), c(fly3$ori*-1))
  fly4$oriradians <- ifelse(fly4$ori >= 0, c((fly4$ori*-1)+(2*pi)), c(fly4$ori*-1))
  fly5$oriradians <- ifelse(fly5$ori >= 0, c((fly5$ori*-1)+(2*pi)), c(fly5$ori*-1))
  fly6$oriradians <- ifelse(fly6$ori >= 0, c((fly6$ori*-1)+(2*pi)), c(fly6$ori*-1))
  fly7$oriradians <- ifelse(fly7$ori >= 0, c((fly7$ori*-1)+(2*pi)), c(fly7$ori*-1))
  fly8$oriradians <- ifelse(fly8$ori >= 0, c((fly8$ori*-1)+(2*pi)), c(fly8$ori*-1))
  fly9$oriradians <- ifelse(fly9$ori >= 0, c((fly9$ori*-1)+(2*pi)), c(fly9$ori*-1))
  fly10$oriradians <- ifelse(fly10$ori >= 0, c((fly10$ori*-1)+(2*pi)), c(fly10$ori*-1))
  fly11$oriradians <- ifelse(fly11$ori >= 0, c((fly11$ori*-1)+(2*pi)), c(fly11$ori*-1))
  fly12$oriradians <- ifelse(fly12$ori >= 0, c((fly12$ori*-1)+(2*pi)), c(fly12$ori*-1))
  fly13$oriradians <- ifelse(fly13$ori >= 0, c((fly13$ori*-1)+(2*pi)), c(fly13$ori*-1))
  fly14$oriradians <- ifelse(fly14$ori >= 0, c((fly14$ori*-1)+(2*pi)), c(fly14$ori*-1))
  fly15$oriradians <- ifelse(fly15$ori >= 0, c((fly15$ori*-1)+(2*pi)), c(fly15$ori*-1))
  fly16$oriradians <- ifelse(fly16$ori >= 0, c((fly16$ori*-1)+(2*pi)), c(fly16$ori*-1))
  fly17$oriradians <- ifelse(fly17$ori >= 0, c((fly17$ori*-1)+(2*pi)), c(fly17$ori*-1))
  fly18$oriradians <- ifelse(fly18$ori >= 0, c((fly18$ori*-1)+(2*pi)), c(fly18$ori*-1))
  fly19$oriradians <- ifelse(fly19$ori >= 0, c((fly19$ori*-1)+(2*pi)), c(fly19$ori*-1))
  fly20$oriradians <- ifelse(fly20$ori >= 0, c((fly20$ori*-1)+(2*pi)), c(fly20$ori*-1))
  
  # CREATE A DATAFRAME THAT CALCULATES ARCSIN BETWEEN ALL PAIRWISE COMBINATION OF FLIES #-------------------------------------
  anglearcsin = data.frame(asin((fly2$pos_y - fly1$pos_y)/(sqrt(((fly1$pos_x - fly2$pos_x)^2) + ((fly1$pos_y - fly2$pos_y)^2)))))
  colnames(anglearcsin) <- c("arcsin1.2")
  anglearcsin$arcsin1.3 <- asin((fly3$pos_y - fly1$pos_y)/(sqrt(((fly1$pos_x - fly3$pos_x)^2) + ((fly1$pos_y - fly3$pos_y)^2))))
  anglearcsin$arcsin1.4 <- asin((fly4$pos_y - fly1$pos_y)/(sqrt(((fly1$pos_x - fly4$pos_x)^2) + ((fly1$pos_y - fly4$pos_y)^2))))
  anglearcsin$arcsin1.5 <- asin((fly5$pos_y - fly1$pos_y)/(sqrt(((fly1$pos_x - fly5$pos_x)^2) + ((fly1$pos_y - fly5$pos_y)^2))))
  anglearcsin$arcsin1.6 <- asin((fly6$pos_y - fly1$pos_y)/(sqrt(((fly1$pos_x - fly6$pos_x)^2) + ((fly1$pos_y - fly6$pos_y)^2))))
  anglearcsin$arcsin1.7 <- asin((fly7$pos_y - fly1$pos_y)/(sqrt(((fly1$pos_x - fly7$pos_x)^2) + ((fly1$pos_y - fly7$pos_y)^2))))
  anglearcsin$arcsin1.8 <- asin((fly8$pos_y - fly1$pos_y)/(sqrt(((fly1$pos_x - fly8$pos_x)^2) + ((fly1$pos_y - fly8$pos_y)^2))))
  anglearcsin$arcsin1.9 <- asin((fly9$pos_y - fly1$pos_y)/(sqrt(((fly1$pos_x - fly9$pos_x)^2) + ((fly1$pos_y - fly9$pos_y)^2))))
  anglearcsin$arcsin1.10 <- asin((fly10$pos_y - fly1$pos_y)/(sqrt(((fly1$pos_x - fly10$pos_x)^2) + ((fly1$pos_y - fly10$pos_y)^2))))
  anglearcsin$arcsin1.11 <- asin((fly11$pos_y - fly1$pos_y)/(sqrt(((fly1$pos_x - fly11$pos_x)^2) + ((fly1$pos_y - fly11$pos_y)^2))))
  anglearcsin$arcsin1.12 <- asin((fly12$pos_y - fly1$pos_y)/(sqrt(((fly1$pos_x - fly12$pos_x)^2) + ((fly1$pos_y - fly12$pos_y)^2))))
  anglearcsin$arcsin1.13 <- asin((fly13$pos_y - fly1$pos_y)/(sqrt(((fly1$pos_x - fly13$pos_x)^2) + ((fly1$pos_y - fly13$pos_y)^2))))
  anglearcsin$arcsin1.14 <- asin((fly14$pos_y - fly1$pos_y)/(sqrt(((fly1$pos_x - fly14$pos_x)^2) + ((fly1$pos_y - fly14$pos_y)^2))))
  anglearcsin$arcsin1.15 <- asin((fly15$pos_y - fly1$pos_y)/(sqrt(((fly1$pos_x - fly15$pos_x)^2) + ((fly1$pos_y - fly15$pos_y)^2))))
  anglearcsin$arcsin1.16 <- asin((fly16$pos_y - fly1$pos_y)/(sqrt(((fly1$pos_x - fly16$pos_x)^2) + ((fly1$pos_y - fly16$pos_y)^2))))
  anglearcsin$arcsin1.17 <- asin((fly17$pos_y - fly1$pos_y)/(sqrt(((fly1$pos_x - fly17$pos_x)^2) + ((fly1$pos_y - fly17$pos_y)^2))))
  anglearcsin$arcsin1.18 <- asin((fly18$pos_y - fly1$pos_y)/(sqrt(((fly1$pos_x - fly18$pos_x)^2) + ((fly1$pos_y - fly18$pos_y)^2))))
  anglearcsin$arcsin1.19 <- asin((fly19$pos_y - fly1$pos_y)/(sqrt(((fly1$pos_x - fly19$pos_x)^2) + ((fly1$pos_y - fly19$pos_y)^2))))
  anglearcsin$arcsin1.20 <- asin((fly20$pos_y - fly1$pos_y)/(sqrt(((fly1$pos_x - fly20$pos_x)^2) + ((fly1$pos_y - fly20$pos_y)^2))))
  anglearcsin$arcsin2.1 <- asin((fly1$pos_y - fly2$pos_y)/(sqrt(((fly2$pos_x - fly1$pos_x)^2) + ((fly2$pos_y - fly1$pos_y)^2))))
  anglearcsin$arcsin2.3 <- asin((fly3$pos_y - fly2$pos_y)/(sqrt(((fly2$pos_x - fly3$pos_x)^2) + ((fly2$pos_y - fly3$pos_y)^2))))
  anglearcsin$arcsin2.4 <- asin((fly4$pos_y - fly2$pos_y)/(sqrt(((fly2$pos_x - fly4$pos_x)^2) + ((fly2$pos_y - fly4$pos_y)^2))))
  anglearcsin$arcsin2.5 <- asin((fly5$pos_y - fly2$pos_y)/(sqrt(((fly2$pos_x - fly5$pos_x)^2) + ((fly2$pos_y - fly5$pos_y)^2))))
  anglearcsin$arcsin2.6 <- asin((fly6$pos_y - fly2$pos_y)/(sqrt(((fly2$pos_x - fly6$pos_x)^2) + ((fly2$pos_y - fly6$pos_y)^2))))
  anglearcsin$arcsin2.7 <- asin((fly7$pos_y - fly2$pos_y)/(sqrt(((fly2$pos_x - fly7$pos_x)^2) + ((fly2$pos_y - fly7$pos_y)^2))))
  anglearcsin$arcsin2.8 <- asin((fly8$pos_y - fly2$pos_y)/(sqrt(((fly2$pos_x - fly8$pos_x)^2) + ((fly2$pos_y - fly8$pos_y)^2))))
  anglearcsin$arcsin2.9 <- asin((fly9$pos_y - fly2$pos_y)/(sqrt(((fly2$pos_x - fly9$pos_x)^2) + ((fly2$pos_y - fly9$pos_y)^2))))
  anglearcsin$arcsin2.10 <- asin((fly10$pos_y - fly2$pos_y)/(sqrt(((fly2$pos_x - fly10$pos_x)^2) + ((fly2$pos_y - fly10$pos_y)^2))))
  anglearcsin$arcsin2.11 <- asin((fly11$pos_y - fly2$pos_y)/(sqrt(((fly2$pos_x - fly11$pos_x)^2) + ((fly2$pos_y - fly11$pos_y)^2))))
  anglearcsin$arcsin2.12 <- asin((fly12$pos_y - fly2$pos_y)/(sqrt(((fly2$pos_x - fly12$pos_x)^2) + ((fly2$pos_y - fly12$pos_y)^2))))
  anglearcsin$arcsin2.13 <- asin((fly13$pos_y - fly2$pos_y)/(sqrt(((fly2$pos_x - fly13$pos_x)^2) + ((fly2$pos_y - fly13$pos_y)^2))))
  anglearcsin$arcsin2.14 <- asin((fly14$pos_y - fly2$pos_y)/(sqrt(((fly2$pos_x - fly14$pos_x)^2) + ((fly2$pos_y - fly14$pos_y)^2))))
  anglearcsin$arcsin2.15 <- asin((fly15$pos_y - fly2$pos_y)/(sqrt(((fly2$pos_x - fly15$pos_x)^2) + ((fly2$pos_y - fly15$pos_y)^2))))
  anglearcsin$arcsin2.16 <- asin((fly16$pos_y - fly2$pos_y)/(sqrt(((fly2$pos_x - fly16$pos_x)^2) + ((fly2$pos_y - fly16$pos_y)^2))))
  anglearcsin$arcsin2.17 <- asin((fly17$pos_y - fly2$pos_y)/(sqrt(((fly2$pos_x - fly17$pos_x)^2) + ((fly2$pos_y - fly17$pos_y)^2))))
  anglearcsin$arcsin2.18 <- asin((fly18$pos_y - fly2$pos_y)/(sqrt(((fly2$pos_x - fly18$pos_x)^2) + ((fly2$pos_y - fly18$pos_y)^2))))
  anglearcsin$arcsin2.19 <- asin((fly19$pos_y - fly2$pos_y)/(sqrt(((fly2$pos_x - fly19$pos_x)^2) + ((fly2$pos_y - fly19$pos_y)^2))))
  anglearcsin$arcsin2.20 <- asin((fly20$pos_y - fly2$pos_y)/(sqrt(((fly2$pos_x - fly20$pos_x)^2) + ((fly2$pos_y - fly20$pos_y)^2))))
  anglearcsin$arcsin3.1 <- asin((fly1$pos_y - fly3$pos_y)/(sqrt(((fly3$pos_x - fly1$pos_x)^2) + ((fly3$pos_y - fly1$pos_y)^2))))
  anglearcsin$arcsin3.2 <- asin((fly2$pos_y - fly3$pos_y)/(sqrt(((fly3$pos_x - fly2$pos_x)^2) + ((fly3$pos_y - fly2$pos_y)^2))))
  anglearcsin$arcsin3.4 <- asin((fly4$pos_y - fly3$pos_y)/(sqrt(((fly3$pos_x - fly4$pos_x)^2) + ((fly3$pos_y - fly4$pos_y)^2))))
  anglearcsin$arcsin3.5 <- asin((fly5$pos_y - fly3$pos_y)/(sqrt(((fly3$pos_x - fly5$pos_x)^2) + ((fly3$pos_y - fly5$pos_y)^2))))
  anglearcsin$arcsin3.6 <- asin((fly6$pos_y - fly3$pos_y)/(sqrt(((fly3$pos_x - fly6$pos_x)^2) + ((fly3$pos_y - fly6$pos_y)^2))))
  anglearcsin$arcsin3.7 <- asin((fly7$pos_y - fly3$pos_y)/(sqrt(((fly3$pos_x - fly7$pos_x)^2) + ((fly3$pos_y - fly7$pos_y)^2))))
  anglearcsin$arcsin3.8 <- asin((fly8$pos_y - fly3$pos_y)/(sqrt(((fly3$pos_x - fly8$pos_x)^2) + ((fly3$pos_y - fly8$pos_y)^2))))
  anglearcsin$arcsin3.9 <- asin((fly9$pos_y - fly3$pos_y)/(sqrt(((fly3$pos_x - fly9$pos_x)^2) + ((fly3$pos_y - fly9$pos_y)^2))))
  anglearcsin$arcsin3.10 <- asin((fly10$pos_y - fly3$pos_y)/(sqrt(((fly3$pos_x - fly10$pos_x)^2) + ((fly3$pos_y - fly10$pos_y)^2))))
  anglearcsin$arcsin3.11 <- asin((fly11$pos_y - fly3$pos_y)/(sqrt(((fly3$pos_x - fly11$pos_x)^2) + ((fly3$pos_y - fly11$pos_y)^2))))
  anglearcsin$arcsin3.12 <- asin((fly12$pos_y - fly3$pos_y)/(sqrt(((fly3$pos_x - fly12$pos_x)^2) + ((fly3$pos_y - fly12$pos_y)^2))))
  anglearcsin$arcsin3.13 <- asin((fly13$pos_y - fly3$pos_y)/(sqrt(((fly3$pos_x - fly13$pos_x)^2) + ((fly3$pos_y - fly13$pos_y)^2))))
  anglearcsin$arcsin3.14 <- asin((fly14$pos_y - fly3$pos_y)/(sqrt(((fly3$pos_x - fly14$pos_x)^2) + ((fly3$pos_y - fly14$pos_y)^2))))
  anglearcsin$arcsin3.15 <- asin((fly15$pos_y - fly3$pos_y)/(sqrt(((fly3$pos_x - fly15$pos_x)^2) + ((fly3$pos_y - fly15$pos_y)^2))))
  anglearcsin$arcsin3.16 <- asin((fly16$pos_y - fly3$pos_y)/(sqrt(((fly3$pos_x - fly16$pos_x)^2) + ((fly3$pos_y - fly16$pos_y)^2))))
  anglearcsin$arcsin3.17 <- asin((fly17$pos_y - fly3$pos_y)/(sqrt(((fly3$pos_x - fly17$pos_x)^2) + ((fly3$pos_y - fly17$pos_y)^2))))
  anglearcsin$arcsin3.18 <- asin((fly18$pos_y - fly3$pos_y)/(sqrt(((fly3$pos_x - fly18$pos_x)^2) + ((fly3$pos_y - fly18$pos_y)^2))))
  anglearcsin$arcsin3.19 <- asin((fly19$pos_y - fly3$pos_y)/(sqrt(((fly3$pos_x - fly19$pos_x)^2) + ((fly3$pos_y - fly19$pos_y)^2))))
  anglearcsin$arcsin3.20 <- asin((fly20$pos_y - fly3$pos_y)/(sqrt(((fly3$pos_x - fly20$pos_x)^2) + ((fly3$pos_y - fly20$pos_y)^2))))
  anglearcsin$arcsin4.1 <- asin((fly1$pos_y - fly4$pos_y)/(sqrt(((fly4$pos_x - fly1$pos_x)^2) + ((fly4$pos_y - fly1$pos_y)^2))))
  anglearcsin$arcsin4.2 <- asin((fly2$pos_y - fly4$pos_y)/(sqrt(((fly4$pos_x - fly2$pos_x)^2) + ((fly4$pos_y - fly2$pos_y)^2))))
  anglearcsin$arcsin4.3 <- asin((fly3$pos_y - fly4$pos_y)/(sqrt(((fly4$pos_x - fly3$pos_x)^2) + ((fly4$pos_y - fly3$pos_y)^2))))
  anglearcsin$arcsin4.5 <- asin((fly5$pos_y - fly4$pos_y)/(sqrt(((fly4$pos_x - fly5$pos_x)^2) + ((fly4$pos_y - fly5$pos_y)^2))))
  anglearcsin$arcsin4.6 <- asin((fly6$pos_y - fly4$pos_y)/(sqrt(((fly4$pos_x - fly6$pos_x)^2) + ((fly4$pos_y - fly6$pos_y)^2))))
  anglearcsin$arcsin4.7 <- asin((fly7$pos_y - fly4$pos_y)/(sqrt(((fly4$pos_x - fly7$pos_x)^2) + ((fly4$pos_y - fly7$pos_y)^2))))
  anglearcsin$arcsin4.8 <- asin((fly8$pos_y - fly4$pos_y)/(sqrt(((fly4$pos_x - fly8$pos_x)^2) + ((fly4$pos_y - fly8$pos_y)^2))))
  anglearcsin$arcsin4.9 <- asin((fly9$pos_y - fly4$pos_y)/(sqrt(((fly4$pos_x - fly9$pos_x)^2) + ((fly4$pos_y - fly9$pos_y)^2))))
  anglearcsin$arcsin4.10 <- asin((fly10$pos_y - fly4$pos_y)/(sqrt(((fly4$pos_x - fly10$pos_x)^2) + ((fly4$pos_y - fly10$pos_y)^2))))
  anglearcsin$arcsin4.11 <- asin((fly11$pos_y - fly4$pos_y)/(sqrt(((fly4$pos_x - fly11$pos_x)^2) + ((fly4$pos_y - fly11$pos_y)^2))))
  anglearcsin$arcsin4.12 <- asin((fly12$pos_y - fly4$pos_y)/(sqrt(((fly4$pos_x - fly12$pos_x)^2) + ((fly4$pos_y - fly12$pos_y)^2))))
  anglearcsin$arcsin4.13 <- asin((fly13$pos_y - fly4$pos_y)/(sqrt(((fly4$pos_x - fly13$pos_x)^2) + ((fly4$pos_y - fly13$pos_y)^2))))
  anglearcsin$arcsin4.14 <- asin((fly14$pos_y - fly4$pos_y)/(sqrt(((fly4$pos_x - fly14$pos_x)^2) + ((fly4$pos_y - fly14$pos_y)^2))))
  anglearcsin$arcsin4.15 <- asin((fly15$pos_y - fly4$pos_y)/(sqrt(((fly4$pos_x - fly15$pos_x)^2) + ((fly4$pos_y - fly15$pos_y)^2))))
  anglearcsin$arcsin4.16 <- asin((fly16$pos_y - fly4$pos_y)/(sqrt(((fly4$pos_x - fly16$pos_x)^2) + ((fly4$pos_y - fly16$pos_y)^2))))
  anglearcsin$arcsin4.17 <- asin((fly17$pos_y - fly4$pos_y)/(sqrt(((fly4$pos_x - fly17$pos_x)^2) + ((fly4$pos_y - fly17$pos_y)^2))))
  anglearcsin$arcsin4.18 <- asin((fly18$pos_y - fly4$pos_y)/(sqrt(((fly4$pos_x - fly18$pos_x)^2) + ((fly4$pos_y - fly18$pos_y)^2))))
  anglearcsin$arcsin4.19 <- asin((fly19$pos_y - fly4$pos_y)/(sqrt(((fly4$pos_x - fly19$pos_x)^2) + ((fly4$pos_y - fly19$pos_y)^2))))
  anglearcsin$arcsin4.20 <- asin((fly20$pos_y - fly4$pos_y)/(sqrt(((fly4$pos_x - fly20$pos_x)^2) + ((fly4$pos_y - fly20$pos_y)^2))))
  anglearcsin$arcsin5.1 <- asin((fly1$pos_y - fly5$pos_y)/(sqrt(((fly5$pos_x - fly1$pos_x)^2) + ((fly5$pos_y - fly1$pos_y)^2))))
  anglearcsin$arcsin5.2 <- asin((fly2$pos_y - fly5$pos_y)/(sqrt(((fly5$pos_x - fly2$pos_x)^2) + ((fly5$pos_y - fly2$pos_y)^2))))
  anglearcsin$arcsin5.3 <- asin((fly3$pos_y - fly5$pos_y)/(sqrt(((fly5$pos_x - fly3$pos_x)^2) + ((fly5$pos_y - fly3$pos_y)^2))))
  anglearcsin$arcsin5.4 <- asin((fly4$pos_y - fly5$pos_y)/(sqrt(((fly5$pos_x - fly4$pos_x)^2) + ((fly5$pos_y - fly4$pos_y)^2))))
  anglearcsin$arcsin5.6 <- asin((fly6$pos_y - fly5$pos_y)/(sqrt(((fly5$pos_x - fly6$pos_x)^2) + ((fly5$pos_y - fly6$pos_y)^2))))
  anglearcsin$arcsin5.7 <- asin((fly7$pos_y - fly5$pos_y)/(sqrt(((fly5$pos_x - fly7$pos_x)^2) + ((fly5$pos_y - fly7$pos_y)^2))))
  anglearcsin$arcsin5.8 <- asin((fly8$pos_y - fly5$pos_y)/(sqrt(((fly5$pos_x - fly8$pos_x)^2) + ((fly5$pos_y - fly8$pos_y)^2))))
  anglearcsin$arcsin5.9 <- asin((fly9$pos_y - fly5$pos_y)/(sqrt(((fly5$pos_x - fly9$pos_x)^2) + ((fly5$pos_y - fly9$pos_y)^2))))
  anglearcsin$arcsin5.10 <- asin((fly10$pos_y - fly5$pos_y)/(sqrt(((fly5$pos_x - fly10$pos_x)^2) + ((fly5$pos_y - fly10$pos_y)^2))))
  anglearcsin$arcsin5.11 <- asin((fly11$pos_y - fly5$pos_y)/(sqrt(((fly5$pos_x - fly11$pos_x)^2) + ((fly5$pos_y - fly11$pos_y)^2))))
  anglearcsin$arcsin5.12 <- asin((fly12$pos_y - fly5$pos_y)/(sqrt(((fly5$pos_x - fly12$pos_x)^2) + ((fly5$pos_y - fly12$pos_y)^2))))
  anglearcsin$arcsin5.13 <- asin((fly13$pos_y - fly5$pos_y)/(sqrt(((fly5$pos_x - fly13$pos_x)^2) + ((fly5$pos_y - fly13$pos_y)^2))))
  anglearcsin$arcsin5.14 <- asin((fly14$pos_y - fly5$pos_y)/(sqrt(((fly5$pos_x - fly14$pos_x)^2) + ((fly5$pos_y - fly14$pos_y)^2))))
  anglearcsin$arcsin5.15 <- asin((fly15$pos_y - fly5$pos_y)/(sqrt(((fly5$pos_x - fly15$pos_x)^2) + ((fly5$pos_y - fly15$pos_y)^2))))
  anglearcsin$arcsin5.16 <- asin((fly16$pos_y - fly5$pos_y)/(sqrt(((fly5$pos_x - fly16$pos_x)^2) + ((fly5$pos_y - fly16$pos_y)^2))))
  anglearcsin$arcsin5.17 <- asin((fly17$pos_y - fly5$pos_y)/(sqrt(((fly5$pos_x - fly17$pos_x)^2) + ((fly5$pos_y - fly17$pos_y)^2))))
  anglearcsin$arcsin5.18 <- asin((fly18$pos_y - fly5$pos_y)/(sqrt(((fly5$pos_x - fly18$pos_x)^2) + ((fly5$pos_y - fly18$pos_y)^2))))
  anglearcsin$arcsin5.19 <- asin((fly19$pos_y - fly5$pos_y)/(sqrt(((fly5$pos_x - fly19$pos_x)^2) + ((fly5$pos_y - fly19$pos_y)^2))))
  anglearcsin$arcsin5.20 <- asin((fly20$pos_y - fly5$pos_y)/(sqrt(((fly5$pos_x - fly20$pos_x)^2) + ((fly5$pos_y - fly20$pos_y)^2))))
  anglearcsin$arcsin6.1 <- asin((fly1$pos_y - fly6$pos_y)/(sqrt(((fly6$pos_x - fly1$pos_x)^2) + ((fly6$pos_y - fly1$pos_y)^2))))
  anglearcsin$arcsin6.2 <- asin((fly2$pos_y - fly6$pos_y)/(sqrt(((fly6$pos_x - fly2$pos_x)^2) + ((fly6$pos_y - fly2$pos_y)^2))))
  anglearcsin$arcsin6.3 <- asin((fly3$pos_y - fly6$pos_y)/(sqrt(((fly6$pos_x - fly3$pos_x)^2) + ((fly6$pos_y - fly3$pos_y)^2))))
  anglearcsin$arcsin6.4 <- asin((fly4$pos_y - fly6$pos_y)/(sqrt(((fly6$pos_x - fly4$pos_x)^2) + ((fly6$pos_y - fly4$pos_y)^2))))
  anglearcsin$arcsin6.5 <- asin((fly5$pos_y - fly6$pos_y)/(sqrt(((fly6$pos_x - fly5$pos_x)^2) + ((fly6$pos_y - fly5$pos_y)^2))))
  anglearcsin$arcsin6.7 <- asin((fly7$pos_y - fly6$pos_y)/(sqrt(((fly6$pos_x - fly7$pos_x)^2) + ((fly6$pos_y - fly7$pos_y)^2))))
  anglearcsin$arcsin6.8 <- asin((fly8$pos_y - fly6$pos_y)/(sqrt(((fly6$pos_x - fly8$pos_x)^2) + ((fly6$pos_y - fly8$pos_y)^2))))
  anglearcsin$arcsin6.9 <- asin((fly9$pos_y - fly6$pos_y)/(sqrt(((fly6$pos_x - fly9$pos_x)^2) + ((fly6$pos_y - fly9$pos_y)^2))))
  anglearcsin$arcsin6.10 <- asin((fly10$pos_y - fly6$pos_y)/(sqrt(((fly6$pos_x - fly10$pos_x)^2) + ((fly6$pos_y - fly10$pos_y)^2))))
  anglearcsin$arcsin6.11 <- asin((fly11$pos_y - fly6$pos_y)/(sqrt(((fly6$pos_x - fly11$pos_x)^2) + ((fly6$pos_y - fly11$pos_y)^2))))
  anglearcsin$arcsin6.12 <- asin((fly12$pos_y - fly6$pos_y)/(sqrt(((fly6$pos_x - fly12$pos_x)^2) + ((fly6$pos_y - fly12$pos_y)^2))))
  anglearcsin$arcsin6.13 <- asin((fly13$pos_y - fly6$pos_y)/(sqrt(((fly6$pos_x - fly13$pos_x)^2) + ((fly6$pos_y - fly13$pos_y)^2))))
  anglearcsin$arcsin6.14 <- asin((fly14$pos_y - fly6$pos_y)/(sqrt(((fly6$pos_x - fly14$pos_x)^2) + ((fly6$pos_y - fly14$pos_y)^2))))
  anglearcsin$arcsin6.15 <- asin((fly15$pos_y - fly6$pos_y)/(sqrt(((fly6$pos_x - fly15$pos_x)^2) + ((fly6$pos_y - fly15$pos_y)^2))))
  anglearcsin$arcsin6.16 <- asin((fly16$pos_y - fly6$pos_y)/(sqrt(((fly6$pos_x - fly16$pos_x)^2) + ((fly6$pos_y - fly16$pos_y)^2))))
  anglearcsin$arcsin6.17 <- asin((fly17$pos_y - fly6$pos_y)/(sqrt(((fly6$pos_x - fly17$pos_x)^2) + ((fly6$pos_y - fly17$pos_y)^2))))
  anglearcsin$arcsin6.18 <- asin((fly18$pos_y - fly6$pos_y)/(sqrt(((fly6$pos_x - fly18$pos_x)^2) + ((fly6$pos_y - fly18$pos_y)^2))))
  anglearcsin$arcsin6.19 <- asin((fly19$pos_y - fly6$pos_y)/(sqrt(((fly6$pos_x - fly19$pos_x)^2) + ((fly6$pos_y - fly19$pos_y)^2))))
  anglearcsin$arcsin6.20 <- asin((fly20$pos_y - fly6$pos_y)/(sqrt(((fly6$pos_x - fly20$pos_x)^2) + ((fly6$pos_y - fly20$pos_y)^2))))
  anglearcsin$arcsin7.1 <- asin((fly1$pos_y - fly7$pos_y)/(sqrt(((fly7$pos_x - fly1$pos_x)^2) + ((fly7$pos_y - fly1$pos_y)^2))))
  anglearcsin$arcsin7.2 <- asin((fly2$pos_y - fly7$pos_y)/(sqrt(((fly7$pos_x - fly2$pos_x)^2) + ((fly7$pos_y - fly2$pos_y)^2))))
  anglearcsin$arcsin7.3 <- asin((fly3$pos_y - fly7$pos_y)/(sqrt(((fly7$pos_x - fly3$pos_x)^2) + ((fly7$pos_y - fly3$pos_y)^2))))
  anglearcsin$arcsin7.4 <- asin((fly4$pos_y - fly7$pos_y)/(sqrt(((fly7$pos_x - fly4$pos_x)^2) + ((fly7$pos_y - fly4$pos_y)^2))))
  anglearcsin$arcsin7.5 <- asin((fly5$pos_y - fly7$pos_y)/(sqrt(((fly7$pos_x - fly5$pos_x)^2) + ((fly7$pos_y - fly5$pos_y)^2))))
  anglearcsin$arcsin7.6 <- asin((fly6$pos_y - fly7$pos_y)/(sqrt(((fly7$pos_x - fly6$pos_x)^2) + ((fly7$pos_y - fly6$pos_y)^2))))
  anglearcsin$arcsin7.8 <- asin((fly8$pos_y - fly7$pos_y)/(sqrt(((fly7$pos_x - fly8$pos_x)^2) + ((fly7$pos_y - fly8$pos_y)^2))))
  anglearcsin$arcsin7.9 <- asin((fly9$pos_y - fly7$pos_y)/(sqrt(((fly7$pos_x - fly9$pos_x)^2) + ((fly7$pos_y - fly9$pos_y)^2))))
  anglearcsin$arcsin7.10 <- asin((fly10$pos_y - fly7$pos_y)/(sqrt(((fly7$pos_x - fly10$pos_x)^2) + ((fly7$pos_y - fly10$pos_y)^2))))
  anglearcsin$arcsin7.11 <- asin((fly11$pos_y - fly7$pos_y)/(sqrt(((fly7$pos_x - fly11$pos_x)^2) + ((fly7$pos_y - fly11$pos_y)^2))))
  anglearcsin$arcsin7.12 <- asin((fly12$pos_y - fly7$pos_y)/(sqrt(((fly7$pos_x - fly12$pos_x)^2) + ((fly7$pos_y - fly12$pos_y)^2))))
  anglearcsin$arcsin7.13 <- asin((fly13$pos_y - fly7$pos_y)/(sqrt(((fly7$pos_x - fly13$pos_x)^2) + ((fly7$pos_y - fly13$pos_y)^2))))
  anglearcsin$arcsin7.14 <- asin((fly14$pos_y - fly7$pos_y)/(sqrt(((fly7$pos_x - fly14$pos_x)^2) + ((fly7$pos_y - fly14$pos_y)^2))))
  anglearcsin$arcsin7.15 <- asin((fly15$pos_y - fly7$pos_y)/(sqrt(((fly7$pos_x - fly15$pos_x)^2) + ((fly7$pos_y - fly15$pos_y)^2))))
  anglearcsin$arcsin7.16 <- asin((fly16$pos_y - fly7$pos_y)/(sqrt(((fly7$pos_x - fly16$pos_x)^2) + ((fly7$pos_y - fly16$pos_y)^2))))
  anglearcsin$arcsin7.17 <- asin((fly17$pos_y - fly7$pos_y)/(sqrt(((fly7$pos_x - fly17$pos_x)^2) + ((fly7$pos_y - fly17$pos_y)^2))))
  anglearcsin$arcsin7.18 <- asin((fly18$pos_y - fly7$pos_y)/(sqrt(((fly7$pos_x - fly18$pos_x)^2) + ((fly7$pos_y - fly18$pos_y)^2))))
  anglearcsin$arcsin7.19 <- asin((fly19$pos_y - fly7$pos_y)/(sqrt(((fly7$pos_x - fly19$pos_x)^2) + ((fly7$pos_y - fly19$pos_y)^2))))
  anglearcsin$arcsin7.20 <- asin((fly20$pos_y - fly7$pos_y)/(sqrt(((fly7$pos_x - fly20$pos_x)^2) + ((fly7$pos_y - fly20$pos_y)^2))))
  anglearcsin$arcsin8.1 <- asin((fly1$pos_y - fly8$pos_y)/(sqrt(((fly8$pos_x - fly1$pos_x)^2) + ((fly8$pos_y - fly1$pos_y)^2))))
  anglearcsin$arcsin8.2 <- asin((fly2$pos_y - fly8$pos_y)/(sqrt(((fly8$pos_x - fly2$pos_x)^2) + ((fly8$pos_y - fly2$pos_y)^2))))
  anglearcsin$arcsin8.3 <- asin((fly3$pos_y - fly8$pos_y)/(sqrt(((fly8$pos_x - fly3$pos_x)^2) + ((fly8$pos_y - fly3$pos_y)^2))))
  anglearcsin$arcsin8.4 <- asin((fly4$pos_y - fly8$pos_y)/(sqrt(((fly8$pos_x - fly4$pos_x)^2) + ((fly8$pos_y - fly4$pos_y)^2))))
  anglearcsin$arcsin8.5 <- asin((fly5$pos_y - fly8$pos_y)/(sqrt(((fly8$pos_x - fly5$pos_x)^2) + ((fly8$pos_y - fly5$pos_y)^2))))
  anglearcsin$arcsin8.6 <- asin((fly6$pos_y - fly8$pos_y)/(sqrt(((fly8$pos_x - fly6$pos_x)^2) + ((fly8$pos_y - fly6$pos_y)^2))))
  anglearcsin$arcsin8.7 <- asin((fly7$pos_y - fly8$pos_y)/(sqrt(((fly8$pos_x - fly7$pos_x)^2) + ((fly8$pos_y - fly7$pos_y)^2))))
  anglearcsin$arcsin8.9 <- asin((fly9$pos_y - fly8$pos_y)/(sqrt(((fly8$pos_x - fly9$pos_x)^2) + ((fly8$pos_y - fly9$pos_y)^2))))
  anglearcsin$arcsin8.10 <- asin((fly10$pos_y - fly8$pos_y)/(sqrt(((fly8$pos_x - fly10$pos_x)^2) + ((fly8$pos_y - fly10$pos_y)^2))))
  anglearcsin$arcsin8.11 <- asin((fly11$pos_y - fly8$pos_y)/(sqrt(((fly8$pos_x - fly11$pos_x)^2) + ((fly8$pos_y - fly11$pos_y)^2))))
  anglearcsin$arcsin8.12 <- asin((fly12$pos_y - fly8$pos_y)/(sqrt(((fly8$pos_x - fly12$pos_x)^2) + ((fly8$pos_y - fly12$pos_y)^2))))
  anglearcsin$arcsin8.13 <- asin((fly13$pos_y - fly8$pos_y)/(sqrt(((fly8$pos_x - fly13$pos_x)^2) + ((fly8$pos_y - fly13$pos_y)^2))))
  anglearcsin$arcsin8.14 <- asin((fly14$pos_y - fly8$pos_y)/(sqrt(((fly8$pos_x - fly14$pos_x)^2) + ((fly8$pos_y - fly14$pos_y)^2))))
  anglearcsin$arcsin8.15 <- asin((fly15$pos_y - fly8$pos_y)/(sqrt(((fly8$pos_x - fly15$pos_x)^2) + ((fly8$pos_y - fly15$pos_y)^2))))
  anglearcsin$arcsin8.16 <- asin((fly16$pos_y - fly8$pos_y)/(sqrt(((fly8$pos_x - fly16$pos_x)^2) + ((fly8$pos_y - fly16$pos_y)^2))))
  anglearcsin$arcsin8.17 <- asin((fly17$pos_y - fly8$pos_y)/(sqrt(((fly8$pos_x - fly17$pos_x)^2) + ((fly8$pos_y - fly17$pos_y)^2))))
  anglearcsin$arcsin8.18 <- asin((fly18$pos_y - fly8$pos_y)/(sqrt(((fly8$pos_x - fly18$pos_x)^2) + ((fly8$pos_y - fly18$pos_y)^2))))
  anglearcsin$arcsin8.19 <- asin((fly19$pos_y - fly8$pos_y)/(sqrt(((fly8$pos_x - fly19$pos_x)^2) + ((fly8$pos_y - fly19$pos_y)^2))))
  anglearcsin$arcsin8.20 <- asin((fly20$pos_y - fly8$pos_y)/(sqrt(((fly8$pos_x - fly20$pos_x)^2) + ((fly8$pos_y - fly20$pos_y)^2))))
  anglearcsin$arcsin9.1 <- asin((fly1$pos_y - fly9$pos_y)/(sqrt(((fly9$pos_x - fly1$pos_x)^2) + ((fly9$pos_y - fly1$pos_y)^2))))
  anglearcsin$arcsin9.2 <- asin((fly2$pos_y - fly9$pos_y)/(sqrt(((fly9$pos_x - fly2$pos_x)^2) + ((fly9$pos_y - fly2$pos_y)^2))))
  anglearcsin$arcsin9.3 <- asin((fly3$pos_y - fly9$pos_y)/(sqrt(((fly9$pos_x - fly3$pos_x)^2) + ((fly9$pos_y - fly3$pos_y)^2))))
  anglearcsin$arcsin9.4 <- asin((fly4$pos_y - fly9$pos_y)/(sqrt(((fly9$pos_x - fly4$pos_x)^2) + ((fly9$pos_y - fly4$pos_y)^2))))
  anglearcsin$arcsin9.5 <- asin((fly5$pos_y - fly9$pos_y)/(sqrt(((fly9$pos_x - fly5$pos_x)^2) + ((fly9$pos_y - fly5$pos_y)^2))))
  anglearcsin$arcsin9.6 <- asin((fly6$pos_y - fly9$pos_y)/(sqrt(((fly9$pos_x - fly6$pos_x)^2) + ((fly9$pos_y - fly6$pos_y)^2))))
  anglearcsin$arcsin9.7 <- asin((fly7$pos_y - fly9$pos_y)/(sqrt(((fly9$pos_x - fly7$pos_x)^2) + ((fly9$pos_y - fly7$pos_y)^2))))
  anglearcsin$arcsin9.8 <- asin((fly8$pos_y - fly9$pos_y)/(sqrt(((fly9$pos_x - fly8$pos_x)^2) + ((fly9$pos_y - fly8$pos_y)^2))))
  anglearcsin$arcsin9.10 <- asin((fly10$pos_y - fly9$pos_y)/(sqrt(((fly9$pos_x - fly10$pos_x)^2) + ((fly9$pos_y - fly10$pos_y)^2))))
  anglearcsin$arcsin9.11 <- asin((fly11$pos_y - fly9$pos_y)/(sqrt(((fly9$pos_x - fly11$pos_x)^2) + ((fly9$pos_y - fly11$pos_y)^2))))
  anglearcsin$arcsin9.12 <- asin((fly12$pos_y - fly9$pos_y)/(sqrt(((fly9$pos_x - fly12$pos_x)^2) + ((fly9$pos_y - fly12$pos_y)^2))))
  anglearcsin$arcsin9.13 <- asin((fly13$pos_y - fly9$pos_y)/(sqrt(((fly9$pos_x - fly13$pos_x)^2) + ((fly9$pos_y - fly13$pos_y)^2))))
  anglearcsin$arcsin9.14 <- asin((fly14$pos_y - fly9$pos_y)/(sqrt(((fly9$pos_x - fly14$pos_x)^2) + ((fly9$pos_y - fly14$pos_y)^2))))
  anglearcsin$arcsin9.15 <- asin((fly15$pos_y - fly9$pos_y)/(sqrt(((fly9$pos_x - fly15$pos_x)^2) + ((fly9$pos_y - fly15$pos_y)^2))))
  anglearcsin$arcsin9.16 <- asin((fly16$pos_y - fly9$pos_y)/(sqrt(((fly9$pos_x - fly16$pos_x)^2) + ((fly9$pos_y - fly16$pos_y)^2))))
  anglearcsin$arcsin9.17 <- asin((fly17$pos_y - fly9$pos_y)/(sqrt(((fly9$pos_x - fly17$pos_x)^2) + ((fly9$pos_y - fly17$pos_y)^2))))
  anglearcsin$arcsin9.18 <- asin((fly18$pos_y - fly9$pos_y)/(sqrt(((fly9$pos_x - fly18$pos_x)^2) + ((fly9$pos_y - fly18$pos_y)^2))))
  anglearcsin$arcsin9.19 <- asin((fly19$pos_y - fly9$pos_y)/(sqrt(((fly9$pos_x - fly19$pos_x)^2) + ((fly9$pos_y - fly19$pos_y)^2))))
  anglearcsin$arcsin9.20 <- asin((fly20$pos_y - fly9$pos_y)/(sqrt(((fly9$pos_x - fly20$pos_x)^2) + ((fly9$pos_y - fly20$pos_y)^2))))
  anglearcsin$arcsin10.1 <- asin((fly1$pos_y - fly10$pos_y)/(sqrt(((fly10$pos_x - fly1$pos_x)^2) + ((fly10$pos_y - fly1$pos_y)^2))))
  anglearcsin$arcsin10.2 <- asin((fly2$pos_y - fly10$pos_y)/(sqrt(((fly10$pos_x - fly2$pos_x)^2) + ((fly10$pos_y - fly2$pos_y)^2))))
  anglearcsin$arcsin10.3 <- asin((fly3$pos_y - fly10$pos_y)/(sqrt(((fly10$pos_x - fly3$pos_x)^2) + ((fly10$pos_y - fly3$pos_y)^2))))
  anglearcsin$arcsin10.4 <- asin((fly4$pos_y - fly10$pos_y)/(sqrt(((fly10$pos_x - fly4$pos_x)^2) + ((fly10$pos_y - fly4$pos_y)^2))))
  anglearcsin$arcsin10.5 <- asin((fly5$pos_y - fly10$pos_y)/(sqrt(((fly10$pos_x - fly5$pos_x)^2) + ((fly10$pos_y - fly5$pos_y)^2))))
  anglearcsin$arcsin10.6 <- asin((fly6$pos_y - fly10$pos_y)/(sqrt(((fly10$pos_x - fly6$pos_x)^2) + ((fly10$pos_y - fly6$pos_y)^2))))
  anglearcsin$arcsin10.7 <- asin((fly7$pos_y - fly10$pos_y)/(sqrt(((fly10$pos_x - fly7$pos_x)^2) + ((fly10$pos_y - fly7$pos_y)^2))))
  anglearcsin$arcsin10.8 <- asin((fly8$pos_y - fly10$pos_y)/(sqrt(((fly10$pos_x - fly8$pos_x)^2) + ((fly10$pos_y - fly8$pos_y)^2))))
  anglearcsin$arcsin10.9 <- asin((fly9$pos_y - fly10$pos_y)/(sqrt(((fly10$pos_x - fly9$pos_x)^2) + ((fly10$pos_y - fly9$pos_y)^2))))
  anglearcsin$arcsin10.11 <- asin((fly11$pos_y - fly10$pos_y)/(sqrt(((fly10$pos_x - fly11$pos_x)^2) + ((fly10$pos_y - fly11$pos_y)^2))))
  anglearcsin$arcsin10.12 <- asin((fly12$pos_y - fly10$pos_y)/(sqrt(((fly10$pos_x - fly12$pos_x)^2) + ((fly10$pos_y - fly12$pos_y)^2))))
  anglearcsin$arcsin10.13 <- asin((fly13$pos_y - fly10$pos_y)/(sqrt(((fly10$pos_x - fly13$pos_x)^2) + ((fly10$pos_y - fly13$pos_y)^2))))
  anglearcsin$arcsin10.14 <- asin((fly14$pos_y - fly10$pos_y)/(sqrt(((fly10$pos_x - fly14$pos_x)^2) + ((fly10$pos_y - fly14$pos_y)^2))))
  anglearcsin$arcsin10.15 <- asin((fly15$pos_y - fly10$pos_y)/(sqrt(((fly10$pos_x - fly15$pos_x)^2) + ((fly10$pos_y - fly15$pos_y)^2))))
  anglearcsin$arcsin10.16 <- asin((fly16$pos_y - fly10$pos_y)/(sqrt(((fly10$pos_x - fly16$pos_x)^2) + ((fly10$pos_y - fly16$pos_y)^2))))
  anglearcsin$arcsin10.17 <- asin((fly17$pos_y - fly10$pos_y)/(sqrt(((fly10$pos_x - fly17$pos_x)^2) + ((fly10$pos_y - fly17$pos_y)^2))))
  anglearcsin$arcsin10.18 <- asin((fly18$pos_y - fly10$pos_y)/(sqrt(((fly10$pos_x - fly18$pos_x)^2) + ((fly10$pos_y - fly18$pos_y)^2))))
  anglearcsin$arcsin10.19 <- asin((fly19$pos_y - fly10$pos_y)/(sqrt(((fly10$pos_x - fly19$pos_x)^2) + ((fly10$pos_y - fly19$pos_y)^2))))
  anglearcsin$arcsin10.20 <- asin((fly20$pos_y - fly10$pos_y)/(sqrt(((fly10$pos_x - fly20$pos_x)^2) + ((fly10$pos_y - fly20$pos_y)^2))))
  anglearcsin$arcsin11.1 <- asin((fly1$pos_y - fly11$pos_y)/(sqrt(((fly11$pos_x - fly1$pos_x)^2) + ((fly11$pos_y - fly1$pos_y)^2))))
  anglearcsin$arcsin11.2 <- asin((fly2$pos_y - fly11$pos_y)/(sqrt(((fly11$pos_x - fly2$pos_x)^2) + ((fly11$pos_y - fly2$pos_y)^2))))
  anglearcsin$arcsin11.3 <- asin((fly3$pos_y - fly11$pos_y)/(sqrt(((fly11$pos_x - fly3$pos_x)^2) + ((fly11$pos_y - fly3$pos_y)^2))))
  anglearcsin$arcsin11.4 <- asin((fly4$pos_y - fly11$pos_y)/(sqrt(((fly11$pos_x - fly4$pos_x)^2) + ((fly11$pos_y - fly4$pos_y)^2))))
  anglearcsin$arcsin11.5 <- asin((fly5$pos_y - fly11$pos_y)/(sqrt(((fly11$pos_x - fly5$pos_x)^2) + ((fly11$pos_y - fly5$pos_y)^2))))
  anglearcsin$arcsin11.6 <- asin((fly6$pos_y - fly11$pos_y)/(sqrt(((fly11$pos_x - fly6$pos_x)^2) + ((fly11$pos_y - fly6$pos_y)^2))))
  anglearcsin$arcsin11.7 <- asin((fly7$pos_y - fly11$pos_y)/(sqrt(((fly11$pos_x - fly7$pos_x)^2) + ((fly11$pos_y - fly7$pos_y)^2))))
  anglearcsin$arcsin11.8 <- asin((fly8$pos_y - fly11$pos_y)/(sqrt(((fly11$pos_x - fly8$pos_x)^2) + ((fly11$pos_y - fly8$pos_y)^2))))
  anglearcsin$arcsin11.9 <- asin((fly9$pos_y - fly11$pos_y)/(sqrt(((fly11$pos_x - fly9$pos_x)^2) + ((fly11$pos_y - fly9$pos_y)^2))))
  anglearcsin$arcsin11.10 <- asin((fly10$pos_y - fly11$pos_y)/(sqrt(((fly11$pos_x - fly10$pos_x)^2) + ((fly11$pos_y - fly10$pos_y)^2))))
  anglearcsin$arcsin11.12 <- asin((fly12$pos_y - fly11$pos_y)/(sqrt(((fly11$pos_x - fly12$pos_x)^2) + ((fly11$pos_y - fly12$pos_y)^2))))
  anglearcsin$arcsin11.13 <- asin((fly13$pos_y - fly11$pos_y)/(sqrt(((fly11$pos_x - fly13$pos_x)^2) + ((fly11$pos_y - fly13$pos_y)^2))))
  anglearcsin$arcsin11.14 <- asin((fly14$pos_y - fly11$pos_y)/(sqrt(((fly11$pos_x - fly14$pos_x)^2) + ((fly11$pos_y - fly14$pos_y)^2))))
  anglearcsin$arcsin11.15 <- asin((fly15$pos_y - fly11$pos_y)/(sqrt(((fly11$pos_x - fly15$pos_x)^2) + ((fly11$pos_y - fly15$pos_y)^2))))
  anglearcsin$arcsin11.16 <- asin((fly16$pos_y - fly11$pos_y)/(sqrt(((fly11$pos_x - fly16$pos_x)^2) + ((fly11$pos_y - fly16$pos_y)^2))))
  anglearcsin$arcsin11.17 <- asin((fly17$pos_y - fly11$pos_y)/(sqrt(((fly11$pos_x - fly17$pos_x)^2) + ((fly11$pos_y - fly17$pos_y)^2))))
  anglearcsin$arcsin11.18 <- asin((fly18$pos_y - fly11$pos_y)/(sqrt(((fly11$pos_x - fly18$pos_x)^2) + ((fly11$pos_y - fly18$pos_y)^2))))
  anglearcsin$arcsin11.19 <- asin((fly19$pos_y - fly11$pos_y)/(sqrt(((fly11$pos_x - fly19$pos_x)^2) + ((fly11$pos_y - fly19$pos_y)^2))))
  anglearcsin$arcsin11.20 <- asin((fly20$pos_y - fly11$pos_y)/(sqrt(((fly11$pos_x - fly20$pos_x)^2) + ((fly11$pos_y - fly20$pos_y)^2))))
  anglearcsin$arcsin12.1 <- asin((fly1$pos_y - fly12$pos_y)/(sqrt(((fly12$pos_x - fly1$pos_x)^2) + ((fly12$pos_y - fly1$pos_y)^2))))
  anglearcsin$arcsin12.2 <- asin((fly2$pos_y - fly12$pos_y)/(sqrt(((fly12$pos_x - fly2$pos_x)^2) + ((fly12$pos_y - fly2$pos_y)^2))))
  anglearcsin$arcsin12.3 <- asin((fly3$pos_y - fly12$pos_y)/(sqrt(((fly12$pos_x - fly3$pos_x)^2) + ((fly12$pos_y - fly3$pos_y)^2))))
  anglearcsin$arcsin12.4 <- asin((fly4$pos_y - fly12$pos_y)/(sqrt(((fly12$pos_x - fly4$pos_x)^2) + ((fly12$pos_y - fly4$pos_y)^2))))
  anglearcsin$arcsin12.5 <- asin((fly5$pos_y - fly12$pos_y)/(sqrt(((fly12$pos_x - fly5$pos_x)^2) + ((fly12$pos_y - fly5$pos_y)^2))))
  anglearcsin$arcsin12.6 <- asin((fly6$pos_y - fly12$pos_y)/(sqrt(((fly12$pos_x - fly6$pos_x)^2) + ((fly12$pos_y - fly6$pos_y)^2))))
  anglearcsin$arcsin12.7 <- asin((fly7$pos_y - fly12$pos_y)/(sqrt(((fly12$pos_x - fly7$pos_x)^2) + ((fly12$pos_y - fly7$pos_y)^2))))
  anglearcsin$arcsin12.8 <- asin((fly8$pos_y - fly12$pos_y)/(sqrt(((fly12$pos_x - fly8$pos_x)^2) + ((fly12$pos_y - fly8$pos_y)^2))))
  anglearcsin$arcsin12.9 <- asin((fly9$pos_y - fly12$pos_y)/(sqrt(((fly12$pos_x - fly9$pos_x)^2) + ((fly12$pos_y - fly9$pos_y)^2))))
  anglearcsin$arcsin12.10 <- asin((fly10$pos_y - fly12$pos_y)/(sqrt(((fly12$pos_x - fly10$pos_x)^2) + ((fly12$pos_y - fly10$pos_y)^2))))
  anglearcsin$arcsin12.11 <- asin((fly11$pos_y - fly12$pos_y)/(sqrt(((fly12$pos_x - fly11$pos_x)^2) + ((fly12$pos_y - fly11$pos_y)^2))))
  anglearcsin$arcsin12.13 <- asin((fly13$pos_y - fly12$pos_y)/(sqrt(((fly12$pos_x - fly13$pos_x)^2) + ((fly12$pos_y - fly13$pos_y)^2))))
  anglearcsin$arcsin12.14 <- asin((fly14$pos_y - fly12$pos_y)/(sqrt(((fly12$pos_x - fly14$pos_x)^2) + ((fly12$pos_y - fly14$pos_y)^2))))
  anglearcsin$arcsin12.15 <- asin((fly15$pos_y - fly12$pos_y)/(sqrt(((fly12$pos_x - fly15$pos_x)^2) + ((fly12$pos_y - fly15$pos_y)^2))))
  anglearcsin$arcsin12.16 <- asin((fly16$pos_y - fly12$pos_y)/(sqrt(((fly12$pos_x - fly16$pos_x)^2) + ((fly12$pos_y - fly16$pos_y)^2))))
  anglearcsin$arcsin12.17 <- asin((fly17$pos_y - fly12$pos_y)/(sqrt(((fly12$pos_x - fly17$pos_x)^2) + ((fly12$pos_y - fly17$pos_y)^2))))
  anglearcsin$arcsin12.18 <- asin((fly18$pos_y - fly12$pos_y)/(sqrt(((fly12$pos_x - fly18$pos_x)^2) + ((fly12$pos_y - fly18$pos_y)^2))))
  anglearcsin$arcsin12.19 <- asin((fly19$pos_y - fly12$pos_y)/(sqrt(((fly12$pos_x - fly19$pos_x)^2) + ((fly12$pos_y - fly19$pos_y)^2))))
  anglearcsin$arcsin12.20 <- asin((fly20$pos_y - fly12$pos_y)/(sqrt(((fly12$pos_x - fly20$pos_x)^2) + ((fly12$pos_y - fly20$pos_y)^2))))
  anglearcsin$arcsin13.1 <- asin((fly1$pos_y - fly13$pos_y)/(sqrt(((fly13$pos_x - fly1$pos_x)^2) + ((fly13$pos_y - fly1$pos_y)^2))))
  anglearcsin$arcsin13.2 <- asin((fly2$pos_y - fly13$pos_y)/(sqrt(((fly13$pos_x - fly2$pos_x)^2) + ((fly13$pos_y - fly2$pos_y)^2))))
  anglearcsin$arcsin13.3 <- asin((fly3$pos_y - fly13$pos_y)/(sqrt(((fly13$pos_x - fly3$pos_x)^2) + ((fly13$pos_y - fly3$pos_y)^2))))
  anglearcsin$arcsin13.4 <- asin((fly4$pos_y - fly13$pos_y)/(sqrt(((fly13$pos_x - fly4$pos_x)^2) + ((fly13$pos_y - fly4$pos_y)^2))))
  anglearcsin$arcsin13.5 <- asin((fly5$pos_y - fly13$pos_y)/(sqrt(((fly13$pos_x - fly5$pos_x)^2) + ((fly13$pos_y - fly5$pos_y)^2))))
  anglearcsin$arcsin13.6 <- asin((fly6$pos_y - fly13$pos_y)/(sqrt(((fly13$pos_x - fly6$pos_x)^2) + ((fly13$pos_y - fly6$pos_y)^2))))
  anglearcsin$arcsin13.7 <- asin((fly7$pos_y - fly13$pos_y)/(sqrt(((fly13$pos_x - fly7$pos_x)^2) + ((fly13$pos_y - fly7$pos_y)^2))))
  anglearcsin$arcsin13.8 <- asin((fly8$pos_y - fly13$pos_y)/(sqrt(((fly13$pos_x - fly8$pos_x)^2) + ((fly13$pos_y - fly8$pos_y)^2))))
  anglearcsin$arcsin13.9 <- asin((fly9$pos_y - fly13$pos_y)/(sqrt(((fly13$pos_x - fly9$pos_x)^2) + ((fly13$pos_y - fly9$pos_y)^2))))
  anglearcsin$arcsin13.10 <- asin((fly10$pos_y - fly13$pos_y)/(sqrt(((fly13$pos_x - fly10$pos_x)^2) + ((fly13$pos_y - fly10$pos_y)^2))))
  anglearcsin$arcsin13.11 <- asin((fly11$pos_y - fly13$pos_y)/(sqrt(((fly13$pos_x - fly11$pos_x)^2) + ((fly13$pos_y - fly11$pos_y)^2))))
  anglearcsin$arcsin13.12 <- asin((fly12$pos_y - fly13$pos_y)/(sqrt(((fly13$pos_x - fly12$pos_x)^2) + ((fly13$pos_y - fly12$pos_y)^2))))
  anglearcsin$arcsin13.14 <- asin((fly14$pos_y - fly13$pos_y)/(sqrt(((fly13$pos_x - fly14$pos_x)^2) + ((fly13$pos_y - fly14$pos_y)^2))))
  anglearcsin$arcsin13.15 <- asin((fly15$pos_y - fly13$pos_y)/(sqrt(((fly13$pos_x - fly15$pos_x)^2) + ((fly13$pos_y - fly15$pos_y)^2))))
  anglearcsin$arcsin13.16 <- asin((fly16$pos_y - fly13$pos_y)/(sqrt(((fly13$pos_x - fly16$pos_x)^2) + ((fly13$pos_y - fly16$pos_y)^2))))
  anglearcsin$arcsin13.17 <- asin((fly17$pos_y - fly13$pos_y)/(sqrt(((fly13$pos_x - fly17$pos_x)^2) + ((fly13$pos_y - fly17$pos_y)^2))))
  anglearcsin$arcsin13.18 <- asin((fly18$pos_y - fly13$pos_y)/(sqrt(((fly13$pos_x - fly18$pos_x)^2) + ((fly13$pos_y - fly18$pos_y)^2))))
  anglearcsin$arcsin13.19 <- asin((fly19$pos_y - fly13$pos_y)/(sqrt(((fly13$pos_x - fly19$pos_x)^2) + ((fly13$pos_y - fly19$pos_y)^2))))
  anglearcsin$arcsin13.20 <- asin((fly20$pos_y - fly13$pos_y)/(sqrt(((fly13$pos_x - fly20$pos_x)^2) + ((fly13$pos_y - fly20$pos_y)^2))))
  anglearcsin$arcsin14.1 <- asin((fly1$pos_y - fly14$pos_y)/(sqrt(((fly14$pos_x - fly1$pos_x)^2) + ((fly14$pos_y - fly1$pos_y)^2))))
  anglearcsin$arcsin14.2 <- asin((fly2$pos_y - fly14$pos_y)/(sqrt(((fly14$pos_x - fly2$pos_x)^2) + ((fly14$pos_y - fly2$pos_y)^2))))
  anglearcsin$arcsin14.3 <- asin((fly3$pos_y - fly14$pos_y)/(sqrt(((fly14$pos_x - fly3$pos_x)^2) + ((fly14$pos_y - fly3$pos_y)^2))))
  anglearcsin$arcsin14.4 <- asin((fly4$pos_y - fly14$pos_y)/(sqrt(((fly14$pos_x - fly4$pos_x)^2) + ((fly14$pos_y - fly4$pos_y)^2))))
  anglearcsin$arcsin14.5 <- asin((fly5$pos_y - fly14$pos_y)/(sqrt(((fly14$pos_x - fly5$pos_x)^2) + ((fly14$pos_y - fly5$pos_y)^2))))
  anglearcsin$arcsin14.6 <- asin((fly6$pos_y - fly14$pos_y)/(sqrt(((fly14$pos_x - fly6$pos_x)^2) + ((fly14$pos_y - fly6$pos_y)^2))))
  anglearcsin$arcsin14.7 <- asin((fly7$pos_y - fly14$pos_y)/(sqrt(((fly14$pos_x - fly7$pos_x)^2) + ((fly14$pos_y - fly7$pos_y)^2))))
  anglearcsin$arcsin14.8 <- asin((fly8$pos_y - fly14$pos_y)/(sqrt(((fly14$pos_x - fly8$pos_x)^2) + ((fly14$pos_y - fly8$pos_y)^2))))
  anglearcsin$arcsin14.9 <- asin((fly9$pos_y - fly14$pos_y)/(sqrt(((fly14$pos_x - fly9$pos_x)^2) + ((fly14$pos_y - fly9$pos_y)^2))))
  anglearcsin$arcsin14.10 <- asin((fly10$pos_y - fly14$pos_y)/(sqrt(((fly14$pos_x - fly10$pos_x)^2) + ((fly14$pos_y - fly10$pos_y)^2))))
  anglearcsin$arcsin14.11 <- asin((fly11$pos_y - fly14$pos_y)/(sqrt(((fly14$pos_x - fly11$pos_x)^2) + ((fly14$pos_y - fly11$pos_y)^2))))
  anglearcsin$arcsin14.12 <- asin((fly12$pos_y - fly14$pos_y)/(sqrt(((fly14$pos_x - fly12$pos_x)^2) + ((fly14$pos_y - fly12$pos_y)^2))))
  anglearcsin$arcsin14.13 <- asin((fly13$pos_y - fly14$pos_y)/(sqrt(((fly14$pos_x - fly13$pos_x)^2) + ((fly14$pos_y - fly13$pos_y)^2))))
  anglearcsin$arcsin14.15 <- asin((fly15$pos_y - fly14$pos_y)/(sqrt(((fly14$pos_x - fly15$pos_x)^2) + ((fly14$pos_y - fly15$pos_y)^2))))
  anglearcsin$arcsin14.16 <- asin((fly16$pos_y - fly14$pos_y)/(sqrt(((fly14$pos_x - fly16$pos_x)^2) + ((fly14$pos_y - fly16$pos_y)^2))))
  anglearcsin$arcsin14.17 <- asin((fly17$pos_y - fly14$pos_y)/(sqrt(((fly14$pos_x - fly17$pos_x)^2) + ((fly14$pos_y - fly17$pos_y)^2))))
  anglearcsin$arcsin14.18 <- asin((fly18$pos_y - fly14$pos_y)/(sqrt(((fly14$pos_x - fly18$pos_x)^2) + ((fly14$pos_y - fly18$pos_y)^2))))
  anglearcsin$arcsin14.19 <- asin((fly19$pos_y - fly14$pos_y)/(sqrt(((fly14$pos_x - fly19$pos_x)^2) + ((fly14$pos_y - fly19$pos_y)^2))))
  anglearcsin$arcsin14.20 <- asin((fly20$pos_y - fly14$pos_y)/(sqrt(((fly14$pos_x - fly20$pos_x)^2) + ((fly14$pos_y - fly20$pos_y)^2))))
  anglearcsin$arcsin15.1 <- asin((fly1$pos_y - fly15$pos_y)/(sqrt(((fly15$pos_x - fly1$pos_x)^2) + ((fly15$pos_y - fly1$pos_y)^2))))
  anglearcsin$arcsin15.2 <- asin((fly2$pos_y - fly15$pos_y)/(sqrt(((fly15$pos_x - fly2$pos_x)^2) + ((fly15$pos_y - fly2$pos_y)^2))))
  anglearcsin$arcsin15.3 <- asin((fly3$pos_y - fly15$pos_y)/(sqrt(((fly15$pos_x - fly3$pos_x)^2) + ((fly15$pos_y - fly3$pos_y)^2))))
  anglearcsin$arcsin15.4 <- asin((fly4$pos_y - fly15$pos_y)/(sqrt(((fly15$pos_x - fly4$pos_x)^2) + ((fly15$pos_y - fly4$pos_y)^2))))
  anglearcsin$arcsin15.5 <- asin((fly5$pos_y - fly15$pos_y)/(sqrt(((fly15$pos_x - fly5$pos_x)^2) + ((fly15$pos_y - fly5$pos_y)^2))))
  anglearcsin$arcsin15.6 <- asin((fly6$pos_y - fly15$pos_y)/(sqrt(((fly15$pos_x - fly6$pos_x)^2) + ((fly15$pos_y - fly6$pos_y)^2))))
  anglearcsin$arcsin15.7 <- asin((fly7$pos_y - fly15$pos_y)/(sqrt(((fly15$pos_x - fly7$pos_x)^2) + ((fly15$pos_y - fly7$pos_y)^2))))
  anglearcsin$arcsin15.8 <- asin((fly8$pos_y - fly15$pos_y)/(sqrt(((fly15$pos_x - fly8$pos_x)^2) + ((fly15$pos_y - fly8$pos_y)^2))))
  anglearcsin$arcsin15.9 <- asin((fly9$pos_y - fly15$pos_y)/(sqrt(((fly15$pos_x - fly9$pos_x)^2) + ((fly15$pos_y - fly9$pos_y)^2))))
  anglearcsin$arcsin15.10 <- asin((fly10$pos_y - fly15$pos_y)/(sqrt(((fly15$pos_x - fly10$pos_x)^2) + ((fly15$pos_y - fly10$pos_y)^2))))
  anglearcsin$arcsin15.11 <- asin((fly11$pos_y - fly15$pos_y)/(sqrt(((fly15$pos_x - fly11$pos_x)^2) + ((fly15$pos_y - fly11$pos_y)^2))))
  anglearcsin$arcsin15.12 <- asin((fly12$pos_y - fly15$pos_y)/(sqrt(((fly15$pos_x - fly12$pos_x)^2) + ((fly15$pos_y - fly12$pos_y)^2))))
  anglearcsin$arcsin15.13 <- asin((fly13$pos_y - fly15$pos_y)/(sqrt(((fly15$pos_x - fly13$pos_x)^2) + ((fly15$pos_y - fly13$pos_y)^2))))
  anglearcsin$arcsin15.14 <- asin((fly14$pos_y - fly15$pos_y)/(sqrt(((fly15$pos_x - fly14$pos_x)^2) + ((fly15$pos_y - fly14$pos_y)^2))))
  anglearcsin$arcsin15.16 <- asin((fly16$pos_y - fly15$pos_y)/(sqrt(((fly15$pos_x - fly16$pos_x)^2) + ((fly15$pos_y - fly16$pos_y)^2))))
  anglearcsin$arcsin15.17 <- asin((fly17$pos_y - fly15$pos_y)/(sqrt(((fly15$pos_x - fly17$pos_x)^2) + ((fly15$pos_y - fly17$pos_y)^2))))
  anglearcsin$arcsin15.18 <- asin((fly18$pos_y - fly15$pos_y)/(sqrt(((fly15$pos_x - fly18$pos_x)^2) + ((fly15$pos_y - fly18$pos_y)^2))))
  anglearcsin$arcsin15.19 <- asin((fly19$pos_y - fly15$pos_y)/(sqrt(((fly15$pos_x - fly19$pos_x)^2) + ((fly15$pos_y - fly19$pos_y)^2))))
  anglearcsin$arcsin15.20 <- asin((fly20$pos_y - fly15$pos_y)/(sqrt(((fly15$pos_x - fly20$pos_x)^2) + ((fly15$pos_y - fly20$pos_y)^2))))
  anglearcsin$arcsin16.1 <- asin((fly1$pos_y - fly16$pos_y)/(sqrt(((fly16$pos_x - fly1$pos_x)^2) + ((fly16$pos_y - fly1$pos_y)^2))))
  anglearcsin$arcsin16.2 <- asin((fly2$pos_y - fly16$pos_y)/(sqrt(((fly16$pos_x - fly2$pos_x)^2) + ((fly16$pos_y - fly2$pos_y)^2))))
  anglearcsin$arcsin16.3 <- asin((fly3$pos_y - fly16$pos_y)/(sqrt(((fly16$pos_x - fly3$pos_x)^2) + ((fly16$pos_y - fly3$pos_y)^2))))
  anglearcsin$arcsin16.4 <- asin((fly4$pos_y - fly16$pos_y)/(sqrt(((fly16$pos_x - fly4$pos_x)^2) + ((fly16$pos_y - fly4$pos_y)^2))))
  anglearcsin$arcsin16.5 <- asin((fly5$pos_y - fly16$pos_y)/(sqrt(((fly16$pos_x - fly5$pos_x)^2) + ((fly16$pos_y - fly5$pos_y)^2))))
  anglearcsin$arcsin16.6 <- asin((fly6$pos_y - fly16$pos_y)/(sqrt(((fly16$pos_x - fly6$pos_x)^2) + ((fly16$pos_y - fly6$pos_y)^2))))
  anglearcsin$arcsin16.7 <- asin((fly7$pos_y - fly16$pos_y)/(sqrt(((fly16$pos_x - fly7$pos_x)^2) + ((fly16$pos_y - fly7$pos_y)^2))))
  anglearcsin$arcsin16.8 <- asin((fly8$pos_y - fly16$pos_y)/(sqrt(((fly16$pos_x - fly8$pos_x)^2) + ((fly16$pos_y - fly8$pos_y)^2))))
  anglearcsin$arcsin16.9 <- asin((fly9$pos_y - fly16$pos_y)/(sqrt(((fly16$pos_x - fly9$pos_x)^2) + ((fly16$pos_y - fly9$pos_y)^2))))
  anglearcsin$arcsin16.10 <- asin((fly10$pos_y - fly16$pos_y)/(sqrt(((fly16$pos_x - fly10$pos_x)^2) + ((fly16$pos_y - fly10$pos_y)^2))))
  anglearcsin$arcsin16.11 <- asin((fly11$pos_y - fly16$pos_y)/(sqrt(((fly16$pos_x - fly11$pos_x)^2) + ((fly16$pos_y - fly11$pos_y)^2))))
  anglearcsin$arcsin16.12 <- asin((fly12$pos_y - fly16$pos_y)/(sqrt(((fly16$pos_x - fly12$pos_x)^2) + ((fly16$pos_y - fly12$pos_y)^2))))
  anglearcsin$arcsin16.13 <- asin((fly13$pos_y - fly16$pos_y)/(sqrt(((fly16$pos_x - fly13$pos_x)^2) + ((fly16$pos_y - fly13$pos_y)^2))))
  anglearcsin$arcsin16.14 <- asin((fly14$pos_y - fly16$pos_y)/(sqrt(((fly16$pos_x - fly14$pos_x)^2) + ((fly16$pos_y - fly14$pos_y)^2))))
  anglearcsin$arcsin16.15 <- asin((fly15$pos_y - fly16$pos_y)/(sqrt(((fly16$pos_x - fly15$pos_x)^2) + ((fly16$pos_y - fly15$pos_y)^2))))
  anglearcsin$arcsin16.17 <- asin((fly17$pos_y - fly16$pos_y)/(sqrt(((fly16$pos_x - fly17$pos_x)^2) + ((fly16$pos_y - fly17$pos_y)^2))))
  anglearcsin$arcsin16.18 <- asin((fly18$pos_y - fly16$pos_y)/(sqrt(((fly16$pos_x - fly18$pos_x)^2) + ((fly16$pos_y - fly18$pos_y)^2))))
  anglearcsin$arcsin16.19 <- asin((fly19$pos_y - fly16$pos_y)/(sqrt(((fly16$pos_x - fly19$pos_x)^2) + ((fly16$pos_y - fly19$pos_y)^2))))
  anglearcsin$arcsin16.20 <- asin((fly20$pos_y - fly16$pos_y)/(sqrt(((fly16$pos_x - fly20$pos_x)^2) + ((fly16$pos_y - fly20$pos_y)^2))))
  anglearcsin$arcsin17.1 <- asin((fly1$pos_y - fly17$pos_y)/(sqrt(((fly17$pos_x - fly1$pos_x)^2) + ((fly17$pos_y - fly1$pos_y)^2))))
  anglearcsin$arcsin17.2 <- asin((fly2$pos_y - fly17$pos_y)/(sqrt(((fly17$pos_x - fly2$pos_x)^2) + ((fly17$pos_y - fly2$pos_y)^2))))
  anglearcsin$arcsin17.3 <- asin((fly3$pos_y - fly17$pos_y)/(sqrt(((fly17$pos_x - fly3$pos_x)^2) + ((fly17$pos_y - fly3$pos_y)^2))))
  anglearcsin$arcsin17.4 <- asin((fly4$pos_y - fly17$pos_y)/(sqrt(((fly17$pos_x - fly4$pos_x)^2) + ((fly17$pos_y - fly4$pos_y)^2))))
  anglearcsin$arcsin17.5 <- asin((fly5$pos_y - fly17$pos_y)/(sqrt(((fly17$pos_x - fly5$pos_x)^2) + ((fly17$pos_y - fly5$pos_y)^2))))
  anglearcsin$arcsin17.6 <- asin((fly6$pos_y - fly17$pos_y)/(sqrt(((fly17$pos_x - fly6$pos_x)^2) + ((fly17$pos_y - fly6$pos_y)^2))))
  anglearcsin$arcsin17.7 <- asin((fly7$pos_y - fly17$pos_y)/(sqrt(((fly17$pos_x - fly7$pos_x)^2) + ((fly17$pos_y - fly7$pos_y)^2))))
  anglearcsin$arcsin17.8 <- asin((fly8$pos_y - fly17$pos_y)/(sqrt(((fly17$pos_x - fly8$pos_x)^2) + ((fly17$pos_y - fly8$pos_y)^2))))
  anglearcsin$arcsin17.9 <- asin((fly9$pos_y - fly17$pos_y)/(sqrt(((fly17$pos_x - fly9$pos_x)^2) + ((fly17$pos_y - fly9$pos_y)^2))))
  anglearcsin$arcsin17.10 <- asin((fly10$pos_y - fly17$pos_y)/(sqrt(((fly17$pos_x - fly10$pos_x)^2) + ((fly17$pos_y - fly10$pos_y)^2))))
  anglearcsin$arcsin17.11 <- asin((fly11$pos_y - fly17$pos_y)/(sqrt(((fly17$pos_x - fly11$pos_x)^2) + ((fly17$pos_y - fly11$pos_y)^2))))
  anglearcsin$arcsin17.12 <- asin((fly12$pos_y - fly17$pos_y)/(sqrt(((fly17$pos_x - fly12$pos_x)^2) + ((fly17$pos_y - fly12$pos_y)^2))))
  anglearcsin$arcsin17.13 <- asin((fly13$pos_y - fly17$pos_y)/(sqrt(((fly17$pos_x - fly13$pos_x)^2) + ((fly17$pos_y - fly13$pos_y)^2))))
  anglearcsin$arcsin17.14 <- asin((fly14$pos_y - fly17$pos_y)/(sqrt(((fly17$pos_x - fly14$pos_x)^2) + ((fly17$pos_y - fly14$pos_y)^2))))
  anglearcsin$arcsin17.15 <- asin((fly15$pos_y - fly17$pos_y)/(sqrt(((fly17$pos_x - fly15$pos_x)^2) + ((fly17$pos_y - fly15$pos_y)^2))))
  anglearcsin$arcsin17.16 <- asin((fly16$pos_y - fly17$pos_y)/(sqrt(((fly17$pos_x - fly16$pos_x)^2) + ((fly17$pos_y - fly16$pos_y)^2))))
  anglearcsin$arcsin17.18 <- asin((fly18$pos_y - fly17$pos_y)/(sqrt(((fly17$pos_x - fly18$pos_x)^2) + ((fly17$pos_y - fly18$pos_y)^2))))
  anglearcsin$arcsin17.19 <- asin((fly19$pos_y - fly17$pos_y)/(sqrt(((fly17$pos_x - fly19$pos_x)^2) + ((fly17$pos_y - fly19$pos_y)^2))))
  anglearcsin$arcsin17.20 <- asin((fly20$pos_y - fly17$pos_y)/(sqrt(((fly17$pos_x - fly20$pos_x)^2) + ((fly17$pos_y - fly20$pos_y)^2))))
  anglearcsin$arcsin18.1 <- asin((fly1$pos_y - fly18$pos_y)/(sqrt(((fly18$pos_x - fly1$pos_x)^2) + ((fly18$pos_y - fly1$pos_y)^2))))
  anglearcsin$arcsin18.2 <- asin((fly2$pos_y - fly18$pos_y)/(sqrt(((fly18$pos_x - fly2$pos_x)^2) + ((fly18$pos_y - fly2$pos_y)^2))))
  anglearcsin$arcsin18.3 <- asin((fly3$pos_y - fly18$pos_y)/(sqrt(((fly18$pos_x - fly3$pos_x)^2) + ((fly18$pos_y - fly3$pos_y)^2))))
  anglearcsin$arcsin18.4 <- asin((fly4$pos_y - fly18$pos_y)/(sqrt(((fly18$pos_x - fly4$pos_x)^2) + ((fly18$pos_y - fly4$pos_y)^2))))
  anglearcsin$arcsin18.5 <- asin((fly5$pos_y - fly18$pos_y)/(sqrt(((fly18$pos_x - fly5$pos_x)^2) + ((fly18$pos_y - fly5$pos_y)^2))))
  anglearcsin$arcsin18.6 <- asin((fly6$pos_y - fly18$pos_y)/(sqrt(((fly18$pos_x - fly6$pos_x)^2) + ((fly18$pos_y - fly6$pos_y)^2))))
  anglearcsin$arcsin18.7 <- asin((fly7$pos_y - fly18$pos_y)/(sqrt(((fly18$pos_x - fly7$pos_x)^2) + ((fly18$pos_y - fly7$pos_y)^2))))
  anglearcsin$arcsin18.8 <- asin((fly8$pos_y - fly18$pos_y)/(sqrt(((fly18$pos_x - fly8$pos_x)^2) + ((fly18$pos_y - fly8$pos_y)^2))))
  anglearcsin$arcsin18.9 <- asin((fly9$pos_y - fly18$pos_y)/(sqrt(((fly18$pos_x - fly9$pos_x)^2) + ((fly18$pos_y - fly9$pos_y)^2))))
  anglearcsin$arcsin18.10 <- asin((fly10$pos_y - fly18$pos_y)/(sqrt(((fly18$pos_x - fly10$pos_x)^2) + ((fly18$pos_y - fly10$pos_y)^2))))
  anglearcsin$arcsin18.11 <- asin((fly11$pos_y - fly18$pos_y)/(sqrt(((fly18$pos_x - fly11$pos_x)^2) + ((fly18$pos_y - fly11$pos_y)^2))))
  anglearcsin$arcsin18.12 <- asin((fly12$pos_y - fly18$pos_y)/(sqrt(((fly18$pos_x - fly12$pos_x)^2) + ((fly18$pos_y - fly12$pos_y)^2))))
  anglearcsin$arcsin18.13 <- asin((fly13$pos_y - fly18$pos_y)/(sqrt(((fly18$pos_x - fly13$pos_x)^2) + ((fly18$pos_y - fly13$pos_y)^2))))
  anglearcsin$arcsin18.14 <- asin((fly14$pos_y - fly18$pos_y)/(sqrt(((fly18$pos_x - fly14$pos_x)^2) + ((fly18$pos_y - fly14$pos_y)^2))))
  anglearcsin$arcsin18.15 <- asin((fly15$pos_y - fly18$pos_y)/(sqrt(((fly18$pos_x - fly15$pos_x)^2) + ((fly18$pos_y - fly15$pos_y)^2))))
  anglearcsin$arcsin18.16 <- asin((fly16$pos_y - fly18$pos_y)/(sqrt(((fly18$pos_x - fly16$pos_x)^2) + ((fly18$pos_y - fly16$pos_y)^2))))
  anglearcsin$arcsin18.17 <- asin((fly17$pos_y - fly18$pos_y)/(sqrt(((fly18$pos_x - fly17$pos_x)^2) + ((fly18$pos_y - fly17$pos_y)^2))))
  anglearcsin$arcsin18.19 <- asin((fly19$pos_y - fly18$pos_y)/(sqrt(((fly18$pos_x - fly19$pos_x)^2) + ((fly18$pos_y - fly19$pos_y)^2))))
  anglearcsin$arcsin18.20 <- asin((fly20$pos_y - fly18$pos_y)/(sqrt(((fly18$pos_x - fly20$pos_x)^2) + ((fly18$pos_y - fly20$pos_y)^2))))
  anglearcsin$arcsin19.1 <- asin((fly1$pos_y - fly19$pos_y)/(sqrt(((fly19$pos_x - fly1$pos_x)^2) + ((fly19$pos_y - fly1$pos_y)^2))))
  anglearcsin$arcsin19.2 <- asin((fly2$pos_y - fly19$pos_y)/(sqrt(((fly19$pos_x - fly2$pos_x)^2) + ((fly19$pos_y - fly2$pos_y)^2))))
  anglearcsin$arcsin19.3 <- asin((fly3$pos_y - fly19$pos_y)/(sqrt(((fly19$pos_x - fly3$pos_x)^2) + ((fly19$pos_y - fly3$pos_y)^2))))
  anglearcsin$arcsin19.4 <- asin((fly4$pos_y - fly19$pos_y)/(sqrt(((fly19$pos_x - fly4$pos_x)^2) + ((fly19$pos_y - fly4$pos_y)^2))))
  anglearcsin$arcsin19.5 <- asin((fly5$pos_y - fly19$pos_y)/(sqrt(((fly19$pos_x - fly5$pos_x)^2) + ((fly19$pos_y - fly5$pos_y)^2))))
  anglearcsin$arcsin19.6 <- asin((fly6$pos_y - fly19$pos_y)/(sqrt(((fly19$pos_x - fly6$pos_x)^2) + ((fly19$pos_y - fly6$pos_y)^2))))
  anglearcsin$arcsin19.7 <- asin((fly7$pos_y - fly19$pos_y)/(sqrt(((fly19$pos_x - fly7$pos_x)^2) + ((fly19$pos_y - fly7$pos_y)^2))))
  anglearcsin$arcsin19.8 <- asin((fly8$pos_y - fly19$pos_y)/(sqrt(((fly19$pos_x - fly8$pos_x)^2) + ((fly19$pos_y - fly8$pos_y)^2))))
  anglearcsin$arcsin19.9 <- asin((fly9$pos_y - fly19$pos_y)/(sqrt(((fly19$pos_x - fly9$pos_x)^2) + ((fly19$pos_y - fly9$pos_y)^2))))
  anglearcsin$arcsin19.10 <- asin((fly10$pos_y - fly19$pos_y)/(sqrt(((fly19$pos_x - fly10$pos_x)^2) + ((fly19$pos_y - fly10$pos_y)^2))))
  anglearcsin$arcsin19.11 <- asin((fly11$pos_y - fly19$pos_y)/(sqrt(((fly19$pos_x - fly11$pos_x)^2) + ((fly19$pos_y - fly11$pos_y)^2))))
  anglearcsin$arcsin19.12 <- asin((fly12$pos_y - fly19$pos_y)/(sqrt(((fly19$pos_x - fly12$pos_x)^2) + ((fly19$pos_y - fly12$pos_y)^2))))
  anglearcsin$arcsin19.13 <- asin((fly13$pos_y - fly19$pos_y)/(sqrt(((fly19$pos_x - fly13$pos_x)^2) + ((fly19$pos_y - fly13$pos_y)^2))))
  anglearcsin$arcsin19.14 <- asin((fly14$pos_y - fly19$pos_y)/(sqrt(((fly19$pos_x - fly14$pos_x)^2) + ((fly19$pos_y - fly14$pos_y)^2))))
  anglearcsin$arcsin19.15 <- asin((fly15$pos_y - fly19$pos_y)/(sqrt(((fly19$pos_x - fly15$pos_x)^2) + ((fly19$pos_y - fly15$pos_y)^2))))
  anglearcsin$arcsin19.16 <- asin((fly16$pos_y - fly19$pos_y)/(sqrt(((fly19$pos_x - fly16$pos_x)^2) + ((fly19$pos_y - fly16$pos_y)^2))))
  anglearcsin$arcsin19.17 <- asin((fly17$pos_y - fly19$pos_y)/(sqrt(((fly19$pos_x - fly17$pos_x)^2) + ((fly19$pos_y - fly17$pos_y)^2))))
  anglearcsin$arcsin19.18 <- asin((fly18$pos_y - fly19$pos_y)/(sqrt(((fly19$pos_x - fly18$pos_x)^2) + ((fly19$pos_y - fly18$pos_y)^2))))
  anglearcsin$arcsin19.20 <- asin((fly20$pos_y - fly19$pos_y)/(sqrt(((fly19$pos_x - fly20$pos_x)^2) + ((fly19$pos_y - fly20$pos_y)^2))))
  anglearcsin$arcsin20.1 <- asin((fly1$pos_y - fly20$pos_y)/(sqrt(((fly20$pos_x - fly1$pos_x)^2) + ((fly20$pos_y - fly1$pos_y)^2))))
  anglearcsin$arcsin20.2 <- asin((fly2$pos_y - fly20$pos_y)/(sqrt(((fly20$pos_x - fly2$pos_x)^2) + ((fly20$pos_y - fly2$pos_y)^2))))
  anglearcsin$arcsin20.3 <- asin((fly3$pos_y - fly20$pos_y)/(sqrt(((fly20$pos_x - fly3$pos_x)^2) + ((fly20$pos_y - fly3$pos_y)^2))))
  anglearcsin$arcsin20.4 <- asin((fly4$pos_y - fly20$pos_y)/(sqrt(((fly20$pos_x - fly4$pos_x)^2) + ((fly20$pos_y - fly4$pos_y)^2))))
  anglearcsin$arcsin20.5 <- asin((fly5$pos_y - fly20$pos_y)/(sqrt(((fly20$pos_x - fly5$pos_x)^2) + ((fly20$pos_y - fly5$pos_y)^2))))
  anglearcsin$arcsin20.6 <- asin((fly6$pos_y - fly20$pos_y)/(sqrt(((fly20$pos_x - fly6$pos_x)^2) + ((fly20$pos_y - fly6$pos_y)^2))))
  anglearcsin$arcsin20.7 <- asin((fly7$pos_y - fly20$pos_y)/(sqrt(((fly20$pos_x - fly7$pos_x)^2) + ((fly20$pos_y - fly7$pos_y)^2))))
  anglearcsin$arcsin20.8 <- asin((fly8$pos_y - fly20$pos_y)/(sqrt(((fly20$pos_x - fly8$pos_x)^2) + ((fly20$pos_y - fly8$pos_y)^2))))
  anglearcsin$arcsin20.9 <- asin((fly9$pos_y - fly20$pos_y)/(sqrt(((fly20$pos_x - fly9$pos_x)^2) + ((fly20$pos_y - fly9$pos_y)^2))))
  anglearcsin$arcsin20.10 <- asin((fly10$pos_y - fly20$pos_y)/(sqrt(((fly20$pos_x - fly10$pos_x)^2) + ((fly20$pos_y - fly10$pos_y)^2))))
  anglearcsin$arcsin20.11 <- asin((fly11$pos_y - fly20$pos_y)/(sqrt(((fly20$pos_x - fly11$pos_x)^2) + ((fly20$pos_y - fly11$pos_y)^2))))
  anglearcsin$arcsin20.12 <- asin((fly12$pos_y - fly20$pos_y)/(sqrt(((fly20$pos_x - fly12$pos_x)^2) + ((fly20$pos_y - fly12$pos_y)^2))))
  anglearcsin$arcsin20.13 <- asin((fly13$pos_y - fly20$pos_y)/(sqrt(((fly20$pos_x - fly13$pos_x)^2) + ((fly20$pos_y - fly13$pos_y)^2))))
  anglearcsin$arcsin20.14 <- asin((fly14$pos_y - fly20$pos_y)/(sqrt(((fly20$pos_x - fly14$pos_x)^2) + ((fly20$pos_y - fly14$pos_y)^2))))
  anglearcsin$arcsin20.15 <- asin((fly15$pos_y - fly20$pos_y)/(sqrt(((fly20$pos_x - fly15$pos_x)^2) + ((fly20$pos_y - fly15$pos_y)^2))))
  anglearcsin$arcsin20.16 <- asin((fly16$pos_y - fly20$pos_y)/(sqrt(((fly20$pos_x - fly16$pos_x)^2) + ((fly20$pos_y - fly16$pos_y)^2))))
  anglearcsin$arcsin20.17 <- asin((fly17$pos_y - fly20$pos_y)/(sqrt(((fly20$pos_x - fly17$pos_x)^2) + ((fly20$pos_y - fly17$pos_y)^2))))
  anglearcsin$arcsin20.18 <- asin((fly18$pos_y - fly20$pos_y)/(sqrt(((fly20$pos_x - fly18$pos_x)^2) + ((fly20$pos_y - fly18$pos_y)^2))))
  anglearcsin$arcsin20.19 <- asin((fly19$pos_y - fly20$pos_y)/(sqrt(((fly20$pos_x - fly19$pos_x)^2) + ((fly20$pos_y - fly19$pos_y)^2))))
  
  # CREATE A DATAFRAME THAT CALCULATES THE ANGLE BETWEEN EVERY PAIRWISE COMBINATION OF FLIES # ---------
  angles <- data.frame(ifelse(fly2$pos_x <= fly1$pos_x, c((anglearcsin$arcsin1.2*-1)+pi), c(ifelse((fly2$pos_x > fly1$pos_x & fly2$pos_y <= fly1$pos_y), c(anglearcsin$arcsin1.2+(2*pi)), c(anglearcsin$arcsin1.2)))))
  colnames(angles) <- c("ang1.2")
  angles$ang1.3 <- ifelse(fly3$pos_x <= fly1$pos_x, c((anglearcsin$arcsin1.3*-1)+pi), c(ifelse((fly3$pos_x > fly1$pos_x & fly3$pos_y <= fly1$pos_y), c(anglearcsin$arcsin1.3+(2*pi)), c(anglearcsin$arcsin1.3))))
  angles$ang1.4 <- ifelse(fly4$pos_x <= fly1$pos_x, c((anglearcsin$arcsin1.4*-1)+pi), c(ifelse((fly4$pos_x > fly1$pos_x & fly4$pos_y <= fly1$pos_y), c(anglearcsin$arcsin1.4+(2*pi)), c(anglearcsin$arcsin1.4))))
  angles$ang1.5 <- ifelse(fly5$pos_x <= fly1$pos_x, c((anglearcsin$arcsin1.5*-1)+pi), c(ifelse((fly5$pos_x > fly1$pos_x & fly5$pos_y <= fly1$pos_y), c(anglearcsin$arcsin1.5+(2*pi)), c(anglearcsin$arcsin1.5))))
  angles$ang1.6 <- ifelse(fly6$pos_x <= fly1$pos_x, c((anglearcsin$arcsin1.6*-1)+pi), c(ifelse((fly6$pos_x > fly1$pos_x & fly6$pos_y <= fly1$pos_y), c(anglearcsin$arcsin1.6+(2*pi)), c(anglearcsin$arcsin1.6))))
  angles$ang1.7 <- ifelse(fly7$pos_x <= fly1$pos_x, c((anglearcsin$arcsin1.7*-1)+pi), c(ifelse((fly7$pos_x > fly1$pos_x & fly7$pos_y <= fly1$pos_y), c(anglearcsin$arcsin1.7+(2*pi)), c(anglearcsin$arcsin1.7))))
  angles$ang1.8 <- ifelse(fly8$pos_x <= fly1$pos_x, c((anglearcsin$arcsin1.8*-1)+pi), c(ifelse((fly8$pos_x > fly1$pos_x & fly8$pos_y <= fly1$pos_y), c(anglearcsin$arcsin1.8+(2*pi)), c(anglearcsin$arcsin1.8))))
  angles$ang1.9 <- ifelse(fly9$pos_x <= fly1$pos_x, c((anglearcsin$arcsin1.9*-1)+pi), c(ifelse((fly9$pos_x > fly1$pos_x & fly9$pos_y <= fly1$pos_y), c(anglearcsin$arcsin1.9+(2*pi)), c(anglearcsin$arcsin1.9))))
  angles$ang1.10 <- ifelse(fly10$pos_x <= fly1$pos_x, c((anglearcsin$arcsin1.10*-1)+pi), c(ifelse((fly10$pos_x > fly1$pos_x & fly10$pos_y <= fly1$pos_y), c(anglearcsin$arcsin1.10+(2*pi)), c(anglearcsin$arcsin1.10))))
  angles$ang1.11 <- ifelse(fly11$pos_x <= fly1$pos_x, c((anglearcsin$arcsin1.11*-1)+pi), c(ifelse((fly11$pos_x > fly1$pos_x & fly11$pos_y <= fly1$pos_y), c(anglearcsin$arcsin1.11+(2*pi)), c(anglearcsin$arcsin1.11))))
  angles$ang1.12 <- ifelse(fly12$pos_x <= fly1$pos_x, c((anglearcsin$arcsin1.12*-1)+pi), c(ifelse((fly12$pos_x > fly1$pos_x & fly12$pos_y <= fly1$pos_y), c(anglearcsin$arcsin1.12+(2*pi)), c(anglearcsin$arcsin1.12))))
  angles$ang1.13 <- ifelse(fly13$pos_x <= fly1$pos_x, c((anglearcsin$arcsin1.13*-1)+pi), c(ifelse((fly13$pos_x > fly1$pos_x & fly13$pos_y <= fly1$pos_y), c(anglearcsin$arcsin1.13+(2*pi)), c(anglearcsin$arcsin1.13))))
  angles$ang1.14 <- ifelse(fly14$pos_x <= fly1$pos_x, c((anglearcsin$arcsin1.14*-1)+pi), c(ifelse((fly14$pos_x > fly1$pos_x & fly14$pos_y <= fly1$pos_y), c(anglearcsin$arcsin1.14+(2*pi)), c(anglearcsin$arcsin1.14))))
  angles$ang1.15 <- ifelse(fly15$pos_x <= fly1$pos_x, c((anglearcsin$arcsin1.15*-1)+pi), c(ifelse((fly15$pos_x > fly1$pos_x & fly15$pos_y <= fly1$pos_y), c(anglearcsin$arcsin1.15+(2*pi)), c(anglearcsin$arcsin1.15))))
  angles$ang1.16 <- ifelse(fly16$pos_x <= fly1$pos_x, c((anglearcsin$arcsin1.16*-1)+pi), c(ifelse((fly16$pos_x > fly1$pos_x & fly16$pos_y <= fly1$pos_y), c(anglearcsin$arcsin1.16+(2*pi)), c(anglearcsin$arcsin1.16))))
  angles$ang1.17 <- ifelse(fly17$pos_x <= fly1$pos_x, c((anglearcsin$arcsin1.17*-1)+pi), c(ifelse((fly17$pos_x > fly1$pos_x & fly17$pos_y <= fly1$pos_y), c(anglearcsin$arcsin1.17+(2*pi)), c(anglearcsin$arcsin1.17))))
  angles$ang1.18 <- ifelse(fly18$pos_x <= fly1$pos_x, c((anglearcsin$arcsin1.18*-1)+pi), c(ifelse((fly18$pos_x > fly1$pos_x & fly18$pos_y <= fly1$pos_y), c(anglearcsin$arcsin1.18+(2*pi)), c(anglearcsin$arcsin1.18))))
  angles$ang1.19 <- ifelse(fly19$pos_x <= fly1$pos_x, c((anglearcsin$arcsin1.19*-1)+pi), c(ifelse((fly19$pos_x > fly1$pos_x & fly19$pos_y <= fly1$pos_y), c(anglearcsin$arcsin1.19+(2*pi)), c(anglearcsin$arcsin1.19))))
  angles$ang1.20 <- ifelse(fly20$pos_x <= fly1$pos_x, c((anglearcsin$arcsin1.20*-1)+pi), c(ifelse((fly20$pos_x > fly1$pos_x & fly20$pos_y <= fly1$pos_y), c(anglearcsin$arcsin1.20+(2*pi)), c(anglearcsin$arcsin1.20))))
  angles$ang2.1 <- ifelse(fly1$pos_x <= fly2$pos_x, c((anglearcsin$arcsin2.1*-1)+pi), c(ifelse((fly1$pos_x > fly2$pos_x & fly1$pos_y <= fly2$pos_y), c(anglearcsin$arcsin2.1+(2*pi)), c(anglearcsin$arcsin2.1))))
  angles$ang2.3 <- ifelse(fly3$pos_x <= fly2$pos_x, c((anglearcsin$arcsin2.3*-1)+pi), c(ifelse((fly3$pos_x > fly2$pos_x & fly3$pos_y <= fly2$pos_y), c(anglearcsin$arcsin2.3+(2*pi)), c(anglearcsin$arcsin2.3))))
  angles$ang2.4 <- ifelse(fly4$pos_x <= fly2$pos_x, c((anglearcsin$arcsin2.4*-1)+pi), c(ifelse((fly4$pos_x > fly2$pos_x & fly4$pos_y <= fly2$pos_y), c(anglearcsin$arcsin2.4+(2*pi)), c(anglearcsin$arcsin2.4))))
  angles$ang2.5 <- ifelse(fly5$pos_x <= fly2$pos_x, c((anglearcsin$arcsin2.5*-1)+pi), c(ifelse((fly5$pos_x > fly2$pos_x & fly5$pos_y <= fly2$pos_y), c(anglearcsin$arcsin2.5+(2*pi)), c(anglearcsin$arcsin2.5))))
  angles$ang2.6 <- ifelse(fly6$pos_x <= fly2$pos_x, c((anglearcsin$arcsin2.6*-1)+pi), c(ifelse((fly6$pos_x > fly2$pos_x & fly6$pos_y <= fly2$pos_y), c(anglearcsin$arcsin2.6+(2*pi)), c(anglearcsin$arcsin2.6))))
  angles$ang2.7 <- ifelse(fly7$pos_x <= fly2$pos_x, c((anglearcsin$arcsin2.7*-1)+pi), c(ifelse((fly7$pos_x > fly2$pos_x & fly7$pos_y <= fly2$pos_y), c(anglearcsin$arcsin2.7+(2*pi)), c(anglearcsin$arcsin2.7))))
  angles$ang2.8 <- ifelse(fly8$pos_x <= fly2$pos_x, c((anglearcsin$arcsin2.8*-1)+pi), c(ifelse((fly8$pos_x > fly2$pos_x & fly8$pos_y <= fly2$pos_y), c(anglearcsin$arcsin2.8+(2*pi)), c(anglearcsin$arcsin2.8))))
  angles$ang2.9 <- ifelse(fly9$pos_x <= fly2$pos_x, c((anglearcsin$arcsin2.9*-1)+pi), c(ifelse((fly9$pos_x > fly2$pos_x & fly9$pos_y <= fly2$pos_y), c(anglearcsin$arcsin2.9+(2*pi)), c(anglearcsin$arcsin2.9))))
  angles$ang2.10 <- ifelse(fly10$pos_x <= fly2$pos_x, c((anglearcsin$arcsin2.10*-1)+pi), c(ifelse((fly10$pos_x > fly2$pos_x & fly10$pos_y <= fly2$pos_y), c(anglearcsin$arcsin2.10+(2*pi)), c(anglearcsin$arcsin2.10))))
  angles$ang2.11 <- ifelse(fly11$pos_x <= fly2$pos_x, c((anglearcsin$arcsin2.11*-1)+pi), c(ifelse((fly11$pos_x > fly2$pos_x & fly11$pos_y <= fly2$pos_y), c(anglearcsin$arcsin2.11+(2*pi)), c(anglearcsin$arcsin2.11))))
  angles$ang2.12 <- ifelse(fly12$pos_x <= fly2$pos_x, c((anglearcsin$arcsin2.12*-1)+pi), c(ifelse((fly12$pos_x > fly2$pos_x & fly12$pos_y <= fly2$pos_y), c(anglearcsin$arcsin2.12+(2*pi)), c(anglearcsin$arcsin2.12))))
  angles$ang2.13 <- ifelse(fly13$pos_x <= fly2$pos_x, c((anglearcsin$arcsin2.13*-1)+pi), c(ifelse((fly13$pos_x > fly2$pos_x & fly13$pos_y <= fly2$pos_y), c(anglearcsin$arcsin2.13+(2*pi)), c(anglearcsin$arcsin2.13))))
  angles$ang2.14 <- ifelse(fly14$pos_x <= fly2$pos_x, c((anglearcsin$arcsin2.14*-1)+pi), c(ifelse((fly14$pos_x > fly2$pos_x & fly14$pos_y <= fly2$pos_y), c(anglearcsin$arcsin2.14+(2*pi)), c(anglearcsin$arcsin2.14))))
  angles$ang2.15 <- ifelse(fly15$pos_x <= fly2$pos_x, c((anglearcsin$arcsin2.15*-1)+pi), c(ifelse((fly15$pos_x > fly2$pos_x & fly15$pos_y <= fly2$pos_y), c(anglearcsin$arcsin2.15+(2*pi)), c(anglearcsin$arcsin2.15))))
  angles$ang2.16 <- ifelse(fly16$pos_x <= fly2$pos_x, c((anglearcsin$arcsin2.16*-1)+pi), c(ifelse((fly16$pos_x > fly2$pos_x & fly16$pos_y <= fly2$pos_y), c(anglearcsin$arcsin2.16+(2*pi)), c(anglearcsin$arcsin2.16))))
  angles$ang2.17 <- ifelse(fly17$pos_x <= fly2$pos_x, c((anglearcsin$arcsin2.17*-1)+pi), c(ifelse((fly17$pos_x > fly2$pos_x & fly17$pos_y <= fly2$pos_y), c(anglearcsin$arcsin2.17+(2*pi)), c(anglearcsin$arcsin2.17))))
  angles$ang2.18 <- ifelse(fly18$pos_x <= fly2$pos_x, c((anglearcsin$arcsin2.18*-1)+pi), c(ifelse((fly18$pos_x > fly2$pos_x & fly18$pos_y <= fly2$pos_y), c(anglearcsin$arcsin2.18+(2*pi)), c(anglearcsin$arcsin2.18))))
  angles$ang2.19 <- ifelse(fly19$pos_x <= fly2$pos_x, c((anglearcsin$arcsin2.19*-1)+pi), c(ifelse((fly19$pos_x > fly2$pos_x & fly19$pos_y <= fly2$pos_y), c(anglearcsin$arcsin2.19+(2*pi)), c(anglearcsin$arcsin2.19))))
  angles$ang2.20 <- ifelse(fly20$pos_x <= fly2$pos_x, c((anglearcsin$arcsin2.20*-1)+pi), c(ifelse((fly20$pos_x > fly2$pos_x & fly20$pos_y <= fly2$pos_y), c(anglearcsin$arcsin2.20+(2*pi)), c(anglearcsin$arcsin2.20))))
  angles$ang3.1 <- ifelse(fly1$pos_x <= fly3$pos_x, c((anglearcsin$arcsin3.1*-1)+pi), c(ifelse((fly1$pos_x > fly3$pos_x & fly1$pos_y <= fly3$pos_y), c(anglearcsin$arcsin3.1+(2*pi)), c(anglearcsin$arcsin3.1))))
  angles$ang3.2 <- ifelse(fly2$pos_x <= fly3$pos_x, c((anglearcsin$arcsin3.2*-1)+pi), c(ifelse((fly2$pos_x > fly3$pos_x & fly2$pos_y <= fly3$pos_y), c(anglearcsin$arcsin3.2+(2*pi)), c(anglearcsin$arcsin3.2))))
  angles$ang3.4 <- ifelse(fly4$pos_x <= fly3$pos_x, c((anglearcsin$arcsin3.4*-1)+pi), c(ifelse((fly4$pos_x > fly3$pos_x & fly4$pos_y <= fly3$pos_y), c(anglearcsin$arcsin3.4+(2*pi)), c(anglearcsin$arcsin3.4))))
  angles$ang3.5 <- ifelse(fly5$pos_x <= fly3$pos_x, c((anglearcsin$arcsin3.5*-1)+pi), c(ifelse((fly5$pos_x > fly3$pos_x & fly5$pos_y <= fly3$pos_y), c(anglearcsin$arcsin3.5+(2*pi)), c(anglearcsin$arcsin3.5))))
  angles$ang3.6 <- ifelse(fly6$pos_x <= fly3$pos_x, c((anglearcsin$arcsin3.6*-1)+pi), c(ifelse((fly6$pos_x > fly3$pos_x & fly6$pos_y <= fly3$pos_y), c(anglearcsin$arcsin3.6+(2*pi)), c(anglearcsin$arcsin3.6))))
  angles$ang3.7 <- ifelse(fly7$pos_x <= fly3$pos_x, c((anglearcsin$arcsin3.7*-1)+pi), c(ifelse((fly7$pos_x > fly3$pos_x & fly7$pos_y <= fly3$pos_y), c(anglearcsin$arcsin3.7+(2*pi)), c(anglearcsin$arcsin3.7))))
  angles$ang3.8 <- ifelse(fly8$pos_x <= fly3$pos_x, c((anglearcsin$arcsin3.8*-1)+pi), c(ifelse((fly8$pos_x > fly3$pos_x & fly8$pos_y <= fly3$pos_y), c(anglearcsin$arcsin3.8+(2*pi)), c(anglearcsin$arcsin3.8))))
  angles$ang3.9 <- ifelse(fly9$pos_x <= fly3$pos_x, c((anglearcsin$arcsin3.9*-1)+pi), c(ifelse((fly9$pos_x > fly3$pos_x & fly9$pos_y <= fly3$pos_y), c(anglearcsin$arcsin3.9+(2*pi)), c(anglearcsin$arcsin3.9))))
  angles$ang3.10 <- ifelse(fly10$pos_x <= fly3$pos_x, c((anglearcsin$arcsin3.10*-1)+pi), c(ifelse((fly10$pos_x > fly3$pos_x & fly10$pos_y <= fly3$pos_y), c(anglearcsin$arcsin3.10+(2*pi)), c(anglearcsin$arcsin3.10))))
  angles$ang3.11 <- ifelse(fly11$pos_x <= fly3$pos_x, c((anglearcsin$arcsin3.11*-1)+pi), c(ifelse((fly11$pos_x > fly3$pos_x & fly11$pos_y <= fly3$pos_y), c(anglearcsin$arcsin3.11+(2*pi)), c(anglearcsin$arcsin3.11))))
  angles$ang3.12 <- ifelse(fly12$pos_x <= fly3$pos_x, c((anglearcsin$arcsin3.12*-1)+pi), c(ifelse((fly12$pos_x > fly3$pos_x & fly12$pos_y <= fly3$pos_y), c(anglearcsin$arcsin3.12+(2*pi)), c(anglearcsin$arcsin3.12))))
  angles$ang3.13 <- ifelse(fly13$pos_x <= fly3$pos_x, c((anglearcsin$arcsin3.13*-1)+pi), c(ifelse((fly13$pos_x > fly3$pos_x & fly13$pos_y <= fly3$pos_y), c(anglearcsin$arcsin3.13+(2*pi)), c(anglearcsin$arcsin3.13))))
  angles$ang3.14 <- ifelse(fly14$pos_x <= fly3$pos_x, c((anglearcsin$arcsin3.14*-1)+pi), c(ifelse((fly14$pos_x > fly3$pos_x & fly14$pos_y <= fly3$pos_y), c(anglearcsin$arcsin3.14+(2*pi)), c(anglearcsin$arcsin3.14))))
  angles$ang3.15 <- ifelse(fly15$pos_x <= fly3$pos_x, c((anglearcsin$arcsin3.15*-1)+pi), c(ifelse((fly15$pos_x > fly3$pos_x & fly15$pos_y <= fly3$pos_y), c(anglearcsin$arcsin3.15+(2*pi)), c(anglearcsin$arcsin3.15))))
  angles$ang3.16 <- ifelse(fly16$pos_x <= fly3$pos_x, c((anglearcsin$arcsin3.16*-1)+pi), c(ifelse((fly16$pos_x > fly3$pos_x & fly16$pos_y <= fly3$pos_y), c(anglearcsin$arcsin3.16+(2*pi)), c(anglearcsin$arcsin3.16))))
  angles$ang3.17 <- ifelse(fly17$pos_x <= fly3$pos_x, c((anglearcsin$arcsin3.17*-1)+pi), c(ifelse((fly17$pos_x > fly3$pos_x & fly17$pos_y <= fly3$pos_y), c(anglearcsin$arcsin3.17+(2*pi)), c(anglearcsin$arcsin3.17))))
  angles$ang3.18 <- ifelse(fly18$pos_x <= fly3$pos_x, c((anglearcsin$arcsin3.18*-1)+pi), c(ifelse((fly18$pos_x > fly3$pos_x & fly18$pos_y <= fly3$pos_y), c(anglearcsin$arcsin3.18+(2*pi)), c(anglearcsin$arcsin3.18))))
  angles$ang3.19 <- ifelse(fly19$pos_x <= fly3$pos_x, c((anglearcsin$arcsin3.19*-1)+pi), c(ifelse((fly19$pos_x > fly3$pos_x & fly19$pos_y <= fly3$pos_y), c(anglearcsin$arcsin3.19+(2*pi)), c(anglearcsin$arcsin3.19))))
  angles$ang3.20 <- ifelse(fly20$pos_x <= fly3$pos_x, c((anglearcsin$arcsin3.20*-1)+pi), c(ifelse((fly20$pos_x > fly3$pos_x & fly20$pos_y <= fly3$pos_y), c(anglearcsin$arcsin3.20+(2*pi)), c(anglearcsin$arcsin3.20))))
  angles$ang4.1 <- ifelse(fly1$pos_x <= fly4$pos_x, c((anglearcsin$arcsin4.1*-1)+pi), c(ifelse((fly1$pos_x > fly4$pos_x & fly1$pos_y <= fly4$pos_y), c(anglearcsin$arcsin4.1+(2*pi)), c(anglearcsin$arcsin4.1))))
  angles$ang4.2 <- ifelse(fly2$pos_x <= fly4$pos_x, c((anglearcsin$arcsin4.2*-1)+pi), c(ifelse((fly2$pos_x > fly4$pos_x & fly2$pos_y <= fly4$pos_y), c(anglearcsin$arcsin4.2+(2*pi)), c(anglearcsin$arcsin4.2))))
  angles$ang4.3 <- ifelse(fly3$pos_x <= fly4$pos_x, c((anglearcsin$arcsin4.3*-1)+pi), c(ifelse((fly3$pos_x > fly4$pos_x & fly3$pos_y <= fly4$pos_y), c(anglearcsin$arcsin4.3+(2*pi)), c(anglearcsin$arcsin4.3))))
  angles$ang4.5 <- ifelse(fly5$pos_x <= fly4$pos_x, c((anglearcsin$arcsin4.5*-1)+pi), c(ifelse((fly5$pos_x > fly4$pos_x & fly5$pos_y <= fly4$pos_y), c(anglearcsin$arcsin4.5+(2*pi)), c(anglearcsin$arcsin4.5))))
  angles$ang4.6 <- ifelse(fly6$pos_x <= fly4$pos_x, c((anglearcsin$arcsin4.6*-1)+pi), c(ifelse((fly6$pos_x > fly4$pos_x & fly6$pos_y <= fly4$pos_y), c(anglearcsin$arcsin4.6+(2*pi)), c(anglearcsin$arcsin4.6))))
  angles$ang4.7 <- ifelse(fly7$pos_x <= fly4$pos_x, c((anglearcsin$arcsin4.7*-1)+pi), c(ifelse((fly7$pos_x > fly4$pos_x & fly7$pos_y <= fly4$pos_y), c(anglearcsin$arcsin4.7+(2*pi)), c(anglearcsin$arcsin4.7))))
  angles$ang4.8 <- ifelse(fly8$pos_x <= fly4$pos_x, c((anglearcsin$arcsin4.8*-1)+pi), c(ifelse((fly8$pos_x > fly4$pos_x & fly8$pos_y <= fly4$pos_y), c(anglearcsin$arcsin4.8+(2*pi)), c(anglearcsin$arcsin4.8))))
  angles$ang4.9 <- ifelse(fly9$pos_x <= fly4$pos_x, c((anglearcsin$arcsin4.9*-1)+pi), c(ifelse((fly9$pos_x > fly4$pos_x & fly9$pos_y <= fly4$pos_y), c(anglearcsin$arcsin4.9+(2*pi)), c(anglearcsin$arcsin4.9))))
  angles$ang4.10 <- ifelse(fly10$pos_x <= fly4$pos_x, c((anglearcsin$arcsin4.10*-1)+pi), c(ifelse((fly10$pos_x > fly4$pos_x & fly10$pos_y <= fly4$pos_y), c(anglearcsin$arcsin4.10+(2*pi)), c(anglearcsin$arcsin4.10))))
  angles$ang4.11 <- ifelse(fly11$pos_x <= fly4$pos_x, c((anglearcsin$arcsin4.11*-1)+pi), c(ifelse((fly11$pos_x > fly4$pos_x & fly11$pos_y <= fly4$pos_y), c(anglearcsin$arcsin4.11+(2*pi)), c(anglearcsin$arcsin4.11))))
  angles$ang4.12 <- ifelse(fly12$pos_x <= fly4$pos_x, c((anglearcsin$arcsin4.12*-1)+pi), c(ifelse((fly12$pos_x > fly4$pos_x & fly12$pos_y <= fly4$pos_y), c(anglearcsin$arcsin4.12+(2*pi)), c(anglearcsin$arcsin4.12))))
  angles$ang4.13 <- ifelse(fly13$pos_x <= fly4$pos_x, c((anglearcsin$arcsin4.13*-1)+pi), c(ifelse((fly13$pos_x > fly4$pos_x & fly13$pos_y <= fly4$pos_y), c(anglearcsin$arcsin4.13+(2*pi)), c(anglearcsin$arcsin4.13))))
  angles$ang4.14 <- ifelse(fly14$pos_x <= fly4$pos_x, c((anglearcsin$arcsin4.14*-1)+pi), c(ifelse((fly14$pos_x > fly4$pos_x & fly14$pos_y <= fly4$pos_y), c(anglearcsin$arcsin4.14+(2*pi)), c(anglearcsin$arcsin4.14))))
  angles$ang4.15 <- ifelse(fly15$pos_x <= fly4$pos_x, c((anglearcsin$arcsin4.15*-1)+pi), c(ifelse((fly15$pos_x > fly4$pos_x & fly15$pos_y <= fly4$pos_y), c(anglearcsin$arcsin4.15+(2*pi)), c(anglearcsin$arcsin4.15))))
  angles$ang4.16 <- ifelse(fly16$pos_x <= fly4$pos_x, c((anglearcsin$arcsin4.16*-1)+pi), c(ifelse((fly16$pos_x > fly4$pos_x & fly16$pos_y <= fly4$pos_y), c(anglearcsin$arcsin4.16+(2*pi)), c(anglearcsin$arcsin4.16))))
  angles$ang4.17 <- ifelse(fly17$pos_x <= fly4$pos_x, c((anglearcsin$arcsin4.17*-1)+pi), c(ifelse((fly17$pos_x > fly4$pos_x & fly17$pos_y <= fly4$pos_y), c(anglearcsin$arcsin4.17+(2*pi)), c(anglearcsin$arcsin4.17))))
  angles$ang4.18 <- ifelse(fly18$pos_x <= fly4$pos_x, c((anglearcsin$arcsin4.18*-1)+pi), c(ifelse((fly18$pos_x > fly4$pos_x & fly18$pos_y <= fly4$pos_y), c(anglearcsin$arcsin4.18+(2*pi)), c(anglearcsin$arcsin4.18))))
  angles$ang4.19 <- ifelse(fly19$pos_x <= fly4$pos_x, c((anglearcsin$arcsin4.19*-1)+pi), c(ifelse((fly19$pos_x > fly4$pos_x & fly19$pos_y <= fly4$pos_y), c(anglearcsin$arcsin4.19+(2*pi)), c(anglearcsin$arcsin4.19))))
  angles$ang4.20 <- ifelse(fly20$pos_x <= fly4$pos_x, c((anglearcsin$arcsin4.20*-1)+pi), c(ifelse((fly20$pos_x > fly4$pos_x & fly20$pos_y <= fly4$pos_y), c(anglearcsin$arcsin4.20+(2*pi)), c(anglearcsin$arcsin4.20))))
  angles$ang5.1 <- ifelse(fly1$pos_x <= fly5$pos_x, c((anglearcsin$arcsin5.1*-1)+pi), c(ifelse((fly1$pos_x > fly5$pos_x & fly1$pos_y <= fly5$pos_y), c(anglearcsin$arcsin5.1+(2*pi)), c(anglearcsin$arcsin5.1))))
  angles$ang5.2 <- ifelse(fly2$pos_x <= fly5$pos_x, c((anglearcsin$arcsin5.2*-1)+pi), c(ifelse((fly2$pos_x > fly5$pos_x & fly2$pos_y <= fly5$pos_y), c(anglearcsin$arcsin5.2+(2*pi)), c(anglearcsin$arcsin5.2))))
  angles$ang5.3 <- ifelse(fly3$pos_x <= fly5$pos_x, c((anglearcsin$arcsin5.3*-1)+pi), c(ifelse((fly3$pos_x > fly5$pos_x & fly3$pos_y <= fly5$pos_y), c(anglearcsin$arcsin5.3+(2*pi)), c(anglearcsin$arcsin5.3))))
  angles$ang5.4 <- ifelse(fly4$pos_x <= fly5$pos_x, c((anglearcsin$arcsin5.4*-1)+pi), c(ifelse((fly4$pos_x > fly5$pos_x & fly4$pos_y <= fly5$pos_y), c(anglearcsin$arcsin5.4+(2*pi)), c(anglearcsin$arcsin5.4))))
  angles$ang5.6 <- ifelse(fly6$pos_x <= fly5$pos_x, c((anglearcsin$arcsin5.6*-1)+pi), c(ifelse((fly6$pos_x > fly5$pos_x & fly6$pos_y <= fly5$pos_y), c(anglearcsin$arcsin5.6+(2*pi)), c(anglearcsin$arcsin5.6))))
  angles$ang5.7 <- ifelse(fly7$pos_x <= fly5$pos_x, c((anglearcsin$arcsin5.7*-1)+pi), c(ifelse((fly7$pos_x > fly5$pos_x & fly7$pos_y <= fly5$pos_y), c(anglearcsin$arcsin5.7+(2*pi)), c(anglearcsin$arcsin5.7))))
  angles$ang5.8 <- ifelse(fly8$pos_x <= fly5$pos_x, c((anglearcsin$arcsin5.8*-1)+pi), c(ifelse((fly8$pos_x > fly5$pos_x & fly8$pos_y <= fly5$pos_y), c(anglearcsin$arcsin5.8+(2*pi)), c(anglearcsin$arcsin5.8))))
  angles$ang5.9 <- ifelse(fly9$pos_x <= fly5$pos_x, c((anglearcsin$arcsin5.9*-1)+pi), c(ifelse((fly9$pos_x > fly5$pos_x & fly9$pos_y <= fly5$pos_y), c(anglearcsin$arcsin5.9+(2*pi)), c(anglearcsin$arcsin5.9))))
  angles$ang5.10 <- ifelse(fly10$pos_x <= fly5$pos_x, c((anglearcsin$arcsin5.10*-1)+pi), c(ifelse((fly10$pos_x > fly5$pos_x & fly10$pos_y <= fly5$pos_y), c(anglearcsin$arcsin5.10+(2*pi)), c(anglearcsin$arcsin5.10))))
  angles$ang5.11 <- ifelse(fly11$pos_x <= fly5$pos_x, c((anglearcsin$arcsin5.11*-1)+pi), c(ifelse((fly11$pos_x > fly5$pos_x & fly11$pos_y <= fly5$pos_y), c(anglearcsin$arcsin5.11+(2*pi)), c(anglearcsin$arcsin5.11))))
  angles$ang5.12 <- ifelse(fly12$pos_x <= fly5$pos_x, c((anglearcsin$arcsin5.12*-1)+pi), c(ifelse((fly12$pos_x > fly5$pos_x & fly12$pos_y <= fly5$pos_y), c(anglearcsin$arcsin5.12+(2*pi)), c(anglearcsin$arcsin5.12))))
  angles$ang5.13 <- ifelse(fly13$pos_x <= fly5$pos_x, c((anglearcsin$arcsin5.13*-1)+pi), c(ifelse((fly13$pos_x > fly5$pos_x & fly13$pos_y <= fly5$pos_y), c(anglearcsin$arcsin5.13+(2*pi)), c(anglearcsin$arcsin5.13))))
  angles$ang5.14 <- ifelse(fly14$pos_x <= fly5$pos_x, c((anglearcsin$arcsin5.14*-1)+pi), c(ifelse((fly14$pos_x > fly5$pos_x & fly14$pos_y <= fly5$pos_y), c(anglearcsin$arcsin5.14+(2*pi)), c(anglearcsin$arcsin5.14))))
  angles$ang5.15 <- ifelse(fly15$pos_x <= fly5$pos_x, c((anglearcsin$arcsin5.15*-1)+pi), c(ifelse((fly15$pos_x > fly5$pos_x & fly15$pos_y <= fly5$pos_y), c(anglearcsin$arcsin5.15+(2*pi)), c(anglearcsin$arcsin5.15))))
  angles$ang5.16 <- ifelse(fly16$pos_x <= fly5$pos_x, c((anglearcsin$arcsin5.16*-1)+pi), c(ifelse((fly16$pos_x > fly5$pos_x & fly16$pos_y <= fly5$pos_y), c(anglearcsin$arcsin5.16+(2*pi)), c(anglearcsin$arcsin5.16))))
  angles$ang5.17 <- ifelse(fly17$pos_x <= fly5$pos_x, c((anglearcsin$arcsin5.17*-1)+pi), c(ifelse((fly17$pos_x > fly5$pos_x & fly17$pos_y <= fly5$pos_y), c(anglearcsin$arcsin5.17+(2*pi)), c(anglearcsin$arcsin5.17))))
  angles$ang5.18 <- ifelse(fly18$pos_x <= fly5$pos_x, c((anglearcsin$arcsin5.18*-1)+pi), c(ifelse((fly18$pos_x > fly5$pos_x & fly18$pos_y <= fly5$pos_y), c(anglearcsin$arcsin5.18+(2*pi)), c(anglearcsin$arcsin5.18))))
  angles$ang5.19 <- ifelse(fly19$pos_x <= fly5$pos_x, c((anglearcsin$arcsin5.19*-1)+pi), c(ifelse((fly19$pos_x > fly5$pos_x & fly19$pos_y <= fly5$pos_y), c(anglearcsin$arcsin5.19+(2*pi)), c(anglearcsin$arcsin5.19))))
  angles$ang5.20 <- ifelse(fly20$pos_x <= fly5$pos_x, c((anglearcsin$arcsin5.20*-1)+pi), c(ifelse((fly20$pos_x > fly5$pos_x & fly20$pos_y <= fly5$pos_y), c(anglearcsin$arcsin5.20+(2*pi)), c(anglearcsin$arcsin5.20))))
  angles$ang6.1 <- ifelse(fly1$pos_x <= fly6$pos_x, c((anglearcsin$arcsin6.1*-1)+pi), c(ifelse((fly1$pos_x > fly6$pos_x & fly1$pos_y <= fly6$pos_y), c(anglearcsin$arcsin6.1+(2*pi)), c(anglearcsin$arcsin6.1))))
  angles$ang6.2 <- ifelse(fly2$pos_x <= fly6$pos_x, c((anglearcsin$arcsin6.2*-1)+pi), c(ifelse((fly2$pos_x > fly6$pos_x & fly2$pos_y <= fly6$pos_y), c(anglearcsin$arcsin6.2+(2*pi)), c(anglearcsin$arcsin6.2))))
  angles$ang6.3 <- ifelse(fly3$pos_x <= fly6$pos_x, c((anglearcsin$arcsin6.3*-1)+pi), c(ifelse((fly3$pos_x > fly6$pos_x & fly3$pos_y <= fly6$pos_y), c(anglearcsin$arcsin6.3+(2*pi)), c(anglearcsin$arcsin6.3))))
  angles$ang6.4 <- ifelse(fly4$pos_x <= fly6$pos_x, c((anglearcsin$arcsin6.4*-1)+pi), c(ifelse((fly4$pos_x > fly6$pos_x & fly4$pos_y <= fly6$pos_y), c(anglearcsin$arcsin6.4+(2*pi)), c(anglearcsin$arcsin6.4))))
  angles$ang6.5 <- ifelse(fly5$pos_x <= fly6$pos_x, c((anglearcsin$arcsin6.5*-1)+pi), c(ifelse((fly5$pos_x > fly6$pos_x & fly5$pos_y <= fly6$pos_y), c(anglearcsin$arcsin6.5+(2*pi)), c(anglearcsin$arcsin6.5))))
  angles$ang6.7 <- ifelse(fly7$pos_x <= fly6$pos_x, c((anglearcsin$arcsin6.7*-1)+pi), c(ifelse((fly7$pos_x > fly6$pos_x & fly7$pos_y <= fly6$pos_y), c(anglearcsin$arcsin6.7+(2*pi)), c(anglearcsin$arcsin6.7))))
  angles$ang6.8 <- ifelse(fly8$pos_x <= fly6$pos_x, c((anglearcsin$arcsin6.8*-1)+pi), c(ifelse((fly8$pos_x > fly6$pos_x & fly8$pos_y <= fly6$pos_y), c(anglearcsin$arcsin6.8+(2*pi)), c(anglearcsin$arcsin6.8))))
  angles$ang6.9 <- ifelse(fly9$pos_x <= fly6$pos_x, c((anglearcsin$arcsin6.9*-1)+pi), c(ifelse((fly9$pos_x > fly6$pos_x & fly9$pos_y <= fly6$pos_y), c(anglearcsin$arcsin6.9+(2*pi)), c(anglearcsin$arcsin6.9))))
  angles$ang6.10 <- ifelse(fly10$pos_x <= fly6$pos_x, c((anglearcsin$arcsin6.10*-1)+pi), c(ifelse((fly10$pos_x > fly6$pos_x & fly10$pos_y <= fly6$pos_y), c(anglearcsin$arcsin6.10+(2*pi)), c(anglearcsin$arcsin6.10))))
  angles$ang6.11 <- ifelse(fly11$pos_x <= fly6$pos_x, c((anglearcsin$arcsin6.11*-1)+pi), c(ifelse((fly11$pos_x > fly6$pos_x & fly11$pos_y <= fly6$pos_y), c(anglearcsin$arcsin6.11+(2*pi)), c(anglearcsin$arcsin6.11))))
  angles$ang6.12 <- ifelse(fly12$pos_x <= fly6$pos_x, c((anglearcsin$arcsin6.12*-1)+pi), c(ifelse((fly12$pos_x > fly6$pos_x & fly12$pos_y <= fly6$pos_y), c(anglearcsin$arcsin6.12+(2*pi)), c(anglearcsin$arcsin6.12))))
  angles$ang6.13 <- ifelse(fly13$pos_x <= fly6$pos_x, c((anglearcsin$arcsin6.13*-1)+pi), c(ifelse((fly13$pos_x > fly6$pos_x & fly13$pos_y <= fly6$pos_y), c(anglearcsin$arcsin6.13+(2*pi)), c(anglearcsin$arcsin6.13))))
  angles$ang6.14 <- ifelse(fly14$pos_x <= fly6$pos_x, c((anglearcsin$arcsin6.14*-1)+pi), c(ifelse((fly14$pos_x > fly6$pos_x & fly14$pos_y <= fly6$pos_y), c(anglearcsin$arcsin6.14+(2*pi)), c(anglearcsin$arcsin6.14))))
  angles$ang6.15 <- ifelse(fly15$pos_x <= fly6$pos_x, c((anglearcsin$arcsin6.15*-1)+pi), c(ifelse((fly15$pos_x > fly6$pos_x & fly15$pos_y <= fly6$pos_y), c(anglearcsin$arcsin6.15+(2*pi)), c(anglearcsin$arcsin6.15))))
  angles$ang6.16 <- ifelse(fly16$pos_x <= fly6$pos_x, c((anglearcsin$arcsin6.16*-1)+pi), c(ifelse((fly16$pos_x > fly6$pos_x & fly16$pos_y <= fly6$pos_y), c(anglearcsin$arcsin6.16+(2*pi)), c(anglearcsin$arcsin6.16))))
  angles$ang6.17 <- ifelse(fly17$pos_x <= fly6$pos_x, c((anglearcsin$arcsin6.17*-1)+pi), c(ifelse((fly17$pos_x > fly6$pos_x & fly17$pos_y <= fly6$pos_y), c(anglearcsin$arcsin6.17+(2*pi)), c(anglearcsin$arcsin6.17))))
  angles$ang6.18 <- ifelse(fly18$pos_x <= fly6$pos_x, c((anglearcsin$arcsin6.18*-1)+pi), c(ifelse((fly18$pos_x > fly6$pos_x & fly18$pos_y <= fly6$pos_y), c(anglearcsin$arcsin6.18+(2*pi)), c(anglearcsin$arcsin6.18))))
  angles$ang6.19 <- ifelse(fly19$pos_x <= fly6$pos_x, c((anglearcsin$arcsin6.19*-1)+pi), c(ifelse((fly19$pos_x > fly6$pos_x & fly19$pos_y <= fly6$pos_y), c(anglearcsin$arcsin6.19+(2*pi)), c(anglearcsin$arcsin6.19))))
  angles$ang6.20 <- ifelse(fly20$pos_x <= fly6$pos_x, c((anglearcsin$arcsin6.20*-1)+pi), c(ifelse((fly20$pos_x > fly6$pos_x & fly20$pos_y <= fly6$pos_y), c(anglearcsin$arcsin6.20+(2*pi)), c(anglearcsin$arcsin6.20))))
  angles$ang7.1 <- ifelse(fly1$pos_x <= fly7$pos_x, c((anglearcsin$arcsin7.1*-1)+pi), c(ifelse((fly1$pos_x > fly7$pos_x & fly1$pos_y <= fly7$pos_y), c(anglearcsin$arcsin7.1+(2*pi)), c(anglearcsin$arcsin7.1))))
  angles$ang7.2 <- ifelse(fly2$pos_x <= fly7$pos_x, c((anglearcsin$arcsin7.2*-1)+pi), c(ifelse((fly2$pos_x > fly7$pos_x & fly2$pos_y <= fly7$pos_y), c(anglearcsin$arcsin7.2+(2*pi)), c(anglearcsin$arcsin7.2))))
  angles$ang7.3 <- ifelse(fly3$pos_x <= fly7$pos_x, c((anglearcsin$arcsin7.3*-1)+pi), c(ifelse((fly3$pos_x > fly7$pos_x & fly3$pos_y <= fly7$pos_y), c(anglearcsin$arcsin7.3+(2*pi)), c(anglearcsin$arcsin7.3))))
  angles$ang7.4 <- ifelse(fly4$pos_x <= fly7$pos_x, c((anglearcsin$arcsin7.4*-1)+pi), c(ifelse((fly4$pos_x > fly7$pos_x & fly4$pos_y <= fly7$pos_y), c(anglearcsin$arcsin7.4+(2*pi)), c(anglearcsin$arcsin7.4))))
  angles$ang7.5 <- ifelse(fly5$pos_x <= fly7$pos_x, c((anglearcsin$arcsin7.5*-1)+pi), c(ifelse((fly5$pos_x > fly7$pos_x & fly5$pos_y <= fly7$pos_y), c(anglearcsin$arcsin7.5+(2*pi)), c(anglearcsin$arcsin7.5))))
  angles$ang7.6 <- ifelse(fly6$pos_x <= fly7$pos_x, c((anglearcsin$arcsin7.6*-1)+pi), c(ifelse((fly6$pos_x > fly7$pos_x & fly6$pos_y <= fly7$pos_y), c(anglearcsin$arcsin7.6+(2*pi)), c(anglearcsin$arcsin7.6))))
  angles$ang7.8 <- ifelse(fly8$pos_x <= fly7$pos_x, c((anglearcsin$arcsin7.8*-1)+pi), c(ifelse((fly8$pos_x > fly7$pos_x & fly8$pos_y <= fly7$pos_y), c(anglearcsin$arcsin7.8+(2*pi)), c(anglearcsin$arcsin7.8))))
  angles$ang7.9 <- ifelse(fly9$pos_x <= fly7$pos_x, c((anglearcsin$arcsin7.9*-1)+pi), c(ifelse((fly9$pos_x > fly7$pos_x & fly9$pos_y <= fly7$pos_y), c(anglearcsin$arcsin7.9+(2*pi)), c(anglearcsin$arcsin7.9))))
  angles$ang7.10 <- ifelse(fly10$pos_x <= fly7$pos_x, c((anglearcsin$arcsin7.10*-1)+pi), c(ifelse((fly10$pos_x > fly7$pos_x & fly10$pos_y <= fly7$pos_y), c(anglearcsin$arcsin7.10+(2*pi)), c(anglearcsin$arcsin7.10))))
  angles$ang7.11 <- ifelse(fly11$pos_x <= fly7$pos_x, c((anglearcsin$arcsin7.11*-1)+pi), c(ifelse((fly11$pos_x > fly7$pos_x & fly11$pos_y <= fly7$pos_y), c(anglearcsin$arcsin7.11+(2*pi)), c(anglearcsin$arcsin7.11))))
  angles$ang7.12 <- ifelse(fly12$pos_x <= fly7$pos_x, c((anglearcsin$arcsin7.12*-1)+pi), c(ifelse((fly12$pos_x > fly7$pos_x & fly12$pos_y <= fly7$pos_y), c(anglearcsin$arcsin7.12+(2*pi)), c(anglearcsin$arcsin7.12))))
  angles$ang7.13 <- ifelse(fly13$pos_x <= fly7$pos_x, c((anglearcsin$arcsin7.13*-1)+pi), c(ifelse((fly13$pos_x > fly7$pos_x & fly13$pos_y <= fly7$pos_y), c(anglearcsin$arcsin7.13+(2*pi)), c(anglearcsin$arcsin7.13))))
  angles$ang7.14 <- ifelse(fly14$pos_x <= fly7$pos_x, c((anglearcsin$arcsin7.14*-1)+pi), c(ifelse((fly14$pos_x > fly7$pos_x & fly14$pos_y <= fly7$pos_y), c(anglearcsin$arcsin7.14+(2*pi)), c(anglearcsin$arcsin7.14))))
  angles$ang7.15 <- ifelse(fly15$pos_x <= fly7$pos_x, c((anglearcsin$arcsin7.15*-1)+pi), c(ifelse((fly15$pos_x > fly7$pos_x & fly15$pos_y <= fly7$pos_y), c(anglearcsin$arcsin7.15+(2*pi)), c(anglearcsin$arcsin7.15))))
  angles$ang7.16 <- ifelse(fly16$pos_x <= fly7$pos_x, c((anglearcsin$arcsin7.16*-1)+pi), c(ifelse((fly16$pos_x > fly7$pos_x & fly16$pos_y <= fly7$pos_y), c(anglearcsin$arcsin7.16+(2*pi)), c(anglearcsin$arcsin7.16))))
  angles$ang7.17 <- ifelse(fly17$pos_x <= fly7$pos_x, c((anglearcsin$arcsin7.17*-1)+pi), c(ifelse((fly17$pos_x > fly7$pos_x & fly17$pos_y <= fly7$pos_y), c(anglearcsin$arcsin7.17+(2*pi)), c(anglearcsin$arcsin7.17))))
  angles$ang7.18 <- ifelse(fly18$pos_x <= fly7$pos_x, c((anglearcsin$arcsin7.18*-1)+pi), c(ifelse((fly18$pos_x > fly7$pos_x & fly18$pos_y <= fly7$pos_y), c(anglearcsin$arcsin7.18+(2*pi)), c(anglearcsin$arcsin7.18))))
  angles$ang7.19 <- ifelse(fly19$pos_x <= fly7$pos_x, c((anglearcsin$arcsin7.19*-1)+pi), c(ifelse((fly19$pos_x > fly7$pos_x & fly19$pos_y <= fly7$pos_y), c(anglearcsin$arcsin7.19+(2*pi)), c(anglearcsin$arcsin7.19))))
  angles$ang7.20 <- ifelse(fly20$pos_x <= fly7$pos_x, c((anglearcsin$arcsin7.20*-1)+pi), c(ifelse((fly20$pos_x > fly7$pos_x & fly20$pos_y <= fly7$pos_y), c(anglearcsin$arcsin7.20+(2*pi)), c(anglearcsin$arcsin7.20))))
  angles$ang8.1 <- ifelse(fly1$pos_x <= fly8$pos_x, c((anglearcsin$arcsin8.1*-1)+pi), c(ifelse((fly1$pos_x > fly8$pos_x & fly1$pos_y <= fly8$pos_y), c(anglearcsin$arcsin8.1+(2*pi)), c(anglearcsin$arcsin8.1))))
  angles$ang8.2 <- ifelse(fly2$pos_x <= fly8$pos_x, c((anglearcsin$arcsin8.2*-1)+pi), c(ifelse((fly2$pos_x > fly8$pos_x & fly2$pos_y <= fly8$pos_y), c(anglearcsin$arcsin8.2+(2*pi)), c(anglearcsin$arcsin8.2))))
  angles$ang8.3 <- ifelse(fly3$pos_x <= fly8$pos_x, c((anglearcsin$arcsin8.3*-1)+pi), c(ifelse((fly3$pos_x > fly8$pos_x & fly3$pos_y <= fly8$pos_y), c(anglearcsin$arcsin8.3+(2*pi)), c(anglearcsin$arcsin8.3))))
  angles$ang8.4 <- ifelse(fly4$pos_x <= fly8$pos_x, c((anglearcsin$arcsin8.4*-1)+pi), c(ifelse((fly4$pos_x > fly8$pos_x & fly4$pos_y <= fly8$pos_y), c(anglearcsin$arcsin8.4+(2*pi)), c(anglearcsin$arcsin8.4))))
  angles$ang8.5 <- ifelse(fly5$pos_x <= fly8$pos_x, c((anglearcsin$arcsin8.5*-1)+pi), c(ifelse((fly5$pos_x > fly8$pos_x & fly5$pos_y <= fly8$pos_y), c(anglearcsin$arcsin8.5+(2*pi)), c(anglearcsin$arcsin8.5))))
  angles$ang8.6 <- ifelse(fly6$pos_x <= fly8$pos_x, c((anglearcsin$arcsin8.6*-1)+pi), c(ifelse((fly6$pos_x > fly8$pos_x & fly6$pos_y <= fly8$pos_y), c(anglearcsin$arcsin8.6+(2*pi)), c(anglearcsin$arcsin8.6))))
  angles$ang8.7 <- ifelse(fly7$pos_x <= fly8$pos_x, c((anglearcsin$arcsin8.7*-1)+pi), c(ifelse((fly7$pos_x > fly8$pos_x & fly7$pos_y <= fly8$pos_y), c(anglearcsin$arcsin8.7+(2*pi)), c(anglearcsin$arcsin8.7))))
  angles$ang8.9 <- ifelse(fly9$pos_x <= fly8$pos_x, c((anglearcsin$arcsin8.9*-1)+pi), c(ifelse((fly9$pos_x > fly8$pos_x & fly9$pos_y <= fly8$pos_y), c(anglearcsin$arcsin8.9+(2*pi)), c(anglearcsin$arcsin8.9))))
  angles$ang8.10 <- ifelse(fly10$pos_x <= fly8$pos_x, c((anglearcsin$arcsin8.10*-1)+pi), c(ifelse((fly10$pos_x > fly8$pos_x & fly10$pos_y <= fly8$pos_y), c(anglearcsin$arcsin8.10+(2*pi)), c(anglearcsin$arcsin8.10))))
  angles$ang8.11 <- ifelse(fly11$pos_x <= fly8$pos_x, c((anglearcsin$arcsin8.11*-1)+pi), c(ifelse((fly11$pos_x > fly8$pos_x & fly11$pos_y <= fly8$pos_y), c(anglearcsin$arcsin8.11+(2*pi)), c(anglearcsin$arcsin8.11))))
  angles$ang8.12 <- ifelse(fly12$pos_x <= fly8$pos_x, c((anglearcsin$arcsin8.12*-1)+pi), c(ifelse((fly12$pos_x > fly8$pos_x & fly12$pos_y <= fly8$pos_y), c(anglearcsin$arcsin8.12+(2*pi)), c(anglearcsin$arcsin8.12))))
  angles$ang8.13 <- ifelse(fly13$pos_x <= fly8$pos_x, c((anglearcsin$arcsin8.13*-1)+pi), c(ifelse((fly13$pos_x > fly8$pos_x & fly13$pos_y <= fly8$pos_y), c(anglearcsin$arcsin8.13+(2*pi)), c(anglearcsin$arcsin8.13))))
  angles$ang8.14 <- ifelse(fly14$pos_x <= fly8$pos_x, c((anglearcsin$arcsin8.14*-1)+pi), c(ifelse((fly14$pos_x > fly8$pos_x & fly14$pos_y <= fly8$pos_y), c(anglearcsin$arcsin8.14+(2*pi)), c(anglearcsin$arcsin8.14))))
  angles$ang8.15 <- ifelse(fly15$pos_x <= fly8$pos_x, c((anglearcsin$arcsin8.15*-1)+pi), c(ifelse((fly15$pos_x > fly8$pos_x & fly15$pos_y <= fly8$pos_y), c(anglearcsin$arcsin8.15+(2*pi)), c(anglearcsin$arcsin8.15))))
  angles$ang8.16 <- ifelse(fly16$pos_x <= fly8$pos_x, c((anglearcsin$arcsin8.16*-1)+pi), c(ifelse((fly16$pos_x > fly8$pos_x & fly16$pos_y <= fly8$pos_y), c(anglearcsin$arcsin8.16+(2*pi)), c(anglearcsin$arcsin8.16))))
  angles$ang8.17 <- ifelse(fly17$pos_x <= fly8$pos_x, c((anglearcsin$arcsin8.17*-1)+pi), c(ifelse((fly17$pos_x > fly8$pos_x & fly17$pos_y <= fly8$pos_y), c(anglearcsin$arcsin8.17+(2*pi)), c(anglearcsin$arcsin8.17))))
  angles$ang8.18 <- ifelse(fly18$pos_x <= fly8$pos_x, c((anglearcsin$arcsin8.18*-1)+pi), c(ifelse((fly18$pos_x > fly8$pos_x & fly18$pos_y <= fly8$pos_y), c(anglearcsin$arcsin8.18+(2*pi)), c(anglearcsin$arcsin8.18))))
  angles$ang8.19 <- ifelse(fly19$pos_x <= fly8$pos_x, c((anglearcsin$arcsin8.19*-1)+pi), c(ifelse((fly19$pos_x > fly8$pos_x & fly19$pos_y <= fly8$pos_y), c(anglearcsin$arcsin8.19+(2*pi)), c(anglearcsin$arcsin8.19))))
  angles$ang8.20 <- ifelse(fly20$pos_x <= fly8$pos_x, c((anglearcsin$arcsin8.20*-1)+pi), c(ifelse((fly20$pos_x > fly8$pos_x & fly20$pos_y <= fly8$pos_y), c(anglearcsin$arcsin8.20+(2*pi)), c(anglearcsin$arcsin8.20))))
  angles$ang9.1 <- ifelse(fly1$pos_x <= fly9$pos_x, c((anglearcsin$arcsin9.1*-1)+pi), c(ifelse((fly1$pos_x > fly9$pos_x & fly1$pos_y <= fly9$pos_y), c(anglearcsin$arcsin9.1+(2*pi)), c(anglearcsin$arcsin9.1))))
  angles$ang9.2 <- ifelse(fly2$pos_x <= fly9$pos_x, c((anglearcsin$arcsin9.2*-1)+pi), c(ifelse((fly2$pos_x > fly9$pos_x & fly2$pos_y <= fly9$pos_y), c(anglearcsin$arcsin9.2+(2*pi)), c(anglearcsin$arcsin9.2))))
  angles$ang9.3 <- ifelse(fly3$pos_x <= fly9$pos_x, c((anglearcsin$arcsin9.3*-1)+pi), c(ifelse((fly3$pos_x > fly9$pos_x & fly3$pos_y <= fly9$pos_y), c(anglearcsin$arcsin9.3+(2*pi)), c(anglearcsin$arcsin9.3))))
  angles$ang9.4 <- ifelse(fly4$pos_x <= fly9$pos_x, c((anglearcsin$arcsin9.4*-1)+pi), c(ifelse((fly4$pos_x > fly9$pos_x & fly4$pos_y <= fly9$pos_y), c(anglearcsin$arcsin9.4+(2*pi)), c(anglearcsin$arcsin9.4))))
  angles$ang9.5 <- ifelse(fly5$pos_x <= fly9$pos_x, c((anglearcsin$arcsin9.5*-1)+pi), c(ifelse((fly5$pos_x > fly9$pos_x & fly5$pos_y <= fly9$pos_y), c(anglearcsin$arcsin9.5+(2*pi)), c(anglearcsin$arcsin9.5))))
  angles$ang9.6 <- ifelse(fly6$pos_x <= fly9$pos_x, c((anglearcsin$arcsin9.6*-1)+pi), c(ifelse((fly6$pos_x > fly9$pos_x & fly6$pos_y <= fly9$pos_y), c(anglearcsin$arcsin9.6+(2*pi)), c(anglearcsin$arcsin9.6))))
  angles$ang9.7 <- ifelse(fly7$pos_x <= fly9$pos_x, c((anglearcsin$arcsin9.7*-1)+pi), c(ifelse((fly7$pos_x > fly9$pos_x & fly7$pos_y <= fly9$pos_y), c(anglearcsin$arcsin9.7+(2*pi)), c(anglearcsin$arcsin9.7))))
  angles$ang9.8 <- ifelse(fly8$pos_x <= fly9$pos_x, c((anglearcsin$arcsin9.8*-1)+pi), c(ifelse((fly8$pos_x > fly9$pos_x & fly8$pos_y <= fly9$pos_y), c(anglearcsin$arcsin9.8+(2*pi)), c(anglearcsin$arcsin9.8))))
  angles$ang9.10 <- ifelse(fly10$pos_x <= fly9$pos_x, c((anglearcsin$arcsin9.10*-1)+pi), c(ifelse((fly10$pos_x > fly9$pos_x & fly10$pos_y <= fly9$pos_y), c(anglearcsin$arcsin9.10+(2*pi)), c(anglearcsin$arcsin9.10))))
  angles$ang9.11 <- ifelse(fly11$pos_x <= fly9$pos_x, c((anglearcsin$arcsin9.11*-1)+pi), c(ifelse((fly11$pos_x > fly9$pos_x & fly11$pos_y <= fly9$pos_y), c(anglearcsin$arcsin9.11+(2*pi)), c(anglearcsin$arcsin9.11))))
  angles$ang9.12 <- ifelse(fly12$pos_x <= fly9$pos_x, c((anglearcsin$arcsin9.12*-1)+pi), c(ifelse((fly12$pos_x > fly9$pos_x & fly12$pos_y <= fly9$pos_y), c(anglearcsin$arcsin9.12+(2*pi)), c(anglearcsin$arcsin9.12))))
  angles$ang9.13 <- ifelse(fly13$pos_x <= fly9$pos_x, c((anglearcsin$arcsin9.13*-1)+pi), c(ifelse((fly13$pos_x > fly9$pos_x & fly13$pos_y <= fly9$pos_y), c(anglearcsin$arcsin9.13+(2*pi)), c(anglearcsin$arcsin9.13))))
  angles$ang9.14 <- ifelse(fly14$pos_x <= fly9$pos_x, c((anglearcsin$arcsin9.14*-1)+pi), c(ifelse((fly14$pos_x > fly9$pos_x & fly14$pos_y <= fly9$pos_y), c(anglearcsin$arcsin9.14+(2*pi)), c(anglearcsin$arcsin9.14))))
  angles$ang9.15 <- ifelse(fly15$pos_x <= fly9$pos_x, c((anglearcsin$arcsin9.15*-1)+pi), c(ifelse((fly15$pos_x > fly9$pos_x & fly15$pos_y <= fly9$pos_y), c(anglearcsin$arcsin9.15+(2*pi)), c(anglearcsin$arcsin9.15))))
  angles$ang9.16 <- ifelse(fly16$pos_x <= fly9$pos_x, c((anglearcsin$arcsin9.16*-1)+pi), c(ifelse((fly16$pos_x > fly9$pos_x & fly16$pos_y <= fly9$pos_y), c(anglearcsin$arcsin9.16+(2*pi)), c(anglearcsin$arcsin9.16))))
  angles$ang9.17 <- ifelse(fly17$pos_x <= fly9$pos_x, c((anglearcsin$arcsin9.17*-1)+pi), c(ifelse((fly17$pos_x > fly9$pos_x & fly17$pos_y <= fly9$pos_y), c(anglearcsin$arcsin9.17+(2*pi)), c(anglearcsin$arcsin9.17))))
  angles$ang9.18 <- ifelse(fly18$pos_x <= fly9$pos_x, c((anglearcsin$arcsin9.18*-1)+pi), c(ifelse((fly18$pos_x > fly9$pos_x & fly18$pos_y <= fly9$pos_y), c(anglearcsin$arcsin9.18+(2*pi)), c(anglearcsin$arcsin9.18))))
  angles$ang9.19 <- ifelse(fly19$pos_x <= fly9$pos_x, c((anglearcsin$arcsin9.19*-1)+pi), c(ifelse((fly19$pos_x > fly9$pos_x & fly19$pos_y <= fly9$pos_y), c(anglearcsin$arcsin9.19+(2*pi)), c(anglearcsin$arcsin9.19))))
  angles$ang9.20 <- ifelse(fly20$pos_x <= fly9$pos_x, c((anglearcsin$arcsin9.20*-1)+pi), c(ifelse((fly20$pos_x > fly9$pos_x & fly20$pos_y <= fly9$pos_y), c(anglearcsin$arcsin9.20+(2*pi)), c(anglearcsin$arcsin9.20))))
  angles$ang10.1 <- ifelse(fly1$pos_x <= fly10$pos_x, c((anglearcsin$arcsin10.1*-1)+pi), c(ifelse((fly1$pos_x > fly10$pos_x & fly1$pos_y <= fly10$pos_y), c(anglearcsin$arcsin10.1+(2*pi)), c(anglearcsin$arcsin10.1))))
  angles$ang10.2 <- ifelse(fly2$pos_x <= fly10$pos_x, c((anglearcsin$arcsin10.2*-1)+pi), c(ifelse((fly2$pos_x > fly10$pos_x & fly2$pos_y <= fly10$pos_y), c(anglearcsin$arcsin10.2+(2*pi)), c(anglearcsin$arcsin10.2))))
  angles$ang10.3 <- ifelse(fly3$pos_x <= fly10$pos_x, c((anglearcsin$arcsin10.3*-1)+pi), c(ifelse((fly3$pos_x > fly10$pos_x & fly3$pos_y <= fly10$pos_y), c(anglearcsin$arcsin10.3+(2*pi)), c(anglearcsin$arcsin10.3))))
  angles$ang10.4 <- ifelse(fly4$pos_x <= fly10$pos_x, c((anglearcsin$arcsin10.4*-1)+pi), c(ifelse((fly4$pos_x > fly10$pos_x & fly4$pos_y <= fly10$pos_y), c(anglearcsin$arcsin10.4+(2*pi)), c(anglearcsin$arcsin10.4))))
  angles$ang10.5 <- ifelse(fly5$pos_x <= fly10$pos_x, c((anglearcsin$arcsin10.5*-1)+pi), c(ifelse((fly5$pos_x > fly10$pos_x & fly5$pos_y <= fly10$pos_y), c(anglearcsin$arcsin10.5+(2*pi)), c(anglearcsin$arcsin10.5))))
  angles$ang10.6 <- ifelse(fly6$pos_x <= fly10$pos_x, c((anglearcsin$arcsin10.6*-1)+pi), c(ifelse((fly6$pos_x > fly10$pos_x & fly6$pos_y <= fly10$pos_y), c(anglearcsin$arcsin10.6+(2*pi)), c(anglearcsin$arcsin10.6))))
  angles$ang10.7 <- ifelse(fly7$pos_x <= fly10$pos_x, c((anglearcsin$arcsin10.7*-1)+pi), c(ifelse((fly7$pos_x > fly10$pos_x & fly7$pos_y <= fly10$pos_y), c(anglearcsin$arcsin10.7+(2*pi)), c(anglearcsin$arcsin10.7))))
  angles$ang10.8 <- ifelse(fly8$pos_x <= fly10$pos_x, c((anglearcsin$arcsin10.8*-1)+pi), c(ifelse((fly8$pos_x > fly10$pos_x & fly8$pos_y <= fly10$pos_y), c(anglearcsin$arcsin10.8+(2*pi)), c(anglearcsin$arcsin10.8))))
  angles$ang10.9 <- ifelse(fly9$pos_x <= fly10$pos_x, c((anglearcsin$arcsin10.9*-1)+pi), c(ifelse((fly9$pos_x > fly10$pos_x & fly9$pos_y <= fly10$pos_y), c(anglearcsin$arcsin10.9+(2*pi)), c(anglearcsin$arcsin10.9))))
  angles$ang10.11 <- ifelse(fly11$pos_x <= fly10$pos_x, c((anglearcsin$arcsin10.11*-1)+pi), c(ifelse((fly11$pos_x > fly10$pos_x & fly11$pos_y <= fly10$pos_y), c(anglearcsin$arcsin10.11+(2*pi)), c(anglearcsin$arcsin10.11))))
  angles$ang10.12 <- ifelse(fly12$pos_x <= fly10$pos_x, c((anglearcsin$arcsin10.12*-1)+pi), c(ifelse((fly12$pos_x > fly10$pos_x & fly12$pos_y <= fly10$pos_y), c(anglearcsin$arcsin10.12+(2*pi)), c(anglearcsin$arcsin10.12))))
  angles$ang10.13 <- ifelse(fly13$pos_x <= fly10$pos_x, c((anglearcsin$arcsin10.13*-1)+pi), c(ifelse((fly13$pos_x > fly10$pos_x & fly13$pos_y <= fly10$pos_y), c(anglearcsin$arcsin10.13+(2*pi)), c(anglearcsin$arcsin10.13))))
  angles$ang10.14 <- ifelse(fly14$pos_x <= fly10$pos_x, c((anglearcsin$arcsin10.14*-1)+pi), c(ifelse((fly14$pos_x > fly10$pos_x & fly14$pos_y <= fly10$pos_y), c(anglearcsin$arcsin10.14+(2*pi)), c(anglearcsin$arcsin10.14))))
  angles$ang10.15 <- ifelse(fly15$pos_x <= fly10$pos_x, c((anglearcsin$arcsin10.15*-1)+pi), c(ifelse((fly15$pos_x > fly10$pos_x & fly15$pos_y <= fly10$pos_y), c(anglearcsin$arcsin10.15+(2*pi)), c(anglearcsin$arcsin10.15))))
  angles$ang10.16 <- ifelse(fly16$pos_x <= fly10$pos_x, c((anglearcsin$arcsin10.16*-1)+pi), c(ifelse((fly16$pos_x > fly10$pos_x & fly16$pos_y <= fly10$pos_y), c(anglearcsin$arcsin10.16+(2*pi)), c(anglearcsin$arcsin10.16))))
  angles$ang10.17 <- ifelse(fly17$pos_x <= fly10$pos_x, c((anglearcsin$arcsin10.17*-1)+pi), c(ifelse((fly17$pos_x > fly10$pos_x & fly17$pos_y <= fly10$pos_y), c(anglearcsin$arcsin10.17+(2*pi)), c(anglearcsin$arcsin10.17))))
  angles$ang10.18 <- ifelse(fly18$pos_x <= fly10$pos_x, c((anglearcsin$arcsin10.18*-1)+pi), c(ifelse((fly18$pos_x > fly10$pos_x & fly18$pos_y <= fly10$pos_y), c(anglearcsin$arcsin10.18+(2*pi)), c(anglearcsin$arcsin10.18))))
  angles$ang10.19 <- ifelse(fly19$pos_x <= fly10$pos_x, c((anglearcsin$arcsin10.19*-1)+pi), c(ifelse((fly19$pos_x > fly10$pos_x & fly19$pos_y <= fly10$pos_y), c(anglearcsin$arcsin10.19+(2*pi)), c(anglearcsin$arcsin10.19))))
  angles$ang10.20 <- ifelse(fly20$pos_x <= fly10$pos_x, c((anglearcsin$arcsin10.20*-1)+pi), c(ifelse((fly20$pos_x > fly10$pos_x & fly20$pos_y <= fly10$pos_y), c(anglearcsin$arcsin10.20+(2*pi)), c(anglearcsin$arcsin10.20))))
  angles$ang11.1 <- ifelse(fly1$pos_x <= fly11$pos_x, c((anglearcsin$arcsin11.1*-1)+pi), c(ifelse((fly1$pos_x > fly11$pos_x & fly1$pos_y <= fly11$pos_y), c(anglearcsin$arcsin11.1+(2*pi)), c(anglearcsin$arcsin11.1))))
  angles$ang11.2 <- ifelse(fly2$pos_x <= fly11$pos_x, c((anglearcsin$arcsin11.2*-1)+pi), c(ifelse((fly2$pos_x > fly11$pos_x & fly2$pos_y <= fly11$pos_y), c(anglearcsin$arcsin11.2+(2*pi)), c(anglearcsin$arcsin11.2))))
  angles$ang11.3 <- ifelse(fly3$pos_x <= fly11$pos_x, c((anglearcsin$arcsin11.3*-1)+pi), c(ifelse((fly3$pos_x > fly11$pos_x & fly3$pos_y <= fly11$pos_y), c(anglearcsin$arcsin11.3+(2*pi)), c(anglearcsin$arcsin11.3))))
  angles$ang11.4 <- ifelse(fly4$pos_x <= fly11$pos_x, c((anglearcsin$arcsin11.4*-1)+pi), c(ifelse((fly4$pos_x > fly11$pos_x & fly4$pos_y <= fly11$pos_y), c(anglearcsin$arcsin11.4+(2*pi)), c(anglearcsin$arcsin11.4))))
  angles$ang11.5 <- ifelse(fly5$pos_x <= fly11$pos_x, c((anglearcsin$arcsin11.5*-1)+pi), c(ifelse((fly5$pos_x > fly11$pos_x & fly5$pos_y <= fly11$pos_y), c(anglearcsin$arcsin11.5+(2*pi)), c(anglearcsin$arcsin11.5))))
  angles$ang11.6 <- ifelse(fly6$pos_x <= fly11$pos_x, c((anglearcsin$arcsin11.6*-1)+pi), c(ifelse((fly6$pos_x > fly11$pos_x & fly6$pos_y <= fly11$pos_y), c(anglearcsin$arcsin11.6+(2*pi)), c(anglearcsin$arcsin11.6))))
  angles$ang11.7 <- ifelse(fly7$pos_x <= fly11$pos_x, c((anglearcsin$arcsin11.7*-1)+pi), c(ifelse((fly7$pos_x > fly11$pos_x & fly7$pos_y <= fly11$pos_y), c(anglearcsin$arcsin11.7+(2*pi)), c(anglearcsin$arcsin11.7))))
  angles$ang11.8 <- ifelse(fly8$pos_x <= fly11$pos_x, c((anglearcsin$arcsin11.8*-1)+pi), c(ifelse((fly8$pos_x > fly11$pos_x & fly8$pos_y <= fly11$pos_y), c(anglearcsin$arcsin11.8+(2*pi)), c(anglearcsin$arcsin11.8))))
  angles$ang11.9 <- ifelse(fly9$pos_x <= fly11$pos_x, c((anglearcsin$arcsin11.9*-1)+pi), c(ifelse((fly9$pos_x > fly11$pos_x & fly9$pos_y <= fly11$pos_y), c(anglearcsin$arcsin11.9+(2*pi)), c(anglearcsin$arcsin11.9))))
  angles$ang11.10 <- ifelse(fly10$pos_x <= fly11$pos_x, c((anglearcsin$arcsin11.10*-1)+pi), c(ifelse((fly10$pos_x > fly11$pos_x & fly10$pos_y <= fly11$pos_y), c(anglearcsin$arcsin11.10+(2*pi)), c(anglearcsin$arcsin11.10))))
  angles$ang11.12 <- ifelse(fly12$pos_x <= fly11$pos_x, c((anglearcsin$arcsin11.12*-1)+pi), c(ifelse((fly12$pos_x > fly11$pos_x & fly12$pos_y <= fly11$pos_y), c(anglearcsin$arcsin11.12+(2*pi)), c(anglearcsin$arcsin11.12))))
  angles$ang11.13 <- ifelse(fly13$pos_x <= fly11$pos_x, c((anglearcsin$arcsin11.13*-1)+pi), c(ifelse((fly13$pos_x > fly11$pos_x & fly13$pos_y <= fly11$pos_y), c(anglearcsin$arcsin11.13+(2*pi)), c(anglearcsin$arcsin11.13))))
  angles$ang11.14 <- ifelse(fly14$pos_x <= fly11$pos_x, c((anglearcsin$arcsin11.14*-1)+pi), c(ifelse((fly14$pos_x > fly11$pos_x & fly14$pos_y <= fly11$pos_y), c(anglearcsin$arcsin11.14+(2*pi)), c(anglearcsin$arcsin11.14))))
  angles$ang11.15 <- ifelse(fly15$pos_x <= fly11$pos_x, c((anglearcsin$arcsin11.15*-1)+pi), c(ifelse((fly15$pos_x > fly11$pos_x & fly15$pos_y <= fly11$pos_y), c(anglearcsin$arcsin11.15+(2*pi)), c(anglearcsin$arcsin11.15))))
  angles$ang11.16 <- ifelse(fly16$pos_x <= fly11$pos_x, c((anglearcsin$arcsin11.16*-1)+pi), c(ifelse((fly16$pos_x > fly11$pos_x & fly16$pos_y <= fly11$pos_y), c(anglearcsin$arcsin11.16+(2*pi)), c(anglearcsin$arcsin11.16))))
  angles$ang11.17 <- ifelse(fly17$pos_x <= fly11$pos_x, c((anglearcsin$arcsin11.17*-1)+pi), c(ifelse((fly17$pos_x > fly11$pos_x & fly17$pos_y <= fly11$pos_y), c(anglearcsin$arcsin11.17+(2*pi)), c(anglearcsin$arcsin11.17))))
  angles$ang11.18 <- ifelse(fly18$pos_x <= fly11$pos_x, c((anglearcsin$arcsin11.18*-1)+pi), c(ifelse((fly18$pos_x > fly11$pos_x & fly18$pos_y <= fly11$pos_y), c(anglearcsin$arcsin11.18+(2*pi)), c(anglearcsin$arcsin11.18))))
  angles$ang11.19 <- ifelse(fly19$pos_x <= fly11$pos_x, c((anglearcsin$arcsin11.19*-1)+pi), c(ifelse((fly19$pos_x > fly11$pos_x & fly19$pos_y <= fly11$pos_y), c(anglearcsin$arcsin11.19+(2*pi)), c(anglearcsin$arcsin11.19))))
  angles$ang11.20 <- ifelse(fly20$pos_x <= fly11$pos_x, c((anglearcsin$arcsin11.20*-1)+pi), c(ifelse((fly20$pos_x > fly11$pos_x & fly20$pos_y <= fly11$pos_y), c(anglearcsin$arcsin11.20+(2*pi)), c(anglearcsin$arcsin11.20))))
  angles$ang12.1 <- ifelse(fly1$pos_x <= fly12$pos_x, c((anglearcsin$arcsin12.1*-1)+pi), c(ifelse((fly1$pos_x > fly12$pos_x & fly1$pos_y <= fly12$pos_y), c(anglearcsin$arcsin12.1+(2*pi)), c(anglearcsin$arcsin12.1))))
  angles$ang12.2 <- ifelse(fly2$pos_x <= fly12$pos_x, c((anglearcsin$arcsin12.2*-1)+pi), c(ifelse((fly2$pos_x > fly12$pos_x & fly2$pos_y <= fly12$pos_y), c(anglearcsin$arcsin12.2+(2*pi)), c(anglearcsin$arcsin12.2))))
  angles$ang12.3 <- ifelse(fly3$pos_x <= fly12$pos_x, c((anglearcsin$arcsin12.3*-1)+pi), c(ifelse((fly3$pos_x > fly12$pos_x & fly3$pos_y <= fly12$pos_y), c(anglearcsin$arcsin12.3+(2*pi)), c(anglearcsin$arcsin12.3))))
  angles$ang12.4 <- ifelse(fly4$pos_x <= fly12$pos_x, c((anglearcsin$arcsin12.4*-1)+pi), c(ifelse((fly4$pos_x > fly12$pos_x & fly4$pos_y <= fly12$pos_y), c(anglearcsin$arcsin12.4+(2*pi)), c(anglearcsin$arcsin12.4))))
  angles$ang12.5 <- ifelse(fly5$pos_x <= fly12$pos_x, c((anglearcsin$arcsin12.5*-1)+pi), c(ifelse((fly5$pos_x > fly12$pos_x & fly5$pos_y <= fly12$pos_y), c(anglearcsin$arcsin12.5+(2*pi)), c(anglearcsin$arcsin12.5))))
  angles$ang12.6 <- ifelse(fly6$pos_x <= fly12$pos_x, c((anglearcsin$arcsin12.6*-1)+pi), c(ifelse((fly6$pos_x > fly12$pos_x & fly6$pos_y <= fly12$pos_y), c(anglearcsin$arcsin12.6+(2*pi)), c(anglearcsin$arcsin12.6))))
  angles$ang12.7 <- ifelse(fly7$pos_x <= fly12$pos_x, c((anglearcsin$arcsin12.7*-1)+pi), c(ifelse((fly7$pos_x > fly12$pos_x & fly7$pos_y <= fly12$pos_y), c(anglearcsin$arcsin12.7+(2*pi)), c(anglearcsin$arcsin12.7))))
  angles$ang12.8 <- ifelse(fly8$pos_x <= fly12$pos_x, c((anglearcsin$arcsin12.8*-1)+pi), c(ifelse((fly8$pos_x > fly12$pos_x & fly8$pos_y <= fly12$pos_y), c(anglearcsin$arcsin12.8+(2*pi)), c(anglearcsin$arcsin12.8))))
  angles$ang12.9 <- ifelse(fly9$pos_x <= fly12$pos_x, c((anglearcsin$arcsin12.9*-1)+pi), c(ifelse((fly9$pos_x > fly12$pos_x & fly9$pos_y <= fly12$pos_y), c(anglearcsin$arcsin12.9+(2*pi)), c(anglearcsin$arcsin12.9))))
  angles$ang12.10 <- ifelse(fly10$pos_x <= fly12$pos_x, c((anglearcsin$arcsin12.10*-1)+pi), c(ifelse((fly10$pos_x > fly12$pos_x & fly10$pos_y <= fly12$pos_y), c(anglearcsin$arcsin12.10+(2*pi)), c(anglearcsin$arcsin12.10))))
  angles$ang12.11 <- ifelse(fly11$pos_x <= fly12$pos_x, c((anglearcsin$arcsin12.11*-1)+pi), c(ifelse((fly11$pos_x > fly12$pos_x & fly11$pos_y <= fly12$pos_y), c(anglearcsin$arcsin12.11+(2*pi)), c(anglearcsin$arcsin12.11))))
  angles$ang12.13 <- ifelse(fly13$pos_x <= fly12$pos_x, c((anglearcsin$arcsin12.13*-1)+pi), c(ifelse((fly13$pos_x > fly12$pos_x & fly13$pos_y <= fly12$pos_y), c(anglearcsin$arcsin12.13+(2*pi)), c(anglearcsin$arcsin12.13))))
  angles$ang12.14 <- ifelse(fly14$pos_x <= fly12$pos_x, c((anglearcsin$arcsin12.14*-1)+pi), c(ifelse((fly14$pos_x > fly12$pos_x & fly14$pos_y <= fly12$pos_y), c(anglearcsin$arcsin12.14+(2*pi)), c(anglearcsin$arcsin12.14))))
  angles$ang12.15 <- ifelse(fly15$pos_x <= fly12$pos_x, c((anglearcsin$arcsin12.15*-1)+pi), c(ifelse((fly15$pos_x > fly12$pos_x & fly15$pos_y <= fly12$pos_y), c(anglearcsin$arcsin12.15+(2*pi)), c(anglearcsin$arcsin12.15))))
  angles$ang12.16 <- ifelse(fly16$pos_x <= fly12$pos_x, c((anglearcsin$arcsin12.16*-1)+pi), c(ifelse((fly16$pos_x > fly12$pos_x & fly16$pos_y <= fly12$pos_y), c(anglearcsin$arcsin12.16+(2*pi)), c(anglearcsin$arcsin12.16))))
  angles$ang12.17 <- ifelse(fly17$pos_x <= fly12$pos_x, c((anglearcsin$arcsin12.17*-1)+pi), c(ifelse((fly17$pos_x > fly12$pos_x & fly17$pos_y <= fly12$pos_y), c(anglearcsin$arcsin12.17+(2*pi)), c(anglearcsin$arcsin12.17))))
  angles$ang12.18 <- ifelse(fly18$pos_x <= fly12$pos_x, c((anglearcsin$arcsin12.18*-1)+pi), c(ifelse((fly18$pos_x > fly12$pos_x & fly18$pos_y <= fly12$pos_y), c(anglearcsin$arcsin12.18+(2*pi)), c(anglearcsin$arcsin12.18))))
  angles$ang12.19 <- ifelse(fly19$pos_x <= fly12$pos_x, c((anglearcsin$arcsin12.19*-1)+pi), c(ifelse((fly19$pos_x > fly12$pos_x & fly19$pos_y <= fly12$pos_y), c(anglearcsin$arcsin12.19+(2*pi)), c(anglearcsin$arcsin12.19))))
  angles$ang12.20 <- ifelse(fly20$pos_x <= fly12$pos_x, c((anglearcsin$arcsin12.20*-1)+pi), c(ifelse((fly20$pos_x > fly12$pos_x & fly20$pos_y <= fly12$pos_y), c(anglearcsin$arcsin12.20+(2*pi)), c(anglearcsin$arcsin12.20))))
  angles$ang13.1 <- ifelse(fly1$pos_x <= fly13$pos_x, c((anglearcsin$arcsin13.1*-1)+pi), c(ifelse((fly1$pos_x > fly13$pos_x & fly1$pos_y <= fly13$pos_y), c(anglearcsin$arcsin13.1+(2*pi)), c(anglearcsin$arcsin13.1))))
  angles$ang13.2 <- ifelse(fly2$pos_x <= fly13$pos_x, c((anglearcsin$arcsin13.2*-1)+pi), c(ifelse((fly2$pos_x > fly13$pos_x & fly2$pos_y <= fly13$pos_y), c(anglearcsin$arcsin13.2+(2*pi)), c(anglearcsin$arcsin13.2))))
  angles$ang13.3 <- ifelse(fly3$pos_x <= fly13$pos_x, c((anglearcsin$arcsin13.3*-1)+pi), c(ifelse((fly3$pos_x > fly13$pos_x & fly3$pos_y <= fly13$pos_y), c(anglearcsin$arcsin13.3+(2*pi)), c(anglearcsin$arcsin13.3))))
  angles$ang13.4 <- ifelse(fly4$pos_x <= fly13$pos_x, c((anglearcsin$arcsin13.4*-1)+pi), c(ifelse((fly4$pos_x > fly13$pos_x & fly4$pos_y <= fly13$pos_y), c(anglearcsin$arcsin13.4+(2*pi)), c(anglearcsin$arcsin13.4))))
  angles$ang13.5 <- ifelse(fly5$pos_x <= fly13$pos_x, c((anglearcsin$arcsin13.5*-1)+pi), c(ifelse((fly5$pos_x > fly13$pos_x & fly5$pos_y <= fly13$pos_y), c(anglearcsin$arcsin13.5+(2*pi)), c(anglearcsin$arcsin13.5))))
  angles$ang13.6 <- ifelse(fly6$pos_x <= fly13$pos_x, c((anglearcsin$arcsin13.6*-1)+pi), c(ifelse((fly6$pos_x > fly13$pos_x & fly6$pos_y <= fly13$pos_y), c(anglearcsin$arcsin13.6+(2*pi)), c(anglearcsin$arcsin13.6))))
  angles$ang13.7 <- ifelse(fly7$pos_x <= fly13$pos_x, c((anglearcsin$arcsin13.7*-1)+pi), c(ifelse((fly7$pos_x > fly13$pos_x & fly7$pos_y <= fly13$pos_y), c(anglearcsin$arcsin13.7+(2*pi)), c(anglearcsin$arcsin13.7))))
  angles$ang13.8 <- ifelse(fly8$pos_x <= fly13$pos_x, c((anglearcsin$arcsin13.8*-1)+pi), c(ifelse((fly8$pos_x > fly13$pos_x & fly8$pos_y <= fly13$pos_y), c(anglearcsin$arcsin13.8+(2*pi)), c(anglearcsin$arcsin13.8))))
  angles$ang13.9 <- ifelse(fly9$pos_x <= fly13$pos_x, c((anglearcsin$arcsin13.9*-1)+pi), c(ifelse((fly9$pos_x > fly13$pos_x & fly9$pos_y <= fly13$pos_y), c(anglearcsin$arcsin13.9+(2*pi)), c(anglearcsin$arcsin13.9))))
  angles$ang13.10 <- ifelse(fly10$pos_x <= fly13$pos_x, c((anglearcsin$arcsin13.10*-1)+pi), c(ifelse((fly10$pos_x > fly13$pos_x & fly10$pos_y <= fly13$pos_y), c(anglearcsin$arcsin13.10+(2*pi)), c(anglearcsin$arcsin13.10))))
  angles$ang13.11 <- ifelse(fly11$pos_x <= fly13$pos_x, c((anglearcsin$arcsin13.11*-1)+pi), c(ifelse((fly11$pos_x > fly13$pos_x & fly11$pos_y <= fly13$pos_y), c(anglearcsin$arcsin13.11+(2*pi)), c(anglearcsin$arcsin13.11))))
  angles$ang13.12 <- ifelse(fly12$pos_x <= fly13$pos_x, c((anglearcsin$arcsin13.12*-1)+pi), c(ifelse((fly12$pos_x > fly13$pos_x & fly12$pos_y <= fly13$pos_y), c(anglearcsin$arcsin13.12+(2*pi)), c(anglearcsin$arcsin13.12))))
  angles$ang13.14 <- ifelse(fly14$pos_x <= fly13$pos_x, c((anglearcsin$arcsin13.14*-1)+pi), c(ifelse((fly14$pos_x > fly13$pos_x & fly14$pos_y <= fly13$pos_y), c(anglearcsin$arcsin13.14+(2*pi)), c(anglearcsin$arcsin13.14))))
  angles$ang13.15 <- ifelse(fly15$pos_x <= fly13$pos_x, c((anglearcsin$arcsin13.15*-1)+pi), c(ifelse((fly15$pos_x > fly13$pos_x & fly15$pos_y <= fly13$pos_y), c(anglearcsin$arcsin13.15+(2*pi)), c(anglearcsin$arcsin13.15))))
  angles$ang13.16 <- ifelse(fly16$pos_x <= fly13$pos_x, c((anglearcsin$arcsin13.16*-1)+pi), c(ifelse((fly16$pos_x > fly13$pos_x & fly16$pos_y <= fly13$pos_y), c(anglearcsin$arcsin13.16+(2*pi)), c(anglearcsin$arcsin13.16))))
  angles$ang13.17 <- ifelse(fly17$pos_x <= fly13$pos_x, c((anglearcsin$arcsin13.17*-1)+pi), c(ifelse((fly17$pos_x > fly13$pos_x & fly17$pos_y <= fly13$pos_y), c(anglearcsin$arcsin13.17+(2*pi)), c(anglearcsin$arcsin13.17))))
  angles$ang13.18 <- ifelse(fly18$pos_x <= fly13$pos_x, c((anglearcsin$arcsin13.18*-1)+pi), c(ifelse((fly18$pos_x > fly13$pos_x & fly18$pos_y <= fly13$pos_y), c(anglearcsin$arcsin13.18+(2*pi)), c(anglearcsin$arcsin13.18))))
  angles$ang13.19 <- ifelse(fly19$pos_x <= fly13$pos_x, c((anglearcsin$arcsin13.19*-1)+pi), c(ifelse((fly19$pos_x > fly13$pos_x & fly19$pos_y <= fly13$pos_y), c(anglearcsin$arcsin13.19+(2*pi)), c(anglearcsin$arcsin13.19))))
  angles$ang13.20 <- ifelse(fly20$pos_x <= fly13$pos_x, c((anglearcsin$arcsin13.20*-1)+pi), c(ifelse((fly20$pos_x > fly13$pos_x & fly20$pos_y <= fly13$pos_y), c(anglearcsin$arcsin13.20+(2*pi)), c(anglearcsin$arcsin13.20))))
  angles$ang14.1 <- ifelse(fly1$pos_x <= fly14$pos_x, c((anglearcsin$arcsin14.1*-1)+pi), c(ifelse((fly1$pos_x > fly14$pos_x & fly1$pos_y <= fly14$pos_y), c(anglearcsin$arcsin14.1+(2*pi)), c(anglearcsin$arcsin14.1))))
  angles$ang14.2 <- ifelse(fly2$pos_x <= fly14$pos_x, c((anglearcsin$arcsin14.2*-1)+pi), c(ifelse((fly2$pos_x > fly14$pos_x & fly2$pos_y <= fly14$pos_y), c(anglearcsin$arcsin14.2+(2*pi)), c(anglearcsin$arcsin14.2))))
  angles$ang14.3 <- ifelse(fly3$pos_x <= fly14$pos_x, c((anglearcsin$arcsin14.3*-1)+pi), c(ifelse((fly3$pos_x > fly14$pos_x & fly3$pos_y <= fly14$pos_y), c(anglearcsin$arcsin14.3+(2*pi)), c(anglearcsin$arcsin14.3))))
  angles$ang14.4 <- ifelse(fly4$pos_x <= fly14$pos_x, c((anglearcsin$arcsin14.4*-1)+pi), c(ifelse((fly4$pos_x > fly14$pos_x & fly4$pos_y <= fly14$pos_y), c(anglearcsin$arcsin14.4+(2*pi)), c(anglearcsin$arcsin14.4))))
  angles$ang14.5 <- ifelse(fly5$pos_x <= fly14$pos_x, c((anglearcsin$arcsin14.5*-1)+pi), c(ifelse((fly5$pos_x > fly14$pos_x & fly5$pos_y <= fly14$pos_y), c(anglearcsin$arcsin14.5+(2*pi)), c(anglearcsin$arcsin14.5))))
  angles$ang14.6 <- ifelse(fly6$pos_x <= fly14$pos_x, c((anglearcsin$arcsin14.6*-1)+pi), c(ifelse((fly6$pos_x > fly14$pos_x & fly6$pos_y <= fly14$pos_y), c(anglearcsin$arcsin14.6+(2*pi)), c(anglearcsin$arcsin14.6))))
  angles$ang14.7 <- ifelse(fly7$pos_x <= fly14$pos_x, c((anglearcsin$arcsin14.7*-1)+pi), c(ifelse((fly7$pos_x > fly14$pos_x & fly7$pos_y <= fly14$pos_y), c(anglearcsin$arcsin14.7+(2*pi)), c(anglearcsin$arcsin14.7))))
  angles$ang14.8 <- ifelse(fly8$pos_x <= fly14$pos_x, c((anglearcsin$arcsin14.8*-1)+pi), c(ifelse((fly8$pos_x > fly14$pos_x & fly8$pos_y <= fly14$pos_y), c(anglearcsin$arcsin14.8+(2*pi)), c(anglearcsin$arcsin14.8))))
  angles$ang14.9 <- ifelse(fly9$pos_x <= fly14$pos_x, c((anglearcsin$arcsin14.9*-1)+pi), c(ifelse((fly9$pos_x > fly14$pos_x & fly9$pos_y <= fly14$pos_y), c(anglearcsin$arcsin14.9+(2*pi)), c(anglearcsin$arcsin14.9))))
  angles$ang14.10 <- ifelse(fly10$pos_x <= fly14$pos_x, c((anglearcsin$arcsin14.10*-1)+pi), c(ifelse((fly10$pos_x > fly14$pos_x & fly10$pos_y <= fly14$pos_y), c(anglearcsin$arcsin14.10+(2*pi)), c(anglearcsin$arcsin14.10))))
  angles$ang14.11 <- ifelse(fly11$pos_x <= fly14$pos_x, c((anglearcsin$arcsin14.11*-1)+pi), c(ifelse((fly11$pos_x > fly14$pos_x & fly11$pos_y <= fly14$pos_y), c(anglearcsin$arcsin14.11+(2*pi)), c(anglearcsin$arcsin14.11))))
  angles$ang14.12 <- ifelse(fly12$pos_x <= fly14$pos_x, c((anglearcsin$arcsin14.12*-1)+pi), c(ifelse((fly12$pos_x > fly14$pos_x & fly12$pos_y <= fly14$pos_y), c(anglearcsin$arcsin14.12+(2*pi)), c(anglearcsin$arcsin14.12))))
  angles$ang14.13 <- ifelse(fly13$pos_x <= fly14$pos_x, c((anglearcsin$arcsin14.13*-1)+pi), c(ifelse((fly13$pos_x > fly14$pos_x & fly13$pos_y <= fly14$pos_y), c(anglearcsin$arcsin14.13+(2*pi)), c(anglearcsin$arcsin14.13))))
  angles$ang14.15 <- ifelse(fly15$pos_x <= fly14$pos_x, c((anglearcsin$arcsin14.15*-1)+pi), c(ifelse((fly15$pos_x > fly14$pos_x & fly15$pos_y <= fly14$pos_y), c(anglearcsin$arcsin14.15+(2*pi)), c(anglearcsin$arcsin14.15))))
  angles$ang14.16 <- ifelse(fly16$pos_x <= fly14$pos_x, c((anglearcsin$arcsin14.16*-1)+pi), c(ifelse((fly16$pos_x > fly14$pos_x & fly16$pos_y <= fly14$pos_y), c(anglearcsin$arcsin14.16+(2*pi)), c(anglearcsin$arcsin14.16))))
  angles$ang14.17 <- ifelse(fly17$pos_x <= fly14$pos_x, c((anglearcsin$arcsin14.17*-1)+pi), c(ifelse((fly17$pos_x > fly14$pos_x & fly17$pos_y <= fly14$pos_y), c(anglearcsin$arcsin14.17+(2*pi)), c(anglearcsin$arcsin14.17))))
  angles$ang14.18 <- ifelse(fly18$pos_x <= fly14$pos_x, c((anglearcsin$arcsin14.18*-1)+pi), c(ifelse((fly18$pos_x > fly14$pos_x & fly18$pos_y <= fly14$pos_y), c(anglearcsin$arcsin14.18+(2*pi)), c(anglearcsin$arcsin14.18))))
  angles$ang14.19 <- ifelse(fly19$pos_x <= fly14$pos_x, c((anglearcsin$arcsin14.19*-1)+pi), c(ifelse((fly19$pos_x > fly14$pos_x & fly19$pos_y <= fly14$pos_y), c(anglearcsin$arcsin14.19+(2*pi)), c(anglearcsin$arcsin14.19))))
  angles$ang14.20 <- ifelse(fly20$pos_x <= fly14$pos_x, c((anglearcsin$arcsin14.20*-1)+pi), c(ifelse((fly20$pos_x > fly14$pos_x & fly20$pos_y <= fly14$pos_y), c(anglearcsin$arcsin14.20+(2*pi)), c(anglearcsin$arcsin14.20))))
  angles$ang15.1 <- ifelse(fly1$pos_x <= fly15$pos_x, c((anglearcsin$arcsin15.1*-1)+pi), c(ifelse((fly1$pos_x > fly15$pos_x & fly1$pos_y <= fly15$pos_y), c(anglearcsin$arcsin15.1+(2*pi)), c(anglearcsin$arcsin15.1))))
  angles$ang15.2 <- ifelse(fly2$pos_x <= fly15$pos_x, c((anglearcsin$arcsin15.2*-1)+pi), c(ifelse((fly2$pos_x > fly15$pos_x & fly2$pos_y <= fly15$pos_y), c(anglearcsin$arcsin15.2+(2*pi)), c(anglearcsin$arcsin15.2))))
  angles$ang15.3 <- ifelse(fly3$pos_x <= fly15$pos_x, c((anglearcsin$arcsin15.3*-1)+pi), c(ifelse((fly3$pos_x > fly15$pos_x & fly3$pos_y <= fly15$pos_y), c(anglearcsin$arcsin15.3+(2*pi)), c(anglearcsin$arcsin15.3))))
  angles$ang15.4 <- ifelse(fly4$pos_x <= fly15$pos_x, c((anglearcsin$arcsin15.4*-1)+pi), c(ifelse((fly4$pos_x > fly15$pos_x & fly4$pos_y <= fly15$pos_y), c(anglearcsin$arcsin15.4+(2*pi)), c(anglearcsin$arcsin15.4))))
  angles$ang15.5 <- ifelse(fly5$pos_x <= fly15$pos_x, c((anglearcsin$arcsin15.5*-1)+pi), c(ifelse((fly5$pos_x > fly15$pos_x & fly5$pos_y <= fly15$pos_y), c(anglearcsin$arcsin15.5+(2*pi)), c(anglearcsin$arcsin15.5))))
  angles$ang15.6 <- ifelse(fly6$pos_x <= fly15$pos_x, c((anglearcsin$arcsin15.6*-1)+pi), c(ifelse((fly6$pos_x > fly15$pos_x & fly6$pos_y <= fly15$pos_y), c(anglearcsin$arcsin15.6+(2*pi)), c(anglearcsin$arcsin15.6))))
  angles$ang15.7 <- ifelse(fly7$pos_x <= fly15$pos_x, c((anglearcsin$arcsin15.7*-1)+pi), c(ifelse((fly7$pos_x > fly15$pos_x & fly7$pos_y <= fly15$pos_y), c(anglearcsin$arcsin15.7+(2*pi)), c(anglearcsin$arcsin15.7))))
  angles$ang15.8 <- ifelse(fly8$pos_x <= fly15$pos_x, c((anglearcsin$arcsin15.8*-1)+pi), c(ifelse((fly8$pos_x > fly15$pos_x & fly8$pos_y <= fly15$pos_y), c(anglearcsin$arcsin15.8+(2*pi)), c(anglearcsin$arcsin15.8))))
  angles$ang15.9 <- ifelse(fly9$pos_x <= fly15$pos_x, c((anglearcsin$arcsin15.9*-1)+pi), c(ifelse((fly9$pos_x > fly15$pos_x & fly9$pos_y <= fly15$pos_y), c(anglearcsin$arcsin15.9+(2*pi)), c(anglearcsin$arcsin15.9))))
  angles$ang15.10 <- ifelse(fly10$pos_x <= fly15$pos_x, c((anglearcsin$arcsin15.10*-1)+pi), c(ifelse((fly10$pos_x > fly15$pos_x & fly10$pos_y <= fly15$pos_y), c(anglearcsin$arcsin15.10+(2*pi)), c(anglearcsin$arcsin15.10))))
  angles$ang15.11 <- ifelse(fly11$pos_x <= fly15$pos_x, c((anglearcsin$arcsin15.11*-1)+pi), c(ifelse((fly11$pos_x > fly15$pos_x & fly11$pos_y <= fly15$pos_y), c(anglearcsin$arcsin15.11+(2*pi)), c(anglearcsin$arcsin15.11))))
  angles$ang15.12 <- ifelse(fly12$pos_x <= fly15$pos_x, c((anglearcsin$arcsin15.12*-1)+pi), c(ifelse((fly12$pos_x > fly15$pos_x & fly12$pos_y <= fly15$pos_y), c(anglearcsin$arcsin15.12+(2*pi)), c(anglearcsin$arcsin15.12))))
  angles$ang15.13 <- ifelse(fly13$pos_x <= fly15$pos_x, c((anglearcsin$arcsin15.13*-1)+pi), c(ifelse((fly13$pos_x > fly15$pos_x & fly13$pos_y <= fly15$pos_y), c(anglearcsin$arcsin15.13+(2*pi)), c(anglearcsin$arcsin15.13))))
  angles$ang15.14 <- ifelse(fly14$pos_x <= fly15$pos_x, c((anglearcsin$arcsin15.14*-1)+pi), c(ifelse((fly14$pos_x > fly15$pos_x & fly14$pos_y <= fly15$pos_y), c(anglearcsin$arcsin15.14+(2*pi)), c(anglearcsin$arcsin15.14))))
  angles$ang15.16 <- ifelse(fly16$pos_x <= fly15$pos_x, c((anglearcsin$arcsin15.16*-1)+pi), c(ifelse((fly16$pos_x > fly15$pos_x & fly16$pos_y <= fly15$pos_y), c(anglearcsin$arcsin15.16+(2*pi)), c(anglearcsin$arcsin15.16))))
  angles$ang15.17 <- ifelse(fly17$pos_x <= fly15$pos_x, c((anglearcsin$arcsin15.17*-1)+pi), c(ifelse((fly17$pos_x > fly15$pos_x & fly17$pos_y <= fly15$pos_y), c(anglearcsin$arcsin15.17+(2*pi)), c(anglearcsin$arcsin15.17))))
  angles$ang15.18 <- ifelse(fly18$pos_x <= fly15$pos_x, c((anglearcsin$arcsin15.18*-1)+pi), c(ifelse((fly18$pos_x > fly15$pos_x & fly18$pos_y <= fly15$pos_y), c(anglearcsin$arcsin15.18+(2*pi)), c(anglearcsin$arcsin15.18))))
  angles$ang15.19 <- ifelse(fly19$pos_x <= fly15$pos_x, c((anglearcsin$arcsin15.19*-1)+pi), c(ifelse((fly19$pos_x > fly15$pos_x & fly19$pos_y <= fly15$pos_y), c(anglearcsin$arcsin15.19+(2*pi)), c(anglearcsin$arcsin15.19))))
  angles$ang15.20 <- ifelse(fly20$pos_x <= fly15$pos_x, c((anglearcsin$arcsin15.20*-1)+pi), c(ifelse((fly20$pos_x > fly15$pos_x & fly20$pos_y <= fly15$pos_y), c(anglearcsin$arcsin15.20+(2*pi)), c(anglearcsin$arcsin15.20))))
  angles$ang16.1 <- ifelse(fly1$pos_x <= fly16$pos_x, c((anglearcsin$arcsin16.1*-1)+pi), c(ifelse((fly1$pos_x > fly16$pos_x & fly1$pos_y <= fly16$pos_y), c(anglearcsin$arcsin16.1+(2*pi)), c(anglearcsin$arcsin16.1))))
  angles$ang16.2 <- ifelse(fly2$pos_x <= fly16$pos_x, c((anglearcsin$arcsin16.2*-1)+pi), c(ifelse((fly2$pos_x > fly16$pos_x & fly2$pos_y <= fly16$pos_y), c(anglearcsin$arcsin16.2+(2*pi)), c(anglearcsin$arcsin16.2))))
  angles$ang16.3 <- ifelse(fly3$pos_x <= fly16$pos_x, c((anglearcsin$arcsin16.3*-1)+pi), c(ifelse((fly3$pos_x > fly16$pos_x & fly3$pos_y <= fly16$pos_y), c(anglearcsin$arcsin16.3+(2*pi)), c(anglearcsin$arcsin16.3))))
  angles$ang16.4 <- ifelse(fly4$pos_x <= fly16$pos_x, c((anglearcsin$arcsin16.4*-1)+pi), c(ifelse((fly4$pos_x > fly16$pos_x & fly4$pos_y <= fly16$pos_y), c(anglearcsin$arcsin16.4+(2*pi)), c(anglearcsin$arcsin16.4))))
  angles$ang16.5 <- ifelse(fly5$pos_x <= fly16$pos_x, c((anglearcsin$arcsin16.5*-1)+pi), c(ifelse((fly5$pos_x > fly16$pos_x & fly5$pos_y <= fly16$pos_y), c(anglearcsin$arcsin16.5+(2*pi)), c(anglearcsin$arcsin16.5))))
  angles$ang16.6 <- ifelse(fly6$pos_x <= fly16$pos_x, c((anglearcsin$arcsin16.6*-1)+pi), c(ifelse((fly6$pos_x > fly16$pos_x & fly6$pos_y <= fly16$pos_y), c(anglearcsin$arcsin16.6+(2*pi)), c(anglearcsin$arcsin16.6))))
  angles$ang16.7 <- ifelse(fly7$pos_x <= fly16$pos_x, c((anglearcsin$arcsin16.7*-1)+pi), c(ifelse((fly7$pos_x > fly16$pos_x & fly7$pos_y <= fly16$pos_y), c(anglearcsin$arcsin16.7+(2*pi)), c(anglearcsin$arcsin16.7))))
  angles$ang16.8 <- ifelse(fly8$pos_x <= fly16$pos_x, c((anglearcsin$arcsin16.8*-1)+pi), c(ifelse((fly8$pos_x > fly16$pos_x & fly8$pos_y <= fly16$pos_y), c(anglearcsin$arcsin16.8+(2*pi)), c(anglearcsin$arcsin16.8))))
  angles$ang16.9 <- ifelse(fly9$pos_x <= fly16$pos_x, c((anglearcsin$arcsin16.9*-1)+pi), c(ifelse((fly9$pos_x > fly16$pos_x & fly9$pos_y <= fly16$pos_y), c(anglearcsin$arcsin16.9+(2*pi)), c(anglearcsin$arcsin16.9))))
  angles$ang16.10 <- ifelse(fly10$pos_x <= fly16$pos_x, c((anglearcsin$arcsin16.10*-1)+pi), c(ifelse((fly10$pos_x > fly16$pos_x & fly10$pos_y <= fly16$pos_y), c(anglearcsin$arcsin16.10+(2*pi)), c(anglearcsin$arcsin16.10))))
  angles$ang16.11 <- ifelse(fly11$pos_x <= fly16$pos_x, c((anglearcsin$arcsin16.11*-1)+pi), c(ifelse((fly11$pos_x > fly16$pos_x & fly11$pos_y <= fly16$pos_y), c(anglearcsin$arcsin16.11+(2*pi)), c(anglearcsin$arcsin16.11))))
  angles$ang16.12 <- ifelse(fly12$pos_x <= fly16$pos_x, c((anglearcsin$arcsin16.12*-1)+pi), c(ifelse((fly12$pos_x > fly16$pos_x & fly12$pos_y <= fly16$pos_y), c(anglearcsin$arcsin16.12+(2*pi)), c(anglearcsin$arcsin16.12))))
  angles$ang16.13 <- ifelse(fly13$pos_x <= fly16$pos_x, c((anglearcsin$arcsin16.13*-1)+pi), c(ifelse((fly13$pos_x > fly16$pos_x & fly13$pos_y <= fly16$pos_y), c(anglearcsin$arcsin16.13+(2*pi)), c(anglearcsin$arcsin16.13))))
  angles$ang16.14 <- ifelse(fly14$pos_x <= fly16$pos_x, c((anglearcsin$arcsin16.14*-1)+pi), c(ifelse((fly14$pos_x > fly16$pos_x & fly14$pos_y <= fly16$pos_y), c(anglearcsin$arcsin16.14+(2*pi)), c(anglearcsin$arcsin16.14))))
  angles$ang16.15 <- ifelse(fly15$pos_x <= fly16$pos_x, c((anglearcsin$arcsin16.15*-1)+pi), c(ifelse((fly15$pos_x > fly16$pos_x & fly15$pos_y <= fly16$pos_y), c(anglearcsin$arcsin16.15+(2*pi)), c(anglearcsin$arcsin16.15))))
  angles$ang16.17 <- ifelse(fly17$pos_x <= fly16$pos_x, c((anglearcsin$arcsin16.17*-1)+pi), c(ifelse((fly17$pos_x > fly16$pos_x & fly17$pos_y <= fly16$pos_y), c(anglearcsin$arcsin16.17+(2*pi)), c(anglearcsin$arcsin16.17))))
  angles$ang16.18 <- ifelse(fly18$pos_x <= fly16$pos_x, c((anglearcsin$arcsin16.18*-1)+pi), c(ifelse((fly18$pos_x > fly16$pos_x & fly18$pos_y <= fly16$pos_y), c(anglearcsin$arcsin16.18+(2*pi)), c(anglearcsin$arcsin16.18))))
  angles$ang16.19 <- ifelse(fly19$pos_x <= fly16$pos_x, c((anglearcsin$arcsin16.19*-1)+pi), c(ifelse((fly19$pos_x > fly16$pos_x & fly19$pos_y <= fly16$pos_y), c(anglearcsin$arcsin16.19+(2*pi)), c(anglearcsin$arcsin16.19))))
  angles$ang16.20 <- ifelse(fly20$pos_x <= fly16$pos_x, c((anglearcsin$arcsin16.20*-1)+pi), c(ifelse((fly20$pos_x > fly16$pos_x & fly20$pos_y <= fly16$pos_y), c(anglearcsin$arcsin16.20+(2*pi)), c(anglearcsin$arcsin16.20))))
  angles$ang17.1 <- ifelse(fly1$pos_x <= fly17$pos_x, c((anglearcsin$arcsin17.1*-1)+pi), c(ifelse((fly1$pos_x > fly17$pos_x & fly1$pos_y <= fly17$pos_y), c(anglearcsin$arcsin17.1+(2*pi)), c(anglearcsin$arcsin17.1))))
  angles$ang17.2 <- ifelse(fly2$pos_x <= fly17$pos_x, c((anglearcsin$arcsin17.2*-1)+pi), c(ifelse((fly2$pos_x > fly17$pos_x & fly2$pos_y <= fly17$pos_y), c(anglearcsin$arcsin17.2+(2*pi)), c(anglearcsin$arcsin17.2))))
  angles$ang17.3 <- ifelse(fly3$pos_x <= fly17$pos_x, c((anglearcsin$arcsin17.3*-1)+pi), c(ifelse((fly3$pos_x > fly17$pos_x & fly3$pos_y <= fly17$pos_y), c(anglearcsin$arcsin17.3+(2*pi)), c(anglearcsin$arcsin17.3))))
  angles$ang17.4 <- ifelse(fly4$pos_x <= fly17$pos_x, c((anglearcsin$arcsin17.4*-1)+pi), c(ifelse((fly4$pos_x > fly17$pos_x & fly4$pos_y <= fly17$pos_y), c(anglearcsin$arcsin17.4+(2*pi)), c(anglearcsin$arcsin17.4))))
  angles$ang17.5 <- ifelse(fly5$pos_x <= fly17$pos_x, c((anglearcsin$arcsin17.5*-1)+pi), c(ifelse((fly5$pos_x > fly17$pos_x & fly5$pos_y <= fly17$pos_y), c(anglearcsin$arcsin17.5+(2*pi)), c(anglearcsin$arcsin17.5))))
  angles$ang17.6 <- ifelse(fly6$pos_x <= fly17$pos_x, c((anglearcsin$arcsin17.6*-1)+pi), c(ifelse((fly6$pos_x > fly17$pos_x & fly6$pos_y <= fly17$pos_y), c(anglearcsin$arcsin17.6+(2*pi)), c(anglearcsin$arcsin17.6))))
  angles$ang17.7 <- ifelse(fly7$pos_x <= fly17$pos_x, c((anglearcsin$arcsin17.7*-1)+pi), c(ifelse((fly7$pos_x > fly17$pos_x & fly7$pos_y <= fly17$pos_y), c(anglearcsin$arcsin17.7+(2*pi)), c(anglearcsin$arcsin17.7))))
  angles$ang17.8 <- ifelse(fly8$pos_x <= fly17$pos_x, c((anglearcsin$arcsin17.8*-1)+pi), c(ifelse((fly8$pos_x > fly17$pos_x & fly8$pos_y <= fly17$pos_y), c(anglearcsin$arcsin17.8+(2*pi)), c(anglearcsin$arcsin17.8))))
  angles$ang17.9 <- ifelse(fly9$pos_x <= fly17$pos_x, c((anglearcsin$arcsin17.9*-1)+pi), c(ifelse((fly9$pos_x > fly17$pos_x & fly9$pos_y <= fly17$pos_y), c(anglearcsin$arcsin17.9+(2*pi)), c(anglearcsin$arcsin17.9))))
  angles$ang17.10 <- ifelse(fly10$pos_x <= fly17$pos_x, c((anglearcsin$arcsin17.10*-1)+pi), c(ifelse((fly10$pos_x > fly17$pos_x & fly10$pos_y <= fly17$pos_y), c(anglearcsin$arcsin17.10+(2*pi)), c(anglearcsin$arcsin17.10))))
  angles$ang17.11 <- ifelse(fly11$pos_x <= fly17$pos_x, c((anglearcsin$arcsin17.11*-1)+pi), c(ifelse((fly11$pos_x > fly17$pos_x & fly11$pos_y <= fly17$pos_y), c(anglearcsin$arcsin17.11+(2*pi)), c(anglearcsin$arcsin17.11))))
  angles$ang17.12 <- ifelse(fly12$pos_x <= fly17$pos_x, c((anglearcsin$arcsin17.12*-1)+pi), c(ifelse((fly12$pos_x > fly17$pos_x & fly12$pos_y <= fly17$pos_y), c(anglearcsin$arcsin17.12+(2*pi)), c(anglearcsin$arcsin17.12))))
  angles$ang17.13 <- ifelse(fly13$pos_x <= fly17$pos_x, c((anglearcsin$arcsin17.13*-1)+pi), c(ifelse((fly13$pos_x > fly17$pos_x & fly13$pos_y <= fly17$pos_y), c(anglearcsin$arcsin17.13+(2*pi)), c(anglearcsin$arcsin17.13))))
  angles$ang17.14 <- ifelse(fly14$pos_x <= fly17$pos_x, c((anglearcsin$arcsin17.14*-1)+pi), c(ifelse((fly14$pos_x > fly17$pos_x & fly14$pos_y <= fly17$pos_y), c(anglearcsin$arcsin17.14+(2*pi)), c(anglearcsin$arcsin17.14))))
  angles$ang17.15 <- ifelse(fly15$pos_x <= fly17$pos_x, c((anglearcsin$arcsin17.15*-1)+pi), c(ifelse((fly15$pos_x > fly17$pos_x & fly15$pos_y <= fly17$pos_y), c(anglearcsin$arcsin17.15+(2*pi)), c(anglearcsin$arcsin17.15))))
  angles$ang17.16 <- ifelse(fly16$pos_x <= fly17$pos_x, c((anglearcsin$arcsin17.16*-1)+pi), c(ifelse((fly16$pos_x > fly17$pos_x & fly16$pos_y <= fly17$pos_y), c(anglearcsin$arcsin17.16+(2*pi)), c(anglearcsin$arcsin17.16))))
  angles$ang17.18 <- ifelse(fly18$pos_x <= fly17$pos_x, c((anglearcsin$arcsin17.18*-1)+pi), c(ifelse((fly18$pos_x > fly17$pos_x & fly18$pos_y <= fly17$pos_y), c(anglearcsin$arcsin17.18+(2*pi)), c(anglearcsin$arcsin17.18))))
  angles$ang17.19 <- ifelse(fly19$pos_x <= fly17$pos_x, c((anglearcsin$arcsin17.19*-1)+pi), c(ifelse((fly19$pos_x > fly17$pos_x & fly19$pos_y <= fly17$pos_y), c(anglearcsin$arcsin17.19+(2*pi)), c(anglearcsin$arcsin17.19))))
  angles$ang17.20 <- ifelse(fly20$pos_x <= fly17$pos_x, c((anglearcsin$arcsin17.20*-1)+pi), c(ifelse((fly20$pos_x > fly17$pos_x & fly20$pos_y <= fly17$pos_y), c(anglearcsin$arcsin17.20+(2*pi)), c(anglearcsin$arcsin17.20))))
  angles$ang18.1 <- ifelse(fly1$pos_x <= fly18$pos_x, c((anglearcsin$arcsin18.1*-1)+pi), c(ifelse((fly1$pos_x > fly18$pos_x & fly1$pos_y <= fly18$pos_y), c(anglearcsin$arcsin18.1+(2*pi)), c(anglearcsin$arcsin18.1))))
  angles$ang18.2 <- ifelse(fly2$pos_x <= fly18$pos_x, c((anglearcsin$arcsin18.2*-1)+pi), c(ifelse((fly2$pos_x > fly18$pos_x & fly2$pos_y <= fly18$pos_y), c(anglearcsin$arcsin18.2+(2*pi)), c(anglearcsin$arcsin18.2))))
  angles$ang18.3 <- ifelse(fly3$pos_x <= fly18$pos_x, c((anglearcsin$arcsin18.3*-1)+pi), c(ifelse((fly3$pos_x > fly18$pos_x & fly3$pos_y <= fly18$pos_y), c(anglearcsin$arcsin18.3+(2*pi)), c(anglearcsin$arcsin18.3))))
  angles$ang18.4 <- ifelse(fly4$pos_x <= fly18$pos_x, c((anglearcsin$arcsin18.4*-1)+pi), c(ifelse((fly4$pos_x > fly18$pos_x & fly4$pos_y <= fly18$pos_y), c(anglearcsin$arcsin18.4+(2*pi)), c(anglearcsin$arcsin18.4))))
  angles$ang18.5 <- ifelse(fly5$pos_x <= fly18$pos_x, c((anglearcsin$arcsin18.5*-1)+pi), c(ifelse((fly5$pos_x > fly18$pos_x & fly5$pos_y <= fly18$pos_y), c(anglearcsin$arcsin18.5+(2*pi)), c(anglearcsin$arcsin18.5))))
  angles$ang18.6 <- ifelse(fly6$pos_x <= fly18$pos_x, c((anglearcsin$arcsin18.6*-1)+pi), c(ifelse((fly6$pos_x > fly18$pos_x & fly6$pos_y <= fly18$pos_y), c(anglearcsin$arcsin18.6+(2*pi)), c(anglearcsin$arcsin18.6))))
  angles$ang18.7 <- ifelse(fly7$pos_x <= fly18$pos_x, c((anglearcsin$arcsin18.7*-1)+pi), c(ifelse((fly7$pos_x > fly18$pos_x & fly7$pos_y <= fly18$pos_y), c(anglearcsin$arcsin18.7+(2*pi)), c(anglearcsin$arcsin18.7))))
  angles$ang18.8 <- ifelse(fly8$pos_x <= fly18$pos_x, c((anglearcsin$arcsin18.8*-1)+pi), c(ifelse((fly8$pos_x > fly18$pos_x & fly8$pos_y <= fly18$pos_y), c(anglearcsin$arcsin18.8+(2*pi)), c(anglearcsin$arcsin18.8))))
  angles$ang18.9 <- ifelse(fly9$pos_x <= fly18$pos_x, c((anglearcsin$arcsin18.9*-1)+pi), c(ifelse((fly9$pos_x > fly18$pos_x & fly9$pos_y <= fly18$pos_y), c(anglearcsin$arcsin18.9+(2*pi)), c(anglearcsin$arcsin18.9))))
  angles$ang18.10 <- ifelse(fly10$pos_x <= fly18$pos_x, c((anglearcsin$arcsin18.10*-1)+pi), c(ifelse((fly10$pos_x > fly18$pos_x & fly10$pos_y <= fly18$pos_y), c(anglearcsin$arcsin18.10+(2*pi)), c(anglearcsin$arcsin18.10))))
  angles$ang18.11 <- ifelse(fly11$pos_x <= fly18$pos_x, c((anglearcsin$arcsin18.11*-1)+pi), c(ifelse((fly11$pos_x > fly18$pos_x & fly11$pos_y <= fly18$pos_y), c(anglearcsin$arcsin18.11+(2*pi)), c(anglearcsin$arcsin18.11))))
  angles$ang18.12 <- ifelse(fly12$pos_x <= fly18$pos_x, c((anglearcsin$arcsin18.12*-1)+pi), c(ifelse((fly12$pos_x > fly18$pos_x & fly12$pos_y <= fly18$pos_y), c(anglearcsin$arcsin18.12+(2*pi)), c(anglearcsin$arcsin18.12))))
  angles$ang18.13 <- ifelse(fly13$pos_x <= fly18$pos_x, c((anglearcsin$arcsin18.13*-1)+pi), c(ifelse((fly13$pos_x > fly18$pos_x & fly13$pos_y <= fly18$pos_y), c(anglearcsin$arcsin18.13+(2*pi)), c(anglearcsin$arcsin18.13))))
  angles$ang18.14 <- ifelse(fly14$pos_x <= fly18$pos_x, c((anglearcsin$arcsin18.14*-1)+pi), c(ifelse((fly14$pos_x > fly18$pos_x & fly14$pos_y <= fly18$pos_y), c(anglearcsin$arcsin18.14+(2*pi)), c(anglearcsin$arcsin18.14))))
  angles$ang18.15 <- ifelse(fly15$pos_x <= fly18$pos_x, c((anglearcsin$arcsin18.15*-1)+pi), c(ifelse((fly15$pos_x > fly18$pos_x & fly15$pos_y <= fly18$pos_y), c(anglearcsin$arcsin18.15+(2*pi)), c(anglearcsin$arcsin18.15))))
  angles$ang18.16 <- ifelse(fly16$pos_x <= fly18$pos_x, c((anglearcsin$arcsin18.16*-1)+pi), c(ifelse((fly16$pos_x > fly18$pos_x & fly16$pos_y <= fly18$pos_y), c(anglearcsin$arcsin18.16+(2*pi)), c(anglearcsin$arcsin18.16))))
  angles$ang18.17 <- ifelse(fly17$pos_x <= fly18$pos_x, c((anglearcsin$arcsin18.17*-1)+pi), c(ifelse((fly17$pos_x > fly18$pos_x & fly17$pos_y <= fly18$pos_y), c(anglearcsin$arcsin18.17+(2*pi)), c(anglearcsin$arcsin18.17))))
  angles$ang18.19 <- ifelse(fly19$pos_x <= fly18$pos_x, c((anglearcsin$arcsin18.19*-1)+pi), c(ifelse((fly19$pos_x > fly18$pos_x & fly19$pos_y <= fly18$pos_y), c(anglearcsin$arcsin18.19+(2*pi)), c(anglearcsin$arcsin18.19))))
  angles$ang18.20 <- ifelse(fly20$pos_x <= fly18$pos_x, c((anglearcsin$arcsin18.20*-1)+pi), c(ifelse((fly20$pos_x > fly18$pos_x & fly20$pos_y <= fly18$pos_y), c(anglearcsin$arcsin18.20+(2*pi)), c(anglearcsin$arcsin18.20))))
  angles$ang19.1 <- ifelse(fly1$pos_x <= fly19$pos_x, c((anglearcsin$arcsin19.1*-1)+pi), c(ifelse((fly1$pos_x > fly19$pos_x & fly1$pos_y <= fly19$pos_y), c(anglearcsin$arcsin19.1+(2*pi)), c(anglearcsin$arcsin19.1))))
  angles$ang19.2 <- ifelse(fly2$pos_x <= fly19$pos_x, c((anglearcsin$arcsin19.2*-1)+pi), c(ifelse((fly2$pos_x > fly19$pos_x & fly2$pos_y <= fly19$pos_y), c(anglearcsin$arcsin19.2+(2*pi)), c(anglearcsin$arcsin19.2))))
  angles$ang19.3 <- ifelse(fly3$pos_x <= fly19$pos_x, c((anglearcsin$arcsin19.3*-1)+pi), c(ifelse((fly3$pos_x > fly19$pos_x & fly3$pos_y <= fly19$pos_y), c(anglearcsin$arcsin19.3+(2*pi)), c(anglearcsin$arcsin19.3))))
  angles$ang19.4 <- ifelse(fly4$pos_x <= fly19$pos_x, c((anglearcsin$arcsin19.4*-1)+pi), c(ifelse((fly4$pos_x > fly19$pos_x & fly4$pos_y <= fly19$pos_y), c(anglearcsin$arcsin19.4+(2*pi)), c(anglearcsin$arcsin19.4))))
  angles$ang19.5 <- ifelse(fly5$pos_x <= fly19$pos_x, c((anglearcsin$arcsin19.5*-1)+pi), c(ifelse((fly5$pos_x > fly19$pos_x & fly5$pos_y <= fly19$pos_y), c(anglearcsin$arcsin19.5+(2*pi)), c(anglearcsin$arcsin19.5))))
  angles$ang19.6 <- ifelse(fly6$pos_x <= fly19$pos_x, c((anglearcsin$arcsin19.6*-1)+pi), c(ifelse((fly6$pos_x > fly19$pos_x & fly6$pos_y <= fly19$pos_y), c(anglearcsin$arcsin19.6+(2*pi)), c(anglearcsin$arcsin19.6))))
  angles$ang19.7 <- ifelse(fly7$pos_x <= fly19$pos_x, c((anglearcsin$arcsin19.7*-1)+pi), c(ifelse((fly7$pos_x > fly19$pos_x & fly7$pos_y <= fly19$pos_y), c(anglearcsin$arcsin19.7+(2*pi)), c(anglearcsin$arcsin19.7))))
  angles$ang19.8 <- ifelse(fly8$pos_x <= fly19$pos_x, c((anglearcsin$arcsin19.8*-1)+pi), c(ifelse((fly8$pos_x > fly19$pos_x & fly8$pos_y <= fly19$pos_y), c(anglearcsin$arcsin19.8+(2*pi)), c(anglearcsin$arcsin19.8))))
  angles$ang19.9 <- ifelse(fly9$pos_x <= fly19$pos_x, c((anglearcsin$arcsin19.9*-1)+pi), c(ifelse((fly9$pos_x > fly19$pos_x & fly9$pos_y <= fly19$pos_y), c(anglearcsin$arcsin19.9+(2*pi)), c(anglearcsin$arcsin19.9))))
  angles$ang19.10 <- ifelse(fly10$pos_x <= fly19$pos_x, c((anglearcsin$arcsin19.10*-1)+pi), c(ifelse((fly10$pos_x > fly19$pos_x & fly10$pos_y <= fly19$pos_y), c(anglearcsin$arcsin19.10+(2*pi)), c(anglearcsin$arcsin19.10))))
  angles$ang19.11 <- ifelse(fly11$pos_x <= fly19$pos_x, c((anglearcsin$arcsin19.11*-1)+pi), c(ifelse((fly11$pos_x > fly19$pos_x & fly11$pos_y <= fly19$pos_y), c(anglearcsin$arcsin19.11+(2*pi)), c(anglearcsin$arcsin19.11))))
  angles$ang19.12 <- ifelse(fly12$pos_x <= fly19$pos_x, c((anglearcsin$arcsin19.12*-1)+pi), c(ifelse((fly12$pos_x > fly19$pos_x & fly12$pos_y <= fly19$pos_y), c(anglearcsin$arcsin19.12+(2*pi)), c(anglearcsin$arcsin19.12))))
  angles$ang19.13 <- ifelse(fly13$pos_x <= fly19$pos_x, c((anglearcsin$arcsin19.13*-1)+pi), c(ifelse((fly13$pos_x > fly19$pos_x & fly13$pos_y <= fly19$pos_y), c(anglearcsin$arcsin19.13+(2*pi)), c(anglearcsin$arcsin19.13))))
  angles$ang19.14 <- ifelse(fly14$pos_x <= fly19$pos_x, c((anglearcsin$arcsin19.14*-1)+pi), c(ifelse((fly14$pos_x > fly19$pos_x & fly14$pos_y <= fly19$pos_y), c(anglearcsin$arcsin19.14+(2*pi)), c(anglearcsin$arcsin19.14))))
  angles$ang19.15 <- ifelse(fly15$pos_x <= fly19$pos_x, c((anglearcsin$arcsin19.15*-1)+pi), c(ifelse((fly15$pos_x > fly19$pos_x & fly15$pos_y <= fly19$pos_y), c(anglearcsin$arcsin19.15+(2*pi)), c(anglearcsin$arcsin19.15))))
  angles$ang19.16 <- ifelse(fly16$pos_x <= fly19$pos_x, c((anglearcsin$arcsin19.16*-1)+pi), c(ifelse((fly16$pos_x > fly19$pos_x & fly16$pos_y <= fly19$pos_y), c(anglearcsin$arcsin19.16+(2*pi)), c(anglearcsin$arcsin19.16))))
  angles$ang19.17 <- ifelse(fly17$pos_x <= fly19$pos_x, c((anglearcsin$arcsin19.17*-1)+pi), c(ifelse((fly17$pos_x > fly19$pos_x & fly17$pos_y <= fly19$pos_y), c(anglearcsin$arcsin19.17+(2*pi)), c(anglearcsin$arcsin19.17))))
  angles$ang19.18 <- ifelse(fly18$pos_x <= fly19$pos_x, c((anglearcsin$arcsin19.18*-1)+pi), c(ifelse((fly18$pos_x > fly19$pos_x & fly18$pos_y <= fly19$pos_y), c(anglearcsin$arcsin19.18+(2*pi)), c(anglearcsin$arcsin19.18))))
  angles$ang19.20 <- ifelse(fly20$pos_x <= fly19$pos_x, c((anglearcsin$arcsin19.20*-1)+pi), c(ifelse((fly20$pos_x > fly19$pos_x & fly20$pos_y <= fly19$pos_y), c(anglearcsin$arcsin19.20+(2*pi)), c(anglearcsin$arcsin19.20))))
  angles$ang20.1 <- ifelse(fly1$pos_x <= fly20$pos_x, c((anglearcsin$arcsin20.1*-1)+pi), c(ifelse((fly1$pos_x > fly20$pos_x & fly1$pos_y <= fly20$pos_y), c(anglearcsin$arcsin20.1+(2*pi)), c(anglearcsin$arcsin20.1))))
  angles$ang20.2 <- ifelse(fly2$pos_x <= fly20$pos_x, c((anglearcsin$arcsin20.2*-1)+pi), c(ifelse((fly2$pos_x > fly20$pos_x & fly2$pos_y <= fly20$pos_y), c(anglearcsin$arcsin20.2+(2*pi)), c(anglearcsin$arcsin20.2))))
  angles$ang20.3 <- ifelse(fly3$pos_x <= fly20$pos_x, c((anglearcsin$arcsin20.3*-1)+pi), c(ifelse((fly3$pos_x > fly20$pos_x & fly3$pos_y <= fly20$pos_y), c(anglearcsin$arcsin20.3+(2*pi)), c(anglearcsin$arcsin20.3))))
  angles$ang20.4 <- ifelse(fly4$pos_x <= fly20$pos_x, c((anglearcsin$arcsin20.4*-1)+pi), c(ifelse((fly4$pos_x > fly20$pos_x & fly4$pos_y <= fly20$pos_y), c(anglearcsin$arcsin20.4+(2*pi)), c(anglearcsin$arcsin20.4))))
  angles$ang20.5 <- ifelse(fly5$pos_x <= fly20$pos_x, c((anglearcsin$arcsin20.5*-1)+pi), c(ifelse((fly5$pos_x > fly20$pos_x & fly5$pos_y <= fly20$pos_y), c(anglearcsin$arcsin20.5+(2*pi)), c(anglearcsin$arcsin20.5))))
  angles$ang20.6 <- ifelse(fly6$pos_x <= fly20$pos_x, c((anglearcsin$arcsin20.6*-1)+pi), c(ifelse((fly6$pos_x > fly20$pos_x & fly6$pos_y <= fly20$pos_y), c(anglearcsin$arcsin20.6+(2*pi)), c(anglearcsin$arcsin20.6))))
  angles$ang20.7 <- ifelse(fly7$pos_x <= fly20$pos_x, c((anglearcsin$arcsin20.7*-1)+pi), c(ifelse((fly7$pos_x > fly20$pos_x & fly7$pos_y <= fly20$pos_y), c(anglearcsin$arcsin20.7+(2*pi)), c(anglearcsin$arcsin20.7))))
  angles$ang20.8 <- ifelse(fly8$pos_x <= fly20$pos_x, c((anglearcsin$arcsin20.8*-1)+pi), c(ifelse((fly8$pos_x > fly20$pos_x & fly8$pos_y <= fly20$pos_y), c(anglearcsin$arcsin20.8+(2*pi)), c(anglearcsin$arcsin20.8))))
  angles$ang20.9 <- ifelse(fly9$pos_x <= fly20$pos_x, c((anglearcsin$arcsin20.9*-1)+pi), c(ifelse((fly9$pos_x > fly20$pos_x & fly9$pos_y <= fly20$pos_y), c(anglearcsin$arcsin20.9+(2*pi)), c(anglearcsin$arcsin20.9))))
  angles$ang20.10 <- ifelse(fly10$pos_x <= fly20$pos_x, c((anglearcsin$arcsin20.10*-1)+pi), c(ifelse((fly10$pos_x > fly20$pos_x & fly10$pos_y <= fly20$pos_y), c(anglearcsin$arcsin20.10+(2*pi)), c(anglearcsin$arcsin20.10))))
  angles$ang20.11 <- ifelse(fly11$pos_x <= fly20$pos_x, c((anglearcsin$arcsin20.11*-1)+pi), c(ifelse((fly11$pos_x > fly20$pos_x & fly11$pos_y <= fly20$pos_y), c(anglearcsin$arcsin20.11+(2*pi)), c(anglearcsin$arcsin20.11))))
  angles$ang20.12 <- ifelse(fly12$pos_x <= fly20$pos_x, c((anglearcsin$arcsin20.12*-1)+pi), c(ifelse((fly12$pos_x > fly20$pos_x & fly12$pos_y <= fly20$pos_y), c(anglearcsin$arcsin20.12+(2*pi)), c(anglearcsin$arcsin20.12))))
  angles$ang20.13 <- ifelse(fly13$pos_x <= fly20$pos_x, c((anglearcsin$arcsin20.13*-1)+pi), c(ifelse((fly13$pos_x > fly20$pos_x & fly13$pos_y <= fly20$pos_y), c(anglearcsin$arcsin20.13+(2*pi)), c(anglearcsin$arcsin20.13))))
  angles$ang20.14 <- ifelse(fly14$pos_x <= fly20$pos_x, c((anglearcsin$arcsin20.14*-1)+pi), c(ifelse((fly14$pos_x > fly20$pos_x & fly14$pos_y <= fly20$pos_y), c(anglearcsin$arcsin20.14+(2*pi)), c(anglearcsin$arcsin20.14))))
  angles$ang20.15 <- ifelse(fly15$pos_x <= fly20$pos_x, c((anglearcsin$arcsin20.15*-1)+pi), c(ifelse((fly15$pos_x > fly20$pos_x & fly15$pos_y <= fly20$pos_y), c(anglearcsin$arcsin20.15+(2*pi)), c(anglearcsin$arcsin20.15))))
  angles$ang20.16 <- ifelse(fly16$pos_x <= fly20$pos_x, c((anglearcsin$arcsin20.16*-1)+pi), c(ifelse((fly16$pos_x > fly20$pos_x & fly16$pos_y <= fly20$pos_y), c(anglearcsin$arcsin20.16+(2*pi)), c(anglearcsin$arcsin20.16))))
  angles$ang20.17 <- ifelse(fly17$pos_x <= fly20$pos_x, c((anglearcsin$arcsin20.17*-1)+pi), c(ifelse((fly17$pos_x > fly20$pos_x & fly17$pos_y <= fly20$pos_y), c(anglearcsin$arcsin20.17+(2*pi)), c(anglearcsin$arcsin20.17))))
  angles$ang20.18 <- ifelse(fly18$pos_x <= fly20$pos_x, c((anglearcsin$arcsin20.18*-1)+pi), c(ifelse((fly18$pos_x > fly20$pos_x & fly18$pos_y <= fly20$pos_y), c(anglearcsin$arcsin20.18+(2*pi)), c(anglearcsin$arcsin20.18))))
  angles$ang20.19 <- ifelse(fly19$pos_x <= fly20$pos_x, c((anglearcsin$arcsin20.19*-1)+pi), c(ifelse((fly19$pos_x > fly20$pos_x & fly19$pos_y <= fly20$pos_y), c(anglearcsin$arcsin20.19+(2*pi)), c(anglearcsin$arcsin20.19))))
  
  # CREATE A DATAFRAME THAT CALCULATES THE DIFFERENCE BETWEEN A FOCAL FLY'S ORIENTATION AND IT'S ANGLE TO EVERY OTHER FLY # --------
  angledifference = data.frame(abs(fly1$oriradians - angles$ang1.2))
  colnames(angledifference) <- c("diff1.2")
  angledifference$diff1.2 <- ifelse(angledifference$diff1.2 > pi, c(abs(angledifference$diff1.2 - (2*pi))), c(angledifference$diff1.2))
  angledifference$diff1.3 <- abs(fly1$oriradians - angles$ang1.3)
  angledifference$diff1.3 <- ifelse(angledifference$diff1.3 > pi, c(abs(angledifference$diff1.3 - (2*pi))), c(angledifference$diff1.3))
  angledifference$diff1.4 <- abs(fly1$oriradians - angles$ang1.4)
  angledifference$diff1.4 <- ifelse(angledifference$diff1.4 > pi, c(abs(angledifference$diff1.4 - (2*pi))), c(angledifference$diff1.4))
  angledifference$diff1.5 <- abs(fly1$oriradians - angles$ang1.5)
  angledifference$diff1.5 <- ifelse(angledifference$diff1.5 > pi, c(abs(angledifference$diff1.5 - (2*pi))), c(angledifference$diff1.5))
  angledifference$diff1.6 <- abs(fly1$oriradians - angles$ang1.6)
  angledifference$diff1.6 <- ifelse(angledifference$diff1.6 > pi, c(abs(angledifference$diff1.6 - (2*pi))), c(angledifference$diff1.6))
  angledifference$diff1.7 <- abs(fly1$oriradians - angles$ang1.7)
  angledifference$diff1.7 <- ifelse(angledifference$diff1.7 > pi, c(abs(angledifference$diff1.7 - (2*pi))), c(angledifference$diff1.7))
  angledifference$diff1.8 <- abs(fly1$oriradians - angles$ang1.8)
  angledifference$diff1.8 <- ifelse(angledifference$diff1.8 > pi, c(abs(angledifference$diff1.8 - (2*pi))), c(angledifference$diff1.8))
  angledifference$diff1.9 <- abs(fly1$oriradians - angles$ang1.9)
  angledifference$diff1.9 <- ifelse(angledifference$diff1.9 > pi, c(abs(angledifference$diff1.9 - (2*pi))), c(angledifference$diff1.9))
  angledifference$diff1.10 <- abs(fly1$oriradians - angles$ang1.10)
  angledifference$diff1.10 <- ifelse(angledifference$diff1.10 > pi, c(abs(angledifference$diff1.10 - (2*pi))), c(angledifference$diff1.10))
  angledifference$diff1.11 <- abs(fly1$oriradians - angles$ang1.11)
  angledifference$diff1.11 <- ifelse(angledifference$diff1.11 > pi, c(abs(angledifference$diff1.11 - (2*pi))), c(angledifference$diff1.11))
  angledifference$diff1.12 <- abs(fly1$oriradians - angles$ang1.12)
  angledifference$diff1.12 <- ifelse(angledifference$diff1.12 > pi, c(abs(angledifference$diff1.12 - (2*pi))), c(angledifference$diff1.12))
  angledifference$diff1.13 <- abs(fly1$oriradians - angles$ang1.13)
  angledifference$diff1.13 <- ifelse(angledifference$diff1.13 > pi, c(abs(angledifference$diff1.13 - (2*pi))), c(angledifference$diff1.13))
  angledifference$diff1.14 <- abs(fly1$oriradians - angles$ang1.14)
  angledifference$diff1.14 <- ifelse(angledifference$diff1.14 > pi, c(abs(angledifference$diff1.14 - (2*pi))), c(angledifference$diff1.14))
  angledifference$diff1.15 <- abs(fly1$oriradians - angles$ang1.15)
  angledifference$diff1.15 <- ifelse(angledifference$diff1.15 > pi, c(abs(angledifference$diff1.15 - (2*pi))), c(angledifference$diff1.15))
  angledifference$diff1.16 <- abs(fly1$oriradians - angles$ang1.16)
  angledifference$diff1.16 <- ifelse(angledifference$diff1.16 > pi, c(abs(angledifference$diff1.16 - (2*pi))), c(angledifference$diff1.16))
  angledifference$diff1.17 <- abs(fly1$oriradians - angles$ang1.17)
  angledifference$diff1.17 <- ifelse(angledifference$diff1.17 > pi, c(abs(angledifference$diff1.17 - (2*pi))), c(angledifference$diff1.17))
  angledifference$diff1.18 <- abs(fly1$oriradians - angles$ang1.18)
  angledifference$diff1.18 <- ifelse(angledifference$diff1.18 > pi, c(abs(angledifference$diff1.18 - (2*pi))), c(angledifference$diff1.18))
  angledifference$diff1.19 <- abs(fly1$oriradians - angles$ang1.19)
  angledifference$diff1.19 <- ifelse(angledifference$diff1.19 > pi, c(abs(angledifference$diff1.19 - (2*pi))), c(angledifference$diff1.19))
  angledifference$diff1.20 <- abs(fly1$oriradians - angles$ang1.20)
  angledifference$diff1.20 <- ifelse(angledifference$diff1.20 > pi, c(abs(angledifference$diff1.20 - (2*pi))), c(angledifference$diff1.20))
  angledifference$diff2.1 <- abs(fly2$oriradians - angles$ang2.1)
  angledifference$diff2.1 <- ifelse(angledifference$diff2.1 > pi, c(abs(angledifference$diff2.1 - (2*pi))), c(angledifference$diff2.1))
  angledifference$diff2.3 <- abs(fly2$oriradians - angles$ang2.3)
  angledifference$diff2.3 <- ifelse(angledifference$diff2.3 > pi, c(abs(angledifference$diff2.3 - (2*pi))), c(angledifference$diff2.3))
  angledifference$diff2.4 <- abs(fly2$oriradians - angles$ang2.4)
  angledifference$diff2.4 <- ifelse(angledifference$diff2.4 > pi, c(abs(angledifference$diff2.4 - (2*pi))), c(angledifference$diff2.4))
  angledifference$diff2.5 <- abs(fly2$oriradians - angles$ang2.5)
  angledifference$diff2.5 <- ifelse(angledifference$diff2.5 > pi, c(abs(angledifference$diff2.5 - (2*pi))), c(angledifference$diff2.5))
  angledifference$diff2.6 <- abs(fly2$oriradians - angles$ang2.6)
  angledifference$diff2.6 <- ifelse(angledifference$diff2.6 > pi, c(abs(angledifference$diff2.6 - (2*pi))), c(angledifference$diff2.6))
  angledifference$diff2.7 <- abs(fly2$oriradians - angles$ang2.7)
  angledifference$diff2.7 <- ifelse(angledifference$diff2.7 > pi, c(abs(angledifference$diff2.7 - (2*pi))), c(angledifference$diff2.7))
  angledifference$diff2.8 <- abs(fly2$oriradians - angles$ang2.8)
  angledifference$diff2.8 <- ifelse(angledifference$diff2.8 > pi, c(abs(angledifference$diff2.8 - (2*pi))), c(angledifference$diff2.8))
  angledifference$diff2.9 <- abs(fly2$oriradians - angles$ang2.9)
  angledifference$diff2.9 <- ifelse(angledifference$diff2.9 > pi, c(abs(angledifference$diff2.9 - (2*pi))), c(angledifference$diff2.9))
  angledifference$diff2.10 <- abs(fly2$oriradians - angles$ang2.10)
  angledifference$diff2.10 <- ifelse(angledifference$diff2.10 > pi, c(abs(angledifference$diff2.10 - (2*pi))), c(angledifference$diff2.10))
  angledifference$diff2.11 <- abs(fly2$oriradians - angles$ang2.11)
  angledifference$diff2.11 <- ifelse(angledifference$diff2.11 > pi, c(abs(angledifference$diff2.11 - (2*pi))), c(angledifference$diff2.11))
  angledifference$diff2.12 <- abs(fly2$oriradians - angles$ang2.12)
  angledifference$diff2.12 <- ifelse(angledifference$diff2.12 > pi, c(abs(angledifference$diff2.12 - (2*pi))), c(angledifference$diff2.12))
  angledifference$diff2.13 <- abs(fly2$oriradians - angles$ang2.13)
  angledifference$diff2.13 <- ifelse(angledifference$diff2.13 > pi, c(abs(angledifference$diff2.13 - (2*pi))), c(angledifference$diff2.13))
  angledifference$diff2.14 <- abs(fly2$oriradians - angles$ang2.14)
  angledifference$diff2.14 <- ifelse(angledifference$diff2.14 > pi, c(abs(angledifference$diff2.14 - (2*pi))), c(angledifference$diff2.14))
  angledifference$diff2.15 <- abs(fly2$oriradians - angles$ang2.15)
  angledifference$diff2.15 <- ifelse(angledifference$diff2.15 > pi, c(abs(angledifference$diff2.15 - (2*pi))), c(angledifference$diff2.15))
  angledifference$diff2.16 <- abs(fly2$oriradians - angles$ang2.16)
  angledifference$diff2.16 <- ifelse(angledifference$diff2.16 > pi, c(abs(angledifference$diff2.16 - (2*pi))), c(angledifference$diff2.16))
  angledifference$diff2.17 <- abs(fly2$oriradians - angles$ang2.17)
  angledifference$diff2.17 <- ifelse(angledifference$diff2.17 > pi, c(abs(angledifference$diff2.17 - (2*pi))), c(angledifference$diff2.17))
  angledifference$diff2.18 <- abs(fly2$oriradians - angles$ang2.18)
  angledifference$diff2.18 <- ifelse(angledifference$diff2.18 > pi, c(abs(angledifference$diff2.18 - (2*pi))), c(angledifference$diff2.18))
  angledifference$diff2.19 <- abs(fly2$oriradians - angles$ang2.19)
  angledifference$diff2.19 <- ifelse(angledifference$diff2.19 > pi, c(abs(angledifference$diff2.19 - (2*pi))), c(angledifference$diff2.19))
  angledifference$diff2.20 <- abs(fly2$oriradians - angles$ang2.20)
  angledifference$diff2.20 <- ifelse(angledifference$diff2.20 > pi, c(abs(angledifference$diff2.20 - (2*pi))), c(angledifference$diff2.20))
  angledifference$diff3.1 <- abs(fly3$oriradians - angles$ang3.1)
  angledifference$diff3.1 <- ifelse(angledifference$diff3.1 > pi, c(abs(angledifference$diff3.1 - (2*pi))), c(angledifference$diff3.1))
  angledifference$diff3.2 <- abs(fly3$oriradians - angles$ang3.2)
  angledifference$diff3.2 <- ifelse(angledifference$diff3.2 > pi, c(abs(angledifference$diff3.2 - (2*pi))), c(angledifference$diff3.2))
  angledifference$diff3.4 <- abs(fly3$oriradians - angles$ang3.4)
  angledifference$diff3.4 <- ifelse(angledifference$diff3.4 > pi, c(abs(angledifference$diff3.4 - (2*pi))), c(angledifference$diff3.4))
  angledifference$diff3.5 <- abs(fly3$oriradians - angles$ang3.5)
  angledifference$diff3.5 <- ifelse(angledifference$diff3.5 > pi, c(abs(angledifference$diff3.5 - (2*pi))), c(angledifference$diff3.5))
  angledifference$diff3.6 <- abs(fly3$oriradians - angles$ang3.6)
  angledifference$diff3.6 <- ifelse(angledifference$diff3.6 > pi, c(abs(angledifference$diff3.6 - (2*pi))), c(angledifference$diff3.6))
  angledifference$diff3.7 <- abs(fly3$oriradians - angles$ang3.7)
  angledifference$diff3.7 <- ifelse(angledifference$diff3.7 > pi, c(abs(angledifference$diff3.7 - (2*pi))), c(angledifference$diff3.7))
  angledifference$diff3.8 <- abs(fly3$oriradians - angles$ang3.8)
  angledifference$diff3.8 <- ifelse(angledifference$diff3.8 > pi, c(abs(angledifference$diff3.8 - (2*pi))), c(angledifference$diff3.8))
  angledifference$diff3.9 <- abs(fly3$oriradians - angles$ang3.9)
  angledifference$diff3.9 <- ifelse(angledifference$diff3.9 > pi, c(abs(angledifference$diff3.9 - (2*pi))), c(angledifference$diff3.9))
  angledifference$diff3.10 <- abs(fly3$oriradians - angles$ang3.10)
  angledifference$diff3.10 <- ifelse(angledifference$diff3.10 > pi, c(abs(angledifference$diff3.10 - (2*pi))), c(angledifference$diff3.10))
  angledifference$diff3.11 <- abs(fly3$oriradians - angles$ang3.11)
  angledifference$diff3.11 <- ifelse(angledifference$diff3.11 > pi, c(abs(angledifference$diff3.11 - (2*pi))), c(angledifference$diff3.11))
  angledifference$diff3.12 <- abs(fly3$oriradians - angles$ang3.12)
  angledifference$diff3.12 <- ifelse(angledifference$diff3.12 > pi, c(abs(angledifference$diff3.12 - (2*pi))), c(angledifference$diff3.12))
  angledifference$diff3.13 <- abs(fly3$oriradians - angles$ang3.13)
  angledifference$diff3.13 <- ifelse(angledifference$diff3.13 > pi, c(abs(angledifference$diff3.13 - (2*pi))), c(angledifference$diff3.13))
  angledifference$diff3.14 <- abs(fly3$oriradians - angles$ang3.14)
  angledifference$diff3.14 <- ifelse(angledifference$diff3.14 > pi, c(abs(angledifference$diff3.14 - (2*pi))), c(angledifference$diff3.14))
  angledifference$diff3.15 <- abs(fly3$oriradians - angles$ang3.15)
  angledifference$diff3.15 <- ifelse(angledifference$diff3.15 > pi, c(abs(angledifference$diff3.15 - (2*pi))), c(angledifference$diff3.15))
  angledifference$diff3.16 <- abs(fly3$oriradians - angles$ang3.16)
  angledifference$diff3.16 <- ifelse(angledifference$diff3.16 > pi, c(abs(angledifference$diff3.16 - (2*pi))), c(angledifference$diff3.16))
  angledifference$diff3.17 <- abs(fly3$oriradians - angles$ang3.17)
  angledifference$diff3.17 <- ifelse(angledifference$diff3.17 > pi, c(abs(angledifference$diff3.17 - (2*pi))), c(angledifference$diff3.17))
  angledifference$diff3.18 <- abs(fly3$oriradians - angles$ang3.18)
  angledifference$diff3.18 <- ifelse(angledifference$diff3.18 > pi, c(abs(angledifference$diff3.18 - (2*pi))), c(angledifference$diff3.18))
  angledifference$diff3.19 <- abs(fly3$oriradians - angles$ang3.19)
  angledifference$diff3.19 <- ifelse(angledifference$diff3.19 > pi, c(abs(angledifference$diff3.19 - (2*pi))), c(angledifference$diff3.19))
  angledifference$diff3.20 <- abs(fly3$oriradians - angles$ang3.20)
  angledifference$diff3.20 <- ifelse(angledifference$diff3.20 > pi, c(abs(angledifference$diff3.20 - (2*pi))), c(angledifference$diff3.20))
  angledifference$diff4.1 <- abs(fly4$oriradians - angles$ang4.1)
  angledifference$diff4.1 <- ifelse(angledifference$diff4.1 > pi, c(abs(angledifference$diff4.1 - (2*pi))), c(angledifference$diff4.1))
  angledifference$diff4.2 <- abs(fly4$oriradians - angles$ang4.2)
  angledifference$diff4.2 <- ifelse(angledifference$diff4.2 > pi, c(abs(angledifference$diff4.2 - (2*pi))), c(angledifference$diff4.2))
  angledifference$diff4.3 <- abs(fly4$oriradians - angles$ang4.3)
  angledifference$diff4.3 <- ifelse(angledifference$diff4.3 > pi, c(abs(angledifference$diff4.3 - (2*pi))), c(angledifference$diff4.3))
  angledifference$diff4.5 <- abs(fly4$oriradians - angles$ang4.5)
  angledifference$diff4.5 <- ifelse(angledifference$diff4.5 > pi, c(abs(angledifference$diff4.5 - (2*pi))), c(angledifference$diff4.5))
  angledifference$diff4.6 <- abs(fly4$oriradians - angles$ang4.6)
  angledifference$diff4.6 <- ifelse(angledifference$diff4.6 > pi, c(abs(angledifference$diff4.6 - (2*pi))), c(angledifference$diff4.6))
  angledifference$diff4.7 <- abs(fly4$oriradians - angles$ang4.7)
  angledifference$diff4.7 <- ifelse(angledifference$diff4.7 > pi, c(abs(angledifference$diff4.7 - (2*pi))), c(angledifference$diff4.7))
  angledifference$diff4.8 <- abs(fly4$oriradians - angles$ang4.8)
  angledifference$diff4.8 <- ifelse(angledifference$diff4.8 > pi, c(abs(angledifference$diff4.8 - (2*pi))), c(angledifference$diff4.8))
  angledifference$diff4.9 <- abs(fly4$oriradians - angles$ang4.9)
  angledifference$diff4.9 <- ifelse(angledifference$diff4.9 > pi, c(abs(angledifference$diff4.9 - (2*pi))), c(angledifference$diff4.9))
  angledifference$diff4.10 <- abs(fly4$oriradians - angles$ang4.10)
  angledifference$diff4.10 <- ifelse(angledifference$diff4.10 > pi, c(abs(angledifference$diff4.10 - (2*pi))), c(angledifference$diff4.10))
  angledifference$diff4.11 <- abs(fly4$oriradians - angles$ang4.11)
  angledifference$diff4.11 <- ifelse(angledifference$diff4.11 > pi, c(abs(angledifference$diff4.11 - (2*pi))), c(angledifference$diff4.11))
  angledifference$diff4.12 <- abs(fly4$oriradians - angles$ang4.12)
  angledifference$diff4.12 <- ifelse(angledifference$diff4.12 > pi, c(abs(angledifference$diff4.12 - (2*pi))), c(angledifference$diff4.12))
  angledifference$diff4.13 <- abs(fly4$oriradians - angles$ang4.13)
  angledifference$diff4.13 <- ifelse(angledifference$diff4.13 > pi, c(abs(angledifference$diff4.13 - (2*pi))), c(angledifference$diff4.13))
  angledifference$diff4.14 <- abs(fly4$oriradians - angles$ang4.14)
  angledifference$diff4.14 <- ifelse(angledifference$diff4.14 > pi, c(abs(angledifference$diff4.14 - (2*pi))), c(angledifference$diff4.14))
  angledifference$diff4.15 <- abs(fly4$oriradians - angles$ang4.15)
  angledifference$diff4.15 <- ifelse(angledifference$diff4.15 > pi, c(abs(angledifference$diff4.15 - (2*pi))), c(angledifference$diff4.15))
  angledifference$diff4.16 <- abs(fly4$oriradians - angles$ang4.16)
  angledifference$diff4.16 <- ifelse(angledifference$diff4.16 > pi, c(abs(angledifference$diff4.16 - (2*pi))), c(angledifference$diff4.16))
  angledifference$diff4.17 <- abs(fly4$oriradians - angles$ang4.17)
  angledifference$diff4.17 <- ifelse(angledifference$diff4.17 > pi, c(abs(angledifference$diff4.17 - (2*pi))), c(angledifference$diff4.17))
  angledifference$diff4.18 <- abs(fly4$oriradians - angles$ang4.18)
  angledifference$diff4.18 <- ifelse(angledifference$diff4.18 > pi, c(abs(angledifference$diff4.18 - (2*pi))), c(angledifference$diff4.18))
  angledifference$diff4.19 <- abs(fly4$oriradians - angles$ang4.19)
  angledifference$diff4.19 <- ifelse(angledifference$diff4.19 > pi, c(abs(angledifference$diff4.19 - (2*pi))), c(angledifference$diff4.19))
  angledifference$diff4.20 <- abs(fly4$oriradians - angles$ang4.20)
  angledifference$diff4.20 <- ifelse(angledifference$diff4.20 > pi, c(abs(angledifference$diff4.20 - (2*pi))), c(angledifference$diff4.20))
  angledifference$diff5.1 <- abs(fly5$oriradians - angles$ang5.1)
  angledifference$diff5.1 <- ifelse(angledifference$diff5.1 > pi, c(abs(angledifference$diff5.1 - (2*pi))), c(angledifference$diff5.1))
  angledifference$diff5.2 <- abs(fly5$oriradians - angles$ang5.2)
  angledifference$diff5.2 <- ifelse(angledifference$diff5.2 > pi, c(abs(angledifference$diff5.2 - (2*pi))), c(angledifference$diff5.2))
  angledifference$diff5.3 <- abs(fly5$oriradians - angles$ang5.3)
  angledifference$diff5.3 <- ifelse(angledifference$diff5.3 > pi, c(abs(angledifference$diff5.3 - (2*pi))), c(angledifference$diff5.3))
  angledifference$diff5.4 <- abs(fly5$oriradians - angles$ang5.4)
  angledifference$diff5.4 <- ifelse(angledifference$diff5.4 > pi, c(abs(angledifference$diff5.4 - (2*pi))), c(angledifference$diff5.4))
  angledifference$diff5.6 <- abs(fly5$oriradians - angles$ang5.6)
  angledifference$diff5.6 <- ifelse(angledifference$diff5.6 > pi, c(abs(angledifference$diff5.6 - (2*pi))), c(angledifference$diff5.6))
  angledifference$diff5.7 <- abs(fly5$oriradians - angles$ang5.7)
  angledifference$diff5.7 <- ifelse(angledifference$diff5.7 > pi, c(abs(angledifference$diff5.7 - (2*pi))), c(angledifference$diff5.7))
  angledifference$diff5.8 <- abs(fly5$oriradians - angles$ang5.8)
  angledifference$diff5.8 <- ifelse(angledifference$diff5.8 > pi, c(abs(angledifference$diff5.8 - (2*pi))), c(angledifference$diff5.8))
  angledifference$diff5.9 <- abs(fly5$oriradians - angles$ang5.9)
  angledifference$diff5.9 <- ifelse(angledifference$diff5.9 > pi, c(abs(angledifference$diff5.9 - (2*pi))), c(angledifference$diff5.9))
  angledifference$diff5.10 <- abs(fly5$oriradians - angles$ang5.10)
  angledifference$diff5.10 <- ifelse(angledifference$diff5.10 > pi, c(abs(angledifference$diff5.10 - (2*pi))), c(angledifference$diff5.10))
  angledifference$diff5.11 <- abs(fly5$oriradians - angles$ang5.11)
  angledifference$diff5.11 <- ifelse(angledifference$diff5.11 > pi, c(abs(angledifference$diff5.11 - (2*pi))), c(angledifference$diff5.11))
  angledifference$diff5.12 <- abs(fly5$oriradians - angles$ang5.12)
  angledifference$diff5.12 <- ifelse(angledifference$diff5.12 > pi, c(abs(angledifference$diff5.12 - (2*pi))), c(angledifference$diff5.12))
  angledifference$diff5.13 <- abs(fly5$oriradians - angles$ang5.13)
  angledifference$diff5.13 <- ifelse(angledifference$diff5.13 > pi, c(abs(angledifference$diff5.13 - (2*pi))), c(angledifference$diff5.13))
  angledifference$diff5.14 <- abs(fly5$oriradians - angles$ang5.14)
  angledifference$diff5.14 <- ifelse(angledifference$diff5.14 > pi, c(abs(angledifference$diff5.14 - (2*pi))), c(angledifference$diff5.14))
  angledifference$diff5.15 <- abs(fly5$oriradians - angles$ang5.15)
  angledifference$diff5.15 <- ifelse(angledifference$diff5.15 > pi, c(abs(angledifference$diff5.15 - (2*pi))), c(angledifference$diff5.15))
  angledifference$diff5.16 <- abs(fly5$oriradians - angles$ang5.16)
  angledifference$diff5.16 <- ifelse(angledifference$diff5.16 > pi, c(abs(angledifference$diff5.16 - (2*pi))), c(angledifference$diff5.16))
  angledifference$diff5.17 <- abs(fly5$oriradians - angles$ang5.17)
  angledifference$diff5.17 <- ifelse(angledifference$diff5.17 > pi, c(abs(angledifference$diff5.17 - (2*pi))), c(angledifference$diff5.17))
  angledifference$diff5.18 <- abs(fly5$oriradians - angles$ang5.18)
  angledifference$diff5.18 <- ifelse(angledifference$diff5.18 > pi, c(abs(angledifference$diff5.18 - (2*pi))), c(angledifference$diff5.18))
  angledifference$diff5.19 <- abs(fly5$oriradians - angles$ang5.19)
  angledifference$diff5.19 <- ifelse(angledifference$diff5.19 > pi, c(abs(angledifference$diff5.19 - (2*pi))), c(angledifference$diff5.19))
  angledifference$diff5.20 <- abs(fly5$oriradians - angles$ang5.20)
  angledifference$diff5.20 <- ifelse(angledifference$diff5.20 > pi, c(abs(angledifference$diff5.20 - (2*pi))), c(angledifference$diff5.20))
  angledifference$diff6.1 <- abs(fly6$oriradians - angles$ang6.1)
  angledifference$diff6.1 <- ifelse(angledifference$diff6.1 > pi, c(abs(angledifference$diff6.1 - (2*pi))), c(angledifference$diff6.1))
  angledifference$diff6.2 <- abs(fly6$oriradians - angles$ang6.2)
  angledifference$diff6.2 <- ifelse(angledifference$diff6.2 > pi, c(abs(angledifference$diff6.2 - (2*pi))), c(angledifference$diff6.2))
  angledifference$diff6.3 <- abs(fly6$oriradians - angles$ang6.3)
  angledifference$diff6.3 <- ifelse(angledifference$diff6.3 > pi, c(abs(angledifference$diff6.3 - (2*pi))), c(angledifference$diff6.3))
  angledifference$diff6.4 <- abs(fly6$oriradians - angles$ang6.4)
  angledifference$diff6.4 <- ifelse(angledifference$diff6.4 > pi, c(abs(angledifference$diff6.4 - (2*pi))), c(angledifference$diff6.4))
  angledifference$diff6.5 <- abs(fly6$oriradians - angles$ang6.5)
  angledifference$diff6.5 <- ifelse(angledifference$diff6.5 > pi, c(abs(angledifference$diff6.5 - (2*pi))), c(angledifference$diff6.5))
  angledifference$diff6.7 <- abs(fly6$oriradians - angles$ang6.7)
  angledifference$diff6.7 <- ifelse(angledifference$diff6.7 > pi, c(abs(angledifference$diff6.7 - (2*pi))), c(angledifference$diff6.7))
  angledifference$diff6.8 <- abs(fly6$oriradians - angles$ang6.8)
  angledifference$diff6.8 <- ifelse(angledifference$diff6.8 > pi, c(abs(angledifference$diff6.8 - (2*pi))), c(angledifference$diff6.8))
  angledifference$diff6.9 <- abs(fly6$oriradians - angles$ang6.9)
  angledifference$diff6.9 <- ifelse(angledifference$diff6.9 > pi, c(abs(angledifference$diff6.9 - (2*pi))), c(angledifference$diff6.9))
  angledifference$diff6.10 <- abs(fly6$oriradians - angles$ang6.10)
  angledifference$diff6.10 <- ifelse(angledifference$diff6.10 > pi, c(abs(angledifference$diff6.10 - (2*pi))), c(angledifference$diff6.10))
  angledifference$diff6.11 <- abs(fly6$oriradians - angles$ang6.11)
  angledifference$diff6.11 <- ifelse(angledifference$diff6.11 > pi, c(abs(angledifference$diff6.11 - (2*pi))), c(angledifference$diff6.11))
  angledifference$diff6.12 <- abs(fly6$oriradians - angles$ang6.12)
  angledifference$diff6.12 <- ifelse(angledifference$diff6.12 > pi, c(abs(angledifference$diff6.12 - (2*pi))), c(angledifference$diff6.12))
  angledifference$diff6.13 <- abs(fly6$oriradians - angles$ang6.13)
  angledifference$diff6.13 <- ifelse(angledifference$diff6.13 > pi, c(abs(angledifference$diff6.13 - (2*pi))), c(angledifference$diff6.13))
  angledifference$diff6.14 <- abs(fly6$oriradians - angles$ang6.14)
  angledifference$diff6.14 <- ifelse(angledifference$diff6.14 > pi, c(abs(angledifference$diff6.14 - (2*pi))), c(angledifference$diff6.14))
  angledifference$diff6.15 <- abs(fly6$oriradians - angles$ang6.15)
  angledifference$diff6.15 <- ifelse(angledifference$diff6.15 > pi, c(abs(angledifference$diff6.15 - (2*pi))), c(angledifference$diff6.15))
  angledifference$diff6.16 <- abs(fly6$oriradians - angles$ang6.16)
  angledifference$diff6.16 <- ifelse(angledifference$diff6.16 > pi, c(abs(angledifference$diff6.16 - (2*pi))), c(angledifference$diff6.16))
  angledifference$diff6.17 <- abs(fly6$oriradians - angles$ang6.17)
  angledifference$diff6.17 <- ifelse(angledifference$diff6.17 > pi, c(abs(angledifference$diff6.17 - (2*pi))), c(angledifference$diff6.17))
  angledifference$diff6.18 <- abs(fly6$oriradians - angles$ang6.18)
  angledifference$diff6.18 <- ifelse(angledifference$diff6.18 > pi, c(abs(angledifference$diff6.18 - (2*pi))), c(angledifference$diff6.18))
  angledifference$diff6.19 <- abs(fly6$oriradians - angles$ang6.19)
  angledifference$diff6.19 <- ifelse(angledifference$diff6.19 > pi, c(abs(angledifference$diff6.19 - (2*pi))), c(angledifference$diff6.19))
  angledifference$diff6.20 <- abs(fly6$oriradians - angles$ang6.20)
  angledifference$diff6.20 <- ifelse(angledifference$diff6.20 > pi, c(abs(angledifference$diff6.20 - (2*pi))), c(angledifference$diff6.20))
  angledifference$diff7.1 <- abs(fly7$oriradians - angles$ang7.1)
  angledifference$diff7.1 <- ifelse(angledifference$diff7.1 > pi, c(abs(angledifference$diff7.1 - (2*pi))), c(angledifference$diff7.1))
  angledifference$diff7.2 <- abs(fly7$oriradians - angles$ang7.2)
  angledifference$diff7.2 <- ifelse(angledifference$diff7.2 > pi, c(abs(angledifference$diff7.2 - (2*pi))), c(angledifference$diff7.2))
  angledifference$diff7.3 <- abs(fly7$oriradians - angles$ang7.3)
  angledifference$diff7.3 <- ifelse(angledifference$diff7.3 > pi, c(abs(angledifference$diff7.3 - (2*pi))), c(angledifference$diff7.3))
  angledifference$diff7.4 <- abs(fly7$oriradians - angles$ang7.4)
  angledifference$diff7.4 <- ifelse(angledifference$diff7.4 > pi, c(abs(angledifference$diff7.4 - (2*pi))), c(angledifference$diff7.4))
  angledifference$diff7.5 <- abs(fly7$oriradians - angles$ang7.5)
  angledifference$diff7.5 <- ifelse(angledifference$diff7.5 > pi, c(abs(angledifference$diff7.5 - (2*pi))), c(angledifference$diff7.5))
  angledifference$diff7.6 <- abs(fly7$oriradians - angles$ang7.6)
  angledifference$diff7.6 <- ifelse(angledifference$diff7.6 > pi, c(abs(angledifference$diff7.6 - (2*pi))), c(angledifference$diff7.6))
  angledifference$diff7.8 <- abs(fly7$oriradians - angles$ang7.8)
  angledifference$diff7.8 <- ifelse(angledifference$diff7.8 > pi, c(abs(angledifference$diff7.8 - (2*pi))), c(angledifference$diff7.8))
  angledifference$diff7.9 <- abs(fly7$oriradians - angles$ang7.9)
  angledifference$diff7.9 <- ifelse(angledifference$diff7.9 > pi, c(abs(angledifference$diff7.9 - (2*pi))), c(angledifference$diff7.9))
  angledifference$diff7.10 <- abs(fly7$oriradians - angles$ang7.10)
  angledifference$diff7.10 <- ifelse(angledifference$diff7.10 > pi, c(abs(angledifference$diff7.10 - (2*pi))), c(angledifference$diff7.10))
  angledifference$diff7.11 <- abs(fly7$oriradians - angles$ang7.11)
  angledifference$diff7.11 <- ifelse(angledifference$diff7.11 > pi, c(abs(angledifference$diff7.11 - (2*pi))), c(angledifference$diff7.11))
  angledifference$diff7.12 <- abs(fly7$oriradians - angles$ang7.12)
  angledifference$diff7.12 <- ifelse(angledifference$diff7.12 > pi, c(abs(angledifference$diff7.12 - (2*pi))), c(angledifference$diff7.12))
  angledifference$diff7.13 <- abs(fly7$oriradians - angles$ang7.13)
  angledifference$diff7.13 <- ifelse(angledifference$diff7.13 > pi, c(abs(angledifference$diff7.13 - (2*pi))), c(angledifference$diff7.13))
  angledifference$diff7.14 <- abs(fly7$oriradians - angles$ang7.14)
  angledifference$diff7.14 <- ifelse(angledifference$diff7.14 > pi, c(abs(angledifference$diff7.14 - (2*pi))), c(angledifference$diff7.14))
  angledifference$diff7.15 <- abs(fly7$oriradians - angles$ang7.15)
  angledifference$diff7.15 <- ifelse(angledifference$diff7.15 > pi, c(abs(angledifference$diff7.15 - (2*pi))), c(angledifference$diff7.15))
  angledifference$diff7.16 <- abs(fly7$oriradians - angles$ang7.16)
  angledifference$diff7.16 <- ifelse(angledifference$diff7.16 > pi, c(abs(angledifference$diff7.16 - (2*pi))), c(angledifference$diff7.16))
  angledifference$diff7.17 <- abs(fly7$oriradians - angles$ang7.17)
  angledifference$diff7.17 <- ifelse(angledifference$diff7.17 > pi, c(abs(angledifference$diff7.17 - (2*pi))), c(angledifference$diff7.17))
  angledifference$diff7.18 <- abs(fly7$oriradians - angles$ang7.18)
  angledifference$diff7.18 <- ifelse(angledifference$diff7.18 > pi, c(abs(angledifference$diff7.18 - (2*pi))), c(angledifference$diff7.18))
  angledifference$diff7.19 <- abs(fly7$oriradians - angles$ang7.19)
  angledifference$diff7.19 <- ifelse(angledifference$diff7.19 > pi, c(abs(angledifference$diff7.19 - (2*pi))), c(angledifference$diff7.19))
  angledifference$diff7.20 <- abs(fly7$oriradians - angles$ang7.20)
  angledifference$diff7.20 <- ifelse(angledifference$diff7.20 > pi, c(abs(angledifference$diff7.20 - (2*pi))), c(angledifference$diff7.20))
  angledifference$diff8.1 <- abs(fly8$oriradians - angles$ang8.1)
  angledifference$diff8.1 <- ifelse(angledifference$diff8.1 > pi, c(abs(angledifference$diff8.1 - (2*pi))), c(angledifference$diff8.1))
  angledifference$diff8.2 <- abs(fly8$oriradians - angles$ang8.2)
  angledifference$diff8.2 <- ifelse(angledifference$diff8.2 > pi, c(abs(angledifference$diff8.2 - (2*pi))), c(angledifference$diff8.2))
  angledifference$diff8.3 <- abs(fly8$oriradians - angles$ang8.3)
  angledifference$diff8.3 <- ifelse(angledifference$diff8.3 > pi, c(abs(angledifference$diff8.3 - (2*pi))), c(angledifference$diff8.3))
  angledifference$diff8.4 <- abs(fly8$oriradians - angles$ang8.4)
  angledifference$diff8.4 <- ifelse(angledifference$diff8.4 > pi, c(abs(angledifference$diff8.4 - (2*pi))), c(angledifference$diff8.4))
  angledifference$diff8.5 <- abs(fly8$oriradians - angles$ang8.5)
  angledifference$diff8.5 <- ifelse(angledifference$diff8.5 > pi, c(abs(angledifference$diff8.5 - (2*pi))), c(angledifference$diff8.5))
  angledifference$diff8.6 <- abs(fly8$oriradians - angles$ang8.6)
  angledifference$diff8.6 <- ifelse(angledifference$diff8.6 > pi, c(abs(angledifference$diff8.6 - (2*pi))), c(angledifference$diff8.6))
  angledifference$diff8.7 <- abs(fly8$oriradians - angles$ang8.7)
  angledifference$diff8.7 <- ifelse(angledifference$diff8.7 > pi, c(abs(angledifference$diff8.7 - (2*pi))), c(angledifference$diff8.7))
  angledifference$diff8.9 <- abs(fly8$oriradians - angles$ang8.9)
  angledifference$diff8.9 <- ifelse(angledifference$diff8.9 > pi, c(abs(angledifference$diff8.9 - (2*pi))), c(angledifference$diff8.9))
  angledifference$diff8.10 <- abs(fly8$oriradians - angles$ang8.10)
  angledifference$diff8.10 <- ifelse(angledifference$diff8.10 > pi, c(abs(angledifference$diff8.10 - (2*pi))), c(angledifference$diff8.10))
  angledifference$diff8.11 <- abs(fly8$oriradians - angles$ang8.11)
  angledifference$diff8.11 <- ifelse(angledifference$diff8.11 > pi, c(abs(angledifference$diff8.11 - (2*pi))), c(angledifference$diff8.11))
  angledifference$diff8.12 <- abs(fly8$oriradians - angles$ang8.12)
  angledifference$diff8.12 <- ifelse(angledifference$diff8.12 > pi, c(abs(angledifference$diff8.12 - (2*pi))), c(angledifference$diff8.12))
  angledifference$diff8.13 <- abs(fly8$oriradians - angles$ang8.13)
  angledifference$diff8.13 <- ifelse(angledifference$diff8.13 > pi, c(abs(angledifference$diff8.13 - (2*pi))), c(angledifference$diff8.13))
  angledifference$diff8.14 <- abs(fly8$oriradians - angles$ang8.14)
  angledifference$diff8.14 <- ifelse(angledifference$diff8.14 > pi, c(abs(angledifference$diff8.14 - (2*pi))), c(angledifference$diff8.14))
  angledifference$diff8.15 <- abs(fly8$oriradians - angles$ang8.15)
  angledifference$diff8.15 <- ifelse(angledifference$diff8.15 > pi, c(abs(angledifference$diff8.15 - (2*pi))), c(angledifference$diff8.15))
  angledifference$diff8.16 <- abs(fly8$oriradians - angles$ang8.16)
  angledifference$diff8.16 <- ifelse(angledifference$diff8.16 > pi, c(abs(angledifference$diff8.16 - (2*pi))), c(angledifference$diff8.16))
  angledifference$diff8.17 <- abs(fly8$oriradians - angles$ang8.17)
  angledifference$diff8.17 <- ifelse(angledifference$diff8.17 > pi, c(abs(angledifference$diff8.17 - (2*pi))), c(angledifference$diff8.17))
  angledifference$diff8.18 <- abs(fly8$oriradians - angles$ang8.18)
  angledifference$diff8.18 <- ifelse(angledifference$diff8.18 > pi, c(abs(angledifference$diff8.18 - (2*pi))), c(angledifference$diff8.18))
  angledifference$diff8.19 <- abs(fly8$oriradians - angles$ang8.19)
  angledifference$diff8.19 <- ifelse(angledifference$diff8.19 > pi, c(abs(angledifference$diff8.19 - (2*pi))), c(angledifference$diff8.19))
  angledifference$diff8.20 <- abs(fly8$oriradians - angles$ang8.20)
  angledifference$diff8.20 <- ifelse(angledifference$diff8.20 > pi, c(abs(angledifference$diff8.20 - (2*pi))), c(angledifference$diff8.20))
  angledifference$diff9.1 <- abs(fly9$oriradians - angles$ang9.1)
  angledifference$diff9.1 <- ifelse(angledifference$diff9.1 > pi, c(abs(angledifference$diff9.1 - (2*pi))), c(angledifference$diff9.1))
  angledifference$diff9.2 <- abs(fly9$oriradians - angles$ang9.2)
  angledifference$diff9.2 <- ifelse(angledifference$diff9.2 > pi, c(abs(angledifference$diff9.2 - (2*pi))), c(angledifference$diff9.2))
  angledifference$diff9.3 <- abs(fly9$oriradians - angles$ang9.3)
  angledifference$diff9.3 <- ifelse(angledifference$diff9.3 > pi, c(abs(angledifference$diff9.3 - (2*pi))), c(angledifference$diff9.3))
  angledifference$diff9.4 <- abs(fly9$oriradians - angles$ang9.4)
  angledifference$diff9.4 <- ifelse(angledifference$diff9.4 > pi, c(abs(angledifference$diff9.4 - (2*pi))), c(angledifference$diff9.4))
  angledifference$diff9.5 <- abs(fly9$oriradians - angles$ang9.5)
  angledifference$diff9.5 <- ifelse(angledifference$diff9.5 > pi, c(abs(angledifference$diff9.5 - (2*pi))), c(angledifference$diff9.5))
  angledifference$diff9.6 <- abs(fly9$oriradians - angles$ang9.6)
  angledifference$diff9.6 <- ifelse(angledifference$diff9.6 > pi, c(abs(angledifference$diff9.6 - (2*pi))), c(angledifference$diff9.6))
  angledifference$diff9.7 <- abs(fly9$oriradians - angles$ang9.7)
  angledifference$diff9.7 <- ifelse(angledifference$diff9.7 > pi, c(abs(angledifference$diff9.7 - (2*pi))), c(angledifference$diff9.7))
  angledifference$diff9.8 <- abs(fly9$oriradians - angles$ang9.8)
  angledifference$diff9.8 <- ifelse(angledifference$diff9.8 > pi, c(abs(angledifference$diff9.8 - (2*pi))), c(angledifference$diff9.8))
  angledifference$diff9.10 <- abs(fly9$oriradians - angles$ang9.10)
  angledifference$diff9.10 <- ifelse(angledifference$diff9.10 > pi, c(abs(angledifference$diff9.10 - (2*pi))), c(angledifference$diff9.10))
  angledifference$diff9.11 <- abs(fly9$oriradians - angles$ang9.11)
  angledifference$diff9.11 <- ifelse(angledifference$diff9.11 > pi, c(abs(angledifference$diff9.11 - (2*pi))), c(angledifference$diff9.11))
  angledifference$diff9.12 <- abs(fly9$oriradians - angles$ang9.12)
  angledifference$diff9.12 <- ifelse(angledifference$diff9.12 > pi, c(abs(angledifference$diff9.12 - (2*pi))), c(angledifference$diff9.12))
  angledifference$diff9.13 <- abs(fly9$oriradians - angles$ang9.13)
  angledifference$diff9.13 <- ifelse(angledifference$diff9.13 > pi, c(abs(angledifference$diff9.13 - (2*pi))), c(angledifference$diff9.13))
  angledifference$diff9.14 <- abs(fly9$oriradians - angles$ang9.14)
  angledifference$diff9.14 <- ifelse(angledifference$diff9.14 > pi, c(abs(angledifference$diff9.14 - (2*pi))), c(angledifference$diff9.14))
  angledifference$diff9.15 <- abs(fly9$oriradians - angles$ang9.15)
  angledifference$diff9.15 <- ifelse(angledifference$diff9.15 > pi, c(abs(angledifference$diff9.15 - (2*pi))), c(angledifference$diff9.15))
  angledifference$diff9.16 <- abs(fly9$oriradians - angles$ang9.16)
  angledifference$diff9.16 <- ifelse(angledifference$diff9.16 > pi, c(abs(angledifference$diff9.16 - (2*pi))), c(angledifference$diff9.16))
  angledifference$diff9.17 <- abs(fly9$oriradians - angles$ang9.17)
  angledifference$diff9.17 <- ifelse(angledifference$diff9.17 > pi, c(abs(angledifference$diff9.17 - (2*pi))), c(angledifference$diff9.17))
  angledifference$diff9.18 <- abs(fly9$oriradians - angles$ang9.18)
  angledifference$diff9.18 <- ifelse(angledifference$diff9.18 > pi, c(abs(angledifference$diff9.18 - (2*pi))), c(angledifference$diff9.18))
  angledifference$diff9.19 <- abs(fly9$oriradians - angles$ang9.19)
  angledifference$diff9.19 <- ifelse(angledifference$diff9.19 > pi, c(abs(angledifference$diff9.19 - (2*pi))), c(angledifference$diff9.19))
  angledifference$diff9.20 <- abs(fly9$oriradians - angles$ang9.20)
  angledifference$diff9.20 <- ifelse(angledifference$diff9.20 > pi, c(abs(angledifference$diff9.20 - (2*pi))), c(angledifference$diff9.20))
  angledifference$diff10.1 <- abs(fly10$oriradians - angles$ang10.1)
  angledifference$diff10.1 <- ifelse(angledifference$diff10.1 > pi, c(abs(angledifference$diff10.1 - (2*pi))), c(angledifference$diff10.1))
  angledifference$diff10.2 <- abs(fly10$oriradians - angles$ang10.2)
  angledifference$diff10.2 <- ifelse(angledifference$diff10.2 > pi, c(abs(angledifference$diff10.2 - (2*pi))), c(angledifference$diff10.2))
  angledifference$diff10.3 <- abs(fly10$oriradians - angles$ang10.3)
  angledifference$diff10.3 <- ifelse(angledifference$diff10.3 > pi, c(abs(angledifference$diff10.3 - (2*pi))), c(angledifference$diff10.3))
  angledifference$diff10.4 <- abs(fly10$oriradians - angles$ang10.4)
  angledifference$diff10.4 <- ifelse(angledifference$diff10.4 > pi, c(abs(angledifference$diff10.4 - (2*pi))), c(angledifference$diff10.4))
  angledifference$diff10.5 <- abs(fly10$oriradians - angles$ang10.5)
  angledifference$diff10.5 <- ifelse(angledifference$diff10.5 > pi, c(abs(angledifference$diff10.5 - (2*pi))), c(angledifference$diff10.5))
  angledifference$diff10.6 <- abs(fly10$oriradians - angles$ang10.6)
  angledifference$diff10.6 <- ifelse(angledifference$diff10.6 > pi, c(abs(angledifference$diff10.6 - (2*pi))), c(angledifference$diff10.6))
  angledifference$diff10.7 <- abs(fly10$oriradians - angles$ang10.7)
  angledifference$diff10.7 <- ifelse(angledifference$diff10.7 > pi, c(abs(angledifference$diff10.7 - (2*pi))), c(angledifference$diff10.7))
  angledifference$diff10.8 <- abs(fly10$oriradians - angles$ang10.8)
  angledifference$diff10.8 <- ifelse(angledifference$diff10.8 > pi, c(abs(angledifference$diff10.8 - (2*pi))), c(angledifference$diff10.8))
  angledifference$diff10.9 <- abs(fly10$oriradians - angles$ang10.9)
  angledifference$diff10.9 <- ifelse(angledifference$diff10.9 > pi, c(abs(angledifference$diff10.9 - (2*pi))), c(angledifference$diff10.9))
  angledifference$diff10.11 <- abs(fly10$oriradians - angles$ang10.11)
  angledifference$diff10.11 <- ifelse(angledifference$diff10.11 > pi, c(abs(angledifference$diff10.11 - (2*pi))), c(angledifference$diff10.11))
  angledifference$diff10.12 <- abs(fly10$oriradians - angles$ang10.12)
  angledifference$diff10.12 <- ifelse(angledifference$diff10.12 > pi, c(abs(angledifference$diff10.12 - (2*pi))), c(angledifference$diff10.12))
  angledifference$diff10.13 <- abs(fly10$oriradians - angles$ang10.13)
  angledifference$diff10.13 <- ifelse(angledifference$diff10.13 > pi, c(abs(angledifference$diff10.13 - (2*pi))), c(angledifference$diff10.13))
  angledifference$diff10.14 <- abs(fly10$oriradians - angles$ang10.14)
  angledifference$diff10.14 <- ifelse(angledifference$diff10.14 > pi, c(abs(angledifference$diff10.14 - (2*pi))), c(angledifference$diff10.14))
  angledifference$diff10.15 <- abs(fly10$oriradians - angles$ang10.15)
  angledifference$diff10.15 <- ifelse(angledifference$diff10.15 > pi, c(abs(angledifference$diff10.15 - (2*pi))), c(angledifference$diff10.15))
  angledifference$diff10.16 <- abs(fly10$oriradians - angles$ang10.16)
  angledifference$diff10.16 <- ifelse(angledifference$diff10.16 > pi, c(abs(angledifference$diff10.16 - (2*pi))), c(angledifference$diff10.16))
  angledifference$diff10.17 <- abs(fly10$oriradians - angles$ang10.17)
  angledifference$diff10.17 <- ifelse(angledifference$diff10.17 > pi, c(abs(angledifference$diff10.17 - (2*pi))), c(angledifference$diff10.17))
  angledifference$diff10.18 <- abs(fly10$oriradians - angles$ang10.18)
  angledifference$diff10.18 <- ifelse(angledifference$diff10.18 > pi, c(abs(angledifference$diff10.18 - (2*pi))), c(angledifference$diff10.18))
  angledifference$diff10.19 <- abs(fly10$oriradians - angles$ang10.19)
  angledifference$diff10.19 <- ifelse(angledifference$diff10.19 > pi, c(abs(angledifference$diff10.19 - (2*pi))), c(angledifference$diff10.19))
  angledifference$diff10.20 <- abs(fly10$oriradians - angles$ang10.20)
  angledifference$diff10.20 <- ifelse(angledifference$diff10.20 > pi, c(abs(angledifference$diff10.20 - (2*pi))), c(angledifference$diff10.20))
  angledifference$diff11.1 <- abs(fly11$oriradians - angles$ang11.1)
  angledifference$diff11.1 <- ifelse(angledifference$diff11.1 > pi, c(abs(angledifference$diff11.1 - (2*pi))), c(angledifference$diff11.1))
  angledifference$diff11.2 <- abs(fly11$oriradians - angles$ang11.2)
  angledifference$diff11.2 <- ifelse(angledifference$diff11.2 > pi, c(abs(angledifference$diff11.2 - (2*pi))), c(angledifference$diff11.2))
  angledifference$diff11.3 <- abs(fly11$oriradians - angles$ang11.3)
  angledifference$diff11.3 <- ifelse(angledifference$diff11.3 > pi, c(abs(angledifference$diff11.3 - (2*pi))), c(angledifference$diff11.3))
  angledifference$diff11.4 <- abs(fly11$oriradians - angles$ang11.4)
  angledifference$diff11.4 <- ifelse(angledifference$diff11.4 > pi, c(abs(angledifference$diff11.4 - (2*pi))), c(angledifference$diff11.4))
  angledifference$diff11.5 <- abs(fly11$oriradians - angles$ang11.5)
  angledifference$diff11.5 <- ifelse(angledifference$diff11.5 > pi, c(abs(angledifference$diff11.5 - (2*pi))), c(angledifference$diff11.5))
  angledifference$diff11.6 <- abs(fly11$oriradians - angles$ang11.6)
  angledifference$diff11.6 <- ifelse(angledifference$diff11.6 > pi, c(abs(angledifference$diff11.6 - (2*pi))), c(angledifference$diff11.6))
  angledifference$diff11.7 <- abs(fly11$oriradians - angles$ang11.7)
  angledifference$diff11.7 <- ifelse(angledifference$diff11.7 > pi, c(abs(angledifference$diff11.7 - (2*pi))), c(angledifference$diff11.7))
  angledifference$diff11.8 <- abs(fly11$oriradians - angles$ang11.8)
  angledifference$diff11.8 <- ifelse(angledifference$diff11.8 > pi, c(abs(angledifference$diff11.8 - (2*pi))), c(angledifference$diff11.8))
  angledifference$diff11.9 <- abs(fly11$oriradians - angles$ang11.9)
  angledifference$diff11.9 <- ifelse(angledifference$diff11.9 > pi, c(abs(angledifference$diff11.9 - (2*pi))), c(angledifference$diff11.9))
  angledifference$diff11.10 <- abs(fly11$oriradians - angles$ang11.10)
  angledifference$diff11.10 <- ifelse(angledifference$diff11.10 > pi, c(abs(angledifference$diff11.10 - (2*pi))), c(angledifference$diff11.10))
  angledifference$diff11.12 <- abs(fly11$oriradians - angles$ang11.12)
  angledifference$diff11.12 <- ifelse(angledifference$diff11.12 > pi, c(abs(angledifference$diff11.12 - (2*pi))), c(angledifference$diff11.12))
  angledifference$diff11.13 <- abs(fly11$oriradians - angles$ang11.13)
  angledifference$diff11.13 <- ifelse(angledifference$diff11.13 > pi, c(abs(angledifference$diff11.13 - (2*pi))), c(angledifference$diff11.13))
  angledifference$diff11.14 <- abs(fly11$oriradians - angles$ang11.14)
  angledifference$diff11.14 <- ifelse(angledifference$diff11.14 > pi, c(abs(angledifference$diff11.14 - (2*pi))), c(angledifference$diff11.14))
  angledifference$diff11.15 <- abs(fly11$oriradians - angles$ang11.15)
  angledifference$diff11.15 <- ifelse(angledifference$diff11.15 > pi, c(abs(angledifference$diff11.15 - (2*pi))), c(angledifference$diff11.15))
  angledifference$diff11.16 <- abs(fly11$oriradians - angles$ang11.16)
  angledifference$diff11.16 <- ifelse(angledifference$diff11.16 > pi, c(abs(angledifference$diff11.16 - (2*pi))), c(angledifference$diff11.16))
  angledifference$diff11.17 <- abs(fly11$oriradians - angles$ang11.17)
  angledifference$diff11.17 <- ifelse(angledifference$diff11.17 > pi, c(abs(angledifference$diff11.17 - (2*pi))), c(angledifference$diff11.17))
  angledifference$diff11.18 <- abs(fly11$oriradians - angles$ang11.18)
  angledifference$diff11.18 <- ifelse(angledifference$diff11.18 > pi, c(abs(angledifference$diff11.18 - (2*pi))), c(angledifference$diff11.18))
  angledifference$diff11.19 <- abs(fly11$oriradians - angles$ang11.19)
  angledifference$diff11.19 <- ifelse(angledifference$diff11.19 > pi, c(abs(angledifference$diff11.19 - (2*pi))), c(angledifference$diff11.19))
  angledifference$diff11.20 <- abs(fly11$oriradians - angles$ang11.20)
  angledifference$diff11.20 <- ifelse(angledifference$diff11.20 > pi, c(abs(angledifference$diff11.20 - (2*pi))), c(angledifference$diff11.20))
  angledifference$diff12.1 <- abs(fly12$oriradians - angles$ang12.1)
  angledifference$diff12.1 <- ifelse(angledifference$diff12.1 > pi, c(abs(angledifference$diff12.1 - (2*pi))), c(angledifference$diff12.1))
  angledifference$diff12.2 <- abs(fly12$oriradians - angles$ang12.2)
  angledifference$diff12.2 <- ifelse(angledifference$diff12.2 > pi, c(abs(angledifference$diff12.2 - (2*pi))), c(angledifference$diff12.2))
  angledifference$diff12.3 <- abs(fly12$oriradians - angles$ang12.3)
  angledifference$diff12.3 <- ifelse(angledifference$diff12.3 > pi, c(abs(angledifference$diff12.3 - (2*pi))), c(angledifference$diff12.3))
  angledifference$diff12.4 <- abs(fly12$oriradians - angles$ang12.4)
  angledifference$diff12.4 <- ifelse(angledifference$diff12.4 > pi, c(abs(angledifference$diff12.4 - (2*pi))), c(angledifference$diff12.4))
  angledifference$diff12.5 <- abs(fly12$oriradians - angles$ang12.5)
  angledifference$diff12.5 <- ifelse(angledifference$diff12.5 > pi, c(abs(angledifference$diff12.5 - (2*pi))), c(angledifference$diff12.5))
  angledifference$diff12.6 <- abs(fly12$oriradians - angles$ang12.6)
  angledifference$diff12.6 <- ifelse(angledifference$diff12.6 > pi, c(abs(angledifference$diff12.6 - (2*pi))), c(angledifference$diff12.6))
  angledifference$diff12.7 <- abs(fly12$oriradians - angles$ang12.7)
  angledifference$diff12.7 <- ifelse(angledifference$diff12.7 > pi, c(abs(angledifference$diff12.7 - (2*pi))), c(angledifference$diff12.7))
  angledifference$diff12.8 <- abs(fly12$oriradians - angles$ang12.8)
  angledifference$diff12.8 <- ifelse(angledifference$diff12.8 > pi, c(abs(angledifference$diff12.8 - (2*pi))), c(angledifference$diff12.8))
  angledifference$diff12.9 <- abs(fly12$oriradians - angles$ang12.9)
  angledifference$diff12.9 <- ifelse(angledifference$diff12.9 > pi, c(abs(angledifference$diff12.9 - (2*pi))), c(angledifference$diff12.9))
  angledifference$diff12.10 <- abs(fly12$oriradians - angles$ang12.10)
  angledifference$diff12.10 <- ifelse(angledifference$diff12.10 > pi, c(abs(angledifference$diff12.10 - (2*pi))), c(angledifference$diff12.10))
  angledifference$diff12.11 <- abs(fly12$oriradians - angles$ang12.11)
  angledifference$diff12.11 <- ifelse(angledifference$diff12.11 > pi, c(abs(angledifference$diff12.11 - (2*pi))), c(angledifference$diff12.11))
  angledifference$diff12.13 <- abs(fly12$oriradians - angles$ang12.13)
  angledifference$diff12.13 <- ifelse(angledifference$diff12.13 > pi, c(abs(angledifference$diff12.13 - (2*pi))), c(angledifference$diff12.13))
  angledifference$diff12.14 <- abs(fly12$oriradians - angles$ang12.14)
  angledifference$diff12.14 <- ifelse(angledifference$diff12.14 > pi, c(abs(angledifference$diff12.14 - (2*pi))), c(angledifference$diff12.14))
  angledifference$diff12.15 <- abs(fly12$oriradians - angles$ang12.15)
  angledifference$diff12.15 <- ifelse(angledifference$diff12.15 > pi, c(abs(angledifference$diff12.15 - (2*pi))), c(angledifference$diff12.15))
  angledifference$diff12.16 <- abs(fly12$oriradians - angles$ang12.16)
  angledifference$diff12.16 <- ifelse(angledifference$diff12.16 > pi, c(abs(angledifference$diff12.16 - (2*pi))), c(angledifference$diff12.16))
  angledifference$diff12.17 <- abs(fly12$oriradians - angles$ang12.17)
  angledifference$diff12.17 <- ifelse(angledifference$diff12.17 > pi, c(abs(angledifference$diff12.17 - (2*pi))), c(angledifference$diff12.17))
  angledifference$diff12.18 <- abs(fly12$oriradians - angles$ang12.18)
  angledifference$diff12.18 <- ifelse(angledifference$diff12.18 > pi, c(abs(angledifference$diff12.18 - (2*pi))), c(angledifference$diff12.18))
  angledifference$diff12.19 <- abs(fly12$oriradians - angles$ang12.19)
  angledifference$diff12.19 <- ifelse(angledifference$diff12.19 > pi, c(abs(angledifference$diff12.19 - (2*pi))), c(angledifference$diff12.19))
  angledifference$diff12.20 <- abs(fly12$oriradians - angles$ang12.20)
  angledifference$diff12.20 <- ifelse(angledifference$diff12.20 > pi, c(abs(angledifference$diff12.20 - (2*pi))), c(angledifference$diff12.20))
  angledifference$diff13.1 <- abs(fly13$oriradians - angles$ang13.1)
  angledifference$diff13.1 <- ifelse(angledifference$diff13.1 > pi, c(abs(angledifference$diff13.1 - (2*pi))), c(angledifference$diff13.1))
  angledifference$diff13.2 <- abs(fly13$oriradians - angles$ang13.2)
  angledifference$diff13.2 <- ifelse(angledifference$diff13.2 > pi, c(abs(angledifference$diff13.2 - (2*pi))), c(angledifference$diff13.2))
  angledifference$diff13.3 <- abs(fly13$oriradians - angles$ang13.3)
  angledifference$diff13.3 <- ifelse(angledifference$diff13.3 > pi, c(abs(angledifference$diff13.3 - (2*pi))), c(angledifference$diff13.3))
  angledifference$diff13.4 <- abs(fly13$oriradians - angles$ang13.4)
  angledifference$diff13.4 <- ifelse(angledifference$diff13.4 > pi, c(abs(angledifference$diff13.4 - (2*pi))), c(angledifference$diff13.4))
  angledifference$diff13.5 <- abs(fly13$oriradians - angles$ang13.5)
  angledifference$diff13.5 <- ifelse(angledifference$diff13.5 > pi, c(abs(angledifference$diff13.5 - (2*pi))), c(angledifference$diff13.5))
  angledifference$diff13.6 <- abs(fly13$oriradians - angles$ang13.6)
  angledifference$diff13.6 <- ifelse(angledifference$diff13.6 > pi, c(abs(angledifference$diff13.6 - (2*pi))), c(angledifference$diff13.6))
  angledifference$diff13.7 <- abs(fly13$oriradians - angles$ang13.7)
  angledifference$diff13.7 <- ifelse(angledifference$diff13.7 > pi, c(abs(angledifference$diff13.7 - (2*pi))), c(angledifference$diff13.7))
  angledifference$diff13.8 <- abs(fly13$oriradians - angles$ang13.8)
  angledifference$diff13.8 <- ifelse(angledifference$diff13.8 > pi, c(abs(angledifference$diff13.8 - (2*pi))), c(angledifference$diff13.8))
  angledifference$diff13.9 <- abs(fly13$oriradians - angles$ang13.9)
  angledifference$diff13.9 <- ifelse(angledifference$diff13.9 > pi, c(abs(angledifference$diff13.9 - (2*pi))), c(angledifference$diff13.9))
  angledifference$diff13.10 <- abs(fly13$oriradians - angles$ang13.10)
  angledifference$diff13.10 <- ifelse(angledifference$diff13.10 > pi, c(abs(angledifference$diff13.10 - (2*pi))), c(angledifference$diff13.10))
  angledifference$diff13.11 <- abs(fly13$oriradians - angles$ang13.11)
  angledifference$diff13.11 <- ifelse(angledifference$diff13.11 > pi, c(abs(angledifference$diff13.11 - (2*pi))), c(angledifference$diff13.11))
  angledifference$diff13.12 <- abs(fly13$oriradians - angles$ang13.12)
  angledifference$diff13.12 <- ifelse(angledifference$diff13.12 > pi, c(abs(angledifference$diff13.12 - (2*pi))), c(angledifference$diff13.12))
  angledifference$diff13.14 <- abs(fly13$oriradians - angles$ang13.14)
  angledifference$diff13.14 <- ifelse(angledifference$diff13.14 > pi, c(abs(angledifference$diff13.14 - (2*pi))), c(angledifference$diff13.14))
  angledifference$diff13.15 <- abs(fly13$oriradians - angles$ang13.15)
  angledifference$diff13.15 <- ifelse(angledifference$diff13.15 > pi, c(abs(angledifference$diff13.15 - (2*pi))), c(angledifference$diff13.15))
  angledifference$diff13.16 <- abs(fly13$oriradians - angles$ang13.16)
  angledifference$diff13.16 <- ifelse(angledifference$diff13.16 > pi, c(abs(angledifference$diff13.16 - (2*pi))), c(angledifference$diff13.16))
  angledifference$diff13.17 <- abs(fly13$oriradians - angles$ang13.17)
  angledifference$diff13.17 <- ifelse(angledifference$diff13.17 > pi, c(abs(angledifference$diff13.17 - (2*pi))), c(angledifference$diff13.17))
  angledifference$diff13.18 <- abs(fly13$oriradians - angles$ang13.18)
  angledifference$diff13.18 <- ifelse(angledifference$diff13.18 > pi, c(abs(angledifference$diff13.18 - (2*pi))), c(angledifference$diff13.18))
  angledifference$diff13.19 <- abs(fly13$oriradians - angles$ang13.19)
  angledifference$diff13.19 <- ifelse(angledifference$diff13.19 > pi, c(abs(angledifference$diff13.19 - (2*pi))), c(angledifference$diff13.19))
  angledifference$diff13.20 <- abs(fly13$oriradians - angles$ang13.20)
  angledifference$diff13.20 <- ifelse(angledifference$diff13.20 > pi, c(abs(angledifference$diff13.20 - (2*pi))), c(angledifference$diff13.20))
  angledifference$diff14.1 <- abs(fly14$oriradians - angles$ang14.1)
  angledifference$diff14.1 <- ifelse(angledifference$diff14.1 > pi, c(abs(angledifference$diff14.1 - (2*pi))), c(angledifference$diff14.1))
  angledifference$diff14.2 <- abs(fly14$oriradians - angles$ang14.2)
  angledifference$diff14.2 <- ifelse(angledifference$diff14.2 > pi, c(abs(angledifference$diff14.2 - (2*pi))), c(angledifference$diff14.2))
  angledifference$diff14.3 <- abs(fly14$oriradians - angles$ang14.3)
  angledifference$diff14.3 <- ifelse(angledifference$diff14.3 > pi, c(abs(angledifference$diff14.3 - (2*pi))), c(angledifference$diff14.3))
  angledifference$diff14.4 <- abs(fly14$oriradians - angles$ang14.4)
  angledifference$diff14.4 <- ifelse(angledifference$diff14.4 > pi, c(abs(angledifference$diff14.4 - (2*pi))), c(angledifference$diff14.4))
  angledifference$diff14.5 <- abs(fly14$oriradians - angles$ang14.5)
  angledifference$diff14.5 <- ifelse(angledifference$diff14.5 > pi, c(abs(angledifference$diff14.5 - (2*pi))), c(angledifference$diff14.5))
  angledifference$diff14.6 <- abs(fly14$oriradians - angles$ang14.6)
  angledifference$diff14.6 <- ifelse(angledifference$diff14.6 > pi, c(abs(angledifference$diff14.6 - (2*pi))), c(angledifference$diff14.6))
  angledifference$diff14.7 <- abs(fly14$oriradians - angles$ang14.7)
  angledifference$diff14.7 <- ifelse(angledifference$diff14.7 > pi, c(abs(angledifference$diff14.7 - (2*pi))), c(angledifference$diff14.7))
  angledifference$diff14.8 <- abs(fly14$oriradians - angles$ang14.8)
  angledifference$diff14.8 <- ifelse(angledifference$diff14.8 > pi, c(abs(angledifference$diff14.8 - (2*pi))), c(angledifference$diff14.8))
  angledifference$diff14.9 <- abs(fly14$oriradians - angles$ang14.9)
  angledifference$diff14.9 <- ifelse(angledifference$diff14.9 > pi, c(abs(angledifference$diff14.9 - (2*pi))), c(angledifference$diff14.9))
  angledifference$diff14.10 <- abs(fly14$oriradians - angles$ang14.10)
  angledifference$diff14.10 <- ifelse(angledifference$diff14.10 > pi, c(abs(angledifference$diff14.10 - (2*pi))), c(angledifference$diff14.10))
  angledifference$diff14.11 <- abs(fly14$oriradians - angles$ang14.11)
  angledifference$diff14.11 <- ifelse(angledifference$diff14.11 > pi, c(abs(angledifference$diff14.11 - (2*pi))), c(angledifference$diff14.11))
  angledifference$diff14.12 <- abs(fly14$oriradians - angles$ang14.12)
  angledifference$diff14.12 <- ifelse(angledifference$diff14.12 > pi, c(abs(angledifference$diff14.12 - (2*pi))), c(angledifference$diff14.12))
  angledifference$diff14.13 <- abs(fly14$oriradians - angles$ang14.13)
  angledifference$diff14.13 <- ifelse(angledifference$diff14.13 > pi, c(abs(angledifference$diff14.13 - (2*pi))), c(angledifference$diff14.13))
  angledifference$diff14.15 <- abs(fly14$oriradians - angles$ang14.15)
  angledifference$diff14.15 <- ifelse(angledifference$diff14.15 > pi, c(abs(angledifference$diff14.15 - (2*pi))), c(angledifference$diff14.15))
  angledifference$diff14.16 <- abs(fly14$oriradians - angles$ang14.16)
  angledifference$diff14.16 <- ifelse(angledifference$diff14.16 > pi, c(abs(angledifference$diff14.16 - (2*pi))), c(angledifference$diff14.16))
  angledifference$diff14.17 <- abs(fly14$oriradians - angles$ang14.17)
  angledifference$diff14.17 <- ifelse(angledifference$diff14.17 > pi, c(abs(angledifference$diff14.17 - (2*pi))), c(angledifference$diff14.17))
  angledifference$diff14.18 <- abs(fly14$oriradians - angles$ang14.18)
  angledifference$diff14.18 <- ifelse(angledifference$diff14.18 > pi, c(abs(angledifference$diff14.18 - (2*pi))), c(angledifference$diff14.18))
  angledifference$diff14.19 <- abs(fly14$oriradians - angles$ang14.19)
  angledifference$diff14.19 <- ifelse(angledifference$diff14.19 > pi, c(abs(angledifference$diff14.19 - (2*pi))), c(angledifference$diff14.19))
  angledifference$diff14.20 <- abs(fly14$oriradians - angles$ang14.20)
  angledifference$diff14.20 <- ifelse(angledifference$diff14.20 > pi, c(abs(angledifference$diff14.20 - (2*pi))), c(angledifference$diff14.20))
  angledifference$diff15.1 <- abs(fly15$oriradians - angles$ang15.1)
  angledifference$diff15.1 <- ifelse(angledifference$diff15.1 > pi, c(abs(angledifference$diff15.1 - (2*pi))), c(angledifference$diff15.1))
  angledifference$diff15.2 <- abs(fly15$oriradians - angles$ang15.2)
  angledifference$diff15.2 <- ifelse(angledifference$diff15.2 > pi, c(abs(angledifference$diff15.2 - (2*pi))), c(angledifference$diff15.2))
  angledifference$diff15.3 <- abs(fly15$oriradians - angles$ang15.3)
  angledifference$diff15.3 <- ifelse(angledifference$diff15.3 > pi, c(abs(angledifference$diff15.3 - (2*pi))), c(angledifference$diff15.3))
  angledifference$diff15.4 <- abs(fly15$oriradians - angles$ang15.4)
  angledifference$diff15.4 <- ifelse(angledifference$diff15.4 > pi, c(abs(angledifference$diff15.4 - (2*pi))), c(angledifference$diff15.4))
  angledifference$diff15.5 <- abs(fly15$oriradians - angles$ang15.5)
  angledifference$diff15.5 <- ifelse(angledifference$diff15.5 > pi, c(abs(angledifference$diff15.5 - (2*pi))), c(angledifference$diff15.5))
  angledifference$diff15.6 <- abs(fly15$oriradians - angles$ang15.6)
  angledifference$diff15.6 <- ifelse(angledifference$diff15.6 > pi, c(abs(angledifference$diff15.6 - (2*pi))), c(angledifference$diff15.6))
  angledifference$diff15.7 <- abs(fly15$oriradians - angles$ang15.7)
  angledifference$diff15.7 <- ifelse(angledifference$diff15.7 > pi, c(abs(angledifference$diff15.7 - (2*pi))), c(angledifference$diff15.7))
  angledifference$diff15.8 <- abs(fly15$oriradians - angles$ang15.8)
  angledifference$diff15.8 <- ifelse(angledifference$diff15.8 > pi, c(abs(angledifference$diff15.8 - (2*pi))), c(angledifference$diff15.8))
  angledifference$diff15.9 <- abs(fly15$oriradians - angles$ang15.9)
  angledifference$diff15.9 <- ifelse(angledifference$diff15.9 > pi, c(abs(angledifference$diff15.9 - (2*pi))), c(angledifference$diff15.9))
  angledifference$diff15.10 <- abs(fly15$oriradians - angles$ang15.10)
  angledifference$diff15.10 <- ifelse(angledifference$diff15.10 > pi, c(abs(angledifference$diff15.10 - (2*pi))), c(angledifference$diff15.10))
  angledifference$diff15.11 <- abs(fly15$oriradians - angles$ang15.11)
  angledifference$diff15.11 <- ifelse(angledifference$diff15.11 > pi, c(abs(angledifference$diff15.11 - (2*pi))), c(angledifference$diff15.11))
  angledifference$diff15.12 <- abs(fly15$oriradians - angles$ang15.12)
  angledifference$diff15.12 <- ifelse(angledifference$diff15.12 > pi, c(abs(angledifference$diff15.12 - (2*pi))), c(angledifference$diff15.12))
  angledifference$diff15.13 <- abs(fly15$oriradians - angles$ang15.13)
  angledifference$diff15.13 <- ifelse(angledifference$diff15.13 > pi, c(abs(angledifference$diff15.13 - (2*pi))), c(angledifference$diff15.13))
  angledifference$diff15.14 <- abs(fly15$oriradians - angles$ang15.14)
  angledifference$diff15.14 <- ifelse(angledifference$diff15.14 > pi, c(abs(angledifference$diff15.14 - (2*pi))), c(angledifference$diff15.14))
  angledifference$diff15.16 <- abs(fly15$oriradians - angles$ang15.16)
  angledifference$diff15.16 <- ifelse(angledifference$diff15.16 > pi, c(abs(angledifference$diff15.16 - (2*pi))), c(angledifference$diff15.16))
  angledifference$diff15.17 <- abs(fly15$oriradians - angles$ang15.17)
  angledifference$diff15.17 <- ifelse(angledifference$diff15.17 > pi, c(abs(angledifference$diff15.17 - (2*pi))), c(angledifference$diff15.17))
  angledifference$diff15.18 <- abs(fly15$oriradians - angles$ang15.18)
  angledifference$diff15.18 <- ifelse(angledifference$diff15.18 > pi, c(abs(angledifference$diff15.18 - (2*pi))), c(angledifference$diff15.18))
  angledifference$diff15.19 <- abs(fly15$oriradians - angles$ang15.19)
  angledifference$diff15.19 <- ifelse(angledifference$diff15.19 > pi, c(abs(angledifference$diff15.19 - (2*pi))), c(angledifference$diff15.19))
  angledifference$diff15.20 <- abs(fly15$oriradians - angles$ang15.20)
  angledifference$diff15.20 <- ifelse(angledifference$diff15.20 > pi, c(abs(angledifference$diff15.20 - (2*pi))), c(angledifference$diff15.20))
  angledifference$diff16.1 <- abs(fly16$oriradians - angles$ang16.1)
  angledifference$diff16.1 <- ifelse(angledifference$diff16.1 > pi, c(abs(angledifference$diff16.1 - (2*pi))), c(angledifference$diff16.1))
  angledifference$diff16.2 <- abs(fly16$oriradians - angles$ang16.2)
  angledifference$diff16.2 <- ifelse(angledifference$diff16.2 > pi, c(abs(angledifference$diff16.2 - (2*pi))), c(angledifference$diff16.2))
  angledifference$diff16.3 <- abs(fly16$oriradians - angles$ang16.3)
  angledifference$diff16.3 <- ifelse(angledifference$diff16.3 > pi, c(abs(angledifference$diff16.3 - (2*pi))), c(angledifference$diff16.3))
  angledifference$diff16.4 <- abs(fly16$oriradians - angles$ang16.4)
  angledifference$diff16.4 <- ifelse(angledifference$diff16.4 > pi, c(abs(angledifference$diff16.4 - (2*pi))), c(angledifference$diff16.4))
  angledifference$diff16.5 <- abs(fly16$oriradians - angles$ang16.5)
  angledifference$diff16.5 <- ifelse(angledifference$diff16.5 > pi, c(abs(angledifference$diff16.5 - (2*pi))), c(angledifference$diff16.5))
  angledifference$diff16.6 <- abs(fly16$oriradians - angles$ang16.6)
  angledifference$diff16.6 <- ifelse(angledifference$diff16.6 > pi, c(abs(angledifference$diff16.6 - (2*pi))), c(angledifference$diff16.6))
  angledifference$diff16.7 <- abs(fly16$oriradians - angles$ang16.7)
  angledifference$diff16.7 <- ifelse(angledifference$diff16.7 > pi, c(abs(angledifference$diff16.7 - (2*pi))), c(angledifference$diff16.7))
  angledifference$diff16.8 <- abs(fly16$oriradians - angles$ang16.8)
  angledifference$diff16.8 <- ifelse(angledifference$diff16.8 > pi, c(abs(angledifference$diff16.8 - (2*pi))), c(angledifference$diff16.8))
  angledifference$diff16.9 <- abs(fly16$oriradians - angles$ang16.9)
  angledifference$diff16.9 <- ifelse(angledifference$diff16.9 > pi, c(abs(angledifference$diff16.9 - (2*pi))), c(angledifference$diff16.9))
  angledifference$diff16.10 <- abs(fly16$oriradians - angles$ang16.10)
  angledifference$diff16.10 <- ifelse(angledifference$diff16.10 > pi, c(abs(angledifference$diff16.10 - (2*pi))), c(angledifference$diff16.10))
  angledifference$diff16.11 <- abs(fly16$oriradians - angles$ang16.11)
  angledifference$diff16.11 <- ifelse(angledifference$diff16.11 > pi, c(abs(angledifference$diff16.11 - (2*pi))), c(angledifference$diff16.11))
  angledifference$diff16.12 <- abs(fly16$oriradians - angles$ang16.12)
  angledifference$diff16.12 <- ifelse(angledifference$diff16.12 > pi, c(abs(angledifference$diff16.12 - (2*pi))), c(angledifference$diff16.12))
  angledifference$diff16.13 <- abs(fly16$oriradians - angles$ang16.13)
  angledifference$diff16.13 <- ifelse(angledifference$diff16.13 > pi, c(abs(angledifference$diff16.13 - (2*pi))), c(angledifference$diff16.13))
  angledifference$diff16.14 <- abs(fly16$oriradians - angles$ang16.14)
  angledifference$diff16.14 <- ifelse(angledifference$diff16.14 > pi, c(abs(angledifference$diff16.14 - (2*pi))), c(angledifference$diff16.14))
  angledifference$diff16.15 <- abs(fly16$oriradians - angles$ang16.15)
  angledifference$diff16.15 <- ifelse(angledifference$diff16.15 > pi, c(abs(angledifference$diff16.15 - (2*pi))), c(angledifference$diff16.15))
  angledifference$diff16.17 <- abs(fly16$oriradians - angles$ang16.17)
  angledifference$diff16.17 <- ifelse(angledifference$diff16.17 > pi, c(abs(angledifference$diff16.17 - (2*pi))), c(angledifference$diff16.17))
  angledifference$diff16.18 <- abs(fly16$oriradians - angles$ang16.18)
  angledifference$diff16.18 <- ifelse(angledifference$diff16.18 > pi, c(abs(angledifference$diff16.18 - (2*pi))), c(angledifference$diff16.18))
  angledifference$diff16.19 <- abs(fly16$oriradians - angles$ang16.19)
  angledifference$diff16.19 <- ifelse(angledifference$diff16.19 > pi, c(abs(angledifference$diff16.19 - (2*pi))), c(angledifference$diff16.19))
  angledifference$diff16.20 <- abs(fly16$oriradians - angles$ang16.20)
  angledifference$diff16.20 <- ifelse(angledifference$diff16.20 > pi, c(abs(angledifference$diff16.20 - (2*pi))), c(angledifference$diff16.20))
  angledifference$diff17.1 <- abs(fly17$oriradians - angles$ang17.1)
  angledifference$diff17.1 <- ifelse(angledifference$diff17.1 > pi, c(abs(angledifference$diff17.1 - (2*pi))), c(angledifference$diff17.1))
  angledifference$diff17.2 <- abs(fly17$oriradians - angles$ang17.2)
  angledifference$diff17.2 <- ifelse(angledifference$diff17.2 > pi, c(abs(angledifference$diff17.2 - (2*pi))), c(angledifference$diff17.2))
  angledifference$diff17.3 <- abs(fly17$oriradians - angles$ang17.3)
  angledifference$diff17.3 <- ifelse(angledifference$diff17.3 > pi, c(abs(angledifference$diff17.3 - (2*pi))), c(angledifference$diff17.3))
  angledifference$diff17.4 <- abs(fly17$oriradians - angles$ang17.4)
  angledifference$diff17.4 <- ifelse(angledifference$diff17.4 > pi, c(abs(angledifference$diff17.4 - (2*pi))), c(angledifference$diff17.4))
  angledifference$diff17.5 <- abs(fly17$oriradians - angles$ang17.5)
  angledifference$diff17.5 <- ifelse(angledifference$diff17.5 > pi, c(abs(angledifference$diff17.5 - (2*pi))), c(angledifference$diff17.5))
  angledifference$diff17.6 <- abs(fly17$oriradians - angles$ang17.6)
  angledifference$diff17.6 <- ifelse(angledifference$diff17.6 > pi, c(abs(angledifference$diff17.6 - (2*pi))), c(angledifference$diff17.6))
  angledifference$diff17.7 <- abs(fly17$oriradians - angles$ang17.7)
  angledifference$diff17.7 <- ifelse(angledifference$diff17.7 > pi, c(abs(angledifference$diff17.7 - (2*pi))), c(angledifference$diff17.7))
  angledifference$diff17.8 <- abs(fly17$oriradians - angles$ang17.8)
  angledifference$diff17.8 <- ifelse(angledifference$diff17.8 > pi, c(abs(angledifference$diff17.8 - (2*pi))), c(angledifference$diff17.8))
  angledifference$diff17.9 <- abs(fly17$oriradians - angles$ang17.9)
  angledifference$diff17.9 <- ifelse(angledifference$diff17.9 > pi, c(abs(angledifference$diff17.9 - (2*pi))), c(angledifference$diff17.9))
  angledifference$diff17.10 <- abs(fly17$oriradians - angles$ang17.10)
  angledifference$diff17.10 <- ifelse(angledifference$diff17.10 > pi, c(abs(angledifference$diff17.10 - (2*pi))), c(angledifference$diff17.10))
  angledifference$diff17.11 <- abs(fly17$oriradians - angles$ang17.11)
  angledifference$diff17.11 <- ifelse(angledifference$diff17.11 > pi, c(abs(angledifference$diff17.11 - (2*pi))), c(angledifference$diff17.11))
  angledifference$diff17.12 <- abs(fly17$oriradians - angles$ang17.12)
  angledifference$diff17.12 <- ifelse(angledifference$diff17.12 > pi, c(abs(angledifference$diff17.12 - (2*pi))), c(angledifference$diff17.12))
  angledifference$diff17.13 <- abs(fly17$oriradians - angles$ang17.13)
  angledifference$diff17.13 <- ifelse(angledifference$diff17.13 > pi, c(abs(angledifference$diff17.13 - (2*pi))), c(angledifference$diff17.13))
  angledifference$diff17.14 <- abs(fly17$oriradians - angles$ang17.14)
  angledifference$diff17.14 <- ifelse(angledifference$diff17.14 > pi, c(abs(angledifference$diff17.14 - (2*pi))), c(angledifference$diff17.14))
  angledifference$diff17.15 <- abs(fly17$oriradians - angles$ang17.15)
  angledifference$diff17.15 <- ifelse(angledifference$diff17.15 > pi, c(abs(angledifference$diff17.15 - (2*pi))), c(angledifference$diff17.15))
  angledifference$diff17.16 <- abs(fly17$oriradians - angles$ang17.16)
  angledifference$diff17.16 <- ifelse(angledifference$diff17.16 > pi, c(abs(angledifference$diff17.16 - (2*pi))), c(angledifference$diff17.16))
  angledifference$diff17.18 <- abs(fly17$oriradians - angles$ang17.18)
  angledifference$diff17.18 <- ifelse(angledifference$diff17.18 > pi, c(abs(angledifference$diff17.18 - (2*pi))), c(angledifference$diff17.18))
  angledifference$diff17.19 <- abs(fly17$oriradians - angles$ang17.19)
  angledifference$diff17.19 <- ifelse(angledifference$diff17.19 > pi, c(abs(angledifference$diff17.19 - (2*pi))), c(angledifference$diff17.19))
  angledifference$diff17.20 <- abs(fly17$oriradians - angles$ang17.20)
  angledifference$diff17.20 <- ifelse(angledifference$diff17.20 > pi, c(abs(angledifference$diff17.20 - (2*pi))), c(angledifference$diff17.20))
  angledifference$diff18.1 <- abs(fly18$oriradians - angles$ang18.1)
  angledifference$diff18.1 <- ifelse(angledifference$diff18.1 > pi, c(abs(angledifference$diff18.1 - (2*pi))), c(angledifference$diff18.1))
  angledifference$diff18.2 <- abs(fly18$oriradians - angles$ang18.2)
  angledifference$diff18.2 <- ifelse(angledifference$diff18.2 > pi, c(abs(angledifference$diff18.2 - (2*pi))), c(angledifference$diff18.2))
  angledifference$diff18.3 <- abs(fly18$oriradians - angles$ang18.3)
  angledifference$diff18.3 <- ifelse(angledifference$diff18.3 > pi, c(abs(angledifference$diff18.3 - (2*pi))), c(angledifference$diff18.3))
  angledifference$diff18.4 <- abs(fly18$oriradians - angles$ang18.4)
  angledifference$diff18.4 <- ifelse(angledifference$diff18.4 > pi, c(abs(angledifference$diff18.4 - (2*pi))), c(angledifference$diff18.4))
  angledifference$diff18.5 <- abs(fly18$oriradians - angles$ang18.5)
  angledifference$diff18.5 <- ifelse(angledifference$diff18.5 > pi, c(abs(angledifference$diff18.5 - (2*pi))), c(angledifference$diff18.5))
  angledifference$diff18.6 <- abs(fly18$oriradians - angles$ang18.6)
  angledifference$diff18.6 <- ifelse(angledifference$diff18.6 > pi, c(abs(angledifference$diff18.6 - (2*pi))), c(angledifference$diff18.6))
  angledifference$diff18.7 <- abs(fly18$oriradians - angles$ang18.7)
  angledifference$diff18.7 <- ifelse(angledifference$diff18.7 > pi, c(abs(angledifference$diff18.7 - (2*pi))), c(angledifference$diff18.7))
  angledifference$diff18.8 <- abs(fly18$oriradians - angles$ang18.8)
  angledifference$diff18.8 <- ifelse(angledifference$diff18.8 > pi, c(abs(angledifference$diff18.8 - (2*pi))), c(angledifference$diff18.8))
  angledifference$diff18.9 <- abs(fly18$oriradians - angles$ang18.9)
  angledifference$diff18.9 <- ifelse(angledifference$diff18.9 > pi, c(abs(angledifference$diff18.9 - (2*pi))), c(angledifference$diff18.9))
  angledifference$diff18.10 <- abs(fly18$oriradians - angles$ang18.10)
  angledifference$diff18.10 <- ifelse(angledifference$diff18.10 > pi, c(abs(angledifference$diff18.10 - (2*pi))), c(angledifference$diff18.10))
  angledifference$diff18.11 <- abs(fly18$oriradians - angles$ang18.11)
  angledifference$diff18.11 <- ifelse(angledifference$diff18.11 > pi, c(abs(angledifference$diff18.11 - (2*pi))), c(angledifference$diff18.11))
  angledifference$diff18.12 <- abs(fly18$oriradians - angles$ang18.12)
  angledifference$diff18.12 <- ifelse(angledifference$diff18.12 > pi, c(abs(angledifference$diff18.12 - (2*pi))), c(angledifference$diff18.12))
  angledifference$diff18.13 <- abs(fly18$oriradians - angles$ang18.13)
  angledifference$diff18.13 <- ifelse(angledifference$diff18.13 > pi, c(abs(angledifference$diff18.13 - (2*pi))), c(angledifference$diff18.13))
  angledifference$diff18.14 <- abs(fly18$oriradians - angles$ang18.14)
  angledifference$diff18.14 <- ifelse(angledifference$diff18.14 > pi, c(abs(angledifference$diff18.14 - (2*pi))), c(angledifference$diff18.14))
  angledifference$diff18.15 <- abs(fly18$oriradians - angles$ang18.15)
  angledifference$diff18.15 <- ifelse(angledifference$diff18.15 > pi, c(abs(angledifference$diff18.15 - (2*pi))), c(angledifference$diff18.15))
  angledifference$diff18.16 <- abs(fly18$oriradians - angles$ang18.16)
  angledifference$diff18.16 <- ifelse(angledifference$diff18.16 > pi, c(abs(angledifference$diff18.16 - (2*pi))), c(angledifference$diff18.16))
  angledifference$diff18.17 <- abs(fly18$oriradians - angles$ang18.17)
  angledifference$diff18.17 <- ifelse(angledifference$diff18.17 > pi, c(abs(angledifference$diff18.17 - (2*pi))), c(angledifference$diff18.17))
  angledifference$diff18.19 <- abs(fly18$oriradians - angles$ang18.19)
  angledifference$diff18.19 <- ifelse(angledifference$diff18.19 > pi, c(abs(angledifference$diff18.19 - (2*pi))), c(angledifference$diff18.19))
  angledifference$diff18.20 <- abs(fly18$oriradians - angles$ang18.20)
  angledifference$diff18.20 <- ifelse(angledifference$diff18.20 > pi, c(abs(angledifference$diff18.20 - (2*pi))), c(angledifference$diff18.20))
  angledifference$diff19.1 <- abs(fly19$oriradians - angles$ang19.1)
  angledifference$diff19.1 <- ifelse(angledifference$diff19.1 > pi, c(abs(angledifference$diff19.1 - (2*pi))), c(angledifference$diff19.1))
  angledifference$diff19.2 <- abs(fly19$oriradians - angles$ang19.2)
  angledifference$diff19.2 <- ifelse(angledifference$diff19.2 > pi, c(abs(angledifference$diff19.2 - (2*pi))), c(angledifference$diff19.2))
  angledifference$diff19.3 <- abs(fly19$oriradians - angles$ang19.3)
  angledifference$diff19.3 <- ifelse(angledifference$diff19.3 > pi, c(abs(angledifference$diff19.3 - (2*pi))), c(angledifference$diff19.3))
  angledifference$diff19.4 <- abs(fly19$oriradians - angles$ang19.4)
  angledifference$diff19.4 <- ifelse(angledifference$diff19.4 > pi, c(abs(angledifference$diff19.4 - (2*pi))), c(angledifference$diff19.4))
  angledifference$diff19.5 <- abs(fly19$oriradians - angles$ang19.5)
  angledifference$diff19.5 <- ifelse(angledifference$diff19.5 > pi, c(abs(angledifference$diff19.5 - (2*pi))), c(angledifference$diff19.5))
  angledifference$diff19.6 <- abs(fly19$oriradians - angles$ang19.6)
  angledifference$diff19.6 <- ifelse(angledifference$diff19.6 > pi, c(abs(angledifference$diff19.6 - (2*pi))), c(angledifference$diff19.6))
  angledifference$diff19.7 <- abs(fly19$oriradians - angles$ang19.7)
  angledifference$diff19.7 <- ifelse(angledifference$diff19.7 > pi, c(abs(angledifference$diff19.7 - (2*pi))), c(angledifference$diff19.7))
  angledifference$diff19.8 <- abs(fly19$oriradians - angles$ang19.8)
  angledifference$diff19.8 <- ifelse(angledifference$diff19.8 > pi, c(abs(angledifference$diff19.8 - (2*pi))), c(angledifference$diff19.8))
  angledifference$diff19.9 <- abs(fly19$oriradians - angles$ang19.9)
  angledifference$diff19.9 <- ifelse(angledifference$diff19.9 > pi, c(abs(angledifference$diff19.9 - (2*pi))), c(angledifference$diff19.9))
  angledifference$diff19.10 <- abs(fly19$oriradians - angles$ang19.10)
  angledifference$diff19.10 <- ifelse(angledifference$diff19.10 > pi, c(abs(angledifference$diff19.10 - (2*pi))), c(angledifference$diff19.10))
  angledifference$diff19.11 <- abs(fly19$oriradians - angles$ang19.11)
  angledifference$diff19.11 <- ifelse(angledifference$diff19.11 > pi, c(abs(angledifference$diff19.11 - (2*pi))), c(angledifference$diff19.11))
  angledifference$diff19.12 <- abs(fly19$oriradians - angles$ang19.12)
  angledifference$diff19.12 <- ifelse(angledifference$diff19.12 > pi, c(abs(angledifference$diff19.12 - (2*pi))), c(angledifference$diff19.12))
  angledifference$diff19.13 <- abs(fly19$oriradians - angles$ang19.13)
  angledifference$diff19.13 <- ifelse(angledifference$diff19.13 > pi, c(abs(angledifference$diff19.13 - (2*pi))), c(angledifference$diff19.13))
  angledifference$diff19.14 <- abs(fly19$oriradians - angles$ang19.14)
  angledifference$diff19.14 <- ifelse(angledifference$diff19.14 > pi, c(abs(angledifference$diff19.14 - (2*pi))), c(angledifference$diff19.14))
  angledifference$diff19.15 <- abs(fly19$oriradians - angles$ang19.15)
  angledifference$diff19.15 <- ifelse(angledifference$diff19.15 > pi, c(abs(angledifference$diff19.15 - (2*pi))), c(angledifference$diff19.15))
  angledifference$diff19.16 <- abs(fly19$oriradians - angles$ang19.16)
  angledifference$diff19.16 <- ifelse(angledifference$diff19.16 > pi, c(abs(angledifference$diff19.16 - (2*pi))), c(angledifference$diff19.16))
  angledifference$diff19.17 <- abs(fly19$oriradians - angles$ang19.17)
  angledifference$diff19.17 <- ifelse(angledifference$diff19.17 > pi, c(abs(angledifference$diff19.17 - (2*pi))), c(angledifference$diff19.17))
  angledifference$diff19.18 <- abs(fly19$oriradians - angles$ang19.18)
  angledifference$diff19.18 <- ifelse(angledifference$diff19.18 > pi, c(abs(angledifference$diff19.18 - (2*pi))), c(angledifference$diff19.18))
  angledifference$diff19.20 <- abs(fly19$oriradians - angles$ang19.20)
  angledifference$diff19.20 <- ifelse(angledifference$diff19.20 > pi, c(abs(angledifference$diff19.20 - (2*pi))), c(angledifference$diff19.20))
  angledifference$diff20.1 <- abs(fly20$oriradians - angles$ang20.1)
  angledifference$diff20.1 <- ifelse(angledifference$diff20.1 > pi, c(abs(angledifference$diff20.1 - (2*pi))), c(angledifference$diff20.1))
  angledifference$diff20.2 <- abs(fly20$oriradians - angles$ang20.2)
  angledifference$diff20.2 <- ifelse(angledifference$diff20.2 > pi, c(abs(angledifference$diff20.2 - (2*pi))), c(angledifference$diff20.2))
  angledifference$diff20.3 <- abs(fly20$oriradians - angles$ang20.3)
  angledifference$diff20.3 <- ifelse(angledifference$diff20.3 > pi, c(abs(angledifference$diff20.3 - (2*pi))), c(angledifference$diff20.3))
  angledifference$diff20.4 <- abs(fly20$oriradians - angles$ang20.4)
  angledifference$diff20.4 <- ifelse(angledifference$diff20.4 > pi, c(abs(angledifference$diff20.4 - (2*pi))), c(angledifference$diff20.4))
  angledifference$diff20.5 <- abs(fly20$oriradians - angles$ang20.5)
  angledifference$diff20.5 <- ifelse(angledifference$diff20.5 > pi, c(abs(angledifference$diff20.5 - (2*pi))), c(angledifference$diff20.5))
  angledifference$diff20.6 <- abs(fly20$oriradians - angles$ang20.6)
  angledifference$diff20.6 <- ifelse(angledifference$diff20.6 > pi, c(abs(angledifference$diff20.6 - (2*pi))), c(angledifference$diff20.6))
  angledifference$diff20.7 <- abs(fly20$oriradians - angles$ang20.7)
  angledifference$diff20.7 <- ifelse(angledifference$diff20.7 > pi, c(abs(angledifference$diff20.7 - (2*pi))), c(angledifference$diff20.7))
  angledifference$diff20.8 <- abs(fly20$oriradians - angles$ang20.8)
  angledifference$diff20.8 <- ifelse(angledifference$diff20.8 > pi, c(abs(angledifference$diff20.8 - (2*pi))), c(angledifference$diff20.8))
  angledifference$diff20.9 <- abs(fly20$oriradians - angles$ang20.9)
  angledifference$diff20.9 <- ifelse(angledifference$diff20.9 > pi, c(abs(angledifference$diff20.9 - (2*pi))), c(angledifference$diff20.9))
  angledifference$diff20.10 <- abs(fly20$oriradians - angles$ang20.10)
  angledifference$diff20.10 <- ifelse(angledifference$diff20.10 > pi, c(abs(angledifference$diff20.10 - (2*pi))), c(angledifference$diff20.10))
  angledifference$diff20.11 <- abs(fly20$oriradians - angles$ang20.11)
  angledifference$diff20.11 <- ifelse(angledifference$diff20.11 > pi, c(abs(angledifference$diff20.11 - (2*pi))), c(angledifference$diff20.11))
  angledifference$diff20.12 <- abs(fly20$oriradians - angles$ang20.12)
  angledifference$diff20.12 <- ifelse(angledifference$diff20.12 > pi, c(abs(angledifference$diff20.12 - (2*pi))), c(angledifference$diff20.12))
  angledifference$diff20.13 <- abs(fly20$oriradians - angles$ang20.13)
  angledifference$diff20.13 <- ifelse(angledifference$diff20.13 > pi, c(abs(angledifference$diff20.13 - (2*pi))), c(angledifference$diff20.13))
  angledifference$diff20.14 <- abs(fly20$oriradians - angles$ang20.14)
  angledifference$diff20.14 <- ifelse(angledifference$diff20.14 > pi, c(abs(angledifference$diff20.14 - (2*pi))), c(angledifference$diff20.14))
  angledifference$diff20.15 <- abs(fly20$oriradians - angles$ang20.15)
  angledifference$diff20.15 <- ifelse(angledifference$diff20.15 > pi, c(abs(angledifference$diff20.15 - (2*pi))), c(angledifference$diff20.15))
  angledifference$diff20.16 <- abs(fly20$oriradians - angles$ang20.16)
  angledifference$diff20.16 <- ifelse(angledifference$diff20.16 > pi, c(abs(angledifference$diff20.16 - (2*pi))), c(angledifference$diff20.16))
  angledifference$diff20.17 <- abs(fly20$oriradians - angles$ang20.17)
  angledifference$diff20.17 <- ifelse(angledifference$diff20.17 > pi, c(abs(angledifference$diff20.17 - (2*pi))), c(angledifference$diff20.17))
  angledifference$diff20.18 <- abs(fly20$oriradians - angles$ang20.18)
  angledifference$diff20.18 <- ifelse(angledifference$diff20.18 > pi, c(abs(angledifference$diff20.18 - (2*pi))), c(angledifference$diff20.18))
  angledifference$diff20.19 <- abs(fly20$oriradians - angles$ang20.19)
  angledifference$diff20.19 <- ifelse(angledifference$diff20.19 > pi, c(abs(angledifference$diff20.19 - (2*pi))), c(angledifference$diff20.19))
  
  # CREATE A DATAFRAME WHERE IF DIFFERENCE BETWEEN FOCAL FLY'S ORIENTATION AND IT'S ANGLE TO OTHER FLY <  (8pi/9), THEN 1, IF NO = 0 # -------------
  angleinteraction = data.frame(ifelse(angledifference$diff1.2 < ((8*pi)/9), c(1), c(0)))
  colnames(angleinteraction) <- c("int1.2")
  angleinteraction$int1.3 <- ifelse(angledifference$diff1.3 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int1.4 <- ifelse(angledifference$diff1.4 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int1.5 <- ifelse(angledifference$diff1.5 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int1.6 <- ifelse(angledifference$diff1.6 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int1.7 <- ifelse(angledifference$diff1.7 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int1.8 <- ifelse(angledifference$diff1.8 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int1.9 <- ifelse(angledifference$diff1.9 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int1.10 <- ifelse(angledifference$diff1.10 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int1.11 <- ifelse(angledifference$diff1.11 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int1.12 <- ifelse(angledifference$diff1.12 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int1.13 <- ifelse(angledifference$diff1.13 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int1.14 <- ifelse(angledifference$diff1.14 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int1.15 <- ifelse(angledifference$diff1.15 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int1.16 <- ifelse(angledifference$diff1.16 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int1.17 <- ifelse(angledifference$diff1.17 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int1.18 <- ifelse(angledifference$diff1.18 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int1.19 <- ifelse(angledifference$diff1.19 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int1.20 <- ifelse(angledifference$diff1.20 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int2.1 <- ifelse(angledifference$diff2.1 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int2.3 <- ifelse(angledifference$diff2.3 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int2.4 <- ifelse(angledifference$diff2.4 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int2.5 <- ifelse(angledifference$diff2.5 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int2.6 <- ifelse(angledifference$diff2.6 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int2.7 <- ifelse(angledifference$diff2.7 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int2.8 <- ifelse(angledifference$diff2.8 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int2.9 <- ifelse(angledifference$diff2.9 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int2.10 <- ifelse(angledifference$diff2.10 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int2.11 <- ifelse(angledifference$diff2.11 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int2.12 <- ifelse(angledifference$diff2.12 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int2.13 <- ifelse(angledifference$diff2.13 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int2.14 <- ifelse(angledifference$diff2.14 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int2.15 <- ifelse(angledifference$diff2.15 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int2.16 <- ifelse(angledifference$diff2.16 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int2.17 <- ifelse(angledifference$diff2.17 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int2.18 <- ifelse(angledifference$diff2.18 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int2.19 <- ifelse(angledifference$diff2.19 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int2.20 <- ifelse(angledifference$diff2.20 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int3.1 <- ifelse(angledifference$diff3.1 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int3.2 <- ifelse(angledifference$diff3.2 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int3.4 <- ifelse(angledifference$diff3.4 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int3.5 <- ifelse(angledifference$diff3.5 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int3.6 <- ifelse(angledifference$diff3.6 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int3.7 <- ifelse(angledifference$diff3.7 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int3.8 <- ifelse(angledifference$diff3.8 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int3.9 <- ifelse(angledifference$diff3.9 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int3.10 <- ifelse(angledifference$diff3.10 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int3.11 <- ifelse(angledifference$diff3.11 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int3.12 <- ifelse(angledifference$diff3.12 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int3.13 <- ifelse(angledifference$diff3.13 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int3.14 <- ifelse(angledifference$diff3.14 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int3.15 <- ifelse(angledifference$diff3.15 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int3.16 <- ifelse(angledifference$diff3.16 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int3.17 <- ifelse(angledifference$diff3.17 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int3.18 <- ifelse(angledifference$diff3.18 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int3.19 <- ifelse(angledifference$diff3.19 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int3.20 <- ifelse(angledifference$diff3.20 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int4.1 <- ifelse(angledifference$diff4.1 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int4.2 <- ifelse(angledifference$diff4.2 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int4.3 <- ifelse(angledifference$diff4.3 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int4.5 <- ifelse(angledifference$diff4.5 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int4.6 <- ifelse(angledifference$diff4.6 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int4.7 <- ifelse(angledifference$diff4.7 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int4.8 <- ifelse(angledifference$diff4.8 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int4.9 <- ifelse(angledifference$diff4.9 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int4.10 <- ifelse(angledifference$diff4.10 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int4.11 <- ifelse(angledifference$diff4.11 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int4.12 <- ifelse(angledifference$diff4.12 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int4.13 <- ifelse(angledifference$diff4.13 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int4.14 <- ifelse(angledifference$diff4.14 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int4.15 <- ifelse(angledifference$diff4.15 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int4.16 <- ifelse(angledifference$diff4.16 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int4.17 <- ifelse(angledifference$diff4.17 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int4.18 <- ifelse(angledifference$diff4.18 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int4.19 <- ifelse(angledifference$diff4.19 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int4.20 <- ifelse(angledifference$diff4.20 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int5.1 <- ifelse(angledifference$diff5.1 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int5.2 <- ifelse(angledifference$diff5.2 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int5.3 <- ifelse(angledifference$diff5.3 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int5.4 <- ifelse(angledifference$diff5.4 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int5.6 <- ifelse(angledifference$diff5.6 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int5.7 <- ifelse(angledifference$diff5.7 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int5.8 <- ifelse(angledifference$diff5.8 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int5.9 <- ifelse(angledifference$diff5.9 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int5.10 <- ifelse(angledifference$diff5.10 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int5.11 <- ifelse(angledifference$diff5.11 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int5.12 <- ifelse(angledifference$diff5.12 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int5.13 <- ifelse(angledifference$diff5.13 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int5.14 <- ifelse(angledifference$diff5.14 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int5.15 <- ifelse(angledifference$diff5.15 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int5.16 <- ifelse(angledifference$diff5.16 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int5.17 <- ifelse(angledifference$diff5.17 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int5.18 <- ifelse(angledifference$diff5.18 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int5.19 <- ifelse(angledifference$diff5.19 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int5.20 <- ifelse(angledifference$diff5.20 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int6.1 <- ifelse(angledifference$diff6.1 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int6.2 <- ifelse(angledifference$diff6.2 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int6.3 <- ifelse(angledifference$diff6.3 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int6.4 <- ifelse(angledifference$diff6.4 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int6.5 <- ifelse(angledifference$diff6.5 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int6.7 <- ifelse(angledifference$diff6.7 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int6.8 <- ifelse(angledifference$diff6.8 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int6.9 <- ifelse(angledifference$diff6.9 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int6.10 <- ifelse(angledifference$diff6.10 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int6.11 <- ifelse(angledifference$diff6.11 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int6.12 <- ifelse(angledifference$diff6.12 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int6.13 <- ifelse(angledifference$diff6.13 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int6.14 <- ifelse(angledifference$diff6.14 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int6.15 <- ifelse(angledifference$diff6.15 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int6.16 <- ifelse(angledifference$diff6.16 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int6.17 <- ifelse(angledifference$diff6.17 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int6.18 <- ifelse(angledifference$diff6.18 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int6.19 <- ifelse(angledifference$diff6.19 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int6.20 <- ifelse(angledifference$diff6.20 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int7.1 <- ifelse(angledifference$diff7.1 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int7.2 <- ifelse(angledifference$diff7.2 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int7.3 <- ifelse(angledifference$diff7.3 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int7.4 <- ifelse(angledifference$diff7.4 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int7.5 <- ifelse(angledifference$diff7.5 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int7.6 <- ifelse(angledifference$diff7.6 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int7.8 <- ifelse(angledifference$diff7.8 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int7.9 <- ifelse(angledifference$diff7.9 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int7.10 <- ifelse(angledifference$diff7.10 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int7.11 <- ifelse(angledifference$diff7.11 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int7.12 <- ifelse(angledifference$diff7.12 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int7.13 <- ifelse(angledifference$diff7.13 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int7.14 <- ifelse(angledifference$diff7.14 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int7.15 <- ifelse(angledifference$diff7.15 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int7.16 <- ifelse(angledifference$diff7.16 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int7.17 <- ifelse(angledifference$diff7.17 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int7.18 <- ifelse(angledifference$diff7.18 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int7.19 <- ifelse(angledifference$diff7.19 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int7.20 <- ifelse(angledifference$diff7.20 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int8.1 <- ifelse(angledifference$diff8.1 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int8.2 <- ifelse(angledifference$diff8.2 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int8.3 <- ifelse(angledifference$diff8.3 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int8.4 <- ifelse(angledifference$diff8.4 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int8.5 <- ifelse(angledifference$diff8.5 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int8.6 <- ifelse(angledifference$diff8.6 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int8.7 <- ifelse(angledifference$diff8.7 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int8.9 <- ifelse(angledifference$diff8.9 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int8.10 <- ifelse(angledifference$diff8.10 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int8.11 <- ifelse(angledifference$diff8.11 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int8.12 <- ifelse(angledifference$diff8.12 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int8.13 <- ifelse(angledifference$diff8.13 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int8.14 <- ifelse(angledifference$diff8.14 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int8.15 <- ifelse(angledifference$diff8.15 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int8.16 <- ifelse(angledifference$diff8.16 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int8.17 <- ifelse(angledifference$diff8.17 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int8.18 <- ifelse(angledifference$diff8.18 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int8.19 <- ifelse(angledifference$diff8.19 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int8.20 <- ifelse(angledifference$diff8.20 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int9.1 <- ifelse(angledifference$diff9.1 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int9.2 <- ifelse(angledifference$diff9.2 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int9.3 <- ifelse(angledifference$diff9.3 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int9.4 <- ifelse(angledifference$diff9.4 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int9.5 <- ifelse(angledifference$diff9.5 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int9.6 <- ifelse(angledifference$diff9.6 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int9.7 <- ifelse(angledifference$diff9.7 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int9.8 <- ifelse(angledifference$diff9.8 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int9.10 <- ifelse(angledifference$diff9.10 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int9.11 <- ifelse(angledifference$diff9.11 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int9.12 <- ifelse(angledifference$diff9.12 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int9.13 <- ifelse(angledifference$diff9.13 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int9.14 <- ifelse(angledifference$diff9.14 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int9.15 <- ifelse(angledifference$diff9.15 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int9.16 <- ifelse(angledifference$diff9.16 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int9.17 <- ifelse(angledifference$diff9.17 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int9.18 <- ifelse(angledifference$diff9.18 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int9.19 <- ifelse(angledifference$diff9.19 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int9.20 <- ifelse(angledifference$diff9.20 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int10.1 <- ifelse(angledifference$diff10.1 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int10.2 <- ifelse(angledifference$diff10.2 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int10.3 <- ifelse(angledifference$diff10.3 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int10.4 <- ifelse(angledifference$diff10.4 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int10.5 <- ifelse(angledifference$diff10.5 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int10.6 <- ifelse(angledifference$diff10.6 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int10.7 <- ifelse(angledifference$diff10.7 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int10.8 <- ifelse(angledifference$diff10.8 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int10.9 <- ifelse(angledifference$diff10.9 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int10.11 <- ifelse(angledifference$diff10.11 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int10.12 <- ifelse(angledifference$diff10.12 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int10.13 <- ifelse(angledifference$diff10.13 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int10.14 <- ifelse(angledifference$diff10.14 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int10.15 <- ifelse(angledifference$diff10.15 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int10.16 <- ifelse(angledifference$diff10.16 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int10.17 <- ifelse(angledifference$diff10.17 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int10.18 <- ifelse(angledifference$diff10.18 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int10.19 <- ifelse(angledifference$diff10.19 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int10.20 <- ifelse(angledifference$diff10.20 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int11.1 <- ifelse(angledifference$diff11.1 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int11.2 <- ifelse(angledifference$diff11.2 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int11.3 <- ifelse(angledifference$diff11.3 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int11.4 <- ifelse(angledifference$diff11.4 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int11.5 <- ifelse(angledifference$diff11.5 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int11.6 <- ifelse(angledifference$diff11.6 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int11.7 <- ifelse(angledifference$diff11.7 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int11.8 <- ifelse(angledifference$diff11.8 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int11.9 <- ifelse(angledifference$diff11.9 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int11.10 <- ifelse(angledifference$diff11.10 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int11.12 <- ifelse(angledifference$diff11.12 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int11.13 <- ifelse(angledifference$diff11.13 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int11.14 <- ifelse(angledifference$diff11.14 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int11.15 <- ifelse(angledifference$diff11.15 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int11.16 <- ifelse(angledifference$diff11.16 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int11.17 <- ifelse(angledifference$diff11.17 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int11.18 <- ifelse(angledifference$diff11.18 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int11.19 <- ifelse(angledifference$diff11.19 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int11.20 <- ifelse(angledifference$diff11.20 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int12.1 <- ifelse(angledifference$diff12.1 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int12.2 <- ifelse(angledifference$diff12.2 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int12.3 <- ifelse(angledifference$diff12.3 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int12.4 <- ifelse(angledifference$diff12.4 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int12.5 <- ifelse(angledifference$diff12.5 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int12.6 <- ifelse(angledifference$diff12.6 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int12.7 <- ifelse(angledifference$diff12.7 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int12.8 <- ifelse(angledifference$diff12.8 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int12.9 <- ifelse(angledifference$diff12.9 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int12.10 <- ifelse(angledifference$diff12.10 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int12.11 <- ifelse(angledifference$diff12.11 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int12.13 <- ifelse(angledifference$diff12.13 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int12.14 <- ifelse(angledifference$diff12.14 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int12.15 <- ifelse(angledifference$diff12.15 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int12.16 <- ifelse(angledifference$diff12.16 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int12.17 <- ifelse(angledifference$diff12.17 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int12.18 <- ifelse(angledifference$diff12.18 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int12.19 <- ifelse(angledifference$diff12.19 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int12.20 <- ifelse(angledifference$diff12.20 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int13.1 <- ifelse(angledifference$diff13.1 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int13.2 <- ifelse(angledifference$diff13.2 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int13.3 <- ifelse(angledifference$diff13.3 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int13.4 <- ifelse(angledifference$diff13.4 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int13.5 <- ifelse(angledifference$diff13.5 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int13.6 <- ifelse(angledifference$diff13.6 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int13.7 <- ifelse(angledifference$diff13.7 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int13.8 <- ifelse(angledifference$diff13.8 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int13.9 <- ifelse(angledifference$diff13.9 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int13.10 <- ifelse(angledifference$diff13.10 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int13.11 <- ifelse(angledifference$diff13.11 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int13.12 <- ifelse(angledifference$diff13.12 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int13.14 <- ifelse(angledifference$diff13.14 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int13.15 <- ifelse(angledifference$diff13.15 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int13.16 <- ifelse(angledifference$diff13.16 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int13.17 <- ifelse(angledifference$diff13.17 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int13.18 <- ifelse(angledifference$diff13.18 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int13.19 <- ifelse(angledifference$diff13.19 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int13.20 <- ifelse(angledifference$diff13.20 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int14.1 <- ifelse(angledifference$diff14.1 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int14.2 <- ifelse(angledifference$diff14.2 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int14.3 <- ifelse(angledifference$diff14.3 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int14.4 <- ifelse(angledifference$diff14.4 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int14.5 <- ifelse(angledifference$diff14.5 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int14.6 <- ifelse(angledifference$diff14.6 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int14.7 <- ifelse(angledifference$diff14.7 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int14.8 <- ifelse(angledifference$diff14.8 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int14.9 <- ifelse(angledifference$diff14.9 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int14.10 <- ifelse(angledifference$diff14.10 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int14.11 <- ifelse(angledifference$diff14.11 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int14.12 <- ifelse(angledifference$diff14.12 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int14.13 <- ifelse(angledifference$diff14.13 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int14.15 <- ifelse(angledifference$diff14.15 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int14.16 <- ifelse(angledifference$diff14.16 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int14.17 <- ifelse(angledifference$diff14.17 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int14.18 <- ifelse(angledifference$diff14.18 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int14.19 <- ifelse(angledifference$diff14.19 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int14.20 <- ifelse(angledifference$diff14.20 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int15.1 <- ifelse(angledifference$diff15.1 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int15.2 <- ifelse(angledifference$diff15.2 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int15.3 <- ifelse(angledifference$diff15.3 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int15.4 <- ifelse(angledifference$diff15.4 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int15.5 <- ifelse(angledifference$diff15.5 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int15.6 <- ifelse(angledifference$diff15.6 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int15.7 <- ifelse(angledifference$diff15.7 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int15.8 <- ifelse(angledifference$diff15.8 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int15.9 <- ifelse(angledifference$diff15.9 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int15.10 <- ifelse(angledifference$diff15.10 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int15.11 <- ifelse(angledifference$diff15.11 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int15.12 <- ifelse(angledifference$diff15.12 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int15.13 <- ifelse(angledifference$diff15.13 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int15.14 <- ifelse(angledifference$diff15.14 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int15.16 <- ifelse(angledifference$diff15.16 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int15.17 <- ifelse(angledifference$diff15.17 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int15.18 <- ifelse(angledifference$diff15.18 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int15.19 <- ifelse(angledifference$diff15.19 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int15.20 <- ifelse(angledifference$diff15.20 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int16.1 <- ifelse(angledifference$diff16.1 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int16.2 <- ifelse(angledifference$diff16.2 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int16.3 <- ifelse(angledifference$diff16.3 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int16.4 <- ifelse(angledifference$diff16.4 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int16.5 <- ifelse(angledifference$diff16.5 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int16.6 <- ifelse(angledifference$diff16.6 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int16.7 <- ifelse(angledifference$diff16.7 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int16.8 <- ifelse(angledifference$diff16.8 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int16.9 <- ifelse(angledifference$diff16.9 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int16.10 <- ifelse(angledifference$diff16.10 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int16.11 <- ifelse(angledifference$diff16.11 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int16.12 <- ifelse(angledifference$diff16.12 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int16.13 <- ifelse(angledifference$diff16.13 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int16.14 <- ifelse(angledifference$diff16.14 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int16.15 <- ifelse(angledifference$diff16.15 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int16.17 <- ifelse(angledifference$diff16.17 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int16.18 <- ifelse(angledifference$diff16.18 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int16.19 <- ifelse(angledifference$diff16.19 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int16.20 <- ifelse(angledifference$diff16.20 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int17.1 <- ifelse(angledifference$diff17.1 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int17.2 <- ifelse(angledifference$diff17.2 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int17.3 <- ifelse(angledifference$diff17.3 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int17.4 <- ifelse(angledifference$diff17.4 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int17.5 <- ifelse(angledifference$diff17.5 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int17.6 <- ifelse(angledifference$diff17.6 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int17.7 <- ifelse(angledifference$diff17.7 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int17.8 <- ifelse(angledifference$diff17.8 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int17.9 <- ifelse(angledifference$diff17.9 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int17.10 <- ifelse(angledifference$diff17.10 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int17.11 <- ifelse(angledifference$diff17.11 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int17.12 <- ifelse(angledifference$diff17.12 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int17.13 <- ifelse(angledifference$diff17.13 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int17.14 <- ifelse(angledifference$diff17.14 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int17.15 <- ifelse(angledifference$diff17.15 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int17.16 <- ifelse(angledifference$diff17.16 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int17.18 <- ifelse(angledifference$diff17.18 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int17.19 <- ifelse(angledifference$diff17.19 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int17.20 <- ifelse(angledifference$diff17.20 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int18.1 <- ifelse(angledifference$diff18.1 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int18.2 <- ifelse(angledifference$diff18.2 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int18.3 <- ifelse(angledifference$diff18.3 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int18.4 <- ifelse(angledifference$diff18.4 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int18.5 <- ifelse(angledifference$diff18.5 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int18.6 <- ifelse(angledifference$diff18.6 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int18.7 <- ifelse(angledifference$diff18.7 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int18.8 <- ifelse(angledifference$diff18.8 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int18.9 <- ifelse(angledifference$diff18.9 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int18.10 <- ifelse(angledifference$diff18.10 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int18.11 <- ifelse(angledifference$diff18.11 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int18.12 <- ifelse(angledifference$diff18.12 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int18.13 <- ifelse(angledifference$diff18.13 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int18.14 <- ifelse(angledifference$diff18.14 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int18.15 <- ifelse(angledifference$diff18.15 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int18.16 <- ifelse(angledifference$diff18.16 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int18.17 <- ifelse(angledifference$diff18.17 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int18.19 <- ifelse(angledifference$diff18.19 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int18.20 <- ifelse(angledifference$diff18.20 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int19.1 <- ifelse(angledifference$diff19.1 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int19.2 <- ifelse(angledifference$diff19.2 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int19.3 <- ifelse(angledifference$diff19.3 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int19.4 <- ifelse(angledifference$diff19.4 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int19.5 <- ifelse(angledifference$diff19.5 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int19.6 <- ifelse(angledifference$diff19.6 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int19.7 <- ifelse(angledifference$diff19.7 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int19.8 <- ifelse(angledifference$diff19.8 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int19.9 <- ifelse(angledifference$diff19.9 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int19.10 <- ifelse(angledifference$diff19.10 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int19.11 <- ifelse(angledifference$diff19.11 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int19.12 <- ifelse(angledifference$diff19.12 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int19.13 <- ifelse(angledifference$diff19.13 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int19.14 <- ifelse(angledifference$diff19.14 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int19.15 <- ifelse(angledifference$diff19.15 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int19.16 <- ifelse(angledifference$diff19.16 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int19.17 <- ifelse(angledifference$diff19.17 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int19.18 <- ifelse(angledifference$diff19.18 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int19.20 <- ifelse(angledifference$diff19.20 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int20.1 <- ifelse(angledifference$diff20.1 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int20.2 <- ifelse(angledifference$diff20.2 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int20.3 <- ifelse(angledifference$diff20.3 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int20.4 <- ifelse(angledifference$diff20.4 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int20.5 <- ifelse(angledifference$diff20.5 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int20.6 <- ifelse(angledifference$diff20.6 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int20.7 <- ifelse(angledifference$diff20.7 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int20.8 <- ifelse(angledifference$diff20.8 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int20.9 <- ifelse(angledifference$diff20.9 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int20.10 <- ifelse(angledifference$diff20.10 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int20.11 <- ifelse(angledifference$diff20.11 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int20.12 <- ifelse(angledifference$diff20.12 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int20.13 <- ifelse(angledifference$diff20.13 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int20.14 <- ifelse(angledifference$diff20.14 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int20.15 <- ifelse(angledifference$diff20.15 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int20.16 <- ifelse(angledifference$diff20.16 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int20.17 <- ifelse(angledifference$diff20.17 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int20.18 <- ifelse(angledifference$diff20.18 < ((8*pi)/9), c(1), c(0))
  angleinteraction$int20.19 <- ifelse(angledifference$diff20.19 < ((8*pi)/9), c(1), c(0))
  
  # CREATE A DATAFRAME WHERE IF DISTANCE AND ANGLE INTERACTION = 1, THEN 1, IF NO = 0 # ------------------------------------------
  distanceangleinteraction = data.frame(ifelse(distanceinteraction$int1.2 == 1, c(ifelse(angleinteraction$int1.2 == 1, c(1), c(0))), c(0)))
  colnames(distanceangleinteraction) <- c("int1.2")
  distanceangleinteraction$int1.3 <- ifelse(distanceinteraction$int1.3 == 1, c(ifelse(angleinteraction$int1.3 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int1.4 <- ifelse(distanceinteraction$int1.4 == 1, c(ifelse(angleinteraction$int1.4 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int1.5 <- ifelse(distanceinteraction$int1.5 == 1, c(ifelse(angleinteraction$int1.5 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int1.6 <- ifelse(distanceinteraction$int1.6 == 1, c(ifelse(angleinteraction$int1.6 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int1.7 <- ifelse(distanceinteraction$int1.7 == 1, c(ifelse(angleinteraction$int1.7 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int1.8 <- ifelse(distanceinteraction$int1.8 == 1, c(ifelse(angleinteraction$int1.8 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int1.9 <- ifelse(distanceinteraction$int1.9 == 1, c(ifelse(angleinteraction$int1.9 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int1.10 <- ifelse(distanceinteraction$int1.10 == 1, c(ifelse(angleinteraction$int1.10 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int1.11 <- ifelse(distanceinteraction$int1.11 == 1, c(ifelse(angleinteraction$int1.11 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int1.12 <- ifelse(distanceinteraction$int1.12 == 1, c(ifelse(angleinteraction$int1.12 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int1.13 <- ifelse(distanceinteraction$int1.13 == 1, c(ifelse(angleinteraction$int1.13 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int1.14 <- ifelse(distanceinteraction$int1.14 == 1, c(ifelse(angleinteraction$int1.14 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int1.15 <- ifelse(distanceinteraction$int1.15 == 1, c(ifelse(angleinteraction$int1.15 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int1.16 <- ifelse(distanceinteraction$int1.16 == 1, c(ifelse(angleinteraction$int1.16 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int1.17 <- ifelse(distanceinteraction$int1.17 == 1, c(ifelse(angleinteraction$int1.17 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int1.18 <- ifelse(distanceinteraction$int1.18 == 1, c(ifelse(angleinteraction$int1.18 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int1.19 <- ifelse(distanceinteraction$int1.19 == 1, c(ifelse(angleinteraction$int1.19 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int1.20 <- ifelse(distanceinteraction$int1.20 == 1, c(ifelse(angleinteraction$int1.20 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int2.1 <- ifelse(distanceinteraction$int2.1 == 1, c(ifelse(angleinteraction$int2.1 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int2.3 <- ifelse(distanceinteraction$int2.3 == 1, c(ifelse(angleinteraction$int2.3 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int2.4 <- ifelse(distanceinteraction$int2.4 == 1, c(ifelse(angleinteraction$int2.4 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int2.5 <- ifelse(distanceinteraction$int2.5 == 1, c(ifelse(angleinteraction$int2.5 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int2.6 <- ifelse(distanceinteraction$int2.6 == 1, c(ifelse(angleinteraction$int2.6 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int2.7 <- ifelse(distanceinteraction$int2.7 == 1, c(ifelse(angleinteraction$int2.7 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int2.8 <- ifelse(distanceinteraction$int2.8 == 1, c(ifelse(angleinteraction$int2.8 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int2.9 <- ifelse(distanceinteraction$int2.9 == 1, c(ifelse(angleinteraction$int2.9 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int2.10 <- ifelse(distanceinteraction$int2.10 == 1, c(ifelse(angleinteraction$int2.10 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int2.11 <- ifelse(distanceinteraction$int2.11 == 1, c(ifelse(angleinteraction$int2.11 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int2.12 <- ifelse(distanceinteraction$int2.12 == 1, c(ifelse(angleinteraction$int2.12 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int2.13 <- ifelse(distanceinteraction$int2.13 == 1, c(ifelse(angleinteraction$int2.13 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int2.14 <- ifelse(distanceinteraction$int2.14 == 1, c(ifelse(angleinteraction$int2.14 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int2.15 <- ifelse(distanceinteraction$int2.15 == 1, c(ifelse(angleinteraction$int2.15 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int2.16 <- ifelse(distanceinteraction$int2.16 == 1, c(ifelse(angleinteraction$int2.16 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int2.17 <- ifelse(distanceinteraction$int2.17 == 1, c(ifelse(angleinteraction$int2.17 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int2.18 <- ifelse(distanceinteraction$int2.18 == 1, c(ifelse(angleinteraction$int2.18 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int2.19 <- ifelse(distanceinteraction$int2.19 == 1, c(ifelse(angleinteraction$int2.19 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int2.20 <- ifelse(distanceinteraction$int2.20 == 1, c(ifelse(angleinteraction$int2.20 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int3.1 <- ifelse(distanceinteraction$int3.1 == 1, c(ifelse(angleinteraction$int3.1 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int3.2 <- ifelse(distanceinteraction$int3.2 == 1, c(ifelse(angleinteraction$int3.2 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int3.4 <- ifelse(distanceinteraction$int3.4 == 1, c(ifelse(angleinteraction$int3.4 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int3.5 <- ifelse(distanceinteraction$int3.5 == 1, c(ifelse(angleinteraction$int3.5 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int3.6 <- ifelse(distanceinteraction$int3.6 == 1, c(ifelse(angleinteraction$int3.6 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int3.7 <- ifelse(distanceinteraction$int3.7 == 1, c(ifelse(angleinteraction$int3.7 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int3.8 <- ifelse(distanceinteraction$int3.8 == 1, c(ifelse(angleinteraction$int3.8 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int3.9 <- ifelse(distanceinteraction$int3.9 == 1, c(ifelse(angleinteraction$int3.9 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int3.10 <- ifelse(distanceinteraction$int3.10 == 1, c(ifelse(angleinteraction$int3.10 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int3.11 <- ifelse(distanceinteraction$int3.11 == 1, c(ifelse(angleinteraction$int3.11 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int3.12 <- ifelse(distanceinteraction$int3.12 == 1, c(ifelse(angleinteraction$int3.12 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int3.13 <- ifelse(distanceinteraction$int3.13 == 1, c(ifelse(angleinteraction$int3.13 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int3.14 <- ifelse(distanceinteraction$int3.14 == 1, c(ifelse(angleinteraction$int3.14 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int3.15 <- ifelse(distanceinteraction$int3.15 == 1, c(ifelse(angleinteraction$int3.15 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int3.16 <- ifelse(distanceinteraction$int3.16 == 1, c(ifelse(angleinteraction$int3.16 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int3.17 <- ifelse(distanceinteraction$int3.17 == 1, c(ifelse(angleinteraction$int3.17 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int3.18 <- ifelse(distanceinteraction$int3.18 == 1, c(ifelse(angleinteraction$int3.18 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int3.19 <- ifelse(distanceinteraction$int3.19 == 1, c(ifelse(angleinteraction$int3.19 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int3.20 <- ifelse(distanceinteraction$int3.20 == 1, c(ifelse(angleinteraction$int3.20 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int4.1 <- ifelse(distanceinteraction$int4.1 == 1, c(ifelse(angleinteraction$int4.1 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int4.2 <- ifelse(distanceinteraction$int4.2 == 1, c(ifelse(angleinteraction$int4.2 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int4.3 <- ifelse(distanceinteraction$int4.3 == 1, c(ifelse(angleinteraction$int4.3 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int4.5 <- ifelse(distanceinteraction$int4.5 == 1, c(ifelse(angleinteraction$int4.5 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int4.6 <- ifelse(distanceinteraction$int4.6 == 1, c(ifelse(angleinteraction$int4.6 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int4.7 <- ifelse(distanceinteraction$int4.7 == 1, c(ifelse(angleinteraction$int4.7 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int4.8 <- ifelse(distanceinteraction$int4.8 == 1, c(ifelse(angleinteraction$int4.8 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int4.9 <- ifelse(distanceinteraction$int4.9 == 1, c(ifelse(angleinteraction$int4.9 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int4.10 <- ifelse(distanceinteraction$int4.10 == 1, c(ifelse(angleinteraction$int4.10 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int4.11 <- ifelse(distanceinteraction$int4.11 == 1, c(ifelse(angleinteraction$int4.11 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int4.12 <- ifelse(distanceinteraction$int4.12 == 1, c(ifelse(angleinteraction$int4.12 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int4.13 <- ifelse(distanceinteraction$int4.13 == 1, c(ifelse(angleinteraction$int4.13 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int4.14 <- ifelse(distanceinteraction$int4.14 == 1, c(ifelse(angleinteraction$int4.14 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int4.15 <- ifelse(distanceinteraction$int4.15 == 1, c(ifelse(angleinteraction$int4.15 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int4.16 <- ifelse(distanceinteraction$int4.16 == 1, c(ifelse(angleinteraction$int4.16 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int4.17 <- ifelse(distanceinteraction$int4.17 == 1, c(ifelse(angleinteraction$int4.17 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int4.18 <- ifelse(distanceinteraction$int4.18 == 1, c(ifelse(angleinteraction$int4.18 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int4.19 <- ifelse(distanceinteraction$int4.19 == 1, c(ifelse(angleinteraction$int4.19 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int4.20 <- ifelse(distanceinteraction$int4.20 == 1, c(ifelse(angleinteraction$int4.20 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int5.1 <- ifelse(distanceinteraction$int5.1 == 1, c(ifelse(angleinteraction$int5.1 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int5.2 <- ifelse(distanceinteraction$int5.2 == 1, c(ifelse(angleinteraction$int5.2 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int5.3 <- ifelse(distanceinteraction$int5.3 == 1, c(ifelse(angleinteraction$int5.3 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int5.4 <- ifelse(distanceinteraction$int5.4 == 1, c(ifelse(angleinteraction$int5.4 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int5.6 <- ifelse(distanceinteraction$int5.6 == 1, c(ifelse(angleinteraction$int5.6 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int5.7 <- ifelse(distanceinteraction$int5.7 == 1, c(ifelse(angleinteraction$int5.7 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int5.8 <- ifelse(distanceinteraction$int5.8 == 1, c(ifelse(angleinteraction$int5.8 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int5.9 <- ifelse(distanceinteraction$int5.9 == 1, c(ifelse(angleinteraction$int5.9 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int5.10 <- ifelse(distanceinteraction$int5.10 == 1, c(ifelse(angleinteraction$int5.10 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int5.11 <- ifelse(distanceinteraction$int5.11 == 1, c(ifelse(angleinteraction$int5.11 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int5.12 <- ifelse(distanceinteraction$int5.12 == 1, c(ifelse(angleinteraction$int5.12 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int5.13 <- ifelse(distanceinteraction$int5.13 == 1, c(ifelse(angleinteraction$int5.13 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int5.14 <- ifelse(distanceinteraction$int5.14 == 1, c(ifelse(angleinteraction$int5.14 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int5.15 <- ifelse(distanceinteraction$int5.15 == 1, c(ifelse(angleinteraction$int5.15 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int5.16 <- ifelse(distanceinteraction$int5.16 == 1, c(ifelse(angleinteraction$int5.16 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int5.17 <- ifelse(distanceinteraction$int5.17 == 1, c(ifelse(angleinteraction$int5.17 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int5.18 <- ifelse(distanceinteraction$int5.18 == 1, c(ifelse(angleinteraction$int5.18 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int5.19 <- ifelse(distanceinteraction$int5.19 == 1, c(ifelse(angleinteraction$int5.19 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int5.20 <- ifelse(distanceinteraction$int5.20 == 1, c(ifelse(angleinteraction$int5.20 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int6.1 <- ifelse(distanceinteraction$int6.1 == 1, c(ifelse(angleinteraction$int6.1 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int6.2 <- ifelse(distanceinteraction$int6.2 == 1, c(ifelse(angleinteraction$int6.2 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int6.3 <- ifelse(distanceinteraction$int6.3 == 1, c(ifelse(angleinteraction$int6.3 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int6.4 <- ifelse(distanceinteraction$int6.4 == 1, c(ifelse(angleinteraction$int6.4 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int6.5 <- ifelse(distanceinteraction$int6.5 == 1, c(ifelse(angleinteraction$int6.5 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int6.7 <- ifelse(distanceinteraction$int6.7 == 1, c(ifelse(angleinteraction$int6.7 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int6.8 <- ifelse(distanceinteraction$int6.8 == 1, c(ifelse(angleinteraction$int6.8 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int6.9 <- ifelse(distanceinteraction$int6.9 == 1, c(ifelse(angleinteraction$int6.9 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int6.10 <- ifelse(distanceinteraction$int6.10 == 1, c(ifelse(angleinteraction$int6.10 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int6.11 <- ifelse(distanceinteraction$int6.11 == 1, c(ifelse(angleinteraction$int6.11 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int6.12 <- ifelse(distanceinteraction$int6.12 == 1, c(ifelse(angleinteraction$int6.12 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int6.13 <- ifelse(distanceinteraction$int6.13 == 1, c(ifelse(angleinteraction$int6.13 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int6.14 <- ifelse(distanceinteraction$int6.14 == 1, c(ifelse(angleinteraction$int6.14 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int6.15 <- ifelse(distanceinteraction$int6.15 == 1, c(ifelse(angleinteraction$int6.15 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int6.16 <- ifelse(distanceinteraction$int6.16 == 1, c(ifelse(angleinteraction$int6.16 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int6.17 <- ifelse(distanceinteraction$int6.17 == 1, c(ifelse(angleinteraction$int6.17 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int6.18 <- ifelse(distanceinteraction$int6.18 == 1, c(ifelse(angleinteraction$int6.18 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int6.19 <- ifelse(distanceinteraction$int6.19 == 1, c(ifelse(angleinteraction$int6.19 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int6.20 <- ifelse(distanceinteraction$int6.20 == 1, c(ifelse(angleinteraction$int6.20 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int7.1 <- ifelse(distanceinteraction$int7.1 == 1, c(ifelse(angleinteraction$int7.1 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int7.2 <- ifelse(distanceinteraction$int7.2 == 1, c(ifelse(angleinteraction$int7.2 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int7.3 <- ifelse(distanceinteraction$int7.3 == 1, c(ifelse(angleinteraction$int7.3 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int7.4 <- ifelse(distanceinteraction$int7.4 == 1, c(ifelse(angleinteraction$int7.4 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int7.5 <- ifelse(distanceinteraction$int7.5 == 1, c(ifelse(angleinteraction$int7.5 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int7.6 <- ifelse(distanceinteraction$int7.6 == 1, c(ifelse(angleinteraction$int7.6 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int7.8 <- ifelse(distanceinteraction$int7.8 == 1, c(ifelse(angleinteraction$int7.8 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int7.9 <- ifelse(distanceinteraction$int7.9 == 1, c(ifelse(angleinteraction$int7.9 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int7.10 <- ifelse(distanceinteraction$int7.10 == 1, c(ifelse(angleinteraction$int7.10 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int7.11 <- ifelse(distanceinteraction$int7.11 == 1, c(ifelse(angleinteraction$int7.11 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int7.12 <- ifelse(distanceinteraction$int7.12 == 1, c(ifelse(angleinteraction$int7.12 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int7.13 <- ifelse(distanceinteraction$int7.13 == 1, c(ifelse(angleinteraction$int7.13 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int7.14 <- ifelse(distanceinteraction$int7.14 == 1, c(ifelse(angleinteraction$int7.14 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int7.15 <- ifelse(distanceinteraction$int7.15 == 1, c(ifelse(angleinteraction$int7.15 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int7.16 <- ifelse(distanceinteraction$int7.16 == 1, c(ifelse(angleinteraction$int7.16 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int7.17 <- ifelse(distanceinteraction$int7.17 == 1, c(ifelse(angleinteraction$int7.17 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int7.18 <- ifelse(distanceinteraction$int7.18 == 1, c(ifelse(angleinteraction$int7.18 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int7.19 <- ifelse(distanceinteraction$int7.19 == 1, c(ifelse(angleinteraction$int7.19 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int7.20 <- ifelse(distanceinteraction$int7.20 == 1, c(ifelse(angleinteraction$int7.20 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int8.1 <- ifelse(distanceinteraction$int8.1 == 1, c(ifelse(angleinteraction$int8.1 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int8.2 <- ifelse(distanceinteraction$int8.2 == 1, c(ifelse(angleinteraction$int8.2 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int8.3 <- ifelse(distanceinteraction$int8.3 == 1, c(ifelse(angleinteraction$int8.3 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int8.4 <- ifelse(distanceinteraction$int8.4 == 1, c(ifelse(angleinteraction$int8.4 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int8.5 <- ifelse(distanceinteraction$int8.5 == 1, c(ifelse(angleinteraction$int8.5 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int8.6 <- ifelse(distanceinteraction$int8.6 == 1, c(ifelse(angleinteraction$int8.6 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int8.7 <- ifelse(distanceinteraction$int8.7 == 1, c(ifelse(angleinteraction$int8.7 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int8.9 <- ifelse(distanceinteraction$int8.9 == 1, c(ifelse(angleinteraction$int8.9 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int8.10 <- ifelse(distanceinteraction$int8.10 == 1, c(ifelse(angleinteraction$int8.10 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int8.11 <- ifelse(distanceinteraction$int8.11 == 1, c(ifelse(angleinteraction$int8.11 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int8.12 <- ifelse(distanceinteraction$int8.12 == 1, c(ifelse(angleinteraction$int8.12 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int8.13 <- ifelse(distanceinteraction$int8.13 == 1, c(ifelse(angleinteraction$int8.13 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int8.14 <- ifelse(distanceinteraction$int8.14 == 1, c(ifelse(angleinteraction$int8.14 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int8.15 <- ifelse(distanceinteraction$int8.15 == 1, c(ifelse(angleinteraction$int8.15 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int8.16 <- ifelse(distanceinteraction$int8.16 == 1, c(ifelse(angleinteraction$int8.16 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int8.17 <- ifelse(distanceinteraction$int8.17 == 1, c(ifelse(angleinteraction$int8.17 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int8.18 <- ifelse(distanceinteraction$int8.18 == 1, c(ifelse(angleinteraction$int8.18 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int8.19 <- ifelse(distanceinteraction$int8.19 == 1, c(ifelse(angleinteraction$int8.19 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int8.20 <- ifelse(distanceinteraction$int8.20 == 1, c(ifelse(angleinteraction$int8.20 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int9.1 <- ifelse(distanceinteraction$int9.1 == 1, c(ifelse(angleinteraction$int9.1 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int9.2 <- ifelse(distanceinteraction$int9.2 == 1, c(ifelse(angleinteraction$int9.2 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int9.3 <- ifelse(distanceinteraction$int9.3 == 1, c(ifelse(angleinteraction$int9.3 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int9.4 <- ifelse(distanceinteraction$int9.4 == 1, c(ifelse(angleinteraction$int9.4 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int9.5 <- ifelse(distanceinteraction$int9.5 == 1, c(ifelse(angleinteraction$int9.5 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int9.6 <- ifelse(distanceinteraction$int9.6 == 1, c(ifelse(angleinteraction$int9.6 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int9.7 <- ifelse(distanceinteraction$int9.7 == 1, c(ifelse(angleinteraction$int9.7 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int9.8 <- ifelse(distanceinteraction$int9.8 == 1, c(ifelse(angleinteraction$int9.8 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int9.10 <- ifelse(distanceinteraction$int9.10 == 1, c(ifelse(angleinteraction$int9.10 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int9.11 <- ifelse(distanceinteraction$int9.11 == 1, c(ifelse(angleinteraction$int9.11 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int9.12 <- ifelse(distanceinteraction$int9.12 == 1, c(ifelse(angleinteraction$int9.12 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int9.13 <- ifelse(distanceinteraction$int9.13 == 1, c(ifelse(angleinteraction$int9.13 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int9.14 <- ifelse(distanceinteraction$int9.14 == 1, c(ifelse(angleinteraction$int9.14 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int9.15 <- ifelse(distanceinteraction$int9.15 == 1, c(ifelse(angleinteraction$int9.15 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int9.16 <- ifelse(distanceinteraction$int9.16 == 1, c(ifelse(angleinteraction$int9.16 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int9.17 <- ifelse(distanceinteraction$int9.17 == 1, c(ifelse(angleinteraction$int9.17 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int9.18 <- ifelse(distanceinteraction$int9.18 == 1, c(ifelse(angleinteraction$int9.18 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int9.19 <- ifelse(distanceinteraction$int9.19 == 1, c(ifelse(angleinteraction$int9.19 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int9.20 <- ifelse(distanceinteraction$int9.20 == 1, c(ifelse(angleinteraction$int9.20 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int10.1 <- ifelse(distanceinteraction$int10.1 == 1, c(ifelse(angleinteraction$int10.1 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int10.2 <- ifelse(distanceinteraction$int10.2 == 1, c(ifelse(angleinteraction$int10.2 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int10.3 <- ifelse(distanceinteraction$int10.3 == 1, c(ifelse(angleinteraction$int10.3 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int10.4 <- ifelse(distanceinteraction$int10.4 == 1, c(ifelse(angleinteraction$int10.4 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int10.5 <- ifelse(distanceinteraction$int10.5 == 1, c(ifelse(angleinteraction$int10.5 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int10.6 <- ifelse(distanceinteraction$int10.6 == 1, c(ifelse(angleinteraction$int10.6 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int10.7 <- ifelse(distanceinteraction$int10.7 == 1, c(ifelse(angleinteraction$int10.7 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int10.8 <- ifelse(distanceinteraction$int10.8 == 1, c(ifelse(angleinteraction$int10.8 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int10.9 <- ifelse(distanceinteraction$int10.9 == 1, c(ifelse(angleinteraction$int10.9 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int10.11 <- ifelse(distanceinteraction$int10.11 == 1, c(ifelse(angleinteraction$int10.11 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int10.12 <- ifelse(distanceinteraction$int10.12 == 1, c(ifelse(angleinteraction$int10.12 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int10.13 <- ifelse(distanceinteraction$int10.13 == 1, c(ifelse(angleinteraction$int10.13 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int10.14 <- ifelse(distanceinteraction$int10.14 == 1, c(ifelse(angleinteraction$int10.14 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int10.15 <- ifelse(distanceinteraction$int10.15 == 1, c(ifelse(angleinteraction$int10.15 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int10.16 <- ifelse(distanceinteraction$int10.16 == 1, c(ifelse(angleinteraction$int10.16 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int10.17 <- ifelse(distanceinteraction$int10.17 == 1, c(ifelse(angleinteraction$int10.17 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int10.18 <- ifelse(distanceinteraction$int10.18 == 1, c(ifelse(angleinteraction$int10.18 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int10.19 <- ifelse(distanceinteraction$int10.19 == 1, c(ifelse(angleinteraction$int10.19 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int10.20 <- ifelse(distanceinteraction$int10.20 == 1, c(ifelse(angleinteraction$int10.20 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int11.1 <- ifelse(distanceinteraction$int11.1 == 1, c(ifelse(angleinteraction$int11.1 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int11.2 <- ifelse(distanceinteraction$int11.2 == 1, c(ifelse(angleinteraction$int11.2 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int11.3 <- ifelse(distanceinteraction$int11.3 == 1, c(ifelse(angleinteraction$int11.3 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int11.4 <- ifelse(distanceinteraction$int11.4 == 1, c(ifelse(angleinteraction$int11.4 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int11.5 <- ifelse(distanceinteraction$int11.5 == 1, c(ifelse(angleinteraction$int11.5 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int11.6 <- ifelse(distanceinteraction$int11.6 == 1, c(ifelse(angleinteraction$int11.6 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int11.7 <- ifelse(distanceinteraction$int11.7 == 1, c(ifelse(angleinteraction$int11.7 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int11.8 <- ifelse(distanceinteraction$int11.8 == 1, c(ifelse(angleinteraction$int11.8 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int11.9 <- ifelse(distanceinteraction$int11.9 == 1, c(ifelse(angleinteraction$int11.9 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int11.10 <- ifelse(distanceinteraction$int11.10 == 1, c(ifelse(angleinteraction$int11.10 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int11.12 <- ifelse(distanceinteraction$int11.12 == 1, c(ifelse(angleinteraction$int11.12 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int11.13 <- ifelse(distanceinteraction$int11.13 == 1, c(ifelse(angleinteraction$int11.13 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int11.14 <- ifelse(distanceinteraction$int11.14 == 1, c(ifelse(angleinteraction$int11.14 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int11.15 <- ifelse(distanceinteraction$int11.15 == 1, c(ifelse(angleinteraction$int11.15 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int11.16 <- ifelse(distanceinteraction$int11.16 == 1, c(ifelse(angleinteraction$int11.16 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int11.17 <- ifelse(distanceinteraction$int11.17 == 1, c(ifelse(angleinteraction$int11.17 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int11.18 <- ifelse(distanceinteraction$int11.18 == 1, c(ifelse(angleinteraction$int11.18 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int11.19 <- ifelse(distanceinteraction$int11.19 == 1, c(ifelse(angleinteraction$int11.19 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int11.20 <- ifelse(distanceinteraction$int11.20 == 1, c(ifelse(angleinteraction$int11.20 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int12.1 <- ifelse(distanceinteraction$int12.1 == 1, c(ifelse(angleinteraction$int12.1 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int12.2 <- ifelse(distanceinteraction$int12.2 == 1, c(ifelse(angleinteraction$int12.2 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int12.3 <- ifelse(distanceinteraction$int12.3 == 1, c(ifelse(angleinteraction$int12.3 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int12.4 <- ifelse(distanceinteraction$int12.4 == 1, c(ifelse(angleinteraction$int12.4 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int12.5 <- ifelse(distanceinteraction$int12.5 == 1, c(ifelse(angleinteraction$int12.5 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int12.6 <- ifelse(distanceinteraction$int12.6 == 1, c(ifelse(angleinteraction$int12.6 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int12.7 <- ifelse(distanceinteraction$int12.7 == 1, c(ifelse(angleinteraction$int12.7 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int12.8 <- ifelse(distanceinteraction$int12.8 == 1, c(ifelse(angleinteraction$int12.8 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int12.9 <- ifelse(distanceinteraction$int12.9 == 1, c(ifelse(angleinteraction$int12.9 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int12.10 <- ifelse(distanceinteraction$int12.10 == 1, c(ifelse(angleinteraction$int12.10 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int12.11 <- ifelse(distanceinteraction$int12.11 == 1, c(ifelse(angleinteraction$int12.11 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int12.13 <- ifelse(distanceinteraction$int12.13 == 1, c(ifelse(angleinteraction$int12.13 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int12.14 <- ifelse(distanceinteraction$int12.14 == 1, c(ifelse(angleinteraction$int12.14 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int12.15 <- ifelse(distanceinteraction$int12.15 == 1, c(ifelse(angleinteraction$int12.15 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int12.16 <- ifelse(distanceinteraction$int12.16 == 1, c(ifelse(angleinteraction$int12.16 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int12.17 <- ifelse(distanceinteraction$int12.17 == 1, c(ifelse(angleinteraction$int12.17 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int12.18 <- ifelse(distanceinteraction$int12.18 == 1, c(ifelse(angleinteraction$int12.18 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int12.19 <- ifelse(distanceinteraction$int12.19 == 1, c(ifelse(angleinteraction$int12.19 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int12.20 <- ifelse(distanceinteraction$int12.20 == 1, c(ifelse(angleinteraction$int12.20 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int13.1 <- ifelse(distanceinteraction$int13.1 == 1, c(ifelse(angleinteraction$int13.1 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int13.2 <- ifelse(distanceinteraction$int13.2 == 1, c(ifelse(angleinteraction$int13.2 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int13.3 <- ifelse(distanceinteraction$int13.3 == 1, c(ifelse(angleinteraction$int13.3 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int13.4 <- ifelse(distanceinteraction$int13.4 == 1, c(ifelse(angleinteraction$int13.4 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int13.5 <- ifelse(distanceinteraction$int13.5 == 1, c(ifelse(angleinteraction$int13.5 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int13.6 <- ifelse(distanceinteraction$int13.6 == 1, c(ifelse(angleinteraction$int13.6 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int13.7 <- ifelse(distanceinteraction$int13.7 == 1, c(ifelse(angleinteraction$int13.7 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int13.8 <- ifelse(distanceinteraction$int13.8 == 1, c(ifelse(angleinteraction$int13.8 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int13.9 <- ifelse(distanceinteraction$int13.9 == 1, c(ifelse(angleinteraction$int13.9 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int13.10 <- ifelse(distanceinteraction$int13.10 == 1, c(ifelse(angleinteraction$int13.10 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int13.11 <- ifelse(distanceinteraction$int13.11 == 1, c(ifelse(angleinteraction$int13.11 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int13.12 <- ifelse(distanceinteraction$int13.12 == 1, c(ifelse(angleinteraction$int13.12 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int13.14 <- ifelse(distanceinteraction$int13.14 == 1, c(ifelse(angleinteraction$int13.14 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int13.15 <- ifelse(distanceinteraction$int13.15 == 1, c(ifelse(angleinteraction$int13.15 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int13.16 <- ifelse(distanceinteraction$int13.16 == 1, c(ifelse(angleinteraction$int13.16 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int13.17 <- ifelse(distanceinteraction$int13.17 == 1, c(ifelse(angleinteraction$int13.17 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int13.18 <- ifelse(distanceinteraction$int13.18 == 1, c(ifelse(angleinteraction$int13.18 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int13.19 <- ifelse(distanceinteraction$int13.19 == 1, c(ifelse(angleinteraction$int13.19 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int13.20 <- ifelse(distanceinteraction$int13.20 == 1, c(ifelse(angleinteraction$int13.20 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int14.1 <- ifelse(distanceinteraction$int14.1 == 1, c(ifelse(angleinteraction$int14.1 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int14.2 <- ifelse(distanceinteraction$int14.2 == 1, c(ifelse(angleinteraction$int14.2 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int14.3 <- ifelse(distanceinteraction$int14.3 == 1, c(ifelse(angleinteraction$int14.3 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int14.4 <- ifelse(distanceinteraction$int14.4 == 1, c(ifelse(angleinteraction$int14.4 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int14.5 <- ifelse(distanceinteraction$int14.5 == 1, c(ifelse(angleinteraction$int14.5 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int14.6 <- ifelse(distanceinteraction$int14.6 == 1, c(ifelse(angleinteraction$int14.6 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int14.7 <- ifelse(distanceinteraction$int14.7 == 1, c(ifelse(angleinteraction$int14.7 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int14.8 <- ifelse(distanceinteraction$int14.8 == 1, c(ifelse(angleinteraction$int14.8 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int14.9 <- ifelse(distanceinteraction$int14.9 == 1, c(ifelse(angleinteraction$int14.9 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int14.10 <- ifelse(distanceinteraction$int14.10 == 1, c(ifelse(angleinteraction$int14.10 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int14.11 <- ifelse(distanceinteraction$int14.11 == 1, c(ifelse(angleinteraction$int14.11 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int14.12 <- ifelse(distanceinteraction$int14.12 == 1, c(ifelse(angleinteraction$int14.12 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int14.13 <- ifelse(distanceinteraction$int14.13 == 1, c(ifelse(angleinteraction$int14.13 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int14.15 <- ifelse(distanceinteraction$int14.15 == 1, c(ifelse(angleinteraction$int14.15 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int14.16 <- ifelse(distanceinteraction$int14.16 == 1, c(ifelse(angleinteraction$int14.16 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int14.17 <- ifelse(distanceinteraction$int14.17 == 1, c(ifelse(angleinteraction$int14.17 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int14.18 <- ifelse(distanceinteraction$int14.18 == 1, c(ifelse(angleinteraction$int14.18 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int14.19 <- ifelse(distanceinteraction$int14.19 == 1, c(ifelse(angleinteraction$int14.19 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int14.20 <- ifelse(distanceinteraction$int14.20 == 1, c(ifelse(angleinteraction$int14.20 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int15.1 <- ifelse(distanceinteraction$int15.1 == 1, c(ifelse(angleinteraction$int15.1 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int15.2 <- ifelse(distanceinteraction$int15.2 == 1, c(ifelse(angleinteraction$int15.2 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int15.3 <- ifelse(distanceinteraction$int15.3 == 1, c(ifelse(angleinteraction$int15.3 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int15.4 <- ifelse(distanceinteraction$int15.4 == 1, c(ifelse(angleinteraction$int15.4 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int15.5 <- ifelse(distanceinteraction$int15.5 == 1, c(ifelse(angleinteraction$int15.5 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int15.6 <- ifelse(distanceinteraction$int15.6 == 1, c(ifelse(angleinteraction$int15.6 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int15.7 <- ifelse(distanceinteraction$int15.7 == 1, c(ifelse(angleinteraction$int15.7 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int15.8 <- ifelse(distanceinteraction$int15.8 == 1, c(ifelse(angleinteraction$int15.8 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int15.9 <- ifelse(distanceinteraction$int15.9 == 1, c(ifelse(angleinteraction$int15.9 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int15.10 <- ifelse(distanceinteraction$int15.10 == 1, c(ifelse(angleinteraction$int15.10 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int15.11 <- ifelse(distanceinteraction$int15.11 == 1, c(ifelse(angleinteraction$int15.11 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int15.12 <- ifelse(distanceinteraction$int15.12 == 1, c(ifelse(angleinteraction$int15.12 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int15.13 <- ifelse(distanceinteraction$int15.13 == 1, c(ifelse(angleinteraction$int15.13 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int15.14 <- ifelse(distanceinteraction$int15.14 == 1, c(ifelse(angleinteraction$int15.14 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int15.16 <- ifelse(distanceinteraction$int15.16 == 1, c(ifelse(angleinteraction$int15.16 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int15.17 <- ifelse(distanceinteraction$int15.17 == 1, c(ifelse(angleinteraction$int15.17 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int15.18 <- ifelse(distanceinteraction$int15.18 == 1, c(ifelse(angleinteraction$int15.18 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int15.19 <- ifelse(distanceinteraction$int15.19 == 1, c(ifelse(angleinteraction$int15.19 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int15.20 <- ifelse(distanceinteraction$int15.20 == 1, c(ifelse(angleinteraction$int15.20 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int16.1 <- ifelse(distanceinteraction$int16.1 == 1, c(ifelse(angleinteraction$int16.1 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int16.2 <- ifelse(distanceinteraction$int16.2 == 1, c(ifelse(angleinteraction$int16.2 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int16.3 <- ifelse(distanceinteraction$int16.3 == 1, c(ifelse(angleinteraction$int16.3 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int16.4 <- ifelse(distanceinteraction$int16.4 == 1, c(ifelse(angleinteraction$int16.4 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int16.5 <- ifelse(distanceinteraction$int16.5 == 1, c(ifelse(angleinteraction$int16.5 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int16.6 <- ifelse(distanceinteraction$int16.6 == 1, c(ifelse(angleinteraction$int16.6 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int16.7 <- ifelse(distanceinteraction$int16.7 == 1, c(ifelse(angleinteraction$int16.7 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int16.8 <- ifelse(distanceinteraction$int16.8 == 1, c(ifelse(angleinteraction$int16.8 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int16.9 <- ifelse(distanceinteraction$int16.9 == 1, c(ifelse(angleinteraction$int16.9 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int16.10 <- ifelse(distanceinteraction$int16.10 == 1, c(ifelse(angleinteraction$int16.10 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int16.11 <- ifelse(distanceinteraction$int16.11 == 1, c(ifelse(angleinteraction$int16.11 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int16.12 <- ifelse(distanceinteraction$int16.12 == 1, c(ifelse(angleinteraction$int16.12 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int16.13 <- ifelse(distanceinteraction$int16.13 == 1, c(ifelse(angleinteraction$int16.13 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int16.14 <- ifelse(distanceinteraction$int16.14 == 1, c(ifelse(angleinteraction$int16.14 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int16.15 <- ifelse(distanceinteraction$int16.15 == 1, c(ifelse(angleinteraction$int16.15 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int16.17 <- ifelse(distanceinteraction$int16.17 == 1, c(ifelse(angleinteraction$int16.17 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int16.18 <- ifelse(distanceinteraction$int16.18 == 1, c(ifelse(angleinteraction$int16.18 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int16.19 <- ifelse(distanceinteraction$int16.19 == 1, c(ifelse(angleinteraction$int16.19 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int16.20 <- ifelse(distanceinteraction$int16.20 == 1, c(ifelse(angleinteraction$int16.20 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int17.1 <- ifelse(distanceinteraction$int17.1 == 1, c(ifelse(angleinteraction$int17.1 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int17.2 <- ifelse(distanceinteraction$int17.2 == 1, c(ifelse(angleinteraction$int17.2 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int17.3 <- ifelse(distanceinteraction$int17.3 == 1, c(ifelse(angleinteraction$int17.3 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int17.4 <- ifelse(distanceinteraction$int17.4 == 1, c(ifelse(angleinteraction$int17.4 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int17.5 <- ifelse(distanceinteraction$int17.5 == 1, c(ifelse(angleinteraction$int17.5 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int17.6 <- ifelse(distanceinteraction$int17.6 == 1, c(ifelse(angleinteraction$int17.6 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int17.7 <- ifelse(distanceinteraction$int17.7 == 1, c(ifelse(angleinteraction$int17.7 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int17.8 <- ifelse(distanceinteraction$int17.8 == 1, c(ifelse(angleinteraction$int17.8 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int17.9 <- ifelse(distanceinteraction$int17.9 == 1, c(ifelse(angleinteraction$int17.9 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int17.10 <- ifelse(distanceinteraction$int17.10 == 1, c(ifelse(angleinteraction$int17.10 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int17.11 <- ifelse(distanceinteraction$int17.11 == 1, c(ifelse(angleinteraction$int17.11 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int17.12 <- ifelse(distanceinteraction$int17.12 == 1, c(ifelse(angleinteraction$int17.12 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int17.13 <- ifelse(distanceinteraction$int17.13 == 1, c(ifelse(angleinteraction$int17.13 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int17.14 <- ifelse(distanceinteraction$int17.14 == 1, c(ifelse(angleinteraction$int17.14 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int17.15 <- ifelse(distanceinteraction$int17.15 == 1, c(ifelse(angleinteraction$int17.15 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int17.16 <- ifelse(distanceinteraction$int17.16 == 1, c(ifelse(angleinteraction$int17.16 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int17.18 <- ifelse(distanceinteraction$int17.18 == 1, c(ifelse(angleinteraction$int17.18 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int17.19 <- ifelse(distanceinteraction$int17.19 == 1, c(ifelse(angleinteraction$int17.19 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int17.20 <- ifelse(distanceinteraction$int17.20 == 1, c(ifelse(angleinteraction$int17.20 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int18.1 <- ifelse(distanceinteraction$int18.1 == 1, c(ifelse(angleinteraction$int18.1 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int18.2 <- ifelse(distanceinteraction$int18.2 == 1, c(ifelse(angleinteraction$int18.2 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int18.3 <- ifelse(distanceinteraction$int18.3 == 1, c(ifelse(angleinteraction$int18.3 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int18.4 <- ifelse(distanceinteraction$int18.4 == 1, c(ifelse(angleinteraction$int18.4 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int18.5 <- ifelse(distanceinteraction$int18.5 == 1, c(ifelse(angleinteraction$int18.5 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int18.6 <- ifelse(distanceinteraction$int18.6 == 1, c(ifelse(angleinteraction$int18.6 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int18.7 <- ifelse(distanceinteraction$int18.7 == 1, c(ifelse(angleinteraction$int18.7 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int18.8 <- ifelse(distanceinteraction$int18.8 == 1, c(ifelse(angleinteraction$int18.8 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int18.9 <- ifelse(distanceinteraction$int18.9 == 1, c(ifelse(angleinteraction$int18.9 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int18.10 <- ifelse(distanceinteraction$int18.10 == 1, c(ifelse(angleinteraction$int18.10 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int18.11 <- ifelse(distanceinteraction$int18.11 == 1, c(ifelse(angleinteraction$int18.11 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int18.12 <- ifelse(distanceinteraction$int18.12 == 1, c(ifelse(angleinteraction$int18.12 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int18.13 <- ifelse(distanceinteraction$int18.13 == 1, c(ifelse(angleinteraction$int18.13 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int18.14 <- ifelse(distanceinteraction$int18.14 == 1, c(ifelse(angleinteraction$int18.14 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int18.15 <- ifelse(distanceinteraction$int18.15 == 1, c(ifelse(angleinteraction$int18.15 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int18.16 <- ifelse(distanceinteraction$int18.16 == 1, c(ifelse(angleinteraction$int18.16 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int18.17 <- ifelse(distanceinteraction$int18.17 == 1, c(ifelse(angleinteraction$int18.17 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int18.19 <- ifelse(distanceinteraction$int18.19 == 1, c(ifelse(angleinteraction$int18.19 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int18.20 <- ifelse(distanceinteraction$int18.20 == 1, c(ifelse(angleinteraction$int18.20 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int19.1 <- ifelse(distanceinteraction$int19.1 == 1, c(ifelse(angleinteraction$int19.1 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int19.2 <- ifelse(distanceinteraction$int19.2 == 1, c(ifelse(angleinteraction$int19.2 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int19.3 <- ifelse(distanceinteraction$int19.3 == 1, c(ifelse(angleinteraction$int19.3 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int19.4 <- ifelse(distanceinteraction$int19.4 == 1, c(ifelse(angleinteraction$int19.4 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int19.5 <- ifelse(distanceinteraction$int19.5 == 1, c(ifelse(angleinteraction$int19.5 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int19.6 <- ifelse(distanceinteraction$int19.6 == 1, c(ifelse(angleinteraction$int19.6 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int19.7 <- ifelse(distanceinteraction$int19.7 == 1, c(ifelse(angleinteraction$int19.7 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int19.8 <- ifelse(distanceinteraction$int19.8 == 1, c(ifelse(angleinteraction$int19.8 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int19.9 <- ifelse(distanceinteraction$int19.9 == 1, c(ifelse(angleinteraction$int19.9 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int19.10 <- ifelse(distanceinteraction$int19.10 == 1, c(ifelse(angleinteraction$int19.10 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int19.11 <- ifelse(distanceinteraction$int19.11 == 1, c(ifelse(angleinteraction$int19.11 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int19.12 <- ifelse(distanceinteraction$int19.12 == 1, c(ifelse(angleinteraction$int19.12 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int19.13 <- ifelse(distanceinteraction$int19.13 == 1, c(ifelse(angleinteraction$int19.13 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int19.14 <- ifelse(distanceinteraction$int19.14 == 1, c(ifelse(angleinteraction$int19.14 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int19.15 <- ifelse(distanceinteraction$int19.15 == 1, c(ifelse(angleinteraction$int19.15 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int19.16 <- ifelse(distanceinteraction$int19.16 == 1, c(ifelse(angleinteraction$int19.16 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int19.17 <- ifelse(distanceinteraction$int19.17 == 1, c(ifelse(angleinteraction$int19.17 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int19.18 <- ifelse(distanceinteraction$int19.18 == 1, c(ifelse(angleinteraction$int19.18 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int19.20 <- ifelse(distanceinteraction$int19.20 == 1, c(ifelse(angleinteraction$int19.20 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int20.1 <- ifelse(distanceinteraction$int20.1 == 1, c(ifelse(angleinteraction$int20.1 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int20.2 <- ifelse(distanceinteraction$int20.2 == 1, c(ifelse(angleinteraction$int20.2 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int20.3 <- ifelse(distanceinteraction$int20.3 == 1, c(ifelse(angleinteraction$int20.3 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int20.4 <- ifelse(distanceinteraction$int20.4 == 1, c(ifelse(angleinteraction$int20.4 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int20.5 <- ifelse(distanceinteraction$int20.5 == 1, c(ifelse(angleinteraction$int20.5 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int20.6 <- ifelse(distanceinteraction$int20.6 == 1, c(ifelse(angleinteraction$int20.6 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int20.7 <- ifelse(distanceinteraction$int20.7 == 1, c(ifelse(angleinteraction$int20.7 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int20.8 <- ifelse(distanceinteraction$int20.8 == 1, c(ifelse(angleinteraction$int20.8 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int20.9 <- ifelse(distanceinteraction$int20.9 == 1, c(ifelse(angleinteraction$int20.9 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int20.10 <- ifelse(distanceinteraction$int20.10 == 1, c(ifelse(angleinteraction$int20.10 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int20.11 <- ifelse(distanceinteraction$int20.11 == 1, c(ifelse(angleinteraction$int20.11 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int20.12 <- ifelse(distanceinteraction$int20.12 == 1, c(ifelse(angleinteraction$int20.12 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int20.13 <- ifelse(distanceinteraction$int20.13 == 1, c(ifelse(angleinteraction$int20.13 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int20.14 <- ifelse(distanceinteraction$int20.14 == 1, c(ifelse(angleinteraction$int20.14 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int20.15 <- ifelse(distanceinteraction$int20.15 == 1, c(ifelse(angleinteraction$int20.15 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int20.16 <- ifelse(distanceinteraction$int20.16 == 1, c(ifelse(angleinteraction$int20.16 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int20.17 <- ifelse(distanceinteraction$int20.17 == 1, c(ifelse(angleinteraction$int20.17 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int20.18 <- ifelse(distanceinteraction$int20.18 == 1, c(ifelse(angleinteraction$int20.18 == 1, c(1), c(0))), c(0))
  distanceangleinteraction$int20.19 <- ifelse(distanceinteraction$int20.19 == 1, c(ifelse(angleinteraction$int20.19 == 1, c(1), c(0))), c(0))
  
  # CREATE AN ASSOCIATION MATRIX WHERE INTERACTIONS LASTING LESS THAN 0.6 SECONDS (18 FRAMES) ARE REMOVED -----------
  associationmatrix <- setnames(setDF(lapply(integer(20), function(...) character(0L))),paste0(1:20))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int1.2), ", ", ""), "01{1,17}0", "0")
  associationmatrix[1,2]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int1.3), ", ", ""), "01{1,17}0", "0")
  associationmatrix[1,3]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int1.4), ", ", ""), "01{1,17}0", "0")
  associationmatrix[1,4]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int1.5), ", ", ""), "01{1,17}0", "0")
  associationmatrix[1,5]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int1.6), ", ", ""), "01{1,17}0", "0")
  associationmatrix[1,6]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int1.7), ", ", ""), "01{1,17}0", "0")
  associationmatrix[1,7]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int1.8), ", ", ""), "01{1,17}0", "0")
  associationmatrix[1,8]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int1.9), ", ", ""), "01{1,17}0", "0")
  associationmatrix[1,9]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int1.10), ", ", ""), "01{1,17}0", "0")
  associationmatrix[1,10]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int1.11), ", ", ""), "01{1,17}0", "0")
  associationmatrix[1,11]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int1.12), ", ", ""), "01{1,17}0", "0")
  associationmatrix[1,12]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int1.13), ", ", ""), "01{1,17}0", "0")
  associationmatrix[1,13]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int1.14), ", ", ""), "01{1,17}0", "0")
  associationmatrix[1,14]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int1.15), ", ", ""), "01{1,17}0", "0")
  associationmatrix[1,15]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int1.16), ", ", ""), "01{1,17}0", "0")
  associationmatrix[1,16]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int1.17), ", ", ""), "01{1,17}0", "0")
  associationmatrix[1,17]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int1.18), ", ", ""), "01{1,17}0", "0")
  associationmatrix[1,18]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int1.19), ", ", ""), "01{1,17}0", "0")
  associationmatrix[1,19]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int1.20), ", ", ""), "01{1,17}0", "0")
  associationmatrix[1,20]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int2.1), ", ", ""), "01{1,17}0", "0")
  associationmatrix[2,1]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int2.3), ", ", ""), "01{1,17}0", "0")
  associationmatrix[2,3]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int2.4), ", ", ""), "01{1,17}0", "0")
  associationmatrix[2,4]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int2.5), ", ", ""), "01{1,17}0", "0")
  associationmatrix[2,5]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int2.6), ", ", ""), "01{1,17}0", "0")
  associationmatrix[2,6]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int2.7), ", ", ""), "01{1,17}0", "0")
  associationmatrix[2,7]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int2.8), ", ", ""), "01{1,17}0", "0")
  associationmatrix[2,8]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int2.9), ", ", ""), "01{1,17}0", "0")
  associationmatrix[2,9]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int2.10), ", ", ""), "01{1,17}0", "0")
  associationmatrix[2,10]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int2.11), ", ", ""), "01{1,17}0", "0")
  associationmatrix[2,11]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int2.12), ", ", ""), "01{1,17}0", "0")
  associationmatrix[2,12]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int2.13), ", ", ""), "01{1,17}0", "0")
  associationmatrix[2,13]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int2.14), ", ", ""), "01{1,17}0", "0")
  associationmatrix[2,14]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int2.15), ", ", ""), "01{1,17}0", "0")
  associationmatrix[2,15]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int2.16), ", ", ""), "01{1,17}0", "0")
  associationmatrix[2,16]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int2.17), ", ", ""), "01{1,17}0", "0")
  associationmatrix[2,17]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int2.18), ", ", ""), "01{1,17}0", "0")
  associationmatrix[2,18]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int2.19), ", ", ""), "01{1,17}0", "0")
  associationmatrix[2,19]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int2.20), ", ", ""), "01{1,17}0", "0")
  associationmatrix[2,20]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int3.1), ", ", ""), "01{1,17}0", "0")
  associationmatrix[3,1]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int3.2), ", ", ""), "01{1,17}0", "0")
  associationmatrix[3,2]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int3.4), ", ", ""), "01{1,17}0", "0")
  associationmatrix[3,4]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int3.5), ", ", ""), "01{1,17}0", "0")
  associationmatrix[3,5]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int3.6), ", ", ""), "01{1,17}0", "0")
  associationmatrix[3,6]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int3.7), ", ", ""), "01{1,17}0", "0")
  associationmatrix[3,7]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int3.8), ", ", ""), "01{1,17}0", "0")
  associationmatrix[3,8]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int3.9), ", ", ""), "01{1,17}0", "0")
  associationmatrix[3,9]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int3.10), ", ", ""), "01{1,17}0", "0")
  associationmatrix[3,10]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int3.11), ", ", ""), "01{1,17}0", "0")
  associationmatrix[3,11]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int3.12), ", ", ""), "01{1,17}0", "0")
  associationmatrix[3,12]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int3.13), ", ", ""), "01{1,17}0", "0")
  associationmatrix[3,13]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int3.14), ", ", ""), "01{1,17}0", "0")
  associationmatrix[3,14]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int3.15), ", ", ""), "01{1,17}0", "0")
  associationmatrix[3,15]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int3.16), ", ", ""), "01{1,17}0", "0")
  associationmatrix[3,16]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int3.17), ", ", ""), "01{1,17}0", "0")
  associationmatrix[3,17]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int3.18), ", ", ""), "01{1,17}0", "0")
  associationmatrix[3,18]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int3.19), ", ", ""), "01{1,17}0", "0")
  associationmatrix[3,19]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int3.20), ", ", ""), "01{1,17}0", "0")
  associationmatrix[3,20]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int4.1), ", ", ""), "01{1,17}0", "0")
  associationmatrix[4,1]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int4.2), ", ", ""), "01{1,17}0", "0")
  associationmatrix[4,2]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int4.3), ", ", ""), "01{1,17}0", "0")
  associationmatrix[4,3]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int4.5), ", ", ""), "01{1,17}0", "0")
  associationmatrix[4,5]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int4.6), ", ", ""), "01{1,17}0", "0")
  associationmatrix[4,6]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int4.7), ", ", ""), "01{1,17}0", "0")
  associationmatrix[4,7]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int4.8), ", ", ""), "01{1,17}0", "0")
  associationmatrix[4,8]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int4.9), ", ", ""), "01{1,17}0", "0")
  associationmatrix[4,9]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int4.10), ", ", ""), "01{1,17}0", "0")
  associationmatrix[4,10]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int4.11), ", ", ""), "01{1,17}0", "0")
  associationmatrix[4,11]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int4.12), ", ", ""), "01{1,17}0", "0")
  associationmatrix[4,12]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int4.13), ", ", ""), "01{1,17}0", "0")
  associationmatrix[4,13]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int4.14), ", ", ""), "01{1,17}0", "0")
  associationmatrix[4,14]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int4.15), ", ", ""), "01{1,17}0", "0")
  associationmatrix[4,15]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int4.16), ", ", ""), "01{1,17}0", "0")
  associationmatrix[4,16]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int4.17), ", ", ""), "01{1,17}0", "0")
  associationmatrix[4,17]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int4.18), ", ", ""), "01{1,17}0", "0")
  associationmatrix[4,18]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int4.19), ", ", ""), "01{1,17}0", "0")
  associationmatrix[4,19]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int4.20), ", ", ""), "01{1,17}0", "0")
  associationmatrix[4,20]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int5.1), ", ", ""), "01{1,17}0", "0")
  associationmatrix[5,1]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int5.2), ", ", ""), "01{1,17}0", "0")
  associationmatrix[5,2]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int5.3), ", ", ""), "01{1,17}0", "0")
  associationmatrix[5,3]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int5.4), ", ", ""), "01{1,17}0", "0")
  associationmatrix[5,4]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int5.6), ", ", ""), "01{1,17}0", "0")
  associationmatrix[5,6]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int5.7), ", ", ""), "01{1,17}0", "0")
  associationmatrix[5,7]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int5.8), ", ", ""), "01{1,17}0", "0")
  associationmatrix[5,8]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int5.9), ", ", ""), "01{1,17}0", "0")
  associationmatrix[5,9]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int5.10), ", ", ""), "01{1,17}0", "0")
  associationmatrix[5,10]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int5.11), ", ", ""), "01{1,17}0", "0")
  associationmatrix[5,11]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int5.12), ", ", ""), "01{1,17}0", "0")
  associationmatrix[5,12]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int5.13), ", ", ""), "01{1,17}0", "0")
  associationmatrix[5,13]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int5.14), ", ", ""), "01{1,17}0", "0")
  associationmatrix[5,14]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int5.15), ", ", ""), "01{1,17}0", "0")
  associationmatrix[5,15]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int5.16), ", ", ""), "01{1,17}0", "0")
  associationmatrix[5,16]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int5.17), ", ", ""), "01{1,17}0", "0")
  associationmatrix[5,17]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int5.18), ", ", ""), "01{1,17}0", "0")
  associationmatrix[5,18]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int5.19), ", ", ""), "01{1,17}0", "0")
  associationmatrix[5,19]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int5.20), ", ", ""), "01{1,17}0", "0")
  associationmatrix[5,20]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int6.1), ", ", ""), "01{1,17}0", "0")
  associationmatrix[6,1]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int6.2), ", ", ""), "01{1,17}0", "0")
  associationmatrix[6,2]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int6.3), ", ", ""), "01{1,17}0", "0")
  associationmatrix[6,3]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int6.4), ", ", ""), "01{1,17}0", "0")
  associationmatrix[6,4]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int6.5), ", ", ""), "01{1,17}0", "0")
  associationmatrix[6,5]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int6.7), ", ", ""), "01{1,17}0", "0")
  associationmatrix[6,7]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int6.8), ", ", ""), "01{1,17}0", "0")
  associationmatrix[6,8]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int6.9), ", ", ""), "01{1,17}0", "0")
  associationmatrix[6,9]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int6.10), ", ", ""), "01{1,17}0", "0")
  associationmatrix[6,10]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int6.11), ", ", ""), "01{1,17}0", "0")
  associationmatrix[6,11]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int6.12), ", ", ""), "01{1,17}0", "0")
  associationmatrix[6,12]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int6.13), ", ", ""), "01{1,17}0", "0")
  associationmatrix[6,13]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int6.14), ", ", ""), "01{1,17}0", "0")
  associationmatrix[6,14]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int6.15), ", ", ""), "01{1,17}0", "0")
  associationmatrix[6,15]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int6.16), ", ", ""), "01{1,17}0", "0")
  associationmatrix[6,16]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int6.17), ", ", ""), "01{1,17}0", "0")
  associationmatrix[6,17]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int6.18), ", ", ""), "01{1,17}0", "0")
  associationmatrix[6,18]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int6.19), ", ", ""), "01{1,17}0", "0")
  associationmatrix[6,19]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int6.20), ", ", ""), "01{1,17}0", "0")
  associationmatrix[6,20]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int7.1), ", ", ""), "01{1,17}0", "0")
  associationmatrix[7,1]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int7.2), ", ", ""), "01{1,17}0", "0")
  associationmatrix[7,2]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int7.3), ", ", ""), "01{1,17}0", "0")
  associationmatrix[7,3]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int7.4), ", ", ""), "01{1,17}0", "0")
  associationmatrix[7,4]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int7.5), ", ", ""), "01{1,17}0", "0")
  associationmatrix[7,5]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int7.6), ", ", ""), "01{1,17}0", "0")
  associationmatrix[7,6]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int7.8), ", ", ""), "01{1,17}0", "0")
  associationmatrix[7,8]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int7.9), ", ", ""), "01{1,17}0", "0")
  associationmatrix[7,9]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int7.10), ", ", ""), "01{1,17}0", "0")
  associationmatrix[7,10]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int7.11), ", ", ""), "01{1,17}0", "0")
  associationmatrix[7,11]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int7.12), ", ", ""), "01{1,17}0", "0")
  associationmatrix[7,12]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int7.13), ", ", ""), "01{1,17}0", "0")
  associationmatrix[7,13]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int7.14), ", ", ""), "01{1,17}0", "0")
  associationmatrix[7,14]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int7.15), ", ", ""), "01{1,17}0", "0")
  associationmatrix[7,15]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int7.16), ", ", ""), "01{1,17}0", "0")
  associationmatrix[7,16]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int7.17), ", ", ""), "01{1,17}0", "0")
  associationmatrix[7,17]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int7.18), ", ", ""), "01{1,17}0", "0")
  associationmatrix[7,18]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int7.19), ", ", ""), "01{1,17}0", "0")
  associationmatrix[7,19]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int7.20), ", ", ""), "01{1,17}0", "0")
  associationmatrix[7,20]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int8.1), ", ", ""), "01{1,17}0", "0")
  associationmatrix[8,1]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int8.2), ", ", ""), "01{1,17}0", "0")
  associationmatrix[8,2]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int8.3), ", ", ""), "01{1,17}0", "0")
  associationmatrix[8,3]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int8.4), ", ", ""), "01{1,17}0", "0")
  associationmatrix[8,4]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int8.5), ", ", ""), "01{1,17}0", "0")
  associationmatrix[8,5]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int8.6), ", ", ""), "01{1,17}0", "0")
  associationmatrix[8,6]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int8.7), ", ", ""), "01{1,17}0", "0")
  associationmatrix[8,7]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int8.9), ", ", ""), "01{1,17}0", "0")
  associationmatrix[8,9]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int8.10), ", ", ""), "01{1,17}0", "0")
  associationmatrix[8,10]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int8.11), ", ", ""), "01{1,17}0", "0")
  associationmatrix[8,11]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int8.12), ", ", ""), "01{1,17}0", "0")
  associationmatrix[8,12]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int8.13), ", ", ""), "01{1,17}0", "0")
  associationmatrix[8,13]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int8.14), ", ", ""), "01{1,17}0", "0")
  associationmatrix[8,14]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int8.15), ", ", ""), "01{1,17}0", "0")
  associationmatrix[8,15]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int8.16), ", ", ""), "01{1,17}0", "0")
  associationmatrix[8,16]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int8.17), ", ", ""), "01{1,17}0", "0")
  associationmatrix[8,17]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int8.18), ", ", ""), "01{1,17}0", "0")
  associationmatrix[8,18]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int8.19), ", ", ""), "01{1,17}0", "0")
  associationmatrix[8,19]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int8.20), ", ", ""), "01{1,17}0", "0")
  associationmatrix[8,20]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int9.1), ", ", ""), "01{1,17}0", "0")
  associationmatrix[9,1]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int9.2), ", ", ""), "01{1,17}0", "0")
  associationmatrix[9,2]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int9.3), ", ", ""), "01{1,17}0", "0")
  associationmatrix[9,3]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int9.4), ", ", ""), "01{1,17}0", "0")
  associationmatrix[9,4]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int9.5), ", ", ""), "01{1,17}0", "0")
  associationmatrix[9,5]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int9.6), ", ", ""), "01{1,17}0", "0")
  associationmatrix[9,6]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int9.7), ", ", ""), "01{1,17}0", "0")
  associationmatrix[9,7]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int9.8), ", ", ""), "01{1,17}0", "0")
  associationmatrix[9,8]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int9.10), ", ", ""), "01{1,17}0", "0")
  associationmatrix[9,10]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int9.11), ", ", ""), "01{1,17}0", "0")
  associationmatrix[9,11]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int9.12), ", ", ""), "01{1,17}0", "0")
  associationmatrix[9,12]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int9.13), ", ", ""), "01{1,17}0", "0")
  associationmatrix[9,13]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int9.14), ", ", ""), "01{1,17}0", "0")
  associationmatrix[9,14]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int9.15), ", ", ""), "01{1,17}0", "0")
  associationmatrix[9,15]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int9.16), ", ", ""), "01{1,17}0", "0")
  associationmatrix[9,16]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int9.17), ", ", ""), "01{1,17}0", "0")
  associationmatrix[9,17]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int9.18), ", ", ""), "01{1,17}0", "0")
  associationmatrix[9,18]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int9.19), ", ", ""), "01{1,17}0", "0")
  associationmatrix[9,19]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int9.20), ", ", ""), "01{1,17}0", "0")
  associationmatrix[9,20]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int10.1), ", ", ""), "01{1,17}0", "0")
  associationmatrix[10,1]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int10.2), ", ", ""), "01{1,17}0", "0")
  associationmatrix[10,2]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int10.3), ", ", ""), "01{1,17}0", "0")
  associationmatrix[10,3]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int10.4), ", ", ""), "01{1,17}0", "0")
  associationmatrix[10,4]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int10.5), ", ", ""), "01{1,17}0", "0")
  associationmatrix[10,5]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int10.6), ", ", ""), "01{1,17}0", "0")
  associationmatrix[10,6]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int10.7), ", ", ""), "01{1,17}0", "0")
  associationmatrix[10,7]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int10.8), ", ", ""), "01{1,17}0", "0")
  associationmatrix[10,8]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int10.9), ", ", ""), "01{1,17}0", "0")
  associationmatrix[10,9]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int10.11), ", ", ""), "01{1,17}0", "0")
  associationmatrix[10,11]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int10.12), ", ", ""), "01{1,17}0", "0")
  associationmatrix[10,12]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int10.13), ", ", ""), "01{1,17}0", "0")
  associationmatrix[10,13]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int10.14), ", ", ""), "01{1,17}0", "0")
  associationmatrix[10,14]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int10.15), ", ", ""), "01{1,17}0", "0")
  associationmatrix[10,15]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int10.16), ", ", ""), "01{1,17}0", "0")
  associationmatrix[10,16]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int10.17), ", ", ""), "01{1,17}0", "0")
  associationmatrix[10,17]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int10.18), ", ", ""), "01{1,17}0", "0")
  associationmatrix[10,18]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int10.19), ", ", ""), "01{1,17}0", "0")
  associationmatrix[10,19]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int10.20), ", ", ""), "01{1,17}0", "0")
  associationmatrix[10,20]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int11.1), ", ", ""), "01{1,17}0", "0")
  associationmatrix[11,1]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int11.2), ", ", ""), "01{1,17}0", "0")
  associationmatrix[11,2]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int11.3), ", ", ""), "01{1,17}0", "0")
  associationmatrix[11,3]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int11.4), ", ", ""), "01{1,17}0", "0")
  associationmatrix[11,4]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int11.5), ", ", ""), "01{1,17}0", "0")
  associationmatrix[11,5]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int11.6), ", ", ""), "01{1,17}0", "0")
  associationmatrix[11,6]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int11.7), ", ", ""), "01{1,17}0", "0")
  associationmatrix[11,7]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int11.8), ", ", ""), "01{1,17}0", "0")
  associationmatrix[11,8]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int11.9), ", ", ""), "01{1,17}0", "0")
  associationmatrix[11,9]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int11.10), ", ", ""), "01{1,17}0", "0")
  associationmatrix[11,10]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int11.12), ", ", ""), "01{1,17}0", "0")
  associationmatrix[11,12]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int11.13), ", ", ""), "01{1,17}0", "0")
  associationmatrix[11,13]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int11.14), ", ", ""), "01{1,17}0", "0")
  associationmatrix[11,14]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int11.15), ", ", ""), "01{1,17}0", "0")
  associationmatrix[11,15]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int11.16), ", ", ""), "01{1,17}0", "0")
  associationmatrix[11,16]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int11.17), ", ", ""), "01{1,17}0", "0")
  associationmatrix[11,17]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int11.18), ", ", ""), "01{1,17}0", "0")
  associationmatrix[11,18]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int11.19), ", ", ""), "01{1,17}0", "0")
  associationmatrix[11,19]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int11.20), ", ", ""), "01{1,17}0", "0")
  associationmatrix[11,20]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int12.1), ", ", ""), "01{1,17}0", "0")
  associationmatrix[12,1]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int12.2), ", ", ""), "01{1,17}0", "0")
  associationmatrix[12,2]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int12.3), ", ", ""), "01{1,17}0", "0")
  associationmatrix[12,3]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int12.4), ", ", ""), "01{1,17}0", "0")
  associationmatrix[12,4]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int12.5), ", ", ""), "01{1,17}0", "0")
  associationmatrix[12,5]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int12.6), ", ", ""), "01{1,17}0", "0")
  associationmatrix[12,6]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int12.7), ", ", ""), "01{1,17}0", "0")
  associationmatrix[12,7]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int12.8), ", ", ""), "01{1,17}0", "0")
  associationmatrix[12,8]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int12.9), ", ", ""), "01{1,17}0", "0")
  associationmatrix[12,9]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int12.10), ", ", ""), "01{1,17}0", "0")
  associationmatrix[12,10]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int12.11), ", ", ""), "01{1,17}0", "0")
  associationmatrix[12,11]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int12.13), ", ", ""), "01{1,17}0", "0")
  associationmatrix[12,13]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int12.14), ", ", ""), "01{1,17}0", "0")
  associationmatrix[12,14]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int12.15), ", ", ""), "01{1,17}0", "0")
  associationmatrix[12,15]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int12.16), ", ", ""), "01{1,17}0", "0")
  associationmatrix[12,16]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int12.17), ", ", ""), "01{1,17}0", "0")
  associationmatrix[12,17]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int12.18), ", ", ""), "01{1,17}0", "0")
  associationmatrix[12,18]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int12.19), ", ", ""), "01{1,17}0", "0")
  associationmatrix[12,19]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int12.20), ", ", ""), "01{1,17}0", "0")
  associationmatrix[12,20]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int13.1), ", ", ""), "01{1,17}0", "0")
  associationmatrix[13,1]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int13.2), ", ", ""), "01{1,17}0", "0")
  associationmatrix[13,2]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int13.3), ", ", ""), "01{1,17}0", "0")
  associationmatrix[13,3]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int13.4), ", ", ""), "01{1,17}0", "0")
  associationmatrix[13,4]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int13.5), ", ", ""), "01{1,17}0", "0")
  associationmatrix[13,5]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int13.6), ", ", ""), "01{1,17}0", "0")
  associationmatrix[13,6]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int13.7), ", ", ""), "01{1,17}0", "0")
  associationmatrix[13,7]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int13.8), ", ", ""), "01{1,17}0", "0")
  associationmatrix[13,8]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int13.9), ", ", ""), "01{1,17}0", "0")
  associationmatrix[13,9]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int13.10), ", ", ""), "01{1,17}0", "0")
  associationmatrix[13,10]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int13.11), ", ", ""), "01{1,17}0", "0")
  associationmatrix[13,11]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int13.12), ", ", ""), "01{1,17}0", "0")
  associationmatrix[13,12]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int13.14), ", ", ""), "01{1,17}0", "0")
  associationmatrix[13,14]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int13.15), ", ", ""), "01{1,17}0", "0")
  associationmatrix[13,15]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int13.16), ", ", ""), "01{1,17}0", "0")
  associationmatrix[13,16]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int13.17), ", ", ""), "01{1,17}0", "0")
  associationmatrix[13,17]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int13.18), ", ", ""), "01{1,17}0", "0")
  associationmatrix[13,18]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int13.19), ", ", ""), "01{1,17}0", "0")
  associationmatrix[13,19]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int13.20), ", ", ""), "01{1,17}0", "0")
  associationmatrix[13,20]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int14.1), ", ", ""), "01{1,17}0", "0")
  associationmatrix[14,1]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int14.2), ", ", ""), "01{1,17}0", "0")
  associationmatrix[14,2]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int14.3), ", ", ""), "01{1,17}0", "0")
  associationmatrix[14,3]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int14.4), ", ", ""), "01{1,17}0", "0")
  associationmatrix[14,4]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int14.5), ", ", ""), "01{1,17}0", "0")
  associationmatrix[14,5]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int14.6), ", ", ""), "01{1,17}0", "0")
  associationmatrix[14,6]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int14.7), ", ", ""), "01{1,17}0", "0")
  associationmatrix[14,7]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int14.8), ", ", ""), "01{1,17}0", "0")
  associationmatrix[14,8]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int14.9), ", ", ""), "01{1,17}0", "0")
  associationmatrix[14,9]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int14.10), ", ", ""), "01{1,17}0", "0")
  associationmatrix[14,10]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int14.11), ", ", ""), "01{1,17}0", "0")
  associationmatrix[14,11]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int14.12), ", ", ""), "01{1,17}0", "0")
  associationmatrix[14,12]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int14.13), ", ", ""), "01{1,17}0", "0")
  associationmatrix[14,13]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int14.15), ", ", ""), "01{1,17}0", "0")
  associationmatrix[14,15]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int14.16), ", ", ""), "01{1,17}0", "0")
  associationmatrix[14,16]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int14.17), ", ", ""), "01{1,17}0", "0")
  associationmatrix[14,17]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int14.18), ", ", ""), "01{1,17}0", "0")
  associationmatrix[14,18]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int14.19), ", ", ""), "01{1,17}0", "0")
  associationmatrix[14,19]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int14.20), ", ", ""), "01{1,17}0", "0")
  associationmatrix[14,20]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int15.1), ", ", ""), "01{1,17}0", "0")
  associationmatrix[15,1]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int15.2), ", ", ""), "01{1,17}0", "0")
  associationmatrix[15,2]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int15.3), ", ", ""), "01{1,17}0", "0")
  associationmatrix[15,3]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int15.4), ", ", ""), "01{1,17}0", "0")
  associationmatrix[15,4]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int15.5), ", ", ""), "01{1,17}0", "0")
  associationmatrix[15,5]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int15.6), ", ", ""), "01{1,17}0", "0")
  associationmatrix[15,6]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int15.7), ", ", ""), "01{1,17}0", "0")
  associationmatrix[15,7]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int15.8), ", ", ""), "01{1,17}0", "0")
  associationmatrix[15,8]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int15.9), ", ", ""), "01{1,17}0", "0")
  associationmatrix[15,9]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int15.10), ", ", ""), "01{1,17}0", "0")
  associationmatrix[15,10]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int15.11), ", ", ""), "01{1,17}0", "0")
  associationmatrix[15,11]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int15.12), ", ", ""), "01{1,17}0", "0")
  associationmatrix[15,12]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int15.13), ", ", ""), "01{1,17}0", "0")
  associationmatrix[15,13]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int15.14), ", ", ""), "01{1,17}0", "0")
  associationmatrix[15,14]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int15.16), ", ", ""), "01{1,17}0", "0")
  associationmatrix[15,16]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int15.17), ", ", ""), "01{1,17}0", "0")
  associationmatrix[15,17]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int15.18), ", ", ""), "01{1,17}0", "0")
  associationmatrix[15,18]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int15.19), ", ", ""), "01{1,17}0", "0")
  associationmatrix[15,19]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int15.20), ", ", ""), "01{1,17}0", "0")
  associationmatrix[15,20]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int16.1), ", ", ""), "01{1,17}0", "0")
  associationmatrix[16,1]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int16.2), ", ", ""), "01{1,17}0", "0")
  associationmatrix[16,2]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int16.3), ", ", ""), "01{1,17}0", "0")
  associationmatrix[16,3]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int16.4), ", ", ""), "01{1,17}0", "0")
  associationmatrix[16,4]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int16.5), ", ", ""), "01{1,17}0", "0")
  associationmatrix[16,5]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int16.6), ", ", ""), "01{1,17}0", "0")
  associationmatrix[16,6]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int16.7), ", ", ""), "01{1,17}0", "0")
  associationmatrix[16,7]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int16.8), ", ", ""), "01{1,17}0", "0")
  associationmatrix[16,8]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int16.9), ", ", ""), "01{1,17}0", "0")
  associationmatrix[16,9]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int16.10), ", ", ""), "01{1,17}0", "0")
  associationmatrix[16,10]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int16.11), ", ", ""), "01{1,17}0", "0")
  associationmatrix[16,11]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int16.12), ", ", ""), "01{1,17}0", "0")
  associationmatrix[16,12]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int16.13), ", ", ""), "01{1,17}0", "0")
  associationmatrix[16,13]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int16.14), ", ", ""), "01{1,17}0", "0")
  associationmatrix[16,14]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int16.15), ", ", ""), "01{1,17}0", "0")
  associationmatrix[16,15]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int16.17), ", ", ""), "01{1,17}0", "0")
  associationmatrix[16,17]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int16.18), ", ", ""), "01{1,17}0", "0")
  associationmatrix[16,18]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int16.19), ", ", ""), "01{1,17}0", "0")
  associationmatrix[16,19]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int16.20), ", ", ""), "01{1,17}0", "0")
  associationmatrix[16,20]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int17.1), ", ", ""), "01{1,17}0", "0")
  associationmatrix[17,1]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int17.2), ", ", ""), "01{1,17}0", "0")
  associationmatrix[17,2]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int17.3), ", ", ""), "01{1,17}0", "0")
  associationmatrix[17,3]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int17.4), ", ", ""), "01{1,17}0", "0")
  associationmatrix[17,4]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int17.5), ", ", ""), "01{1,17}0", "0")
  associationmatrix[17,5]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int17.6), ", ", ""), "01{1,17}0", "0")
  associationmatrix[17,6]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int17.7), ", ", ""), "01{1,17}0", "0")
  associationmatrix[17,7]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int17.8), ", ", ""), "01{1,17}0", "0")
  associationmatrix[17,8]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int17.9), ", ", ""), "01{1,17}0", "0")
  associationmatrix[17,9]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int17.10), ", ", ""), "01{1,17}0", "0")
  associationmatrix[17,10]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int17.11), ", ", ""), "01{1,17}0", "0")
  associationmatrix[17,11]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int17.12), ", ", ""), "01{1,17}0", "0")
  associationmatrix[17,12]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int17.13), ", ", ""), "01{1,17}0", "0")
  associationmatrix[17,13]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int17.14), ", ", ""), "01{1,17}0", "0")
  associationmatrix[17,14]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int17.15), ", ", ""), "01{1,17}0", "0")
  associationmatrix[17,15]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int17.16), ", ", ""), "01{1,17}0", "0")
  associationmatrix[17,16]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int17.18), ", ", ""), "01{1,17}0", "0")
  associationmatrix[17,18]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int17.19), ", ", ""), "01{1,17}0", "0")
  associationmatrix[17,19]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int17.20), ", ", ""), "01{1,17}0", "0")
  associationmatrix[17,20]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int18.1), ", ", ""), "01{1,17}0", "0")
  associationmatrix[18,1]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int18.2), ", ", ""), "01{1,17}0", "0")
  associationmatrix[18,2]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int18.3), ", ", ""), "01{1,17}0", "0")
  associationmatrix[18,3]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int18.4), ", ", ""), "01{1,17}0", "0")
  associationmatrix[18,4]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int18.5), ", ", ""), "01{1,17}0", "0")
  associationmatrix[18,5]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int18.6), ", ", ""), "01{1,17}0", "0")
  associationmatrix[18,6]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int18.7), ", ", ""), "01{1,17}0", "0")
  associationmatrix[18,7]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int18.8), ", ", ""), "01{1,17}0", "0")
  associationmatrix[18,8]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int18.9), ", ", ""), "01{1,17}0", "0")
  associationmatrix[18,9]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int18.10), ", ", ""), "01{1,17}0", "0")
  associationmatrix[18,10]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int18.11), ", ", ""), "01{1,17}0", "0")
  associationmatrix[18,11]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int18.12), ", ", ""), "01{1,17}0", "0")
  associationmatrix[18,12]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int18.13), ", ", ""), "01{1,17}0", "0")
  associationmatrix[18,13]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int18.14), ", ", ""), "01{1,17}0", "0")
  associationmatrix[18,14]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int18.15), ", ", ""), "01{1,17}0", "0")
  associationmatrix[18,15]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int18.16), ", ", ""), "01{1,17}0", "0")
  associationmatrix[18,16]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int18.17), ", ", ""), "01{1,17}0", "0")
  associationmatrix[18,17]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int18.19), ", ", ""), "01{1,17}0", "0")
  associationmatrix[18,19]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int18.20), ", ", ""), "01{1,17}0", "0")
  associationmatrix[18,20]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int19.1), ", ", ""), "01{1,17}0", "0")
  associationmatrix[19,1]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int19.2), ", ", ""), "01{1,17}0", "0")
  associationmatrix[19,2]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int19.3), ", ", ""), "01{1,17}0", "0")
  associationmatrix[19,3]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int19.4), ", ", ""), "01{1,17}0", "0")
  associationmatrix[19,4]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int19.5), ", ", ""), "01{1,17}0", "0")
  associationmatrix[19,5]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int19.6), ", ", ""), "01{1,17}0", "0")
  associationmatrix[19,6]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int19.7), ", ", ""), "01{1,17}0", "0")
  associationmatrix[19,7]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int19.8), ", ", ""), "01{1,17}0", "0")
  associationmatrix[19,8]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int19.9), ", ", ""), "01{1,17}0", "0")
  associationmatrix[19,9]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int19.10), ", ", ""), "01{1,17}0", "0")
  associationmatrix[19,10]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int19.11), ", ", ""), "01{1,17}0", "0")
  associationmatrix[19,11]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int19.12), ", ", ""), "01{1,17}0", "0")
  associationmatrix[19,12]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int19.13), ", ", ""), "01{1,17}0", "0")
  associationmatrix[19,13]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int19.14), ", ", ""), "01{1,17}0", "0")
  associationmatrix[19,14]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int19.15), ", ", ""), "01{1,17}0", "0")
  associationmatrix[19,15]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int19.16), ", ", ""), "01{1,17}0", "0")
  associationmatrix[19,16]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int19.17), ", ", ""), "01{1,17}0", "0")
  associationmatrix[19,17]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int19.18), ", ", ""), "01{1,17}0", "0")
  associationmatrix[19,18]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int19.20), ", ", ""), "01{1,17}0", "0")
  associationmatrix[19,20]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int20.1), ", ", ""), "01{1,17}0", "0")
  associationmatrix[20,1]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int20.2), ", ", ""), "01{1,17}0", "0")
  associationmatrix[20,2]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int20.3), ", ", ""), "01{1,17}0", "0")
  associationmatrix[20,3]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int20.4), ", ", ""), "01{1,17}0", "0")
  associationmatrix[20,4]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int20.5), ", ", ""), "01{1,17}0", "0")
  associationmatrix[20,5]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int20.6), ", ", ""), "01{1,17}0", "0")
  associationmatrix[20,6]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int20.7), ", ", ""), "01{1,17}0", "0")
  associationmatrix[20,7]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int20.8), ", ", ""), "01{1,17}0", "0")
  associationmatrix[20,8]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int20.9), ", ", ""), "01{1,17}0", "0")
  associationmatrix[20,9]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int20.10), ", ", ""), "01{1,17}0", "0")
  associationmatrix[20,10]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int20.11), ", ", ""), "01{1,17}0", "0")
  associationmatrix[20,11]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int20.12), ", ", ""), "01{1,17}0", "0")
  associationmatrix[20,12]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int20.13), ", ", ""), "01{1,17}0", "0")
  associationmatrix[20,13]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int20.14), ", ", ""), "01{1,17}0", "0")
  associationmatrix[20,14]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int20.15), ", ", ""), "01{1,17}0", "0")
  associationmatrix[20,15]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int20.16), ", ", ""), "01{1,17}0", "0")
  associationmatrix[20,16]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int20.17), ", ", ""), "01{1,17}0", "0")
  associationmatrix[20,17]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int20.18), ", ", ""), "01{1,17}0", "0")
  associationmatrix[20,18]=sum(str_count(distanceangletimeinteraction, "1"))
  distanceangletimeinteraction <- str_replace_all(str_replace_all(toString(distanceangleinteraction$int20.19), ", ", ""), "01{1,17}0", "0")
  associationmatrix[20,19]=sum(str_count(distanceangletimeinteraction, "1"))
  
  # SAVE ADJACENCY MATRIX AS CSV # ------------------------------------------
  currentmatrix = paste(currentvideo,".csv", sep = "")
  write.table(associationmatrix, file=currentmatrix, row.names=FALSE, col.names=FALSE, sep=",")
  
  # COMPUTE ACTIVITY FOR EACH FLY # -------------
  # DELETE ROWS WHERE pos_x OR pos_y ARE NAs (FOR CALCULATING DISTANCE MOVED FROM FRAME TO FRAME, IF THERE'S MISSING DATA IN ONE ROW, WE CAN STILL INFER HOW FAR THE FLY MOVED OVER THAT UNKNOWN TIME) #
  fly1 <- fly1[which(!is.na(fly1$pos_x)),]
  fly1 <- fly1[which(!is.na(fly1$pos_y)),]
  fly2 <- fly2[which(!is.na(fly2$pos_x)),]
  fly2 <- fly2[which(!is.na(fly2$pos_y)),]
  fly3 <- fly3[which(!is.na(fly3$pos_x)),]
  fly3 <- fly3[which(!is.na(fly3$pos_y)),]
  fly4 <- fly4[which(!is.na(fly4$pos_x)),]
  fly4 <- fly4[which(!is.na(fly4$pos_y)),]
  fly5 <- fly5[which(!is.na(fly5$pos_x)),]
  fly5 <- fly5[which(!is.na(fly5$pos_y)),]
  fly6 <- fly6[which(!is.na(fly6$pos_x)),]
  fly6 <- fly6[which(!is.na(fly6$pos_y)),]
  fly7 <- fly7[which(!is.na(fly7$pos_x)),]
  fly7 <- fly7[which(!is.na(fly7$pos_y)),]
  fly8 <- fly8[which(!is.na(fly8$pos_x)),]
  fly8 <- fly8[which(!is.na(fly8$pos_y)),]
  fly9 <- fly9[which(!is.na(fly9$pos_x)),]
  fly9 <- fly9[which(!is.na(fly9$pos_y)),]
  fly10 <- fly10[which(!is.na(fly10$pos_x)),]
  fly10 <- fly10[which(!is.na(fly10$pos_y)),]
  fly11 <- fly11[which(!is.na(fly11$pos_x)),]
  fly11 <- fly11[which(!is.na(fly11$pos_y)),]
  fly12 <- fly12[which(!is.na(fly12$pos_x)),]
  fly12 <- fly12[which(!is.na(fly12$pos_y)),]
  fly13 <- fly13[which(!is.na(fly13$pos_x)),]
  fly13 <- fly13[which(!is.na(fly13$pos_y)),]
  fly14 <- fly14[which(!is.na(fly14$pos_x)),]
  fly14 <- fly14[which(!is.na(fly14$pos_y)),]
  fly15 <- fly15[which(!is.na(fly15$pos_x)),]
  fly15 <- fly15[which(!is.na(fly15$pos_y)),]
  fly16 <- fly16[which(!is.na(fly16$pos_x)),]
  fly16 <- fly16[which(!is.na(fly16$pos_y)),]
  fly17 <- fly17[which(!is.na(fly17$pos_x)),]
  fly17 <- fly17[which(!is.na(fly17$pos_y)),]
  fly18 <- fly18[which(!is.na(fly18$pos_x)),]
  fly18 <- fly18[which(!is.na(fly18$pos_y)),]
  fly19 <- fly19[which(!is.na(fly19$pos_x)),]
  fly19 <- fly19[which(!is.na(fly19$pos_y)),]
  fly20 <- fly20[which(!is.na(fly20$pos_x)),]
  fly20 <- fly20[which(!is.na(fly20$pos_y)),]
  # CREATE A NEW COLUMN IN THE FLY# DATAFRAMES THAT SHOWS HOW MANY MMs THEY MOVED FROM CURRENT FRAME TO THE NEXT FRAME #
  for (currentFrame in c(1:length(fly1$pos_x))){
    nextFrame = currentFrame +1
    fly1$distmoved[currentFrame] <- sqrt(((fly1$pos_x[currentFrame] - fly1$pos_x[nextFrame])^2) + ((fly1$pos_y[currentFrame] - fly1$pos_y[nextFrame])^2))}
  for (currentFrame in c(1:length(fly2$pos_x))){
    nextFrame = currentFrame +1
    fly2$distmoved[currentFrame] <- sqrt(((fly2$pos_x[currentFrame] - fly2$pos_x[nextFrame])^2) + ((fly2$pos_y[currentFrame] - fly2$pos_y[nextFrame])^2))}
  for (currentFrame in c(1:length(fly3$pos_x))){
    nextFrame = currentFrame +1
    fly3$distmoved[currentFrame] <- sqrt(((fly3$pos_x[currentFrame] - fly3$pos_x[nextFrame])^2) + ((fly3$pos_y[currentFrame] - fly3$pos_y[nextFrame])^2))}
  for (currentFrame in c(1:length(fly4$pos_x))){
    nextFrame = currentFrame +1
    fly4$distmoved[currentFrame] <- sqrt(((fly4$pos_x[currentFrame] - fly4$pos_x[nextFrame])^2) + ((fly4$pos_y[currentFrame] - fly4$pos_y[nextFrame])^2))}
  for (currentFrame in c(1:length(fly5$pos_x))){
    nextFrame = currentFrame +1
    fly5$distmoved[currentFrame] <- sqrt(((fly5$pos_x[currentFrame] - fly5$pos_x[nextFrame])^2) + ((fly5$pos_y[currentFrame] - fly5$pos_y[nextFrame])^2))}
  for (currentFrame in c(1:length(fly6$pos_x))){
    nextFrame = currentFrame +1
    fly6$distmoved[currentFrame] <- sqrt(((fly6$pos_x[currentFrame] - fly6$pos_x[nextFrame])^2) + ((fly6$pos_y[currentFrame] - fly6$pos_y[nextFrame])^2))}
  for (currentFrame in c(1:length(fly7$pos_x))){
    nextFrame = currentFrame +1
    fly7$distmoved[currentFrame] <- sqrt(((fly7$pos_x[currentFrame] - fly7$pos_x[nextFrame])^2) + ((fly7$pos_y[currentFrame] - fly7$pos_y[nextFrame])^2))}
  for (currentFrame in c(1:length(fly8$pos_x))){
    nextFrame = currentFrame +1
    fly8$distmoved[currentFrame] <- sqrt(((fly8$pos_x[currentFrame] - fly8$pos_x[nextFrame])^2) + ((fly8$pos_y[currentFrame] - fly8$pos_y[nextFrame])^2))}
  for (currentFrame in c(1:length(fly9$pos_x))){
    nextFrame = currentFrame +1
    fly9$distmoved[currentFrame] <- sqrt(((fly9$pos_x[currentFrame] - fly9$pos_x[nextFrame])^2) + ((fly9$pos_y[currentFrame] - fly9$pos_y[nextFrame])^2))}
  for (currentFrame in c(1:length(fly10$pos_x))){
    nextFrame = currentFrame +1
    fly10$distmoved[currentFrame] <- sqrt(((fly10$pos_x[currentFrame] - fly10$pos_x[nextFrame])^2) + ((fly10$pos_y[currentFrame] - fly10$pos_y[nextFrame])^2))}
  for (currentFrame in c(1:length(fly11$pos_x))){
    nextFrame = currentFrame +1
    fly11$distmoved[currentFrame] <- sqrt(((fly11$pos_x[currentFrame] - fly11$pos_x[nextFrame])^2) + ((fly11$pos_y[currentFrame] - fly11$pos_y[nextFrame])^2))}
  for (currentFrame in c(1:length(fly12$pos_x))){
    nextFrame = currentFrame +1
    fly12$distmoved[currentFrame] <- sqrt(((fly12$pos_x[currentFrame] - fly12$pos_x[nextFrame])^2) + ((fly12$pos_y[currentFrame] - fly12$pos_y[nextFrame])^2))}
  for (currentFrame in c(1:length(fly13$pos_x))){
    nextFrame = currentFrame +1
    fly13$distmoved[currentFrame] <- sqrt(((fly13$pos_x[currentFrame] - fly13$pos_x[nextFrame])^2) + ((fly13$pos_y[currentFrame] - fly13$pos_y[nextFrame])^2))}
  for (currentFrame in c(1:length(fly14$pos_x))){
    nextFrame = currentFrame +1
    fly14$distmoved[currentFrame] <- sqrt(((fly14$pos_x[currentFrame] - fly14$pos_x[nextFrame])^2) + ((fly14$pos_y[currentFrame] - fly14$pos_y[nextFrame])^2))}
  for (currentFrame in c(1:length(fly15$pos_x))){
    nextFrame = currentFrame +1
    fly15$distmoved[currentFrame] <- sqrt(((fly15$pos_x[currentFrame] - fly15$pos_x[nextFrame])^2) + ((fly15$pos_y[currentFrame] - fly15$pos_y[nextFrame])^2))}
  for (currentFrame in c(1:length(fly16$pos_x))){
    nextFrame = currentFrame +1
    fly16$distmoved[currentFrame] <- sqrt(((fly16$pos_x[currentFrame] - fly16$pos_x[nextFrame])^2) + ((fly16$pos_y[currentFrame] - fly16$pos_y[nextFrame])^2))}
  for (currentFrame in c(1:length(fly17$pos_x))){
    nextFrame = currentFrame +1
    fly17$distmoved[currentFrame] <- sqrt(((fly17$pos_x[currentFrame] - fly17$pos_x[nextFrame])^2) + ((fly17$pos_y[currentFrame] - fly17$pos_y[nextFrame])^2))}
  for (currentFrame in c(1:length(fly18$pos_x))){
    nextFrame = currentFrame +1
    fly18$distmoved[currentFrame] <- sqrt(((fly18$pos_x[currentFrame] - fly18$pos_x[nextFrame])^2) + ((fly18$pos_y[currentFrame] - fly18$pos_y[nextFrame])^2))}
  for (currentFrame in c(1:length(fly19$pos_x))){
    nextFrame = currentFrame +1
    fly19$distmoved[currentFrame] <- sqrt(((fly19$pos_x[currentFrame] - fly19$pos_x[nextFrame])^2) + ((fly19$pos_y[currentFrame] - fly19$pos_y[nextFrame])^2))}
  for (currentFrame in c(1:length(fly20$pos_x))){
    nextFrame = currentFrame +1
    fly20$distmoved[currentFrame] <- sqrt(((fly20$pos_x[currentFrame] - fly20$pos_x[nextFrame])^2) + ((fly20$pos_y[currentFrame] - fly20$pos_y[nextFrame])^2))}
  
  # COMPUTE FLY ACTIVITIES #
  currentflyactivities <- data.frame(matrix(ncol = 1, nrow = 20))
  colnames(currentflyactivities) <- c("activity")
  currentflyactivities$activity[1] <- sum(fly1$distmoved, na.rm = T)
  currentflyactivities$activity[2] <- sum(fly2$distmoved, na.rm = T)
  currentflyactivities$activity[3] <- sum(fly3$distmoved, na.rm = T)
  currentflyactivities$activity[4] <- sum(fly4$distmoved, na.rm = T)
  currentflyactivities$activity[5] <- sum(fly5$distmoved, na.rm = T)
  currentflyactivities$activity[6] <- sum(fly6$distmoved, na.rm = T)
  currentflyactivities$activity[7] <- sum(fly7$distmoved, na.rm = T)
  currentflyactivities$activity[8] <- sum(fly8$distmoved, na.rm = T)
  currentflyactivities$activity[9] <- sum(fly9$distmoved, na.rm = T)
  currentflyactivities$activity[10] <- sum(fly10$distmoved, na.rm = T)
  currentflyactivities$activity[11] <- sum(fly11$distmoved, na.rm = T)
  currentflyactivities$activity[12] <- sum(fly12$distmoved, na.rm = T)
  currentflyactivities$activity[13] <- sum(fly13$distmoved, na.rm = T)
  currentflyactivities$activity[14] <- sum(fly14$distmoved, na.rm = T)
  currentflyactivities$activity[15] <- sum(fly15$distmoved, na.rm = T)
  currentflyactivities$activity[16] <- sum(fly16$distmoved, na.rm = T)
  currentflyactivities$activity[17] <- sum(fly17$distmoved, na.rm = T)
  currentflyactivities$activity[18] <- sum(fly18$distmoved, na.rm = T)
  currentflyactivities$activity[19] <- sum(fly19$distmoved, na.rm = T)
  currentflyactivities$activity[20] <- sum(fly20$distmoved, na.rm = T)
  currentflyactivities$fly <- 1:nrow(currentflyactivities)
  currentflyactivities$video <- currentvideo
  # MERGE ACTIVITY DATAFRAME #
  activity <- rbind(activity, currentflyactivities)
  rm(currentflyactivities)
}
rm(anglearcsin,angledifference,angleinteraction,angles,associationmatrix,distanceangleinteraction,distanceinteraction,distances,fly1,fly2,fly3,fly4,fly5,fly6,fly7,fly8,fly9,fly10,fly11,fly12,fly13,fly14,fly15,fly16,fly17,fly18,fly19,fly20,currentFrame,currentmatrix,currentvideo,currentvideodata,currentvideopxmm,distanceangletimeinteraction,headers,meansize,nextFrame)

# PULL OUT NODE AND NETWORK METRICS FROM ADJACENCY MATRICES #
nodemetrics <- data.frame(matrix(ncol = 1, nrow = 0))
networkmetrics <- data.frame(matrix(ncol = 1, nrow = 0))
for (currentvideo in videos) {
  # LOAD IN ADJACENCY MATRIX #
  currentmatrix = paste(currentvideo, ".csv", sep = "")
  associationmatrix <- read_csv(currentmatrix, col_names = FALSE)
  associationmatrix <- as.matrix(associationmatrix)
  colnames(associationmatrix) <- flies
  rownames(associationmatrix) <- flies
  associationmatrix <- graph.adjacency(associationmatrix, mode = "directed", weighted = TRUE, diag = FALSE)
  # CALCULATE NODE METRICS #
  currentnodemetrics <- strength(associationmatrix, mode=c("in"))
  currentnodemetrics <- data.frame(currentnodemetrics)
  colnames(currentnodemetrics) <- c("indegree")
  currentnodemetrics$outdegree <- strength(associationmatrix, mode=c("out"))
  eigencentrality <- eigen_centrality(associationmatrix, directed=TRUE, scale = FALSE, weights=NULL)
  currentnodemetrics$eigencentrality <- eigencentrality$vector
  rm(eigencentrality)
  currentnodemetrics$betweenness <- betweenness(associationmatrix, nobigint = FALSE)
  currentnodemetrics$clusteringcoeff <- transitivity(associationmatrix, type=c("weighted"))
  setDT(currentnodemetrics, keep.rownames = TRUE)[]
  colnames(currentnodemetrics)[1] <- c("flyID")
  currentnodemetrics$video <- currentvideo
  nodemetrics <- rbind(nodemetrics, currentnodemetrics, fill = TRUE)
  # CALCULATE NETWORK METRICS #
  currentnetworkmetrics <- edge_density(associationmatrix, loops = FALSE)
  currentnetworkmetrics <- data.frame(currentnetworkmetrics)
  colnames(currentnetworkmetrics) <- c("density")
  matrixcluster <- cluster_edge_betweenness(associationmatrix)
  currentnetworkmetrics$modularity <- modularity(associationmatrix, membership(matrixcluster))
  currentnetworkmetrics$transitivity <- transitivity(associationmatrix, type=c("global"))
  currentnetworkmetrics$diameter <- diameter(associationmatrix, directed = TRUE)
  currentnetworkmetrics$video <- currentvideo
  networkmetrics <- rbind(networkmetrics, currentnetworkmetrics)}
rm(currentnetworkmetrics,currentnodemetrics,associationmatrix,currentmatrix,currentvideo,flies,matrixcluster,videos)
nodemetrics$matrix.ncol...1..nrow...0. <- NULL

# LOAD IN INDIVIDUAL FLY DATA, TRIAL DATA, MATING DATA, AND LIFETIME FITNESS DATA # --------------
trialdata <- read_csv("trialdata.csv", 
                      col_types = cols(trial = col_character(),
                                       food = col_character(),
                                       genotype = col_character(),
                                       sex = col_character(),
                                       color = col_character(),
                                       emerged = col_date(format = "%m/%d/%y"), 
                                       lightslefton = col_number(), 
                                       lifestatus = col_number()))
matings <- read_csv("firstmatings.csv", 
                    col_types = cols(trial = col_character(), 
                                     male.color = col_character(), 
                                     female.color = col_character(), 
                                     setuptime = col_time(format = "%H:%M"), 
                                     matingstart = col_time(format = "%H:%M"), 
                                     matingend = col_time(format = "%H:%M"), 
                                     setupdate = col_date(format = "%m/%d/%y")))
fitness <- read_csv("fitness.csv", 
                    col_types = cols(trial = col_character(), 
                                     genotype = col_character(), 
                                     emerged = col_date(format = "%m/%d/%y"), 
                                     deathdate = col_date(format = "%m/%d/%y"), 
                                     causeofdeath = col_character(), 
                                     prematuremale = col_number(), 
                                     d0c1 = col_number(), 
                                     d0c2 = col_number(), 
                                     d7c1 = col_number(), 
                                     d7c2 = col_number(), 
                                     d14c1 = col_number(), 
                                     d14c2 = col_number(), 
                                     d21c1 = col_number(), 
                                     d21c2 = col_number(), 
                                     d28c1 = col_number(), 
                                     d28c2 = col_number(), 
                                     d35c1 = col_number(), 
                                     d35c2 = col_number(), 
                                     d42c1 = col_number(), 
                                     d42c2 = col_number(), 
                                     d49c1 = col_number(), 
                                     d49c2 = col_number(), 
                                     d56c1 = col_number(), 
                                     d56c2 = col_number(), 
                                     d63c1 = col_number(), 
                                     d63c2 = col_number(), 
                                     d70c1 = col_number(), 
                                     d70c2 = col_number(), 
                                     d77c1 = col_number(), 
                                     d77c2 = col_number(), 
                                     d84c1 = col_number(), 
                                     d84c2 = col_number(), 
                                     d91c1 = col_number(), 
                                     d91c2 = col_number()))

# TURN FLY COLORS AND FLY IDS INTO A MERGEABLE DATASET # ----------------------
videoflycolors <- data.frame(videoinformation$video)
colnames(videoflycolors) <- c("video")
videoflycolors$femaleNR <- videoinformation$femaleNR
videoflycolors$femaleCO <- videoinformation$femaleCO
videoflycolors$femaleCY <- videoinformation$femaleCY
videoflycolors$femaleLG <- videoinformation$femaleLG
videoflycolors$femaleCG <- videoinformation$femaleCG
videoflycolors$femaleUB <- videoinformation$femaleUB
videoflycolors$femaleBB <- videoinformation$femaleBB
videoflycolors$femaleBP <- videoinformation$femaleBP
videoflycolors$femalePP <- videoinformation$femalePP
videoflycolors$femaleTW <- videoinformation$femaleTW
videoflycolors$maleNR <- videoinformation$maleNR
videoflycolors$maleCO <- videoinformation$maleCO
videoflycolors$maleCY <- videoinformation$maleCY
videoflycolors$maleLG <- videoinformation$maleLG
videoflycolors$maleCG <- videoinformation$maleCG
videoflycolors$maleUB <- videoinformation$maleUB
videoflycolors$maleBB <- videoinformation$maleBB
videoflycolors$maleBP <- videoinformation$maleBP
videoflycolors$malePP <- videoinformation$malePP
videoflycolors$maleTW <- videoinformation$maleTW
videoflycolors <- melt(videoflycolors)
videoflycolors$video <- as.character(videoflycolors$video)
videoflycolors$sex <- str_split_fixed(videoflycolors$variable, "ale", 2)
videoflycolors$color <- str_sub(videoflycolors$variable, -2, -1)
videoflycolors$sex <- ifelse(videoflycolors$sex == "fem", "F", "M")
videoflycolors$variable <- NULL
colnames(videoflycolors) <- c("video", "flyID", "sex", "color")
videoinformation$femaleNR <- NULL
videoinformation$femaleCO <- NULL
videoinformation$femaleCY <- NULL
videoinformation$femaleLG <- NULL
videoinformation$femaleCG <- NULL
videoinformation$femaleUB <- NULL
videoinformation$femaleBB <- NULL
videoinformation$femaleBP <- NULL
videoinformation$femalePP <- NULL
videoinformation$femaleTW <- NULL
videoinformation$maleNR <- NULL
videoinformation$maleCO <- NULL
videoinformation$maleCY <- NULL
videoinformation$maleLG <- NULL
videoinformation$maleCG <- NULL
videoinformation$maleUB <- NULL
videoinformation$maleBB <- NULL
videoinformation$maleBP <- NULL
videoinformation$malePP <- NULL
videoinformation$maleTW <- NULL

# MERGE FLY COLOR INFORMATION, NODE METRICS, ACTIVITY, VIDEO INFORMATION, AND TRIAL INFORMATION # ---------------
nodemetrics$flyID <- sub("fly", "", nodemetrics$flyID)
nodemetrics$flyID <- as.numeric(nodemetrics$flyID)
activity$flyID <- activity$fly
activity$fly <- NULL
flydata <- merge(nodemetrics, activity, by=c("video", "flyID"))
rm(nodemetrics, activity)
flydata <- merge(flydata, videoflycolors, by=c("video", "flyID"))
rm(videoflycolors)
flydata <- merge(flydata, videoinformation, by=c("video"))
rm(videoinformation)
flydata <- merge(flydata, trialdata, by=c("trial", "sex", "color"), all = TRUE)
flydata$time <- ifelse(flydata$time == "NA", NA, flydata$time)
rm(trialdata)

# MERGE LIFETIME FITNESS DATA AND MATING DATA INTO THE FLY DATASET # --------
flydata <- merge(flydata, fitness, by=c("trial", "genotype"), all=TRUE)
rm(fitness)
flydata$lifespan <- as.numeric(flydata$deathdate - flydata$emerged.x)
flydata$lifespan <- ifelse(flydata$causeofdeath == "Escape", NA, flydata$lifespan)
flydata$lifespan <- ifelse(flydata$causeofdeath == "Cold Knockout", NA, flydata$lifespan)
flydata$lifespan <- ifelse(flydata$causeofdeath == "NA", NA, flydata$lifespan)
flydata$lifespan <- ifelse(flydata$prematuremale == 1, NA, flydata$lifespan)
flydata$offspring <- rowSums(flydata[,c(27:54)], na.rm = TRUE)
flydata$offspring <- ifelse(flydata$causeofdeath == "Escape", NA, flydata$offspring)
flydata$offspring <- ifelse(flydata$causeofdeath == "Cold Knockout", NA, flydata$offspring)
flydata$offspring <- ifelse(flydata$causeofdeath == "NA", NA, flydata$offspring)
flydata$offspring <- ifelse(flydata$prematuremale == 1, NA, flydata$offspring)
flydata$reprorate <- flydata$offspring / flydata$lifespan

# MERGING MATING DATA INTO THE FLY DATASET # --------------
malematings <- matings[c(1:2,4:7)]
malematings$sex <- "M"
names(malematings)[2] <- "color"
femalematings <- matings[c(1,3:7)]
femalematings$sex <- "F"
names(femalematings)[2] <- "color"
matedata <- rbind(malematings, femalematings)
rm(malematings,femalematings)
matedata$latency <- as.numeric(difftime(matedata$matingstart,matedata$setuptime))
matenumber <- data.frame(table(matedata$trial,matedata$color,matedata$sex))
colnames(matenumber) <- c("trial", "color", "sex", "matings")
latency <- matedata %>% 
  group_by(trial,color,sex) %>% 
  slice(which.min(latency))
matedata <- merge(matenumber, latency, by = c("trial","color","sex"), all=TRUE)
matedata$matingstart <- NULL
matedata$matingend <- NULL
rm(matenumber,latency)
datetime <- unique(matedata[,c('trial','setuptime','setupdate')])
datetime <- datetime[complete.cases(datetime), ]
matedata <- merge(matedata, datetime, by = c("trial"))
rm(datetime)
matedata$setuptime.x <- NULL
matedata$setupdate.x <- NULL
names(matedata)[6] <- "setuptime"
names(matedata)[7] <- "setupdate"
flydata <- merge(flydata, matedata, by=c("trial", "color", "sex"), all=TRUE)
rm(matedata, matings)

# MERGE NETWORK-WIDE DATA INTO FLY DATASET # ---------------
flydata <- merge(flydata, networkmetrics, by=c("video"), all = TRUE)
rm(networkmetrics)

# CLEAN UP DATA # -------------------
cleanflydata <- data.frame(flydata$trial)
colnames(cleanflydata) <- c("trial")
cleanflydata$trial <- as.character(cleanflydata$trial)
cleanflydata$food <- as.character(flydata$food)
cleanflydata$pc <- as.numeric(ifelse(flydata$food == 3, 1, (ifelse(flydata$food == 2, 2, 4))))
cleanflydata$pc.rescaled <- as.numeric(rescale(cleanflydata$pc, to=c(0,1)))
cleanflydata$dilution <- as.numeric(ifelse(flydata$food == 5, 1, (ifelse(flydata$food == 4, 2, 4))))
cleanflydata$dilution.rescaled <- as.numeric(rescale(cleanflydata$dilution, to=c(0,1)))
cleanflydata$sex <- as.character(flydata$sex)
cleanflydata$genotype <- as.character(flydata$genotype)
cleanflydata$color <- as.character(flydata$color)
cleanflydata$emerged <- as.Date(flydata$emerged.x, "%Y-%M-%D")
cleanflydata$setuptime <- as.hms(flydata$setuptime)
cleanflydata$setuptime.rescaled <- as.numeric(cleanflydata$setuptime)
cleanflydata$setuptime.rescaled <- as.numeric(rescale(cleanflydata$setuptime.rescaled, to=c(0,1)))
cleanflydata$latency <- as.numeric(flydata$latency)
cleanflydata$matings <- as.numeric(flydata$matings)
cleanflydata$lifestatus <- as.numeric(flydata$lifestatus)
cleanflydata$video <- as.character(flydata$video)
cleanflydata$videodate <- as.Date(flydata$videodate, "%Y-%M-%D")
cleanflydata$videoday <- as.numeric(flydata$day)
cleanflydata$videoday.rescaled <- as.numeric(rescale(cleanflydata$videoday, to=c(0,1)))
cleanflydata$videotime <- as.hms(flydata$time)
cleanflydata$videotime.rescaled <- as.numeric(cleanflydata$videotime)
cleanflydata$videotime.rescaled <- as.numeric(rescale(cleanflydata$videotime.rescaled, to=c(0,1)))
cleanflydata$activity <- as.numeric(flydata$activity)
cleanflydata$activity.rescaled <- as.numeric(rescale(cleanflydata$activity, to=c(0,1)))
cleanflydata$indegree <- as.numeric(flydata$indegree)
cleanflydata$indegree.rescaled <- as.numeric(rescale(cleanflydata$indegree, to=c(0,1)))
cleanflydata$outdegree <- as.numeric(flydata$outdegree)
cleanflydata$outdegree.rescaled <- as.numeric(rescale(cleanflydata$outdegree, to=c(0,1)))
cleanflydata$eigencentrality <- as.numeric(flydata$eigencentrality)
cleanflydata$eigencentrality.rescaled <- as.numeric(rescale(cleanflydata$eigencentrality, to=c(0,1)))
cleanflydata$betweenness <- as.numeric(flydata$betweenness)
cleanflydata$betweenness.rescaled <- as.numeric(rescale(cleanflydata$betweenness, to=c(0,1)))
cleanflydata$clusteringcoeff <- as.numeric(flydata$clusteringcoeff)
cleanflydata$clusteringcoeff.rescaled <- as.numeric(rescale(cleanflydata$clusteringcoeff, to=c(0,1)))
cleanflydata$density <- as.numeric(flydata$density)
cleanflydata$modularity <- as.numeric(flydata$modularity)
cleanflydata$transitivity <- as.numeric(flydata$transitivity)
cleanflydata$diameter <- as.numeric(flydata$diameter)
cleanflydata$offspring <- as.numeric(flydata$offspring)
cleanflydata$deathdate <- as.Date(flydata$deathdate, "%Y-%M-%D")
cleanflydata$causeofdeath <- as.character(flydata$causeofdeath)
cleanflydata$lifespan <- as.numeric(flydata$lifespan)
cleanflydata$reprorate <- as.numeric(flydata$reprorate)
flydata <- cleanflydata[order(cleanflydata$trial, cleanflydata$sex, cleanflydata$genotype),]
rm(cleanflydata)

# EFFECTS ON NODE-LEVEL NETWORK POSITION # --------------------
# Reduce dataset to only cases where network position was measured.
network  <- flydata[complete.cases(flydata[, 16]),]
# Make individual and observation effect columns
network$individual <- paste(network$trial, network$genotype, sep = "_")
network$observation <- paste(network$trial, network$genotype, network$videoday, sep = "_")
# Make random effects factors to work with glmmadmb (doesn't change output of lme4 models)
network$trial <- as.factor(network$trial)
network$food <- as.factor(network$food)
network$genotype <- as.factor(network$genotype)
network$individual <- as.factor(network$individual)
network$observation <- as.factor(network$observation)

# OBSERVED SEX*E EFFECTS ON NETWORK POSITIONS AND ACTIVITY # ----------------
# Model fit testing
# Instrength
indegree_G <- glmer(indegree ~  (1|genotype) + sex*food + (1|trial), data = network, family = poisson)
indegree_G_sim <- simulateResiduals(fittedModel = indegree_G, n = 1000)
testResiduals(indegree_G_sim) # Overdispersion
indegree_G <- glmer(indegree ~  (1|genotype) + sex*food + (1|trial) + (1|observation), data = network, family = poisson)
summary(indegree_G)
# Outstrength
outdegree_G <- glmer(outdegree ~  (1|genotype) + sex*food + (1|trial), data = network, family = poisson)
outdegree_G_sim <- simulateResiduals(fittedModel = outdegree_G, n = 1000)
testResiduals(outdegree_G_sim) # Overdispersion
outdegree_G <- glmer(outdegree ~  (1|genotype) + sex*food + (1|trial) + (1|observation), data = network, family = poisson)
summary(outdegree_G)
# Eigenvector Centrality
eigencentrality_G <- lmer(eigencentrality ~  (1|genotype) + sex*food + (1|individual), data = network)
eigencentrality_G_sim <- simulateResiduals(fittedModel = eigencentrality_G, n = 1000)
testResiduals(eigencentrality_G_sim) # Adequate model fit
summary(eigencentrality_G)
# Betweenness Centrality
betweenness_G <- glmer(betweenness ~  (1|genotype) + sex*food + (1|trial), data = network, family = poisson)
betweenness_G_sim <- simulateResiduals(fittedModel = betweenness_G, n = 1000)
testResiduals(betweenness_G_sim) # 0 inflation
betweenness_G <- glmmTMB(betweenness ~  (1|genotype) + sex*food + (1|trial), data = network, family = "poisson", ziformula = ~1)
summary(betweenness_G)
# Clustering Coefficient
clusteringcoeff_G <- lmer(clusteringcoeff ~  (1|genotype) + sex*food + (1|trial), data = network)
clusteringcoeff_G_sim <- simulateResiduals(fittedModel = clusteringcoeff_G, n = 1000)
testResiduals(clusteringcoeff_G_sim) # Adequate model fit
summary(clusteringcoeff_G)
# Activity
activity_G <- lmer(activity ~  (1|genotype) + sex*food + (1|trial), data = network)
activity_G_sim <- simulateResiduals(fittedModel = activity_G, n = 1000)
testResiduals(activity_G_sim) # Adequate model fit
summary(activity_G)
# Observed significance of Sex*E
Anova(indegree_G, type = "3")          # Observed significance of Sex*E:   Chisq = 4.1841,  Df = 4,  p = 0.38166
Anova(outdegree_G, type = "3")         # Observed significance of Sex*E:   Chisq = 13.6913, Df = 4,  p = 0.008348
Anova(eigencentrality_G, type = "3")   # Observed significance of Sex*E:   Chisq = 4.7778,  Df = 4,  p = 0.31087
Anova(betweenness_G, type = "3")       # Observed significance of Sex*E:   Chisq = 55.9562, Df = 4,  p = 2.048e-11
Anova(clusteringcoeff_G, type = "3")   # Observed significance of Sex*E:   LR = 38.1301,    Df = 4,  p = 1.053e-07
Anova(activity_G, type = "3")          # Observed significance of Sex*E:   LR = 12.6146,    Df = 4,  p = 0.01332
# Save observed test statistics as values
indegree_SexxE_obs        <- Anova(indegree_G, type = "3")$Chisq[4]
outdegree_SexxE_obs       <- Anova(outdegree_G, type = "3")$Chisq[4]
eigencentrality_SexxE_obs <- Anova(eigencentrality_G, type = "3")$Chisq[4]
betweenness_SexxE_obs     <- Anova(betweenness_G, type = "3")$Chisq[4]
clusteringcoeff_SexxE_obs <- Anova(clusteringcoeff_G, type = "3")$Chisq[4]
activity_SexxE_obs        <- Anova(activity_G, type = "3")$Chisq[4]

# PERMUTED SEX*E EFFECTS ON NETWORK POSITIONS AND ACTIVITY # ----------------
indegree_SexxE_null        = numeric(1000)
outdegree_SexxE_null       = numeric(1000)
eigencentrality_SexxE_null = numeric(1000)
betweenness_SexxE_null     = numeric(1000)
clusteringcoeff_SexxE_null = numeric(1000)
activity_SexxE_null        = numeric(1000)
for (i in 1:1000) {
  network_rand = network %>% 
    group_by(video) %>%
    mutate(sex_rand = sample(sex, replace=F))
  indegree_SexxE_rand        <- glmer(indegree       ~ (1|genotype) + sex_rand*food + (1|trial) + (1|observation), data = network_rand, family = poisson)
  outdegree_SexxE_rand       <- glmer(outdegree      ~ (1|genotype) + sex_rand*food + (1|trial) + (1|observation), data = network_rand, family = poisson)
  eigencentrality_SexxE_rand <- lmer(eigencentrality ~ (1|genotype) + sex_rand*food + (1|individual), data = network_rand)
  betweenness_SexxE_rand     <- glmmTMB(betweenness  ~ (1|genotype) + sex_rand*food + (1|trial), data = network_rand, family = "poisson", ziformula = ~1)
  clusteringcoeff_SexxE_rand <- lmer(clusteringcoeff ~ (1|genotype) + sex_rand*food + (1|trial), data = network_rand)
  activity_SexxE_rand        <- lmer(activity        ~ (1|genotype) + sex_rand*food + (1|trial), data = network_rand)
  indegree_SexxE_null[i]        <- Anova(indegree_SexxE_rand,  type = "3")$Chisq[4]
  outdegree_SexxE_null[i]       <- Anova(outdegree_SexxE_rand, type = "3")$Chisq[4]
  eigencentrality_SexxE_null[i] <- Anova(eigencentrality_SexxE_rand, type = "3")$Chisq[4]
  betweenness_SexxE_null[i]     <- Anova(betweenness_SexxE_rand, type = "3")$Chisq[4]
  clusteringcoeff_SexxE_null[i] <- Anova(clusteringcoeff_SexxE_rand, type = "3")$Chisq[4]
  activity_SexxE_null[i] <- Anova(activity_SexxE_rand, type = "3")$Chisq[4]}
sum(indegree_SexxE_obs<indegree_SexxE_null)/1000               # Permuted significance of Sex*E: p = 0.287
sum(outdegree_SexxE_obs<outdegree_SexxE_null)/1000             # Permuted significance of Sex*E: p = 0.002
sum(eigencentrality_SexxE_obs<eigencentrality_SexxE_null)/1000 # Permuted significance of Sex*E: p = 0.319
sum(betweenness_SexxE_obs<betweenness_SexxE_null)/1000         # Permuted significance of Sex*E: p = 0.885
sum(clusteringcoeff_SexxE_obs<clusteringcoeff_SexxE_null)/1000 # Permuted significance of Sex*E: p < 0.001
sum(activity_SexxE_obs<activity_SexxE_null)/1000               # Permuted significance of Sex*E: p = 0.003

# OBSERVED G*E, G, E, AND SEX EFFECTS ON NETWORK POSITIONS AND ACTIVITY # ----------------
indegree_GxE        <- glmer(indegree       ~ (food|genotype) + food + sex + (1|trial) + (1|observation), data = network, family = poisson)
indegree_G          <- glmer(indegree       ~ (1|genotype)    + food + sex + (1|trial) + (1|observation), data = network, family = poisson)
indegree_NA         <- glmer(indegree       ~                   food + sex + (1|trial) + (1|observation), data = network, family = poisson)
outdegree_GxE       <- glmer(outdegree      ~ (food|genotype) + sex*food   + (1|trial) + (1|observation), data = network, family = poisson)
outdegree_G         <- glmer(outdegree      ~ (1|genotype)    + sex*food   + (1|trial) + (1|observation), data = network, family = poisson)
outdegree_NA        <- glmer(outdegree      ~                   sex*food   + (1|trial) + (1|observation), data = network, family = poisson)
eigencentrality_GxE <- lmer(eigencentrality ~ (food|genotype) + sex        + (1|individual),              data = network)
eigencentrality_G   <- lmer(eigencentrality ~ (1|genotype)    + sex        + (1|individual),              data = network)
eigencentrality_NA  <- lmer(eigencentrality ~                   sex        + (1|individual),              data = network)
betweenness_GxE     <- glmmTMB(betweenness  ~ (food|genotype) + food + sex + (1|trial),                   data = network, family = "poisson", ziformula = ~1)
betweenness_G       <- glmmTMB(betweenness  ~ (1|genotype)    + food + sex + (1|trial),                   data = network, family = "poisson", ziformula = ~1)
betweenness_NA      <- glmmTMB(betweenness  ~                   food + sex + (1|trial),                   data = network, family = "poisson", ziformula = ~1)
clusteringcoeff_GxE <- lmer(clusteringcoeff ~ (food|genotype) + sex*food   + (1|trial),                   data = network)
clusteringcoeff_G   <- lmer(clusteringcoeff ~ (1|genotype)    + sex*food   + (1|trial),                   data = network)
clusteringcoeff_NA  <- lmer(clusteringcoeff ~                   sex*food   + (1|trial),                   data = network)
activity_GxE        <- lmer(activity        ~ (food|genotype) + sex*food   + (1|trial),                   data = network)
activity_G          <- lmer(activity        ~ (1|genotype)    + sex*food   + (1|trial),                   data = network)
activity_NA         <- lmer(activity        ~                   sex*food   + (1|trial),                   data = network)
# Observed significance of GxE
anova(indegree_GxE, indegree_G)                         # Observed significance of GxE:            LR = 4.284,       Df = 14, p = 0.9935
anova(outdegree_GxE, outdegree_G)                       # Observed significance of GxE:            LR = 1.0264,      Df = 14, p = 1
anova(eigencentrality_GxE, eigencentrality_G)           # Observed significance of GxE:            LR = 0,           Df = 14, p = 1
anova(betweenness_GxE, betweenness_G)                   # Observed significance of GxE:            LR = 2133.3,      Df = 14, p < 2.2e-16
anova(clusteringcoeff_GxE, clusteringcoeff_G)           # Observed significance of GxE:            LR = 5.0002,      Df = 14, p = 0.9858
anova(activity_GxE, activity_G)                         # Observed significance of GxE:            LR = 5.3016,      Df = 14, p = 0.9812
# Observed significance of G
anova(indegree_G, indegree_NA)                          # Observed significance of G:              LR = 22.713,      Df = 1,  p = 1.881e-06
anova(outdegree_G, outdegree_NA)                        # Observed significance of G:              LR = 24.13,       Df = 1,  p = 9.005e-07
anova(eigencentrality_G, eigencentrality_NA)            # Observed significance of G:              LR = 24.348,      Df = 1,  p = 8.039e-07
anova(betweenness_G, betweenness_NA)                    # Observed significance of G:              LR = 1086,        Df = 1,  p < 2.2e-16
anova(clusteringcoeff_G, clusteringcoeff_NA)            # Observed significance of G:              LR = 54.726,      Df = 1,  p = 1.385e-13
anova(activity_G, activity_NA)                          # Observed significance of G:              LR = 118.87,      Df = 1,  p < 2.2e-16
# Observed significance of E
Anova(indegree_G, type = "3")                           # Observed significance of E:              Chisq = 6.306,    Df = 4,  p = 0.177430
Anova(outdegree_G, type = "3")                          # Observed significance of E:              Chisq = 9.4918,   Df = 4,  p = 0.049916
fligner.test(eigencentrality ~ food, data = network)    # Observed significance of E on variation: Chisq = 8.3337,   Df = 4,  p = 0.08009
Anova(betweenness_G, type = "3")                        # Observed significance of E:              Chisq = 2.2667,   Df = 4,  p = 0.68684
Anova(clusteringcoeff_G, type = "3")                    # Observed significance of E:              Chisq = 9.8733,   Df = 4,  p = 0.04262
Anova(activity_G, type = "3")                           # Observed significance of E:              Chisq = 4.5636,   Df = 4,  p = 0.33508
# Observed significance of Sex
Anova(indegree_G, type = "3")                           # Observed significance of Sex:            Chisq = 8.072,    Df = 1,  p = 0.004495
Anova(outdegree_G, type = "3")                          # Observed significance of Sex:            Chisq = 10.1084,  Df = 1,  p = 0.001476
Anova(eigencentrality_G, type = "3")                    # Observed significance of Sex:            Chisq = 9.2482,   Df = 1,  p = 0.002357
Anova(betweenness_G, type = "3")                        # Observed significance of Sex:            Chisq = 6.2083,   Df = 1,  p = 0.01272
Anova(clusteringcoeff_G, type = "3")                    # Observed significance of Sex:            Chisq = 106.8320, Df = 1,  p < 2e-16
Anova(activity_G, type = "3")                           # Observed significance of Sex:            Chisq = 15.6791,  Df = 1,  p = 7.505e-05
# Save observed test statistics as values
indegree_GxE_obs        <- anova(indegree_GxE, indegree_G)$Chisq[2]
outdegree_GxE_obs       <- anova(outdegree_GxE, outdegree_G)$Chisq[2]
eigencentrality_GxE_obs <- anova(eigencentrality_GxE, eigencentrality_G)$Chisq[2]
betweenness_GxE_obs     <- anova(betweenness_GxE, betweenness_G)$Chisq[2]
clusteringcoeff_GxE_obs <- anova(clusteringcoeff_GxE, clusteringcoeff_G)$Chisq[2]
activity_GxE_obs        <- anova(activity_GxE, activity_G)$Chisq[2]
indegree_G_obs          <- anova(indegree_G, indegree_NA)$Chisq[2]
outdegree_G_obs         <- anova(outdegree_G, outdegree_NA)$Chisq[2]
eigencentrality_G_obs   <- anova(eigencentrality_G, eigencentrality_NA)$Chisq[2]
betweenness_G_obs       <- anova(betweenness_G, betweenness_NA)$Chisq[2]
clusteringcoeff_G_obs   <- anova(clusteringcoeff_G, clusteringcoeff_NA)$Chisq[2]
activity_G_obs          <- anova(activity_G, activity_NA)$Chisq[2]
indegree_Sex_obs        <- Anova(indegree_G, type = "3")$Chisq[3]
outdegree_Sex_obs       <- Anova(outdegree_G, type = "3")$Chisq[2]
eigencentrality_Sex_obs <- Anova(eigencentrality_G, type = "3")$Chisq[2]
betweenness_Sex_obs     <- Anova(betweenness_G, type = "3")$Chisq[3]
clusteringcoeff_Sex_obs <- Anova(clusteringcoeff_G, type = "3")$Chisq[2]
activity_Sex_obs        <- Anova(activity_G, type = "3")$Chisq[2]

# PERMUTED G*E, G, AND SEX EFFECTS ON NETWORK POSITIONS AND ACTIVITY # ----------------
indegree_GxE_null        = numeric(1000)
indegree_G_null          = numeric(1000)
indegree_Sex_null        = numeric(1000)
outdegree_GxE_null       = numeric(1000)
outdegree_G_null         = numeric(1000)
outdegree_Sex_null       = numeric(1000)
eigencentrality_GxE_null = numeric(1000)
eigencentrality_G_null   = numeric(1000)
eigencentrality_Sex_null = numeric(1000)
betweenness_GxE_null     = numeric(1000)
betweenness_G_null       = numeric(1000)
betweenness_Sex_null     = numeric(1000)
clusteringcoeff_GxE_null = numeric(1000)
clusteringcoeff_G_null   = numeric(1000)
clusteringcoeff_Sex_null = numeric(1000)
activity_GxE_null        = numeric(1000)
activity_G_null          = numeric(1000)
activity_Sex_null        = numeric(1000)
for (i in 1:1000) {
  network_rand = network %>% 
    group_by(video) %>%
    mutate(sex_rand = sample(sex, replace=F)) %>%
    group_by(video, sex) %>%
    mutate(genotype_rand = sample(genotype, replace=F))
  indegree_GxE_rand        <- glmer(indegree       ~ (food|genotype_rand) + food + sex      + (1|trial) + (1|observation), data = network_rand, family = poisson)
  indegree_G_rand          <- glmer(indegree       ~ (1|genotype_rand)    + food + sex      + (1|trial) + (1|observation), data = network_rand, family = poisson)
  indegree_NA_rand         <- glmer(indegree       ~                        food + sex      + (1|trial) + (1|observation), data = network_rand, family = poisson)
  indegree_Sex_rand        <- glmer(indegree       ~ (1|genotype)         + food + sex_rand + (1|trial) + (1|observation), data = network_rand, family = poisson)
  outdegree_GxE_rand       <- glmer(outdegree      ~ (food|genotype_rand) + sex*food        + (1|trial) + (1|observation), data = network_rand, family = poisson)
  outdegree_G_rand         <- glmer(outdegree      ~ (1|genotype_rand)    + sex*food        + (1|trial) + (1|observation), data = network_rand, family = poisson)
  outdegree_NA_rand        <- glmer(outdegree      ~                        sex*food        + (1|trial) + (1|observation), data = network_rand, family = poisson)
  outdegree_Sex_rand       <- glmer(outdegree      ~ (1|genotype)         + sex_rand*food   + (1|trial) + (1|observation), data = network_rand, family = poisson)
  eigencentrality_GxE_rand <- lmer(eigencentrality ~ (food|genotype_rand) + sex             + (1|individual),              data = network_rand)
  eigencentrality_G_rand   <- lmer(eigencentrality ~ (1|genotype_rand)    + sex             + (1|individual),              data = network_rand)
  eigencentrality_NA_rand  <- lmer(eigencentrality ~                        sex             + (1|individual),              data = network_rand)
  eigencentrality_Sex_rand <- lmer(eigencentrality ~ (1|genotype)         + sex_rand        + (1|individual),              data = network_rand)
  betweenness_GxE_rand     <- glmmTMB(betweenness  ~ (food|genotype_rand) + food + sex      + (1|trial),                   data = network_rand, family = "poisson", ziformula = ~1)
  betweenness_G_rand       <- glmmTMB(betweenness  ~ (1|genotype_rand)    + food + sex      + (1|trial),                   data = network_rand, family = "poisson", ziformula = ~1)
  betweenness_NA_rand      <- glmmTMB(betweenness  ~                        food + sex      + (1|trial),                   data = network_rand, family = "poisson", ziformula = ~1)
  betweenness_Sex_rand     <- glmmTMB(betweenness  ~ (1|genotype)         + food + sex_rand + (1|trial),                   data = network_rand, family = "poisson", ziformula = ~1)
  clusteringcoeff_GxE_rand <- lmer(clusteringcoeff ~ (food|genotype_rand) + sex*food        + (1|trial),                   data = network_rand)
  clusteringcoeff_G_rand   <- lmer(clusteringcoeff ~ (1|genotype_rand)    + sex*food        + (1|trial),                   data = network_rand)
  clusteringcoeff_NA_rand  <- lmer(clusteringcoeff ~                        sex*food        + (1|trial),                   data = network_rand)
  clusteringcoeff_Sex_rand <- lmer(clusteringcoeff ~ (1|genotype)         + sex_rand*food   + (1|trial),                   data = network_rand)
  activity_GxE_rand        <- lmer(activity        ~ (food|genotype_rand) + sex*food        + (1|trial),                   data = network_rand)
  activity_G_rand          <- lmer(activity        ~ (1|genotype_rand)    + sex*food        + (1|trial),                   data = network_rand)
  activity_NA_rand         <- lmer(activity        ~                        sex*food        + (1|trial),                   data = network_rand)
  activity_Sex_rand        <- lmer(activity        ~ (1|genotype)         + sex_rand*food   + (1|trial),                   data = network_rand)
  indegree_GxE_null[i]        <- anova(indegree_GxE_rand, indegree_G_rand)$Chisq[2]
  outdegree_GxE_null[i]       <- anova(outdegree_GxE_rand, outdegree_G_rand)$Chisq[2]
  eigencentrality_GxE_null[i] <- anova(eigencentrality_GxE_rand, eigencentrality_G_rand)$Chisq[2]
  betweenness_GxE_null[i]     <- anova(betweenness_GxE_rand, betweenness_G_rand)$Chisq[2]
  clusteringcoeff_GxE_null[i] <- anova(clusteringcoeff_GxE_rand, clusteringcoeff_G_rand)$Chisq[2]
  activity_GxE_null[i]        <- anova(activity_GxE_rand, activity_G_rand)$Chisq[2]
  indegree_G_null[i]          <- anova(indegree_G_rand, indegree_NA_rand)$Chisq[2]
  outdegree_G_null[i]         <- anova(outdegree_G_rand, outdegree_NA_rand)$Chisq[2]
  eigencentrality_G_null[i]   <- anova(eigencentrality_G_rand, eigencentrality_NA_rand)$Chisq[2]
  betweenness_G_null[i]       <- anova(betweenness_G_rand, betweenness_NA_rand)$Chisq[2]
  clusteringcoeff_G_null[i]   <- anova(clusteringcoeff_G_rand, clusteringcoeff_NA_rand)$Chisq[2]
  activity_G_null[i]          <- anova(activity_G_rand, activity_NA_rand)$Chisq[2]
  indegree_Sex_null[i]        <- Anova(indegree_Sex_rand, type = "3")$Chisq[3]
  outdegree_Sex_null[i]       <- Anova(outdegree_Sex_rand, type = "3")$Chisq[2]
  eigencentrality_Sex_null[i] <- Anova(eigencentrality_Sex_rand, type = "3")$Chisq[2]
  betweenness_Sex_null[i]     <- Anova(betweenness_Sex_rand, type = "3")$Chisq[3]
  clusteringcoeff_Sex_null[i] <- Anova(clusteringcoeff_Sex_rand, type = "3")$Chisq[2]
  activity_Sex_null[i]        <- Anova(activity_Sex_rand, type = "3")$Chisq[2]}
# Permuted significance of GxE
sum(indegree_GxE_obs<indegree_GxE_null)/1000                      # Permuted significance of GxE: p = 0.195
sum(outdegree_GxE_obs<outdegree_GxE_null)/1000                    # Permuted significance of GxE: p = 0.543
sum(eigencentrality_GxE_obs<eigencentrality_GxE_null)/1000        # Permuted significance of GxE: p = 0.997
sum(betweenness_GxE_obs<betweenness_GxE_null, na.rm = T)/1000     # Permuted significance of GxE: p = 0.853
sum(clusteringcoeff_GxE_obs<clusteringcoeff_GxE_null)/1000        # Permuted significance of GxE: p = 0.007
sum(activity_GxE_obs<activity_GxE_null)/1000                      # Permuted significance of GxE: p = 0.104
# Permuted significance of G
sum(indegree_G_obs<indegree_G_null)/1000                          # Permuted significance of G:   p < 0.001
sum(outdegree_G_obs<outdegree_G_null)/1000                        # Permuted significance of G:   p < 0.001
sum(eigencentrality_G_obs<eigencentrality_G_null)/1000            # Permuted significance of G:   p < 0.001 
sum(betweenness_G_obs<betweenness_G_null)/1000                    # Permuted significance of G:   p = 0.022
sum(clusteringcoeff_G_obs<clusteringcoeff_G_null)/1000            # Permuted significance of G:   p < 0.001
sum(activity_G_obs<activity_G_null)/1000                          # Permuted significance of G:   p < 0.001
# Permuted significance of Sex
sum(indegree_Sex_obs<indegree_Sex_null)/1000                      # Permuted significance of Sex: p = 0.002
sum(outdegree_Sex_obs<outdegree_Sex_null)/1000                    # Permuted significance of Sex: p = 0.001
sum(eigencentrality_Sex_obs<eigencentrality_Sex_null)/1000        # Permuted significance of Sex: p = 0.006
sum(betweenness_Sex_obs<betweenness_Sex_null)/1000                # Permuted significance of Sex: p = 0.683
sum(clusteringcoeff_Sex_obs<clusteringcoeff_Sex_null)/1000        # Permuted significance of Sex: p < 0.001
sum(activity_Sex_obs<activity_Sex_null)/1000                      # Permuted significance of Sex: p = 0.001

# OBSERVED SEX*E EFFECTS ON NETWORK POSITIONS (WITH ACTIVITY COVARIATE) # ----------------
indegree_G_act        <- glmer(indegree       ~  (1|genotype) + sex*food + (1|trial) + (1|observation) + activity.rescaled, data = network, family = poisson)
outdegree_G_act       <- glmer(outdegree      ~  (1|genotype) + sex*food + (1|trial) + (1|observation) + activity.rescaled, data = network, family = poisson)
eigencentrality_G_act <- lmer(eigencentrality ~  (1|genotype) + sex*food + (1|individual)              + activity.rescaled, data = network)
betweenness_G_act     <- glmmTMB(betweenness  ~  (1|genotype) + sex*food + (1|trial)                   + activity.rescaled, data = network, family = "poisson", ziformula = ~1)
clusteringcoeff_G_act <- lmer(clusteringcoeff ~  (1|genotype) + sex*food + (1|trial)                   + activity.rescaled, data = network)
# Observed significance of Sex*E
Anova(indegree_G_act, type = "3")         # Observed significance of Sex*E: Chisq = 3.5035,   Df = 4,  p = 0.47735
Anova(outdegree_G_act, type = "3")        # Observed significance of Sex*E: Chisq = 12.5983,  Df = 4,  p = 0.0134146
Anova(eigencentrality_G_act, type = "3")  # Observed significance of Sex*E: Chisq = 4.7520,   Df = 4,  p = 0.31370
Anova(betweenness_G_act, type = "3")      # Observed significance of Sex*E: Chisq = 56.0581,  Df = 4,  p = 1.95e-11
Anova(clusteringcoeff_G_act, type = "3")  # Observed significance of Sex*E: Chisq = 36.7106,  Df = 4,  p = 2.066e-07
# Save observed test statistics as values
indegree_SexxE_Act_obs <- Anova(indegree_G_act, type = "3")$Chisq[5]
outdegree_SexxE_Act_obs <- Anova(outdegree_G_act, type = "3")$Chisq[5]
eigencentrality_SexxE_Act_obs <- Anova(eigencentrality_G_act, type = "3")$Chisq[5]
betweenness_SexxE_Act_obs <- Anova(betweenness_G_act, type = "3")$Chisq[5]
clusteringcoeff_SexxE_Act_obs <- Anova(clusteringcoeff_G_act, type = "3")$Chisq[5]

# PERMUTED SEX*E EFFECTS ON NETWORK POSITIONS (WITH ACTIVITY COVARIATE) # ----------------
indegree_SexxE_Act_null        = numeric(1000)
outdegree_SexxE_Act_null       = numeric(1000)
eigencentrality_SexxE_Act_null = numeric(1000)
betweenness_SexxE_Act_null     = numeric(1000)
clusteringcoeff_SexxE_Act_null = numeric(1000)
for (i in 1:1000) {
  network_rand = network %>% 
    group_by(video) %>%
    mutate(sex_rand = sample(sex, replace=F))
  indegree_SexxE_Act_rand        <- glmer(indegree       ~ (1|genotype) + sex_rand*food + (1|trial) + (1|observation) + activity.rescaled, data = network_rand, family = poisson)
  outdegree_SexxE_Act_rand       <- glmer(outdegree      ~ (1|genotype) + sex_rand*food + (1|trial) + (1|observation) + activity.rescaled, data = network_rand, family = poisson)
  eigencentrality_SexxE_Act_rand <- lmer(eigencentrality ~ (1|genotype) + sex_rand*food + (1|individual)              + activity.rescaled, data = network_rand)
  betweenness_SexxE_Act_rand     <- glmmTMB(betweenness  ~ (1|genotype) + sex_rand*food + (1|trial)                   + activity.rescaled, data = network_rand, family = "poisson", ziformula = ~1)
  clusteringcoeff_SexxE_Act_rand <- lmer(clusteringcoeff ~ (1|genotype) + sex_rand*food + (1|trial)                   + activity.rescaled, data = network_rand)
  indegree_SexxE_Act_null[i]        <- Anova(indegree_SexxE_Act_rand, type = "3")$Chisq[5]
  outdegree_SexxE_Act_null[i]       <- Anova(outdegree_SexxE_Act_rand, type = "3")$Chisq[5]
  eigencentrality_SexxE_Act_null[i] <- Anova(eigencentrality_SexxE_Act_rand, type = "3")$Chisq[5]
  betweenness_SexxE_Act_null[i]     <- Anova(betweenness_SexxE_Act_rand, type = "3")$Chisq[5]
  clusteringcoeff_SexxE_Act_null[i] <- Anova(clusteringcoeff_SexxE_Act_rand, type = "3")$Chisq[5]}
sum(indegree_SexxE_Act_obs<indegree_SexxE_Act_null)/1000                 # Permuted significance of Sex*E: p = 0.380
sum(outdegree_SexxE_Act_obs<outdegree_SexxE_Act_null)/1000               # Permuted significance of Sex*E: p = 0.003
sum(eigencentrality_SexxE_Act_obs<eigencentrality_SexxE_Act_null)/1000   # Permuted significance of Sex*E: p = 0.291
sum(betweenness_SexxE_Act_obs<betweenness_SexxE_Act_null)/1000           # Permuted significance of Sex*E: p = 0.798
sum(clusteringcoeff_SexxE_Act_obs<clusteringcoeff_SexxE_Act_null)/1000   # Permuted significance of Sex*E: p < 0.001

# OBSERVED G*E, G, E, SEX, AND ACTIVITY EFFECTS ON NETWORK POSITIONS (WITH ACTIVITY COVARIATE) # ------------------
indegree_GxE_act        <- glmer(indegree       ~  (food|genotype) + food + sex + (1|trial) + (1|observation) + activity.rescaled, data = network, family = poisson)
indegree_G_act          <- glmer(indegree       ~  (1|genotype)    + food + sex + (1|trial) + (1|observation) + activity.rescaled, data = network, family = poisson)
indegree_NA_act         <- glmer(indegree       ~                    food + sex + (1|trial) + (1|observation) + activity.rescaled, data = network, family = poisson)
outdegree_GxE_act       <- glmer(outdegree      ~  (food|genotype) + sex*food   + (1|trial) + (1|observation) + activity.rescaled, data = network, family = poisson)
outdegree_G_act         <- glmer(outdegree      ~  (1|genotype)    + sex*food   + (1|trial) + (1|observation) + activity.rescaled, data = network, family = poisson)
outdegree_NA_act        <- glmer(outdegree      ~                    sex*food   + (1|trial) + (1|observation) + activity.rescaled, data = network, family = poisson)
eigencentrality_GxE_act <- lmer(eigencentrality ~  (food|genotype) + sex        + (1|individual)              + activity.rescaled, data = network)
eigencentrality_G_act   <- lmer(eigencentrality ~  (1|genotype)    + sex        + (1|individual)              + activity.rescaled, data = network)
eigencentrality_NA_act  <- lmer(eigencentrality ~                    sex        + (1|individual)              + activity.rescaled, data = network)
betweenness_GxE_act     <- glmmTMB(betweenness  ~  (food|genotype) + food + sex + (1|trial)                   + activity.rescaled, data = network, family = "poisson", ziformula = ~1)
betweenness_G_act       <- glmmTMB(betweenness  ~  (1|genotype)    + food + sex + (1|trial)                   + activity.rescaled, data = network, family = "poisson", ziformula = ~1)
betweenness_NA_act      <- glmmTMB(betweenness  ~                    food + sex + (1|trial)                   + activity.rescaled, data = network, family = "poisson", ziformula = ~1)
clusteringcoeff_GxE_act <- lmer(clusteringcoeff ~  (food|genotype) + sex*food   + (1|trial)                   + activity.rescaled, data = network)
clusteringcoeff_G_act   <- lmer(clusteringcoeff ~  (1|genotype)    + sex*food   + (1|trial)                   + activity.rescaled, data = network)
clusteringcoeff_NA_act  <- lmer(clusteringcoeff ~                    sex*food   + (1|trial)                   + activity.rescaled, data = network)
# Observed significance of GxE
anova(indegree_GxE_act, indegree_G_act)                         # Observed significance of GxE:      LR    = 4.2308,     Df = 14, p = 0.9939
anova(outdegree_GxE_act, outdegree_G_act)                       # Observed significance of GxE:      LR    = 1.0519,     Df = 14, p = 1
anova(eigencentrality_GxE_act, eigencentrality_G_act)           # Observed significance of GxE:      LR    = 0,          Df = 14, p = 1
anova(betweenness_GxE_act, betweenness_G_act)                   # Observed significance of GxE:      LR    = 1775.9,     Df = 14, p < 2.2e-16
anova(clusteringcoeff_GxE_act, clusteringcoeff_G_act)           # Observed significance of GxE:      LR    = 3.9583,     Df = 14, p = 0.9957
# Observed significance of G
anova(indegree_G_act, indegree_NA_act)                          # Observed significance of G:        LR    = 25.339,     Df = 1,  p = 4.808e-07
anova(outdegree_G_act, outdegree_NA_act)                        # Observed significance of G:        LR    = 26.814,     Df = 1,  p = 2.24e-07
anova(eigencentrality_G_act, eigencentrality_NA_act)            # Observed significance of G:        LR    = 24.236,     Df = 1,  p = 8.521e-07
anova(betweenness_G_act, betweenness_NA_act)                    # Observed significance of G:        LR    = 984.26,     Df = 1,  p < 2.2e-16 
anova(clusteringcoeff_G_act, clusteringcoeff_NA_act)            # Observed significance of G:        LR    = 57.175,     Df = 1,  p = 3.987e-14
# Observed significance of E
Anova(indegree_G_act, type = "3")                               # Observed significance of E:        Chisq = 6.1658,     Df = 4,  p = 0.18710
Anova(outdegree_G_act, type = "3")                              # Observed significance of E:        Chisq = 9.3142,     Df = 4,  p = 0.0537091
Anova(betweenness_G_act, type = "3")                            # Observed significance of E:        Chisq = 12.8286,    Df = 4,  p = 0.01214
Anova(clusteringcoeff_G_act, type = "3")                        # Observed significance of E:        Chisq = 9.9298,     Df = 4,  p = 0.0416270
# Observed significance of Sex
Anova(indegree_G_act, type = "3")                               # Observed significance of Sex:      Chisq = 4.8146,     Df = 1,  p = 0.02822
Anova(outdegree_G_act, type = "3")                              # Observed significance of Sex:      Chisq = 11.9477,    Df = 1,  p = 0.0005472
Anova(eigencentrality_G_act, type = "3")                        # Observed significance of Sex:      Chisq = 8.7982,     Df = 1,  p = 0.003015
Anova(betweenness_G_act, type = "3")                            # Observed significance of Sex:      Chisq = 0.6535,     Df = 1,  p = 0.41887
Anova(clusteringcoeff_G_act, type = "3")                        # Observed significance of Sex:      Chisq = 114.1249,   Df = 1,  p < 2.2e-16
# Observed significance of Activity
Anova(indegree_G_act, type = "3")                               # Observed significance of Activity: Chisq = 4.6241,     Df = 1,  p = 0.03153
Anova(outdegree_G_act, type = "3")                              # Observed significance of Activity: Chisq = 5.0013,     Df = 1,  p = 0.0253277
Anova(eigencentrality_G_act, type = "3")                        # Observed significance of Activity: Chisq = 0.0060,     Df = 1,  p = 0.938493
Anova(betweenness_G_act, type = "3")                            # Observed significance of Activity: Chisq = 2774.8466,  Df = 1,  p < 2e-16 
Anova(clusteringcoeff_G_act, type = "3")                        # Observed significance of Activity: Chisq = 13.6008,    Df = 1,  p = 0.0002261
# Save observed test statistics as values
indegree_GxE_Act_obs        <- anova(indegree_GxE_act, indegree_G_act)$Chisq[2]
outdegree_GxE_Act_obs       <- anova(outdegree_GxE_act, outdegree_G_act)$Chisq[2]
eigencentrality_GxE_Act_obs <- anova(eigencentrality_GxE_act, eigencentrality_G_act)$Chisq[2]
betweenness_GxE_Act_obs     <- anova(betweenness_GxE_act, betweenness_G_act)$Chisq[2]
clusteringcoeff_GxE_Act_obs <- anova(clusteringcoeff_GxE_act, clusteringcoeff_G_act)$Chisq[2]
indegree_G_Act_obs          <- anova(indegree_G_act, indegree_NA_act)$Chisq[2]
outdegree_G_Act_obs         <- anova(outdegree_G_act, outdegree_NA_act)$Chisq[2]
eigencentrality_G_Act_obs   <- anova(eigencentrality_G_act, eigencentrality_NA_act)$Chisq[2]
betweenness_G_Act_obs       <- anova(betweenness_G_act, betweenness_NA_act)$Chisq[2]
clusteringcoeff_G_Act_obs   <- anova(clusteringcoeff_G_act, clusteringcoeff_NA_act)$Chisq[2]
indegree_Sex_Act_obs        <- Anova(indegree_G_act, type = "3")$Chisq[3]
outdegree_Sex_Act_obs       <- Anova(outdegree_G_act, type = "3")$Chisq[2]
eigencentrality_Sex_Act_obs <- Anova(eigencentrality_G_act, type = "3")$Chisq[2]
betweenness_Sex_Act_obs     <- Anova(betweenness_G_act, type = "3")$Chisq[3]
clusteringcoeff_Sex_Act_obs <- Anova(clusteringcoeff_G_act, type = "3")$Chisq[2]
indegree_Act_obs            <- Anova(indegree_G_act, type = "3")$Chisq[4]
outdegree_Act_obs           <- Anova(outdegree_G_act, type = "3")$Chisq[4]
eigencentrality_Act_obs     <- Anova(eigencentrality_G_act, type = "3")$Chisq[3]
betweenness_Act_obs         <- Anova(betweenness_G_act, type = "3")$Chisq[4]
clusteringcoeff_Act_obs     <- Anova(clusteringcoeff_G_act, type = "3")$Chisq[4]

# PERMUTED G*E, G, SEX, AND ACTIVITY EFFECTS ON NETWORK POSITIONS (WITH ACTIVITY COVARIATE) # ------------------
indegree_GxE_Act_null        = numeric(1000)
indegree_G_Act_null          = numeric(1000)
indegree_Sex_Act_null        = numeric(1000)
indegree_Act_null            = numeric(1000)
outdegree_GxE_Act_null       = numeric(1000)
outdegree_G_Act_null         = numeric(1000)
outdegree_Sex_Act_null       = numeric(1000)
outdegree_Act_null           = numeric(1000)
eigencentrality_GxE_Act_null = numeric(1000)
eigencentrality_G_Act_null   = numeric(1000)
eigencentrality_Sex_Act_null = numeric(1000)
eigencentrality_Act_null     = numeric(1000)
betweenness_GxE_Act_null     = numeric(1000)
betweenness_G_Act_null       = numeric(1000)
betweenness_Sex_Act_null     = numeric(1000)
betweenness_Act_null         = numeric(1000)
clusteringcoeff_GxE_Act_null = numeric(1000)
clusteringcoeff_G_Act_null   = numeric(1000)
clusteringcoeff_Sex_Act_null = numeric(1000)
clusteringcoeff_Act_null     = numeric(1000)
for (i in 1:1000) {
  network_rand = network %>% 
    group_by(video) %>%
    mutate(sex_rand = sample(sex, replace=F)) %>%
    group_by(video, sex) %>%
    mutate(genotype_rand = sample(genotype, replace=F)) %>%
    mutate(activity.rescaled_rand = sample(activity.rescaled, replace=F))
  indegree_GxE_act_rand        <- glmer(indegree       ~  (food|genotype_rand) + food + sex      + (1|trial) + (1|observation) + activity.rescaled,      data = network_rand, family = poisson)
  indegree_G_act_rand          <- glmer(indegree       ~  (1|genotype_rand)    + food + sex      + (1|trial) + (1|observation) + activity.rescaled,      data = network_rand, family = poisson)
  indegree_NA_act_rand         <- glmer(indegree       ~                         food + sex      + (1|trial) + (1|observation) + activity.rescaled,      data = network_rand, family = poisson)
  indegree_Sex_Act_rand        <- glmer(indegree       ~  (1|genotype)         + food + sex_rand + (1|trial) + (1|observation) + activity.rescaled,      data = network_rand, family = poisson)
  indegree_Act_rand            <- glmer(indegree       ~  (1|genotype)         + food + sex      + (1|trial) + (1|observation) + activity.rescaled_rand, data = network_rand, family = poisson)
  outdegree_GxE_act_rand       <- glmer(outdegree      ~  (food|genotype_rand) + sex*food        + (1|trial) + (1|observation) + activity.rescaled,      data = network_rand, family = poisson)
  outdegree_G_act_rand         <- glmer(outdegree      ~  (1|genotype_rand)    + sex*food        + (1|trial) + (1|observation) + activity.rescaled,      data = network_rand, family = poisson)
  outdegree_NA_act_rand        <- glmer(outdegree      ~                         sex*food        + (1|trial) + (1|observation) + activity.rescaled,      data = network_rand, family = poisson)
  outdegree_Sex_Act_rand       <- glmer(outdegree      ~  (1|genotype)         + sex_rand*food   + (1|trial) + (1|observation) + activity.rescaled,      data = network_rand, family = poisson)
  outdegree_Act_rand           <- glmer(outdegree      ~  (1|genotype)         + sex*food        + (1|trial) + (1|observation) + activity.rescaled_rand, data = network_rand, family = poisson)
  eigencentrality_GxE_act_rand <- lmer(eigencentrality ~  (food|genotype_rand) + sex             + (1|individual)              + activity.rescaled,      data = network_rand)
  eigencentrality_G_act_rand   <- lmer(eigencentrality ~  (1|genotype_rand)    + sex             + (1|individual)              + activity.rescaled,      data = network_rand)
  eigencentrality_NA_act_rand  <- lmer(eigencentrality ~                         sex             + (1|individual)              + activity.rescaled,      data = network_rand)
  eigencentrality_Sex_Act_rand <- lmer(eigencentrality ~  (1|genotype)         + sex_rand        + (1|individual)              + activity.rescaled,      data = network_rand)
  eigencentrality_Act_rand     <- lmer(eigencentrality ~  (1|genotype)         + sex             + (1|individual)              + activity.rescaled_rand, data = network_rand)
  betweenness_GxE_act_rand     <- glmmTMB(betweenness  ~  (food|genotype_rand) + food + sex      + (1|trial)                   + activity.rescaled,      data = network_rand, family = "poisson", ziformula = ~1)
  betweenness_G_act_rand       <- glmmTMB(betweenness  ~  (1|genotype_rand)    + food + sex      + (1|trial)                   + activity.rescaled,      data = network_rand, family = "poisson", ziformula = ~1)
  betweenness_NA_act_rand      <- glmmTMB(betweenness  ~                         food + sex      + (1|trial)                   + activity.rescaled,      data = network_rand, family = "poisson", ziformula = ~1)
  betweenness_Sex_Act_rand     <- glmmTMB(betweenness  ~  (1|genotype)         + food + sex_rand + (1|trial)                   + activity.rescaled,      data = network_rand, family = "poisson", ziformula = ~1)
  betweenness_Act_rand         <- glmmTMB(betweenness  ~  (1|genotype)         + food + sex      + (1|trial)                   + activity.rescaled_rand, data = network_rand, family = "poisson", ziformula = ~1)
  clusteringcoeff_GxE_act_rand <- lmer(clusteringcoeff ~  (food|genotype_rand) + sex*food        + (1|trial)                   + activity.rescaled,      data = network_rand)
  clusteringcoeff_G_act_rand   <- lmer(clusteringcoeff ~  (1|genotype_rand)    + sex*food        + (1|trial)                   + activity.rescaled,      data = network_rand)
  clusteringcoeff_NA_act_rand  <- lmer(clusteringcoeff ~                         sex*food        + (1|trial)                   + activity.rescaled,      data = network_rand)
  clusteringcoeff_Sex_Act_rand <- lmer(clusteringcoeff ~  (1|genotype)         + sex_rand*food   + (1|trial)                   + activity.rescaled,      data = network_rand)
  clusteringcoeff_Act_rand     <- lmer(clusteringcoeff ~  (1|genotype)         + sex*food        + (1|trial)                   + activity.rescaled_rand, data = network_rand)
  indegree_GxE_Act_null[i]        <- anova(indegree_GxE_act_rand, indegree_G_act_rand)$Chisq[2]
  outdegree_GxE_Act_null[i]       <- anova(outdegree_GxE_act_rand, outdegree_G_act_rand)$Chisq[2]
  eigencentrality_GxE_Act_null[i] <- anova(eigencentrality_GxE_act_rand, eigencentrality_G_act_rand)$Chisq[2]
  betweenness_GxE_Act_null[i]     <- anova(betweenness_GxE_act_rand, betweenness_G_act_rand)$Chisq[2]
  clusteringcoeff_GxE_Act_null[i] <- anova(clusteringcoeff_GxE_act_rand, clusteringcoeff_G_act_rand)$Chisq[2]
  indegree_G_Act_null[i]          <- anova(indegree_G_act_rand, indegree_NA_act_rand)$Chisq[2]
  outdegree_G_Act_null[i]         <- anova(outdegree_G_act_rand, outdegree_NA_act_rand)$Chisq[2]
  eigencentrality_G_Act_null[i]   <- anova(eigencentrality_G_act_rand, eigencentrality_NA_act_rand)$Chisq[2]
  betweenness_G_Act_null[i]       <- anova(betweenness_G_act_rand, betweenness_NA_act_rand)$Chisq[2]
  clusteringcoeff_G_Act_null[i]   <- anova(clusteringcoeff_G_act_rand, clusteringcoeff_NA_act_rand)$Chisq[2]
  indegree_Sex_Act_null[i]        <- Anova(indegree_Sex_Act_rand, type = "3")$Chisq[3]
  outdegree_Sex_Act_null[i]       <- Anova(outdegree_Sex_Act_rand, type = "3")$Chisq[2]
  eigencentrality_Sex_Act_null[i] <- Anova(eigencentrality_Sex_Act_rand, type = "3")$Chisq[2]
  betweenness_Sex_Act_null[i]     <- Anova(betweenness_Sex_Act_rand, type = "3")$Chisq[3]
  clusteringcoeff_Sex_Act_null[i] <- Anova(clusteringcoeff_Sex_Act_rand, type = "3")$Chisq[2]
  indegree_Act_null[i]            <- Anova(indegree_Act_rand, type = "3")$Chisq[4]
  outdegree_Act_null[i]           <- Anova(outdegree_Act_rand, type = "3")$Chisq[4]
  eigencentrality_Act_null[i]     <- Anova(eigencentrality_Act_rand, type = "3")$Chisq[3]
  betweenness_Act_null[i]         <- Anova(betweenness_Act_rand, type = "3")$Chisq[4]
  clusteringcoeff_Act_null[i]     <- Anova(clusteringcoeff_Act_rand, type = "3")$Chisq[4]}
# Permuted significance of GxE
sum(indegree_GxE_Act_obs<indegree_GxE_Act_null)/1000                       # Permuted significance of GxE:      p = 0.206
sum(outdegree_GxE_Act_obs<outdegree_GxE_Act_null)/1000                     # Permuted significance of GxE:      p = 0.555
sum(eigencentrality_GxE_Act_obs<eigencentrality_GxE_Act_null)/1000         # Permuted significance of GxE:      p = 0.997
sum(betweenness_GxE_Act_obs<betweenness_GxE_Act_null, na.rm = T)/1000      # Permuted significance of GxE:      p = 0.819
sum(clusteringcoeff_GxE_Act_obs<clusteringcoeff_GxE_Act_null)/1000         # Permuted significance of GxE:      p = 0.016
# Permuted significance of G
sum(indegree_G_Act_obs<indegree_G_Act_null)/1000                           # Permuted significance of G:        p < 0.001 
sum(outdegree_G_Act_obs<outdegree_G_Act_null)/1000                         # Permuted significance of G:        p < 0.001
sum(eigencentrality_G_Act_obs<eigencentrality_G_Act_null)/1000             # Permuted significance of G:        p < 0.001
sum(betweenness_G_Act_obs<betweenness_G_Act_null)/1000                     # Permuted significance of G:        p = 0.012
sum(clusteringcoeff_G_Act_obs<clusteringcoeff_G_Act_null)/1000             # Permuted significance of G:        p < 0.001
# Permuted significance of Sex
sum(indegree_Sex_Act_obs<indegree_Sex_Act_null)/1000                       # Permuted significance of Sex:      p = 0.015
sum(outdegree_Sex_Act_obs<outdegree_Sex_Act_null)/1000                     # Permuted significance of Sex:      p = 0.001
sum(eigencentrality_Sex_Act_obs<eigencentrality_Sex_Act_null)/1000         # Permuted significance of Sex:      p = 0.002
sum(betweenness_Sex_Act_obs<betweenness_Sex_Act_null)/1000                 # Permuted significance of Sex:      p = 0.888
sum(clusteringcoeff_Sex_Act_obs<clusteringcoeff_Sex_Act_null)/1000         # Permuted significance of Sex:      p < 0.001
# Permuted significance of Activity
sum(indegree_Act_obs<indegree_Act_null)/1000                               # Permuted significance of Activity: p = 0.570
sum(outdegree_Act_obs<outdegree_Act_null)/1000                             # Permuted significance of Activity: p = 0.440
sum(eigencentrality_Act_obs<eigencentrality_Act_null)/1000                 # Permuted significance of Activity: p = 0.961
sum(betweenness_Act_obs<betweenness_Act_null)/1000                         # Permuted significance of Activity: p < 0.001
sum(clusteringcoeff_Act_obs<clusteringcoeff_Act_null)/1000                 # Permuted significance of Activity: p = 0.348

# OBSERVED SEX*E EFFECTS ON NETWORK POSITIONS (FOR P:C ONLY) # ----------------
network_pc <- network[which(network$food=='1' | network$food=='2' | network$food=='3'),]
indegree_G_pc        <- glmer(indegree       ~  (1|genotype) + sex*pc.rescaled + (1|trial) + (1|observation), data = network_pc, family = poisson)
outdegree_G_pc       <- glmer(outdegree      ~  (1|genotype) + sex*pc.rescaled + (1|trial) + (1|observation), data = network_pc, family = poisson)
eigencentrality_G_pc <- lmer(eigencentrality ~  (1|genotype) + sex*pc.rescaled + (1|individual),              data = network_pc)
betweenness_G_pc     <- glmmTMB(betweenness  ~  (1|genotype) + sex*pc.rescaled + (1|trial),                   data = network_pc, family = "poisson", ziformula = ~1)
clusteringcoeff_G_pc <- lmer(clusteringcoeff ~  (1|genotype) + sex*pc.rescaled + (1|trial),                   data = network_pc)
# Observed significance of Sex*PC
Anova(indegree_G_pc, type = "3")             # Observed significance of Sex*PC: Chisq = 0.9037,  Df = 1,  p = 0.3418021
Anova(outdegree_G_pc, type = "3")            # Observed significance of Sex*PC: Chisq = 9.5991,  Df = 1,  p = 0.001947
Anova(eigencentrality_G_pc, type = "3")      # Observed significance of Sex*PC: Chisq = 0.7726,  Df = 1,  p = 0.3794074
Anova(betweenness_G_pc, type = "3")          # Observed significance of Sex*PC: Chisq = 18.7612, Df = 1,  p = 1.481e-05
Anova(clusteringcoeff_G_pc, type = "3")      # Observed significance of Sex*PC: Chisq = 25.8897, Df = 1,  p = 3.615e-07
# Save observed test statistics as values
indegree_SexxE_pc_obs        <- Anova(indegree_G_pc, type = "3")$Chisq[4]
outdegree_SexxE_pc_obs       <- Anova(outdegree_G_pc, type = "3")$Chisq[4]
eigencentrality_SexxE_pc_obs <- Anova(eigencentrality_G_pc, type = "3")$Chisq[4]
betweenness_SexxE_pc_obs     <- Anova(betweenness_G_pc, type = "3")$Chisq[4]
clusteringcoeff_SexxE_pc_obs <- Anova(clusteringcoeff_G_pc, type = "3")$Chisq[4]

# PERMUTED SEX*E EFFECTS ON NETWORK POSITIONS (FOR P:C ONLY) # ----------------
indegree_SexxE_pc_null        <- numeric(1000)
outdegree_SexxE_pc_null       <- numeric(1000)
eigencentrality_SexxE_pc_null <- numeric(1000)
betweenness_SexxE_pc_null     <- numeric(1000)
clusteringcoeff_SexxE_pc_null <- numeric(1000)
for (i in 1:1000) {
  network_pc_rand = network_pc %>% 
    group_by(video) %>%
    mutate(sex_rand = sample(sex, replace=F))
  indegree_SexxE_pc_rand        <- glmer(indegree       ~ (1|genotype) + sex_rand*pc.rescaled + (1|trial) + (1|observation), data = network_pc_rand, family = poisson)
  outdegree_SexxE_pc_rand       <- glmer(outdegree      ~ (1|genotype) + sex_rand*pc.rescaled + (1|trial) + (1|observation), data = network_pc_rand, family = poisson)
  eigencentrality_SexxE_pc_rand <- lmer(eigencentrality ~ (1|genotype) + sex_rand*pc.rescaled + (1|individual),              data = network_pc_rand)
  betweenness_SexxE_pc_rand     <- glmmTMB(betweenness  ~ (1|genotype) + sex_rand*pc.rescaled + (1|trial),                   data = network_pc_rand, family = "poisson", ziformula = ~1)
  clusteringcoeff_SexxE_pc_rand <- lmer(clusteringcoeff ~ (1|genotype) + sex_rand*pc.rescaled + (1|trial),                   data = network_pc_rand)
  indegree_SexxE_pc_null[i]        <- Anova(indegree_SexxE_pc_rand, type = "3")$Chisq[4]
  outdegree_SexxE_pc_null[i]       <- Anova(outdegree_SexxE_pc_rand, type = "3")$Chisq[4]
  eigencentrality_SexxE_pc_null[i] <- Anova(eigencentrality_SexxE_pc_rand, type = "3")$Chisq[4]
  betweenness_SexxE_pc_null[i]     <- Anova(betweenness_SexxE_pc_rand, type = "3")$Chisq[4]
  clusteringcoeff_SexxE_pc_null[i] <- Anova(clusteringcoeff_SexxE_pc_rand, type = "3")$Chisq[4]}
sum(indegree_SexxE_pc_obs<indegree_SexxE_pc_null)/1000                 # Permuted significance of Sex*PC: p = 0.301
sum(outdegree_SexxE_pc_obs<outdegree_SexxE_pc_null)/1000               # Permuted significance of Sex*PC: p < 0.001
sum(eigencentrality_SexxE_pc_obs<eigencentrality_SexxE_pc_null)/1000   # Permuted significance of Sex*PC: p = 0.390
sum(betweenness_SexxE_pc_obs<betweenness_SexxE_pc_null)/1000           # Permuted significance of Sex*PC: p = 0.451
sum(clusteringcoeff_SexxE_pc_obs<clusteringcoeff_SexxE_pc_null)/1000   # Permuted significance of Sex*PC: p < 0.001

# OBSERVED G*E, G, E, AND SEX EFFECTS ON NETWORK POSITIONS (FOR P:C ONLY) # ------------------
indegree_GxE_pc        <- glmer(indegree       ~ (pc.rescaled|genotype) + pc.rescaled + sex + (1|trial) + (1|observation), data = network_pc, family = poisson)
indegree_G_pc          <- glmer(indegree       ~ (1|genotype)           + pc.rescaled + sex + (1|trial) + (1|observation), data = network_pc, family = poisson)
indegree_NA_pc         <- glmer(indegree       ~                          pc.rescaled + sex + (1|trial) + (1|observation), data = network_pc, family = poisson)
outdegree_GxE_pc       <- glmer(outdegree      ~ (pc.rescaled|genotype) + sex*pc.rescaled   + (1|trial) + (1|observation), data = network_pc, family = poisson)
outdegree_G_pc         <- glmer(outdegree      ~ (1|genotype)           + sex*pc.rescaled   + (1|trial) + (1|observation), data = network_pc, family = poisson)
outdegree_NA_pc        <- glmer(outdegree      ~                          sex*pc.rescaled   + (1|trial) + (1|observation), data = network_pc, family = poisson)
eigencentrality_GxE_pc <- lmer(eigencentrality ~ (pc.rescaled|genotype) + sex               + (1|individual),              data = network_pc)
eigencentrality_G_pc   <- lmer(eigencentrality ~ (1|genotype)           + sex               + (1|individual),              data = network_pc)
eigencentrality_NA_pc  <- lmer(eigencentrality ~                          sex               + (1|individual),              data = network_pc)
betweenness_GxE_pc     <- glmmTMB(betweenness  ~ (pc.rescaled|genotype) + pc.rescaled + sex + (1|trial),                   data = network_pc, family = "poisson", ziformula = ~1)
betweenness_G_pc       <- glmmTMB(betweenness  ~ (1|genotype)           + pc.rescaled + sex + (1|trial),                   data = network_pc, family = "poisson", ziformula = ~1)
betweenness_NA_pc      <- glmmTMB(betweenness  ~                          pc.rescaled + sex + (1|trial),                   data = network_pc, family = "poisson", ziformula = ~1)
clusteringcoeff_GxE_pc <- lmer(clusteringcoeff ~ (pc.rescaled|genotype) + sex*pc.rescaled   + (1|trial),                   data = network_pc)
clusteringcoeff_G_pc   <- lmer(clusteringcoeff ~ (1|genotype)           + sex*pc.rescaled   + (1|trial),                   data = network_pc)
clusteringcoeff_NA_pc  <- lmer(clusteringcoeff ~                          sex*pc.rescaled   + (1|trial),                   data = network_pc)
# Observed significance of GxE
anova(indegree_GxE_pc, indegree_G_pc)                             # Observed significance of GxE:             LR = 0,           Df = 2,  p = 1
anova(outdegree_GxE_pc, outdegree_G_pc)                           # Observed significance of GxE:             LR = 0,           Df = 2,  p = 1
anova(eigencentrality_GxE_pc, eigencentrality_G_pc)               # Observed significance of GxE:             LR = 0.2362,      Df = 2,  p = 0.8886
anova(betweenness_GxE_pc, betweenness_G_pc)                       # Observed significance of GxE:             LR = 625.07,      Df = 2,  p < 2.2e-16
anova(clusteringcoeff_GxE_pc, clusteringcoeff_G_pc)               # Observed significance of GxE:             LR = 1.8429,      Df = 2,  p = 0.3979
# Observed significance of G
anova(indegree_G_pc, indegree_NA_pc)                              # Observed significance of G:               LR = 7.1597,      Df = 1,  p = 0.007456
anova(outdegree_G_pc, outdegree_NA_pc)                            # Observed significance of G:               LR = 11.131,      Df = 1,  p = 0.000849
anova(eigencentrality_G_pc, eigencentrality_NA_pc)                # Observed significance of G:               LR = 8.6954,      Df = 1,  p = 0.00319
anova(betweenness_G_pc, betweenness_NA_pc)                        # Observed significance of G:               LR = 804.04,      Df = 1,  p < 2.2e-16
anova(clusteringcoeff_G_pc, clusteringcoeff_NA_pc)                # Observed significance of G:               LR = 28.988,      Df = 1,  p = 7.283e-08
# Observed significance of PC
Anova(indegree_G_pc, type = "3")                                  # Observed significance of PC:              Chisq = 5.1746,   Df = 1,  p = 0.0229188
Anova(outdegree_G_pc, type = "3")                                 # Observed significance of PC:              Chisq = 8.1219,   Df = 1,  p = 0.004373
fligner.test(eigencentrality ~ pc.rescaled, data = network_pc)    # Observed significance of PC on variation: Chisq = 0.14684,  Df = 2,  p = 0.9292
Anova(betweenness_G_pc, type = "3")                               # Observed significance of PC:              Chisq = 0.5818,   Df = 1,  p = 0.44562
Anova(clusteringcoeff_G_pc, type = "3")                           # Observed significance of PC:              Chisq = 7.0574,   Df = 1,  p = 0.007894
# Observed significance of Sex
Anova(indegree_G_pc, type = "3")                                  # Observed significance of Sex:             Chisq = 11.5939,  Df = 1,  p = 0.0006617
Anova(outdegree_G_pc, type = "3")                                 # Observed significance of Sex:             Chisq = 0.0496,   Df = 1,  p = 0.823723
Anova(eigencentrality_G_pc, type = "3")                           # Observed significance of Sex:             Chisq = 13.48,    Df = 1,  p = 0.0002411
Anova(betweenness_G_pc, type = "3")                               # Observed significance of Sex:             Chisq = 7.6270,   Df = 1,  p = 0.00575
Anova(clusteringcoeff_G_pc, type = "3")                           # Observed significance of Sex:             Chisq = 36.1203,  Df = 1,  p = 1.855e-09
# Save observed test statistics as values
indegree_GxE_pc_obs        <- anova(indegree_GxE_pc, indegree_G_pc)$Chisq[2]
outdegree_GxE_pc_obs       <- anova(outdegree_GxE_pc, outdegree_G_pc)$Chisq[2]
eigencentrality_GxE_pc_obs <- anova(eigencentrality_GxE_pc, eigencentrality_G_pc)$Chisq[2]
betweenness_GxE_pc_obs     <- anova(betweenness_GxE_pc, betweenness_G_pc)$Chisq[2]
clusteringcoeff_GxE_pc_obs <- anova(clusteringcoeff_GxE_pc, clusteringcoeff_G_pc)$Chisq[2]
indegree_G_pc_obs          <- anova(indegree_G_pc, indegree_NA_pc)$Chisq[2]
outdegree_G_pc_obs         <- anova(outdegree_G_pc, outdegree_NA_pc)$Chisq[2]
eigencentrality_G_pc_obs   <- anova(eigencentrality_G_pc, eigencentrality_NA_pc)$Chisq[2]
betweenness_G_pc_obs       <- anova(betweenness_G_pc, betweenness_NA_pc)$Chisq[2]
clusteringcoeff_G_pc_obs   <- anova(clusteringcoeff_G_pc, clusteringcoeff_NA_pc)$Chisq[2]
indegree_Sex_pc_obs        <- Anova(indegree_G_pc, type = "3")$Chisq[3]
outdegree_Sex_pc_obs       <- Anova(outdegree_G_pc, type = "3")$Chisq[2]
eigencentrality_Sex_pc_obs <- Anova(eigencentrality_G_pc, type = "3")$Chisq[2]
betweenness_Sex_pc_obs     <- Anova(betweenness_G_pc, type = "3")$Chisq[3]
clusteringcoeff_Sex_pc_obs <- Anova(clusteringcoeff_G_pc, type = "3")$Chisq[2]

# PERMUTED G*E, G, AND SEX EFFECTS ON NETWORK POSITIONS (FOR P:C ONLY) # ------------------
indegree_GxE_pc_null        = numeric(1000)
indegree_G_pc_null          = numeric(1000)
indegree_Sex_pc_null        = numeric(1000)
outdegree_GxE_pc_null       = numeric(1000)
outdegree_G_pc_null         = numeric(1000)
outdegree_Sex_pc_null       = numeric(1000)
eigencentrality_GxE_pc_null = numeric(1000)
eigencentrality_G_pc_null   = numeric(1000)
eigencentrality_Sex_pc_null = numeric(1000)
betweenness_GxE_pc_null     = numeric(1000)
betweenness_G_pc_null       = numeric(1000)
betweenness_Sex_pc_null     = numeric(1000)
clusteringcoeff_GxE_pc_null = numeric(1000)
clusteringcoeff_G_pc_null   = numeric(1000)
clusteringcoeff_Sex_pc_null = numeric(1000)
for (i in 1:1000) {
  network_pc_rand = network_pc %>% 
    group_by(video) %>%
    mutate(sex_rand = sample(sex, replace=F)) %>%
    group_by(video, sex) %>%
    mutate(genotype_rand = sample(genotype, replace=F))
  indegree_GxE_pc_rand        <- glmer(indegree       ~ (pc.rescaled|genotype_rand) + pc.rescaled + sex      + (1|trial) + (1|observation), data = network_pc_rand, family = poisson)
  indegree_G_pc_rand          <- glmer(indegree       ~ (1|genotype_rand)           + pc.rescaled + sex      + (1|trial) + (1|observation), data = network_pc_rand, family = poisson)
  indegree_NA_pc_rand         <- glmer(indegree       ~                               pc.rescaled + sex      + (1|trial) + (1|observation), data = network_pc_rand, family = poisson)
  indegree_Sex_pc_rand        <- glmer(indegree       ~ (1|genotype)                + pc.rescaled + sex_rand + (1|trial) + (1|observation), data = network_pc_rand, family = poisson)
  outdegree_GxE_pc_rand       <- glmer(outdegree      ~ (pc.rescaled|genotype_rand) + sex*pc.rescaled        + (1|trial) + (1|observation), data = network_pc_rand, family = poisson)
  outdegree_G_pc_rand         <- glmer(outdegree      ~ (1|genotype_rand)           + sex*pc.rescaled        + (1|trial) + (1|observation), data = network_pc_rand, family = poisson)
  outdegree_NA_pc_rand        <- glmer(outdegree      ~                               sex*pc.rescaled        + (1|trial) + (1|observation), data = network_pc_rand, family = poisson)
  outdegree_Sex_pc_rand       <- glmer(outdegree      ~ (1|genotype)                + sex_rand*pc.rescaled   + (1|trial) + (1|observation), data = network_pc_rand, family = poisson)
  eigencentrality_GxE_pc_rand <- lmer(eigencentrality ~ (pc.rescaled|genotype_rand) + sex                    + (1|individual),              data = network_pc_rand)
  eigencentrality_G_pc_rand   <- lmer(eigencentrality ~ (1|genotype_rand)           + sex                    + (1|individual),              data = network_pc_rand)
  eigencentrality_NA_pc_rand  <- lmer(eigencentrality ~                               sex                    + (1|individual),              data = network_pc_rand)
  eigencentrality_Sex_pc_rand <- lmer(eigencentrality ~ (1|genotype)                + sex_rand               + (1|individual),              data = network_pc_rand)
  betweenness_GxE_pc_rand     <- glmmTMB(betweenness  ~ (pc.rescaled|genotype_rand) + pc.rescaled + sex      + (1|trial),                   data = network_pc_rand, family = "poisson", ziformula = ~1)
  betweenness_G_pc_rand       <- glmmTMB(betweenness  ~ (1|genotype_rand)           + pc.rescaled + sex      + (1|trial),                   data = network_pc_rand, family = "poisson", ziformula = ~1)
  betweenness_NA_pc_rand      <- glmmTMB(betweenness  ~                               pc.rescaled + sex      + (1|trial),                   data = network_pc_rand, family = "poisson", ziformula = ~1)
  betweenness_Sex_pc_rand     <- glmmTMB(betweenness  ~ (1|genotype)                + pc.rescaled + sex_rand + (1|trial),                   data = network_pc_rand, family = "poisson", ziformula = ~1)
  clusteringcoeff_GxE_pc_rand <- lmer(clusteringcoeff ~ (pc.rescaled|genotype_rand) + sex*pc.rescaled        + (1|trial),                   data = network_pc_rand)
  clusteringcoeff_G_pc_rand   <- lmer(clusteringcoeff ~ (1|genotype_rand)           + sex*pc.rescaled        + (1|trial),                   data = network_pc_rand)
  clusteringcoeff_NA_pc_rand  <- lmer(clusteringcoeff ~                               sex*pc.rescaled        + (1|trial),                   data = network_pc_rand)
  clusteringcoeff_Sex_pc_rand <- lmer(clusteringcoeff ~ (1|genotype)                + sex_rand*pc.rescaled   + (1|trial),                   data = network_pc_rand)
  indegree_GxE_pc_null[i]        <- anova(indegree_GxE_pc_rand, indegree_G_pc_rand)$Chisq[2]
  outdegree_GxE_pc_null[i]       <- anova(outdegree_GxE_pc_rand, outdegree_G_pc_rand)$Chisq[2]
  eigencentrality_GxE_pc_null[i] <- anova(eigencentrality_GxE_pc_rand, eigencentrality_G_pc_rand)$Chisq[2]
  betweenness_GxE_pc_null[i]     <- anova(betweenness_GxE_pc_rand, betweenness_G_pc_rand)$Chisq[2]
  clusteringcoeff_GxE_pc_null[i] <- anova(clusteringcoeff_GxE_pc_rand, clusteringcoeff_G_pc_rand)$Chisq[2]
  indegree_G_pc_null[i]          <- anova(indegree_G_pc_rand, indegree_NA_pc_rand)$Chisq[2]
  outdegree_G_pc_null[i]         <- anova(outdegree_G_pc_rand, outdegree_NA_pc_rand)$Chisq[2]
  eigencentrality_G_pc_null[i]   <- anova(eigencentrality_G_pc_rand, eigencentrality_NA_pc_rand)$Chisq[2]
  betweenness_G_pc_null[i]       <- anova(betweenness_G_pc_rand, betweenness_NA_pc_rand)$Chisq[2]
  clusteringcoeff_G_pc_null[i]   <- anova(clusteringcoeff_G_pc_rand, clusteringcoeff_NA_pc_rand)$Chisq[2]
  indegree_Sex_pc_null[i]        <- Anova(indegree_Sex_pc_rand, type = "3")$Chisq[3]
  outdegree_Sex_pc_null[i]       <- Anova(outdegree_Sex_pc_rand, type = "3")$Chisq[2]
  eigencentrality_Sex_pc_null[i] <- Anova(eigencentrality_Sex_pc_rand, type = "3")$Chisq[2]
  betweenness_Sex_pc_null[i]     <- Anova(betweenness_Sex_pc_rand, type = "3")$Chisq[3]
  clusteringcoeff_Sex_pc_null[i] <- Anova(clusteringcoeff_Sex_pc_rand, type = "3")$Chisq[2]}
# Permuted significance of GxPC
sum(indegree_GxE_pc_obs<indegree_GxE_pc_null)/1000                      # Permuted significance of GxPC: p = 0.591
sum(outdegree_GxE_pc_obs<outdegree_GxE_pc_null)/1000                    # Permuted significance of GxPC: p = 0.512
sum(eigencentrality_GxE_pc_obs<eigencentrality_GxE_pc_null)/1000        # Permuted significance of GxPC: p = 0.447
sum(betweenness_GxE_pc_obs<betweenness_GxE_pc_null)/1000                # Permuted significance of GxPC: p = 0.368
sum(clusteringcoeff_GxE_pc_obs<clusteringcoeff_GxE_pc_null)/1000        # Permuted significance of GxPC: p = 0.001
# Permuted significance of G
sum(indegree_G_pc_obs<indegree_G_pc_null)/1000                          # Permuted significance of G:    p < 0.001
sum(outdegree_G_pc_obs<outdegree_G_pc_null)/1000                        # Permuted significance of G:    p < 0.001
sum(eigencentrality_G_pc_obs<eigencentrality_G_pc_null)/1000            # Permuted significance of G:    p = 0.002
sum(betweenness_G_pc_obs<betweenness_G_pc_null)/1000                    # Permuted significance of G:    p = 0.100
sum(clusteringcoeff_G_pc_obs<clusteringcoeff_G_pc_null)/1000            # Permuted significance of G:    p < 0.001
# Permuted significance of Sex
sum(indegree_Sex_pc_obs<indegree_Sex_pc_null)/1000                      # Permuted significance of Sex:  p < 0.001
sum(outdegree_Sex_pc_obs<outdegree_Sex_pc_null)/1000                    # Permuted significance of Sex:  p = 0.799
sum(eigencentrality_Sex_pc_obs<eigencentrality_Sex_pc_null)/1000        # Permuted significance of Sex:  p < 0.001
sum(betweenness_Sex_pc_obs<betweenness_Sex_pc_null)/1000                # Permuted significance of Sex:  p = 0.632
sum(clusteringcoeff_Sex_pc_obs<clusteringcoeff_Sex_pc_null)/1000        # Permuted significance of Sex:  p < 0.001

# OBSERVED SEX*E EFFECTS ON NETWORK POSITIONS (FOR DILUTION ONLY) # ----------------
network_dilution <- network[which(network$food=='1' | network$food=='4' | network$food=='5'),]
indegree_G_dilution        <- glmer(indegree       ~  (1|genotype) + sex*dilution.rescaled + (1|trial) + (1|observation), data = network_dilution, family = poisson)
outdegree_G_dilution       <- glmer(outdegree      ~  (1|genotype) + sex*dilution.rescaled + (1|trial) + (1|observation), data = network_dilution, family = poisson)
eigencentrality_G_dilution <- lmer(eigencentrality ~  (1|genotype) + sex*dilution.rescaled + (1|individual),              data = network_dilution)
betweenness_G_dilution     <- glmmTMB(betweenness  ~  (1|genotype) + sex*dilution.rescaled + (1|trial),                   data = network_dilution, family = "poisson", ziformula = ~1)
clusteringcoeff_G_dilution <- lmer(clusteringcoeff ~  (1|genotype) + sex*dilution.rescaled + (1|trial),                   data = network_dilution)
# Observed significance of Sex*PC
Anova(indegree_G_dilution, type = "3")             # Observed significance of Sex*Dilution: Chisq = 0.0134,  Df = 1,  p = 0.9080
Anova(outdegree_G_dilution, type = "3")            # Observed significance of Sex*Dilution: Chisq = 4.5914,  Df = 1,  p = 0.03213
Anova(eigencentrality_G_dilution, type = "3")      # Observed significance of Sex*Dilution: Chisq = 0.2079,  Df = 1,  p = 0.6484
Anova(betweenness_G_dilution, type = "3")          # Observed significance of Sex*Dilution: Chisq = 11.4318, Df = 1,  p = 0.000722
Anova(clusteringcoeff_G_dilution, type = "3")      # Observed significance of Sex*Dilution: Chisq = 33.250,  Df = 1,  p = 8.104e-09
# Save observed test statistics as values
indegree_SexxE_dilution_obs <- Anova(indegree_G_dilution, type = "3")$Chisq[4]
outdegree_SexxE_dilution_obs <- Anova(outdegree_G_dilution, type = "3")$Chisq[4]
eigencentrality_SexxE_dilution_obs <- Anova(eigencentrality_G_dilution, type = "3")$Chisq[4]
betweenness_SexxE_dilution_obs <- Anova(betweenness_G_dilution, type = "3")$Chisq[4]
clusteringcoeff_SexxE_dilution_obs <- Anova(clusteringcoeff_G_dilution, type = "3")$Chisq[4]

# PERMUTED SEX*E EFFECTS ON NETWORK POSITIONS (FOR DILUTION ONLY) # ----------------
indegree_SexxE_dilution_null        <- numeric(1000)
outdegree_SexxE_dilution_null       <- numeric(1000)
eigencentrality_SexxE_dilution_null <- numeric(1000)
betweenness_SexxE_dilution_null     <- numeric(1000)
clusteringcoeff_SexxE_dilution_null <- numeric(1000)
for (i in 1:1000) {
  network_dilution_rand = network_dilution %>% 
    group_by(video) %>%
    mutate(sex_rand = sample(sex, replace=F))
  indegree_SexxE_dilution_rand        <- glmer(indegree       ~ (1|genotype) + sex_rand*dilution.rescaled + (1|trial) + (1|observation), data = network_dilution_rand, family = poisson)
  outdegree_SexxE_dilution_rand       <- glmer(outdegree      ~ (1|genotype) + sex_rand*dilution.rescaled + (1|trial) + (1|observation), data = network_dilution_rand, family = poisson)
  eigencentrality_SexxE_dilution_rand <- lmer(eigencentrality ~ (1|genotype) + sex_rand*dilution.rescaled + (1|individual),              data = network_dilution_rand)
  betweenness_SexxE_dilution_rand     <- glmmTMB(betweenness  ~ (1|genotype) + sex_rand*dilution.rescaled + (1|trial),                   data = network_dilution_rand, family = "poisson", ziformula = ~1)
  clusteringcoeff_SexxE_dilution_rand <- lmer(clusteringcoeff ~ (1|genotype) + sex_rand*dilution.rescaled + (1|trial),                   data = network_dilution_rand)
  indegree_SexxE_dilution_null[i]        <- Anova(indegree_SexxE_dilution_rand, type = "3")$Chisq[4]
  outdegree_SexxE_dilution_null[i]       <- Anova(outdegree_SexxE_dilution_rand, type = "3")$Chisq[4]
  eigencentrality_SexxE_dilution_null[i] <- Anova(eigencentrality_SexxE_dilution_rand, type = "3")$Chisq[4]
  betweenness_SexxE_dilution_null[i]     <- Anova(betweenness_SexxE_dilution_rand, type = "3")$Chisq[4]
  clusteringcoeff_SexxE_dilution_null[i] <- Anova(clusteringcoeff_SexxE_dilution_rand, type = "3")$Chisq[4]}
sum(indegree_SexxE_dilution_obs<indegree_SexxE_dilution_null)/1000                 # Permuted significance of Sex*Dilution: p = 0.907
sum(outdegree_SexxE_dilution_obs<outdegree_SexxE_dilution_null)/1000               # Permuted significance of Sex*Dilution: p = 0.033
sum(eigencentrality_SexxE_dilution_obs<eigencentrality_SexxE_dilution_null)/1000   # Permuted significance of Sex*Dilution: p = 0.621
sum(betweenness_SexxE_dilution_obs<betweenness_SexxE_dilution_null)/1000           # Permuted significance of Sex*Dilution: p = 0.626
sum(clusteringcoeff_SexxE_dilution_obs<clusteringcoeff_SexxE_dilution_null)/1000   # Permuted significance of Sex*Dilution: p < 0.001

# OBSERVED G*E, G, E, AND SEX EFFECTS ON NETWORK POSITIONS (FOR DILUTION ONLY) # ------------------
indegree_GxE_dilution        <- glmer(indegree       ~ (dilution.rescaled|genotype) + dilution.rescaled + sex + (1|trial) + (1|observation), data = network_dilution, family = poisson)
indegree_G_dilution          <- glmer(indegree       ~ (1|genotype)                 + dilution.rescaled + sex + (1|trial) + (1|observation), data = network_dilution, family = poisson)
indegree_NA_dilution         <- glmer(indegree       ~                                dilution.rescaled + sex + (1|trial) + (1|observation), data = network_dilution, family = poisson)
outdegree_GxE_dilution       <- glmer(outdegree      ~ (dilution.rescaled|genotype) + sex*dilution.rescaled   + (1|trial) + (1|observation), data = network_dilution, family = poisson)
outdegree_G_dilution         <- glmer(outdegree      ~ (1|genotype)                 + sex*dilution.rescaled   + (1|trial) + (1|observation), data = network_dilution, family = poisson)
outdegree_NA_dilution        <- glmer(outdegree      ~                                sex*dilution.rescaled   + (1|trial) + (1|observation), data = network_dilution, family = poisson)
eigencentrality_GxE_dilution <- lmer(eigencentrality ~ (dilution.rescaled|genotype) + sex                     + (1|individual),              data = network_dilution)
eigencentrality_G_dilution   <- lmer(eigencentrality ~ (1|genotype)                 + sex                     + (1|individual),              data = network_dilution)
eigencentrality_NA_dilution  <- lmer(eigencentrality ~                                sex                     + (1|individual),              data = network_dilution)
betweenness_GxE_dilution     <- glmmTMB(betweenness  ~ (dilution.rescaled|genotype) + dilution.rescaled + sex + (1|trial),                   data = network_dilution, family = "poisson", ziformula = ~1)
betweenness_G_dilution       <- glmmTMB(betweenness  ~ (1|genotype)                 + dilution.rescaled + sex + (1|trial),                   data = network_dilution, family = "poisson", ziformula = ~1)
betweenness_NA_dilution      <- glmmTMB(betweenness  ~                                dilution.rescaled + sex + (1|trial),                   data = network_dilution, family = "poisson", ziformula = ~1)
clusteringcoeff_GxE_dilution <- lmer(clusteringcoeff ~ (dilution.rescaled|genotype) + sex*dilution.rescaled   + (1|trial),                   data = network_dilution)
clusteringcoeff_G_dilution   <- lmer(clusteringcoeff ~ (1|genotype)                 + sex*dilution.rescaled   + (1|trial),                   data = network_dilution)
clusteringcoeff_NA_dilution  <- lmer(clusteringcoeff ~                                sex*dilution.rescaled   + (1|trial),                   data = network_dilution)
# Observed significance of GxDilution
anova(indegree_GxE_dilution, indegree_G_dilution)                             # Observed significance of GxDilution:            LR = 1.2601,     Df = 2,  p = 0.5326
anova(outdegree_GxE_dilution, outdegree_G_dilution)                           # Observed significance of GxDilution:            LR = 0.8204,     Df = 2,  p = 0.6635
anova(eigencentrality_GxE_dilution, eigencentrality_G_dilution)               # Observed significance of GxDilution:            LR = 0.5369,     Df = 2,  p = 0.7645
anova(betweenness_GxE_dilution, betweenness_G_dilution)                       # Observed significance of GxDilution:            LR = 565.97,     Df = 2,  p < 2.2e-16
anova(clusteringcoeff_GxE_dilution, clusteringcoeff_G_dilution)               # Observed significance of GxDilution:            LR = 1.0472,     Df = 2,  p = 0.5924
# Observed significance of G
anova(indegree_G_dilution, indegree_NA_dilution)                              # Observed significance of G:                     LR = 15.992,     Df = 1,  p = 6.36e-05
anova(outdegree_G_dilution, outdegree_NA_dilution)                            # Observed significance of G:                     LR = 11.284,     Df = 1,  p = 0.0007816
anova(eigencentrality_G_dilution, eigencentrality_NA_dilution)                # Observed significance of G:                     LR = 12.069,     Df = 1,  p = 0.0005128
anova(betweenness_G_dilution, betweenness_NA_dilution)                        # Observed significance of G:                     LR = 611.21,     Df = 1,  p < 2.2e-16
anova(clusteringcoeff_G_dilution, clusteringcoeff_NA_dilution)                # Observed significance of G:                     LR = 32.614,     Df = 1,  p = 1.124e-08
# Observed significance of Dilution
Anova(indegree_G_dilution, type = "3")                                        # Observed significance of Dilution:              Chisq = 1.1299,  Df = 1,  p = 0.28781
Anova(outdegree_G_dilution, type = "3")                                       # Observed significance of Dilution:              Chisq = 2.4621,  Df = 1,  p = 0.11662
fligner.test(eigencentrality ~ dilution.rescaled, data = network_dilution)    # Observed significance of Dilution on variation: Chisq = 8.0973,  Df = 2,  p = 0.01745
Anova(betweenness_G_dilution, type = "3")                                     # Observed significance of Dilution:              Chisq = 0.3414,  Df = 1,  p = 0.55904
Anova(clusteringcoeff_G_dilution, type = "3")                                 # Observed significance of Dilution:              Chisq = 11.618,  Df = 1,  p = 0.0006532
# Observed significance of Sex
Anova(indegree_G_dilution, type = "3")                                        # Observed significance of Sex:                   Chisq = 2.8539,  Df = 1,  p = 0.09115
Anova(outdegree_G_dilution, type = "3")                                       # Observed significance of Sex:                   Chisq = 0.4647,  Df = 1,  p = 0.49542
Anova(eigencentrality_G_dilution, type = "3")                                 # Observed significance of Sex:                   Chisq = 3.4682,  Df = 1,  p = 0.06256
Anova(betweenness_G_dilution, type = "3")                                     # Observed significance of Sex:                   Chisq = 6.0741,  Df = 1,  p = 0.01372
Anova(clusteringcoeff_G_dilution, type = "3")                                 # Observed significance of Sex:                   Chisq = 16.483,  Df = 1,  p = 4.908e-05
# Save observed test statistics as values
indegree_GxE_dilution_obs        <- anova(indegree_GxE_dilution, indegree_G_dilution)$Chisq[2]
outdegree_GxE_dilution_obs       <- anova(outdegree_GxE_dilution, outdegree_G_dilution)$Chisq[2]
eigencentrality_GxE_dilution_obs <- anova(eigencentrality_GxE_dilution, eigencentrality_G_dilution)$Chisq[2]
betweenness_GxE_dilution_obs     <- anova(betweenness_GxE_dilution, betweenness_G_dilution)$Chisq[2]
clusteringcoeff_GxE_dilution_obs <- anova(clusteringcoeff_GxE_dilution, clusteringcoeff_G_dilution)$Chisq[2]
indegree_G_dilution_obs          <- anova(indegree_G_dilution, indegree_NA_dilution)$Chisq[2]
outdegree_G_dilution_obs         <- anova(outdegree_G_dilution, outdegree_NA_dilution)$Chisq[2]
eigencentrality_G_dilution_obs   <- anova(eigencentrality_G_dilution, eigencentrality_NA_dilution)$Chisq[2]
betweenness_G_dilution_obs       <- anova(betweenness_G_dilution, betweenness_NA_dilution)$Chisq[2]
clusteringcoeff_G_dilution_obs   <- anova(clusteringcoeff_G_dilution, clusteringcoeff_NA_dilution)$Chisq[2]
indegree_Sex_dilution_obs        <- Anova(indegree_G_dilution, type = "3")$Chisq[3]
outdegree_Sex_dilution_obs       <- Anova(outdegree_G_dilution, type = "3")$Chisq[2]
eigencentrality_Sex_dilution_obs <- Anova(eigencentrality_G_dilution, type = "3")$Chisq[2]
betweenness_Sex_dilution_obs     <- Anova(betweenness_G_dilution, type = "3")$Chisq[3]
clusteringcoeff_Sex_dilution_obs <- Anova(clusteringcoeff_G_dilution, type = "3")$Chisq[2]

# PERMUTED G*E, G, AND SEX EFFECTS ON NETWORK POSITIONS (FOR DILUTION ONLY) # ------------------
indegree_GxE_dilution_null        = numeric(1000)
indegree_G_dilution_null          = numeric(1000)
indegree_Sex_dilution_null        = numeric(1000)
outdegree_GxE_dilution_null       = numeric(1000)
outdegree_G_dilution_null         = numeric(1000)
outdegree_Sex_dilution_null       = numeric(1000)
eigencentrality_GxE_dilution_null = numeric(1000)
eigencentrality_G_dilution_null   = numeric(1000)
eigencentrality_Sex_dilution_null = numeric(1000)
betweenness_GxE_dilution_null     = numeric(1000)
betweenness_G_dilution_null       = numeric(1000)
betweenness_Sex_dilution_null     = numeric(1000)
clusteringcoeff_GxE_dilution_null = numeric(1000)
clusteringcoeff_G_dilution_null   = numeric(1000)
clusteringcoeff_Sex_dilution_null = numeric(1000)
for (i in 1:1000) {
  network_dilution_rand = network_dilution %>% 
    group_by(video) %>%
    mutate(sex_rand = sample(sex, replace=F)) %>%
    group_by(video, sex) %>%
    mutate(genotype_rand = sample(genotype, replace=F))
  indegree_GxE_dilution_rand        <- glmer(indegree       ~ (dilution.rescaled|genotype_rand) + dilution.rescaled + sex      + (1|trial) + (1|observation), data = network_dilution_rand, family = poisson)
  indegree_G_dilution_rand          <- glmer(indegree       ~ (1|genotype_rand)                 + dilution.rescaled + sex      + (1|trial) + (1|observation), data = network_dilution_rand, family = poisson)
  indegree_NA_dilution_rand         <- glmer(indegree       ~                                     dilution.rescaled + sex      + (1|trial) + (1|observation), data = network_dilution_rand, family = poisson)
  indegree_Sex_dilution_rand        <- glmer(indegree       ~ (1|genotype)                      + dilution.rescaled + sex_rand + (1|trial) + (1|observation), data = network_dilution_rand, family = poisson)
  outdegree_GxE_dilution_rand       <- glmer(outdegree      ~ (dilution.rescaled|genotype_rand) + sex*dilution.rescaled        + (1|trial) + (1|observation), data = network_dilution_rand, family = poisson)
  outdegree_G_dilution_rand         <- glmer(outdegree      ~ (1|genotype_rand)                 + sex*dilution.rescaled        + (1|trial) + (1|observation), data = network_dilution_rand, family = poisson)
  outdegree_NA_dilution_rand        <- glmer(outdegree      ~                                     sex*dilution.rescaled        + (1|trial) + (1|observation), data = network_dilution_rand, family = poisson)
  outdegree_Sex_dilution_rand       <- glmer(outdegree      ~ (1|genotype)                      + sex_rand*dilution.rescaled   + (1|trial) + (1|observation), data = network_dilution_rand, family = poisson)
  eigencentrality_GxE_dilution_rand <- lmer(eigencentrality ~ (dilution.rescaled|genotype_rand) + sex                          + (1|individual),              data = network_dilution_rand)
  eigencentrality_G_dilution_rand   <- lmer(eigencentrality ~ (1|genotype_rand)                 + sex                          + (1|individual),              data = network_dilution_rand)
  eigencentrality_NA_dilution_rand  <- lmer(eigencentrality ~                                     sex                          + (1|individual),              data = network_dilution_rand)
  eigencentrality_Sex_dilution_rand <- lmer(eigencentrality ~ (1|genotype)                      + sex_rand                     + (1|individual),              data = network_dilution_rand)
  betweenness_GxE_dilution_rand     <- glmmTMB(betweenness  ~ (dilution.rescaled|genotype_rand) + dilution.rescaled + sex      + (1|trial),                   data = network_dilution_rand, family = "poisson", ziformula = ~1)
  betweenness_G_dilution_rand       <- glmmTMB(betweenness  ~ (1|genotype_rand)                 + dilution.rescaled + sex      + (1|trial),                   data = network_dilution_rand, family = "poisson", ziformula = ~1)
  betweenness_NA_dilution_rand      <- glmmTMB(betweenness  ~                                     dilution.rescaled + sex      + (1|trial),                   data = network_dilution_rand, family = "poisson", ziformula = ~1)
  betweenness_Sex_dilution_rand     <- glmmTMB(betweenness  ~ (1|genotype)                      + dilution.rescaled + sex_rand + (1|trial),                   data = network_dilution_rand, family = "poisson", ziformula = ~1)
  clusteringcoeff_GxE_dilution_rand <- lmer(clusteringcoeff ~ (dilution.rescaled|genotype_rand) + sex*dilution.rescaled        + (1|trial),                   data = network_dilution_rand)
  clusteringcoeff_G_dilution_rand   <- lmer(clusteringcoeff ~ (1|genotype_rand)                 + sex*dilution.rescaled        + (1|trial),                   data = network_dilution_rand)
  clusteringcoeff_NA_dilution_rand  <- lmer(clusteringcoeff ~                                     sex*dilution.rescaled        + (1|trial),                   data = network_dilution_rand)
  clusteringcoeff_Sex_dilution_rand <- lmer(clusteringcoeff ~ (1|genotype)                      + sex_rand*dilution.rescaled   + (1|trial),                   data = network_dilution_rand)
  indegree_GxE_dilution_null[i]        <- anova(indegree_GxE_dilution_rand, indegree_G_dilution_rand)$Chisq[2]
  outdegree_GxE_dilution_null[i]       <- anova(outdegree_GxE_dilution_rand, outdegree_G_dilution_rand)$Chisq[2]
  eigencentrality_GxE_dilution_null[i] <- anova(eigencentrality_GxE_dilution_rand, eigencentrality_G_dilution_rand)$Chisq[2]
  betweenness_GxE_dilution_null[i]     <- anova(betweenness_GxE_dilution_rand, betweenness_G_dilution_rand)$Chisq[2]
  clusteringcoeff_GxE_dilution_null[i] <- anova(clusteringcoeff_GxE_dilution_rand, clusteringcoeff_G_dilution_rand)$Chisq[2]
  indegree_G_dilution_null[i]          <- anova(indegree_G_dilution_rand, indegree_NA_dilution_rand)$Chisq[2]
  outdegree_G_dilution_null[i]         <- anova(outdegree_G_dilution_rand, outdegree_NA_dilution_rand)$Chisq[2]
  eigencentrality_G_dilution_null[i]   <- anova(eigencentrality_G_dilution_rand, eigencentrality_NA_dilution_rand)$Chisq[2]
  betweenness_G_dilution_null[i]       <- anova(betweenness_G_dilution_rand, betweenness_NA_dilution_rand)$Chisq[2]
  clusteringcoeff_G_dilution_null[i]   <- anova(clusteringcoeff_G_dilution_rand, clusteringcoeff_NA_dilution_rand)$Chisq[2]
  indegree_Sex_dilution_null[i]        <- Anova(indegree_Sex_dilution_rand, type = "3")$Chisq[3]
  outdegree_Sex_dilution_null[i]       <- Anova(outdegree_Sex_dilution_rand, type = "3")$Chisq[2]
  eigencentrality_Sex_dilution_null[i] <- Anova(eigencentrality_Sex_dilution_rand, type = "3")$Chisq[2]
  betweenness_Sex_dilution_null[i]     <- Anova(betweenness_Sex_dilution_rand, type = "3")$Chisq[3]
  clusteringcoeff_Sex_dilution_null[i] <- Anova(clusteringcoeff_Sex_dilution_rand, type = "3")$Chisq[2]}
# Permuted significance of GxDilution
sum(indegree_GxE_dilution_obs<indegree_GxE_dilution_null)/1000                      # Permuted significance of GxDilution: p = 0.082
sum(outdegree_GxE_dilution_obs<outdegree_GxE_dilution_null)/1000                    # Permuted significance of GxDilution: p = 0.109
sum(eigencentrality_GxE_dilution_obs<eigencentrality_GxE_dilution_null)/1000        # Permuted significance of GxDilution: p = 0.251
sum(betweenness_GxE_dilution_obs<betweenness_GxE_dilution_null)/1000                # Permuted significance of GxDilution: p = 0.782
sum(clusteringcoeff_GxE_dilution_obs<clusteringcoeff_GxE_dilution_null)/1000        # Permuted significance of GxDilution: p = 0.028
# Permuted significance of G
sum(indegree_G_dilution_obs<indegree_G_dilution_null)/1000                          # Permuted significance of G:          p < 0.001
sum(outdegree_G_dilution_obs<outdegree_G_dilution_null)/1000                        # Permuted significance of G:          p < 0.001
sum(eigencentrality_G_dilution_obs<eigencentrality_G_dilution_null)/1000            # Permuted significance of G:          p < 0.001
sum(betweenness_G_dilution_obs<betweenness_G_dilution_null)/1000                    # Permuted significance of G:          p = 0.554
sum(clusteringcoeff_G_dilution_obs<clusteringcoeff_G_dilution_null)/1000            # Permuted significance of G:          p < 0.001
# Permuted significance of Sex
sum(indegree_Sex_dilution_obs<indegree_Sex_dilution_null)/1000                      # Permuted significance of Sex:        p = 0.087
sum(outdegree_Sex_dilution_obs<outdegree_Sex_dilution_null)/1000                    # Permuted significance of Sex:        p = 0.435
sum(eigencentrality_Sex_dilution_obs<eigencentrality_Sex_dilution_null)/1000        # Permuted significance of Sex:        p = 0.072
sum(betweenness_Sex_dilution_obs<betweenness_Sex_dilution_null)/1000                # Permuted significance of Sex:        p = 0.687
sum(clusteringcoeff_Sex_dilution_obs<clusteringcoeff_Sex_dilution_null)/1000        # Permuted significance of Sex:        p < 0.001

# BROAD-SENSE HERITABILITY ESTIMATES OF NODE-LEVEL NETWORK POSITIONS # --------------
# Estimates the total variance attributable to the random effect "genotype" using Nakagawa, Johnson, and Schielzeth's MuMIn package
# Note that as of June 2020, the r.squaredGLMM function in MuMIn cannot handle zero-inflated models (as is the case for betweenness).
# As such, betweenness was modeled with a poisson distribution and OLRE to account for overdispersion.
r.squaredGLMM(indegree_G, pj2014 = TRUE)[5] - r.squaredGLMM(indegree_NA, pj2014 = TRUE)[5]                                          # 0.02362979
r.squaredGLMM(outdegree_G, pj2014 = TRUE)[5] - r.squaredGLMM(outdegree_NA, pj2014 = TRUE)[5]                                        # 0.02504123
r.squaredGLMM(eigencentrality_G)[2] - r.squaredGLMM(eigencentrality_NA)[2]                                                          # 0.04191849
r.squaredGLMM(glmer(betweenness~(1|genotype)+food+sex+(1|trial)+(1|observation),data=network,family=poisson), pj2014 = TRUE)[5] +
  - r.squaredGLMM(glmer(betweenness~food+sex+(1|trial)+(1|observation),data=network,family=poisson), pj2014 = TRUE)[5]              # 0.1659185
r.squaredGLMM(clusteringcoeff_G)[2] - r.squaredGLMM(clusteringcoeff_NA)[2]                                                          # 0.04999238

# CONSISTENCY OF NETWORK POSITIONS ACROSS DAYS # --------------
bothdays <- network[network$individual %in% network$individual[duplicated(network$individual)],] # Reduce dataframe to only cases where both days 1 and 2 were measured for a social group
bothdays_indegree <- dcast(bothdays, trial + genotype ~ videoday, value.var = "indegree") # Cast dataframe so that each day is its own column for each network position measure
bothdays_outdegree <- dcast(bothdays, trial + genotype ~ videoday, value.var = "outdegree")
bothdays_eigencentrality <- dcast(bothdays, trial + genotype ~ videoday, value.var = "eigencentrality")
bothdays_betweenness <- dcast(bothdays, trial + genotype ~ videoday, value.var = "betweenness")
bothdays_clusteringcoeff <- dcast(bothdays, trial + genotype ~ videoday, value.var = "clusteringcoeff")
bothdays_indegree <- bothdays_indegree %>% rename(indegree_day1 = "1", indegree_day2 = "2")
bothdays_outdegree <- bothdays_outdegree %>% rename(outdegree_day1 = "1", outdegree_day2 = "2")
bothdays_eigencentrality <- bothdays_eigencentrality %>% rename(eigencentrality_day1 = "1", eigencentrality_day2 = "2")
bothdays_betweenness <- bothdays_betweenness %>% rename(betweenness_day1 = "1", betweenness_day2 = "2")
bothdays_clusteringcoeff <- bothdays_clusteringcoeff %>% rename(clusteringcoeff_day1 = "1", clusteringcoeff_day2 = "2")
bothdays1 <- merge(bothdays_indegree, bothdays_outdegree, by = c("trial", "genotype"))
bothdays2 <- merge(bothdays_eigencentrality, bothdays_betweenness, by = c("trial", "genotype"))
bothdays3 <- merge(bothdays2, bothdays_clusteringcoeff, by = c("trial", "genotype"))
bothdays <- merge(bothdays1, bothdays3, by = c("trial", "genotype"))
rm(bothdays_indegree, bothdays_outdegree, bothdays_eigencentrality, bothdays_betweenness, bothdays_clusteringcoeff, bothdays1, bothdays2, bothdays3)
plot(bothdays$indegree_day1, bothdays$indegree_day2)
plot(bothdays$outdegree_day1, bothdays$outdegree_day2)
plot(bothdays$eigencentrality_day1, bothdays$eigencentrality_day2)
plot(bothdays$betweenness_day1, bothdays$betweenness_day2)
plot(bothdays$clusteringcoeff_day1, bothdays$clusteringcoeff_day2)
cor.test(bothdays$indegree_day1, bothdays$indegree_day2, method = "kendall")                  # tau = 0.1672161,  z = 4.0172,  p-value = 5.888e-05
cor.test(bothdays$outdegree_day1, bothdays$outdegree_day2, method = "kendall")                # tau = 0.2221694,  z = 5.3373,  p-value = 9.434e-08
cor.test(bothdays$eigencentrality_day1, bothdays$eigencentrality_day2, method = "kendall")    # tau = 0.09070389, z = 2.1792,  p-value = 0.02932
cor.test(bothdays$betweenness_day1, bothdays$betweenness_day2, method = "kendall")            # tau = 0.02381852, z = 0.52447, p-value = 0.6
cor.test(bothdays$clusteringcoeff_day1, bothdays$clusteringcoeff_day2, method = "kendall")    # tau = 0.4581527,  z = 11.007,  p-value < 2.2e-16

# EFFECTS ON FITNESS # ---------------
# Make random effects factors to work with glmmadmb (doesn't change output of lme4 models)
flydata$trial <- as.factor(flydata$trial)
flydata$genotype <- as.factor(flydata$genotype)
# Create dataframes for each fitness variable
# Female offspring and lifespan
fitnessdata                 <- subset(flydata, offspring >= 0)
# Male matings and latency (includes all males, including those that did not mate)
matingsdata                 <- subset(flydata, matings >= 0)
matingsdata                 <- subset(matingsdata, sex == "M")
# Create dataframes where network data was also collected. 
# Individuals could only have one measure for any given fitness metric, but could have two network position measures taken (in cases where videos of the group were taken on both days 1 and 2)
# Because an individual's network position was consistent over time, we used network measures taken on day 1 in our fitness models
fitnessdata_day1            <- subset(fitnessdata, videoday == 1)
matingsdata_day1            <- subset(matingsdata, videoday == 1)
# Node-level network parameters are non-independent in that they are calculated from the same underlying interaction matrix. 
# As such, they are commonly highly correlated, which presents colinearity problems in multiple regression analyses.
# Test how correlated these node-level network parameters are, and if it poses a problem for analyses.
# Test normality of each node-level network parameter to determine whether to use Pearson, Kendall, or Spearman corrleations
shapiro.test(flydata$indegree)            # Non-normal
shapiro.test(flydata$outdegree)           # Non-normal
shapiro.test(flydata$eigencentrality)     # Non-normal
shapiro.test(flydata$betweenness)         # Non-normal
shapiro.test(flydata$clusteringcoeff)     # Non-normal
# Shapiro tests can be prone to over-identifying non-normality when sample sizes are large
# Visualize node-level network parameter distributions
hist(flydata$indegree, breaks = 100)            # Long right tail
hist(flydata$outdegree, breaks = 100)           # Long right tail
hist(flydata$eigencentrality, breaks = 100)     # Appears normal
hist(flydata$betweenness, breaks = 100)         # Zero-inflated
hist(flydata$clusteringcoeff, breaks = 100)     # Appears normal
# Test for correlations using Kendall rank correlations
cor.test(flydata$indegree, flydata$outdegree, method = "kendall")                 # Correlated
cor.test(flydata$indegree, flydata$eigencentrality, method = "kendall")           # Correlated
cor.test(flydata$indegree, flydata$betweenness, method = "kendall")               # Correlated
cor.test(flydata$indegree, flydata$clusteringcoeff, method = "kendall")           # Correlated
cor.test(flydata$outdegree, flydata$eigencentrality, method = "kendall")          # Correlated
cor.test(flydata$outdegree, flydata$betweenness, method = "kendall")              # Correlated
cor.test(flydata$outdegree, flydata$clusteringcoeff, method = "kendall")          # Correlated
cor.test(flydata$eigencentrality, flydata$betweenness, method = "kendall")        # Correlated
cor.test(flydata$eigencentrality, flydata$clusteringcoeff, method = "kendall")    # Correlated
cor.test(flydata$betweenness, flydata$clusteringcoeff, method = "kendall")        # Correlated
# Test whether correlations between node-level network parameters lead to egregious variance inflation factors (>5-10) for a sample model
vif(glmer(offspring ~ (1|genotype) + (1|trial) + indegree.rescaled*food + outdegree.rescaled*food + eigencentrality.rescaled*food + betweenness.rescaled*food + clusteringcoeff.rescaled*food, data = fitnessdata_day1, family = poisson))
vif(lmer(lifespan   ~ (1|genotype) + (1|trial) + indegree.rescaled*food + outdegree.rescaled*food + eigencentrality.rescaled*food + betweenness.rescaled*food + clusteringcoeff.rescaled*food, data = fitnessdata_day1))
vif(glmer(matings   ~ (1|genotype) + (1|trial) + indegree.rescaled*food + outdegree.rescaled*food + eigencentrality.rescaled*food + betweenness.rescaled*food + clusteringcoeff.rescaled*food, data = matingsdata_day1, family = poisson))
vif(lmer(latency    ~ (1|genotype) + (1|trial) + indegree.rescaled*food + outdegree.rescaled*food + eigencentrality.rescaled*food + betweenness.rescaled*food + clusteringcoeff.rescaled*food, data = matingsdata_day1))
# Variance inflation factors for node-level network parameters as high as...
# offspring: 18.376
# lifespan:  31.163
# matings:   48.246
# latency:   52.981
# Model each fitness response separately for each node-level network parameter. Correct for multiple testing with Bonferroni, though this is likely conservative as the multiply tested node-level network parameters are correlated and not truly independent tests.
# Model the effects of genotype and food on fitness responses in models on their own, as testing the effects of food and genotype in every separate node-level network parameter model is redundant.
# Also, many measures of fitness don't have a corresponding network predictor (in cases where individuals died before videos could be taken, or videos were unusable due to various reasons).
# Modeling the effects of genotype and food on fitness separately allows for a more robust sample size.
# Before proceeding, remove duplicate values from the fitness dataframes, in cases where individuals have network data for both days 1 and 2.
fitnessdata <- fitnessdata %>% distinct(trial, genotype, .keep_all = TRUE)
matingsdata <- matingsdata %>% distinct(trial, genotype, .keep_all = TRUE)

# OFFSPRING # --------------
offspring_G         <- glmer(offspring ~ (1|genotype) + (1|trial) + food,        data = fitnessdata, family = poisson)
offspring_NA        <- glmer(offspring ~                (1|trial) + food,        data = fitnessdata, family = poisson)
offspring_G_sim     <- simulateResiduals(fittedModel = offspring_G, n = 1000)
offspring_NA_sim    <- simulateResiduals(fittedModel = offspring_NA, n = 1000)
testResiduals(offspring_G_sim)  # Overdispersion
testResiduals(offspring_NA_sim) # Overdispersion
offspring_G         <- glmmadmb(offspring ~ (1|genotype) + (1|trial) + food,        data = fitnessdata, family = "nbinom")
offspring_NA        <- glmmadmb(offspring ~                (1|trial) + food,        data = fitnessdata, family = "nbinom")
anova(offspring_G, offspring_NA)  # Observed significance of G:   LR = 0.04,     Df = 1,  p = 0.8415
Anova(offspring_G, type = "3")    # Observed significance of E:   LR = 6.84,     Df = 4,  p = 0.1446
offspring_day1_indegreeXfood         <- glmmadmb(offspring ~ (1|genotype) + (1|trial) + indegree.rescaled*food,        data = fitnessdata_day1, family = "nbinom")
offspring_day1_outdegreeXfood        <- glmmadmb(offspring ~ (1|genotype) + (1|trial) + outdegree.rescaled*food,       data = fitnessdata_day1, family = "nbinom")
offspring_day1_eigencentralityXfood  <- glmmadmb(offspring ~ (1|genotype) + (1|trial) + eigencentrality.rescaled*food, data = fitnessdata_day1, family = "nbinom")
offspring_day1_betweennessXfood      <- glmmadmb(offspring ~ (1|genotype) + (1|trial) + betweenness.rescaled*food,     data = fitnessdata_day1, family = "nbinom")
offspring_day1_clusteringcoeffXfood  <- glmmadmb(offspring ~ (1|genotype) + (1|trial) + clusteringcoeff.rescaled*food, data = fitnessdata_day1, family = "nbinom")
Anova(offspring_day1_indegreeXfood, type = "3")          # Observed significance of Indegree*Food on Day 1:          Chisq = 7.0738, Df = 4, p = 0.13204
Anova(offspring_day1_outdegreeXfood, type = "3")         # Observed significance of Outdegree*Food on Day 1:         Chisq = 4.1296, Df = 4, p = 0.38875
Anova(offspring_day1_eigencentralityXfood, type = "3")   # Observed significance of Eigencentrality*Food on Day 1:   Chisq = 2.2688, Df = 4, p = 0.6865
Anova(offspring_day1_betweennessXfood, type = "3")       # Observed significance of Betweenness*Food on Day 1:       Chisq = 2.9432, Df = 4, p = 0.5674
Anova(offspring_day1_clusteringcoeffXfood, type = "3")   # Observed significance of Clustering Coeff*Food on Day 1:  Chisq = 0.8829, Df = 4, p = 0.9270
offspring_day1_indegree         <- glmmadmb(offspring ~ (1|genotype) + (1|trial) + indegree.rescaled + food,        data = fitnessdata_day1, family = "nbinom")
offspring_day1_outdegree        <- glmmadmb(offspring ~ (1|genotype) + (1|trial) + outdegree.rescaled + food,       data = fitnessdata_day1, family = "nbinom")
offspring_day1_eigencentrality  <- glmmadmb(offspring ~ (1|genotype) + (1|trial) + eigencentrality.rescaled + food, data = fitnessdata_day1, family = "nbinom")
offspring_day1_betweenness      <- glmmadmb(offspring ~ (1|genotype) + (1|trial) + betweenness.rescaled + food,     data = fitnessdata_day1, family = "nbinom")
offspring_day1_clusteringcoeff  <- glmmadmb(offspring ~ (1|genotype) + (1|trial) + clusteringcoeff.rescaled + food, data = fitnessdata_day1, family = "nbinom")
Anova(offspring_day1_indegree, type = "3")          # Observed significance of Indegree on Day 1:          Chisq = 0.8571, Df = 1, p = 0.3546
Anova(offspring_day1_outdegree, type = "3")         # Observed significance of Outdegree on Day 1:         Chisq = 0.5697, Df = 1, p = 0.4504
Anova(offspring_day1_eigencentrality, type = "3")   # Observed significance of Eigencentrality on Day 1:   Chisq = 1.4058, Df = 1, p = 0.4387
Anova(offspring_day1_betweenness, type = "3")       # Observed significance of Betweenness on Day 1:       Chisq = 1.5934, Df = 1, p = 0.2068
Anova(offspring_day1_clusteringcoeff, type = "3")   # Observed significance of Clustering Coeff on Day 1:  Chisq = 0.7173, Df = 1, p = 0.3970

# LIFESPAN # --------------
lifespan_G         <- glmer(lifespan ~ (1|genotype) + (1|trial) + food,        data = fitnessdata, family = poisson)
lifespan_NA        <- glmer(lifespan ~                (1|trial) + food,        data = fitnessdata, family = poisson)
lifespan_G_sim     <- simulateResiduals(fittedModel = lifespan_G, n = 1000)
lifespan_NA_sim    <- simulateResiduals(fittedModel = lifespan_NA, n = 1000)
testResiduals(lifespan_G_sim)  # Overdispersion
testResiduals(lifespan_NA_sim) # Overdispersion
lifespan_G         <- glmmadmb(lifespan ~ (1|genotype) + (1|trial) + food,        data = fitnessdata, family = "nbinom")
lifespan_NA        <- glmmadmb(lifespan ~                (1|trial) + food,        data = fitnessdata, family = "nbinom")
anova(lifespan_G, lifespan_NA)   # Observed significance of G:   LR = 0.16,     Df = 1,  p = 0.6892
Anova(lifespan_G, type = "3")    # Observed significance of E:   LR = 27.277,   Df = 4,  p = 1.747e-05
lifespan_day1_indegreeXfood         <- glmmadmb(lifespan ~ (1|genotype) + (1|trial) + indegree.rescaled*food,        data = fitnessdata_day1, family = "nbinom")
lifespan_day1_outdegreeXfood        <- glmmadmb(lifespan ~ (1|genotype) + (1|trial) + outdegree.rescaled*food,       data = fitnessdata_day1, family = "nbinom")
lifespan_day1_eigencentralityXfood  <- glmmadmb(lifespan ~ (1|genotype) + (1|trial) + eigencentrality.rescaled*food, data = fitnessdata_day1, family = "nbinom")
lifespan_day1_betweennessXfood      <- glmmadmb(lifespan ~ (1|genotype) + (1|trial) + betweenness.rescaled*food,     data = fitnessdata_day1, family = "nbinom")
lifespan_day1_clusteringcoeffXfood  <- glmmadmb(lifespan ~ (1|genotype) + (1|trial) + clusteringcoeff.rescaled*food, data = fitnessdata_day1, family = "nbinom")
Anova(lifespan_day1_indegreeXfood, type = "3")          # Observed significance of Indegree*Food on Day 1:          Chisq = 1.0173, Df = 4, p = 0.9072
Anova(lifespan_day1_outdegreeXfood, type = "3")         # Observed significance of Outdegree*Food on Day 1:         Chisq = 2.1840, Df = 4, p = 0.7020
Anova(lifespan_day1_eigencentralityXfood, type = "3")   # Observed significance of Eigencentrality*Food on Day 1:   Chisq = 0.3940, Df = 4, p = 0.9830
Anova(lifespan_day1_betweennessXfood, type = "3")       # Observed significance of Betweenness*Food on Day 1:       Chisq = 5.3711, Df = 4, p = 0.2512982
Anova(lifespan_day1_clusteringcoeffXfood, type = "3")   # Observed significance of Clustering Coeff*Food on Day 1:  Chisq = 1.2195, Df = 4, p = 0.8749
lifespan_day1_indegree         <- glmmadmb(lifespan ~ (1|genotype) + (1|trial) + indegree.rescaled + food,        data = fitnessdata_day1, family = "nbinom")
lifespan_day1_outdegree        <- glmmadmb(lifespan ~ (1|genotype) + (1|trial) + outdegree.rescaled + food,       data = fitnessdata_day1, family = "nbinom")
lifespan_day1_eigencentrality  <- glmmadmb(lifespan ~ (1|genotype) + (1|trial) + eigencentrality.rescaled + food, data = fitnessdata_day1, family = "nbinom")
lifespan_day1_betweenness      <- glmmadmb(lifespan ~ (1|genotype) + (1|trial) + betweenness.rescaled + food,     data = fitnessdata_day1, family = "nbinom")
lifespan_day1_clusteringcoeff  <- glmmadmb(lifespan ~ (1|genotype) + (1|trial) + clusteringcoeff.rescaled + food, data = fitnessdata_day1, family = "nbinom")
Anova(lifespan_day1_indegree, type = "3")          # Observed significance of Indegree on Day 1:          Chisq = 0.0065, Df = 1, p = 0.935522
Anova(lifespan_day1_outdegree, type = "3")         # Observed significance of Outdegree on Day 1:         Chisq = 0.0147, Df = 1, p = 0.90339
Anova(lifespan_day1_eigencentrality, type = "3")   # Observed significance of Eigencentrality on Day 1:   Chisq = 1.1499, Df = 1, p = 0.283576
Anova(lifespan_day1_betweenness, type = "3")       # Observed significance of Betweenness on Day 1:       Chisq = 2.1379, Df = 1, p = 0.143702
Anova(lifespan_day1_clusteringcoeff, type = "3")   # Observed significance of Clustering Coeff on Day 1:  Chisq = 0.000,  Df = 1, p = 0.99558

# MATINGS # --------------
matings_G         <- glmer(matings ~ (1|genotype) + (1|trial) + food,        data = matingsdata, family = poisson)
matings_NA        <- glmer(matings ~                (1|trial) + food,        data = matingsdata, family = poisson)
matings_G_sim     <- simulateResiduals(fittedModel = matings_G, n = 1000)
matings_NA_sim    <- simulateResiduals(fittedModel = matings_NA, n = 1000)
testResiduals(matings_G_sim)  # 0 inflation
testResiduals(matings_NA_sim) # 0 inflation
matings_G         <- glmmadmb(matings ~ (1|genotype) + (1|trial) + food,        data = matingsdata, family = "poisson", zeroInflation = TRUE)
matings_NA        <- glmmadmb(matings ~                (1|trial) + food,        data = matingsdata, family = "poisson", zeroInflation = TRUE)
anova(matings_G, matings_NA)    # Observed significance of G:   LR = 99.76,    Df = 1,  p < 2.2e-16
Anova(matings_G, type = "3")    # Observed significance of E:   LR = 1.0208,   Df = 4,  p = 0.9066
matings_day1_indegreeXfood         <- glmmadmb(matings ~ (1|genotype) + (1|trial) + indegree.rescaled*food,        data = matingsdata_day1, family = "poisson", zeroInflation = TRUE)
matings_day1_outdegreeXfood        <- glmmadmb(matings ~ (1|genotype) + (1|trial) + outdegree.rescaled*food,       data = matingsdata_day1, family = "poisson", zeroInflation = TRUE)
matings_day1_eigencentralityXfood  <- glmmadmb(matings ~ (1|genotype) + (1|trial) + eigencentrality.rescaled*food, data = matingsdata_day1, family = "poisson", zeroInflation = TRUE)
matings_day1_betweennessXfood      <- glmmadmb(matings ~ (1|genotype) + (1|trial) + betweenness.rescaled*food,     data = matingsdata_day1, family = "poisson", zeroInflation = TRUE)
matings_day1_clusteringcoeffXfood  <- glmmadmb(matings ~ (1|genotype) + (1|trial) + clusteringcoeff.rescaled*food, data = matingsdata_day1, family = "poisson", zeroInflation = TRUE)
Anova(matings_day1_indegreeXfood, type = "3")          # Observed significance of Indegree*Food on Day 1:          Chisq = 52.0140, Df = 4, p = 1.370e-10
Anova(matings_day1_outdegreeXfood, type = "3")         # Observed significance of Outdegree*Food on Day 1:         Chisq = 20.4397, Df = 4, p = 0.0004088
Anova(matings_day1_eigencentralityXfood, type = "3")   # Observed significance of Eigencentrality*Food on Day 1:   Chisq = 49.6988, Df = 4, p = 4.174e-10
Anova(matings_day1_betweennessXfood, type = "3")       # Observed significance of Betweenness*Food on Day 1:       Chisq = 7.4725,  Df = 4, p = 0.11293
Anova(matings_day1_clusteringcoeffXfood, type = "3")   # Observed significance of Clustering Coeff*Food on Day 1:  Chisq = 5.8754,  Df = 4, p = 0.2087
matings_day1_indegree         <- glmmadmb(matings ~ (1|genotype) + (1|trial) + indegree.rescaled + food,        data = matingsdata_day1, family = "poisson", zeroInflation = TRUE)
matings_day1_outdegree        <- glmmadmb(matings ~ (1|genotype) + (1|trial) + outdegree.rescaled + food,       data = matingsdata_day1, family = "poisson", zeroInflation = TRUE)
matings_day1_eigencentrality  <- glmmadmb(matings ~ (1|genotype) + (1|trial) + eigencentrality.rescaled + food, data = matingsdata_day1, family = "poisson", zeroInflation = TRUE)
matings_day1_betweenness      <- glmmadmb(matings ~ (1|genotype) + (1|trial) + betweenness.rescaled + food,     data = matingsdata_day1, family = "poisson", zeroInflation = TRUE)
matings_day1_clusteringcoeff  <- glmmadmb(matings ~ (1|genotype) + (1|trial) + clusteringcoeff.rescaled + food, data = matingsdata_day1, family = "poisson", zeroInflation = TRUE)
Anova(matings_day1_indegreeXfood, type = "3")          # Observed significance of Indegree on Day 1:          Chisq = 8.7949, Df = 1, p = 0.003021
Anova(matings_day1_outdegreeXfood, type = "3")         # Observed significance of Outdegree on Day 1:         Chisq = 4.0745, Df = 1, p = 0.0435342
Anova(matings_day1_eigencentralityXfood, type = "3")   # Observed significance of Eigencentrality on Day 1:   Chisq = 9.9375, Df = 1, p = 0.001619
Anova(matings_day1_betweenness, type = "3")            # Observed significance of Betweenness on Day 1:       Chisq = 2.0795, Df = 1, p = 0.1493
Anova(matings_day1_clusteringcoeff, type = "3")        # Observed significance of Clustering Coeff on Day 1:  Chisq = 2.4098, Df = 1, p = 0.1206

# LATENCY # --------------
# Create latency variables to indicate censoring (because matings were only observed up until 120 minutes after social groups were established)
matingsdata$latency.censored <- matingsdata$latency
matingsdata$latency.censored[is.na(matingsdata$latency.censored)] <- 120
matingsdata$is.censored <- ifelse(matingsdata$latency.censored == 120, 0, 1)
# Test proportional hazards assumptions
cox.zph(coxph(Surv(latency.censored, is.censored) ~ food + frailty(trial), data = matingsdata))
# Global p = 0.8191:  Model meets proportional hazards assumptions
latency_G  <- coxme(Surv(latency.censored, is.censored) ~ (1|genotype) + (1|trial) + food, data = matingsdata)
latency_NA <- coxme(Surv(latency.censored, is.censored) ~                (1|trial) + food, data = matingsdata)
anova(latency_G, latency_NA)    # Observed significance of G:   LR = 87.104,  Df = 1,  p < 2.2e-16
Anova(latency_G, type = "3")    # Observed significance of E:   LR = 2.418,   Df = 4,  p = 0.6594
matingsdata_day1$latency.censored <- matingsdata_day1$latency
matingsdata_day1$latency.censored[is.na(matingsdata_day1$latency.censored)] <- 120
matingsdata_day1$is.censored <- ifelse(matingsdata_day1$latency.censored == 120, 0, 1)
latency_day1_indegreeXfood        <- coxme(Surv(latency.censored, is.censored) ~ food*indegree.rescaled        + (1|genotype) + (1|trial), data=matingsdata_day1)
latency_day1_outdegreeXfood       <- coxme(Surv(latency.censored, is.censored) ~ food*outdegree.rescaled       + (1|genotype) + (1|trial), data=matingsdata_day1)
latency_day1_eigencentralityXfood <- coxme(Surv(latency.censored, is.censored) ~ food*eigencentrality.rescaled + (1|genotype) + (1|trial), data=matingsdata_day1)
latency_day1_betweennessXfood     <- coxme(Surv(latency.censored, is.censored) ~ food*betweenness.rescaled     + (1|genotype) + (1|trial), data=matingsdata_day1)
latency_day1_clusteringcoeffXfood <- coxme(Surv(latency.censored, is.censored) ~ food*clusteringcoeff.rescaled + (1|genotype) + (1|trial), data=matingsdata_day1)
Anova(latency_day1_indegreeXfood, type = "3")          # Observed significance of Indegree*Food on Day 1:          Chisq = 13.7200,  Df = 4, p = 0.008244
Anova(latency_day1_outdegreeXfood, type = "3")         # Observed significance of Outdegree*Food on Day 1:         Chisq = 10.6016,  Df = 4, p = 0.03143
Anova(latency_day1_eigencentralityXfood, type = "3")   # Observed significance of Eigencentrality*Food on Day 1:   Chisq = 12.8030,  Df = 4, p = 0.012280
Anova(latency_day1_betweennessXfood, type = "3")       # Observed significance of Betweenness*Food on Day 1:       Chisq = 5.8568,   Df = 4, p = 0.2101
Anova(latency_day1_clusteringcoeffXfood, type = "3")   # Observed significance of Clustering Coeff*Food on Day 1:  Chisq = 4.1280,   Df = 4, p = 0.3890
latency_day1_indegree        <- coxme(Surv(latency.censored, is.censored) ~ food + indegree.rescaled        + (1|genotype) + (1|trial), data=matingsdata_day1)
latency_day1_outdegree       <- coxme(Surv(latency.censored, is.censored) ~ food + outdegree.rescaled       + (1|genotype) + (1|trial), data=matingsdata_day1)
latency_day1_eigencentrality <- coxme(Surv(latency.censored, is.censored) ~ food + eigencentrality.rescaled + (1|genotype) + (1|trial), data=matingsdata_day1)
latency_day1_betweenness     <- coxme(Surv(latency.censored, is.censored) ~ food + betweenness.rescaled     + (1|genotype) + (1|trial), data=matingsdata_day1)
latency_day1_clusteringcoeff <- coxme(Surv(latency.censored, is.censored) ~ food + clusteringcoeff.rescaled + (1|genotype) + (1|trial), data=matingsdata_day1)
Anova(latency_day1_indegreeXfood, type = "3")     # Observed significance of Indegree on Day 1:          Chisq = 9.6871, Df = 1, p = 0.001856
Anova(latency_day1_outdegree, type = "3")         # Observed significance of Outdegree on Day 1:         Chisq = 3.0696, Df = 1, p = 0.07977
Anova(latency_day1_eigencentrality, type = "3")   # Observed significance of Eigencentrality on Day 1:   Chisq = 0.0089, Df = 1, p = 0.9248
Anova(latency_day1_betweenness, type = "3")       # Observed significance of Betweenness on Day 1:       Chisq = 0.5046, Df = 1, p = 0.4775
Anova(latency_day1_clusteringcoeff, type = "3")   # Observed significance of Clustering Coeff on Day 1:  Chisq = 3.8138, Df = 1, p = 0.05083

# OFFSPRING (P:C ONLY) # --------------
fitnessdata_pc <- fitnessdata[which(fitnessdata$food=='1' | fitnessdata$food=='2' | fitnessdata$food=='3'),]
fitnessdata_day1_pc <- fitnessdata_day1[which(fitnessdata_day1$food=='1' | fitnessdata_day1$food=='2' | fitnessdata_day1$food=='3'),]
offspring_G_pc         <- glmer(offspring ~ (1|genotype) + (1|trial) + pc.rescaled,        data = fitnessdata_pc, family = poisson)
offspring_NA_pc        <- glmer(offspring ~                (1|trial) + pc.rescaled,        data = fitnessdata_pc, family = poisson)
offspring_G_pc_sim     <- simulateResiduals(fittedModel = offspring_G_pc, n = 1000)
offspring_NA_pc_sim    <- simulateResiduals(fittedModel = offspring_NA_pc, n = 1000)
testResiduals(offspring_G_pc_sim)  # Overdispersion
testResiduals(offspring_NA_pc_sim) # Overdispersion
offspring_G_pc         <- glmmadmb(offspring ~ (1|genotype) + (1|trial) + pc.rescaled,        data = fitnessdata_pc, family = "nbinom")
offspring_NA_pc        <- glmmadmb(offspring ~                (1|trial) + pc.rescaled,        data = fitnessdata_pc, family = "nbinom")
anova(offspring_G_pc, offspring_NA_pc)  # Observed significance of G:    LR = 0.000,     Df = 1,  p = 1
Anova(offspring_G_pc, type = "3")       # Observed significance of PC:   LR = 4.8603,    Df = 1,  p = 0.02748
offspring_day1_pc_indegreeXpc         <- glmmadmb(offspring ~ (1|genotype) + (1|trial) + indegree.rescaled*pc.rescaled,        data = fitnessdata_day1_pc, family = "nbinom")
offspring_day1_pc_outdegreeXpc        <- glmmadmb(offspring ~ (1|genotype) + (1|trial) + outdegree.rescaled*pc.rescaled,       data = fitnessdata_day1_pc, family = "nbinom")
offspring_day1_pc_eigencentralityXpc  <- glmmadmb(offspring ~ (1|genotype) + (1|trial) + eigencentrality.rescaled*pc.rescaled, data = fitnessdata_day1_pc, family = "nbinom")
offspring_day1_pc_betweennessXpc      <- glmmadmb(offspring ~ (1|genotype) + (1|trial) + betweenness.rescaled*pc.rescaled,     data = fitnessdata_day1_pc, family = "nbinom")
offspring_day1_pc_clusteringcoeffXpc  <- glmmadmb(offspring ~ (1|genotype) + (1|trial) + clusteringcoeff.rescaled*pc.rescaled, data = fitnessdata_day1_pc, family = "nbinom")
Anova(offspring_day1_pc_indegreeXpc, type = "3")          # Observed significance of Indegree*PC on Day 1:          Chisq = 4.1304, Df = 1, p = 0.04212
Anova(offspring_day1_pc_outdegreeXpc, type = "3")         # Observed significance of Outdegree*PC on Day 1:         Chisq = 2.3412, Df = 1, p = 0.1260
Anova(offspring_day1_pc_eigencentralityXpc, type = "3")   # Observed significance of Eigencentrality*PC on Day 1:   Chisq = 0.7323, Df = 1, p = 0.3921
Anova(offspring_day1_pc_betweennessXpc, type = "3")       # Observed significance of Betweenness*PC on Day 1:       Chisq = 1.3586, Df = 1, p = 0.24377
Anova(offspring_day1_pc_clusteringcoeffXpc, type = "3")   # Observed significance of Clustering Coeff*PC on Day 1:  Chisq = 0.5017, Df = 1, p = 0.47877
offspring_day1_pc_indegree         <- glmmadmb(offspring ~ (1|genotype) + (1|trial) + indegree.rescaled + pc.rescaled,        data = fitnessdata_day1_pc, family = "nbinom")
offspring_day1_pc_outdegree        <- glmmadmb(offspring ~ (1|genotype) + (1|trial) + outdegree.rescaled + pc.rescaled,       data = fitnessdata_day1_pc, family = "nbinom")
offspring_day1_pc_eigencentrality  <- glmmadmb(offspring ~ (1|genotype) + (1|trial) + eigencentrality.rescaled + pc.rescaled, data = fitnessdata_day1_pc, family = "nbinom")
offspring_day1_pc_betweenness      <- glmmadmb(offspring ~ (1|genotype) + (1|trial) + betweenness.rescaled + pc.rescaled,     data = fitnessdata_day1_pc, family = "nbinom")
offspring_day1_pc_clusteringcoeff  <- glmmadmb(offspring ~ (1|genotype) + (1|trial) + clusteringcoeff.rescaled + pc.rescaled, data = fitnessdata_day1_pc, family = "nbinom")
Anova(offspring_day1_pc_indegree, type = "3")          # Observed significance of Indegree on Day 1:          Chisq = 0.8131, Df = 1, p = 0.3672
Anova(offspring_day1_pc_outdegree, type = "3")         # Observed significance of Outdegree on Day 1:         Chisq = 0.5821, Df = 1, p = 0.4455
Anova(offspring_day1_pc_eigencentrality, type = "3")   # Observed significance of Eigencentrality on Day 1:   Chisq = 1.9316, Df = 1, p = 0.1646
Anova(offspring_day1_pc_betweenness, type = "3")       # Observed significance of Betweenness on Day 1:       Chisq = 1.3546, Df = 1, p = 0.2445
Anova(offspring_day1_pc_clusteringcoeff, type = "3")   # Observed significance of Clustering Coeff on Day 1:  Chisq = 0.4087, Df = 1, p = 0.5227

# LIFESPAN (P:C ONLY) # --------------
lifespan_G_pc         <- glmer(lifespan ~ (1|genotype) + (1|trial) + pc.rescaled,        data = fitnessdata_pc, family = poisson)
lifespan_NA_pc        <- glmer(lifespan ~                (1|trial) + pc.rescaled,        data = fitnessdata_pc, family = poisson)
lifespan_G_pc_sim     <- simulateResiduals(fittedModel = lifespan_G_pc, n = 1000)
lifespan_NA_pc_sim    <- simulateResiduals(fittedModel = lifespan_NA_pc, n = 1000)
testResiduals(lifespan_G_pc_sim)  # Overdispersion
testResiduals(lifespan_NA_pc_sim) # Overdispersion
lifespan_G_pc         <- glmmadmb(lifespan ~ (1|genotype) + (1|trial) + pc.rescaled,        data = fitnessdata_pc, family = "nbinom")
lifespan_NA_pc        <- glmmadmb(lifespan ~                (1|trial) + pc.rescaled,        data = fitnessdata_pc, family = "nbinom")
anova(lifespan_G_pc, lifespan_NA_pc)   # Observed significance of G:    LR = 0.4,      Df = 1,  p = 0.5271
Anova(lifespan_G_pc, type = "3")       # Observed significance of PC:   LR = 7.4431,   Df = 1,  p = 0.006368
lifespan_day1_pc_indegreeXpc         <- glmmadmb(lifespan ~ (1|genotype) + (1|trial) + indegree.rescaled*pc.rescaled,        data = fitnessdata_day1_pc, family = "nbinom")
lifespan_day1_pc_outdegreeXpc        <- glmmadmb(lifespan ~ (1|genotype) + (1|trial) + outdegree.rescaled*pc.rescaled,       data = fitnessdata_day1_pc, family = "nbinom") # Fails to converge. Model negative binomial distribution another way
lifespan_day1_pc_outdegreeXpc        <- glmer.nb(lifespan ~ (1|genotype) + (1|trial) + outdegree.rescaled*pc.rescaled,       data = fitnessdata_day1_pc)
lifespan_day1_pc_eigencentralityXpc  <- glmmadmb(lifespan ~ (1|genotype) + (1|trial) + eigencentrality.rescaled*pc.rescaled, data = fitnessdata_day1_pc, family = "nbinom")
lifespan_day1_pc_betweennessXpc      <- glmmadmb(lifespan ~ (1|genotype) + (1|trial) + betweenness.rescaled*pc.rescaled,     data = fitnessdata_day1_pc, family = "nbinom")
lifespan_day1_pc_clusteringcoeffXpc  <- glmmadmb(lifespan ~ (1|genotype) + (1|trial) + clusteringcoeff.rescaled*pc.rescaled, data = fitnessdata_day1_pc, family = "nbinom")
Anova(lifespan_day1_pc_indegreeXpc, type = "3")          # Observed significance of Indegree*PC on Day 1:          Chisq = 0.8157, Df = 1, p = 0.3664
Anova(lifespan_day1_pc_outdegreeXpc, type = "3")         # Observed significance of Outdegree*PC on Day 1:         Chisq = 1.5397, Df = 1, p = 0.2147
Anova(lifespan_day1_pc_eigencentralityXpc, type = "3")   # Observed significance of Eigencentrality*PC on Day 1:   Chisq = 0.0228, Df = 1, p = 0.8800
Anova(lifespan_day1_pc_betweennessXpc, type = "3")       # Observed significance of Betweenness*PC on Day 1:       Chisq = 3.4709, Df = 1, p = 0.062456
Anova(lifespan_day1_pc_clusteringcoeffXpc, type = "3")   # Observed significance of Clustering Coeff*PC on Day 1:  Chisq = 0.4020, Df = 1, p = 0.5261
lifespan_day1_pc_indegree         <- glmmadmb(lifespan ~ (1|genotype) + (1|trial) + indegree.rescaled + pc.rescaled,        data = fitnessdata_day1_pc, family = "nbinom")
lifespan_day1_pc_outdegree        <- glmmadmb(lifespan ~ (1|genotype) + (1|trial) + outdegree.rescaled + pc.rescaled,       data = fitnessdata_day1_pc, family = "nbinom")
lifespan_day1_pc_eigencentrality  <- glmmadmb(lifespan ~ (1|genotype) + (1|trial) + eigencentrality.rescaled + pc.rescaled, data = fitnessdata_day1_pc, family = "nbinom")
lifespan_day1_pc_betweenness      <- glmmadmb(lifespan ~ (1|genotype) + (1|trial) + betweenness.rescaled + pc.rescaled,     data = fitnessdata_day1_pc, family = "nbinom")
lifespan_day1_pc_clusteringcoeff  <- glmmadmb(lifespan ~ (1|genotype) + (1|trial) + clusteringcoeff.rescaled + pc.rescaled, data = fitnessdata_day1_pc, family = "nbinom")
Anova(lifespan_day1_pc_indegree, type = "3")          # Observed significance of Indegree on Day 1:          Chisq = 0.0121, Df = 1, p = 0.91232
Anova(lifespan_day1_pc_outdegree, type = "3")         # Observed significance of Outdegree on Day 1:         Chisq = 0.0011, Df = 1, p = 0.9730
Anova(lifespan_day1_pc_eigencentrality, type = "3")   # Observed significance of Eigencentrality on Day 1:   Chisq = 1.0001, Df = 1, p = 0.31729
Anova(lifespan_day1_pc_betweenness, type = "3")       # Observed significance of Betweenness on Day 1:       Chisq = 1.4239, Df = 1, p = 0.23276
Anova(lifespan_day1_pc_clusteringcoeff, type = "3")   # Observed significance of Clustering Coeff on Day 1:  Chisq = 0.0626, Df = 1, p = 0.8024

# MATINGS (P:C ONLY) # --------------
matingsdata_pc <- matingsdata[which(matingsdata$food=='1' | matingsdata$food=='2' | matingsdata$food=='3'),]
matingsdata_day1_pc <- matingsdata_day1[which(matingsdata_day1$food=='1' | matingsdata_day1$food=='2' | matingsdata_day1$food=='3'),]
matings_G_pc         <- glmer(matings ~ (1|genotype) + (1|trial) + pc.rescaled,        data = matingsdata_pc, family = poisson)
matings_NA_pc        <- glmer(matings ~                (1|trial) + pc.rescaled,        data = matingsdata_pc, family = poisson)
matings_G_pc_sim     <- simulateResiduals(fittedModel = matings_G_pc, n = 1000)
matings_NA_pc_sim    <- simulateResiduals(fittedModel = matings_NA_pc, n = 1000)
testResiduals(matings_G_pc_sim)  # 0 inflation
testResiduals(matings_NA_pc_sim) # 0 inflation
matings_G_pc         <- glmmadmb(matings ~ (1|genotype) + (1|trial) + pc.rescaled,        data = matingsdata_pc, family = "poisson", zeroInflation = TRUE)
matings_NA_pc        <- glmmadmb(matings ~                (1|trial) + pc.rescaled,        data = matingsdata_pc, family = "poisson", zeroInflation = TRUE)
anova(matings_G_pc, matings_NA_pc)    # Observed significance of G:   LR = 28.416,   Df = 1,  p = 9.785e-08
Anova(matings_G_pc, type = "3")       # Observed significance of PC:  LR = 0.2734,   Df = 1,  p = 0.60109
matings_day1_pc_indegreeXpc         <- glmmadmb(matings ~ (1|genotype) + (1|trial) + indegree.rescaled*pc.rescaled,        data = matingsdata_day1_pc, family = "poisson", zeroInflation = TRUE)
matings_day1_pc_outdegreeXpc        <- glmmadmb(matings ~ (1|genotype) + (1|trial) + outdegree.rescaled*pc.rescaled,       data = matingsdata_day1_pc, family = "poisson", zeroInflation = TRUE)
matings_day1_pc_eigencentralityXpc  <- glmmadmb(matings ~ (1|genotype) + (1|trial) + eigencentrality.rescaled*pc.rescaled, data = matingsdata_day1_pc, family = "poisson", zeroInflation = TRUE)
matings_day1_pc_betweennessXpc      <- glmmadmb(matings ~ (1|genotype) + (1|trial) + betweenness.rescaled*pc.rescaled,     data = matingsdata_day1_pc, family = "poisson", zeroInflation = TRUE)
matings_day1_pc_clusteringcoeffXpc  <- glmmadmb(matings ~ (1|genotype) + (1|trial) + clusteringcoeff.rescaled*pc.rescaled, data = matingsdata_day1_pc, family = "poisson", zeroInflation = TRUE)
Anova(matings_day1_pc_indegreeXpc, type = "3")          # Observed significance of Indegree*PC on Day 1:          Chisq = 6.4778, Df = 1, p = 0.01092
Anova(matings_day1_pc_outdegreeXpc, type = "3")         # Observed significance of Outdegree*PC on Day 1:         Chisq = 5.7333, Df = 1, p = 0.01665
Anova(matings_day1_pc_eigencentralityXpc, type = "3")   # Observed significance of Eigencentrality*PC on Day 1:   Chisq = 6.6775, Df = 1, p = 0.009764
Anova(matings_day1_pc_betweennessXpc, type = "3")       # Observed significance of Betweenness*PC on Day 1:       Chisq = 2.1166, Df = 1, p = 0.14571
Anova(matings_day1_pc_clusteringcoeffXpc, type = "3")   # Observed significance of Clustering Coeff*PC on Day 1:  Chisq = 0.0190, Df = 1, p = 0.8903
matings_day1_pc_indegree         <- glmmadmb(matings ~ (1|genotype) + (1|trial) + indegree.rescaled + pc.rescaled,        data = matingsdata_day1_pc, family = "poisson", zeroInflation = TRUE)
matings_day1_pc_outdegree        <- glmmadmb(matings ~ (1|genotype) + (1|trial) + outdegree.rescaled + pc.rescaled,       data = matingsdata_day1_pc, family = "poisson", zeroInflation = TRUE)
matings_day1_pc_eigencentrality  <- glmmadmb(matings ~ (1|genotype) + (1|trial) + eigencentrality.rescaled + pc.rescaled, data = matingsdata_day1_pc, family = "poisson", zeroInflation = TRUE)
matings_day1_pc_betweenness      <- glmmadmb(matings ~ (1|genotype) + (1|trial) + betweenness.rescaled + pc.rescaled,     data = matingsdata_day1_pc, family = "poisson", zeroInflation = TRUE)
matings_day1_pc_clusteringcoeff  <- glmmadmb(matings ~ (1|genotype) + (1|trial) + clusteringcoeff.rescaled + pc.rescaled, data = matingsdata_day1_pc, family = "poisson", zeroInflation = TRUE)
Anova(matings_day1_pc_indegreeXpc, type = "3")            # Observed significance of Indegree on Day 1:          Chisq = 2.8080, Df = 1, p = 0.09380
Anova(matings_day1_pc_outdegreeXpc, type = "3")           # Observed significance of Outdegree on Day 1:         Chisq = 2.2196, Df = 1, p = 0.13627
Anova(matings_day1_pc_eigencentralityXpc, type = "3")     # Observed significance of Eigencentrality on Day 1:   Chisq = 1.3243, Df = 1, p = 0.249824
Anova(matings_day1_pc_betweenness, type = "3")            # Observed significance of Betweenness on Day 1:       Chisq = 0.0794, Df = 1, p = 0.77811
Anova(matings_day1_pc_clusteringcoeff, type = "3")        # Observed significance of Clustering Coeff on Day 1:  Chisq = 0.0630, Df = 1, p = 0.8019

# LATENCY (P:C ONLY) # --------------
latency_G_pc  <- coxme(Surv(latency.censored, is.censored) ~ (1|genotype) + (1|trial) + pc.rescaled, data = matingsdata_pc)
latency_NA_pc <- coxme(Surv(latency.censored, is.censored) ~                (1|trial) + pc.rescaled, data = matingsdata_pc)
anova(latency_G_pc, latency_NA_pc)    # Observed significance of G:   LR = 24.579,  Df = 1,  p = 7.133e-07
Anova(latency_G_pc, type = "3")       # Observed significance of E:   LR = 0.661,   Df = 1,  p = 0.4162
latency_day1_pc_indegreeXpc        <- coxme(Surv(latency.censored, is.censored) ~ pc.rescaled*indegree.rescaled        + (1|genotype) + (1|trial), data=matingsdata_day1_pc)
latency_day1_pc_outdegreeXpc       <- coxme(Surv(latency.censored, is.censored) ~ pc.rescaled*outdegree.rescaled       + (1|genotype) + (1|trial), data=matingsdata_day1_pc)
latency_day1_pc_eigencentralityXpc <- coxme(Surv(latency.censored, is.censored) ~ pc.rescaled*eigencentrality.rescaled + (1|genotype) + (1|trial), data=matingsdata_day1_pc)
latency_day1_pc_betweennessXpc     <- coxme(Surv(latency.censored, is.censored) ~ pc.rescaled*betweenness.rescaled     + (1|genotype) + (1|trial), data=matingsdata_day1_pc)
latency_day1_pc_clusteringcoeffXpc <- coxme(Surv(latency.censored, is.censored) ~ pc.rescaled*clusteringcoeff.rescaled + (1|genotype) + (1|trial), data=matingsdata_day1_pc)
Anova(latency_day1_pc_indegreeXpc, type = "3")          # Observed significance of Indegree*PC on Day 1:          Chisq = 7.3174,   Df = 1, p = 0.006829
Anova(latency_day1_pc_outdegreeXpc, type = "3")         # Observed significance of Outdegree*PC on Day 1:         Chisq = 8.1828,   Df = 1, p = 0.004229
Anova(latency_day1_pc_eigencentralityXpc, type = "3")   # Observed significance of Eigencentrality*PC on Day 1:   Chisq = 4.3298,   Df = 1, p = 0.03745
Anova(latency_day1_pc_betweennessXpc, type = "3")       # Observed significance of Betweenness*PC on Day 1:       Chisq = 0.3946,   Df = 1, p = 0.5299
Anova(latency_day1_pc_clusteringcoeffXpc, type = "3")   # Observed significance of Clustering Coeff*PC on Day 1:  Chisq = 0.0572,   Df = 1, p = 0.8110
latency_day1_pc_indegree        <- coxme(Surv(latency.censored, is.censored) ~ pc.rescaled + indegree.rescaled        + (1|genotype) + (1|trial), data=matingsdata_day1_pc)
latency_day1_pc_outdegree       <- coxme(Surv(latency.censored, is.censored) ~ pc.rescaled + outdegree.rescaled       + (1|genotype) + (1|trial), data=matingsdata_day1_pc)
latency_day1_pc_eigencentrality <- coxme(Surv(latency.censored, is.censored) ~ pc.rescaled + eigencentrality.rescaled + (1|genotype) + (1|trial), data=matingsdata_day1_pc)
latency_day1_pc_betweenness     <- coxme(Surv(latency.censored, is.censored) ~ pc.rescaled + betweenness.rescaled     + (1|genotype) + (1|trial), data=matingsdata_day1_pc)
latency_day1_pc_clusteringcoeff <- coxme(Surv(latency.censored, is.censored) ~ pc.rescaled + clusteringcoeff.rescaled + (1|genotype) + (1|trial), data=matingsdata_day1_pc)
Anova(latency_day1_pc_indegreeXpc, type = "3")       # Observed significance of Indegree on Day 1:          Chisq = 5.0059, Df = 1, p = 0.025261
Anova(latency_day1_pc_outdegreeXpc, type = "3")      # Observed significance of Outdegree on Day 1:         Chisq = 4.7971, Df = 1, p = 0.028508
Anova(latency_day1_pc_eigencentrality, type = "3")   # Observed significance of Eigencentrality on Day 1:   Chisq = 2.0039, Df = 1, p = 0.1569
Anova(latency_day1_pc_betweenness, type = "3")       # Observed significance of Betweenness on Day 1:       Chisq = 0.0511, Df = 1, p = 0.8212
Anova(latency_day1_pc_clusteringcoeff, type = "3")   # Observed significance of Clustering Coeff on Day 1:  Chisq = 0.1904, Df = 1, p = 0.6626

# OFFSPRING (DILUTION ONLY) # --------------
fitnessdata_dilution <- fitnessdata[which(fitnessdata$food=='1' | fitnessdata$food=='4' | fitnessdata$food=='5'),]
fitnessdata_day1_dilution <- fitnessdata_day1[which(fitnessdata_day1$food=='1' | fitnessdata_day1$food=='4' | fitnessdata_day1$food=='5'),]
offspring_G_dilution         <- glmer(offspring ~ (1|genotype) + (1|trial) + dilution.rescaled,        data = fitnessdata_dilution, family = poisson)
offspring_NA_dilution        <- glmer(offspring ~                (1|trial) + dilution.rescaled,        data = fitnessdata_dilution, family = poisson)
offspring_G_dilution_sim     <- simulateResiduals(fittedModel = offspring_G_dilution, n = 1000)
offspring_NA_dilution_sim    <- simulateResiduals(fittedModel = offspring_NA_dilution, n = 1000)
testResiduals(offspring_G_dilution_sim)  # Overdispersion
testResiduals(offspring_NA_dilution_sim) # Overdispersion
offspring_G_dilution         <- glmmadmb(offspring ~ (1|genotype) + (1|trial) + dilution.rescaled,        data = fitnessdata_dilution, family = "nbinom")
offspring_NA_dilution        <- glmmadmb(offspring ~                (1|trial) + dilution.rescaled,        data = fitnessdata_dilution, family = "nbinom")
anova(offspring_G_dilution, offspring_NA_dilution)  # Observed significance of G:          LR = 0.000,     Df = 1,  p = 1
Anova(offspring_G_dilution, type = "3")             # Observed significance of Dilution:   LR = 0.0002,    Df = 1,  p = 0.9887
offspring_day1_dilution_indegreeXdilution         <- glmmadmb(offspring ~ (1|genotype) + (1|trial) + indegree.rescaled*dilution.rescaled,        data = fitnessdata_day1_dilution, family = "nbinom")
offspring_day1_dilution_outdegreeXdilution        <- glmmadmb(offspring ~ (1|genotype) + (1|trial) + outdegree.rescaled*dilution.rescaled,       data = fitnessdata_day1_dilution, family = "nbinom")
offspring_day1_dilution_eigencentralityXdilution  <- glmmadmb(offspring ~ (1|genotype) + (1|trial) + eigencentrality.rescaled*dilution.rescaled, data = fitnessdata_day1_dilution, family = "nbinom")
offspring_day1_dilution_betweennessXdilution      <- glmmadmb(offspring ~ (1|genotype) + (1|trial) + betweenness.rescaled*dilution.rescaled,     data = fitnessdata_day1_dilution, family = "nbinom")
offspring_day1_dilution_clusteringcoeffXdilution  <- glmmadmb(offspring ~ (1|genotype) + (1|trial) + clusteringcoeff.rescaled*dilution.rescaled, data = fitnessdata_day1_dilution, family = "nbinom")
Anova(offspring_day1_dilution_indegreeXdilution, type = "3")          # Observed significance of Indegree*Dilution on Day 1:          Chisq = 2.0018, Df = 1, p = 0.1571
Anova(offspring_day1_dilution_outdegreeXdilution, type = "3")         # Observed significance of Outdegree*Dilution on Day 1:         Chisq = 1.3107, Df = 1, p = 0.2523
Anova(offspring_day1_dilution_eigencentralityXdilution, type = "3")   # Observed significance of Eigencentrality*Dilution on Day 1:   Chisq = 1.0516, Df = 1, p = 0.3051
Anova(offspring_day1_dilution_betweennessXdilution, type = "3")       # Observed significance of Betweenness*Dilution on Day 1:       Chisq = 0.7988, Df = 1, p = 0.3714
Anova(offspring_day1_dilution_clusteringcoeffXdilution, type = "3")   # Observed significance of Clustering Coeff*Dilution on Day 1:  Chisq = 0.0001, Df = 1, p = 0.993107
offspring_day1_dilution_indegree         <- glmmadmb(offspring ~ (1|genotype) + (1|trial) + indegree.rescaled + dilution.rescaled,        data = fitnessdata_day1_dilution, family = "nbinom")
offspring_day1_dilution_outdegree        <- glmmadmb(offspring ~ (1|genotype) + (1|trial) + outdegree.rescaled + dilution.rescaled,       data = fitnessdata_day1_dilution, family = "nbinom")
offspring_day1_dilution_eigencentrality  <- glmmadmb(offspring ~ (1|genotype) + (1|trial) + eigencentrality.rescaled + dilution.rescaled, data = fitnessdata_day1_dilution, family = "nbinom")
offspring_day1_dilution_betweenness      <- glmmadmb(offspring ~ (1|genotype) + (1|trial) + betweenness.rescaled + dilution.rescaled,     data = fitnessdata_day1_dilution, family = "nbinom")
offspring_day1_dilution_clusteringcoeff  <- glmmadmb(offspring ~ (1|genotype) + (1|trial) + clusteringcoeff.rescaled + dilution.rescaled, data = fitnessdata_day1_dilution, family = "nbinom")
Anova(offspring_day1_dilution_indegree, type = "3")          # Observed significance of Indegree on Day 1:          Chisq = 1.8522, Df = 1, p = 0.1735
Anova(offspring_day1_dilution_outdegree, type = "3")         # Observed significance of Outdegree on Day 1:         Chisq = 0.9327, Df = 1, p = 0.3342
Anova(offspring_day1_dilution_eigencentrality, type = "3")   # Observed significance of Eigencentrality on Day 1:   Chisq = 0.9678, Df = 1, p = 0.3252
Anova(offspring_day1_dilution_betweenness, type = "3")       # Observed significance of Betweenness on Day 1:       Chisq = 1.6170, Df = 1, p = 0.2035
Anova(offspring_day1_dilution_clusteringcoeff, type = "3")   # Observed significance of Clustering Coeff on Day 1:  Chisq = 0.9990, Df = 1, p = 0.3176

# LIFESPAN (DILUTION ONLY) # --------------
lifespan_G_dilution         <- glmer(lifespan ~ (1|genotype) + (1|trial) + dilution.rescaled,        data = fitnessdata_dilution, family = poisson)
lifespan_NA_dilution        <- glmer(lifespan ~                (1|trial) + dilution.rescaled,        data = fitnessdata_dilution, family = poisson)
lifespan_G_dilution_sim     <- simulateResiduals(fittedModel = lifespan_G_dilution, n = 1000)
lifespan_NA_dilution_sim    <- simulateResiduals(fittedModel = lifespan_NA_dilution, n = 1000)
testResiduals(lifespan_G_dilution_sim)  # Overdispersion
testResiduals(lifespan_NA_dilution_sim) # Overdispersion
lifespan_G_dilution         <- glmmadmb(lifespan ~ (1|genotype) + (1|trial) + dilution.rescaled,        data = fitnessdata_dilution, family = "nbinom")
lifespan_NA_dilution        <- glmmadmb(lifespan ~                (1|trial) + dilution.rescaled,        data = fitnessdata_dilution, family = "nbinom")
anova(lifespan_G_dilution, lifespan_NA_dilution)   # Observed significance of G:          LR = 0.00,     Df = 1,  p = 1
Anova(lifespan_G_dilution, type = "3")             # Observed significance of Dilution:   LR = 10.815,   Df = 1,  p = 0.001007
lifespan_day1_dilution_indegreeXdilution         <- glmmadmb(lifespan ~ (1|genotype) + (1|trial) + indegree.rescaled*dilution.rescaled,        data = fitnessdata_day1_dilution, family = "nbinom")
lifespan_day1_dilution_outdegreeXdilution        <- glmmadmb(lifespan ~ (1|genotype) + (1|trial) + outdegree.rescaled*dilution.rescaled,       data = fitnessdata_day1_dilution, family = "nbinom")
lifespan_day1_dilution_eigencentralityXdilution  <- glmmadmb(lifespan ~ (1|genotype) + (1|trial) + eigencentrality.rescaled*dilution.rescaled, data = fitnessdata_day1_dilution, family = "nbinom")
lifespan_day1_dilution_betweennessXdilution      <- glmmadmb(lifespan ~ (1|genotype) + (1|trial) + betweenness.rescaled*dilution.rescaled,     data = fitnessdata_day1_dilution, family = "nbinom")
lifespan_day1_dilution_clusteringcoeffXdilution  <- glmmadmb(lifespan ~ (1|genotype) + (1|trial) + clusteringcoeff.rescaled*dilution.rescaled, data = fitnessdata_day1_dilution, family = "nbinom")
Anova(lifespan_day1_dilution_indegreeXdilution, type = "3")          # Observed significance of Indegree*Dilution on Day 1:          Chisq = 0.3453, Df = 1, p = 0.5568
Anova(lifespan_day1_dilution_outdegreeXdilution, type = "3")         # Observed significance of Outdegree*Dilution on Day 1:         Chisq = 0.8531, Df = 1, p = 0.3557
Anova(lifespan_day1_dilution_eigencentralityXdilution, type = "3")   # Observed significance of Eigencentrality*Dilution on Day 1:   Chisq = 0.1883, Df = 1, p = 0.6644
Anova(lifespan_day1_dilution_betweennessXdilution, type = "3")       # Observed significance of Betweenness*Dilution on Day 1:       Chisq = 0.9195, Df = 1, p = 0.337618
Anova(lifespan_day1_dilution_clusteringcoeffXdilution, type = "3")   # Observed significance of Clustering Coeff*Dilution on Day 1:  Chisq = 0.8134, Df = 1, p = 0.3671
lifespan_day1_dilution_indegree         <- glmmadmb(lifespan ~ (1|genotype) + (1|trial) + indegree.rescaled + dilution.rescaled,        data = fitnessdata_day1_dilution, family = "nbinom")
lifespan_day1_dilution_outdegree        <- glmmadmb(lifespan ~ (1|genotype) + (1|trial) + outdegree.rescaled + dilution.rescaled,       data = fitnessdata_day1_dilution, family = "nbinom")
lifespan_day1_dilution_eigencentrality  <- glmmadmb(lifespan ~ (1|genotype) + (1|trial) + eigencentrality.rescaled + dilution.rescaled, data = fitnessdata_day1_dilution, family = "nbinom")
lifespan_day1_dilution_betweenness      <- glmmadmb(lifespan ~ (1|genotype) + (1|trial) + betweenness.rescaled + dilution.rescaled,     data = fitnessdata_day1_dilution, family = "nbinom")
lifespan_day1_dilution_clusteringcoeff  <- glmmadmb(lifespan ~ (1|genotype) + (1|trial) + clusteringcoeff.rescaled + dilution.rescaled, data = fitnessdata_day1_dilution, family = "nbinom")
Anova(lifespan_day1_dilution_indegree, type = "3")          # Observed significance of Indegree on Day 1:          Chisq = 0.2160, Df = 1, p = 0.642105
Anova(lifespan_day1_dilution_outdegree, type = "3")         # Observed significance of Outdegree on Day 1:         Chisq = 0.4471, Df = 1, p = 0.503719
Anova(lifespan_day1_dilution_eigencentrality, type = "3")   # Observed significance of Eigencentrality on Day 1:   Chisq = 0.3469, Df = 1, p = 0.555884
Anova(lifespan_day1_dilution_betweenness, type = "3")       # Observed significance of Betweenness on Day 1:       Chisq = 3.1883, Df = 1, p = 0.074166
Anova(lifespan_day1_dilution_clusteringcoeff, type = "3")   # Observed significance of Clustering Coeff on Day 1:  Chisq = 0.0497, Df = 1, p = 0.82361

# MATINGS (DILUTION ONLY) # --------------
matingsdata_dilution <- matingsdata[which(matingsdata$food=='1' | matingsdata$food=='4' | matingsdata$food=='5'),]
matingsdata_day1_dilution <- matingsdata_day1[which(matingsdata_day1$food=='1' | matingsdata_day1$food=='4' | matingsdata_day1$food=='5'),]
matings_G_dilution         <- glmer(matings ~ (1|genotype) + (1|trial) + dilution.rescaled,        data = matingsdata_dilution, family = poisson)
matings_NA_dilution        <- glmer(matings ~                (1|trial) + dilution.rescaled,        data = matingsdata_dilution, family = poisson)
matings_G_dilution_sim     <- simulateResiduals(fittedModel = matings_G_dilution, n = 1000)
matings_NA_dilution_sim    <- simulateResiduals(fittedModel = matings_NA_dilution, n = 1000)
testResiduals(matings_G_dilution_sim)  # 0 inflation
testResiduals(matings_NA_dilution_sim) # 0 inflation
matings_G_dilution         <- glmmadmb(matings ~ (1|genotype) + (1|trial) + dilution.rescaled,        data = matingsdata_dilution, family = "poisson", zeroInflation = TRUE)
matings_NA_dilution        <- glmmadmb(matings ~                (1|trial) + dilution.rescaled,        data = matingsdata_dilution, family = "poisson", zeroInflation = TRUE)
anova(matings_G_dilution, matings_NA_dilution)    # Observed significance of G:         LR = 52.814,   Df = 1,  p = 3.667e-13
Anova(matings_G_dilution, type = "3")             # Observed significance of Dilution:  LR = 0.1126,   Df = 1,  p = 0.73716
matings_day1_dilution_indegreeXdilution         <- glmmadmb(matings ~ (1|genotype) + (1|trial) + indegree.rescaled*dilution.rescaled,        data = matingsdata_day1_dilution, family = "poisson", zeroInflation = TRUE)
matings_day1_dilution_outdegreeXdilution        <- glmmadmb(matings ~ (1|genotype) + (1|trial) + outdegree.rescaled*dilution.rescaled,       data = matingsdata_day1_dilution, family = "poisson", zeroInflation = TRUE)
matings_day1_dilution_eigencentralityXdilution  <- glmmadmb(matings ~ (1|genotype) + (1|trial) + eigencentrality.rescaled*dilution.rescaled, data = matingsdata_day1_dilution, family = "poisson", zeroInflation = TRUE)
matings_day1_dilution_betweennessXdilution      <- glmmadmb(matings ~ (1|genotype) + (1|trial) + betweenness.rescaled*dilution.rescaled,     data = matingsdata_day1_dilution, family = "poisson", zeroInflation = TRUE)
matings_day1_dilution_clusteringcoeffXdilution  <- glmmadmb(matings ~ (1|genotype) + (1|trial) + clusteringcoeff.rescaled*dilution.rescaled, data = matingsdata_day1_dilution, family = "poisson", zeroInflation = TRUE)
Anova(matings_day1_dilution_indegreeXdilution, type = "3")          # Observed significance of Indegree*Dilution on Day 1:          Chisq = 7.0982, Df = 1, p = 0.007716
Anova(matings_day1_dilution_outdegreeXdilution, type = "3")         # Observed significance of Outdegree*Dilution on Day 1:         Chisq = 2.6985, Df = 1, p = 0.1004
Anova(matings_day1_dilution_eigencentralityXdilution, type = "3")   # Observed significance of Eigencentrality*Dilution on Day 1:   Chisq = 9.8587, Df = 1, p = 0.001690
Anova(matings_day1_dilution_betweennessXdilution, type = "3")       # Observed significance of Betweenness*Dilution on Day 1:       Chisq = 0.5159, Df = 1, p = 0.4726
Anova(matings_day1_dilution_clusteringcoeffXdilution, type = "3")   # Observed significance of Clustering Coeff*Dilution on Day 1:  Chisq = 1.6757, Df = 1, p = 0.19550
matings_day1_dilution_indegree         <- glmmadmb(matings ~ (1|genotype) + (1|trial) + indegree.rescaled + dilution.rescaled,        data = matingsdata_day1_dilution, family = "poisson", zeroInflation = TRUE)
matings_day1_dilution_outdegree        <- glmmadmb(matings ~ (1|genotype) + (1|trial) + outdegree.rescaled + dilution.rescaled,       data = matingsdata_day1_dilution, family = "poisson", zeroInflation = TRUE)
matings_day1_dilution_eigencentrality  <- glmmadmb(matings ~ (1|genotype) + (1|trial) + eigencentrality.rescaled + dilution.rescaled, data = matingsdata_day1_dilution, family = "poisson", zeroInflation = TRUE)
matings_day1_dilution_betweenness      <- glmmadmb(matings ~ (1|genotype) + (1|trial) + betweenness.rescaled + dilution.rescaled,     data = matingsdata_day1_dilution, family = "poisson", zeroInflation = TRUE)
matings_day1_dilution_clusteringcoeff  <- glmmadmb(matings ~ (1|genotype) + (1|trial) + clusteringcoeff.rescaled + dilution.rescaled, data = matingsdata_day1_dilution, family = "poisson", zeroInflation = TRUE)
Anova(matings_day1_dilution_indegreeXdilution, type = "3")            # Observed significance of Indegree on Day 1:          Chisq = 6.0898, Df = 1, p = 0.013597
Anova(matings_day1_dilution_outdegree, type = "3")                    # Observed significance of Outdegree on Day 1:         Chisq = 0.0032, Df = 1, p = 0.9547
Anova(matings_day1_dilution_eigencentralityXdilution, type = "3")     # Observed significance of Eigencentrality on Day 1:   Chisq = 7.0589, Df = 1, p = 0.007887
Anova(matings_day1_dilution_betweenness, type = "3")                  # Observed significance of Betweenness on Day 1:       Chisq = 0.8184, Df = 1, p = 0.3656
Anova(matings_day1_dilution_clusteringcoeff, type = "3")              # Observed significance of Clustering Coeff on Day 1:  Chisq = 1.0516, Df = 1, p = 0.3051

# LATENCY (DILUTION ONLY) # --------------
latency_G_dilution  <- coxme(Surv(latency.censored, is.censored) ~ (1|genotype) + (1|trial) + dilution.rescaled, data = matingsdata_dilution)
latency_NA_dilution <- coxme(Surv(latency.censored, is.censored) ~                (1|trial) + dilution.rescaled, data = matingsdata_dilution)
anova(latency_G_dilution, latency_NA_dilution)    # Observed significance of G:   LR = 48.182,  Df = 1,  p = 3.885e-12
Anova(latency_G_dilution, type = "3")             # Observed significance of E:   LR = 0.8908,  Df = 1,  p = 0.3453
latency_day1_dilution_indegreeXdilution        <- coxme(Surv(latency.censored, is.censored) ~ dilution.rescaled*indegree.rescaled        + (1|genotype) + (1|trial), data=matingsdata_day1_dilution)
latency_day1_dilution_outdegreeXdilution       <- coxme(Surv(latency.censored, is.censored) ~ dilution.rescaled*outdegree.rescaled       + (1|genotype) + (1|trial), data=matingsdata_day1_dilution)
latency_day1_dilution_eigencentralityXdilution <- coxme(Surv(latency.censored, is.censored) ~ dilution.rescaled*eigencentrality.rescaled + (1|genotype) + (1|trial), data=matingsdata_day1_dilution)
latency_day1_dilution_betweennessXdilution     <- coxme(Surv(latency.censored, is.censored) ~ dilution.rescaled*betweenness.rescaled     + (1|genotype) + (1|trial), data=matingsdata_day1_dilution)
latency_day1_dilution_clusteringcoeffXdilution <- coxme(Surv(latency.censored, is.censored) ~ dilution.rescaled*clusteringcoeff.rescaled + (1|genotype) + (1|trial), data=matingsdata_day1_dilution)
Anova(latency_day1_dilution_indegreeXdilution, type = "3")          # Observed significance of Indegree*dilution on Day 1:          Chisq = 10.3672,   Df = 1, p = 0.001283
Anova(latency_day1_dilution_outdegreeXdilution, type = "3")         # Observed significance of Outdegree*dilution on Day 1:         Chisq = 6.4710,    Df = 1, p = 0.01096
Anova(latency_day1_dilution_eigencentralityXdilution, type = "3")   # Observed significance of Eigencentrality*dilution on Day 1:   Chisq = 10.6279,   Df = 1, p = 0.001114
Anova(latency_day1_dilution_betweennessXdilution, type = "3")       # Observed significance of Betweenness*dilution on Day 1:       Chisq = 2.9097,    Df = 1, p = 0.08805
Anova(latency_day1_dilution_clusteringcoeffXdilution, type = "3")   # Observed significance of Clustering Coeff*dilution on Day 1:  Chisq = 2.2930,    Df = 1, p = 0.12996
latency_day1_dilution_indegree        <- coxme(Surv(latency.censored, is.censored) ~ dilution.rescaled + indegree.rescaled        + (1|genotype) + (1|trial), data=matingsdata_day1_dilution)
latency_day1_dilution_outdegree       <- coxme(Surv(latency.censored, is.censored) ~ dilution.rescaled + outdegree.rescaled       + (1|genotype) + (1|trial), data=matingsdata_day1_dilution)
latency_day1_dilution_eigencentrality <- coxme(Surv(latency.censored, is.censored) ~ dilution.rescaled + eigencentrality.rescaled + (1|genotype) + (1|trial), data=matingsdata_day1_dilution)
latency_day1_dilution_betweenness     <- coxme(Surv(latency.censored, is.censored) ~ dilution.rescaled + betweenness.rescaled     + (1|genotype) + (1|trial), data=matingsdata_day1_dilution)
latency_day1_dilution_clusteringcoeff <- coxme(Surv(latency.censored, is.censored) ~ dilution.rescaled + clusteringcoeff.rescaled + (1|genotype) + (1|trial), data=matingsdata_day1_dilution)
Anova(latency_day1_dilution_indegreeXdilution, type = "3")          # Observed significance of Indegree on Day 1:          Chisq = 10.1084, Df = 1, p = 0.001476
Anova(latency_day1_dilution_outdegree, type = "3")                  # Observed significance of Outdegree on Day 1:         Chisq = 0.3987,  Df = 1, p = 0.5278
Anova(latency_day1_dilution_eigencentralityXdilution, type = "3")   # Observed significance of Eigencentrality on Day 1:   Chisq = 7.2122,  Df = 1, p = 0.007241
Anova(latency_day1_dilution_betweenness, type = "3")                # Observed significance of Betweenness on Day 1:       Chisq = 0.4724,  Df = 1, p = 0.4919
Anova(latency_day1_dilution_clusteringcoeff, type = "3")            # Observed significance of Clustering Coeff on Day 1:  Chisq = 2.9407,  Df = 1, p = 0.08638

# FIGURE 1 # --------------------
# Indegree
network$indegree.seconds <- network$indegree / 30     # Convert indegree into seconds of interactions
ggplot(data = network, aes(x = fct_reorder(genotype, indegree.seconds, .fun = median, .desc = TRUE), y = indegree.seconds, fill = sex)) +
  geom_boxplot(outlier.shape = NA, lwd = 0.75) +
  scale_fill_manual(values = c("lightgoldenrod2", "darkslategrey")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.75), axis.ticks = element_line(size = 0.75), legend.position='none', 
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1), axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  coord_cartesian(ylim = c(0, 1050)) +
  xlab("Genotypes") +
  ylab("Instrength (s)")
# Outdegree
network$outdegree.seconds <- network$outdegree / 30     # Convert outdegree into seconds of interactions
ggplot(data = network, aes(x = fct_reorder(genotype, outdegree.seconds, .fun = median, .desc = TRUE), y = outdegree.seconds, fill = sex)) +
  geom_boxplot(outlier.shape = NA, lwd = 0.75) +
  scale_fill_manual(values = c("lightgoldenrod2", "darkslategrey")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.75), axis.ticks = element_line(size = 0.75), legend.position='none', 
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1), axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  coord_cartesian(ylim = c(0, 1050)) +
  xlab("Genotypes") +
  ylab("Outstrength (s)")
# Clustering Coefficient
ggplot(data = network, aes(x = fct_reorder(genotype, clusteringcoeff, .fun = median, .desc = TRUE), y = clusteringcoeff, fill = sex)) +
  geom_boxplot(outlier.shape = NA, lwd = 0.75) +
  scale_fill_manual(values = c("lightgoldenrod2", "darkslategrey")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.75), axis.ticks = element_line(size = 0.75), legend.position='none', 
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1), axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  coord_cartesian(ylim = c(0.85, 1.05)) +
  xlab("Genotypes") +
  ylab("Clustering Coefficient")
# Betweenness
ggplot(data = network, aes(x = fct_reorder(genotype, betweenness, .fun = median, .desc = TRUE), y = betweenness, fill = sex)) +
  geom_boxplot(outlier.shape = NA, lwd = 0.75) +
  scale_fill_manual(values = c("lightgoldenrod2", "darkslategrey")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.75), axis.ticks = element_line(size = 0.75), legend.position='none', 
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1), axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  coord_cartesian(ylim = c(0, 150)) +
  xlab("Genotypes") +
  ylab("Betweenness Centrality")
# Eigencentrality
ggplot(data = network, aes(x = fct_reorder(genotype, eigencentrality, .fun = median, .desc = TRUE), y = eigencentrality, fill = sex)) +
  geom_boxplot(outlier.shape = NA, lwd = 0.75) +
  scale_fill_manual(values = c("lightgoldenrod2", "darkslategrey")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.75), axis.ticks = element_line(size = 0.75), legend.position='none', 
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1), axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  coord_cartesian(ylim = c(0.05, 0.4)) +
  xlab("Genotypes") +
  ylab("Eigenvector Centrality")

# FIGURE 2 # ------------------
# Sex*P:C effect on Outdegree
network_pc$outdegree.seconds <- network_pc$outdegree / 30     # Convert outdegree into seconds of interactions
network_pc$pc.labels <- ifelse(network_pc$pc == "1", "1:1", ifelse(network_pc$pc == "2", "1:2", "1:4"))
ggplot(data = network_pc, aes(x = pc.labels, y = outdegree.seconds, fill = sex)) +
  geom_boxplot(outlier.shape = NA, lwd = 0.75) + 
  scale_fill_manual(values = c("lightgoldenrod2", "darkslategrey")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.75), axis.ticks = element_line(size = 0.75), legend.position='none',
        axis.text.x = element_text(size = 20), axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  coord_cartesian(ylim = c(0, 1200)) +
  xlab("Protein:Carbohydrate") +
  ylab("Outstrength (s)")
# Sex*Caloric Concentration effect on Outdegree
network_dilution$outdegree.seconds <- network_dilution$outdegree / 30     # Convert outdegree into seconds of interactions
network_dilution$dilution.labels <- ifelse(network_dilution$dilution == "1", "1x", ifelse(network_dilution$dilution == "2", "2x", "4x"))
ggplot(data = network_dilution, aes(x = dilution.labels, y = outdegree.seconds, fill = sex)) +
  geom_boxplot(outlier.shape = NA, lwd = 0.75) + 
  scale_fill_manual(values = c("lightgoldenrod2", "darkslategrey")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.75), axis.ticks = element_line(size = 0.75), legend.position='none',
        axis.text.x = element_text(size = 20), axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  coord_cartesian(ylim = c(0, 1200)) +
  xlab("Caloric Concentration") +
  ylab("Outstrength (s)")
# Sex*P:C effect on Clustering Coefficient
ggplot(data = network_pc, aes(x = pc.labels, y = clusteringcoeff, fill = sex)) +
  geom_boxplot(outlier.shape = NA, lwd = 0.75) + 
  scale_fill_manual(values = c("lightgoldenrod2", "darkslategrey")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.75), axis.ticks = element_line(size = 0.75), legend.position='none',
        axis.text.x = element_text(size = 20), axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  coord_cartesian(ylim = c(0.85, 1.05)) +
  xlab("Protein:Carbohydrate") +
  ylab("Clustering Coefficient")
# Sex*Caloric Concentration effect on Clustering Coefficient
ggplot(data = network_dilution, aes(x = dilution.labels, y = clusteringcoeff, fill = sex)) +
  geom_boxplot(outlier.shape = NA, lwd = 0.75) + 
  scale_fill_manual(values = c("lightgoldenrod2", "darkslategrey")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.75), axis.ticks = element_line(size = 0.75), legend.position='none',
        axis.text.x = element_text(size = 20), axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  coord_cartesian(ylim = c(0.85, 1.05)) +
  xlab("Caloric Concentration") +
  ylab("Clustering Coefficient")

# FIGURE 3 # ---------------
# Specify color schemes for P:C and Dilution
matingsdata_day1_pc$pc.color <- ifelse(matingsdata_day1_pc$pc == "1", "salmon", ifelse(matingsdata_day1_pc$pc == "2", "firebrick4", "darkslategrey"))
matingsdata_day1_dilution$dilution.color <- ifelse(matingsdata_day1_dilution$dilution == "1", "darkslategrey1", ifelse(matingsdata_day1_dilution$dilution == "2", "darkslategrey4", "darkslategrey"))
# Change indegree and outdegree to seconds
matingsdata_day1_pc$indegree.seconds <- matingsdata_day1_pc$indegree / 30
matingsdata_day1_pc$outdegree.seconds <- matingsdata_day1_pc$outdegree / 30
matingsdata_day1_dilution$indegree.seconds <- matingsdata_day1_dilution$indegree / 30
matingsdata_day1_dilution$outdegree.seconds <- matingsdata_day1_dilution$outdegree / 30
# Indegree (Day 1) and P:C effects on matings
ggplot(data = matingsdata_day1_pc, aes(x = indegree.seconds, y = matings, color = pc.color)) + 
  geom_point(size = 3, shape = 16, alpha = 0.6) +
  geom_smooth(method="lm", se=TRUE, lwd = 0) +
  geom_smooth(method="lm", se=FALSE, lwd = 2) +
  scale_color_manual(values = c("darkslategrey", "salmon", "firebrick4")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.75), axis.ticks = element_line(size = 0.75), legend.position='none', 
        axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),
        axis.text.y = element_text(size = 25), axis.title.y = element_text(size = 25)) +
  scale_y_continuous(breaks=c(0,1,2,3,4)) +
  scale_x_continuous(breaks=c(400,800,1200,1600)) +
  coord_cartesian(ylim = c(0.0, 4.0)) +
  coord_cartesian(xlim = c(200, 1650)) +
  xlab("Instrength (s)") +
  ylab("Number of Matings")
# Outdegree (Day 1) and P:C effects on matings
ggplot(data = matingsdata_day1_pc, aes(x = outdegree.seconds, y = matings, color = pc.color)) + 
  geom_point(size = 3, shape = 16, alpha = 0.6) +
  geom_smooth(method="lm", se=TRUE, lwd = 0) +
  geom_smooth(method="lm", se=FALSE, lwd = 2) +
  scale_color_manual(values = c("darkslategrey", "salmon", "firebrick4")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.75), axis.ticks = element_line(size = 0.75), legend.position='none', 
        axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),
        axis.text.y = element_text(size = 25), axis.title.y = element_text(size = 25)) +
  scale_y_continuous(breaks=c(0,1,2,3,4)) +
  scale_x_continuous(breaks=c(400,800,1200,1600)) +
  coord_cartesian(ylim = c(0.0, 4.0)) +
  coord_cartesian(xlim = c(200, 1650)) +
  xlab("Outstrength (s)") +
  ylab("Number of Matings")
# Eigencentrality (Day 1) and P:C effects on matings
ggplot(data = matingsdata_day1_pc, aes(x = eigencentrality, y = matings, color = pc.color)) + 
  geom_point(size = 3, shape = 16, alpha = 0.6) +
  geom_smooth(method="lm", se=TRUE, lwd = 0) +
  geom_smooth(method="lm", se=FALSE, lwd = 2) +
  scale_color_manual(values = c("darkslategrey", "salmon", "firebrick4")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.75), axis.ticks = element_line(size = 0.75), legend.position='none', 
        axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),
        axis.text.y = element_text(size = 25), axis.title.y = element_text(size = 25)) +
  scale_y_continuous(breaks=c(0,1,2,3,4)) +
  coord_cartesian(ylim = c(-0.3, 4.0)) +
  xlab("Eigenvector Centrality") +
  ylab("Number of Matings")
# Indegree (Day 1) and Dilution effects on matings
ggplot(data = matingsdata_day1_dilution, aes(x = indegree.seconds, y = matings, color = dilution.color)) + 
  geom_point(size = 3, shape = 16, alpha = 0.6) +
  geom_smooth(method="lm", se=TRUE, lwd = 0) +
  geom_smooth(method="lm", se=FALSE, lwd = 2) +
  scale_color_manual(values = c("darkslategrey", "darkslategray1", "darkslategray4")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.75), axis.ticks = element_line(size = 0.75), legend.position='none', 
        axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),
        axis.text.y = element_text(size = 25), axis.title.y = element_text(size = 25)) +
  scale_y_continuous(breaks=c(0,1,2,3,4)) +
  scale_x_continuous(breaks=c(400,800,1200,1600)) +
  coord_cartesian(ylim = c(0.0, 4.0)) +
  coord_cartesian(xlim = c(200, 1850)) +
  xlab("Instrength (s)") +
  ylab("Number of Matings")
# Outdegree (Day 1) and Dilution effects on matings
ggplot(data = matingsdata_day1_dilution, aes(x = outdegree.seconds, y = matings, color = dilution.color)) + 
  geom_point(size = 3, shape = 16, alpha = 0.6) +
  geom_smooth(method="lm", se=TRUE, lwd = 0) +
  geom_smooth(method="lm", se=FALSE, lwd = 2) +
  scale_color_manual(values = c("darkslategrey", "darkslategray1", "darkslategray4")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.75), axis.ticks = element_line(size = 0.75), legend.position='none', 
        axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),
        axis.text.y = element_text(size = 25), axis.title.y = element_text(size = 25)) +
  scale_y_continuous(breaks=c(0,1,2,3,4)) +
  scale_x_continuous(breaks=c(400,800,1200,1600)) +
  coord_cartesian(ylim = c(0.0, 4.0)) +
  coord_cartesian(xlim = c(200, 1800)) +
  xlab("Outstrength (s)") +
  ylab("Number of Matings")
# Eigencentrality (Day 1) and Dilution effects on matings
ggplot(data = matingsdata_day1_dilution, aes(x = eigencentrality, y = matings, color = dilution.color)) + 
  geom_point(size = 3, shape = 16, alpha = 0.6) +
  geom_smooth(method="lm", se=TRUE, lwd = 0) +
  geom_smooth(method="lm", se=FALSE, lwd = 2) +
  scale_color_manual(values = c("darkslategrey", "darkslategray1", "darkslategray4")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.75), axis.ticks = element_line(size = 0.75), legend.position='none', 
        axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),
        axis.text.y = element_text(size = 25), axis.title.y = element_text(size = 25)) +
  scale_y_continuous(breaks=c(0,1,2,3,4)) +
  scale_x_continuous(breaks = c(0.1,0.2,0.3,0.4)) +
  coord_cartesian(ylim = c(0.0, 4.0)) +
  coord_cartesian(xlim = c(0.075,0.4)) +
  xlab("Eigenvector Centrality") +
  ylab("Number of Matings")

# EXTENDED DATA FIGURE 2 # -----------------
network_figure <- network %>%
  rename(
    "Instrength" = indegree.seconds,
    "Outstrength" = outdegree.seconds,
    "Clustering Coefficient" = clusteringcoeff,
    "Betweenness" = betweenness,
    "Eigenvector Centrality" = eigencentrality)
ggpairs(network_figure[, c("Instrength", "Outstrength", "Clustering Coefficient", "Betweenness", "Eigenvector Centrality")], mapping = aes(alpha = 1/10))
























