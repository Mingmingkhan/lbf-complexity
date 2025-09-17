library(spatstat)
library(dplyr)
library(plyr)
library(metacom)
library(cooccur)
source("code/palaeoFunctions.R")

#feed in data and allocate to KLB according to zones

data.pri<-read.csv(file = "data/Seymour_LBF_Fossils.csv")

data.klb7<-data.pri[data.pri$KLB==7,]
data.klb8 <- data.pri[data.pri$KLB==8,]
data.klb9 <- data.pri[data.pri$KLB==9,]

df7.8.9 <- rbind.fill(data.klb7, data.klb8, data.klb9)
rownames(df7.8.9) <- df7.8.9$station
mt.data <- colnames(df7.8.9[1:4])

#Taxa for analyses 
taxa.keep <- readLines("data/taxa.keep.txt")

df7.8.9 <- df7.8.9[, c(mt.data, taxa.keep)]
df7.8.9[is.na(df7.8.9)] <- 0 

#read in line Shapes
line3<-read.delim("data/Line3.txt")
line4<-read.delim("data/Line4.txt")

line3<-line3[,c(4,5)]
line4<-line4[,c(4,5)]

#make KLB 8 window --------
rev.line4 <-cbind(rev(line4[,1]),rev(line4[,2]))
colnames(rev.line4)<-colnames(line4[,1:2])

coords.klb8<-rbind(rev.line4,(line3[,1:2]))
coords.klb8.closed <- rbind(coords.klb8, line4[nrow(line4),])

# get stations using spatstat

win.klb8<-owin(poly=list(x=rev(coords.klb8$Long),y=rev(coords.klb8$Lat)))
ppp.klb8<-ppp(data.pri$Longitude, data.pri$Latitude,
              marks= data.pri$station, window = win.klb8)

klb8.st.ppp <- ppp.klb8$marks


#treatment 1 - move KLB 8 left ----

#move bottom of klb 
res<-move.line(line3,0.007) 
left.shift.bottom <- get.new.line(ori.line = line3, moved.line = res[[2]])

res<-move.line(line4, 0.007)
left.shift.top <- get.new.line(ori.line = line4, moved.line = res[[2]])


rev.left.shift.top<-cbind(rev(left.shift.top[,1]),rev(left.shift.top[,2]))
colnames(rev.left.shift.top)<-colnames(left.shift.top[,1:2])

coords.klb8.treat1 <-rbind(rev.left.shift.top,(left.shift.bottom[,1:2]))
coords.klb8.treat1.closed <- rbind(coords.klb8.treat1, left.shift.top[nrow(left.shift.top),])
coords.klb8.treat1 <- as.data.frame(coords.klb8.treat1)
colnames(coords.klb8.treat1) <- c("Long", "Lat")

win.klb8.treat1 <- owin(poly=list(x=rev(coords.klb8.treat1$Long),
                                  y=rev(coords.klb8.treat1$Lat)))
ppp.klb8.treat1 <- ppp(data.pri$Longitude, data.pri$Latitude,
                       marks= data.pri$station, window = win.klb8.treat1)

klb8.st.treat1 <- ppp.klb8.treat1$marks

#treatment 2 - move KLB 7 right ----

#move bottom of klb 
res<-move.line(line3,0.007) #have dist shows 579m
right.shift.bottom <- get.new.line(ori.line = line3, moved.line = res[[1]])

res<-move.line(line4, 0.007)
right.shift.top <- get.new.line(ori.line = line4, moved.line = res[[1]])


rev.right.shift.top<-cbind(rev(right.shift.top[,1]),rev(right.shift.top[,2]))
colnames(rev.right.shift.top)<-colnames(right.shift.top[,1:2])

coords.klb8.treat2 <-rbind(rev.right.shift.top,(right.shift.bottom[,1:2]))
coords.klb8.treat2.closed <- rbind(coords.klb8.treat2, right.shift.top[nrow(right.shift.top),])

coords.klb8.treat2 <- as.data.frame(coords.klb8.treat2)
colnames(coords.klb8.treat2) <- c("Long", "Lat")

win.klb8.treat2 <- owin(poly=list(x=rev(coords.klb8.treat2$Long),
                                  y=rev(coords.klb8.treat2$Lat)))
ppp.klb8.treat2 <- ppp(data.pri$Longitude, data.pri$Latitude,
                       marks= data.pri$station, window = win.klb8.treat2)

klb8.st.treat2 <- ppp.klb8.treat2$marks

#Treatment 3 - Contract KLB -----

#need right shift bottom, left shift top 
rev.left.shift.top<-cbind(rev(left.shift.top[,1]),rev(left.shift.top[,2]))
colnames(rev.left.shift.top)<-colnames(left.shift.top[,1:2])

coords.klb8.treat3 <-rbind(rev.left.shift.top,(right.shift.bottom[,1:2]))
coords.klb8.treat3.closed <- rbind(coords.klb8.treat3, left.shift.top[nrow(left.shift.top),])

coords.klb8.treat3 <- as.data.frame(coords.klb8.treat3)
colnames(coords.klb8.treat3) <- c("Long", "Lat")

win.klb8.treat3 <- owin(poly=list(x=rev(coords.klb8.treat3$Long),
                                  y=rev(coords.klb8.treat3$Lat)))
ppp.klb8.treat3 <- ppp(data.pri$Longitude, data.pri$Latitude,
                       marks= data.pri$station, window = win.klb8.treat3)

klb8.st.treat3 <- ppp.klb8.treat3$marks

# create sensitivity data frames ------------------------------------------

#raw data 
data.klb8 <- df7.8.9[df7.8.9$KLB == 8,]
rownames(data.klb8) <-data.klb8$station
data.klb8 <- data.klb8[, colnames(data.klb8) %in% taxa.keep]
data.klb8 <- remove.zeroes(data.klb8)

#treatment 1 - left shift
data.klb8.treat1 <- data.pri[data.pri$station %in% klb8.st.treat1,]
data.klb8.treat1 <- data.klb8.treat1[, colnames(data.klb8.treat1) %in% taxa.keep]
data.klb8.treat1 <- remove.zeroes(data.klb8.treat1)

#treatment 2 - right shift 
data.klb8.treat2 <- data.pri[data.pri$station %in% klb8.st.treat2,]
data.klb8.treat2 <- data.klb8.treat2[, colnames(data.klb8.treat2) %in% taxa.keep]
data.klb8.treat2 <- remove.zeroes(data.klb8.treat2)

#treatment 3 - contract
data.klb8.treat3 <- data.pri[data.pri$station %in% klb8.st.treat3,]
data.klb8.treat3 <- data.klb8.treat3[, colnames(data.klb8.treat3) %in% taxa.keep]
data.klb8.treat3 <- remove.zeroes(data.klb8.treat3)

# ANALYSES -------

#make list of KLB 8 res 

klb8.metacom.res <- vector(mode = "list", length = 4)
names(klb8.metacom.res) <- c("raw", "treat1", "treat2", 
                             "treat3")

klb8.metacom.res$raw <- get.metacom.res(data.klb8)
klb8.metacom.res$treat1 <- get.metacom.res(data.klb8.treat1)
klb8.metacom.res$treat2 <- get.metacom.res(data.klb8.treat2)
klb8.metacom.res$treat3 <- get.metacom.res(data.klb8.treat3)

klb8.metacom.res.df <- do.call(cbind, klb8.metacom.res)
klb8.metacom.res.df2 <- do.call(rbind, klb8.metacom.res)

#write.csv(klb8.metacom.res.df, 
#          file = "results/klb8_metacom_res.csv",
#          row.names = FALSE)

klb8.metacom.outputs <- metacom.output(klb8.metacom.res)

save(klb8.metacom.res, klb8.metacom.outputs, klb8.metacom.res.df,
     file = "results/klb8.metacom.results.RData")

# do co-occur -------------------------------------------------------------

library(cooccur)

klb8.cooc.res <- vector("list", length = 4)
names(klb8.cooc.res) <- c("raw", "left.shift",
                          "right.shift",
                          "contract")

data.klb8 <- remove.sing(data.klb8)
data.klb8 <- as.data.frame(t(data.klb8)) #transpose the matrix

klb8.mat <- cooccur(data.klb8, type = "spp_site", spp_names = TRUE, 
                    thresh = TRUE)
summary(klb8.mat)

klb8.cooc.res[[1]] <- klb8.mat

summary(klb8.mat)
plot(klb8.mat)

data.klb8.treat1 <- remove.sing(data.klb8.treat1)
data.klb8.treat1 <- as.data.frame(t(data.klb8.treat1)) #transpose the matrix
klb8.treat1.mat <- cooccur(data.klb8.treat1, type = "spp_site", spp_names = TRUE, 
                           thresh = TRUE)
plot(klb8.treat1.mat)
summary(klb8.treat1.mat)
klb8.cooc.res[[2]] <- klb8.treat1.mat

data.klb8.treat2 <- remove.sing(data.klb8.treat2)
data.klb8.treat2 <- as.data.frame(t(data.klb8.treat2)) #transpose the matrix
klb8.treat2.mat <- cooccur(data.klb8.treat2, type = "spp_site", spp_names = TRUE, 
                           thresh = TRUE)
plot(klb8.treat2.mat)
summary(klb8.treat2.mat)
klb8.cooc.res[[3]] <- klb8.treat2.mat


data.klb8.treat3 <- remove.sing(data.klb8.treat3)
data.klb8.treat3 <- as.data.frame(t(data.klb8.treat3)) #transpose the matrix
klb8.treat3.mat <- cooccur(data.klb8.treat3, type = "spp_site", spp_names = TRUE, 
                           thresh = TRUE)
plot(klb8.treat3.mat)
summary(klb8.treat3.mat)
klb8.cooc.res[[4]] <- klb8.treat3.mat

save(klb8.cooc.res, file = "results/klb8.cooccur.results.RData")




