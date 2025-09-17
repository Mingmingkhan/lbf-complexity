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


#read in Line Shapes
line4<-read.delim("data/Line4.txt")
line5<-read.delim("data/Line5.txt")


#make KLB 9 window --------
rev.line5 <-cbind(rev(line5[,1]),rev(line5[,2]))
colnames(rev.line5)<-colnames(line5[,1:2])

coords.klb9<-rbind(rev.line5,(line4[,1:2]))
coords.klb9.closed <- rbind(coords.klb9, line5[nrow(line5),])

plot(coords.klb9, type = "l")
# get stations using spatstat

win.klb9<-owin(poly=list(x=rev(coords.klb9$Long),y=rev(coords.klb9$Lat)))
ppp.klb9<-ppp(data.pri$Longitude, data.pri$Latitude,
              marks= data.pri$station, window = win.klb9)

klb9.st.ppp <- ppp.klb9$marks

#treatment 2 - contract KLB ----

#move bottom of klb 
res<-move.line(line4,0.007) #have dist shows 579m
right.shift.bottom <- get.new.line(ori.line = line4, moved.line = res[[1]])

coords.klb9.treat2 <-rbind(rev.line5,(right.shift.bottom[,1:2]))
coords.klb9.treat2.closed <- rbind(coords.klb9.treat2, line5[nrow(line5),])

coords.klb9.treat2 <- as.data.frame(coords.klb9.treat2)
colnames(coords.klb9.treat2) <- c("Long", "Lat")

win.klb9.treat2 <- owin(poly=list(x=rev(coords.klb9.treat2$Long),
                                  y=rev(coords.klb9.treat2$Lat)))
ppp.klb9.treat2 <- ppp(data.pri$Longitude, data.pri$Latitude,
                       marks= data.pri$station, window = win.klb9.treat2)

klb9.st.treat2 <- ppp.klb9.treat2$marks

klb9.st.treat1 <- klb9.st.treat2


# create sensitivity data frames ------------------------------------------

#raw data 
data.klb9 <- df7.8.9[df7.8.9$KLB == 9,]
rownames(data.klb9) <-data.klb9$station
data.klb9 <- data.klb9[, colnames(data.klb9) %in% taxa.keep]
data.klb9 <- remove.zeroes(data.klb9)

#treatment 1 - contract KLB
data.klb9.treat1 <- data.pri[data.pri$station %in% klb9.st.treat1,]
data.klb9.treat1 <- data.klb9.treat1[, colnames(data.klb9.treat1) %in% taxa.keep]
data.klb9.treat1 <- remove.zeroes(data.klb9.treat1)


# ANALYSES ------

#res list ----

klb9.metacom.res <- vector(mode = "list", length = 2)
names(klb9.metacom.res) <- c("raw", "treat1")

klb9.metacom.res$raw <- get.metacom.res(data.klb9)
klb9.metacom.res$treat1 <- get.metacom.res(data.klb9.treat1)


klb9.metacom.res.df <- do.call(cbind, klb9.metacom.res)
klb9.metacom.res.df2 <- do.call(rbind, klb9.metacom.res)

#write.csv(klb9.metacom.res.df, 
#          file = "results/klb9_metacom_res.csv",
#          row.names = FALSE)

klb9.metacom.outputs <- metacom.output(klb9.metacom.res)

save(klb9.metacom.res, klb9.metacom.outputs, klb9.metacom.res.df,
     file = "results/klb9.metacom.results.RData")


# do co-occur -------------------------------------------------------------

library(cooccur)
klb9.cooc.res <- vector("list", length = 2)
names(klb9.cooc.res) <- c("raw", "contract")

data.klb9 <- remove.sing(data.klb9)

data.klb9 <- as.data.frame(t(data.klb9)) #transpose the matrix

klb9.mat <- cooccur(data.klb9, type = "spp_site", spp_names = TRUE, 
                    thresh = TRUE)
summary(klb9.mat)
plot(klb9.mat)

klb9.cooc.res[[1]] <- klb9.mat

data.klb9.treat1 <- remove.sing(data.klb9.treat1)
data.klb9.treat1 <- as.data.frame(t(data.klb9.treat1)) #transpose the matrix
klb9.treat1.mat <- cooccur(data.klb9.treat1, type = "spp_site", spp_names = TRUE, 
                           thresh = TRUE)
plot(klb9.treat1.mat)
summary(klb9.treat1.mat)
klb9.cooc.res[[2]] <- klb9.treat1.mat

save(klb9.cooc.res, file = "results/klb9.cooccur.res.RData")



