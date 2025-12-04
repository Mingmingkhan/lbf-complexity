# Load Libraries -----

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


#read in KLB7 shapes Shapes
line2<-read.delim("data/Line2.txt")
line3<-read.delim("data/Line3.txt")

line2<-line2[,c(4,5)]
line3<-line3[,c(4,5)]

# make KLB 7 window --------
rev.line3 <-cbind(rev(line3[,1]),rev(line3[,2]))
colnames(rev.line3)<-colnames(line3[,1:2])

coords.klb7<-rbind(rev.line3,(line2[,1:2]))
coords.klb7.closed <- rbind(coords.klb7, line3[nrow(line3),])

# get stations using spatstat

win.klb7<-owin(poly=list(x=rev(coords.klb7$Long),y=rev(coords.klb7$Lat)))
ppp.klb7<-ppp(data.pri$Longitude, data.pri$Latitude,
              marks= data.pri$station, window = win.klb7)

klb7.st.ppp <- ppp.klb7$marks


#treatment 1 - move KLB 7 left ----

#move bottom of klb 
res<-move.line(line2,0.007) #have dist shows 579m
left.shift.bottom <- get.new.line(ori.line = line2, moved.line = res[[2]])

res<-move.line(line3, 0.007)
left.shift.top <- get.new.line(ori.line = line3, moved.line = res[[2]])


rev.left.shift.top<-cbind(rev(left.shift.top[,1]),rev(left.shift.top[,2]))
colnames(rev.left.shift.top)<-colnames(left.shift.top[,1:2])

coords.klb7.treat1 <-rbind(rev.left.shift.top,(left.shift.bottom[,1:2]))
coords.klb7.treat1.closed <- rbind(coords.klb7.treat1, left.shift.top[nrow(left.shift.top),])
coords.klb7.treat1 <- as.data.frame(coords.klb7.treat1)
colnames(coords.klb7.treat1) <- c("Long", "Lat")

win.klb7.treat1 <- owin(poly=list(x=rev(coords.klb7.treat1$Long),
                                  y=rev(coords.klb7.treat1$Lat)))
ppp.klb7.treat1 <- ppp(data.pri$Longitude, data.pri$Latitude,
                marks=data.pri$station, window = win.klb7.treat1)

klb7.st.treat1 <- ppp.klb7.treat1$marks

#treatment 2 - move KLB 7 right ----

#move bottom of klb 
res<-move.line(line2,0.007) #have dist shows 579m
right.shift.bottom <- get.new.line(ori.line = line2, moved.line = res[[1]])

res<-move.line(line3, 0.007)
right.shift.top <- get.new.line(ori.line = line3, moved.line = res[[1]])


rev.right.shift.top<-cbind(rev(right.shift.top[,1]),rev(right.shift.top[,2]))
colnames(rev.right.shift.top)<-colnames(right.shift.top[,1:2])

coords.klb7.treat2 <-rbind(rev.right.shift.top,(right.shift.bottom[,1:2]))
coords.klb7.treat2.closed <- rbind(coords.klb7.treat2, right.shift.top[nrow(right.shift.top),])

coords.klb7.treat2 <- as.data.frame(coords.klb7.treat2)
colnames(coords.klb7.treat2) <- c("Long", "Lat")

win.klb7.treat2 <- owin(poly=list(x=rev(coords.klb7.treat2$Long),
                                  y=rev(coords.klb7.treat2$Lat)))
ppp.klb7.treat2 <- ppp(data.pri$Longitude, data.pri$Latitude,
                       marks= data.pri$station, window = win.klb7.treat2)

klb7.st.treat2 <- ppp.klb7.treat2$marks


#Treatment 3 - Contract KLB ----
#need right shift bottom, left shift top 

rev.left.shift.top<-cbind(rev(left.shift.top[,1]),rev(left.shift.top[,2]))
colnames(rev.left.shift.top)<-colnames(left.shift.top[,1:2])

coords.klb7.treat3 <-rbind(rev.left.shift.top,(right.shift.bottom[,1:2]))
coords.klb7.treat3.closed <- rbind(coords.klb7.treat3, left.shift.top[nrow(left.shift.top),])

coords.klb7.treat3 <- as.data.frame(coords.klb7.treat3)
colnames(coords.klb7.treat3) <- c("Long", "Lat")

win.klb7.treat3 <- owin(poly=list(x=rev(coords.klb7.treat3$Long),
                                  y=rev(coords.klb7.treat3$Lat)))
ppp.klb7.treat3 <- ppp(data.pri$Longitude, data.pri$Latitude,
                       marks= data.pri$station, window = win.klb7.treat3)

klb7.st.treat3 <- ppp.klb7.treat3$marks


#plot the treatments (Figure 2) -------

klb7.col <- rgb(0, 0.65098, 0.313725, alpha = 0.2)

par(mfrow = c(1,3))
plot(coords.klb7, type = "l", 
     xlim = c(-56.83, -56.73), ylim = c(-64.33, -64.24),
     lwd = 1)
grid()
title("(a) Shift KLB west")
polygon(coords.klb7.closed, col = "white")
points(data.klb7$Longitude, data.klb7$Latitude, 
       col = "black", pch = 1)
polygon(coords.klb7.treat1.closed, col = klb7.col, lty = 2, lwd = 2)
points(ppp.klb7.treat1$x, ppp.klb7.treat1$y, 
       col = "black", pch = 16)

plot(coords.klb7, type = "l", 
     xlim = c(-56.83, -56.73), ylim = c(-64.33, -64.24),
     lwd = 1)
grid()
title("(b) Shift KLB east")
polygon(coords.klb7.closed, col = "white")
points(data.klb7$Longitude, data.klb7$Latitude, 
       col = "black", pch = 1)
polygon(coords.klb7.treat2.closed, col = klb7.col, lty = 2, lwd = 2)
points(ppp.klb7.treat2$x, ppp.klb7.treat2$y, 
       col = "black", pch = 16)

plot(coords.klb7, type = "l", 
     xlim = c(-56.83, -56.73), ylim = c(-64.33, -64.24),
     lwd = 1)
grid()
title("(c) Contract KLB")
polygon(coords.klb7.closed, col = "white")
points(data.klb7$Longitude, data.klb7$Latitude, 
       col = "black", pch = 1)
polygon(coords.klb7.treat3.closed, col = klb7.col, lty = 2, lwd = 2)
points(ppp.klb7.treat3$x, ppp.klb7.treat3$y, 
       col = "black", pch = 16)

# create sensitivity data frames ------------------------------------------

#raw data 
data.klb7 <- df7.8.9[df7.8.9$KLB == 7,]
rownames(data.klb7) <-data.klb7$station
data.klb7 <- data.klb7[, colnames(data.klb7) %in% taxa.keep]
data.klb7 <- remove.zeroes(data.klb7)

#treatment 1 - left shift
data.klb7.treat1 <- data.pri[data.pri$station %in% klb7.st.treat1,]
data.klb7.treat1 <- data.klb7.treat1[, colnames(data.klb7.treat1) %in% taxa.keep]
data.klb7.treat1 <- remove.zeroes(data.klb7.treat1)

#treatment 2 - right shift 
data.klb7.treat2 <- data.pri[data.pri$station %in% klb7.st.treat2,]
data.klb7.treat2 <- data.klb7.treat2[, colnames(data.klb7.treat2) %in% taxa.keep]
data.klb7.treat2 <- remove.zeroes(data.klb7.treat2)

#treatment 3 - contract
data.klb7.treat3 <- data.pri[data.pri$station %in% klb7.st.treat3,]
data.klb7.treat3 <- data.klb7.treat3[, colnames(data.klb7.treat3) %in% taxa.keep]
data.klb7.treat3 <- remove.zeroes(data.klb7.treat3)


# ANALYSIS ----------------------------------------------------------------

#make list of KLB 7 res -------

klb7.metacom.res <- vector(mode = "list", length = 4)
names(klb7.metacom.res) <- c("raw", "treat1", "treat2", 
                             "treat3")

klb7.metacom.res$raw <- get.metacom.res(data.klb7)
klb7.metacom.res$treat1 <- get.metacom.res(data.klb7.treat1)
klb7.metacom.res$treat2 <- get.metacom.res(data.klb7.treat2)
klb7.metacom.res$treat3 <- get.metacom.res(data.klb7.treat3)

klb7.metacom.res.df <- do.call(cbind, klb7.metacom.res)
klb7.metacom.res.df2 <- do.call(rbind, klb7.metacom.res)

#write.csv(klb7.metacom.res.df, 
#          file = "results/klb7_metacom_res.csv",
#          row.names = FALSE)

klb7.metacom.outputs <- metacom.output(klb7.metacom.res)

save(klb7.metacom.res, klb7.metacom.outputs, klb7.metacom.res.df,
     file = "results/klb7.metacom.results.RData")

# do co-occur -------------------------------------------------------------

klb7.cooc.res <- vector("list", 4)
names(klb7.cooc.res) <- c("raw", "left.shift",
                          "right.shift",
                          "contract")
klb7.df <- data.klb7
klb7.df <- remove.sing(klb7.df)

klb7.df <- as.data.frame(t(klb7.df)) #transpose the matrix


klb7.mat <- cooccur(klb7.df, type = "spp_site", spp_names = TRUE, 
                    thresh = TRUE)
summary(klb7.mat)
plot(klb7.mat)

klb7.cooc.res[[1]] <- klb7.mat

data.klb7.treat1 <- remove.sing(data.klb7.treat1)
data.klb7.treat1 <- as.data.frame(t(data.klb7.treat1)) #transpose the matrix
klb7.treat1.mat <- cooccur(data.klb7.treat1, type = "spp_site", spp_names = TRUE, 
                           thresh = TRUE)
plot(klb7.treat1.mat)
summary(klb7.treat1.mat)
klb7.cooc.res[[2]] <- klb7.treat1.mat


data.klb7.treat2 <- as.data.frame(t(data.klb7.treat2)) #transpose the matrix
klb7.treat2.mat <- cooccur(data.klb7.treat2, type = "spp_site", spp_names = TRUE, 
                           thresh = TRUE)
plot(klb7.treat2.mat)
summary(klb7.treat2.mat)
klb7.cooc.res[[3]] <- klb7.treat2.mat


data.klb7.treat3 <- as.data.frame(t(data.klb7.treat3)) #transpose the matrix
klb7.treat3.mat <- cooccur(data.klb7.treat3, type = "spp_site", spp_names = TRUE, 
                           thresh = TRUE)
plot(klb7.treat3.mat)
summary(klb7.treat3.mat)
klb7.cooc.res[[4]] <- klb7.treat3.mat


save(klb7.cooc.res, file = "results/klb7.cooccur.results.RData")



  
  
  