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

df7.8.9 <- df7.8.9[, c(taxa.keep)]
df7.8.9[is.na(df7.8.9)] <- 0 

df <- remove.zeroes(df7.8.9)

svg(file = "imagine_ordinated.svg", w = 12, h = 12, dpi = 300)
df2 <- OrderMatrix(df, scores =1, outputScores = FALSE)
metacom::Imagine(df2, col = c(0, "grey"),fill=FALSE,
                 xlab = "", ylab = "", speciesnames = FALSE,
                 sitenames = FALSE, order = TRUE)
#add taxa names
axis(3, at = c(1:ncol(df2)), labels = FALSE)
axis(side = 3, las = 2, mgp = c(3, 0.75, 0), labels = FALSE)
text(x = 1:ncol(df2), y = par("usr")[2] + 160, 
     labels = colnames(df2), 
     srt = -35, cex = 0.8, xpd = NA, adj = 1)
axis(side = 2, at = c(1:nrow(df2)), labels= FALSE)

#colour part of the labels 
all.sites <- rownames(df2) #ordinated station names

#klb9 
klb9.st <- data.klb9$station
ind.9 <- match(klb9.st, all.sites)
axis(side = 2, at = ind.9[1:length(ind.9)], labels = FALSE, 
     col.ticks = "orange", lwd = 2)
axis(side = 2, las = 2, mgp = c(3, 0.75, 0), labels = FALSE, 
     cex.axis = 0.5)
text(y = ind.9[1:length(ind.9)], x = par("usr")[1] - 1, 
     labels = klb9.st, xpd = NA, cex = 0.5, col = "orange")


#klb7
klb7.st <- data.klb7$station
ind.7 <- match(klb7.st, all.sites)
axis(side = 2, at = ind.7[1:length(ind.7)], labels = FALSE, 
     col.ticks = "forestgreen", lwd = 2)
axis(side = 2, las = 2, mgp = c(3, 0.75, 0), labels = FALSE, 
     cex.axis = 0.5)
text(y = ind.7[1:length(ind.7)], x = par("usr")[1] - 1, 
     labels = klb7.st, xpd = NA, cex = 0.5, col = "forestgreen")

#klb8
klb8.st <- data.klb8$station
ind.8 <- match(klb8.st, all.sites)
axis(side = 2, at = ind.8[1:length(ind.8)], labels = FALSE, 
     col.ticks = "blue", lwd = 2)
axis(side = 2, las = 2, mgp = c(3, 0.75, 0), labels = FALSE, 
     cex.axis = 0.5)
text(y = ind.8[1:length(ind.8)], x = par("usr")[1] - 1, 
     labels = klb8.st, xpd = NA, cex = 0.5, col = "blue")


legend("bottomright", legend = c("ORDINATED", "klb7", 
                                 "klb8", "klb9"),col = c("black", "forestgreen",
                                                         "blue", "orange"),
       bg = "white",lwd = 2, lty = c(0, 1, 1, 1))

dev.off()

# NMDS --------------------------------------------------------------------

data.klb7<-data.pri[data.pri$KLB==7,]
data.klb8 <- data.pri[data.pri$KLB==8,]
data.klb9 <- data.pri[data.pri$KLB==9,]

df7.8.9 <- rbind.fill(data.klb7, data.klb8, data.klb9)
rownames(df7.8.9) <- df7.8.9$station
mt.data <- colnames(df7.8.9[1:4])
df7.8.9 <- df7.8.9[, c(mt.data, taxa.keep)]

df7.8.9[is.na(df7.8.9)] <- 0 

#modify remove.zeroes
remove.zeroes <- function(df, data_cols = 5:ncol(df)) {
        # Separate metadata and data
        metadata <- df[, -data_cols]
        data <- df[, data_cols]

        # Remove taxa (columns) with all zeros
        col_sums <- colSums(data)
        data <- data[, col_sums > 0, drop = FALSE]

        # Remove rows (sites) with all zeros
        row_sums <- rowSums(data)
        keep_rows <- row_sums > 0
        data <- data[keep_rows, , drop = FALSE]
        metadata <- metadata[keep_rows, , drop = FALSE]

        # Convert to presence-absence
        data[data >= 1] <- 1

        # Recombine metadata and cleaned data
        cleaned_df <- cbind(metadata, data)

        return(cleaned_df)
}

df7.8.9.cleaned <- remove.zeroes(df7.8.9, data_cols = 5:ncol(df7.8.9))
pa_data <- df7.8.9.cleaned[5:ncol(df7.8.9.cleaned)]

dist_matrix <- vegdist(pa_data, method = "jaccard", binary = TRUE)

# 
nmds_result <- metaMDS(dist_matrix, k = 2, trymax = 100)
nmds_result

nmds_scores <- as.data.frame(scores(nmds_result, display = "sites"))

# Add the KLB group back into the NMDS scores
nmds_scores$KLB <- df7.8.9.cleaned$KLB
nmds_scores$KLB <- factor(nmds_scores$KLB, levels = c("9", "7", "8"))

# Plot using ggplot2
x <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = factor(KLB))) +
        geom_point(size = 3, alpha = 0.8) +
        scale_color_manual(values = c("7" = "forestgreen", "8" = "blue", "9" = "orange")) +
        labs(title = "NMDS Ordination (Jaccard, Presence/Absence)",
             x = "NMDS1", y = "NMDS2", color = "KLB Group") +
        theme_minimal() +
        xlim(-2, 2) + ylim(-2,2)
x

