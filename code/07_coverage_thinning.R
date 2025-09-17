source("code/palaeoFunctions.R")
library(dplyr)
library(plyr)

modified.remove.zeroes <- function(df) {
  # Split metadata and taxa based on externally defined colnames
  meta_df <- df[, mt.data, drop = FALSE]
  taxa_df <- df[, taxa.keep, drop = FALSE]
  
  # Remove zero-only columns (taxa)
  taxa_matrix <- as.matrix(taxa_df)
  col_sums <- colSums(taxa_matrix, na.rm = TRUE)
  taxa_matrix <- taxa_matrix[, col_sums > 0, drop = FALSE]
  
  # Remove zero-only rows (samples)
  row_sums <- rowSums(taxa_matrix, na.rm = TRUE)
  taxa_matrix <- taxa_matrix[row_sums > 0, , drop = FALSE]
  
  # Convert back to data frame
  taxa_df_clean <- as.data.frame(taxa_matrix)
  
  # Subset metadata to match remaining rows
  meta_df_clean <- meta_df[rownames(taxa_df_clean), , drop = FALSE]
  
  # Combine metadata and cleaned taxa data
  final_df <- cbind(meta_df_clean, taxa_df_clean)
  
  return(final_df)
}

remove.sing <- function(df) {
  # Safe taxa selection
  valid_taxa <- intersect(taxa.keep, colnames(df))
  valid_meta <- intersect(mt.data, colnames(df))
  
  # Subset metadata and taxa
  meta_df <- df[, valid_meta, drop = FALSE]
  taxa_df <- df[, valid_taxa, drop = FALSE]
  
  # Convert to presence/absence
  taxa_df[taxa_df >= 1] <- 1
  
  # Remove singleton taxa
  col_sums <- colSums(taxa_df, na.rm = TRUE)
  taxa_df <- taxa_df[, col_sums > 1, drop = FALSE]
  
  # Remove samples with no remaining taxa
  row_sums <- rowSums(taxa_df, na.rm = TRUE)
  taxa_df <- taxa_df[row_sums > 0, , drop = FALSE]
  
  # Subset metadata to remaining samples
  meta_df <- meta_df[rownames(taxa_df), , drop = FALSE]
  
  # Combine
  cleaned_df <- cbind(meta_df, taxa_df)
  return(cleaned_df)
}

# Supplementary Figure 1: Rarefaction curves ----


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

data.klb7 <- df7.8.9[df7.8.9$KLB == 7,]
data.klb8 <- df7.8.9[df7.8.9$KLB == 8,]
data.klb9 <- df7.8.9[df7.8.9$KLB == 9,]

data.789 <- modified.remove.zeroes(df7.8.9)
data.7 <- modified.remove.zeroes(data.klb7)
data.8 <- modified.remove.zeroes(data.klb8)
data.9 <- modified.remove.zeroes(data.klb9)

#do the rarefaction (Fig s1)
library(vegan)

sp.789 <- specaccum(data.789[,colnames(data.789) %in% taxa.keep], method = "rarefaction")
sp.7 <- specaccum(data.7[,colnames(data.7) %in% taxa.keep], method = "rarefaction")

sp.8 <- specaccum(data.8[,colnames(data.8) %in% taxa.keep], method = "rarefaction")
sp.9 <- specaccum(data.9[,colnames(data.9) %in% taxa.keep], method = "rarefaction")
sl <- diff(sp.9$richness) / diff(sp.9$sites)
last_sl <- tail(sl, 1)

#plot 

plot(sp.789, main = "Accumulation curves", xlab = "number of sites", 
     ylab = "observed richness (number of families)")
grid(col = "lightgray")
plot(sp.9, add = TRUE, col = "orange", lwd = 1.5)
plot(sp.7, add = TRUE, col = "cornflowerblue", lwd = 1.5)
plot(sp.8, add = TRUE, col = "forestgreen", lwd = 1.5)

#coverage rarefaction 
library(iNEXT)
library(ggplot2)

df7 <- data.7[,colnames(data.7) %in% taxa.keep]
df8 <- data.8[,colnames(data.8) %in% taxa.keep]
df9 <- data.9[,colnames(data.9) %in% taxa.keep]

#Sum abundances for each unit
abund_KLB7 <- colSums(df7)
abund_KLB8 <- colSums(df8)
abund_KLB9 <- colSums(df9)


abund.list <- list(KLB7 = abund_KLB7, KLB8 = abund_KLB8, KLB9 = abund_KLB9)
sp.789.out <- iNEXT(abund.list, q = 0, datatype = "abundance", knots = 40, se = TRUE)
ggiNEXT(sp.789.out, type = 2) + ggtitle("Sample Coverage")

# Align species columns across units
# Get union of all species
all_species <- union(names(abund_KLB7), names(abund_KLB8)) %>%
  union(names(abund_KLB9))

# Function to fill in missing species with 0s
fill_missing_species <- function(abund, all_species) {
  missing <- setdiff(all_species, names(abund))
  abund[missing] <- 0
  abund <- abund[all_species]  # reorder to keep consistent order
  return(abund)
}

abund_KLB7 <- fill_missing_species(abund_KLB7, all_species)
abund_KLB8 <- fill_missing_species(abund_KLB8, all_species)
abund_KLB9 <- fill_missing_species(abund_KLB9, all_species)

# Create abundance list for iNEXT
abund_list <- list(
  KLB7 = abund_KLB7,
  KLB8 = abund_KLB8,
  KLB9 = abund_KLB9
)

# Run iNEXT at 95% sample coverage
estimates_95 <- estimateD(
  x = abund_list,
  datatype = "abundance",
  q = 0,                 # species richness
  base = "coverage",
  level = 0.95           # desired coverage level
)

# View results
print(estimates_95)

out <- iNEXT(
  x = abund_list,
  q = 0,
  datatype = "abundance"
)

cov.plot <- ggiNEXT(out, type = 3) + ggtitle("Coverage-Based Rarefaction") + theme_light() +
  theme(legend.position = "bottom") + labs(y = "Number of Families") + ylim(0, 50) +
  scale_color_manual(values = c("KLB9" = "orange",
                                "KLB7" = "cornflowerblue",
                                "KLB8" = "forestgreen")) +
  scale_fill_manual(values = c("KLB9" = "orange",
                               "KLB7" = "cornflowerblue",
                               "KLB8" = "forestgreen")) + 
  geom_vline(xintercept = 0.95, linetype = 3, color = "black")
cov.plot


# spatial thinning -------------------------------------------------------

modified.remove.zeroes <- function(df) {
  # Split metadata and taxa based on externally defined colnames
  meta_df <- df[, mt.data, drop = FALSE]
  taxa_df <- df[, taxa.keep, drop = FALSE]
  #turn to present/absence
  taxa_df[taxa_df >= 1] <- 1
  
  # Remove zero-only columns (taxa)
  taxa_matrix <- as.matrix(taxa_df)
  col_sums <- colSums(taxa_matrix, na.rm = TRUE)
  taxa_matrix <- taxa_matrix[, col_sums > 0, drop = FALSE]
  
  # Remove zero-only rows (samples)
  row_sums <- rowSums(taxa_matrix, na.rm = TRUE)
  taxa_matrix <- taxa_matrix[row_sums > 0, , drop = FALSE]
  
  # Convert back to data frame
  taxa_df_clean <- as.data.frame(taxa_matrix)
  
  # Subset metadata to match remaining rows
  meta_df_clean <- meta_df[rownames(taxa_df_clean), , drop = FALSE]
  
  # Combine metadata and cleaned taxa data
  final_df <- cbind(meta_df_clean, taxa_df_clean)
  
  return(final_df)
}



#first KLB 7 ----

klb7.st <- data.7[,colnames(data.7) %in% mt.data]

#turn coords into an sf object 
library(sf)
library(tidysdm)

stations <- rownames(data.7)
stations <- as.numeric(stations)
klb7.st <- klb7.st[klb7.st$station %in% stations,]

data_coords <- st_as_sf(klb7.st, coords = c("Longitude", "Latitude"))
st_crs(data_coords) <- 4326

#klb7 thinning 30 ----

klb7_thinning30 <- vector(mode = "list", length = 10)
for (i in 1:10){
  names(klb7_thinning30)[i] <- paste("trial_",i)
  klb7_thinning30[[i]]$thinned_st <- thin_by_dist(data_coords, 30)
  klb7_thinning30[[i]]$thinned_dat <- data.klb7[data.klb7$station %in% 
                                klb7_thinning30[[i]]$thinned_st$station,]
  klb7_thinning30[[i]]$metacom_dat <- modified.remove.zeroes(klb7_thinning30[[i]]$thinned_dat)
  
  #do metacom
  klb7_thinning30[[i]]$coh.value <- 
    metacom::Coherence(klb7_thinning30[[i]]$metacom_dat
                       [,colnames(klb7_thinning30[[i]]$metacom_dat) %in% taxa.keep])
  klb7_thinning30[[i]]$coh.z <- klb7_thinning30[[i]]$coh.value[2, 2]
  klb7_thinning30[[i]]$turn.value <- 
    metacom::Turnover(klb7_thinning30[[i]]$metacom_dat
                      [,colnames(klb7_thinning30[[i]]$metacom_dat) %in% taxa.keep])
    
  klb7_thinning30[[i]]$turn.z <- klb7_thinning30[[i]]$turn.value[2,2]
  klb7_thinning30[[i]]$turn.p <- klb7_thinning30[[i]]$turn.value[3,2]
  
  klb7_thinning30[[i]]$bc.value <- 
    metacom::BoundaryClump(klb7_thinning30[[i]]$metacom_dat
                       [,colnames(klb7_thinning30[[i]]$metacom_dat) %in% taxa.keep])
  klb7_thinning30[[i]]$bc.i <- klb7_thinning30[[i]]$bc.value[1,2]
  klb7_thinning30[[i]]$bc.p <- klb7_thinning30[[i]]$bc.value[2,2]
  
  
  #do cooccur
  klb7_thinning30[[i]]$cooccur_dat <- remove.sing(klb7_thinning30[[i]]$metacom_dat)
  cooccur_dat <- klb7_thinning30[[i]]$cooccur_dat
  cooccur_dat <- cooccur_dat[,colnames(cooccur_dat) %in% taxa.keep]
  
  temp2 <- as.data.frame(t(cooccur_dat))
  klb7_thinning30[[i]]$cooccur_dat <- temp2
  cooc.mat <- cooccur::cooccur(temp2, type = "spp_site", spp_names = TRUE,
                               thresh = TRUE)
  temp.sum <- summary(cooc.mat)
  temp.res <- temp.sum[7]
  klb7_thinning30[[i]]$cooc.mat <- cooc.mat
  klb7_thinning30[[i]]$cooc.percent <- temp.res
  cat("trial completed", i, "\n")
}


klb7_thinned.site.number <- lapply(klb7_thinning30, `[[`,"thinned_st")
x <- lapply(klb7_thinned.site.number, `[[`,1)
klb7_thinned.site.number <- matrix(sapply(x, length), ncol = 1)
rownames(klb7_thinned.site.number) <- names(x)

coh.res.all <- lapply(klb7_thinning30, `[[`,"coh.z")
coh.res.all <- do.call(rbind, coh.res.all)

turn.res.all <- lapply(klb7_thinning30, `[[`,"turn.z")
turn.res.all <- do.call(rbind, turn.res.all)

turn.p.all <- lapply(klb7_thinning30, `[[`,"turn.p")
turn.p.all <- do.call(rbind, turn.p.all)

bc.i.all <- lapply(klb7_thinning30, `[[`,"bc.i")
bc.i.all <- do.call(rbind, bc.i.all)

bc.p.all <- lapply(klb7_thinning30, `[[`,"bc.p")
bc.p.all <- do.call(rbind, bc.p.all)

cooc.res.all <- lapply(klb7_thinning30, `[[`,"cooc.percent")
cooc.res.all <- do.call(rbind, cooc.res.all)
klb7_thinned_summary_res30 <- as.data.frame(cbind(klb7_thinned.site.number, 
                                                  coh.res.all, 
                                                  turn.res.all, 
                                                  turn.p.all,
                                                  bc.i.all, 
                                                  bc.p.all,
                                                  cooc.res.all))
colnames(klb7_thinned_summary_res30) <- c("sampled_sites", "coherence",
                                          "turnover", "turnover.p",
                                          "bc.i", "bc.p", "cooccur")

# klb7 thinning 50 --------------------------------------------------------


klb7_thinning50 <- vector(mode = "list", length = 10)
for (i in 1:10){
  names(klb7_thinning50)[i] <- paste("trial_",i)
  klb7_thinning50[[i]]$thinned_st <- thin_by_dist(data_coords, 50)
  klb7_thinning50[[i]]$thinned_dat <- data.klb7[data.klb7$station %in% 
                                                  klb7_thinning50[[i]]$thinned_st$station,]
  klb7_thinning50[[i]]$metacom_dat <- modified.remove.zeroes(klb7_thinning50[[i]]$thinned_dat)
  
  #do metacom
  klb7_thinning50[[i]]$coh.value <- 
    metacom::Coherence(klb7_thinning50[[i]]$metacom_dat
                       [,colnames(klb7_thinning50[[i]]$metacom_dat) %in% taxa.keep])
  klb7_thinning50[[i]]$coh.z <- klb7_thinning50[[i]]$coh.value[2, 2]
  klb7_thinning50[[i]]$turn.value <- 
    metacom::Turnover(klb7_thinning50[[i]]$metacom_dat
                      [,colnames(klb7_thinning50[[i]]$metacom_dat) %in% taxa.keep])
  
  klb7_thinning50[[i]]$turn.z <- klb7_thinning50[[i]]$turn.value[2,2]
  klb7_thinning50[[i]]$turn.p <- klb7_thinning50[[i]]$turn.value[3,2]
  
  klb7_thinning50[[i]]$bc.value <- 
    metacom::BoundaryClump(klb7_thinning50[[i]]$metacom_dat
                           [,colnames(klb7_thinning50[[i]]$metacom_dat) %in% taxa.keep])
  klb7_thinning50[[i]]$bc.i <- klb7_thinning50[[i]]$bc.value[1,2]
  klb7_thinning50[[i]]$bc.p <- klb7_thinning50[[i]]$bc.value[2,2]
  
  
  #do cooccur
  klb7_thinning50[[i]]$cooccur_dat <- remove.sing(klb7_thinning50[[i]]$metacom_dat)
  cooccur_dat <- klb7_thinning50[[i]]$cooccur_dat
  cooccur_dat <- cooccur_dat[,colnames(cooccur_dat) %in% taxa.keep]
  
  temp2 <- as.data.frame(t(cooccur_dat))
  klb7_thinning50[[i]]$cooccur_dat <- temp2
  cooc.mat <- cooccur::cooccur(temp2, type = "spp_site", spp_names = TRUE,
                               thresh = TRUE)
  temp.sum <- summary(cooc.mat)
  temp.res <- temp.sum[7]
  klb7_thinning50[[i]]$cooc.mat <- cooc.mat
  klb7_thinning50[[i]]$cooc.percent <- temp.res
  cat("trial completed", i, "\n")
}


klb7_thinned.site.number <- lapply(klb7_thinning50, `[[`,"thinned_st")
x <- lapply(klb7_thinned.site.number, `[[`,1)
klb7_thinned.site.number <- matrix(sapply(x, length), ncol = 1)
rownames(klb7_thinned.site.number) <- names(x)

coh.res.all <- lapply(klb7_thinning50, `[[`,"coh.z")
coh.res.all <- do.call(rbind, coh.res.all)

turn.res.all <- lapply(klb7_thinning50, `[[`,"turn.z")
turn.res.all <- do.call(rbind, turn.res.all)

turn.p.all <- lapply(klb7_thinning50, `[[`,"turn.p")
turn.p.all <- do.call(rbind, turn.p.all)

bc.i.all <- lapply(klb7_thinning50, `[[`,"bc.i")
bc.i.all <- do.call(rbind, bc.i.all)

bc.p.all <- lapply(klb7_thinning50, `[[`,"bc.p")
bc.p.all <- do.call(rbind, bc.p.all)

cooc.res.all <- lapply(klb7_thinning50, `[[`,"cooc.percent")
cooc.res.all <- do.call(rbind, cooc.res.all)
klb7_thinned_summary_res50 <- as.data.frame(cbind(klb7_thinned.site.number, 
                                                  coh.res.all, 
                                                  turn.res.all, 
                                                  turn.p.all,
                                                  bc.i.all, 
                                                  bc.p.all,
                                                  cooc.res.all))
colnames(klb7_thinned_summary_res50) <- c("sampled_sites", "coherence",
                                          "turnover", "turnover.p",
                                          "bc.i", "bc.p", "cooccur")





# klb7 thinning 70 --------------------------------------------------------


klb7_thinning70 <- vector(mode = "list", length = 10)
for (i in 1:10){
  names(klb7_thinning70)[i] <- paste("trial_",i)
  klb7_thinning70[[i]]$thinned_st <- thin_by_dist(data_coords, 70)
  klb7_thinning70[[i]]$thinned_dat <- data.klb7[data.klb7$station %in% 
                                                  klb7_thinning70[[i]]$thinned_st$station,]
  klb7_thinning70[[i]]$metacom_dat <- modified.remove.zeroes(klb7_thinning70[[i]]$thinned_dat)
  
  #do metacom
  klb7_thinning70[[i]]$coh.value <- 
    metacom::Coherence(klb7_thinning70[[i]]$metacom_dat
                       [,colnames(klb7_thinning70[[i]]$metacom_dat) %in% taxa.keep])
  klb7_thinning70[[i]]$coh.z <- klb7_thinning70[[i]]$coh.value[2, 2]
  klb7_thinning70[[i]]$turn.value <- 
    metacom::Turnover(klb7_thinning70[[i]]$metacom_dat
                      [,colnames(klb7_thinning70[[i]]$metacom_dat) %in% taxa.keep])
  
  klb7_thinning70[[i]]$turn.z <- klb7_thinning70[[i]]$turn.value[2,2]
  klb7_thinning70[[i]]$turn.p <- klb7_thinning70[[i]]$turn.value[3,2]
  
  klb7_thinning70[[i]]$bc.value <- 
    metacom::BoundaryClump(klb7_thinning70[[i]]$metacom_dat
                           [,colnames(klb7_thinning70[[i]]$metacom_dat) %in% taxa.keep])
  klb7_thinning70[[i]]$bc.i <- klb7_thinning70[[i]]$bc.value[1,2]
  klb7_thinning70[[i]]$bc.p <- klb7_thinning70[[i]]$bc.value[2,2]
  
  
  #do cooccur
  klb7_thinning70[[i]]$cooccur_dat <- remove.sing(klb7_thinning70[[i]]$metacom_dat)
  cooccur_dat <- klb7_thinning70[[i]]$cooccur_dat
  cooccur_dat <- cooccur_dat[,colnames(cooccur_dat) %in% taxa.keep]
  
  temp2 <- as.data.frame(t(cooccur_dat))
  klb7_thinning70[[i]]$cooccur_dat <- temp2
  cooc.mat <- cooccur::cooccur(temp2, type = "spp_site", spp_names = TRUE,
                               thresh = TRUE)
  temp.sum <- summary(cooc.mat)
  temp.res <- temp.sum[7]
  klb7_thinning70[[i]]$cooc.mat <- cooc.mat
  klb7_thinning70[[i]]$cooc.percent <- temp.res
  cat("trial completed", i, "\n")
}


klb7_thinned.site.number <- lapply(klb7_thinning70, `[[`,"thinned_st")
x <- lapply(klb7_thinned.site.number, `[[`,1)
klb7_thinned.site.number <- matrix(sapply(x, length), ncol = 1)
rownames(klb7_thinned.site.number) <- names(x)

coh.res.all <- lapply(klb7_thinning70, `[[`,"coh.z")
coh.res.all <- do.call(rbind, coh.res.all)

turn.res.all <- lapply(klb7_thinning70, `[[`,"turn.z")
turn.res.all <- do.call(rbind, turn.res.all)

turn.p.all <- lapply(klb7_thinning70, `[[`,"turn.p")
turn.p.all <- do.call(rbind, turn.p.all)

bc.i.all <- lapply(klb7_thinning70, `[[`,"bc.i")
bc.i.all <- do.call(rbind, bc.i.all)

bc.p.all <- lapply(klb7_thinning70, `[[`,"bc.p")
bc.p.all <- do.call(rbind, bc.p.all)

cooc.res.all <- lapply(klb7_thinning70, `[[`,"cooc.percent")
cooc.res.all <- do.call(rbind, cooc.res.all)
klb7_thinned_summary_res70 <- as.data.frame(cbind(klb7_thinned.site.number, 
                                                  coh.res.all, 
                                                  turn.res.all, 
                                                  turn.p.all,
                                                  bc.i.all, 
                                                  bc.p.all,
                                                  cooc.res.all))
colnames(klb7_thinned_summary_res70) <- c("sampled_sites", "coherence",
                                          "turnover", "turnover.p",
                                          "bc.i", "bc.p", "cooccur")


# NOW KLB 8 ---------------------------------------------------------------

klb8.st <- data.8[,colnames(data.8) %in% mt.data]
stations <- rownames(data.8)
stations <- as.numeric(stations)
klb8.st <- klb8.st[klb8.st$station %in% stations,]

data_coords <- st_as_sf(klb8.st, coords = c("Longitude", "Latitude"))
st_crs(data_coords) <- 4326

#klb8 thinning 30 ----

klb8_thinning30 <- vector(mode = "list", length = 10)
for (i in 1:10){
  names(klb8_thinning30)[i] <- paste("trial_",i)
  klb8_thinning30[[i]]$thinned_st <- thin_by_dist(data_coords, 30)
  klb8_thinning30[[i]]$thinned_dat <- data.klb8[data.klb8$station %in% 
                                                  klb8_thinning30[[i]]$thinned_st$station,]
  klb8_thinning30[[i]]$metacom_dat <- modified.remove.zeroes(klb8_thinning30[[i]]$thinned_dat)
  
  #do metacom
  klb8_thinning30[[i]]$coh.value <- 
    metacom::Coherence(klb8_thinning30[[i]]$metacom_dat
                       [,colnames(klb8_thinning30[[i]]$metacom_dat) %in% taxa.keep])
  klb8_thinning30[[i]]$coh.z <- klb8_thinning30[[i]]$coh.value[2, 2]
  klb8_thinning30[[i]]$turn.value <- 
    metacom::Turnover(klb8_thinning30[[i]]$metacom_dat
                      [,colnames(klb8_thinning30[[i]]$metacom_dat) %in% taxa.keep])
  
  klb8_thinning30[[i]]$turn.z <- klb8_thinning30[[i]]$turn.value[2,2]
  klb8_thinning30[[i]]$turn.p <- klb8_thinning30[[i]]$turn.value[3,2]
  
  klb8_thinning30[[i]]$bc.value <- 
    metacom::BoundaryClump(klb8_thinning30[[i]]$metacom_dat
                           [,colnames(klb8_thinning30[[i]]$metacom_dat) %in% taxa.keep])
  klb8_thinning30[[i]]$bc.i <- klb8_thinning30[[i]]$bc.value[1,2]
  klb8_thinning30[[i]]$bc.p <- klb8_thinning30[[i]]$bc.value[2,2]
  
  
  #do cooccur
  klb8_thinning30[[i]]$cooccur_dat <- remove.sing(klb8_thinning30[[i]]$metacom_dat)
  cooccur_dat <- klb8_thinning30[[i]]$cooccur_dat
  cooccur_dat <- cooccur_dat[,colnames(cooccur_dat) %in% taxa.keep]
  
  temp2 <- as.data.frame(t(cooccur_dat))
  klb8_thinning30[[i]]$cooccur_dat <- temp2
  cooc.mat <- cooccur::cooccur(temp2, type = "spp_site", spp_names = TRUE,
                               thresh = TRUE)
  temp.sum <- summary(cooc.mat)
  temp.res <- temp.sum[7]
  klb8_thinning30[[i]]$cooc.mat <- cooc.mat
  klb8_thinning30[[i]]$cooc.percent <- temp.res
  cat("trial completed", i, "\n")
}


klb8_thinned.site.number <- lapply(klb8_thinning30, `[[`,"thinned_st")
x <- lapply(klb8_thinned.site.number, `[[`,1)
klb8_thinned.site.number <- matrix(sapply(x, length), ncol = 1)
rownames(klb8_thinned.site.number) <- names(x)

coh.res.all <- lapply(klb8_thinning30, `[[`,"coh.z")
coh.res.all <- do.call(rbind, coh.res.all)

turn.res.all <- lapply(klb8_thinning30, `[[`,"turn.z")
turn.res.all <- do.call(rbind, turn.res.all)

turn.p.all <- lapply(klb8_thinning30, `[[`,"turn.p")
turn.p.all <- do.call(rbind, turn.p.all)

bc.i.all <- lapply(klb8_thinning30, `[[`,"bc.i")
bc.i.all <- do.call(rbind, bc.i.all)

bc.p.all <- lapply(klb8_thinning30, `[[`,"bc.p")
bc.p.all <- do.call(rbind, bc.p.all)

cooc.res.all <- lapply(klb8_thinning30, `[[`,"cooc.percent")
cooc.res.all <- do.call(rbind, cooc.res.all)
klb8_thinned_summary_res30 <- as.data.frame(cbind(klb8_thinned.site.number, 
                                                  coh.res.all, 
                                                  turn.res.all, 
                                                  turn.p.all,
                                                  bc.i.all, 
                                                  bc.p.all,
                                                  cooc.res.all))
colnames(klb8_thinned_summary_res30) <- c("sampled_sites", "coherence",
                                          "turnover", "turnover.p",
                                          "bc.i", "bc.p", "cooccur")

# klb8 thinning 50 --------------------------------------------------------


klb8_thinning50 <- vector(mode = "list", length = 10)
for (i in 1:10){
  names(klb8_thinning50)[i] <- paste("trial_",i)
  klb8_thinning50[[i]]$thinned_st <- thin_by_dist(data_coords, 50)
  klb8_thinning50[[i]]$thinned_dat <- data.klb8[data.klb8$station %in% 
                                                  klb8_thinning50[[i]]$thinned_st$station,]
  klb8_thinning50[[i]]$metacom_dat <- modified.remove.zeroes(klb8_thinning50[[i]]$thinned_dat)
  
  #do metacom
  klb8_thinning50[[i]]$coh.value <- 
    metacom::Coherence(klb8_thinning50[[i]]$metacom_dat
                       [,colnames(klb8_thinning50[[i]]$metacom_dat) %in% taxa.keep])
  klb8_thinning50[[i]]$coh.z <- klb8_thinning50[[i]]$coh.value[2, 2]
  klb8_thinning50[[i]]$turn.value <- 
    metacom::Turnover(klb8_thinning50[[i]]$metacom_dat
                      [,colnames(klb8_thinning50[[i]]$metacom_dat) %in% taxa.keep])
  
  klb8_thinning50[[i]]$turn.z <- klb8_thinning50[[i]]$turn.value[2,2]
  klb8_thinning50[[i]]$turn.p <- klb8_thinning50[[i]]$turn.value[3,2]
  
  klb8_thinning50[[i]]$bc.value <- 
    metacom::BoundaryClump(klb8_thinning50[[i]]$metacom_dat
                           [,colnames(klb8_thinning50[[i]]$metacom_dat) %in% taxa.keep])
  klb8_thinning50[[i]]$bc.i <- klb8_thinning50[[i]]$bc.value[1,2]
  klb8_thinning50[[i]]$bc.p <- klb8_thinning50[[i]]$bc.value[2,2]
  
  
  #do cooccur
  klb8_thinning50[[i]]$cooccur_dat <- remove.sing(klb8_thinning50[[i]]$metacom_dat)
  cooccur_dat <- klb8_thinning50[[i]]$cooccur_dat
  cooccur_dat <- cooccur_dat[,colnames(cooccur_dat) %in% taxa.keep]
  
  temp2 <- as.data.frame(t(cooccur_dat))
  klb8_thinning50[[i]]$cooccur_dat <- temp2
  cooc.mat <- cooccur::cooccur(temp2, type = "spp_site", spp_names = TRUE,
                               thresh = TRUE)
  temp.sum <- summary(cooc.mat)
  temp.res <- temp.sum[7]
  klb8_thinning50[[i]]$cooc.mat <- cooc.mat
  klb8_thinning50[[i]]$cooc.percent <- temp.res
  cat("trial completed", i, "\n")
}


klb8_thinned.site.number <- lapply(klb8_thinning50, `[[`,"thinned_st")
x <- lapply(klb8_thinned.site.number, `[[`,1)
klb8_thinned.site.number <- matrix(sapply(x, length), ncol = 1)
rownames(klb8_thinned.site.number) <- names(x)

coh.res.all <- lapply(klb8_thinning50, `[[`,"coh.z")
coh.res.all <- do.call(rbind, coh.res.all)

turn.res.all <- lapply(klb8_thinning50, `[[`,"turn.z")
turn.res.all <- do.call(rbind, turn.res.all)

turn.p.all <- lapply(klb8_thinning50, `[[`,"turn.p")
turn.p.all <- do.call(rbind, turn.p.all)

bc.i.all <- lapply(klb8_thinning50, `[[`,"bc.i")
bc.i.all <- do.call(rbind, bc.i.all)

bc.p.all <- lapply(klb8_thinning50, `[[`,"bc.p")
bc.p.all <- do.call(rbind, bc.p.all)

cooc.res.all <- lapply(klb8_thinning50, `[[`,"cooc.percent")
cooc.res.all <- do.call(rbind, cooc.res.all)
klb8_thinned_summary_res50 <- as.data.frame(cbind(klb8_thinned.site.number, 
                                                  coh.res.all, 
                                                  turn.res.all, 
                                                  turn.p.all,
                                                  bc.i.all, 
                                                  bc.p.all,
                                                  cooc.res.all))
colnames(klb8_thinned_summary_res50) <- c("sampled_sites", "coherence",
                                          "turnover", "turnover.p",
                                          "bc.i", "bc.p", "cooccur")





# klb8 thinning 70 --------------------------------------------------------


klb8_thinning70 <- vector(mode = "list", length = 10)
for (i in 1:10){
  names(klb8_thinning70)[i] <- paste("trial_",i)
  klb8_thinning70[[i]]$thinned_st <- thin_by_dist(data_coords, 70)
  klb8_thinning70[[i]]$thinned_dat <- data.klb8[data.klb8$station %in% 
                                                  klb8_thinning70[[i]]$thinned_st$station,]
  klb8_thinning70[[i]]$metacom_dat <- modified.remove.zeroes(klb8_thinning70[[i]]$thinned_dat)
  
  #do metacom
  klb8_thinning70[[i]]$coh.value <- 
    metacom::Coherence(klb8_thinning70[[i]]$metacom_dat
                       [,colnames(klb8_thinning70[[i]]$metacom_dat) %in% taxa.keep])
  klb8_thinning70[[i]]$coh.z <- klb8_thinning70[[i]]$coh.value[2, 2]
  klb8_thinning70[[i]]$turn.value <- 
    metacom::Turnover(klb8_thinning70[[i]]$metacom_dat
                      [,colnames(klb8_thinning70[[i]]$metacom_dat) %in% taxa.keep])
  
  klb8_thinning70[[i]]$turn.z <- klb8_thinning70[[i]]$turn.value[2,2]
  klb8_thinning70[[i]]$turn.p <- klb8_thinning70[[i]]$turn.value[3,2]
  
  klb8_thinning70[[i]]$bc.value <- 
    metacom::BoundaryClump(klb8_thinning70[[i]]$metacom_dat
                           [,colnames(klb8_thinning70[[i]]$metacom_dat) %in% taxa.keep])
  klb8_thinning70[[i]]$bc.i <- klb8_thinning70[[i]]$bc.value[1,2]
  klb8_thinning70[[i]]$bc.p <- klb8_thinning70[[i]]$bc.value[2,2]
  
  
  #do cooccur
  klb8_thinning70[[i]]$cooccur_dat <- remove.sing(klb8_thinning70[[i]]$metacom_dat)
  cooccur_dat <- klb8_thinning70[[i]]$cooccur_dat
  cooccur_dat <- cooccur_dat[,colnames(cooccur_dat) %in% taxa.keep]
  
  temp2 <- as.data.frame(t(cooccur_dat))
  klb8_thinning70[[i]]$cooccur_dat <- temp2
  cooc.mat <- cooccur::cooccur(temp2, type = "spp_site", spp_names = TRUE,
                               thresh = TRUE)
  temp.sum <- summary(cooc.mat)
  temp.res <- temp.sum[7]
  klb8_thinning70[[i]]$cooc.mat <- cooc.mat
  klb8_thinning70[[i]]$cooc.percent <- temp.res
  cat("trial completed", i, "\n")
}


klb8_thinned.site.number <- lapply(klb8_thinning70, `[[`,"thinned_st")
x <- lapply(klb8_thinned.site.number, `[[`,1)
klb8_thinned.site.number <- matrix(sapply(x, length), ncol = 1)
rownames(klb8_thinned.site.number) <- names(x)

coh.res.all <- lapply(klb8_thinning70, `[[`,"coh.z")
coh.res.all <- do.call(rbind, coh.res.all)

turn.res.all <- lapply(klb8_thinning70, `[[`,"turn.z")
turn.res.all <- do.call(rbind, turn.res.all)

turn.p.all <- lapply(klb8_thinning70, `[[`,"turn.p")
turn.p.all <- do.call(rbind, turn.p.all)

bc.i.all <- lapply(klb8_thinning70, `[[`,"bc.i")
bc.i.all <- do.call(rbind, bc.i.all)

bc.p.all <- lapply(klb8_thinning70, `[[`,"bc.p")
bc.p.all <- do.call(rbind, bc.p.all)

cooc.res.all <- lapply(klb8_thinning70, `[[`,"cooc.percent")
cooc.res.all <- do.call(rbind, cooc.res.all)
klb8_thinned_summary_res70 <- as.data.frame(cbind(klb8_thinned.site.number, 
                                                  coh.res.all, 
                                                  turn.res.all, 
                                                  turn.p.all,
                                                  bc.i.all, 
                                                  bc.p.all,
                                                  cooc.res.all))
colnames(klb8_thinned_summary_res70) <- c("sampled_sites", "coherence",
                                          "turnover", "turnover.p",
                                          "bc.i", "bc.p", "cooccur")

# KLB 9  ------------------------------------------------------------------


klb9.st <- data.9[,colnames(data.9) %in% mt.data]
stations <- rownames(data.9)
stations <- as.numeric(stations)
klb9.st <- klb9.st[klb9.st$station %in% stations,]

data_coords <- st_as_sf(klb9.st, coords = c("Longitude", "Latitude"))
st_crs(data_coords) <- 4326

#klb9 thinning 30 ----

klb9_thinning30 <- vector(mode = "list", length = 10)
for (i in 1:10){
  names(klb9_thinning30)[i] <- paste("trial_",i)
  klb9_thinning30[[i]]$thinned_st <- thin_by_dist(data_coords, 30)
  klb9_thinning30[[i]]$thinned_dat <- data.klb9[data.klb9$station %in% 
                                                  klb9_thinning30[[i]]$thinned_st$station,]
  klb9_thinning30[[i]]$metacom_dat <- modified.remove.zeroes(klb9_thinning30[[i]]$thinned_dat)
  
  #do metacom
  klb9_thinning30[[i]]$coh.value <- 
    metacom::Coherence(klb9_thinning30[[i]]$metacom_dat
                       [,colnames(klb9_thinning30[[i]]$metacom_dat) %in% taxa.keep])
  klb9_thinning30[[i]]$coh.z <- klb9_thinning30[[i]]$coh.value[2, 2]
  klb9_thinning30[[i]]$turn.value <- 
    metacom::Turnover(klb9_thinning30[[i]]$metacom_dat
                      [,colnames(klb9_thinning30[[i]]$metacom_dat) %in% taxa.keep])
  
  klb9_thinning30[[i]]$turn.z <- klb9_thinning30[[i]]$turn.value[2,2]
  klb9_thinning30[[i]]$turn.p <- klb9_thinning30[[i]]$turn.value[3,2]
  
  klb9_thinning30[[i]]$bc.value <- 
    metacom::BoundaryClump(klb9_thinning30[[i]]$metacom_dat
                           [,colnames(klb9_thinning30[[i]]$metacom_dat) %in% taxa.keep])
  klb9_thinning30[[i]]$bc.i <- klb9_thinning30[[i]]$bc.value[1,2]
  klb9_thinning30[[i]]$bc.p <- klb9_thinning30[[i]]$bc.value[2,2]
  
  
  #do cooccur
  klb9_thinning30[[i]]$cooccur_dat <- remove.sing(klb9_thinning30[[i]]$metacom_dat)
  cooccur_dat <- klb9_thinning30[[i]]$cooccur_dat
  cooccur_dat <- cooccur_dat[,colnames(cooccur_dat) %in% taxa.keep]
  
  temp2 <- as.data.frame(t(cooccur_dat))
  klb9_thinning30[[i]]$cooccur_dat <- temp2
  cooc.mat <- cooccur::cooccur(temp2, type = "spp_site", spp_names = TRUE,
                               thresh = TRUE)
  temp.sum <- summary(cooc.mat)
  temp.res <- temp.sum[7]
  klb9_thinning30[[i]]$cooc.mat <- cooc.mat
  klb9_thinning30[[i]]$cooc.percent <- temp.res
  cat("trial completed", i, "\n")
}


klb9_thinned.site.number <- lapply(klb9_thinning30, `[[`,"thinned_st")
x <- lapply(klb9_thinned.site.number, `[[`,1)
klb9_thinned.site.number <- matrix(sapply(x, length), ncol = 1)
rownames(klb9_thinned.site.number) <- names(x)

coh.res.all <- lapply(klb9_thinning30, `[[`,"coh.z")
coh.res.all <- do.call(rbind, coh.res.all)

turn.res.all <- lapply(klb9_thinning30, `[[`,"turn.z")
turn.res.all <- do.call(rbind, turn.res.all)

turn.p.all <- lapply(klb9_thinning30, `[[`,"turn.p")
turn.p.all <- do.call(rbind, turn.p.all)

bc.i.all <- lapply(klb9_thinning30, `[[`,"bc.i")
bc.i.all <- do.call(rbind, bc.i.all)

bc.p.all <- lapply(klb9_thinning30, `[[`,"bc.p")
bc.p.all <- do.call(rbind, bc.p.all)

cooc.res.all <- lapply(klb9_thinning30, `[[`,"cooc.percent")
cooc.res.all <- do.call(rbind, cooc.res.all)
klb9_thinned_summary_res30 <- as.data.frame(cbind(klb9_thinned.site.number, 
                                                  coh.res.all, 
                                                  turn.res.all, 
                                                  turn.p.all,
                                                  bc.i.all, 
                                                  bc.p.all,
                                                  cooc.res.all))
colnames(klb9_thinned_summary_res30) <- c("sampled_sites", "coherence",
                                          "turnover", "turnover.p",
                                          "bc.i", "bc.p", "cooccur")

# klb9 thinning 50 --------------------------------------------------------


klb9_thinning50 <- vector(mode = "list", length = 10)
for (i in 1:10){
  names(klb9_thinning50)[i] <- paste("trial_",i)
  klb9_thinning50[[i]]$thinned_st <- thin_by_dist(data_coords, 50)
  klb9_thinning50[[i]]$thinned_dat <- data.klb9[data.klb9$station %in% 
                                                  klb9_thinning50[[i]]$thinned_st$station,]
  klb9_thinning50[[i]]$metacom_dat <- modified.remove.zeroes(klb9_thinning50[[i]]$thinned_dat)
  
  #do metacom
  klb9_thinning50[[i]]$coh.value <- 
    metacom::Coherence(klb9_thinning50[[i]]$metacom_dat
                       [,colnames(klb9_thinning50[[i]]$metacom_dat) %in% taxa.keep])
  klb9_thinning50[[i]]$coh.z <- klb9_thinning50[[i]]$coh.value[2, 2]
  klb9_thinning50[[i]]$turn.value <- 
    metacom::Turnover(klb9_thinning50[[i]]$metacom_dat
                      [,colnames(klb9_thinning50[[i]]$metacom_dat) %in% taxa.keep])
  
  klb9_thinning50[[i]]$turn.z <- klb9_thinning50[[i]]$turn.value[2,2]
  klb9_thinning50[[i]]$turn.p <- klb9_thinning50[[i]]$turn.value[3,2]
  
  klb9_thinning50[[i]]$bc.value <- 
    metacom::BoundaryClump(klb9_thinning50[[i]]$metacom_dat
                           [,colnames(klb9_thinning50[[i]]$metacom_dat) %in% taxa.keep])
  klb9_thinning50[[i]]$bc.i <- klb9_thinning50[[i]]$bc.value[1,2]
  klb9_thinning50[[i]]$bc.p <- klb9_thinning50[[i]]$bc.value[2,2]
  
  
  #do cooccur
  klb9_thinning50[[i]]$cooccur_dat <- remove.sing(klb9_thinning50[[i]]$metacom_dat)
  cooccur_dat <- klb9_thinning50[[i]]$cooccur_dat
  cooccur_dat <- cooccur_dat[,colnames(cooccur_dat) %in% taxa.keep]
  
  temp2 <- as.data.frame(t(cooccur_dat))
  klb9_thinning50[[i]]$cooccur_dat <- temp2
  cooc.mat <- cooccur::cooccur(temp2, type = "spp_site", spp_names = TRUE,
                               thresh = TRUE)
  temp.sum <- summary(cooc.mat)
  temp.res <- temp.sum[7]
  klb9_thinning50[[i]]$cooc.mat <- cooc.mat
  klb9_thinning50[[i]]$cooc.percent <- temp.res
  cat("trial completed", i, "\n")
}


klb9_thinned.site.number <- lapply(klb9_thinning50, `[[`,"thinned_st")
x <- lapply(klb9_thinned.site.number, `[[`,1)
klb9_thinned.site.number <- matrix(sapply(x, length), ncol = 1)
rownames(klb9_thinned.site.number) <- names(x)

coh.res.all <- lapply(klb9_thinning50, `[[`,"coh.z")
coh.res.all <- do.call(rbind, coh.res.all)

turn.res.all <- lapply(klb9_thinning50, `[[`,"turn.z")
turn.res.all <- do.call(rbind, turn.res.all)

turn.p.all <- lapply(klb9_thinning50, `[[`,"turn.p")
turn.p.all <- do.call(rbind, turn.p.all)

bc.i.all <- lapply(klb9_thinning50, `[[`,"bc.i")
bc.i.all <- do.call(rbind, bc.i.all)

bc.p.all <- lapply(klb9_thinning50, `[[`,"bc.p")
bc.p.all <- do.call(rbind, bc.p.all)

cooc.res.all <- lapply(klb9_thinning50, `[[`,"cooc.percent")
cooc.res.all <- do.call(rbind, cooc.res.all)
klb9_thinned_summary_res50 <- as.data.frame(cbind(klb9_thinned.site.number, 
                                                  coh.res.all, 
                                                  turn.res.all, 
                                                  turn.p.all,
                                                  bc.i.all, 
                                                  bc.p.all,
                                                  cooc.res.all))
colnames(klb9_thinned_summary_res50) <- c("sampled_sites", "coherence",
                                          "turnover", "turnover.p",
                                          "bc.i", "bc.p", "cooccur")





# klb9 thinning 70 --------------------------------------------------------


klb9_thinning70 <- vector(mode = "list", length = 10)
for (i in 1:10){
  names(klb9_thinning70)[i] <- paste("trial_",i)
  klb9_thinning70[[i]]$thinned_st <- thin_by_dist(data_coords, 70)
  klb9_thinning70[[i]]$thinned_dat <- data.klb9[data.klb9$station %in% 
                                                  klb9_thinning70[[i]]$thinned_st$station,]
  klb9_thinning70[[i]]$metacom_dat <- modified.remove.zeroes(klb9_thinning70[[i]]$thinned_dat)
  
  #do metacom
  klb9_thinning70[[i]]$coh.value <- 
    metacom::Coherence(klb9_thinning70[[i]]$metacom_dat
                       [,colnames(klb9_thinning70[[i]]$metacom_dat) %in% taxa.keep])
  klb9_thinning70[[i]]$coh.z <- klb9_thinning70[[i]]$coh.value[2, 2]
  klb9_thinning70[[i]]$turn.value <- 
    metacom::Turnover(klb9_thinning70[[i]]$metacom_dat
                      [,colnames(klb9_thinning70[[i]]$metacom_dat) %in% taxa.keep])
  
  klb9_thinning70[[i]]$turn.z <- klb9_thinning70[[i]]$turn.value[2,2]
  klb9_thinning70[[i]]$turn.p <- klb9_thinning70[[i]]$turn.value[3,2]
  
  klb9_thinning70[[i]]$bc.value <- 
    metacom::BoundaryClump(klb9_thinning70[[i]]$metacom_dat
                           [,colnames(klb9_thinning70[[i]]$metacom_dat) %in% taxa.keep])
  klb9_thinning70[[i]]$bc.i <- klb9_thinning70[[i]]$bc.value[1,2]
  klb9_thinning70[[i]]$bc.p <- klb9_thinning70[[i]]$bc.value[2,2]
  
  
  #do cooccur
  klb9_thinning70[[i]]$cooccur_dat <- remove.sing(klb9_thinning70[[i]]$metacom_dat)
  cooccur_dat <- klb9_thinning70[[i]]$cooccur_dat
  cooccur_dat <- cooccur_dat[,colnames(cooccur_dat) %in% taxa.keep]
  
  temp2 <- as.data.frame(t(cooccur_dat))
  klb9_thinning70[[i]]$cooccur_dat <- temp2
  cooc.mat <- cooccur::cooccur(temp2, type = "spp_site", spp_names = TRUE,
                               thresh = TRUE)
  temp.sum <- summary(cooc.mat)
  temp.res <- temp.sum[7]
  klb9_thinning70[[i]]$cooc.mat <- cooc.mat
  klb9_thinning70[[i]]$cooc.percent <- temp.res
  cat("trial completed", i, "\n")
}


klb9_thinned.site.number <- lapply(klb9_thinning70, `[[`,"thinned_st")
x <- lapply(klb9_thinned.site.number, `[[`,1)
klb9_thinned.site.number <- matrix(sapply(x, length), ncol = 1)
rownames(klb9_thinned.site.number) <- names(x)

coh.res.all <- lapply(klb9_thinning70, `[[`,"coh.z")
coh.res.all <- do.call(rbind, coh.res.all)

turn.res.all <- lapply(klb9_thinning70, `[[`,"turn.z")
turn.res.all <- do.call(rbind, turn.res.all)

turn.p.all <- lapply(klb9_thinning70, `[[`,"turn.p")
turn.p.all <- do.call(rbind, turn.p.all)

bc.i.all <- lapply(klb9_thinning70, `[[`,"bc.i")
bc.i.all <- do.call(rbind, bc.i.all)

bc.p.all <- lapply(klb9_thinning70, `[[`,"bc.p")
bc.p.all <- do.call(rbind, bc.p.all)

cooc.res.all <- lapply(klb9_thinning70, `[[`,"cooc.percent")
cooc.res.all <- do.call(rbind, cooc.res.all)
klb9_thinned_summary_res70 <- as.data.frame(cbind(klb9_thinned.site.number, 
                                                  coh.res.all, 
                                                  turn.res.all, 
                                                  turn.p.all,
                                                  bc.i.all, 
                                                  bc.p.all,
                                                  cooc.res.all))
colnames(klb9_thinned_summary_res70) <- c("sampled_sites", "coherence",
                                          "turnover", "turnover.p",
                                          "bc.i", "bc.p", "cooccur")

# save all results --------------------------------------------------------

save(klb7_thinned_summary_res30, klb7_thinned_summary_res50, 
     klb7_thinned_summary_res70, klb7_thinning30, klb7_thinning70,
     klb7_thinning50, 
     klb8_thinned_summary_res30, klb8_thinned_summary_res50, 
     klb8_thinned_summary_res70, klb8_thinning30, klb8_thinning70,
     klb8_thinning50, 
     klb9_thinned_summary_res30, klb9_thinned_summary_res50, 
     klb9_thinned_summary_res70, klb9_thinning30, klb9_thinning70,
     klb9_thinning50,
     file = "results/all_thinning_res_Jul24.RData")


library(colorspace)

#--- species accum ------- goes to end ---------

par(mfrow = c(3,3))

plot(sp.7, main = "Species accumulation KLB 7, thinning by 30m", 
     xlab = "Number of sites", 
     ylab = "observed richness (number of families)", 
     xlim = c(0, 150),
     ylim = c(0,45), lwd = 2)
grid(col = "lightgrey")
for (i in 1:10){
  temp <- klb7_thinning30[[i]]$thinned_dat
  temp.sp <- specaccum(temp, method = "rarefaction")
  plot(temp.sp, add = TRUE, 
       col = adjustcolor("cornflowerblue", alpha.f = 0.2), 
       lwd = 1.5)
}

plot(sp.7, main = "Species accumulation KLB 7, thinning by 50m", 
     xlab = "Number of sites", 
     ylab = "observed richness (number of families)", xlim = c(0, 150),
     ylim = c(0,45), lwd = 2)
grid(col = "lightgrey")
for (i in 1:10){
  temp <- klb7_thinning50[[i]]$thinned_dat
  temp.sp <- specaccum(temp, method = "rarefaction")
  plot(temp.sp, add = TRUE, 
       col = adjustcolor("cornflowerblue", alpha.f = 0.2), 
       lwd = 1.5)
}

plot(sp.7, main = "Species accumulation KLB 7, thinning by 70m", 
     xlab = "Number of sites", 
     ylab = "observed richness (number of families)", xlim = c(0, 150),
     ylim = c(0,45), lwd = 2)
grid(col = "lightgrey")
for (i in 1:10){
  temp <- klb7_thinning70[[i]]$thinned_dat
  temp.sp <- specaccum(temp, method = "rarefaction")
  plot(temp.sp, add = TRUE,
       col = adjustcolor("cornflowerblue", alpha.f = 0.2),
       lwd = 1.5)
}


plot(sp.8, main = "Species accumulation KLB 8, thinning by 30m", 
     xlab = "Number of sites", 
     ylab = "observed richness (number of families)", xlim = c(0, 150),
     ylim = c(0,45), lwd = 2)
grid(col = "lightgrey")
for (i in 1:10){
  temp <- klb8_thinning30[[i]]$thinned_dat
  temp.sp <- specaccum(temp, method = "rarefaction")
  plot(temp.sp, add = TRUE, lwd = 1.5, 
       col = adjustcolor("forestgreen", alpha.f = 0.2))
}

plot(sp.8, main = "Species accumulation KLB 8, thinning by 50m", 
     xlab = "Number of sites", 
     ylab = "observed richness (number of families)", xlim = c(0, 150),
     ylim = c(0,45), lwd = 2)
grid(col = "lightgrey")
for (i in 1:10){
  temp <- klb8_thinning50[[i]]$thinned_dat
  temp.sp <- specaccum(temp, method = "rarefaction")
  plot(temp.sp, add = TRUE, lwd = 1.5, 
       col = adjustcolor("forestgreen", alpha.f = 0.2))
}

plot(sp.8, main = "Species accumulation KLB 8, thinning by 70m", 
     xlab = "Number of sites", 
     ylab = "observed richness (number of families)", xlim = c(0, 150),
     ylim = c(0,45), lwd = 2)
grid(col = "lightgrey")
for (i in 1:10){
  temp <- klb8_thinning70[[i]]$thinned_dat
  temp.sp <- specaccum(temp, method = "rarefaction")
  plot(temp.sp, add = TRUE, lwd = 1.5, 
       col = adjustcolor("forestgreen", alpha.f = 0.2))
}


plot(sp.9, main = "Species accumulation KLB 9, thinning by 30m", 
     xlab = "Number of sites", 
     ylab = "observed richness (number of families)", xlim = c(0, 150),
     ylim = c(0,45), lwd = 2)
grid(col = "lightgrey")
for (i in 1:10){
  temp <- klb9_thinning30[[i]]$thinned_dat
  temp.sp <- specaccum(temp, method = "rarefaction")
  plot(temp.sp, add = TRUE, lwd = 1.5, 
       col = adjustcolor("darkorange", alpha.f = 0.2))
}

plot(sp.9, main = "Species accumulation KLB 9, thinning by 50m", 
     xlab = "Number of sites", 
     ylab = "observed richness (number of families)", xlim = c(0, 150),
     ylim = c(0,45), lwd = 2)
grid(col = "lightgrey")
for (i in 1:10){
  temp <- klb9_thinning50[[i]]$thinned_dat
  temp.sp <- specaccum(temp, method = "rarefaction")
  plot(temp.sp, add = TRUE, lwd = 1.5, 
       col = adjustcolor("darkorange", alpha.f = 0.2))
}

plot(sp.9, main = "Species accumulation KLB 9, thinning by 70m", 
     xlab = "Number of sites", 
     ylab = "observed richness (number of families)", xlim = c(0, 150),
     ylim = c(0,45), lwd = 2)
grid(col = "lightgrey")
for (i in 1:10){
  temp <- klb9_thinning70[[i]]$thinned_dat
  temp.sp <- specaccum(temp, method = "rarefaction")
  plot(temp.sp, add = TRUE, lwd = 1.5, 
       col = adjustcolor("darkorange", alpha.f = 0.2))
}


# thinning coocur boxplots ------------------------------------------------

#load(file = "results/all_thinning_res_Jul24.RData")
load(file = "results/cooccur.sig.res.RData")


#klb7
klb7_thin30.cooc <- as.data.frame(klb7_thinned_summary_res30$cooccur)
klb7_thin30.cooc$klb <- rep("klb7_thin30", 10)
klb7_thin30.cooc$treatment <- rep("trials", 10)
colnames(klb7_thin30.cooc) <- c("sig.res", "klb", "treatment")

klb7_thin50.cooc <- as.data.frame(klb7_thinned_summary_res50$cooccur)
klb7_thin50.cooc$klb <- rep("klb7_thin50", 10)
klb7_thin50.cooc$treatment <- rep("trials", 10)
colnames(klb7_thin50.cooc) <- c("sig.res", "klb", "treatment")

klb7_thin70.cooc <- as.data.frame(klb7_thinned_summary_res70$cooccur)
klb7_thin70.cooc$klb <- rep("klb7_thin70", 10)
klb7_thin70.cooc$treatment <- rep("trials", 10)
colnames(klb7_thin70.cooc) <- c("sig.res", "klb", "treatment")

#klb8
klb8_thin30.cooc <- as.data.frame(klb8_thinned_summary_res30$cooccur)
klb8_thin30.cooc$klb <- rep("klb8_thin30", 10)
klb8_thin30.cooc$treatment <- rep("trials", 10)
colnames(klb8_thin30.cooc) <- c("sig.res", "klb", "treatment")

klb8_thin50.cooc <- as.data.frame(klb8_thinned_summary_res50$cooccur)
klb8_thin50.cooc$klb <- rep("klb8_thin50", 10)
klb8_thin50.cooc$treatment <- rep("trials", 10)
colnames(klb8_thin50.cooc) <- c("sig.res", "klb", "treatment")

klb8_thin70.cooc <- as.data.frame(klb8_thinned_summary_res70$cooccur)
klb8_thin70.cooc$klb <- rep("klb8_thin70", 10)
klb8_thin70.cooc$treatment <- rep("trials", 10)
colnames(klb8_thin70.cooc) <- c("sig.res", "klb", "treatment")



#klb9
thin30.cooc <- as.data.frame(klb9_thinned_summary_res30$cooccur)
thin30.cooc$klb <- rep("klb9_thin30", 10)
thin30.cooc$treatment <- rep("trials", 10)
colnames(thin30.cooc) <- c("sig.res", "klb", "treatment")

thin50.cooc <- as.data.frame(klb9_thinned_summary_res50$cooccur)
thin50.cooc$klb <- rep("klb9_thin50", 10)
thin50.cooc$treatment <- rep("trials", 10)
colnames(thin50.cooc) <- c("sig.res", "klb", "treatment")

thin70.cooc <- as.data.frame(klb9_thinned_summary_res70$cooccur)
thin70.cooc$klb <- rep("klb9_thin70", 10)
thin70.cooc$treatment <- rep("trials", 10)
colnames(thin70.cooc) <- c("sig.res", "klb", "treatment")


cooccur.sig.res <- rbind(cooccur.sig.res[1:5,], klb7_thin30.cooc, klb7_thin50.cooc,
                         klb7_thin70.cooc, 
                         cooccur.sig.res[6:10,],
                         klb8_thin30.cooc, klb8_thin50.cooc,
                         klb8_thin70.cooc, 
                         cooccur.sig.res[11:13,],
                         thin30.cooc, thin50.cooc,
                         thin70.cooc)
rownames(cooccur.sig.res) <- NULL
cooccur.sig.res$klb <- 
  factor(cooccur.sig.res$klb, 
      levels = c("klb7", "klb7_thin30", 
        "klb7_thin50", "klb7_thin70", 
          "klb8", "klb8_thin30", 
            "klb8_thin50", "klb8_thin70",
                  "klb9", "klb9_thin30", 
                      "klb9_thin50", "klb9_thin70"
      ))

#plot the data

library(ggplot2)

#colors -------


coocc.box <- ggplot(cooccur.sig.res, aes(x = klb, y = sig.res, fill = klb)) +
  geom_boxplot(width = 0.5, position = position_dodge(0.6)) + 
  geom_point(
    position = position_jitterdodge(jitter.width = 0.02),
    aes(shape = treatment),
    size = 5  # control the size of the shape
  ) + 
  #theme(legend.position = "none") + 
  scale_fill_manual(values =  c("darkblue", "royalblue", "cornflowerblue",  
                                "lightblue", 
                                "darkgreen", "forestgreen", "olivedrab","yellowgreen",
                                "darkorange", "orange", "salmon", "pink")) +
  scale_shape_manual(
    values = c("trials" = 8, "raw" = 1, "contract" = 0, 
               "left.shift" = 2, "right.shift" = 6) # 4 = cross; try 8 for star
  ) +
  ggtitle("Cooccur results") +
  guides(fill = "none", shape = "none", color = "none", linetype = "none") + 

  theme_bw() + ylim(0, 25)
  

coocc.box


#draw the metacom plots -----------------

#import klb 7 and 8 and 9 metacom results ----

load(file = "results/klb7.metacom.results.RData")
load(file = "results/klb8.metacom.results.RData")
load(file = "results/klb9.metacom.results.RData")

#plot metacom Res 
coh.res.7 <- (klb7.metacom.res.df[2, c(2, 8, 14, 20)])
turn.res.7 <- (klb7.metacom.res.df[2, c(4, 10, 16, 22)])
bc.res.7 <- (klb7.metacom.res.df[1, c(6, 12, 18, 24)])

coh.res.8 <- (klb8.metacom.res.df[2, c(2, 8, 14, 20)])
turn.res.8 <- (klb8.metacom.res.df[2, c(4, 10, 16, 22)])
bc.res.8 <- (klb8.metacom.res.df[1, c(6, 12, 18, 24)])

coh.res.9 <- (klb9.metacom.res.df[2, c(2, 8)])
turn.res.9 <- (klb9.metacom.res.df[2, c(4, 10)])
bc.res.9 <- (klb9.metacom.res.df[1, c(6, 12)])

plot(NA, xlim = c(-5, 5), ylim = c(-11, -1),
     xlab = "Turnover Z-score", ylab = "Coherence Z-score",
     main = "All Metacommunity Results", xaxt = "n", yaxt = "n")
axis(1, at = seq(-5, 5, by = 1))     
axis(2, at = seq(-11, -1, by = 1))   
abline(v = 0, lwd = 2)
box()

#klb9


#klb9 thinned points 

points(klb9_thinned_summary_res30$turnover, klb9_thinned_summary_res30$coherence, 
       pch = 16, cex = klb9_thinned_summary_res30$bc.i, 
       col = adjustcolor("darkorange", alpha.f = 0.8) )

points(klb9_thinned_summary_res50$turnover, klb9_thinned_summary_res50$coherence, 
       pch = 16, cex = klb9_thinned_summary_res50$bc.i, 
       col = adjustcolor("salmon", alpha.f = 0.5) )

points(klb9_thinned_summary_res70$turnover, klb9_thinned_summary_res70$coherence, 
       pch = 16, cex = klb9_thinned_summary_res50$bc.i, 
       col = adjustcolor("pink", alpha.f = 0.5) )

points(turn.res.9$raw.turnover_res, coh.res.9$raw.coherence_res, 
       pch = 22, cex = bc.res.9$raw.boundary_clumping_res, 
       bg = "darkorange", lwd = 1, col = "black")



#klb8


#thinned points
points(klb8_thinned_summary_res30$turnover, klb8_thinned_summary_res30$coherence, 
       pch = 16, cex = klb8_thinned_summary_res30$bc.i, 
       col = adjustcolor("forestgreen", alpha.f = 0.8) )

points(klb8_thinned_summary_res50$turnover, klb8_thinned_summary_res50$coherence, 
       pch = 16, cex = klb8_thinned_summary_res50$bc.i, 
       col = adjustcolor("olivedrab", alpha.f = 0.5) )

points(klb8_thinned_summary_res70$turnover, klb8_thinned_summary_res70$coherence, 
       pch = 16, cex = klb8_thinned_summary_res50$bc.i, 
       col = adjustcolor("yellowgreen", alpha.f = 0.5) )

points(turn.res.8$raw.turnover_res, coh.res.8$raw.coherence_res, 
       pch = 22, cex = bc.res.8$raw.boundary_clumping_res, 
       bg = "forestgreen", lwd = 1, col = "black")



#klb 7


points(klb7_thinned_summary_res30$turnover, klb7_thinned_summary_res30$coherence, 
       pch = 16, cex = klb7_thinned_summary_res30$bc.i, 
       col = adjustcolor("royalblue", alpha.f = 0.8) )

points(klb7_thinned_summary_res50$turnover, klb7_thinned_summary_res50$coherence, 
       pch = 16, cex = klb7_thinned_summary_res50$bc.i, 
       col = adjustcolor("cornflowerblue", alpha.f = 0.5) )

points(klb7_thinned_summary_res70$turnover, klb7_thinned_summary_res70$coherence, 
       pch = 16, cex = klb7_thinned_summary_res50$bc.i, 
       col = adjustcolor("lightblue", alpha.f = 0.5) )

points(turn.res.7$raw.turnover_res, coh.res.7$raw.coherence_res, 
       pch = 22, cex = bc.res.7$raw.boundary_clumping_res, 
       bg = "darkblue", lwd = 1, col = "black")










