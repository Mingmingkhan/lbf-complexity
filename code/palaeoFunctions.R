
move.line<-function(line1,d1)
{#change line by buffer zone given as distance d1 from the 
  #perpendicular line, line1
  
  #function to move line
  x1<-line1[1,1]
  x2<-line1[nrow(line1),1]
  y1<-line1[1,2]
  y2<-line1[nrow(line1),2]
  
  grad<-(y2-y1)/(x2-x1) #start and end points of the line
  perp.grad<- -1/grad
  
  res<-list()
  res[[1]]<-matrix(NA,nrow(line1),2)
  res[[2]]<-matrix(NA,nrow(line1),2)
  
  for(i in 1:nrow(line1))
  {
    x.coord<-line1[i,1]	
    y.coord<-line1[i,2]	
    
    res[[1]][i,1]<-x.coord+(d1/sqrt(perp.grad^2+1))
    res[[1]][i,2]<-y.coord+perp.grad*(res[[1]][i,1]-x.coord)
    
    res[[2]][i,1]<-x.coord-(d1/sqrt(perp.grad^2+1))
    res[[2]][i,2]<-y.coord+perp.grad*(res[[2]][i,1]-x.coord)
    
    #check its working ok, should be =d1
    print(sqrt((res[[1]][i,2]-line1[i,2])^2+(res[[1]][i,1]-line1[i,1])^2))
  }
  
  
  print(sqrt((res[[1]][1,2]-y1)^2+(res[[1]][1,1]-x1)^2))
  
  return(res)
}

get.new.line <- function(ori.line, moved.line){
  #distance between ori.line and moved.line 
  library(geosphere)
  hav.dist <- matrix(NA, nrow(moved.line))
  for (i in 1:nrow(ori.line)){
    hav.dist[i] <- distHaversine(c(ori.line[i,1],ori.line[i,2]),
                                 c(moved.line[i,1],moved.line[i,2]))
  }
  print(hav.dist)
  
  #shift the bottom coord
  del.y <- ori.line[nrow(ori.line),2] - moved.line[nrow(moved.line), 2]
  new.line <- matrix(NA, nrow(moved.line), 2)
  new.line[,1] <- moved.line[,1] #long
  new.line[,2] <- moved.line[,2]+del.y
  
  #test Haversine again 
  hav.dist <- matrix(NA, nrow(moved.line))
  for (i in 1:nrow(ori.line)){
    hav.dist[i] <- distHaversine(c(ori.line[i,1],ori.line[i,2]),
                                 c(new.line[i,1],new.line[i,2]))
  }
  print(hav.dist)
  
  return(new.line)
  
}

remove.zeroes <- function(df) {
  df <- as.matrix(df)
  
  # Remove columns (taxa) with all zeros
  col_sums <- colSums(df, na.rm = TRUE)
  df <- df[, col_sums > 0, drop = FALSE]
  
  # Remove rows (sites) with all zeros
  row_sums <- rowSums(df, na.rm = TRUE)
  df <- df[row_sums > 0, , drop = FALSE]
  
  # Convert to presence-absence (0/1)
  df[df >= 1] <- 1
  
  # Convert back to data frame (optional)
  df <- as.data.frame(df)
  
  return(df)
}

remove.sing.zeroes <- function(df){
  #remove taxa with 0 counts, change to presence absence
  df[df >=1] <- 1 #change to presence absence 
  #remove singletons 
  x <- as.numeric(colSums(df))
  x.ind <- which(x ==1)
  df2 <- df[,-x.ind]
  c <- as.numeric(colSums(df2))
  c.ind <- which(c == 0)
  #now remove zeroes 
  if (length(c.ind > 0)){
    temp <- df2[,-c.ind ]
  }
  
  r <- as.numeric(rowSums(temp))
  r.ind <- which(r==0)
  if (length(r.ind > 0)){
    temp2 <- temp[-r.ind,]
  }
  else temp2 <- temp
  return(temp2)
}





remove.sing <- function(df) {
  df <- as.matrix(df)
  # Get occupancy
  col_sums <- colSums(df, na.rm = TRUE)
 # Identify singleton columns (occur only once)
  keep_cols <- col_sums > 1
  # If all columns are removed, return warning
  if (sum(keep_cols) == 0) {
    warning("All columns were singletons — returning an empty data frame.")
    return(data.frame())
  }

  df_cleaned <- as.data.frame(df[, keep_cols, drop = FALSE])
  return(df_cleaned)
}

get.metacom.res <- function(df){

  #run metacom
  library(metacom)
  df.coh <- Coherence(df)
  df.turn <- Turnover(df)
  df.bc <- BoundaryClump(df)
  
  #combine results
  res <- cbind(df.coh, df.turn, df.bc)
  colnames(res) <- c("coherence", "coherence_res", "turnover", "turnover_res",
                  "boundary_clumping", "boundary_clumping_res")
  return(res)
  
}

metacom.output <- function(metacom.res){
  ## input is the list of results 
  ##get the named output based on metacom scores 
  metacom.named.output <- vector(mode = "character", length =length(metacom.res) )
  
  for (i in 1:length(metacom.res)) {
    temp <- metacom.res[[i]]
    #extract the parameters
    coh.z <- temp[2, 2]
    coh.p <- temp[3,2]
    turn.z <- temp[2, 4]
    turn.p <- temp[3, 4]
    bc <- temp[1, 6]
    bc.p <- temp[2,6]
    
    #allocate the named outputs 
    #1) Random 
    if (coh.p > 0.05) {
      x <- "Random"
      metacom.named.output[i] <- x
      print(x)}

    # 2) Checkerboard
    if ((coh.p <= 0.05) & (coh.z > 0)) {
      x <- "Checkerboard"
      metacom.named.output[i] <- x
      print(x)}
    
    #3) Clementsian 
    
    
    #3) Clementsian 
    if ((coh.z < 0) & (turn.z > 0) & (turn.p <= 0.05)
            & (bc.p <= 0.05) & (bc > 1)) {
      x <- "Clementsian"
      metacom.named.output[i] <- x
      print(x)}
    
    #4) quasi Clementsian
    if ((coh.z < 0) & (turn.z > 0) & (turn.p > 0.05)
             & (bc.p <= 0.05) & (bc > 1)) {
      x <- "quasi Clementsian"
      metacom.named.output[i] <- x
      print(x)}
    
    #5) Gleasonian 
    if ((coh.z < 0) & (turn.z > 0) & (turn.p <= 0.05) &
             (bc.p > 0.05)) {
      x <- "Gleasonian"
      metacom.named.output[i] <- x
      print(x)}
    
    #6) quasi Gleasonian 
    if ((coh.z < 0) & (turn.z > 0) & (turn.p > 0.05) &
             (bc.p > 0.05)) {
      x <- "quasi Gleasonian"
      metacom.named.output[i] <- x
      print(x)}
    
    #7) Evenly spaced  
    if ((coh.z < 0) & (turn.z > 0) & (turn.p <= 0.05)
             & (bc.p <= 0.05) & (bc < 1)) {
      x <- "Evenly spaced"
      metacom.named.output[i] <- x
      print(x)}
    
    #8) quasi Evenly spaced  
    if ((coh.z < 0) & (turn.z > 0) & (turn.p > 0.05)
             & (bc.p <= 0.05) & (bc < 1)) {
      x <- "quasi Evenly spaced"
      metacom.named.output[i] <- x
      print(x)}
    
    #9) Nested clumped 
    if ((coh.z < 0) & (turn.z < 0) & (turn.p <= 0.05)
             & (bc.p <= 0.05) & (bc > 1)) {
      x <- "Nested clumped"
      metacom.named.output[i] <- x
      print(x)}
    
    #10) quasi Nested clumped
    if ((coh.z < 0) & (turn.z < 0) & (turn.p > 0.05)
             & (bc.p <= 0.05) & (bc > 1)) {
      x <- "quasi Nested clumped"
      metacom.named.output[i] <- x
      print(x)}
    
    #11) Nested random
    if ((coh.z < 0) & (turn.z < 0) & (turn.p <= 0.05)
          & (bc.p > 0.05)) {
      x <- "Nested random"
      metacom.named.output[i] <- x
      print(x)}
    
    #12) quasi nested random 
    if ((coh.z < 0) & (turn.z < 0) & (turn.p > 0.05)
             & (bc.p > 0.05)) {
      x <- "quasi Nested random"
      metacom.named.output[i] <- x
      print(x)}
    
    #13) nested hyper dispered
    if ((coh.z < 0) & (turn.z < 0) & (turn.p <= 0.05)
             & (bc.p < 0.05) & (bc < 1))  {
      x <- "Nested hyper dispersed"
      metacom.named.output[i] <- x
      print(x)}
    
    #14) quasi nested hyper dispered
    if ((coh.z < 0) & (turn.z < 0) & (turn.p > 0.05)
             & (bc.p < 0.05) & (bc < 1))  {
      x <- "quasi Nested hyper dispersed"
      metacom.named.output[i] <- x
      print(x)}
    
  } #end of if statements
  
return(metacom.named.output)
#print(metacom.named.output)
  
}
  
plot.interact <- function(df, int.col, z){
  #windows()
  #plot the imagine function more nicely
  #takes the df of site by taxa, and the color of the treatment
  #z = where the species names plots
  df <- OrderMatrix(df, scores =1, outputScores = FALSE)
  metacom::Imagine(df, col = c(0, int.col),fill=FALSE,
          xlab = "", ylab = "", speciesnames = FALSE,
          sitenames = FALSE, order = TRUE)
  axis(3, at = c(1:ncol(df)), labels = FALSE)
  axis(side = 3, las = 2, mgp = c(3, 0.75, 0), labels = FALSE)
  text(x = 1:ncol(df), y = par("usr")[2] + z, 
       labels = colnames(df), 
       srt = -35, cex = 1, xpd = NA, adj = 1)
  axis(side = 2, at = c(1:nrow(df)), labels= FALSE)
  axis(side = 2, las = 2, mgp = c(3, 0.75, 0), labels = FALSE, 
       cex.axis = 0.8)
  text(y = 1:nrow(df), x = par("usr")[1] - 1, 
       labels = rownames(df), xpd = NA, cex =0.8)  

  
}


plot.interact.klb789 <- function(df){
  #windows()
  #plot the imagine function more nicely
  #takes the df of site by taxa, and the color of the treatment
  #z = where the species names plots
  df2 <- OrderMatrix(df, scores =1, outputScores = FALSE)
  metacom::Imagine(df2, col = c(0, int.col),fill=FALSE,
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
  
  #klb7
  klb7.st <- data.klb7$station
  ind.7 <- match(klb7.st, all.sites)
  axis(side = 2, at = ind.7[1:length(ind.7)], labels = FALSE, 
       col.ticks = "forestgreen", lwd = 2)
  
  #klb8
  klb8.st <- data.klb8$station
  ind.8 <- match(klb8.st, all.sites)
  axis(side = 2, at = ind.8[1:length(ind.8)], labels = FALSE, 
       col.ticks = "blue", lwd = 2)
  
  #klb9 
  klb9.st <- data.klb9$station
  ind.9 <- match(klb9.st, all.sites)
  axis(side = 2, at = ind.9[1:length(ind.9)], labels = FALSE, 
       col.ticks = "orange", lwd = 2)
  
}
