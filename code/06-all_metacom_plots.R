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

x.pol <- c(0, 0, 5, 5)
y.pol <- c(-11, 0, 0, -11)

plot(NA, xlim = c(-5, 5), ylim = c(-11, -1),
     xlab = "Turnover Z-score", ylab = "Coherence Z-score",
     main = "All Metacommunity Results", xaxt = "n", yaxt = "n")
axis(1, at = seq(-5, 5, by = 1))     
axis(2, at = seq(-11, -1, by = 1))   
abline(v = 0, lwd = 2)
box()

#klb9
points(turn.res.9$raw.turnover_res, coh.res.9$raw.coherence_res, 
       pch = 21, cex = bc.res.9$raw.boundary_clumping_res, 
       bg = "#c9f235ff", lwd = 2, col = "black")

#contract
points(turn.res.9$treat1.turnover_res, coh.res.9$treat1.coherence_res, 
       pch = 22, cex = bc.res.9$treat1.boundary_clumping_res,
       bg = "#c9f235ff", lwd = 2)
#text(3.5, -8.0, klb9.metacom.outputs[3])

#klb8
points(turn.res.8$raw.turnover_res, coh.res.8$raw.coherence_res, 
       pch = 21, cex = bc.res.8$raw.boundary_clumping_res, 
       bg = "white", lwd = 2, col = "#7fb62bff")

#left
points(turn.res.8$treat1.turnover_res, coh.res.8$treat1.coherence_res, 
       pch = 6, cex = bc.res.8$treat1.boundary_clumping_res,
       col = "#7fb62bff", lwd = 2)

#right
points(turn.res.8$treat2.turnover_res, coh.res.8$treat2.coherence_res, 
       pch = 2, cex = bc.res.8$treat2.boundary_clumping_res,
       col = "#7fb62bff", lwd = 2)

#contract
points(turn.res.8$treat3.turnover_res, coh.res.8$treat3.coherence_res, 
       pch = 0, cex = bc.res.8$treat3.boundary_clumping_res,
       col = "#7fb62bff", lwd = 2)

#klb7 
points(turn.res.7$raw.turnover_res, coh.res.7$raw.coherence_res, 
       pch = 21, cex = bc.res.7$raw.boundary_clumping_res, 
       bg = "#00602dff", lwd = 2, col = "black")

#left
points(turn.res.7$treat1.turnover_res, coh.res.7$treat1.coherence_res, 
       pch = 6, cex = bc.res.7$treat1.boundary_clumping_res,
       col = "#00602dff", lwd = 2)
#right
points(turn.res.7$treat2.turnover_res, coh.res.7$treat2.coherence_res, 
       pch = 2, cex = bc.res.7$treat2.boundary_clumping_res,
       col = "#00602dff", lwd = 2)

#contract
points(turn.res.7$treat3.turnover_res, coh.res.7$treat3.coherence_res, 
       pch = 22, cex = bc.res.7$treat3.boundary_clumping_res,
       bg = "#00602dff", lwd = 2)

