#Open all cooccur res
# make boxplot of significant aggregations 

load(file = "results/klb7.cooccur.results.RData")
load(file = "results/klb8.cooccur.results.RData")
load(file = "results/klb9.cooccur.res.RData")

temp <- matrix(NA, nrow = length(klb7.cooc.res))
for (i in 1:length(klb7.cooc.res)){
  temp2 <- klb7.cooc.res[[i]]["percent_sig"]
  temp[i] <- as.numeric(temp2)
}


klb7.sig.res <- matrix(NA, nrow = length(klb7.cooc.res), ncol = 3)
klb7.sig.res[,1] <- "klb7"
klb7.sig.res[,2] <- names(klb7.cooc.res)
klb7.sig.res[,3] <- temp
klb7.sig.res <- as.data.frame(klb7.sig.res)

temp <- matrix(NA, nrow = length(klb8.cooc.res))
for (i in 1:length(klb8.cooc.res)){
  temp2 <- klb8.cooc.res[[i]]["percent_sig"]
  temp[i] <- as.numeric(temp2)
}

klb8.sig.res <- matrix(NA, nrow = length(klb8.cooc.res), ncol = 3)
klb8.sig.res[,1] <- "klb8"
klb8.sig.res[,2] <- names(klb8.cooc.res)
klb8.sig.res[,3] <- temp
klb8.sig.res <- as.data.frame(klb8.sig.res)

temp <- matrix(NA, nrow = length(klb9.cooc.res))
for (i in 1:length(klb9.cooc.res)){
  temp2 <- klb9.cooc.res[[i]]["percent_sig"]
  temp[i] <- as.numeric(temp2)
}

klb9.sig.res <- matrix(NA, nrow = length(klb9.cooc.res), ncol = 3)
klb9.sig.res[,1] <- "klb9"
klb9.sig.res[,2] <- names(klb9.cooc.res)
klb9.sig.res[,3] <- temp
klb9.sig.res <- as.data.frame(klb9.sig.res)


cooccur.sig.res <- rbind(klb7.sig.res, klb8.sig.res, klb9.sig.res)
colnames(cooccur.sig.res) <- c("klb", "treatment", "sig.res")
cooccur.sig.res$klb <- as.factor(cooccur.sig.res$klb)
cooccur.sig.res$treatment <- as.factor(cooccur.sig.res$treatment)
cooccur.sig.res$sig.res <- as.numeric(cooccur.sig.res$sig.res)
save(cooccur.sig.res, file = "results/cooccur.sig.res.RData")

library(ggplot2)

coocc.box <- ggplot(cooccur.sig.res, aes(x = klb, y = sig.res, 
                                         fill = klb)) +
  geom_boxplot(width = 0.5, position = position_dodge(0.6)) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.02),
             aes(shape = treatment, size = 0.001)) + 
  theme(legend.position = "none") + 
  scale_fill_manual(values=c("#00a650ff","#a6da5aff", "#cce86eff")) +
  ggtitle("test") + 
  #coord_flip() +
  theme_bw()

coocc.box






