library(JUMP)
library(radjust)
library(STAREG)
library(adaFilter)
library(JMirror)
library(CARLIS)
library(CoCoNuT)
library(ReAD)
source("./R/methods.R")

## Primary
# SCZ2022
# PGC3_SCZ_wave3.eur 
scz_eur <- read.table("./data analysis/BD-SCZ/scz2022/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv",
                      quote = "", fill = TRUE, comment.char = "!", header = TRUE)
scz_eur$SNP_ID <- paste0(scz_eur$CHROM, ":", scz_eur$POS, ":", scz_eur$A1, ":", scz_eur$A2)
rownames(scz_eur) <- scz_eur$SNP_ID

# BIP2024
# bip2024_eur_no23andMe.gz
bip_eur <- read.table("./data analysis/BD-SCZ/bip2024/bip2024_eur_no23andMe.gz",
                      quote = "", fill = TRUE, comment.char = "!", header = TRUE)
bip_eur$SNP_ID <- paste0(bip_eur$CHR, ":", bip_eur$BP, ":", bip_eur$A1, ":", bip_eur$A2)
rownames(bip_eur) <- bip_eur$SNP_ID


## Auxiliary
# PGC3_SCZ_wave3.latino 
scz_latino <- read.table("./data analysis/BD-SCZ/scz2022/PGC3_SCZ_wave3.latino.autosome.public.v3.vcf.tsv",
                         quote = "", fill = TRUE, comment.char = "!", header = TRUE)
scz_latino$SNP_ID <- paste0(scz_latino$CHROM, ":", scz_latino$POS, ":", scz_latino$A1, ":", scz_latino$A2)
rownames(scz_latino) <- scz_latino$SNP_ID

# bip2024_eas_no23andMe.gz
bip_eas <- read.table("./data analysis/BD-SCZ/bip2024/bip2024_eas_no23andMe.gz",
                      quote = "", fill = TRUE, comment.char = "!", header = TRUE)
bip_eas$SNP_ID <- paste0(bip_eas$CHR, ":", bip_eas$BP, ":", bip_eas$A1, ":", bip_eas$A2)
rownames(bip_eas) <- bip_eas$SNP_ID

# PGC3_SCZ_wave3.afram 
scz_afram <- read.table("./data analysis/BD-SCZ/scz2022/PGC3_SCZ_wave3.afram.autosome.public.v3.vcf.tsv",
                        quote = "", fill = TRUE, comment.char = "!", header = TRUE)
scz_afram$SNP_ID <- paste0(scz_afram$CHROM, ":", scz_afram$POS, ":", scz_afram$A1, ":", scz_afram$A2)
rownames(scz_afram) <- scz_afram$SNP_ID

# intersect
snp_id <- intersect(scz_eur$SNP_ID, bip_eur$SNP_ID)
snp_id <- intersect(snp_id, scz_latino$SNP_ID)
snp_id <- intersect(snp_id, bip_eas$SNP_ID)
snp_id <- intersect(snp_id, scz_afram$SNP_ID)


pa <- scz_eur[snp_id,]$PVAL
pb <- bip_eur[snp_id,]$P
x1 <- scz_latino[snp_id,]$PVAL
x2 <- bip_eas[snp_id,]$P
x3 <- scz_afram[snp_id,]$PVAL

x12 <- Cauchy(rbind(x1, x2))
x23 <- Cauchy(rbind(x2, x3))
x13 <- Cauchy(rbind(x1, x3))
x <- Cauchy(rbind(x1, x2, x3))

alphas = c(5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1)

# STAREG
res.eb <- stareg(pa, pb)
# ReAD
res.hmm <- RepLIS(pa, pb)
# AdaFilter
res.adafilter <- adaFilter(cbind(pa, pb), r = 2, type.I.err = "FDR")
# coconut
res.coconut <- coconut(pa, pb, x)
res.coconut1 <- coconut(pa, pb, x1)
res.coconut2 <- coconut(pa, pb, x2)
res.coconut3 <- coconut(pa, pb, x3)
res.coconut12 <- coconut(pa, pb, x12)
res.coconut23 <- coconut(pa, pb, x23)
res.coconut13 <- coconut(pa, pb, x13)

# carlis
res.carlis <- CARLIS(pa, pb, x)
res.carlis1 <- CARLIS(pa, pb, x1)
res.carlis2 <- CARLIS(pa, pb, x2)
res.carlis3 <- CARLIS(pa, pb, x3)
res.carlis12 <- CARLIS(pa, pb, x12)
res.carlis23 <- CARLIS(pa, pb, x23)
res.carlis13 <- CARLIS(pa, pb, x13)

snp_jump <- list()
snp_adafilter <- list()
snp_carlis <- list()
snp_carlis1 <- list()
snp_carlis2 <- list()
snp_carlis3 <- list()
snp_carlis12 <- list()
snp_carlis23 <- list()
snp_carlis13 <- list()
snp_coconut <- list()
snp_coconut1 <- list()
snp_coconut2 <- list()
snp_coconut3 <- list()
snp_coconut12 <- list()
snp_coconut23 <- list()
snp_coconut13 <- list()
snp_read <- list()
snp_stareg <- list()

i = 0
for(q in alphas){
  i = i + 1
  # JUMP
  jump.obj <- JUMP(pa, pb, q)
  jump.thr <- jump.obj$jump.thr
  p.max <- jump.obj$p.max
  snp_jump[[i]] <- snp_id[p.max <= jump.thr]
  
  # STAREG
  padj.eb <- res.eb$fdr
  snp_stareg[[i]] <- snp_id[padj.eb <= q]
  
  # ReAD
  padj.hmm <- res.hmm$fdr
  snp_read[[i]] <- snp_id[padj.hmm <= q]
  
  # AdaFilter
  snp_adafilter[[i]] <- snp_id[res.adafilter$adjusted.p <= q]
  
  # # Joint mirror
  # res.jm <- JointMirror.R(cbind(pa, pb), rank.Mode = "EmptyPoset",
  #                         trgt.fdr.level = q)
  # rej.jm <- rep(0, m)
  # rej.jm[res.jm$selected] = 1
  # snp_jm <- snp_id[rej.jm]
  
  # coconut
  # x <- Cauchy(rbind(x1, x2))
  snp_coconut[[i]] <- snp_id[res.coconut$radj <= q]
  snp_coconut1[[i]] <- snp_id[res.coconut1$radj <= q]
  snp_coconut2[[i]] <- snp_id[res.coconut2$radj <= q]
  snp_coconut3[[i]] <- snp_id[res.coconut3$radj <= q]
  snp_coconut12[[i]] <- snp_id[res.coconut12$radj <= q]
  snp_coconut23[[i]] <- snp_id[res.coconut23$radj <= q]
  snp_coconut13[[i]] <- snp_id[res.coconut13$radj <= q]
  
  # carlis
  snp_carlis[[i]] <- snp_id[res.carlis$fdr <= q]
  snp_carlis1[[i]] <- snp_id[res.carlis1$fdr <= q]
  snp_carlis2[[i]] <- snp_id[res.carlis2$fdr <= q]
  snp_carlis3[[i]] <- snp_id[res.carlis3$fdr <= q]
  snp_carlis12[[i]] <- snp_id[res.carlis12$fdr <= q]
  snp_carlis23[[i]] <- snp_id[res.carlis23$fdr <= q]
  snp_carlis13[[i]] <- snp_id[res.carlis13$fdr <= q]
}


#--------------------------------------------------------------------------------
# Plot the results
#--------------------------------------------------------------------------------

library(ggplot2)
result <- data.frame(alpha = alphas, num.rej = unlist(lapply(snp_adafilter, length)), method = "AdaFilter")
result <- rbind(result, data.frame(alpha = alphas, num.rej = unlist(lapply(snp_jump, length)), method = "JUMP"))
result <- rbind(result, data.frame(alpha = alphas, num.rej = unlist(lapply(snp_stareg, length)), method = "STAREG"))
result <- rbind(result, data.frame(alpha = alphas, num.rej = unlist(lapply(snp_coconut, length)), method = "CoCoNuT-Case7"))
result <- rbind(result, data.frame(alpha = alphas, num.rej = unlist(lapply(snp_read, length)), method = "ReAD"))
result <- rbind(result, data.frame(alpha = alphas, num.rej = unlist(lapply(snp_carlis, length)), method = "CARLIS-Case7"))

result$method <- factor(result$method, levels = c("JUMP", "AdaFilter","STAREG", "ReAD", "CoCoNuT-Case7", "CARLIS-Case7"))
# 6*9
rej.plot <- ggplot(data=result, aes(x=alpha,y=num.rej, group=method, colour=method))+
  # geom_point(size=3)+
  labs(x="Nominal FDR", y="Number of rejections")+
  # geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#696969", size = 1.25)+
  geom_line(linewidth = 1.25, aes(linetype = method))+
  # scale_y_continuous(limits = c(0, 0.12), breaks = seq(0,0.1,0.02)) +
  # scale_x_continuous(trans = "log2") + #scale_y_continuous(trans = "log10") +
  # scale_linetype_manual(values = c("dashed","dotted", "dotdash", "longdash",
  #                                 "twodash", c(4, 2), c(2, 1, 2, 1), "solid"))+
  scale_colour_manual(name="", values = c("JUMP"="#FA8C5F", "AdaFilter"="#6496D2", "STAREG"="#0AB4B4",
                                          "ReAD"="#FAC800", "CoCoNuT-Case7"="#00FFFF", "CARLIS-Case7"="#FF0000"))+
  # + scale_linetype(guide=FALSE) +
  labs(color = "", linetype = "") +  #combined legend
  theme(legend.text=element_text(size=rel(1.2)),
        legend.title = element_blank(),
        legend.position = "right",
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        # panel.background = element_rect(fill = "transparent"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        panel.grid.minor = element_blank(),
        axis.line = element_line(linewidth = 0.6, linetype = "solid"))

library(UpSetR)
# Compare different methods, q = 0.001
listInput <- list(JUMP = snp_jump[[4]], AdaFilter = snp_adafilter[[4]], STAREG = snp_stareg[[4]], 
                  CoCoNuT = snp_coconut[[4]], ReAD = snp_read[[4]], CARLIS = snp_carlis[[4]])
# names(listInput1)[1] <- "MaxP-BH"
# names(listInput1)[8] <- "CoCoNuT_Case 7"
pdf(file="data analysis/BD-SCZ/upset_new.pdf",onefile = FALSE,width=7.5,height=5)
upset(fromList(listInput), nsets=6, sets.x.label = "Number of replicable SNPs", 
      order.by = "freq", text.scale = 1.2, set_size.show = TRUE, number.angles = 0, 
      set_size.scale_max = 4500, point.size = 2.5, 
      queries = list(list(query = intersects, 
                          params = list("CARLIS"),
                          active = T,
                          # color = "sandybrown",
                          query.name = "CARLIS only"),
                     list(query = intersects,
                          params = list("JUMP", "AdaFilter", "STAREG", "CoCoNuT", "ReAD", 
                                        "CARLIS"),
                          active = T,
                          query.name = "All")))
dev.off()

# Compare CARLIS with different auxiliary covariates
listInput3 <- list(JUMP = snp_jump[[4]], AdaFilter = snp_adafilter[[4]], STAREG = snp_stareg[[4]], 
                   CoCoNuT.Case7 = snp_coconut[[4]], ReAD = snp_read[[4]], 
                   CARLIS.Case1 = snp_carlis1[[4]], CARLIS.Case2 = snp_carlis2[[4]],
                   CARLIS.Case3 = snp_carlis3[[4]], CARLIS.Case4 = snp_carlis12[[4]],
                   CARLIS.Case5 = snp_carlis23[[4]], CARLIS.Case6 = snp_carlis13[[4]],
                   CARLIS.Case7 = snp_carlis[[4]])
names(listInput3) <- c("JUMP", "AdaFilter", "STAREG", "CoCoNuT-Case7", "ReAD", "CARLIS-Case1", "CARLIS-Case2", 
                       "CARLIS-Case3", "CARLIS-Case4", "CARLIS-Case5", "CARLIS-Case6", "CARLIS-Case7")
pdf(file="data analysis/BD-SCZ/upset_all.pdf",onefile = FALSE,width=12,height=6)
upset(fromList(listInput3), nsets=12, sets.x.label = "Number of pleiotropic SNPs",
      order.by = "freq", text.scale = 1.2, set_size.show = TRUE,
      set_size.scale_max = 4500, point.size = 2.5,
      queries = list(list(query = intersects, 
                          params = list("CARLIS-Case1","CARLIS-Case2", "CARLIS-Case3", "CARLIS-Case4",
                                        "CARLIS-Case5", "CARLIS-Case6", "CARLIS-Case7"),
                          active = T,
                          # color = "sandybrown",
                          query.name = "CARLIS only"),
                     list(query = intersects,
                          params = list("JUMP", "AdaFilter", "STAREG", "CoCoNuT-Case7", "ReAD", 
                                        "CARLIS-Case1", "CARLIS-Case2", "CARLIS-Case3", "CARLIS-Case4",
                                        "CARLIS-Case5", "CARLIS-Case6", "CARLIS-Case7"),
                          active = T,
                          query.name = "All")
                     ))
dev.off()


##------------------------------------------------------------------
##  GWAS results (Manhattan plot)
##------------------------------------------------------------------
gwasRes <- bip_eur[snp_id,c(1:3, 9)]
colnames(gwasRes)[4] <- "pvalue_bip"
colnames(gwasRes)[3] <- "position"
colnames(gwasRes)[2] <- "chromosome"
gwasRes$pvalue_scz <- scz_eur[snp_id, "PVAL"]
gwasRes$maxP <- pmax(gwasRes$pvalue_bip, gwasRes$pvalue_scz)
gwasRes$Lfdr <- res.eb$Lfdr
gwasRes$covLfdr <- res.coconut$Lfdr
gwasRes$rLIS <- res.hmm$repLIS
gwasRes$carlis <- res.carlis$carLIS
gwasRes$chromosome <- as.numeric(gwasRes$chromosome)
gwasRes$position <- as.numeric(gwasRes$position)

library(dplyr)
don <- gwasRes %>% 
  
  # Compute chromosome size
  group_by(chromosome) %>% 
  summarise(chr_len=max(position)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(total=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwasRes, ., by = "chromosome") %>%
  
  # Add a cumulative position of each SNP
  arrange(chromosome, position) %>%
  mutate(BPcum = position + total)

## Add highlight and annotation information
# don <- don %>% mutate( is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
#   mutate( is_annotate=ifelse(-log10(P)>4, "yes", "no")) 

axisdf = don %>% group_by(chromosome) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

q = 0.001
Rcpp::sourceCpp("src/cutoff.cpp")
cutoff.res <- cutoff(gwasRes$maxP, gwasRes$rLIS, gwasRes$Lfdr, gwasRes$covLfdr, gwasRes$carlis, q)
maxpBH.thr <- bh_cutoff(gwasRes$maxP, q)
rLIS.thr <- lfdr_cutoff(gwasRes$rLIS, q)
Lfdr.thr <- lfdr_cutoff(res.eb$Lfdr, q)
covLfdr.thr <- lfdr_cutoff(gwasRes$covLfdr, q)
carLIS.thr <- lfdr_cutoff(gwasRes$carlis, q)
jump.obj <- JUMP(pa, pb, q)
jump.thr <- jump.obj$jump.thr


library(ggplot2)
# maxpBH.thr <- cutoff.res$pBH_thr
maxP <- ggplot(don, aes(x=BPcum, y=-log10(maxP))) +
  geom_point( aes(color=as.factor(chromosome)), alpha=0.8, size=1) +
  scale_color_manual(values = rep(c("skyblue", "grey"), 22 )) +
  scale_x_continuous(label = axisdf$chromosome[c(1:9,11,13,16,21)], breaks= axisdf$center[c(1:9,11,13,16,21)]) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,20)) +     
  geom_hline(yintercept = -log10(maxpBH.thr), color = "sandybrown", size = 1.2, linetype = "dashed") +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_text(size = 13.5),
    axis.title.x = element_text(size = 13.5),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )
maxP <- maxP + labs(x = "Chromosome", y = expression(-log[10](P[max])))
ggsave(maxP, file = "./data analysis/BD-SCZ/maxP_manhattan.pdf", width = 5.5, height = 4.5)

# rLIS.thr <- cutoff.res$rLIS_thr
rLIS <- ggplot(don, aes(x=BPcum, y=-log10(rLIS))) +
  geom_point( aes(color=as.factor(chromosome)), alpha=0.8, size=1) +
  scale_color_manual(values = rep(c("skyblue", "grey"), 22 )) +
  scale_x_continuous( label = axisdf$chromosome[c(1:9,11,13,16,21)], breaks= axisdf$center[c(1:9,11,13,16,21)] ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,20)) +     
  geom_hline(yintercept = -log10(rLIS.thr), color = "sandybrown", size = 1.2, linetype = "dashed") +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_text(size = 13.5),
    axis.title.x = element_text(size = 13.5),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )
rLIS <- rLIS + labs(x = "Chromosome", y = expression(-log[10](rLIS)))
ggsave(rLIS, file = "./data analysis/BD-SCZ/rLIS_manhattan.pdf", width = 5.5, height = 4.5)


# Lfdr.thr <- cutoff.res$Lfdr_thr
Lfdr <- ggplot(don, aes(x=BPcum, y=-log10(Lfdr))) +
  geom_point( aes(color=as.factor(chromosome)), alpha=0.8, size=1) +
  scale_color_manual(values = rep(c("skyblue", "grey"), 22 )) +
  scale_x_continuous( label = axisdf$chromosome[c(1:9,11,13,16,21)], breaks= axisdf$center[c(1:9,11,13,16,21)] ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,20)) +     
  geom_hline(yintercept = -log10(Lfdr.thr), color = "sandybrown", size = 1.2, linetype = "dashed") +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_text(size = 13.5),
    axis.title.x = element_text(size = 13.5),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )
Lfdr <- Lfdr + labs(x = "Chromosome", y = expression(-log[10](Lfdr)))
ggsave(Lfdr, file = "./data analysis/BD-SCZ/Lfdr_manhattan.pdf", width = 5.5, height = 4.5)

# covLfdr.thr <- cutoff.res$covLfdr_thr
covLfdr <- ggplot(don, aes(x=BPcum, y=-log10(covLfdr))) +
  geom_point( aes(color=as.factor(chromosome)), alpha=0.8, size=1) +
  scale_color_manual(values = rep(c("skyblue", "grey"), 22 )) +
  scale_x_continuous( label = axisdf$chromosome[c(1:9,11,13,16,21)], breaks= axisdf$center[c(1:9,11,13,16,21)] ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,20)) +     
  geom_hline(yintercept = -log10(covLfdr.thr), color = "sandybrown", size = 1.2, linetype = "dashed") +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_text(size = 13.5),
    axis.title.x = element_text(size = 13.5),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )
covLfdr <- covLfdr + labs(x = "Chromosome", y = expression(-log[10](covLfdr)))
ggsave(covLfdr, file = "./data analysis/BD-SCZ/covLfdr_manhattan.pdf", width = 5.5, height = 4.5)

# carLIS.thr <- cutoff.res$carLIS_thr
carLIS <- ggplot(don, aes(x=BPcum, y=-log10(carlis))) +
  geom_point( aes(color=as.factor(chromosome)), alpha=0.8, size=1) +
  scale_color_manual(values = rep(c("skyblue", "grey"), 22 )) +
  scale_x_continuous( label = axisdf$chromosome[c(1:9,11,13,16,21)], breaks= axisdf$center[c(1:9,11,13,16,21)] ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,20)) +     
  geom_hline(yintercept = -log10(carLIS.thr), color = "sandybrown", size = 1.2, linetype = "dashed") +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_text(size = 13.5),
    axis.title.x = element_text(size = 13.5),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )
carLIS <- carLIS + labs(x = "Chromosome", y = expression(-log[10](CARLIS)))
ggsave(carLIS, file = "./data analysis/BD-SCZ/carLIS_manhattan.pdf", width = 5.5, height = 4.5)



