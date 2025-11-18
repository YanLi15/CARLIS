

# NDKP GP2 Parkinson's disease 2025 GWAS European non-Finnish ancestry
pd_eur <- read.table("./data analysis/Parkinson's disease/NDKP-GP2 Parkinson's disease 2025 GWAS European ancestry/GP2_euro_ancestry_meta_analysis_2024/GP2_EUR_ONLY_HG38_12162024.txt.gz",
                               quote = "", fill = TRUE, comment.char = "!", header = TRUE)
rownames(pd_eur) <- pd_eur$SNP_ID

# NDKP GP2 Parkinson's disease 2025 GWAS FIN ancestry
pd_fin <- read.table("./data analysis/Parkinson's disease/NDKP-GP2 Parkinson's disease 2025 GWAS European ancestry/GP2_euro_ancestry_meta_analysis_2024/GP2_FIN_ONLY_HG38_12162024.txt.gz",
                             quote = "", fill = TRUE, comment.char = "!", header = TRUE)
rownames(pd_fin) <- pd_fin$SNP_ID

# NDKP GP2 Parkinson's disease 2025 GWAS Ashkenazi Jewish ancestry
pd_aj <- read.table("./data analysis/Parkinson's disease/NDKP-GP2 Parkinson's disease 2025 GWAS European ancestry/GP2_euro_ancestry_meta_analysis_2024/GP2_AJ_ONLY_HG38_12162024.txt.gz",
                     quote = "", fill = TRUE, comment.char = "!", header = TRUE)
rownames(pd_aj) <- pd_aj$SNP_ID


snp_id <- intersect(pd_eur$SNP_ID, pd_fin$SNP_ID)
snp_id <- intersect(snp_id, pd_aj$SNP)

pa <- pd_eur[snp_id,]$p_value
pb <- pd_fin[snp_id,]$p_value
x <- pd_aj[snp_id,]$p_value


library(JUMP)
library(radjust)
library(STAREG)
library(adaFilter)
library(JMirror)
library(AdaPTGMM)
library(CARLIS)
library(ReAD)
library(CoCoNuT)
source("./R/methods.R")

alphas = c(1e-10, 5e-10, 1e-9, 5e-9, 1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3)

# STAREG
res.eb <- stareg(pa, pb)
# ReAD
res.hmm <- RepLIS(pa, pb)
# AdaFilter
res.adafilter <- adaFilter(cbind(pa, pb), r = 2, type.I.err = "FDR")
# coconut
res.coconut <- coconut(pa, pb, x)
# carlis
res.carlis <- CARLIS(pa, pb, x)

snp_jump <- list()
snp_adafilter <- list()
snp_carlis <- list()
snp_coconut <- list()
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
  snp_coconut[[i]] <- snp_id[res.coconut$radj <= q]
  
  # carlis
  snp_carlis[[i]] <- snp_id[res.carlis$fdr <= q]
}

library(ggplot2)
result <- data.frame(alpha = alphas[6:16], num.rej = unlist(lapply(snp_adafilter, length))[6:16], method = "AdaFilter")
result <- rbind(result, data.frame(alpha = alphas[6:16], num.rej = unlist(lapply(snp_jump, length))[6:16], method = "JUMP"))
result <- rbind(result, data.frame(alpha = alphas[6:16], num.rej = unlist(lapply(snp_stareg, length))[6:16], method = "STAREG"))
result <- rbind(result, data.frame(alpha = alphas[6:16], num.rej = unlist(lapply(snp_coconut, length))[6:16], method = "CoCoNuT"))
result <- rbind(result, data.frame(alpha = alphas[6:16], num.rej = unlist(lapply(snp_read, length))[6:16], method = "ReAD"))
result <- rbind(result, data.frame(alpha = alphas[6:16], num.rej = unlist(lapply(snp_carlis, length))[6:16], method = "CARLIS"))
# 5e-8,1e-7,5e-7,1e-6,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3

result$method <- factor(result$method, levels = c("JUMP", "AdaFilter","STAREG", "ReAD", "CoCoNuT", "CARLIS"))
rej.plot <- ggplot(data=result, aes(x=alpha,y=num.rej, group=method, colour=method))+
  # geom_point(size=3)+
  labs(x="Nominal FDR", y="Number of rejections")+
  # geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#696969", size = 1.25)+
  geom_line(linewidth = 1.25, aes(linetype = method))+
  #scale_y_continuous(limits = c(0, 0.12), breaks = seq(0,0.1,0.02)) +
  scale_x_continuous(trans = "log10") + #scale_y_continuous(trans = "log10") +
  # scale_linetype_manual(values = c("dashed","dotted", "dotdash", "longdash",
  #                                 "twodash", c(4, 2), c(2, 1, 2, 1), "solid"))+
  scale_colour_manual(name="", values = c("JUMP"="#FA8C5F", "AdaFilter"="#6496D2", "STAREG"="#0AB4B4", 
                                          "ReAD"="#FAC800", "CoCoNuT"="#00FFFF", "CARLIS"="#FF0000"))+
  # + scale_linetype(guide=FALSE) +
  labs(color = "", linetype = "") +  #combined legend
  theme(legend.text=element_text(size=rel(1.2)),
        legend.title = element_blank(),
        legend.position = "bottom",
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        # panel.background = element_rect(fill = "transparent"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        panel.grid.minor = element_blank(),
        axis.line = element_line(linewidth = 0.6, linetype = "solid"))


#--------------------------------------------------------------------------------
# Plot the results
#--------------------------------------------------------------------------------
library(UpSetR)
# Compare different methods, q = 1e-7
listInput <- list(JUMP = snp_jump[[7]], AdaFilter = snp_adafilter[[7]], STAREG = snp_stareg[[7]], 
                   CoCoNuT = snp_coconut[[7]], ReAD = snp_read[[7]], CARLIS = snp_carlis[[7]])
# names(listInput1)[1] <- "MaxP-BH"
# names(listInput1)[8] <- "CoCoNuT_Case 7"
pdf(file="data analysis/Parkinson's disease/upset.pdf",onefile = FALSE,width=7.5,height=5)
upset(fromList(listInput), nsets=6, sets.x.label = "Number of replicable SNPs", 
      order.by = "freq", text.scale = 1.2, set_size.show = TRUE, 
      set_size.scale_max = 6000, point.size = 2.5, 
      queries = list(list(query = intersects, 
                          params = list("CARLIS"),
                          active = T,
                          query.name = "CARLIS only"),
                     list(query = intersects,
                          params = list("JUMP", "AdaFilter", "STAREG", "CoCoNuT", "ReAD", "CARLIS"),
                          active = T,
                          query.name = "All")))
dev.off()


##------------------------------------------------------------------
##  GWAS results (Manhattan plot)
##------------------------------------------------------------------
gwasRes.pd <- pd_eur[snp_id,c(1:3, 8)]
colnames(gwasRes.pd)[4] <- "pvalue_eur"
colnames(gwasRes.pd)[2] <- "position"
gwasRes.pd$pvalue_fin <- pd_fin[snp_id, "p_value"]
gwasRes.pd$maxP <- pmax(gwasRes.pd$pvalue_eur, gwasRes.pd$pvalue_fin)
gwasRes.pd$Lfdr <- res.eb$Lfdr
gwasRes.pd$covLfdr <- res.coconut$Lfdr
gwasRes.pd$rLIS <- res.hmm$repLIS
gwasRes.pd$carlis <- res.carlis$carLIS
gwasRes.pd$chromosome <- as.numeric(gwasRes.pd$chromosome)
gwasRes.pd$position <- as.numeric(gwasRes.pd$position)
gwasRes.pd$pa <- pa
gwasRes.pd$pb <- pb

library(dplyr)
don <- gwasRes.pd %>% 
  
  # Compute chromosome size
  group_by(chromosome) %>% 
  summarise(chr_len=max(position)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(total=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwasRes.pd, ., by = "chromosome") %>%
  
  # Add a cumulative position of each SNP
  arrange(chromosome, position) %>%
  mutate(BPcum = position + total)

## Add highlight and annotation information
# don <- don %>% mutate( is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
#   mutate( is_annotate=ifelse(-log10(P)>4, "yes", "no")) 

axisdf = don %>% group_by(chromosome) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

q = 1e-7
Rcpp::sourceCpp("src/cutoff.cpp")
cutoff.res <- cutoff(gwasRes.pd$maxP, gwasRes.pd$rLIS, gwasRes.pd$Lfdr, gwasRes.pd$covLfdr, gwasRes.pd$carlis, q)
rLIS.thr <- lfdr_cutoff(gwasRes.pd$rLIS, q)
Lfdr.thr <- lfdr_cutoff(res.eb$Lfdr, q)
covLfdr.thr <- lfdr_cutoff(gwasRes.pd$covLfdr, q)
carLIS.thr <- lfdr_cutoff(gwasRes.pd$carlis, q)
jump.obj <- JUMP(pa, pb, q)
jump.thr <- jump.obj$jump.thr

library(ggplot2)
# maxpBH.thr <- cutoff.res$pBH_thr
maxP <- ggplot(don, aes(x=BPcum, y=-log10(maxP))) +
  geom_point( aes(color=as.factor(chromosome)), alpha=0.8, size=1) +
  scale_color_manual(values = rep(c("skyblue", "grey"), 22 )) +
  scale_x_continuous(label = axisdf$chromosome[c(1:9,11,13,16,21)], breaks= axisdf$center[c(1:9,11,13,16,21)]) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,25)) +     # remove space between plot area and x axis
  #ylab(expression(-log[10])(P[max])) + xlab("Chromosome")
  # # Add highlighted points
  # geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
  # # Add label using ggrepel to avoid overlapping
  # geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=SNP), size=2) +
  geom_hline(yintercept = -log10(jump.thr), color = "sandybrown", size = 1.2, linetype = "dashed") +
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
ggsave(maxP, file = "./data analysis/Parkinson's disease/maxP_manhattan.pdf", width = 5.5, height = 4.5)

# rLIS.thr <- cutoff.res$rLIS_thr
rLIS <- ggplot(don, aes(x=BPcum, y=-log10(rLIS))) +
  geom_point( aes(color=as.factor(chromosome)), alpha=0.8, size=1) +
  scale_color_manual(values = rep(c("skyblue", "grey"), 22 )) +
  scale_x_continuous( label = axisdf$chromosome[c(1:9,11,13,16,21)], breaks= axisdf$center[c(1:9,11,13,16,21)] ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,25)) +     # remove space between plot area and x axis
  # # Add highlighted points
  # geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
  # # Add label using ggrepel to avoid overlapping
  # geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=SNP), size=2) +
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
# repLIS <- rLIS + labs(x = "Chromosome", y = expression(-log[10](repLIS)))
ggsave(rLIS, file = "./data analysis/Parkinson's disease/rLIS_manhattan.pdf", width = 5.5, height = 4.5)
# ggsave(repLIS, file = "./data analysis/Asthma/repLIS_manhattan.pdf", width = 5.5, height = 4.5)


# Lfdr.thr <- cutoff.res$Lfdr_thr
Lfdr <- ggplot(don, aes(x=BPcum, y=-log10(Lfdr))) +
  geom_point( aes(color=as.factor(chromosome)), alpha=0.8, size=1) +
  scale_color_manual(values = rep(c("skyblue", "grey"), 22 )) +
  scale_x_continuous( label = axisdf$chromosome[c(1:9,11,13,16,21)], breaks= axisdf$center[c(1:9,11,13,16,21)] ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,25)) +     # remove space between plot area and x axis
  # # Add highlighted points
  # geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
  # # Add label using ggrepel to avoid overlapping
  # geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=SNP), size=2) +
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
ggsave(Lfdr, file = "./data analysis/Parkinson's disease/Lfdr_manhattan.pdf", width = 5.5, height = 4.5)

# covLfdr.thr <- cutoff.res$covLfdr_thr
covLfdr <- ggplot(don, aes(x=BPcum, y=-log10(covLfdr))) +
  geom_point( aes(color=as.factor(chromosome)), alpha=0.8, size=1) +
  scale_color_manual(values = rep(c("skyblue", "grey"), 22 )) +
  scale_x_continuous( label = axisdf$chromosome[c(1:9,11,13,16,21)], breaks= axisdf$center[c(1:9,11,13,16,21)] ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,25)) +     # remove space between plot area and x axis
  # # Add highlighted points
  # geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
  # # Add label using ggrepel to avoid overlapping
  # geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=SNP), size=2) +
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
ggsave(covLfdr, file = "./data analysis/Parkinson's disease/covLfdr_manhattan.pdf", width = 5.5, height = 4.5)

# carLIS.thr <- cutoff.res$carLIS_thr
carLIS <- ggplot(don, aes(x=BPcum, y=-log10(carlis))) +
  geom_point( aes(color=as.factor(chromosome)), alpha=0.8, size=1) +
  scale_color_manual(values = rep(c("skyblue", "grey"), 22 )) +
  scale_x_continuous( label = axisdf$chromosome[c(1:9,11,13,16,21)], breaks= axisdf$center[c(1:9,11,13,16,21)] ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,25)) +     # remove space between plot area and x axis
  # # Add highlighted points
  # geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
  # # Add label using ggrepel to avoid overlapping
  # geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=SNP), size=2) +
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
ggsave(carLIS, file = "./data analysis/Parkinson's disease/carLIS_manhattan.pdf", width = 5.5, height = 4.5)
