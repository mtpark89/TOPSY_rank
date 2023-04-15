library(readr)
library(magrittr)
library(dplyr)
library(ggplot2)
library(ggridges)
library(RMINC)
library(matrixStats)
library(Rankcluster)
library(ggpubr)
library(reshape2)

################################################################
###Read in demographics

master <- read.csv('/projects/TOPSY_NBM/Analysis_final/TOPSY_master_NBM_20201122.csv', sep="@")

master$DX <-gsub('Healthy Control', 'HC',
	gsub('First Episode Patient', 'FEP',
	gsub('3\\+ Year Patient', '3+year', master$PatientCat)))

gf <- master

gf$DX<- relevel(factor(gf$DX), ref="HC")
gf <- subset(gf, DX!="NA")

mrs <- read_tsv("TOPSY_MRS_20210530.csv")
gf <- left_join(gf, mrs %>% select(-ConsensusDx, ), by="ID")

#NBM file already has QC data
#QC_CIVET <- read.csv("/projects/TOPSY_final_20200504/QC/QC_CIVET_20200707.csv")
#gf <- left_join(gf, QC_CIVET, by="Subject")

#Subset for comparison analysis with other datasets
master$Dataset <- "TOPSY"

master %<>% 
    mutate(
	PANSS_total = select(., starts_with("PANSS8")) %>% rowSums(),
        PANSS_pos = select(., starts_with("PANSS8P")) %>% rowSums(),
        PANSS_neg = select(., starts_with("PANSS8N")) %>% rowSums(),
	PANSS_g = select(., starts_with("PANSS8G")) %>% rowSums(),
    )

master$CIVET_ID <- paste0(master$Dataset, "_", master$Subject, sep="")

write_tsv(master %>% select(Dataset, Subject, CIVET_ID, Age, Sex, DX, QC_CIVET, PANSS_total, PANSS_pos, PANSS_neg, PANSS_g), file="TOPSY_subset_20230327.tsv")

################################################################
###CIVET vertex-wise files setup

gf$left_ct <- paste("/projects/Datasets/CIVET_TOPSY/thickness/TOPSY_",gf$Subject,"_native_rms_rsl_tlaplace_30mm_left.txt", sep="")
gf$right_ct <- paste("/projects/Datasets/CIVET_TOPSY/thickness/TOPSY_",gf$Subject,"_native_rms_rsl_tlaplace_30mm_right.txt", sep="")

gf_master <- gf

###Filters
gf<- subset(gf, QC_CIVET=="Pass" | QC_CIVET=="Partial")
gf <- left_join(gf, euler, by="Subject")

###Read parcellations map
stat_maps <- read_delim("stat_maps.txt", delim=",")
stat_maps_left <- head(stat_maps, n=40962)
stat_maps_right <- tail(stat_maps, n=40962)

###Subsets
gf_HC <- subset(gf, DX=="HC")

yeo7_left <-  c("Unmatched_left", "Visual_left", "Somatomotor_left", "DorsalAtt_left", "VentralAtt_left", "Limbic_yeo_left", "Frontoparietal_left", "Default_left")
yeo7_right <-  c("Unmatched_right", "Visual_right", "Somatomotor_right", "DorsalAtt_right", "VentralAtt_right", "Limbic_yeo_right", "Frontoparietal_right", "Default_right")

vonE_left <-  c("Unmatched_left", "Motor_left", "Association-FP_left", "Association-FT_left", "Sensory-2_left", "Sensory-1_left", "Limbic_left", "Insular_left")
vonE_right <-  c("Unmatched_right", "Motor_right", "Association-FP_right", "Association-FT_right", "Sensory-2_right", "Sensory-1_right", "Limbic_right", "Insular_right")


###########################################################################################################
######ANALYSIS#####

###1. Visual inspection of cortical thickness & surface area ranks

left_ct_rank <- t(vertexTableRank(gf$left_ct))
right_ct_rank <- t(vertexTableRank(gf$right_ct))

left_ct <- vertexTable(gf$left_ct)
right_ct <- vertexTable(gf$right_ct)

writeVertex(left_ct_rank, "results/left_ct_rank.vertstats")
writeVertex(right_ct_rank, "results/right_ct_rank.vertstats")

writeVertex(cbind(stat_maps, rbind(left_ct_rank, right_ct_rank)), "results/ct_rank.vertstats")
writeVertex(cbind(stat_maps, rbind(left_ct, right_ct)), "results/ct_raw.vertstats")

###Sample subject & group mean, raw and ranked CT
sample_subject_left <- cbind(left_ct[,1], left_ct_rank[,1])
sample_subject_right <- cbind(right_ct[,1], right_ct_rank[,1])

colnames(sample_subject_left) <- c("sub_left_ct", "sub_left_ct_rank")
colnames(sample_subject_right) <- c("sub_right_ct", "sub_right_ct_rank")

#Mean CT, and mean CT then ranked
group_mean_left <- cbind(rowMeans(left_ct), rank(rowMeans(left_ct), ties.method="first"))
group_mean_right <- cbind(rowMeans(right_ct), rank(rowMeans(right_ct), ties.method="first"))

colnames(group_mean_left) <- c("left_ct", "left_ct_rank")
colnames(group_mean_right) <- c("right_ct", "right_ct_rank")

writeVertex(cbind(sample_subject_left, group_mean_left, stat_maps_left), "results/subject_group_mean_left.vertstats")
writeVertex(cbind(sample_subject_right, group_mean_right, stat_maps_right), "results/subject_group_mean_right.vertstats")

###Rank standard deviation: all, HC, FEP.

gf_FEP <- subset(gf, DX=="FEP")
gf_HC <- subset(gf, DX=="HC")

left_ct_rank <- t(vertexTableRank(gf$left_ct))
right_ct_rank <- t(vertexTableRank(gf$right_ct))

left_ct_rank_HC <- t(vertexTableRank(gf_HC$left_ct))
right_ct_rank_HC <- t(vertexTableRank(gf_HC$right_ct))

left_ct_rank_FEP <- t(vertexTableRank(gf_FEP$left_ct))
right_ct_rank_FEP <- t(vertexTableRank(gf_FEP$right_ct))

writeVertex(cbind(rowSds(left_ct_rank), rowSds(left_ct_rank_HC), rowSds(left_ct_rank_FEP), stat_maps_left), "results/left_ct_rank_sd.vertstats")
writeVertex(cbind(rowSds(right_ct_rank), rowSds(right_ct_rank_HC), rowSds(right_ct_rank_FEP), stat_maps_right), "results/right_ct_rank_sd.vertstats")

###2. Cortical rank vs. subregions

group_left_ct_final <- cbind(group_mean_left, rank(rowMedians(left_ct_rank), ties.method="first"), stat_maps_left)
colnames(group_left_ct_final)[1:3] <- c("mean_ct", "mean_ct_rank", "median_ct_rank")

writeVertex(group_left_ct_final, file="results/group_left_ct_final.vertstats")

###Ridges plot

plot1 <- ggplot(group_left_ct_final, aes(x=median_ct_rank, y=as.factor(yeo7))) + geom_density_ridges(alpha=1, size=1) + stat_density_ridges(quantile_lines=TRUE, quantiles=0.5) + 
	theme_ridges(grid=FALSE, center_axis_labels=TRUE) +
	xlab("") + ylab("Yeo Networks") +
	scale_y_discrete(labels=c("Unmatched", "Visual", "Somatomotor", "DorsalAtt", "VentralAtt", "Limbic", "Frontoparietal", "Default")) +
	theme_minimal() + theme(text=element_text(size=22, family="Helvetica"), axis.text.x=element_text(angle=45, hjust=0.35), axis.title.x=element_blank(), axis.text=element_text(size=16))

plot2 <- ggplot(group_left_ct_final, aes(x=median_ct_rank, y=as.factor(vonEconomo))) + geom_density_ridges(alpha=1, size=1) + stat_density_ridges(quantile_lines=TRUE, quantiles=0.5) + 
	theme_ridges(grid=FALSE, center_axis_labels=TRUE) +
	xlab("") + ylab("von Economo class") +
	scale_y_discrete(labels=c("Unmatched", "Primary Motor", "Association", "Association", "Secondary sensory", "Primary sensory", "Limbic", "Insular")) +
	theme_minimal() + theme(text=element_text(size=22, family="Helvetica"), axis.text.x=element_text(angle=45, hjust=0.35), axis.title.x=element_blank(), axis.text=element_text(size=16))

ggsave("Plot_networks_ranks.png", ggarrange(plot1, plot2, ncol=2, label.x = "Rank"), width=10, height=5, bg="white")

###3. Comparing models: raw vs. ranked

gf_compare <- gf %>% filter(DX=="FEP" | DX=="HC")
gf_compare$DX %<>% droplevels()
dataset <- gf_compare

left_ct_compare <- compare_models(gf_compare, "left_ct", "DX")
right_ct_compare <- compare_models(gf_compare, "right_ct", "DX")

writeVertex(left_ct_compare, file="results/left_ct_compare.vertstats")
writeVertex(right_ct_compare, file="results/right_ct_compare.vertstats")

###4. Group differences in average rank

#Using permutation testing
left_ct_rank_diff <- average_rank_diff_permute(gf_HC$left_ct, gf_FEP$left_ct, 10000)
right_ct_rank_diff <- average_rank_diff_permute(gf_HC$right_ct, gf_FEP$right_ct, 10000)

left_ct_rank_diff$logpval <- -log10(left_ct_rank_diff$pvals) * sign(left_ct_rank_diff$diff)
right_ct_rank_diff$logpval <- -log10(right_ct_rank_diff$pvals) * sign(right_ct_rank_diff$diff)

writeVertex(left_ct_rank_diff, "results/left_ct_rank_diff.vertstats")
writeVertex(right_ct_rank_diff, "results/right_ct_rank_diff.vertstats")

#Plotting left & right hemisphere rank differences histogram
left_ct_rank_diff$hemi <- "Left CT"
right_ct_rank_diff$hemi <- "Right CT"

ct_rank_diff <- rbind(left_ct_rank_diff, right_ct_rank_diff)

sd(ct_rank_diff$diff) #SD is 1648.211, *2 is 3296.

plot_histo <- ggplot(ct_rank_diff, aes(x=diff, group=hemi, fill=hemi)) + geom_density(alpha=0.5, size=1) + 
	xlab("Mean rank difference") + ylab("Density") +  theme_minimal() + scale_fill_discrete() +
	theme(text=element_text(size=22, family="Helvetica"), axis.text.x=element_text(angle=45), axis.text=element_text(size=16), legend.title=element_blank()) +
	geom_vline(xintercept=3296, linetype="dashed", colour="darkgray", size=2) + geom_vline(xintercept=-3296, linetype="dashed", colour="darkgray", size=2)

ggsave("Plot_histo_mean_rankdiff.png", plot_histo, width=5, height=5, bg="white")

#Comparing to linear model with raw CT
vertexFDR(vertexLm(left_ct ~ DX, data=gf_compare))
vertexFDR(vertexLm(left_ct ~ Age+ Sex + DX, data=gf_compare))
vertexFDR(vertexLm(right_ct ~ DX, data=gf_compare))
vertexFDR(vertexLm(right_ct ~ Age+ Sex + DX, data=gf_compare))

left_ct <- vertexTable(gf_compare$left_ct)
right_ct <- vertexTable(gf_compare$right_ct)

left_ct_raw_lm <- data.frame()
right_ct_raw_lm <- data.frame()

for (i in 1:nrow(left_ct)){
	model1 <- parseLm(lm(left_ct[i,] ~ gf_compare$DX))
	model2 <- parseLm(lm(right_ct[i,] ~ gf_compare$DX))
	left_ct_raw_lm <- rbind(left_ct_raw_lm, model1)
	right_ct_raw_lm <- rbind(right_ct_raw_lm, model2)
}

left_ct_raw_lm$qvalue <- p.adjust(left_ct_raw_lm$`p_gf_compare$DXFEP`, method="fdr")
right_ct_raw_lm$qvalue <- p.adjust(right_ct_raw_lm$`p_gf_compare$DXFEP`, method="fdr")

left_ct_raw_lm$logpval <- -log10(left_ct_raw_lm$`p_gf_compare$DXFEP`) * sign(left_ct_raw_lm$`t_gf_compare$DXFEP`)
right_ct_raw_lm$logpval <- -log10(right_ct_raw_lm$`p_gf_compare$DXFEP`) * sign(right_ct_raw_lm$`t_gf_compare$DXFEP`)

colnames(left_ct_raw_lm) %<>% gsub("\\$DXFEP", "", .)
colnames(right_ct_raw_lm) %<>% gsub("\\$DXFEP", "", .)
writeVertex(left_ct_raw_lm, "results/left_ct_raw_lm.vertstats")
writeVertex(right_ct_raw_lm, "results/right_ct_raw_lm.vertstats")

###5. Median rank testing

left_ct_rank <- t(vertexTableRank(gf$left_ct))
right_ct_rank <- t(vertexTableRank(gf$right_ct))

left_ct_yeo <- getMedians(left_ct_rank, stat_maps_left$yeo7)
right_ct_yeo <- getMedians(right_ct_rank, stat_maps_right$yeo7)

left_ct_vonEconomo <- getMedians(left_ct_rank, stat_maps_left$vonEconomo)
right_ct_vonEconomo <- getMedians(right_ct_rank, stat_maps_right$vonEconomo)

colnames(left_ct_yeo) <-  c("Unmatched_left", "Visual_left", "Somatomotor_left", "DorsalAtt_left", "VentralAtt_left", "Limbic_yeo_left", "Frontoparietal_left", "Default_left")
colnames(right_ct_yeo) <-  c("Unmatched_right", "Visual_right", "Somatomotor_right", "DorsalAtt_right", "VentralAtt_right", "Limbic_yeo_right", "Frontoparietal_right", "Default_right")

colnames(left_ct_vonEconomo) <-  c("Unmatched_left", "Motor_left", "Association-FP_left", "Association-FT_left", "Sensory-2_left", "Sensory-1_left", "Limbic_left", "Insular_left")
colnames(right_ct_vonEconomo) <-  c("Unmatched_right", "Motor_right", "Association-FP_right", "Association-FT_right", "Sensory-2_right", "Sensory-1_right", "Limbic_right", "Insular_right")

gf_medians <- cbind(gf, left_ct_yeo, right_ct_yeo)
gf_medians_FEP <- subset(gf_medians, DX=="FEP")

##################################################

getCors(gf_medians_FEP$totalP, gf_medians_FEP %>% select(ends_with(c("_left", "_right")))) %>% write_tsv(., "results/results_correlations_median_totalP.tsv")
getCors(gf_medians_FEP$totalN, gf_medians_FEP %>% select(ends_with(c("_left", "_right")))) %>% write_tsv(., "results/results_correlations_median_totalN.tsv")
getCors(gf_medians_FEP$WksTo50P, gf_medians_FEP %>% select(ends_with(c("_left", "_right")))) %>% write_tsv(., "results/results_correlations_median_WksTo50P.tsv")
getCors(gf_medians_FEP$WksToCGI2, gf_medians_FEP %>% select(ends_with(c("_left", "_right")))) %>% write_tsv(., "results/results_correlations_median_WksToCGI2.tsv")

getLms(gf_medians_FEP %>% select(ends_with(c("_left", "_right"))), "Age + Sex + euler_total + totalP", gf_medians_FEP)
getLms(gf_medians_FEP %>% select(ends_with(c("_left", "_right"))), "Age + Sex + euler_total + totalN", gf_medians_FEP)
getLms(gf_medians_FEP %>% select(ends_with(c("_left", "_right"))), "Age + Sex + euler_total + WksTo50P", gf_medians_FEP)
getLms(gf_medians_FEP %>% select(ends_with(c("_left", "_right"))), "Age + Sex + euler_total + WksToCGI2", gf_medians_FEP)

summary(lm(gf_medians_FEP$DorsalAtt_left ~ Age + Sex + euler_total + WksToCGI2, data=gf_medians_FEP))
summary(lm(gf_medians_FEP$DorsalAtt_right ~ Age + Sex + euler_total + WksToCGI2, data=gf_medians_FEP))

summary(lm(gf_medians_FEP$DorsalAtt_left ~ Age + Sex + euler_total + WksTo50P, data=gf_medians_FEP))
summary(lm(gf_medians_FEP$DorsalAtt_right ~ Age + Sex + euler_total + WksTo50P, data=gf_medians_FEP))

summary(lm(gf_medians_FEP$Default_left ~ Age + Sex + euler_total + totalP, data=gf_medians_FEP))

melted <- reshape2::melt(gf_medians_FEP %>% select(DorsalAtt_left, DorsalAtt_right, WksToCGI2), id="WksToCGI2", variable.name="hemisphere")
melted$hemisphere <-gsub('DorsalAtt_left', 'Left',
		gsub('DorsalAtt_right', 'Right', melted$hemisphere))

plot_cgi <- ggplot(melted %>% filter(!is.na(WksToCGI2)), aes(x=value, y=WksToCGI2, group=hemisphere, colour=hemisphere)) + 
    geom_point(size=2) + theme_minimal() + geom_smooth(method=lm, size=2) +
    xlab("Median rank: Dorsal Attention network") + ylab("Time to reach CGI-S score 2 (weeks)") +
    theme(text=element_text(size=16, family="Helvetica"), axis.text.x=element_text(hjust=0.35), axis.text=element_text(size=14), legend.title=element_blank(), legend.position = "NONE")  + scale_colour_Publication()

ggsave("figures/Plot_cgi_dorsal.png", plot_cgi, width=6, height=5, bg="white")

melted <- reshape2::melt(gf_medians_FEP %>% select(DorsalAtt_left, DorsalAtt_right, WksTo50P), id="WksTo50P", variable.name="hemisphere")
melted$hemisphere <-gsub('DorsalAtt_left', 'Left',
		gsub('DorsalAtt_right', 'Right', melted$hemisphere))

plot_50p <- ggplot(melted %>% filter(!is.na(WksTo50P)), aes(x=value, y=WksTo50P, group=hemisphere, colour=hemisphere)) + 
    geom_point(size=2) + theme_minimal() + geom_smooth(method=lm, size=2) +
    xlab("Median rank: Dorsal Attention network") + ylab("Time to reach 50% PANSS (weeks)") +
    theme(text=element_text(size=16, family="Helvetica"), axis.text.x=element_text(hjust=0.35), axis.text=element_text(size=14), legend.title=element_blank(), legend.position = "NONE") + scale_colour_Publication()

ggsave("figures/Plot_50p_dorsal.png", plot_50p, width=6, height=5, bg="white")

plot_totalP1 <- ggplot(gf_medians_FEP, aes(x=Default_left, y=totalP)) + 
	geom_point(size=2) + theme_minimal() + geom_smooth(method=lm, size=2) +
	xlab("Median rank: Default network") + ylab("Total positive symptoms") +
	theme(text=element_text(size=16, family="Helvetica"), axis.text.x=element_text(hjust=0.35), axis.text=element_text(size=14), legend.title=element_blank())

plot_totalP2 <- ggplot(gf_medians_FEP, aes(x=VentralAtt_left, y=totalP)) + 
	geom_point(size=2) + theme_minimal() + geom_smooth(method=lm, size=2) +
	xlab("Median rank: Ventral attention network") + ylab("Total positive symptoms") +
	theme(text=element_text(size=16, family="Helvetica"), axis.text.x=element_text(hjust=0.35), axis.text=element_text(size=14), legend.title=element_blank())

ggsave("figures/Plot_totalP_default.png", plot_totalP1, width=6, height=5, bg="white")
ggsave("figures/Plot_totalP_ventralatt.png", plot_totalP2, width=6, height=5, bg="white")

plot_TOPSY_medians <- ggarrange(plot_totalP1, plot_cgi, plot_50p, ncol=2, nrow=2)

ggsave("figures/Plot_TOPSY_medians.png", plot_TOPSY_medians, width=12, height=10, bg="white")
####################################
###Comparison: raw cortical thickness & network medians

left_ct <- vertexTable(gf$left_ct)
right_ct <- vertexTable(gf$right_ct)

left_ct_yeo <- getMedians(left_ct, stat_maps_left$yeo7)
right_ct_yeo <- getMedians(right_ct, stat_maps_right$yeo7)

gf_medians_raw <- cbind(gf, left_ct_yeo, right_ct_yeo)
gf_medians_raw_FEP <- subset(gf_medians_raw, DX=="FEP")

getCors(gf_medians_raw_FEP$totalP, gf_medians_raw_FEP %>% select(ends_with(c("_left", "_right"))) ) %>% filter(pvalue < 0.05)
getCors(gf_medians_raw_FEP$totalN, gf_medians_raw_FEP %>% select(ends_with(c("_left", "_right"))) ) %>% filter(pvalue < 0.05)
getCors(gf_medians_raw_FEP$WksTo50P, gf_medians_raw_FEP %>% select(ends_with(c("_left", "_right"))) ) %>% filter(pvalue < 0.05)
getCors(gf_medians_raw_FEP$WksToCGI2, gf_medians_raw_FEP %>% select(ends_with(c("_left", "_right"))) ) %>% filter(pvalue < 0.05)

#####################################################

###6. Rank-based normative modelling: using median rank of HCs to get cayley distances between individual & median network ranks

cayley_ct_left_hcmed_yeo7 <- cayley_dist_ref(gf_HC$left_ct, gf$left_ct, stat_maps_left$yeo7)
cayley_ct_right_hcmed_yeo7 <- cayley_dist_ref(gf_HC$right_ct, gf$right_ct, stat_maps_right$yeo7)

###Rank-based normative modelling: using all reference subjects, time consuming, but better results.

cayley_ct_left_yeo7 <- indiv_distance_cayley_sub(gf_HC$left_ct, gf$left_ct, stat_maps_left$yeo7)
cayley_ct_left_yeo17 <- indiv_distance_cayley_sub(gf_HC$left_ct, gf$left_ct, stat_maps_left$yeo17)
cayley_ct_left_vonEconomo <- indiv_distance_cayley_sub(gf_HC$left_ct, gf$left_ct, stat_maps_left$vonEconomo)
cayley_ct_right_yeo7 <- indiv_distance_cayley_sub(gf_HC$right_ct, gf$right_ct, stat_maps_right$yeo7)
cayley_ct_right_yeo17 <- indiv_distance_cayley_sub(gf_HC$right_ct, gf$right_ct, stat_maps_right$yeo17)
cayley_ct_right_vonEconomo <- indiv_distance_cayley_sub(gf_HC$right_ct, gf$right_ct, stat_maps_right$vonEconomo)

write.table(cayley_ct_left_yeo7, "results/cayley_ct_left_yeo7.csv", sep="\t", row.names=F)
write.table(cayley_ct_left_yeo17, "results/cayley_ct_left_yeo17.csv", sep="\t", row.names=F)
write.table(cayley_ct_left_vonEconomo, "results/cayley_ct_left_vonEconomo.csv", sep="\t", row.names=F)

write.table(cayley_ct_right_yeo7, "results/cayley_ct_right_yeo7.csv", sep="\t", row.names=F)
write.table(cayley_ct_right_yeo17, "results/cayley_ct_right_yeo17.csv", sep="\t", row.names=F)
write.table(cayley_ct_right_vonEconomo, "results/cayley_ct_right_vonEconomo.csv", sep="\t", row.names=F)

#######################################

colnames(cayley_ct_left_yeo7) <- yeo7_left
colnames(cayley_ct_right_yeo7) <- yeo7_right

gf_cayley <- cbind(gf, cayley_ct_left_yeo7, cayley_ct_right_yeo7)
gf_cayley_FEP <- subset(gf_cayley, DX=="FEP")

###Comparing all groups against HC
getLms(gf_cayley %>% select(ends_with(c("_left", "_right"))), "Age + Sex + euler_total + DX", gf_cayley) %>%
	filter(!grepl("Unmatched", Variable)) %>%
	write_tsv(., file="results/results_cayley_groupdiff.tsv")

###Comparing only FEP against HC
getLms(subset(gf_cayley, DX=="FEP" | DX=="HC") %>% select(ends_with(c("_left", "_right"))), "Age + Sex + euler_total + DX", subset(gf_cayley, DX=="FEP" | DX=="HC")) %>% 
	filter(!grepl("Unmatched", Variable)) %>%
	write_tsv(., file="results/results_cayley_groupdiff_FEPonly.tsv")

###Comparing raw CT median over networks, then Lms

median_ct_left_yeo7 <- vertexTable(gf_cayley$left_ct) %>% getMedians(., stat_maps_left$yeo7)
median_ct_right_yeo7 <- vertexTable(gf_cayley$right_ct) %>% getMedians(., stat_maps_right$yeo7)

colnames(median_ct_left_yeo7) <- yeo7_left
colnames(median_ct_right_yeo7) <- yeo7_right

getLms(median_ct_left_yeo7, "Age + Sex + euler_total + DX", gf_cayley)
getLms(median_ct_right_yeo7, "Age + Sex + euler_total + DX", gf_cayley)

write_tsv( rbind(getLms(median_ct_left_yeo7, "Age + Sex + euler_total + DX", gf_cayley), getLms(median_ct_right_yeo7, "Age + Sex + euler_total + DX", gf_cayley)) %>% filter(!grepl("Unmatched", Variable)), file="results/results_medianct_groupdiff.tsv")

###Cayley distance vs. symptom correlations: none significant
getCors(gf_cayley_FEP$totalP, gf_cayley_FEP %>% select(ends_with(c("_left", "_right"))) ) %>% filter(pvalue < 0.05)
getCors(gf_cayley_FEP$totalN, gf_cayley_FEP %>% select(ends_with(c("_left", "_right"))) ) %>% filter(pvalue < 0.05)
getCors(gf_cayley_FEP$WksTo50P, gf_cayley_FEP %>% select(ends_with(c("_left", "_right"))) ) %>% filter(pvalue < 0.05)
getCors(gf_cayley_FEP$WksToCGI2, gf_cayley_FEP %>% select(ends_with(c("_left", "_right"))) ) %>% filter(pvalue < 0.05)

###Plotting

results_cayley_groupdiff <- getLms(gf_cayley %>% select(ends_with(c("_left", "_right"))), "Age + Sex + euler_total + DX", gf_cayley) %>% filter(!grepl("Unmatched", Variable))

melted <- melt(results_cayley_groupdiff %>% select(Variable, starts_with("t_DX")), id="Variable", variable.name="network")
melted %<>%  mutate(Hemisphere = ifelse(grepl("_left", Variable), "Left", ifelse(grepl("_right", Variable), "Right", NA)))

melted$Variable <- gsub("_left|_right|_yeo", "", melted$Variable)
melted$Variable <- factor(melted$Variable, levels = unique(melted$Variable))
networks_recode <- c("Visual" = "VIS",
		"Somatomotor" = "SM",
		"DorsalAtt" = "DAN",
		"VentralAtt" = "VAN",
		"Limbic" = "LMB",
		"Frontoparietal" = "FP",
		"Default" = "DMN")
melted %<>% mutate(Variable = recode(Variable, !!!networks_recode))

melted$network <- gsub("t_DX", "", melted$network)
melted$network <- factor(melted$network, levels = c("FEP", "CHR", "3+year"))

plot_groups <- ggplot(melted, aes(x=Variable, y=value, fill=network)) + geom_bar(stat="identity") + facet_grid(Hemisphere ~ network) + scale_fill_Publication2() + theme_Publication() + theme(legend.position="none", axis.text.x=element_text(size=8, angle=0)) + xlab("Networks") + ylab("t-statistic")

ggsave("figures/Plot_cayley_dx_groups.png", plot_groups, width=8, height=5, bg="white")

melted <- melt(gf_cayley %>% select(DX, ends_with(c("left", "right"))), id="DX", variable.name="network")
melted %<>% filter(network=="Somatomotor_left" | network=="DorsalAtt_left")
melted$value %<>% as.numeric()
melted$network <-gsub('Somatomotor_left', 'Left Somatomotor',
		gsub('DorsalAtt_left', 'Left Dorsal Attention', melted$network))

plot_cayley <- ggplot(melted, aes(x=network, y=value, colour=DX)) + geom_boxplot() + geom_point(aes(colour=DX), size=3, position=position_jitterdodge()) + theme_Publication() + facet_wrap( ~ network, scales="free") + xlab("") + ylab("Cayley distance") +     theme(axis.text.x=element_text(hjust=0.35, angle=0), axis.text=element_text(size=14), legend.title=element_blank(), legend.position = "bottom", strip.text.x=element_blank(), axis.title.x=element_blank()) + scale_colour_Publication()

ggsave("figures/Plot_cayley.png", plot_cayley, width=10, height=5, bg="white")

###Rank-rank comparisons of lowest vs. highest distance
###Left dorsal attention

cayley_vis <- gf_cayley %>% filter(DorsalAtt_left<=3353 | DorsalAtt_left >= 3359)
cayley_vis_ranks <- vertexTableRank(cayley_vis$left_ct)
cayley_vis_ranks <- cayley_vis_ranks[,stat_maps_left$yeo7==3] %>% rowRanks()

cayley_vis_ranks <- t(vertexTable(cayley_vis$left_ct))
cayley_vis_ranks <- cayley_vis_ranks[,stat_maps_left$yeo7==3] %>% rowRanks()

qplot(cayley_vis_ranks[which.min(cayley_vis$DorsalAtt_left),], cayley_vis_ranks[which.max(cayley_vis$DorsalAtt_left),], )
qplot(cayley_vis_ranks[which.min(cayley_vis$DorsalAtt_left),], cayley_vis_ranks[cayley_vis$DorsalAtt_left==3353,], )

cor.test(cayley_vis_ranks[which.min(cayley_vis$DorsalAtt_left),], cayley_vis_ranks[which.max(cayley_vis$DorsalAtt_left),], method="spearman")
cor.test(cayley_vis_ranks[which.min(cayley_vis$DorsalAtt_left),], cayley_vis_ranks[cayley_vis$DorsalAtt_left==3353,], method="spearman")

########################################################
###Plotting correlations between networks: medians and cayley

pheatmap::pheatmap(cor(gf_cayley %>% select(ends_with(c("_left", "_right"))) %>% select(!starts_with("Unmatched")) ))
pheatmap::pheatmap(cor(gf_medians_raw_FEP %>% select(ends_with(c("_left", "_right")))) %>% select(!starts_with("Unmatched")) ))
