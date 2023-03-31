library(readr)
library(magrittr)
library(dplyr)
library(ggplot2)
library(ggridges)
library(RMINC)
library(matrixStats)
library(Rankcluster)
library(ggpubr)

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

###########################################################################################################
######ANALYSIS#####

###1. Visual inspection of cortical thickness & surface area ranks

left_ct_rank <- t(vertexTableRank(gf$left_ct))
right_ct_rank <- t(vertexTableRank(gf$right_ct))

left_ct <- vertexTable(gf$left_ct)
right_ct <- vertexTable(gf$right_ct)

###Sample individuals
sample_subject_left <- cbind(left_ct[,1], left_ct_rank[,1], left_sa[,1], left_sa_rank[,1], left_vol[,1], left_vol_rank[,1], left_mc[,1], left_mc_rank[,1])
sample_subject_right <- cbind(right_ct[,1], right_ct_rank[,1], right_sa[,1], right_sa_rank[,1], right_vol[,1], right_vol_rank[,1], right_mc[,1], right_mc_rank[,1])

colnames(sample_subject_left) <- c("left_ct", "left_ct_rank", "left_sa", "left_sa_rank", "left_vol", "left_vol_rank", "left_mc", "left_mc_rank")
colnames(sample_subject_right) <- c("right_ct", "right_ct_rank", "right_sa", "right_sa_rank", "right_vol", "right_vol_rank", "right_mc", "right_mc_rank")

sample_subject_left %<>% cbind(., stat_maps_left)
sample_subject_right %<>% cbind(., stat_maps_right)

writeVertex(sample_subject_left, "results/sample_subject_left.vertstats")
writeVertex(sample_subject_right, "results/sample_subject_right.vertstats")

writeVertex(left_ct_rank, "results/left_ct_rank.vertstats")
writeVertex(right_ct_rank, "results/right_ct_rank.vertstats")

writeVertex(cbind(stat_maps, rbind(left_ct_rank, right_ct_rank)), "results/ct_rank.vertstats")
writeVertex(cbind(stat_maps, rbind(left_ct, right_ct)), "results/ct_raw.vertstats")

###Group mean CT & ranked
group_mean_left <- cbind(rowMeans(left_ct), rank(rowMeans(left_ct), ties.method="first"))
group_mean_right <- cbind(rowMeans(right_ct), rank(rowMeans(right_ct), ties.method="first"))

colnames(group_mean_left) <- c("left_ct", "left_ct_rank")
colnames(group_mean_right) <- c("right_ct", "right_ct_rank")

group_mean_left %<>% cbind(., stat_maps_left)
group_mean_right %<>% cbind(., stat_maps_right)

writeVertex(group_mean_left, "results/group_mean_left.vertstats")
writeVertex(group_mean_right, "results/group_mean_right.vertstats")

###All combined: sample subject & group mean, CT only for final figures
sample_subject_left <- cbind(left_ct[,1], left_ct_rank[,1])
sample_subject_right <- cbind(right_ct[,1], right_ct_rank[,1])

colnames(sample_subject_left) <- c("sub_left_ct", "sub_left_ct_rank")
colnames(sample_subject_right) <- c("sub_right_ct", "sub_right_ct_rank")

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

group_left_ct_final <-
	cbind(group_mean_left, rank(rowMedians(left_ct_rank), ties.method="first"), stat_maps_left)
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

left_sa_compare <- compare_models(gf_compare, "left_sa", "DX")
right_sa_compare <- compare_models(gf_compare, "right_sa", "DX")

writeVertex(left_ct_compare, file="left_ct_compare.vertstats")
writeVertex(right_ct_compare, file="right_ct_compare.vertstats")

writeVertex(left_sa_compare, file="left_sa_compare.vertstats")
writeVertex(right_sa_compare, file="right_sa_compare.vertstats")

###4. Group differences in average rank

left_ct_rank_diff <- average_rank_diff(gf_HC$left_ct, gf_FEP$left_ct)
right_ct_rank_diff <- average_rank_diff(gf_HC$right_ct, gf_FEP$right_ct)

#Using permutation testing
left_ct_rank_diff <- average_rank_diff_permute(gf_HC$left_ct, gf_FEP$left_ct, 10000)
right_ct_rank_diff <- average_rank_diff_permute(gf_HC$right_ct, gf_FEP$right_ct, 10000)

writeVertex(left_ct_rank_diff, "results/left_ct_rank_diff.vertstats")
writeVertex(right_ct_rank_diff, "results/right_ct_rank_diff.vertstats")

ggqqplot(left_ct_rank_diff$diff)

left_ct_rank_diff$hemi <- "Left CT"
right_ct_rank_diff$hemi <- "Right CT"

#Plotting left & right hemisphere rank differences histogram
ct_rank_diff <- rbind(left_ct_rank_diff, right_ct_rank_diff)

sd(ct_rank_diff$diff) #SD is 1648.211, *2 is 3296.

plot_histo <- ggplot(ct_rank_diff, aes(x=diff, group=hemi, fill=hemi)) + geom_density(alpha=0.5, size=1) + 
	xlab("Mean rank difference") + ylab("Density") +  theme_minimal() + scale_fill_discrete() +
	theme(text=element_text(size=22, family="Helvetica"), axis.text.x=element_text(angle=45), axis.text=element_text(size=16), legend.title=element_blank()) +
	geom_vline(xintercept=3296, linetype="dashed", colour="darkgray", size=2) + geom_vline(xintercept=-3296, linetype="dashed", colour="darkgray", size=2)

ggsave("Plot_histo_mean_rankdiff.png", plot_histo, width=5, height=5, bg="white")

###5. Imaging-transcriptomics & comparison to raw CT-based testing

###Compare number of significant genes
###Enrichment using significant genes only vs. 10%, and compare
###DAVID: write out Entrez and

###AHBA: ranked CT
left_ct_rank_diff_ahba <- ahba_CIVET(left_ct_rank_diff$diff)
ahba_CIVET_write("left_ct_rank_diff_ahba", left_ct_rank_diff_ahba, "left_ct_rank_diff_ahba/")

###AHBA: raw CT linear model
vs <- vertexLm(left_ct ~ DX, data=subset(gf, DX=="FEP" | DX=="HC"))
vs %<>% data.frame()

left_ct_raw_lm_ahba <- ahba_CIVET(vs$tvalue.DXFEP)
ahba_CIVET_write("left_ct_raw_lm_ahba", left_ct_raw_lm_ahba, "left_ct_raw_lm_ahba/")

###AHBA: raw CT using mean CT & subtraction

left_ct_raw_meandiff <- average_vertex_diff(gf_HC$left_ct, gf_FEP$left_ct)
writeVertex(left_ct_raw_meandiff, file="results/left_ct_raw_meandiff.vertstats")

left_ct_raw_meandiff_ahba <- ahba_CIVET(left_ct_raw_meandiff$diff)
ahba_CIVET_write("left_ct_raw_meandiff_ahba", left_ct_raw_meandiff_ahba, "left_ct_raw_meandiff_ahba/")

###Toppgene results per decile (1 and 10): DisGeNEt curated

toppgene_rank_dec10 <- read_tsv("left_ct_rank_diff_ahba/ToppGene_rank_decile10_tstat.txt")
toppgene_rawlm_dec10 <- read_tsv("left_ct_raw_lm_ahba/ToppGene_rawlm_decile10_tstat.txt")
toppgene_meandiff_dec10 <- read_tsv("left_ct_raw_meandiff_ahba/ToppGene_meandiff_decile10_tstat.txt")

toppgene_rank_dec1 <- read_tsv("left_ct_rank_diff_ahba/ToppGene_rank_decile1_tstat.txt")
toppgene_rawlm_dec1 <- read_tsv("left_ct_raw_lm_ahba/ToppGene_rawlm_decile1_tstat.txt")
toppgene_meandiff_dec1 <- read_tsv("left_ct_raw_meandiff_ahba/ToppGene_meandiff_decile1_tstat.txt")

toppgene_rank_dec1$rank <- 1:nrow(toppgene_rank_dec1)
toppgene_rawlm_dec1$rank <- 1:nrow(toppgene_rawlm_dec1)
toppgene_meandiff_dec1$rank <- 1:nrow(toppgene_meandiff_dec1)

toppgene_rank_dec10$rank <- 1:nrow(toppgene_rank_dec10)
toppgene_rawlm_dec10$rank <- 1:nrow(toppgene_rawlm_dec10)
toppgene_meandiff_dec10$rank <- 1:nrow(toppgene_meandiff_dec10)

toppgene_rank_dec1$Method <- "Rank"
toppgene_rawlm_dec1$Method <- "CT-linear model"
toppgene_meandiff_dec1$Method <- "CT-difference"

toppgene_rank_dec10$Method <- "Rank"
toppgene_rawlm_dec10$Method <- "CT-linear model"
toppgene_meandiff_dec10$Method <- "CT-difference"

toppgene_dec1 <- rbind(toppgene_rank_dec1, toppgene_rawlm_dec1)
toppgene_dec10 <- rbind(toppgene_rank_dec10, toppgene_rawlm_dec10)

toppgene_dec1 %<>% mutate(scz = case_when(Name=="Schizophrenia" ~ "Schizophrenia"))
toppgene_dec10 %<>% mutate(scz = case_when(Name=="Schizophrenia" ~ "Schizophrenia"))

toppgene_dec1 %<>% mutate(pointsize = case_when(Name=="Schizophrenia" ~ 2, Name!="Schizophrenia" ~ 1))
toppgene_dec10 %<>% mutate(pointsize = case_when(Name=="Schizophrenia" ~ 2, Name!="Schizophrenia" ~ 1))

plot_dec1 <- ggplot(subset(toppgene_dec1, rank < 50), aes(x=rank, y=-log10(`p-value`), group=Method, colour=Method)) + 
	geom_point(size=2) + theme_minimal() +
	xlab("Rank") + ylab("-log10(p value)") +
	geom_hline(yintercept=-log10(0.05/2600), linetype="dashed", colour="darkgray", size=2) +
	theme(text=element_text(size=22, family="Helvetica"), axis.text.x=element_text(hjust=0.35), axis.text=element_text(size=16))

plot_dec10 <- ggplot(subset(toppgene_dec10, rank < 50), aes(x=rank, y=-log10(`p-value`), group=Method, colour=Method)) + 
	geom_point(size=2) + theme_minimal() +
	xlab("Rank") + ylab("-log10(p value)") +
	geom_hline(yintercept=-log10(0.05/2600), linetype="dashed", colour="darkgray", size=2) +
	theme(text=element_text(size=22, family="Helvetica"), axis.text.x=element_text(hjust=0.35), axis.text=element_text(size=16))

write.table(subset(toppgene_dec1, Name=="Schizophrenia"), file="figures/toppgene_dec1.csv", sep="\t", row.names=F)
write.table(subset(toppgene_dec10, Name=="Schizophrenia"), file="figures/toppgene_dec10.csv", sep="\t", row.names=F)

write.table(toppgene_dec1, file="results/toppgene_dec1.csv", sep="\t", row.names=F)
write.table(toppgene_dec10, file="results/toppgene_dec10.csv", sep="\t", row.names=F)

###2600 annotations for Bonferroni

ggsave("figures/Plot_toppgene_deciles_comparison.png", ggarrange(plot_dec1, plot_dec10, ncol=2, common.legend = TRUE, legend="bottom"), width=10, height=5, bg="white")

###6. Median rank testing

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

getCors(gf_medians_FEP$totalP, gf_medians_FEP[(ncol(gf_medians_FEP)-15):ncol(gf_medians_FEP)])
getCors(gf_medians_FEP$totalN, gf_medians_FEP[(ncol(gf_medians_FEP)-15):ncol(gf_medians_FEP)])
getCors(gf_medians_FEP$WksTo50P, gf_medians_FEP[(ncol(gf_medians_FEP)-15):ncol(gf_medians_FEP)])
getCors(gf_medians_FEP$WksToCGI2, gf_medians_FEP[(ncol(gf_medians_FEP)-15):ncol(gf_medians_FEP)])

summary(lm(gf_medians_FEP$DorsalAtt_left ~ Age + Sex + euler_total + WksToCGI2, data=gf_medians_FEP))
summary(lm(gf_medians_FEP$DorsalAtt_right ~ Age + Sex + euler_total + WksToCGI2, data=gf_medians_FEP))

melted <- melt(gf_medians_FEP %>% select(DorsalAtt_left, DorsalAtt_right, WksToCGI2), id="WksToCGI2", variable.name="hemisphere")
melted$hemisphere <-gsub('DorsalAtt_left', 'Left',
		gsub('DorsalAtt_right', 'Right', melted$hemisphere))

plot_cgi <-ggplot(melted, aes(x=value, y=WksToCGI2, group=hemisphere, colour=hemisphere)) + 
    geom_point(size=2) + theme_minimal() + geom_smooth(method=lm, size=2) +
    xlab("Median rank: Dorsal Attention network") + ylab("Weeks to reach CGI-S score 2") +
    theme(text=element_text(size=18, family="Helvetica"), axis.text.x=element_text(hjust=0.35), axis.text=element_text(size=14), legend.title=element_blank(), legend.position = "bottom")

ggsave("Plot_cgi_dorsal.png", plot_cgi, width=6, height=5, bg="white")

summary(lm(gf_medians_FEP$Default_left ~ Age + Sex + euler_total + totalP, data=gf_medians_FEP))

plot_totalP <- ggplot(gf_medians_FEP, aes(x=Default_left, y=totalP)) + 
	geom_point(size=2) + theme_minimal() + geom_smooth(method=lm, size=2) +
	xlab("Median rank: Default network") + ylab("Total positive symptoms") +
	theme(text=element_text(size=18, family="Helvetica"), axis.text.x=element_text(hjust=0.35), axis.text=element_text(size=14), legend.title=element_blank())

ggsave("Plot_totalP.png", plot_totalP, width=6, height=5, bg="white")

####################################
###Comparison: raw cortical thickness & network medians

left_ct <- vertexTable(gf$left_ct)
right_ct <- vertexTable(gf$right_ct)

left_ct_yeo <- getMedians(left_ct, stat_maps_left$yeo7)
right_ct_yeo <- getMedians(right_ct, stat_maps_right$yeo7)

gf_medians_raw <- cbind(gf, left_ct_yeo, right_ct_yeo)
gf_medians_raw_FEP <- subset(gf_medians_raw, DX=="FEP")

getCors(gf_medians_raw_FEP$totalP, gf_medians_raw_FEP[(ncol(gf_medians_raw_FEP)-15):ncol(gf_medians_raw_FEP)])
getCors(gf_medians_raw_FEP$totalN, gf_medians_raw_FEP[(ncol(gf_medians_raw_FEP)-15):ncol(gf_medians_raw_FEP)])
getCors(gf_medians_raw_FEP$WksTo50P, gf_medians_raw_FEP[(ncol(gf_medians_raw_FEP)-15):ncol(gf_medians_raw_FEP)])
getCors(gf_medians_raw_FEP$WksToCGI2, gf_medians_raw_FEP[(ncol(gf_medians_raw_FEP)-15):ncol(gf_medians_raw_FEP)])

###Plotting for comparison
gf_medians_FEP$Method <- "CT-ranked"
gf_medians_raw_FEP$Method <- "CT"

gf_medians_combined <- rbind(gf_medians_FEP, gf_medians_raw_FEP)

#####################################################
###7. Rank-based normative modelling: using median rank of HCs to get cayley distances between individual & median network ranks

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

cayley_ct_left_yeo7 %<>% append_colnames(., "_left")
cayley_ct_right_yeo7 %<>% append_colnames(., "_right")

gf_cayley <- cbind(gf, cayley_ct_left_yeo7, cayley_ct_right_yeo7)
gf_cayley_FEP <- subset(gf_cayley, DX=="FEP")

summary(lm(cayley_1_left ~ Age + Sex + DX, data=gf_cayley))
summary(lm(cayley_2_left ~ Age + Sex + DX, data=gf_cayley))
summary(lm(cayley_3_left ~ Age + Sex + DX, data=gf_cayley))
summary(lm(cayley_4_left ~ Age + Sex + DX, data=gf_cayley))
summary(lm(cayley_5_left ~ Age + Sex + DX, data=gf_cayley))
summary(lm(cayley_6_left ~ Age + Sex + DX, data=gf_cayley))
summary(lm(cayley_7_left ~ Age + Sex + DX, data=gf_cayley))

summary(lm(cayley_1_left ~ Age + Sex + euler_total + DX, data=gf_cayley))
summary(lm(cayley_2_left ~ Age + Sex + euler_total + DX, data=gf_cayley))
summary(lm(cayley_3_left ~ Age + Sex + euler_total + DX, data=gf_cayley))
summary(lm(cayley_4_left ~ Age + Sex + euler_total + DX, data=gf_cayley))
summary(lm(cayley_5_left ~ Age + Sex + euler_total + DX, data=gf_cayley))
summary(lm(cayley_6_left ~ Age + Sex + euler_total + DX, data=gf_cayley))
summary(lm(cayley_7_left ~ Age + Sex + euler_total + DX, data=gf_cayley))

summary(lm(cayley_1_left ~ Age + Sex + euler_total + DX, data=subset(gf_cayley, DX=="FEP" | DX=="HC")))
summary(lm(cayley_2_left ~ Age + Sex + euler_total + DX, data=subset(gf_cayley, DX=="FEP" | DX=="HC")))
summary(lm(cayley_3_left ~ Age + Sex + euler_total + DX, data=subset(gf_cayley, DX=="FEP" | DX=="HC")))
summary(lm(cayley_4_left ~ Age + Sex + euler_total + DX, data=subset(gf_cayley, DX=="FEP" | DX=="HC")))
summary(lm(cayley_5_left ~ Age + Sex + euler_total + DX, data=subset(gf_cayley, DX=="FEP" | DX=="HC")))
summary(lm(cayley_6_left ~ Age + Sex + euler_total + DX, data=subset(gf_cayley, DX=="FEP" | DX=="HC")))
summary(lm(cayley_7_left ~ Age + Sex + euler_total + DX, data=subset(gf_cayley, DX=="FEP" | DX=="HC")))

summary(lm(cayley_1_right ~ Age + Sex + DX, data=gf_cayley))
summary(lm(cayley_2_right ~ Age + Sex + DX, data=gf_cayley))
summary(lm(cayley_3_right ~ Age + Sex + DX, data=gf_cayley))
summary(lm(cayley_4_right ~ Age + Sex + DX, data=gf_cayley))
summary(lm(cayley_5_right ~ Age + Sex + DX, data=gf_cayley))
summary(lm(cayley_6_right ~ Age + Sex + DX, data=gf_cayley))
summary(lm(cayley_7_right ~ Age + Sex + DX, data=gf_cayley))

summary(lm(cayley_1_right ~ Age + Sex + euler_total + DX, data=gf_cayley))
summary(lm(cayley_2_right ~ Age + Sex + euler_total + DX, data=gf_cayley))
summary(lm(cayley_3_right ~ Age + Sex + euler_total + DX, data=gf_cayley))
summary(lm(cayley_4_right ~ Age + Sex + euler_total + DX, data=gf_cayley))
summary(lm(cayley_5_right ~ Age + Sex + euler_total + DX, data=gf_cayley))
summary(lm(cayley_6_right ~ Age + Sex + euler_total + DX, data=gf_cayley))
summary(lm(cayley_7_right ~ Age + Sex + euler_total + DX, data=gf_cayley))

getCors(gf_cayley_FEP$totalP, gf_cayley_FEP[(ncol(gf_cayley_FEP)-15):ncol(gf_cayley_FEP)])
getCors(gf_cayley_FEP$totalN, gf_cayley_FEP[(ncol(gf_cayley_FEP)-15):ncol(gf_cayley_FEP)])
getCors(gf_cayley_FEP$WksTo50P, gf_cayley_FEP[(ncol(gf_cayley_FEP)-15):ncol(gf_cayley_FEP)])
getCors(gf_cayley_FEP$WksToCGI2, gf_cayley_FEP[(ncol(gf_cayley_FEP)-15):ncol(gf_cayley_FEP)])

summary(lm(test$`cayley_3` ~ Age + Sex + euler_total + DX, data=test))

colnames(left_ct_yeo) <-  c("Unmatched_left", "Visual_left", "Somatomotor_left", "DorsalAtt_left", "VentralAtt_left", "Limbic_yeo_left", "Frontoparietal_left", "Default_left")
colnames(right_ct_yeo) <-  c("Unmatched_right", "Visual_right", "Somatomotor_right", "DorsalAtt_right", "VentralAtt_right", "Limbic_yeo_right", "Frontoparietal_right", "Default_right")

###Plotting

melted <- melt(gf_cayley %>% select(DX, contains(c("cayley", "left"))), id="DX", variable.name="network")
melted %<>% filter(network=="cayley_2_left" | network=="cayley_3_left")
melted$value %<>% as.numeric()
melted$network <-gsub('cayley_2_left', 'Left Somatomotor',
		gsub('cayley_3_left', 'Left Dorsal Attention', melted$network))

plot_cayley <- ggplot(melted, aes(x=network, y=value)) + geom_boxplot(aes(colour=DX)) + geom_point(aes(colour=DX), size=3, position=position_jitterdodge()) + theme_minimal() + facet_wrap( ~ network, scales="free") +
	xlab("") + ylab("Cayley distance") +
	theme(text=element_text(size=18, family="Helvetica"), axis.text.x=element_text(hjust=0.35), axis.text=element_text(size=14), legend.title=element_blank(), legend.position = "bottom", strip.text.x=element_blank(), axis.title.x=element_blank())

ggsave("Plot_cayley.png", plot_cayley, width=8, height=5, bg="white")

########################################################
































####AHBA

expressionMatrix_MNI <- read_tsv("/projects/TOPSY_NBM_Allen/data/processed/Cortical_expressionMatrix_qc_15120genes.tsv")
sampleAnnot_MNI <- read_tsv("/projects/TOPSY_NBM_Allen/data/processed/Cortical_expressionMatrix_annotations.tsv")
qc_pass <- read_tsv("/projects/TOPSY_NBM_Allen/data/processed/QCpass_genes_15120.tsv")

identical(colnames(expressionMatrix_MNI[2:ncol(expressionMatrix_MNI)]), sampleAnnot_MNI$uniqueID)
identical(expressionMatrix_MNI$probe_name, qc_pass$probe_name)

rownames(expressionMatrix_MNI) <- qc_pass$gene_symbol
expressionMatrix_MNI %<>% select(-probe_name)

###Correlating imaging to genetic
vs <- as.data.frame(vs)
vs$Vertex <- 1:nrow(vs)

#cortical <- read_tsv(file = "data/Left_cortical_cors.csv")
sampleAnnot_MNI <- left_join(sampleAnnot_MNI, cortical, by="Vertex")

test<- as.data.frame(cor(t(expressionMatrix_MNI), sampleAnnot_MNI$`tvalue-Glutamate`))
test$Gene <- rownames(test)
test %<>% arrange(desc(V1))
write.table(test[1:1512,] %>% select(Gene), file="test1.txt", sep="\t", row.names=F, col.names=F, quote=FALSE)
write.table(test[(15120-1512):15120,] %>% select(Gene), file="test2.txt", sep="\t", row.names=F, col.names=F, quote=FALSE)

library(lmerTest)

ahba_mat <- matrix(data=NA, nrow=nrow(expressionMatrix_MNI), ncol=3)

for (i in 1:nrow(expressionMatrix_MNI)){
	model <- summary(lmer(sampleAnnot_MNI$NBM_cor ~ t(expressionMatrix_MNI[i,]) + (1|sampleAnnot_MNI$donorID), verbose=0))

	ahba_mat[i,1] <- rownames(expressionMatrix_MNI)[i]
	ahba_mat[i,2:3] <- model$coefficients[2,4:5] ###T-statistic and p-value
}


for (i in 1:nrow(expressionMatrix_MNI)){
	model <- summary(lmer(sampleAnnot_MNI$NBM_cor ~ t(expressionMatrix_MNI[i,]) + (1|sampleAnnot_MNI$donorID), verbose=0))

	ahba_mat[i,1] <- rownames(expressionMatrix_MNI)[i]
	ahba_mat[i,2:3] <- model$coefficients[2,4:5] ###T-statistic and p-value
}

ahba_mat2 <- matrix(data=NA, nrow=nrow(expressionMatrix_MNI), ncol=3)

for (i in 1:nrow(expressionMatrix_MNI)){
	model <- summary(lmer(sampleAnnot_MNI$Ch123_cor ~ t(expressionMatrix_MNI[i,]) + (1|sampleAnnot_MNI$donorID), verbose=0))
	ahba_mat2[i,1] <- rownames(expressionMatrix_MNI)[i]
	ahba_mat2[i,2:3] <- model$coefficients[2,4:5] ###T-statistic and p-value
}

colnames(ahba_mat) <- c("Gene", "tstat", "pval")
colnames(ahba_mat2) <- c("Gene", "tstat", "pval")

write.table(ahba_mat, file="results/AHBA_NBM-cortical_backup_20210409.csv", sep=",", row.names=F)
write.table(ahba_mat2, file="results/AHBA_DB-cortical_backup_20210409.csv", sep=",", row.names=F)


###Correlate with gandal transcriptome results?

gandal <- read_tsv(file="../TOPSY_NBM_Allen/data/Science_2018.csv")

test <- left_join(left_ct_rank_diff_ahba, gandal, by=c("Gene"="gene_name") )



###DKT
dkt_labels <- read.csv('/projects/Templates/CIVET_resources/icbm/DKT/DKTatlas40.labels', sep=" ", header=FALSE)
dkt_left <- read.csv('/projects/Templates/CIVET_resources/icbm/DKT/icbm_avg_mid_mc_dkt40_left_40962.txt', sep=" ", header=FALSE)
dkt_right <- read.csv('/projects/Templates/CIVET_resources/icbm/DKT/icbm_avg_mid_mc_dkt40_right_40962.txt', sep=" ", header=FALSE)
colnames(dkt_labels)[2]<- "subregion"

results_dkt <- join(dkt_labels, stats_1)

results_dkt <- join(dkt_right, results_dkt)

results_dkt$subregion<- NULL

mni.write.vertex.stats(results_dkt, file="DKT_results_right_20200926.vertstats", headers=TRUE)


results_dkt <- join(dkt_labels, stats_1)

results_dkt <- join(dkt_left, results_dkt)

results_dkt$subregion<- NULL

mni.write.vertex.stats(results_dkt, file="DKT_results_left_20200926.vertstats", headers=TRUE)


