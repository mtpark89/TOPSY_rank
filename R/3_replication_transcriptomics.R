library(readr)
library(magrittr)
library(dplyr)
library(ggplot2)
library(ggridges)
library(RMINC)
library(matrixStats)
library(Rankcluster)
library(ggpubr)

###5. Imaging-transcriptomics & comparison to raw CT-based testing

ct_compare_ahba <- data.frame()
for (dataset in unique(combined_df$Dataset)){
	
	data <- subset(combined_df, Dataset==dataset)
	data %<>% filter(DX=="HC" | DX=="FEP" | DX=="Scz")
	data$DX %<>% droplevels()

	data_split <- split(data, data$DX)

	#Read in files, raw ct and ranked, and for groups.	
	left_ct <- vertexTable(data$left_ct)
	right_ct <- vertexTable(data$right_ct)

	left_ct_rank <- t(vertexTableRank(data$left_ct))
	right_ct_rank <- t(vertexTableRank(data$right_ct))

	left_ct_rank_HC <- t(vertexTableRank(data_split[[1]]$left_ct))
	right_ct_rank_HC <- t(vertexTableRank(data_split[[1]]$right_ct))

	left_ct_rank_Scz <- t(vertexTableRank(data_split[[2]]$left_ct))
	right_ct_rank_Scz <- t(vertexTableRank(data_split[[2]]$right_ct))
	
	#Group means: mean CT across dataset, and ranked mean CT
	group_mean_left <- cbind(rowMeans(left_ct), rank(rowMeans(left_ct), ties.method="first"))
	group_mean_right <- cbind(rowMeans(right_ct), rank(rowMeans(right_ct), ties.method="first"))

	colnames(group_mean_left) <- c("left_ct", "left_ct_rank")
	colnames(group_mean_right) <- c("right_ct", "right_ct_rank")
	
	#Get group-wise median rank of individual ranked CT (which is ranked again), and combine with cortical maps
	group_left_ct_final <- cbind(group_mean_left, rank(rowMedians(left_ct_rank), ties.method="first"), stat_maps_left)
	colnames(group_left_ct_final)[1:3] <- c("mean_ct", "mean_ct_rank", "median_ct_rank")
	
	#Print out plots for each dataset
	plot1 <- ggplot(group_left_ct_final, aes(x=median_ct_rank, y=as.factor(yeo7))) +
		geom_density_ridges(alpha=1, size=1) + stat_density_ridges(quantile_lines=TRUE, quantiles=0.5) + 
		theme_ridges(grid=FALSE, center_axis_labels=TRUE) +
		labs(x = "", y = "Yeo Networks", title = dataset) +
		scale_y_discrete(labels=c("Unmatched", "Visual", "Somatomotor", "DorsalAtt", "VentralAtt", "Limbic", "Frontoparietal", "Default")) +
		theme_minimal() + 
		theme(text=element_text(size=22, family="Helvetica"),
			axis.text.x=element_text(angle=45, hjust=0.35),
			axis.title.x=element_blank(),
			axis.text=element_text(size=16),
			plot.title = element_text(size=20, face="bold"))

	plot2 <- ggplot(group_left_ct_final, aes(x=median_ct_rank, y=as.factor(vonEconomo))) +
		geom_density_ridges(alpha=1, size=1) + stat_density_ridges(quantile_lines=TRUE, quantiles=0.5) + 
		theme_ridges(grid=FALSE, center_axis_labels=TRUE) +
		labs(x = "", y = "von Economo class") +
		scale_y_discrete(labels=c("Unmatched", "Visual", "Somatomotor", "DorsalAtt", "VentralAtt", "Limbic", "Frontoparietal", "Default")) +
		theme_minimal() + 
		theme(text=element_text(size=22, family="Helvetica"),
			axis.text.x=element_text(angle=45, hjust=0.35),
			axis.title.x=element_blank(),
			axis.text=element_text(size=16),
			plot.title = element_text(size=20, face="bold"))

	ggsave(paste("results/Plot_networks_ranks_",dataset,".png",sep=""), ggarrange(plot1, plot2, ncol=2), width=10, height=5, bg="white")

	#Model comparisons: raw vs. ranked
	
	left_ct_compare <- compare_models(data, "left_ct", "DX")
	right_ct_compare <- compare_models(data, "right_ct", "DX")
	combined_compare <- rbind(left_ct_compare, right_ct_compare)

	assign(paste0("left_ct_compare_", dataset), left_ct_compare)
	assign(paste0("right_ct_compare_", dataset), right_ct_compare)

	#t.tests for combined, and left/right, and add results to df for all datasets
	t_comb <- t.test(combined_compare$AIC_raw, combined_compare$AIC_ranked)
	t_left <- t.test(left_ct_compare$AIC_raw, left_ct_compare$AIC_ranked)
	t_right <- t.test(right_ct_compare$AIC_raw, right_ct_compare$AIC_ranked)

	t_aic <- cbind(dataset, parse_ttest(t_comb), parse_ttest(t_left), parse_ttest(t_right))
	colnames(t_aic) <- c("dataset", "t_combined", "p_combined", "t_left", "p_left", "t_right", "p_right")
		
	datasets_t_aic <- rbind(datasets_t_aic, t_aic)

	#Group differences in average rank using permutation testing, and write out results
	assign(paste0("left_ct_rank_diff_", dataset), average_rank_diff_permute(data_split[[1]]$left_ct, data_split[[2]]$left_ct, 10000))
	assign(paste0("right_ct_rank_diff_", dataset), average_rank_diff_permute(data_split[[1]]$right_ct, data_split[[2]]$right_ct, 10000))

	writeVertex(get(paste0("left_ct_rank_diff_", dataset)), paste0("results/left_ct_rank_diff_",dataset,".vertstats"))
	writeVertex(get(paste0("right_ct_rank_diff_", dataset)), paste0("results/right_ct_rank_diff_",dataset,".vertstats"))

	#Plotting mean rank differences
	left_ct_rank_diff <- get(paste0("left_ct_rank_diff_", dataset))
	right_ct_rank_diff <- get(paste0("right_ct_rank_diff_", dataset))
	
	left_ct_rank_diff$hemi <- "Left CT"
	right_ct_rank_diff$hemi <- "Right CT"
	ct_rank_diff <- rbind(left_ct_rank_diff, right_ct_rank_diff)

	plot_histo <- 
	ggplot(ct_rank_diff, aes(x=diff, group=hemi, fill=hemi)) + 
		geom_density(alpha=0.5, size=1) + 
		labs(x = "Mean rank difference", y = "Density", title = dataset) +
		theme_minimal() + 
		theme(text=element_text(size=22, family="Helvetica"), 
			axis.text.x=element_text(angle=45), 
			axis.text=element_text(size=16),
			legend.position="none", 
			plot.title=element_text(hjust=0.5)) +
		geom_vline(xintercept=sd(ct_rank_diff$diff)*2, linetype="dashed", colour="darkgray", size=2) + 
		geom_vline(xintercept=-sd(ct_rank_diff$diff)*2, linetype="dashed", colour="darkgray", size=2)

	ggsave(paste0("results/Plot_histo_mean_rankdiff_",dataset,".png"), plot_histo, width=5, height=5, bg="white")
	
	###DONE###
}


###AHBA: ranked CT
left_ct_rank_diff_ahba <- ahba_CIVET(left_ct_rank_diff$diff)
ahba_CIVET_write("left_ct_rank_diff_ahba", left_ct_rank_diff_ahba, "left_ct_rank_diff_ahba/")

###AHBA: raw CT linear model
vs <- vertexLm(left_ct ~ DX, data=subset(hcp, DX=="FEP" | DX=="HC"))

vs <- vertexLm(left_ct ~ DX, data=subset(combined_df_TOPSY, Dataset=="COBRE"))
vs %<>% data.frame()

left_ct_raw_lm_ahba <- ahba_CIVET(vs$tvalue.DXScz)

left_ct_raw_lm_ahba <- ahba_CIVET(vs$tvalue.DXFEP)
ahba_CIVET_write("left_ct_raw_lm_ahba", left_ct_raw_lm_ahba, "left_ct_raw_lm_ahba/")

###AHBA: raw CT using mean CT & subtraction

left_ct_raw_meandiff <- average_vertex_diff(hcp_HC$left_ct, hcp_FEP$left_ct)
writeVertex(left_ct_raw_meandiff, file="results/left_ct_raw_meandiff.vertstats")

left_ct_raw_meandiff_ahba <- ahba_CIVET(left_ct_raw_meandiff$diff)
ahba_CIVET_write("left_ct_raw_meandiff_ahba", left_ct_raw_meandiff_ahba, "left_ct_raw_meandiff_ahba/")
