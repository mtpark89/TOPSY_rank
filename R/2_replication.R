library(readr)
library(magrittr)
library(dplyr)
library(ggplot2)
library(ggridges)
library(RMINC)
library(matrixStats)
library(Rankcluster)
library(ggpubr)

source("../TOPSY_rank/1_functions.R")
###

path <- "/projects/Datasets_phenotypic/Datasets_scz/Final_2023/"

# Get the list of all TSV files in the folder
tsv_files <- list.files(path, pattern = "\\.tsv$", full.names = TRUE)

# Loop through each TSV file
df_list <- list()
for (file in tsv_files) {
  # Read in the TSV file
  df <- read_tsv(file)
  
  # Extract the unique value from the "Dataset" column
  VARIABLE <- unique(df$Dataset)

  df$left_ct <- paste("/projects/Datasets/CIVET_scz/thickness/",df$CIVET_ID,"_native_rms_rsl_tlaplace_30mm_left.txt", sep="") 
  df$right_ct <- paste("/projects/Datasets/CIVET_scz/thickness/",df$CIVET_ID,"_native_rms_rsl_tlaplace_30mm_right.txt", sep="") 
  
  df$left_vol <- paste("/projects/Datasets/CIVET_scz/surface/",df$CIVET_ID,"_surface_rsl_left_native_volume_40mm.txt", sep="")
  df$right_vol <- paste("/projects/Datasets/CIVET_scz/surface/",df$CIVET_ID,"_surface_rsl_right_native_volume_40mm.txt", sep="")

  df$left_mc <- paste("/projects/Datasets/CIVET_scz/thickness/",df$CIVET_ID,"_native_mc_rsl_30mm_mid_left.txt", sep="")
  df$right_mc <- paste("/projects/Datasets/CIVET_scz/thickness/",df$CIVET_ID,"_native_mc_rsl_30mm_mid_right.txt", sep="")

  df$left_sa <- paste("/projects/Datasets/CIVET_scz/surface/",df$CIVET_ID,"_mid_surface_rsl_left_native_area_40mm.txt", sep="")
  df$right_sa <- paste("/projects/Datasets/CIVET_scz/surface/",df$CIVET_ID,"_mid_surface_rsl_right_native_area_40mm.txt", sep="")

  df <- mutate(df, Subject = as.character(Subject))

  # Assign the data frame object with a name containing the unique value
  assign(paste0("df_", VARIABLE), df)
  
  # Add the data frame object to a list
  df_list[[VARIABLE]] <- df
}

combined_df <- bind_rows(df_list)

qc_CIVET <- read_tsv("/projects/Datasets/CIVET_scz/QC_20230326_pass.csv")
combined_df %<>% left_join(., qc_CIVET, by="CIVET_ID")
write_tsv(combined_df, "Datasets_combined_20230327.tsv")

combined_df %<>% filter(QC_CIVET=="Pass")
combined_df$DX <- relevel(factor(combined_df$DX), ref="HC")

###Read parcellations map
stat_maps <- read_delim("../TOPSY_rank/stat_maps.txt", delim=",")
stat_maps_left <- head(stat_maps, n=40962)
stat_maps_right <- tail(stat_maps, n=40962)

###########################################################################################################

datasets_t_aic <- data.frame()
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

###Between-dataset similarity

df_TOPSY <- read_tsv("../TOPSY_rank/TOPSY_subset_20230327.tsv")

  df_TOPSY$left_ct <- paste("/projects/Datasets/CIVET_scz/thickness/",df_TOPSY$CIVET_ID,"_native_rms_rsl_tlaplace_30mm_left.txt", sep="") 
  df_TOPSY$right_ct <- paste("/projects/Datasets/CIVET_scz/thickness/",df_TOPSY$CIVET_ID,"_native_rms_rsl_tlaplace_30mm_right.txt", sep="") 
  
  df_TOPSY$left_vol <- paste("/projects/Datasets/CIVET_scz/surface/",df_TOPSY$CIVET_ID,"_surface_rsl_left_native_volume_40mm.txt", sep="")
  df_TOPSY$right_vol <- paste("/projects/Datasets/CIVET_scz/surface/",df_TOPSY$CIVET_ID,"_surface_rsl_right_native_volume_40mm.txt", sep="")

  df_TOPSY$left_mc <- paste("/projects/Datasets/CIVET_scz/thickness/",df_TOPSY$CIVET_ID,"_native_mc_rsl_30mm_mid_left.txt", sep="")
  df_TOPSY$right_mc <- paste("/projects/Datasets/CIVET_scz/thickness/",df_TOPSY$CIVET_ID,"_native_mc_rsl_30mm_mid_right.txt", sep="")

  df_TOPSY$left_sa <- paste("/projects/Datasets/CIVET_scz/surface/",df_TOPSY$CIVET_ID,"_mid_surface_rsl_left_native_area_40mm.txt", sep="")
  df_TOPSY$right_sa <- paste("/projects/Datasets/CIVET_scz/surface/",df_TOPSY$CIVET_ID,"_mid_surface_rsl_right_native_area_40mm.txt", sep="")

df_TOPSY %<>% filter(QC_CIVET=="Pass" | QC_CIVET=="Partial")
df_TOPSY %<>% filter(DX=="HC" | DX=="FEP" | DX=="Scz")
df_TOPSY$DX <- relevel(factor(df_TOPSY$DX), ref="HC")
df_TOPSY$DX %<>% droplevels()

left_ct_compare_TOPSY <- compare_models(df_TOPSY, "left_ct", "DX")
right_ct_compare_TOPSY <- compare_models(df_TOPSY, "right_ct", "DX")

combined_df_TOPSY <- bind_rows(combined_df, df_TOPSY)

###Collecting p-values of logistic regression, DX ~ CT (raw or ranked)
ct_compare_datasets_raw <- data.frame()
ct_compare_datasets_ranked <- data.frame()

for (dataset in unique(combined_df_TOPSY$Dataset)){

	left_ct_compare <- get(paste0("left_ct_compare_",dataset))
	right_ct_compare <- get(paste0("right_ct_compare_",dataset))

	col_raw <- rbind(left_ct_compare, right_ct_compare)$p_raw %>% data.frame()
	colnames(col_raw) <- paste0(dataset, "_p_raw")

	col_ranked <- rbind(left_ct_compare, right_ct_compare)$p_ranked %>% data.frame()
	colnames(col_ranked) <- paste0(dataset, "_p_ranked")

	if (length(ct_compare_datasets_raw) == 0) {
		ct_compare_datasets_raw <- data.frame(col_raw, check.names = FALSE)
		ct_compare_datasets_ranked <- data.frame(col_ranked, check.names = FALSE)

	} else {
		ct_compare_datasets_raw %<>% cbind(., col_raw)
		ct_compare_datasets_ranked %<>% cbind(., col_ranked)
	}
}

for (dataset in unique(combined_df_TOPSY$Dataset)){

	left_ct_compare <- get(paste0("left_ct_compare_",dataset))
	right_ct_compare <- get(paste0("right_ct_compare_",dataset))

	col_raw <- rbind(left_ct_compare, right_ct_compare)$z_raw %>% data.frame()
	colnames(col_raw) <- paste0(dataset, "_z_raw")

	col_ranked <- rbind(left_ct_compare, right_ct_compare)$z_ranked %>% data.frame()
	colnames(col_ranked) <- paste0(dataset, "_z_ranked")

	if (length(ct_compare_datasets_raw) == 0) {
		ct_compare_datasets_raw <- data.frame(col_raw, check.names = FALSE)
		ct_compare_datasets_ranked <- data.frame(col_ranked, check.names = FALSE)

	} else {
		ct_compare_datasets_raw %<>% cbind(., col_raw)
		ct_compare_datasets_ranked %<>% cbind(., col_ranked)
	}
}


###Thresholding over range, then comparing overlapping vertices with chi-square test
thresholds <- seq(0.01, 0.20, 0.01)

for (i in thresholds){
	apply(ct_compare_datasets_raw, 2, threshold_column, percentile=i)
}

test1 <- apply(ct_compare_datasets_raw, 2, threshold_column, percentile=i) %>% as.data.frame
test2 <- apply(ct_compare_datasets_ranked, 2, threshold_column, percentile=i) %>% as.data.frame
test1 %<>% pairwise_similarity()
test2 %<>% pairwise_similarity()

t.test(test1$Chisq, test2$Chisq)
t.test(test1$Overlap, test2$Overlap)

summary(test1)
summary(test2)

###Cayley distance testing per dataset, network distance ~ Age + Sex + DX

cayley_dx <- data.frame()

for (dataset in unique(combined_df_TOPSY$Dataset)){
	
	data <- subset(combined_df_TOPSY, Dataset==dataset)

	cat("\n")
	cat(dataset)

	data %<>% filter(DX=="HC" | DX=="FEP" | DX=="Scz")
	data$DX <- replace(data$DX, data$DX=="FEP", "Scz")
	data$DX %<>% droplevels()

	if (dataset=="HCPEP") {
		data %<>% filter(DX_2!="Other")
	}

	data_split <- split(data, data$DX)
	
	if (!exists(paste0("cayley_ct_left_yeo7_", dataset))) {
		assign(paste0("cayley_ct_left_yeo7_", dataset), indiv_distance_cayley_sub(data_split$HC$left_ct, data$left_ct, stat_maps_left$yeo7))
	}

	if (!exists(paste0("cayley_ct_right_yeo7_", dataset))) {
		assign(paste0("cayley_ct_right_yeo7_", dataset), indiv_distance_cayley_sub(data_split$HC$right_ct, data$right_ct, stat_maps_right$yeo7))
	}
	

	left_dx <- getLms( get(paste0("cayley_ct_left_yeo7_",dataset)), "Age + Sex + DX", data) %>% mutate(Variable=paste0(Variable, "_left"))
	right_dx <- getLms( get(paste0("cayley_ct_right_yeo7_",dataset)), "Age + Sex + DX", data) %>% mutate(Variable=paste0(Variable, "_right"))

	cayley_dx <- rbind(cayley_dx, cbind(Dataset=dataset, Label="Left Yeo", left_dx))
	cayley_dx <- rbind(cayley_dx, cbind(Dataset=dataset, Label="Right Yeo", right_dx))

	temp_left <- get(paste0("cayley_ct_left_yeo7_",dataset))
	temp_right <- get(paste0("cayley_ct_right_yeo7_",dataset))
	rownames(temp_left) <- data$Subject
	rownames(temp_right) <- data$Subject

	assign(paste0("cayley_ct_left_yeo7_", dataset), temp_left)
	assign(paste0("cayley_ct_right_yeo7_", dataset), temp_right)
}


###Symptom correlations in FEP/SCZ only

cayley_cors <- data.frame()

for (dataset in unique(combined_df_TOPSY$Dataset)){
	
	data <- subset(combined_df_TOPSY, Dataset==dataset)
	data %<>% filter(DX=="HC" | DX=="FEP" | DX=="Scz")
	data$DX <- replace(data$DX, data$DX=="FEP", "Scz")
	data$DX %<>% droplevels()

	if (dataset=="HCPEP") {
		data %<>% filter(DX_2!="Other")
	}


	left <- get(paste0("cayley_ct_left_yeo7_",dataset)) %>% append_colnames(., "_left")
	right <- get(paste0("cayley_ct_right_yeo7_",dataset)) %>% append_colnames(., "_right")

	data <- cbind(data, left, right)

	data %<>% filter(DX!="HC")

	if (dataset %in% c("BrainGluSchi", "COBRE", "HCPEP", "TOPSY")) {
		scores <- data %>% select(starts_with(c("PANSS_")))
	} else {
		scores <- data %>% select(ends_with(c("_global")))
	}


	for(i in 1:ncol(scores)){
		yeo7 <- cbind(Score=names(scores)[i], Dataset=dataset, Label="Yeo7", getCors(scores[[i]], data %>% select(starts_with("cayley_"))))
		#left_von <- cbind(Score=names(scores)[i], dataset, Label="Left von Economo", getCors(scores[[i]], left_ct_vonEconomo))
		#right_von <- cbind(Score=names(scores)[i], dataset, Label="Right von Economo", getCors(scores[[i]], right_ct_vonEconomo))

		cayley_cors <- rbind(cayley_cors, yeo7)
	}	

}

cayley_cors %<>% mutate_at(vars(starts_with(c("Correlation", "pvalue"))), as.numeric)
cayley_cors %>% arrange(., pvalue) %>% filter(!grepl(0, Variable))

######################################
###Repeat with representative ranks of HCs in every dataset, and used as refernece for all subjects in every dataset.
###25, median, 75 percentiles.

reference_hc_left <- data.frame()
reference_hc_right <- data.frame()

for (dataset in unique(combined_df_TOPSY$Dataset)){
	
	data <- subset(combined_df_TOPSY, Dataset==dataset)

	data %<>% filter(DX=="HC")
	
	probs<-c(0.25, 0.5, 0.75)

	left <- vertexTableRank(data$left_ct) %>% colQuantiles(probs=probs) %>% colRanks(ties.method = "first")
	right <- vertexTableRank(data$right_ct) %>% colQuantiles(probs=probs) %>% colRanks(ties.method = "first")
	reference_hc_left <- rbind(reference_hc_left, left)
	reference_hc_right <- rbind(reference_hc_right, right)
}

cayley_ct_left_yeo7_quant <- indiv_distance_cayley_sub(reference_hc_left, combined_df_TOPSY$left_ct, stat_maps_left$yeo7)
cayley_ct_right_yeo7_quant <- indiv_distance_cayley_sub(reference_hc_right, combined_df_TOPSY$right_ct, stat_maps_right$yeo7)

cayley_ct_left_yeo17_quant <- indiv_distance_cayley_sub(reference_hc_left, combined_df_TOPSY$left_ct, stat_maps_left$yeo17)
cayley_ct_right_yeo17_quant <- indiv_distance_cayley_sub(reference_hc_right, combined_df_TOPSY$right_ct, stat_maps_right$yeo17)

cayley_ct_left_vonE_quant <- indiv_distance_cayley_sub(reference_hc_left, combined_df_TOPSY$left_ct, stat_maps_left$vonEconomo)
cayley_ct_right_vonE_quant <- indiv_distance_cayley_sub(reference_hc_right, combined_df_TOPSY$right_ct, stat_maps_right$vonEconomo)


###
combined_df_TOPSY_quant <- cbind(combined_df_TOPSY, cayley_ct_left_yeo7_quant %>% append_colnames(., "_yeo7_left"), cayley_ct_right_yeo7_quant %>% append_colnames(., "_yeo7_right"))

###Age correlations?
getCors(combined_df_TOPSY_quant$Age, combined_df_TOPSY_quant %>% select(starts_with("cayley_")))

###Cayley distance (using quatile ranks across all datasets) testing per dataset, network distance ~ Age + Sex + DX
cayley_dx_quant <- data.frame()

for (dataset in unique(combined_df_TOPSY_quant$Dataset)){
	
	data <- subset(combined_df_TOPSY_quant, Dataset==dataset)

	cat("\n")
	cat(dataset)

	data %<>% filter(DX=="HC" | DX=="FEP" | DX=="Scz")
	data$DX <- replace(data$DX, data$DX=="FEP", "Scz")
	data$DX %<>% droplevels()

	if (dataset=="HCPEP") {
		data %<>% filter(DX_2!="Other")
	}

	left_dx <- getLms(data %>% select(starts_with("cayley_")) %>% select(ends_with("left")), "Age + Sex + DX", data) #%>% mutate(Variable=paste0(Variable, "_left"))
	right_dx <- getLms(data %>% select(starts_with("cayley_")) %>% select(ends_with("right")), "Age + Sex + DX", data) #%>% mutate(Variable=paste0(Variable, "_right"))

	cayley_dx_quant <- rbind(cayley_dx_quant, cbind(Dataset=dataset, Label="Left Yeo", left_dx))
	cayley_dx_quant <- rbind(cayley_dx_quant, cbind(Dataset=dataset, Label="Right Yeo", right_dx))
}

combined_df_TOPSY_quant_yeo17 <- cbind(combined_df_TOPSY, cayley_ct_left_yeo17_quant %>% append_colnames(., "_yeo7_left"), cayley_ct_right_yeo17_quant %>% append_colnames(., "_yeo7_right"))

cayley_dx_quant_yeo17 <- data.frame()

for (dataset in unique(combined_df_TOPSY_quant_yeo17$Dataset)){
	
	data <- subset(combined_df_TOPSY_quant_yeo17, Dataset==dataset)

	cat("\n")
	cat(dataset)

	data %<>% filter(DX=="HC" | DX=="FEP" | DX=="Scz")
	data$DX <- replace(data$DX, data$DX=="FEP", "Scz")
	data$DX %<>% droplevels()

	if (dataset=="HCPEP") {
		data %<>% filter(DX_2!="Other")
	}

	left_dx <- getLms(data %>% select(starts_with("cayley_")) %>% select(ends_with("left")), "Age + Sex + DX", data) #%>% mutate(Variable=paste0(Variable, "_left"))
	right_dx <- getLms(data %>% select(starts_with("cayley_")) %>% select(ends_with("right")), "Age + Sex + DX", data) #%>% mutate(Variable=paste0(Variable, "_right"))

	cayley_dx_quant_yeo17 <- rbind(cayley_dx_quant_yeo17, cbind(Dataset=dataset, Label="Left Yeo", left_dx))
	cayley_dx_quant_yeo17 <- rbind(cayley_dx_quant_yeo17, cbind(Dataset=dataset, Label="Right Yeo", right_dx))
}

###Cross-dataset comparisons of effect sizes

reshape2::dcast(cayley_dx_quant %>% select(Dataset, Variable, t_DXScz),  Dataset ~ Variable)

test <- reshape2::dcast(cayley_dx_quant %>% select(Dataset, Variable, t_DXScz),  Variable ~ Dataset)
test_mat <- as.matrix(test %>% select(-Variable))
rownames(test_mat) <- test$Variable
pheatmap::pheatmap(test_mat)


rownames(test) <- test[,1]
test <- test[,-1]


test %<>% select(-Variable)
test %<>% filter(., !grepl("0", Variable))
summary(cor(test)[lower.tri(cor(test))])

pairwise_cors(test)


test <- reshape2::dcast(cayley_dx_quant_yeo17 %>% select(Dataset, Variable, t_DXScz),  Variable ~ Dataset)
test_mat <- as.matrix(test %>% select(-Variable))
rownames(test_mat) <- test$Variable
pheatmap::pheatmap(test_mat)


###Median rank correlation comparisons

ct_compare_correlations <- data.frame()
ct_compare_correlations_yeo17 <- data.frame()
for (dataset in unique(combined_df_TOPSY$Dataset)){
	
	df <- subset(combined_df_TOPSY, Dataset==dataset)
	df %<>% filter(DX=="FEP" | DX=="Scz")
	df$DX %<>% droplevels()
	
	#df %<>% filter(Age < 41)
	#data_split <- split(data, data$DX)
	if (dataset=="HCPEP") {
		df %<>% filter(DX_2!="Other")
	}
	
	left_ct_rank <- t(vertexTableRank(df$left_ct))
	right_ct_rank <- t(vertexTableRank(df$right_ct))

	left_ct_yeo <- getMedians(left_ct_rank, stat_maps_left$yeo7)
	right_ct_yeo <- getMedians(right_ct_rank, stat_maps_right$yeo7)

	left_ct_yeo17 <- getMedians(left_ct_rank, stat_maps_left$yeo17)
	right_ct_yeo17 <- getMedians(right_ct_rank, stat_maps_right$yeo17)

	#left_ct_vonEconomo <- getMedians(left_ct_rank, stat_maps_left$vonEconomo)
	#right_ct_vonEconomo <- getMedians(right_ct_rank, stat_maps_right$vonEconomo)

	colnames(left_ct_yeo) <-  c("Unmatched_left", "Visual_left", "Somatomotor_left", "DorsalAtt_left", "VentralAtt_left", "Limbic_yeo_left", "Frontoparietal_left", "Default_left")
	colnames(right_ct_yeo) <-  c("Unmatched_right", "Visual_right", "Somatomotor_right", "DorsalAtt_right", "VentralAtt_right", "Limbic_yeo_right", "Frontoparietal_right", "Default_right")

	left_ct <- vertexTable(df$left_ct)
	right_ct <- vertexTable(df$right_ct)

	left_ct_yeo_raw <- getMedians(left_ct, stat_maps_left$yeo7)
	right_ct_yeo_raw <- getMedians(right_ct, stat_maps_right$yeo7)

	left_ct_yeo17_raw <- getMedians(left_ct, stat_maps_left$yeo17)
	right_ct_yeo17_raw <- getMedians(right_ct, stat_maps_right$yeo17)


	#colnames(left_ct_vonEconomo) <-  c("Unmatched_left", "Motor_left", "Association-FP_left", "Association-FT_left", "Sensory-2_left", "Sensory-1_left", "Limbic_left", "Insular_left")
	#colnames(right_ct_vonEconomo) <-  c("Unmatched_right", "Motor_right", "Association-FP_right", "Association-FT_right", "Sensory-2_right", "Sensory-1_right", "Limbic_right", "Insular_right")

	if (dataset %in% c("BrainGluSchi", "COBRE", "HCPEP", "TOPSY")) {
		scores <- df %>% select(starts_with(c("PANSS_")))
	} else {
		scores <- df %>% select(ends_with(c("_global")))
	}
	
	
	for(i in 1:ncol(scores)){
		left_yeo <- cbind(Score=names(scores)[i], dataset, Label="Left Yeo", getCors(scores[[i]], left_ct_yeo), getCors(scores[[i]], left_ct_yeo_raw) %>% select(-Variable) %>% append_colnames(., "_raw"))
		right_yeo <- cbind(Score=names(scores)[i], dataset, Label="Right Yeo", getCors(scores[[i]], right_ct_yeo), getCors(scores[[i]], right_ct_yeo_raw) %>% select(-Variable) %>% append_colnames(., "_raw"))

		#left_von <- cbind(Score=names(scores)[i], dataset, Label="Left von Economo", getCors(scores[[i]], left_ct_vonEconomo))
		#right_von <- cbind(Score=names(scores)[i], dataset, Label="Right von Economo", getCors(scores[[i]], right_ct_vonEconomo))

		ct_compare_correlations <- rbind(ct_compare_correlations, rbind(left_yeo, right_yeo))
	}

	for(i in 1:ncol(scores)){
		left_yeo <- cbind(Score=names(scores)[i], dataset, Label="Left Yeo17", getCors(scores[[i]], left_ct_yeo17), getCors(scores[[i]], left_ct_yeo17_raw) %>% select(-Variable) %>% append_colnames(., "_raw"))
		right_yeo <- cbind(Score=names(scores)[i], dataset, Label="Right Yeo17", getCors(scores[[i]], right_ct_yeo17), getCors(scores[[i]], right_ct_yeo17_raw) %>% select(-Variable) %>% append_colnames(., "_raw"))

		#left_von <- cbind(Score=names(scores)[i], dataset, Label="Left von Economo", getCors(scores[[i]], left_ct_vonEconomo))
		#right_von <- cbind(Score=names(scores)[i], dataset, Label="Right von Economo", getCors(scores[[i]], right_ct_vonEconomo))

		ct_compare_correlations_yeo17 <- rbind(ct_compare_correlations_yeo17, rbind(left_yeo, right_yeo))
	}
}

ct_compare_correlations %<>% mutate_at(vars(starts_with(c("Correlation", "pvalue"))), as.numeric)
ct_compare_correlations_yeo17 %<>% mutate_at(vars(starts_with(c("Correlation", "pvalue"))), as.numeric)

test <- reshape2::dcast(ct_compare_correlations %>% filter(Score=="PANSS_pos" | Score=="SAPS_global") %>% select(dataset, Variable, Correlation),  dataset ~ Variable)
test_mat <- test %>% select(-dataset) %>% as.matrix()
rownames(test_mat) <- test$dataset
pheatmap::pheatmap(test_mat)

test <- reshape2::dcast(ct_compare_correlations %>% filter(Score=="PANSS_neg" | Score=="SANS_global") %>% select(dataset, Variable, Correlation),  dataset ~ Variable)
test_mat <- test %>% select(-dataset) %>% as.matrix()
rownames(test_mat) <- test$dataset
pheatmap::pheatmap(test_mat)

test <- ct_compare_correlations %>% select(Variable, pvalue, Score, pvalue_raw) %>% melt(., id.vars=c("Score", "Variable")
test %<>% mutate(Variable=gsub("_left|_right", "", Variable))%>% filter(Variable!="Unmatched")  %>% filter(!grepl("PANSS_total|PANSS_g", Score))

ggplot(test, aes(x=Variable, y=-log10(value), color=Score)) + geom_point(position=position_dodge(0.5), size=2) + facet_wrap(~variable, ncol=2)

ggplot(test, aes(x=Variable, y=-log10(value), color=variable)) + geom_point(position=position_dodge(0.5), size=2) + facet_wrap(~variable, ncol=2)


###HCP-specific testing (DX_2) comparisons

####################################
###Comparison: raw cortical thickness & network medians

left_ct <- vertexTable(hcp$left_ct)
right_ct <- vertexTable(hcp$right_ct)

left_ct_yeo <- getMedians(left_ct, stat_maps_left$yeo7)
right_ct_yeo <- getMedians(right_ct, stat_maps_right$yeo7)

hcp_medians_raw <- cbind(hcp, left_ct_yeo, right_ct_yeo)
hcp_medians_raw_FEP <- subset(hcp_medians_raw, DX=="FEP")

getCors(hcp_medians_raw_FEP$totalP, hcp_medians_raw_FEP[(ncol(hcp_medians_raw_FEP)-15):ncol(hcp_medians_raw_FEP)])
getCors(hcp_medians_raw_FEP$totalN, hcp_medians_raw_FEP[(ncol(hcp_medians_raw_FEP)-15):ncol(hcp_medians_raw_FEP)])
getCors(hcp_medians_raw_FEP$WksTo50P, hcp_medians_raw_FEP[(ncol(hcp_medians_raw_FEP)-15):ncol(hcp_medians_raw_FEP)])
getCors(hcp_medians_raw_FEP$WksToCGI2, hcp_medians_raw_FEP[(ncol(hcp_medians_raw_FEP)-15):ncol(hcp_medians_raw_FEP)])

###Plotting for comparison
hcp_medians_FEP$Method <- "CT-ranked"
hcp_medians_raw_FEP$Method <- "CT"

hcp_medians_combined <- rbind(hcp_medians_FEP, hcp_medians_raw_FEP)


#####

test <-ct_compare_correlations %>% filter(Score %in% c("SAPS_global", "SANS_global", "PANSS_pos", "PANSS_neg")) %>% filter(., grepl("Yeo", Label))
test <-ct_compare_correlations %>% filter(Score %in% c("SAPS_global", "PANSS_pos")) %>% filter(., grepl("Yeo", Label))

test %<>% mutate_at(vars("Correlation", "pvalue"), as.numeric)

###Plotting: dataset (organized by median age), network-symptom correlations
test <-ct_compare_correlations %>% filter(Score %in% c("SAPS_global", "PANSS_pos")) %>% filter(., grepl("Yeo", Label))

test %<>% select(dataset, Correlation, Variable) %>% pivot_wider(names_from = Variable, values_from=Correlation)

test2 <- test %>% select(-dataset) %>% mutate_all(as.numeric) %>% as.matrix()
rownames(test2) <- test$dataset


test <-ct_compare_correlations %>% filter(Score %in% c("SANS_global", "PANSS_neg")) %>% filter(., grepl("Yeo", Label))

test %<>% select(dataset, Correlation, Variable) %>% pivot_wider(names_from = Variable, values_from=Correlation)

test2 <- test %>% select(-dataset) %>% mutate_all(as.numeric) %>% as.matrix()
rownames(test2) <- test$dataset


#####################################################


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
