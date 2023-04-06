library(readr)
library(magrittr)
library(dplyr)
library(ggplot2)
library(RMINC)
library(matrixStats)
library(Rankcluster)

############
###Functions

parseLm <- function (model) 
{
    modelcoef <- summary(model)[[4]]
    rsquared <- summary(model)[[9]]
    fstatistic <- as.numeric(summary(model)[[10]][1])
    DF<- model$df.residual
    rows <- nrow(modelcoef)
    cols <- ncol(modelcoef)
    
    betacolnames= paste("beta_", rownames(modelcoef)[2:rows], sep="")
    secolnames= paste("se_", rownames(modelcoef)[2:rows], sep="")
    tcolnames= paste("t_", rownames(modelcoef)[2:rows], sep="")
    pcolnames= paste("p_", rownames(modelcoef)[2:rows], sep="")
    
    output=cbind(fstatistic, DF, rsquared)
    
    a<-c(1,2,3,4)
    for (i in a){
        for (j in 2:rows) {
            output= cbind(output, modelcoef[j,i])
        }
    }
    
    colnames(output) <- c("F_statistic","DF","R_squared", betacolnames, secolnames, tcolnames, pcolnames)
    return(as.data.frame(output))
}

extract_column <- function (filename, column = 1) #RMINC internal function, called via RMINC:::extract_column
{
    ext = tools::file_ext(filename)
    if (ext %in% c("gz", "xz", "bz2", "GZ", "XZ", "BZ2")) {
        filename_ = tools::file_path_sans_ext(filename)
        ext = tools::file_ext(filename_)
    }
    if (ext %in% c("csv", "CSV")) {
        data_ = readr::read_csv(filename, col_names = T, col_types = readr::cols(.default = readr::col_double()))
    }
    else {
        data_ = readr::read_table(filename, col_names = F, col_types = readr::cols(.default = readr::col_double()))
    }
    return(as.matrix(data_[, column]))
}

vertexTableRank<- function (filenames, column = 1) 
{
	#vertexTable function, added ranking to get ranks instead of values
	#Rows as subjects, columns are vertices--transposed from the usual vertexTable output
	#rowRanks from MatrixStats package
	if (is.factor(filenames)) 
		filenames <- as.character(filenames)
	if (!is.character(filenames)) 
		stop("filenames must be either a character or factor vector")
	nSubjects <- length(filenames)
	nvertices <- nrow(extract_column(filenames[1], column = column))
	output <- matrix(nrow = nvertices, ncol = nSubjects)
	for (i in 1:nSubjects) {
		output[, i] <- as.matrix(extract_column(filenames[i], column = column))
	}
	return(output %<>% t() %>% rowRanks(ties.method="average")) ###TIES ARE AVERAGED
}

vertexTableRank<- function (filenames, column = 1) 
{
	#vertexTable function, added ranking to get ranks instead of values
	#Rows as subjects, columns are vertices--transposed from the usual vertexTable output
	#rowRanks from MatrixStats package

	###NO TIES in this.

	if (is.factor(filenames)) 
		filenames <- as.character(filenames)
	if (!is.character(filenames)) 
		stop("filenames must be either a character or factor vector")
	nSubjects <- length(filenames)
	nvertices <- nrow(extract_column(filenames[1], column = column))
	output <- matrix(nrow = nvertices, ncol = nSubjects)
	for (i in 1:nSubjects) {
		output[, i] <- as.matrix(extract_column(filenames[i], column = column))
	}
	return(output %<>% t() %>% rowRanks(ties.method="first")) ###NO TIES, first element gets rank.
}


compare_models <- function(dataset, column, model) {
	#Function to compare raw vs. ranked values in predicting
	raw <- vertexTable(dataset[[column]])
	ranked <- t(vertexTableRank(dataset[[column]]))

	df <- data.frame(NULL)

	for (i in 1:nrow(raw)) {
		model_raw <- glm(as.formula(paste(model,"raw[i,]", sep="~")), data=dataset, family="binomial")
		model_ranked <- glm(as.formula(paste(model,"ranked[i,]", sep="~")), data=dataset, family="binomial")
		diff <- abs(model_raw$aic) - abs(model_ranked$aic) #If positive, AIC is greater in raw

		#pval_raw <- coef(summary(model_raw))[,'Pr(>|z|)'][-1]
		#pval_ranked <- coef(summary(model_raw))[,'Pr(>|z|)'][-1]


		df <- rbind(df,
		cbind(append_colnames(parseGlm(model_raw), "_raw"), append_colnames(parseGlm(model_ranked), "_ranked"), diff) 
		)
	}
	#colnames(df) <- c("aic_raw", "aic_ranked", "diff")
	colnames(df) %<>% gsub("_raw\\[i, \\]", "", .) %>% gsub("_ranked\\[i, \\]", "", .)
	
	df %<>% mutate(q_raw = p.adjust(p_raw, method="fdr"), q_ranked= p.adjust(p_ranked, method="fdr"))
	return(df)
}

parse_ttest <- function(ttest) {
  t_stat <- ttest$statistic
  p_value <- ttest$p.value
  results_df <- data.frame(t_stat, p_value)
  return(results_df)
}

cor_columns <- function(df1, df2) {
  # Ensure both data frames have same dimensions
  if (!all(dim(df1) == dim(df2))) {
    stop("Input data frames must have the same dimensions.")
  }
  
  # Initialize matrix to store correlations
  n <- ncol(df1)
  cor_matrix <- matrix(0, nrow = 1, ncol = n)
  
  # Compute correlations between same columns in both data frames
  for (i in 1:n) {
    cor_matrix[1, i] <- cor(df1[, i], df2[, i])
  }
  
  # Return correlation matrix
  return(t(data.frame(cor_matrix)))
}

threshold_column <- function(column, percentile) {
  threshold_value <- quantile(column, probs = percentile)
  return(ifelse(column <= threshold_value, 1, 0))
}


pairwise_similarity <- function(df) {
	
	result <- data.frame()

	for (i in 1:(ncol(df)-1)) {
		colA <- colnames(df)[i] #Main column to be comparing against
		
		for (j in (i+1):ncol(df)){
			colB <- colnames(df)[j] #Second column being compared
			
			table <- table(df[,i], df[,j])
			chitest <- chisq.test(table)

			result <- rbind(result, cbind(colA, colB, table[4], chitest$statistic))
		}
	}
	
	colnames(result) <- c("ColA", "ColB", "Overlap", "Chisq")
	result %<>% mutate_at(vars("Overlap", "Chisq"), as.numeric)
	return(result)
}


indiv_rankdiff <- function(reference, sample) {
	#Individual rank differences based on a reference dataset
	#indiv_rankdiff(gf_HC$left_CT, gf$left_CT) for example.
	
	refset<- vertexTableRank(reference)
	testset<- vertexTableRank(sample)

	rankdiff <- matrix(nrow=nrow(testset), ncol=ncol(testset))
	
	#For each row in testset (each subject), subtract it from every row in refset (i.e. healthy population)
	#Then get the column means, and store in new matrix rankdiff
	for (i in 1:nrow(testset)) {
		newrow<- colMeans(sweep(refset, 2, testset[i,], FUN="-")) ###Mean of rank differences
		rankdiff[i,] <-  newrow
	}
	return(rankdiff)
}

indiv_rankdiff <- function(reference, sample) {
	#Individual rank differences based on a reference dataset
	#indiv_rankdiff(gf_HC$left_CT, gf$left_CT) for example.
	
	refset<- vertexTableRank(reference)
	testset<- vertexTableRank(sample)

	rankdiff <- matrix(nrow=nrow(testset), ncol=ncol(testset))
	
	#For each row in testset (each subject), subtract it from every row in refset (i.e. healthy population)
	#Then get the column means, and store in new matrix rankdiff
	for (i in 1:nrow(testset)) {
		newrow<- colMedians(sweep(refset, 2, testset[i,], FUN="-")) ###Median of rank differences
		rankdiff[i,] <-  newrow
	}
	return(rankdiff)
}


cor_vecmat <- function(y, x, type ="pearson") { #y is vector, x is matrix
  res <- matrix(nrow=nrow(x))

  for (i in 1:nrow(x)) res[i] <- cor(y, x[i,], method = type)  
  res
}

indiv_rankdiff_kendall <- function(reference, sample) {
	#Individual rank differences based on a reference dataset
	#indiv_rankdiff(gf_HC$left_CT, gf$left_CT) for example.
	
	refset<- vertexTableRank(reference)
	testset<- vertexTableRank(sample)

	rankdiff <- matrix(nrow=nrow(testset))
	
	#For each row in testset (each subject), subtract it from every row in refset (i.e. healthy population)
	#Then get the column means, and store in new matrix rankdiff
	for (i in 1:nrow(testset)) {
		cors <- cor_vecmat(testset[i,], refset, type="kendall") 
		rankdiff[i] <- median(cors[cors < 1,]) ###Exclude perfect cors (self-correlations)
	}
	return(rankdiff)
}

indiv_distance_kendall_sub <- function(reference, sample, map) {
	#Individual rank differences based on a reference dataset, PER DEFINED SUBREGION
	#indiv_distance_kendall_sub(gf_HC$left_CT, gf$left_CT, stat_maps_left$yeo7) for example.

	refset<- vertexTableRank(reference)
	testset<- vertexTableRank(sample)

	rankdiff <- matrix(nrow=nrow(testset), ncol=length(unique(map)))

	total=length(rankdiff)
	modulo=10
	count=0
	
	#For each row in testset (each subject), subtract it from every row in refset (i.e. healthy population)
	#Then get the column means, and store in new matrix rankdiff
	
	for (i in 1: length(unique(map))) {
		
		label <- sort(unique(map))[i]
		
		#Subset total vertices into label-specific vertices, then re-rank them.
		refset_label <- refset[,map==label] %>% rowRanks()
		testset_label <- testset[,map==label] %>% rowRanks()
		
		for (j in 1:nrow(testset_label)) {
			cors <- cor_vecmat(testset_label[j,], refset_label, type="kendall") 
			rankdiff[j,i] <- median(cors[cors < 0.999,]) ###Exclude perfect cors (self-correlations)

			count=count+1
		
			if (count%%modulo == 0) {
				cat(format((count/total) * 100, digits = 3))
				cat("%  ")
			}
		}
	}

	colnames(rankdiff) <- paste("kendall_",sort(unique(map)), sep="")
	return(data.frame(rankdiff))
}


cayley_vecmat <- function(y, x) { #y is vector, x is matrix
  res <- matrix(nrow=nrow(x))

  for (i in 1:nrow(x)) res[i] <- distCayley(y, x[i,])  
  res
}

indiv_distance_cayley_sub <- function(reference, sample, map) {
	#Individual rank differences based on a reference dataset, PER DEFINED SUBREGION
	#indiv_distance_kendall_sub(gf_HC$left_CT, gf$left_CT, stat_maps_left$yeo7) for example.

	#Checks reference input type
	#If character (df$left_ct), then read in as ranks.
	#Rows are subjects, and columns are vertices
	if (typeof(reference)=="character") {
		refset<- vertexTableRank(reference)
	}
	
	#If integer (already ranked: output of row/colRanks functions), then assign to refset.
	if (typeof(reference)=="integer") {
		refset<- reference
	}
	
	testset<- vertexTableRank(sample)

	rankdiff <- matrix(nrow=nrow(testset), ncol=length(unique(map)))

	total=length(rankdiff)
	count=0
		
	for (i in 1: length(unique(map))) {
		
		label <- sort(unique(map))[i]
		
		#Subset total vertices into label-specific vertices, then re-rank them.
		refset_label <- refset[,map==label] %>% rowRanks()
		testset_label <- testset[,map==label] %>% rowRanks()
		
		for (j in 1:nrow(testset_label)) {
			cors <- cayley_vecmat(testset_label[j,], refset_label) 
			rankdiff[j,i] <- median(cors[cors > 0,]) ###Exclude 0 distance (self-distance)
		
			count=count+1
		
			if (count > 0) {
				cat("\r")
				cat("Progress: ", round(count / total * 100), "%")
			}
		
		}

	}

	colnames(rankdiff) <- paste("cayley_",sort(unique(map)), sep="")
	return(data.frame(rankdiff))
}




cayley_dist_ref <- function(reference, sample, map) {
	#Individual rank differences based on a reference, PER DEFINED SUBREGION
	#cayley_dist_ref(gf_HC$left_CT, gf$left_CT, stat_maps_left$yeo7) for example.
	#The median ranks of the reference set is computed & ranked again to preserve rankings

	refset<- vertexTableRank(reference) %>% colMedians() %>% rank(ties.method="first") %>% t()
	testset<- vertexTableRank(sample)

	rankdiff <- matrix(nrow=nrow(testset), ncol=length(unique(map)))

	for (i in 1: length(unique(map))) {
		
		label <- sort(unique(map))[i]
		
		#Subset total vertices into label-specific vertices, then re-rank them.
		refset_label <- refset[,map==label] %>% rank(ties.method="first")
		testset_label <- testset[,map==label] %>% rowRanks()
		
		rankdiff[,i] <- cayley_vecmat(t(refset_label), testset_label)
	}

	colnames(rankdiff) <- paste("cayley_",sort(unique(map)), sep="")
	return(data.frame(rankdiff))
}

hamming_vecmat <- function(y, x) { #y is vector, x is matrix
  res <- matrix(nrow=nrow(x))

  for (i in 1:nrow(x)) res[i] <- distHamming(y, x[i,])  
  res
}

hamming_dist_ref <- function(reference, sample, map) {
	#Individual rank differences based on a reference, PER DEFINED SUBREGION
	#hamming_dist_ref(gf_HC$left_CT, gf$left_CT, stat_maps_left$yeo7) for example.
	#The median ranks of the reference set is computed & ranked again to preserve rankings

	refset<- vertexTableRank(reference) %>% colMedians() %>% rank(ties.method="first") %>% t()
	testset<- vertexTableRank(sample)

	rankdiff <- matrix(nrow=nrow(testset), ncol=length(unique(map)))

	for (i in 1: length(unique(map))) {
		
		label <- sort(unique(map))[i]
		
		#Subset total vertices into label-specific vertices, then re-rank them.
		refset_label <- refset[,map==label] %>% rank(ties.method="first")
		testset_label <- testset[,map==label] %>% rowRanks()
		
		rankdiff[,i] <- hamming_vecmat(t(refset_label), testset_label)
	}

	colnames(rankdiff) <- paste("hamming_",sort(unique(map)), sep="")
	return(data.frame(rankdiff))
}

getMedians <- function(vtable, map) 
{
	datalist = list()
	for (i in 1: length(unique(map))) {
		label <- sort(unique(map))[i]
		means <- colMedians(vtable[map==label,])
		datalist[[i]] <- means
	}
	output <- datalist %<>% as.data.frame()
	colnames(output) <-  sort(unique(map))
	
	return(output)
}

append_colnames <- function(dataframe, add) {
	colnames(dataframe) <- paste(colnames(dataframe), add, sep="")
	return(dataframe)
}

threshold_groupdiff <- function(dataset, rankdiff, group, model, min, max, steps) {
	#Thresholding testing function with specified model
	#Testing 3 different outcome measures--total number of thresholded vertices, and neg/pos vertices separately
	df <- data.frame(NULL)

	for (i in seq(min, max, by=steps)) { #Range chosen based on min/max of left ct rank differences
		threshold=i
		dataset$total <- rowSums(abs(rankdiff) > threshold)
		dataset$positive <- rowSums(rankdiff > threshold)
		dataset$negative <- rowSums(rankdiff < -threshold)
		
		model1 <- lm(as.formula(paste("total", model, sep="~")), dataset) %>% parseLm %>% append_colnames("_total")
		model2 <- lm(as.formula(paste("positive", model, sep="~")), dataset) %>% parseLm %>% append_colnames("_pos")
		model3 <- lm(as.formula(paste("negative", model, sep="~")), dataset) %>% parseLm %>% append_colnames("_neg")

		#Get average number per group
		mean_total <- tapply(dataset$total, group, mean) %>% t() %>% append_colnames("_total")
		mean_pos <- tapply(dataset$positive, group, mean) %>% t() %>% append_colnames("_pos")
		mean_neg <- tapply(dataset$negative, group, mean) %>% t() %>% append_colnames("_neg")

		df <- rbind(df, cbind(threshold, mean_total, mean_pos, mean_neg, model1, model2, model3))
	}
	return(df)
}

threshold_groupdiff_window <- function(dataset, rankdiff, group, model) {
	#Thresholding testing function with specified model

	#Sliding window approach--min is minimum window size, max is maximum window size
	#Minimum is 250, maximum 5000, scaled by 250?
	#Sliding window by 100. 0 to 100, 100 to 200.

	windowEnd=35000
	windowSize=250

	df <- data.frame(NULL)

	while ( windowSize < 5001) {
		windowStart=0

		while (windowStart + windowSize < windowEnd) {
			
			lower= windowStart
			upper= windowStart + windowSize

			dataset$total <- rowSums(lower < abs(rankdiff) & abs(rankdiff) < upper)
			dataset$positive <- rowSums(lower < rankdiff & rankdiff < upper)
			dataset$negative <- rowSums(-upper < rankdiff & rankdiff < -lower)
		
			model1 <- lm(as.formula(paste("total", model, sep="~")), dataset) %>% parseLm %>% append_colnames("_total")
			model2 <- lm(as.formula(paste("positive", model, sep="~")), dataset) %>% parseLm %>% append_colnames("_pos")
			model3 <- lm(as.formula(paste("negative", model, sep="~")), dataset) %>% parseLm %>% append_colnames("_neg")

			#Get average number per group
			mean_total <- tapply(dataset$total, group, mean) %>% t() %>% append_colnames("_total")
			mean_pos <- tapply(dataset$positive, group, mean) %>% t() %>% append_colnames("_pos")
			mean_neg <- tapply(dataset$negative, group, mean) %>% t() %>% append_colnames("_neg")

			df <- rbind(df, cbind(windowStart, windowStart+windowSize, windowSize, mean_total, mean_pos, mean_neg, model1, model2, model3))

			windowStart <- windowStart + 100

		}
		
		windowSize <- windowSize + 250
	}	
		
	return(df)
}

threshold_reg <- function(dataset, rankdiff, model, min, max, steps) {
	#Thresholding testing function with specified model
	#Testing 3 different outcome measures--total number of thresholded vertices, and neg/pos vertices separately
	#For regressions, not group differences
	df <- data.frame(NULL)

	for (i in seq(min, max, by=steps)) { #Range chosen based on min/max of left ct rank differences
		threshold=i
		dataset$total <- rowSums(abs(rankdiff) > threshold)
		dataset$positive <- rowSums(rankdiff > threshold)
		dataset$negative <- rowSums(rankdiff < -threshold)
		
		model1 <- lm(as.formula(paste("total", model, sep="~")), dataset) %>% parseLm %>% append_colnames("_total")
		model2 <- lm(as.formula(paste("positive", model, sep="~")), dataset) %>% parseLm %>% append_colnames("_pos")
		model3 <- lm(as.formula(paste("negative", model, sep="~")), dataset) %>% parseLm %>% append_colnames("_neg")

		df <- rbind(df, cbind(threshold, model1, model2, model3))
	}
	return(df)
}

average_rank_diff <- function(reference, sample) {
	#Average rank differences based on a reference dataset
	#average_rank_diff(gf_HC$left_CT, gf_FEP$left_CT) for example.
	
	refset<- vertexTableRank(reference) %>% colMeans()
	testset<- vertexTableRank(sample) %>% colMeans()
	output <- as.data.frame(cbind(refset, testset, testset-refset))
	colnames(output) <- c("reference", "test", "diff")
	return(output)
}

average_rank_diff_permute <- function(reference, sample, n_permute) {
	#Average rank differences based on a reference dataset
	#average_rank_diff(gf_HC$left_CT, gf_FEP$left_CT) for example.
	
	refset<- vertexTableRank(reference) 
	testset<- vertexTableRank(sample)

	#Combined dataset for permutation testing
	combined <- rbind(refset, testset)

	n_reference<-nrow(refset)
	n_sample <- nrow(sample)

	perm_matrix <- matrix(data=NA, nrow= ncol(refset), ncol=n_permute)
	
	for (i in 1:n_permute) {
		#Sample
		perm_ref <- sample(1:nrow(combined), size=n_reference, replace=FALSE)

		perm_ref <- combined[perm_ref,] %>% colMeans()
		perm_sample <- combined[-perm_ref,] %>% colMeans()
		perm_matrix[,i] <- perm_ref - perm_sample

		# Print percentage of permutations completed
		if (i %% (n_permute / 20) == 0) {
			message(paste0(round(i / n_permute * 100), "% of permutations completed."))
		}
	}

	refset %<>% colMeans()
	testset %<>% colMeans()
	diff <- data.frame(testset-refset)

	pvals <- matrix(data=NA, nrow=nrow(diff), ncol=1)
	
	for (j in 1:nrow(pvals)) {
		
		if(diff[j,] > 0) {
			pval <- (sum(perm_matrix[j,] > diff[j,]) + 1) / (ncol(perm_matrix) + 1)
		} else {
			pval <- (sum(perm_matrix[j,] < diff[j,]) + 1) / (ncol(perm_matrix) + 1)
		}

		pvals[j,] <- pval
	}

	output <- data.frame(cbind(refset, testset, diff, pvals))
	
	colnames(output) <- c("reference", "test", "diff", "pvals")
	return(output)
}


average_vertex_diff <- function(reference, sample) {
    #Differences in average values based on a reference dataset
    #average_vertex_diff(gf_HC$left_ct, gf_FEP$left_ct) for example.
    
    refset<- t(vertexTable(reference)) %>% colMeans()
    testset<- t(vertexTable(sample)) %>% colMeans()
    output <- as.data.frame(cbind(refset, testset, testset-refset))
    colnames(output) <- c("reference", "test", "diff")
    return(output)
}



getSums <- function(vtable, map) 
{
	datalist = list()
	for (i in 1: length(unique(map))) {
		label <- sort(unique(map))[i]
		means <- colSums(vtable[map==label,])
		datalist[[i]] <- means
	}
	output <- datalist %<>% as.data.frame()
	colnames(output) <-  sort(unique(map))
	
	return(output)
}

getMedians <- function(vtable, map) 
{
	datalist = list()
	for (i in 1: length(unique(map))) {
		label <- sort(unique(map))[i]
		means <- colMedians(vtable[map==label,])
		datalist[[i]] <- means
	}
	output <- datalist %<>% as.data.frame()
	colnames(output) <-  sort(unique(map))
	
	return(output)
}


getCors <- function(var, columns) {
	datalist = data.frame()

	for (i in 1:ncol(columns)) {
		test <- cor.test(var, columns[,i])

		datalist <- rbind(datalist, cbind(colnames(columns)[i], test$estimate, test$p.value))
	}
	colnames(datalist) <- c("Variable", "Correlation", "pvalue")
	return(datalist)
}

getLms <- function(columns, model, df) {
	datalist = data.frame()

	for (i in 1:ncol(columns)) {
		result <- lm(as.formula(paste(columns[i], model, sep="~")), data=df) %>% parseLm()
		datalist <- rbind(datalist, cbind(names(columns)[i], result))
	}
	colnames(datalist)[1] <- "Variable"
	return(datalist)
}

parseGlm <- function(model) {
    modelcoef <- coef(summary(model))
    #rsquared <- summary(model)[[9]]
    #fstatistic <- as.numeric(summary(model)[[10]][1])
    #DF<- model$df.residual
    
    rows <- nrow(modelcoef)
    cols <- ncol(modelcoef)
    
    betacolnames= paste("odds_", rownames(modelcoef)[2:rows], sep="")
    secolnames= paste("se_", rownames(modelcoef)[2:rows], sep="")
    tcolnames= paste("z_", rownames(modelcoef)[2:rows], sep="")
    pcolnames= paste("p_", rownames(modelcoef)[2:rows], sep="")
    
    AIC<- model$aic
    dev<- model$deviance
    df<- model$df.residual
    output=cbind(AIC, dev, df)
    
    a<-c(1,2,3,4)
    for (i in a){
        for (j in 2:rows) {
            output= cbind(output, modelcoef[j,i])
        }
    }
    
    colnames(output) <- c("AIC","Dev","DF", betacolnames, secolnames, tcolnames, pcolnames)
    return(as.data.frame(output))
}

#####AHBA analysis

ahba_CIVET <- function(data) {
	#Data is left cortex data in MNI space (CIVET)--for example, results$t_stat as input.
	data %<>% as.data.frame %>% set_colnames("data")
	data$Vertex <- 1:nrow(data)

	ahba <- read.csv('/projects/Allen_to_CIVET/Analysis_mapping/Donors_cortex_20200823.csv', sep="\t")

	merged <- merge(data, ahba, select="Vertex", all.y=TRUE)

	nrow(merged) #1236 vertices to samples

	library(lmerTest)

	ahba_mat <- matrix(data=NA, nrow=15126, ncol=3)

	for (i in 25:ncol(merged)){
		model <- summary(lmer(data ~ merged[,i] + (1|Donor), data=merged, verbose=0))
		j=i-24
		ahba_mat[j,1] <- colnames(merged)[i]
		ahba_mat[j,2:3] <- model$coefficients[2,4:5] ###t-statistic and p-value
	}

	ahba_mat <- as.data.frame(ahba_mat)
	colnames(ahba_mat) <- c("Gene", "tstat", "pval")
	ahba_mat %<>% mutate_at(c('tstat', 'pval'), as.numeric)
	ahba_mat$qval <- p.adjust(ahba_mat$pval, method="fdr")

	ahba_mat <- ahba_mat[order(ahba_mat$tstat, decreasing=TRUE),]

	ahba_mat %<>% mutate(decile_tstat=ntile(tstat, 10)) %>% mutate(decile_pval=ntile(-log10(pval), 10))

	return(ahba_mat)
}

ahba_CIVET_write <- function(name, mat, output) {

	if (!dir.exists(output)){
		dir.create(output)
	} else {
		print("Dir already exists!")
	}
	write.table(mat, file=paste(output, paste(name, Sys.Date(), "results.csv", sep="_"), sep=""), sep="\t", row.names=F, col.names=T, quote=FALSE)
	write.table(mat %>% select(Gene) %>% head(n=nrow(mat)*0.10), file=paste(output, paste(name, Sys.Date(), "results_top10genes.csv", sep="_"), sep=""), sep="\t", row.names=F, col.names=F, quote=FALSE)
	write.table(mat %>% select(Gene) %>% tail(n=nrow(mat)*0.10), file=paste(output, paste(name, Sys.Date(), "results_bottom10genes.csv", sep="_"), sep=""), sep="\t", row.names=F, col.names=F, quote=FALSE)
	write.table(mat %>% select(Gene, tstat), file=paste(output, paste(name, Sys.Date(), "results_tstats.csv", sep="_"), sep=""), sep="\t", row.names=F, col.names=F, quote=FALSE)
	write.table(mat %>% select(Gene),file=paste(output, paste(name, Sys.Date(), "results_rankedgenes.csv", sep="_"), sep=""), sep="\t", row.names=F, col.names=F, quote=FALSE)

	for (i in 1:length(unique(mat$decile_tstat))) {
		decile <- unique(mat$decile_tstat)[i]
		sub <- subset(mat, decile_tstat==decile)
		write.table(sub$Gene, file=paste(output, paste(name, "decile_tstat", decile, "genes.txt", sep="_"), paste=""), sep="\t", row.names=F, col.names=F, quote=FALSE)
	}

	for (i in 1:length(unique(mat$decile_pval))) {
		decile <- unique(mat$decile_pval)[i]
		sub <- subset(mat, decile_pval==decile)
		write.table(sub$Gene, file=paste(output, paste(name, "decile_pval", decile, "genes.txt", sep="_"), paste=""), sep="\t", row.names=F, col.names=F, quote=FALSE)
	}
}
