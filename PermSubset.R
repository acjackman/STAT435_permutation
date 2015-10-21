library(combinat)
library(plyr)
library(doMC)
doMC::registerDoMC(cores=8)
grads <- read.csv("GradCleanup.csv")

# Subset data to look at only three years (5!)^3 permutations
grads <- subset(grads, Cohort > 2003)

# Block by Cohort, group by Borough
real_test <- friedman.test(GradRate ~ Borough | Cohort, data=grads)

# Create permutations across groups within blocks
gen_perms <- function(dataset){
	permn(dataset$GradRate)
}
block_perms <- dlply(grads, c("Cohort"),
	gen_perms
)

# Compute all possible combinations of different block permutations
num_perms <- factorial(length(unique(grads$Borough)))
combinations <- as.matrix(expand.grid(1:120,1:120,1:120))
str(combinations)
# Run the Permutation Test
get_perm <- function(block, iteration){
	block_perms[[ block ]][[iteration]]
}
block_names <- as.character(attr(block_perms, "split_labels")[,1])
borough_names <- unique(as.character(grads$Borough))
run_test <- function(block_perm_select){
	blk_sel <- as.numeric(block_perm_select)
	perm_grad_rate <- c(get_perm(block_names[1], blk_sel[1]),
		                get_perm(block_names[2], blk_sel[2]),
		                get_perm(block_names[3], blk_sel[3])
		              )

	s_perm <- data.frame(Borough = rep(borough_names,3),
		                 Cohort = rep(block_names, each=length(borough_names)),
					     GradRate = perm_grad_rate,
					     stringsAsFactors=TRUE)

	test <- friedman.test(GradRate ~ Borough | Cohort, data=s_perm)
	pval <- test$p.value
	names(pval) <- "p-value"
	c(test$statistic, pval)
}
results <- adply(combinations, 1, run_test, .parallel = TRUE)

keep_object_list <- c("block_perms","combinations", "get_perm", "results")
rm(list=ls()[!ls() %in% keep_object_list])
save.image("PermSubset.RData")
