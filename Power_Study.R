library(combinat)
library(plyr)
library(reshape2)
library(ggplot2)
library(doMC)
doMC::registerDoMC(cores=8)

run_perm_test <- function(c_dist, t_dist){
	### Take generated data
	control <- eval(c_dist)
	treatment <- eval(t_dist)

	group_labels <- c(rep("control", length(control)),
	                  rep("treatment", length(treatment)))
	fake_data <- c(control,treatment)

	fake_dataset <- data.frame(group = group_labels, value = fake_data)

	### Calculate all different possible combinations of data,
	# taken the number of treatments at a time
	# Using indices just in case of ties and how the other
	# observations are selected
	i_data <- seq(1,length(fake_data))
	group2_perms <- t(combn(i_data,length(treatment)))

	pick_other_index <- function(x){i_data[!(i_data %in% x)]}
	group1_perms <- adply(group2_perms, 1, pick_other_index)[,-1]
	total_num_perm <- nrow(group1_perms)


	### Function to get the kth permutation out of the total_num_perm
	group_labels <- c(rep("control", ncol(group1_perms)),
	                  rep("treatment", ncol(group2_perms)))
	get_k_perm_data <- function(k){
	  data.frame(group = group_labels,
	             value = c(as.numeric(group1_perms[k,]),
	                       as.numeric(group2_perms[k,]))
	  )
	}

	# Function to run the test for a given dataset,
	# this function is used to test each possible permuted dataset
	test_dataset <- function(dataset){
	  test <- t.test(value ~ group, data=dataset)
	  c(test$statistic, "p-value" = test$p.value)
	}

	# Run the Permutation Test
	run_kth_test <- function(k){
	  test_dataset(get_k_perm_data(k))
	}
	possible_perms <- 1:total_num_perm
	# Use a plyr function instead of a loop
	perm_test_dist <- aaply(possible_perms, 1, run_kth_test)

	# plot the p-values from one set of data
	# these are the t-statistics and p-values for each of the 35 permutations
	# 7 choose 4 = 35
	# hist(perm_test_dist[,"p-value"])

	# compute our p-value to be the fraction of t-statistics that are greater than
	  # the first observation
	# Calculate P-value for permutation test
	sum(perm_test_dist[,"p-value"] <= test_dataset(fake_dataset)["p-value"])/nrow(perm_test_dist)
}


run_sims <- function(n_sim, c_dist, t_dist){
	# rdply repeats the function n_sim times and then aggregates the results
	rdply(n_sim, run_perm_test(c_dist, t_dist), .parallel = TRUE)[,2]
}


# expression function allows random number generation to be done at each
# repetition of of the permutation test
# run_sims(n_sim = 100, c_dist = expression(rnorm(3,5,1)),
# 	                  t_dist = expression(rnorm(4,5,1)))
# run_sims(n_sim = 100, c_dist = expression(rnorm(3,5,1)),
# 	                  t_dist = expression(rnorm(4,7,1)))


### Power Studies
# Setup variables
set.seed(1140)
step <- c(1:10, seq(20,100,10))
n_simulations <- 1000
alpha <- .05
ps_sims <- list()
power_study_results <- data.frame(
	step = step,
	normal = numeric(length=length(step)),
	t = numeric(length=length(step)),
	exponential = numeric(length=length(step)),
	contaminated = numeric(length=length(step))
)

# Random number generator for contaminated normal
rcnorm <- function(n, o_mu, c_mu, c_pct, sd=1){
	o_n <- n*(1-c_pct)
	c_n <- n*(c_pct)
	# If we don't get an integer for contamination percent
	# randomly assign to either group
	if(o_n %% 1 !=0){
		ifelse(runif(1)>=.5,
			   o_n <- o_n + 1,
			   c_n <- c_n + 1)
	}

	c(rnorm(o_n, o_mu, sd), rnorm(c_n, c_mu, sd))
}

# Loop through steps
st_power <- proc.time()
for(i in step){
	res <- list()
	res$norm <- run_sims(n_sim = n_simulations,
		                 c_dist = expression(rnorm(3,5,1)),
	                     t_dist = expression(rnorm(4,5+i,1)))
	res$t <- run_sims(n_sim = n_simulations,
		              c_dist = expression(rt(3,df=3,ncp=5)),
	                  t_dist = expression(rt(4,df=3,ncp=5+i)))
	res$exp <- run_sims(n_sim = n_simulations,
                        c_dist = expression(rexp(3,1/5)),
                        t_dist = expression(rexp(4,1/(5+i))))
	res$cnrm <- run_sims(n_sim = n_simulations,
		 c_dist = expression(rcnorm(3,o_mu = 5, c_mu = 8, c_pct= .15)),
	     t_dist = expression(rcnorm(3,o_mu = 5+i, c_mu = 8+i, c_pct= .15)))

    # Calculate the number of p-values less than alpha
    power_study_results[i,-1] <- laply(res,function(x){sum(x<alpha)})
}
et_power <- proc.time()
et_power - st_power
save.image("Power_Study.RData")

# Generate Plot
results_plot <- melt(power_study_results, id.vars = c("step"))
names(results_plot)<- c("step","dist","cnt_sig")
results_plot$power <- results_plot$cnt_sig/n_simulations
power_p <- ggplot(results_plot, aes(x=step, y=power, fill=dist))
power_p <- power_p + geom_bar(stat="identity", position="dodge")
# power_p

# Save Data for further review
save.image("Power_Study.RData")
