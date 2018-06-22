nb <- 10000

# loading the parallel package:
library(parallel) # detectCores(), makeCluster(), clusterEvalQ(),
# clusterExport(), parSapply(), stopCluster()

# loading the required data:
result_k_means <- readRDS("result_k_means.rds")

# defining the function we need in the master workspace:
pairwise_distance <- function(df) {
  df %>%
    `[`(, c("longitude", "latitude")) %>%
    unlist() %>%
    matrix(, 2) %>%
    geosphere::distm() %>%
    {.[lower.tri(.)]}
}

# and same for the other function:
compute_statistic <- function(df) {
  df %>%
    split(df$group) %>%
    lapply(pairwise_distance) %>%
    unlist() %>%
    mean()
}

# create a cluster, with all the available cores but one:
cl <- makeCluster(detectCores() - 1)

# loading the needed packages on each core:
clusterEvalQ(cl, library(dplyr, geosphere))

# sending needed object from master workspace to each core's workspace:
clusterExport(cl, c("result_k_means", "pairwise_distance", "compute_statistic"))

# parallelized replication:
stat_null_distr <- parSapply(cl, 1:nb,
                             function(x) compute_statistic(mutate(result_k_means, group = sample(group))))

# stopping the cluster:
stopCluster(cl)

# saving the data to disk:
saveRDS(stat_null_distr, "output.rds")