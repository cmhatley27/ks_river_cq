cucconi.teststat <- function(x, y, m = length(x), n = length(y)){
  
  # Calculates the test statistic for the Cucconi two-sample location-scale test
  
  N <- m + n
  S <- rank(c(x, y))[(m + 1):N]
  denom <- sqrt(m * n * (N + 1) * (2 * N + 1) * (8 * N + 11) / 5)
  U <- (6 * sum(S^2) - n * (N + 1) * (2 * N + 1)) / denom
  V <- (6 * sum((N + 1 - S)^2) - n * (N + 1) * (2 * N + 1)) / denom
  rho <- (2 * (N^2 - 4)) / ((2 * N + 1) * (8 * N + 11)) - 1
  C <- (U^2 + V^2 - 2 * rho * U * V) / (2 * (1 - rho^2))
  return(C)
}

cucconi.dist.perm <- function(x, y, reps = 1000){
  
  # Computes the distribution of the Cucconi test statistic using random permutations
  
  m <- length(x)
  n <- length(y)
  N <- m + n
  alldata <- c(x, y)
  
  bootFunc <- function(){
    permdata <- alldata[sample(1:N, size = N, replace = FALSE)]
    xperm <- permdata[1:m]
    yperm <- permdata[(m + 1):N]
    return(cucconi.teststat(x = xperm, y = yperm, m = m, n = n))    
  }
  permvals <- replicate(reps, expr = bootFunc())
  
  #  permvals<-rep(NA,times=reps)
  #  for(r in 1:reps){
  #    permdata<-alldata[sample(1:N,size=N,replace=FALSE)]
  #    xperm<-permdata[1:m]
  #    yperm<-permdata[(m+1):N]
  #    permvals[r]<-cucconi.teststat(x=xperm,y=yperm, m=m, n=n)
  #  }
  return(permvals)
}

cucconi.test <- function(x, y, method = c("permutation", "bootstrap")){
  
  # Implementation of the Cucconi test for the two-sample location-scale problem
  # A permutation/bootstrap distribution of the test statistic (C) under the
  # null hypothesis is used to calculate the p-value.
  # Reference: Marozzi (2013), p. 1302-1303
  
  m <- length(x)
  n <- length(y)
  C <- cucconi.teststat(x = x, y = y, m = m, n = n)
  
  if(method[1] == "permutation"){
    h0dist <- cucconi.dist.perm(x = x, y = y)
  }
  
  if(method[1] == "bootstrap"){
    h0dist <- cucconi.dist.boot(x = x, y = y)
  }
  
  p.value <- length(h0dist[h0dist >= C]) / length(h0dist)
  
  cat("\nCucconi two-sample location-scale test\n")
  cat("\nNull hypothesis: The locations and scales of the two population distributions are equal.\n")
  cat("Alternative hypothesis: The locations and/or scales of the two population distributions differ.\n")
  cat(paste("\nC = ", round(C, 3), ", p-value = ", round(p.value, 4), "\n\n", sep=""))
  
  return(list(C = C,
              method = method[1],
              p.value = p.value))
}

c_dist <- cucconi.dist.perm(event_summary_cluster$reservoir_precip_ratio[event_summary_cluster$cluster == 5], event_summary_cluster$reservoir_precip_ratio[event_summary_cluster$cluster == 6])
c_stat <- cucconi.teststat(event_summary_cluster$reservoir_precip_ratio[event_summary_cluster$cluster == 5], event_summary_cluster$reservoir_precip_ratio[event_summary_cluster$cluster == 6])
ggplot() +
  geom_histogram(data = tibble(c_dist), aes(x = c_dist)) +
  geom_vline(xintercept = c_stat)

combs <- combn(k, 2)

cluster_cuc <- tibble(
  var = sort(unique(event_cdfs$var))
)

for(comb in seq(ncol(combs))){
  comb_cuc <- event_cdfs %>%
    group_by(var) %>%
    summarize(ks = cucconi.test(value[cluster == combs[1,comb]], value[cluster == combs[2,comb]])$C,
              sig = cucconi.test(value[cluster == combs[1,comb]], value[cluster == combs[2,comb]])$p.value)
  
  cluster_cuc[,comb+1] <- comb_cuc[,3]
  colnames(cluster_cuc)[comb+1] <- paste0('sig_',combs[1,comb],'_',combs[2,comb])
}

#All significant pairwise differences
diff_vars_cuc <- cluster_cuc %>%
  pivot_longer(!var, names_to = 'combo', values_to = 'sig') %>%
  filter(sig <= 0.05)

vars_obs_cuc <- diff_vars_cuc %>%
  count(var) %>%
  arrange(desc(n)) %>%
  left_join(vars_obs, by = 'var', suffix = c('.cuc', '.ks')) %>%
  replace_na(list(n.ks = 0))

ggplot(data = vars_obs_cuc) +
  geom_col(aes(x = fct_rev(fct_reorder(var, n.cuc - n.ks)), y = n.cuc - n.ks)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


