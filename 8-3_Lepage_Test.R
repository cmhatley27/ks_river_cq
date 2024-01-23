
### LEPAGE TEST FUNCTIONS (FROM NSM3 PACKAGE)

lepage.test <- function(x, y = NA, g = NA, method = NA, n.mc = 10000){
  
  ##Adapted from kruskal.test()##
  if (is.list(x)) {
    if (length(x) < 2L) 
      stop("'x' must be a list with at least 2 elements")
    y <- x[[2]]
    x <- x[[1]]
  }
  else {
    if(min(is.na(y)) != 0){
      k <- length(unique(g))
      if (length(x) != length(g)) 
        stop("'x' and 'g' must have the same length")
      if (k < 2) 
        stop("all observations are in the same group")
      y <- x[g == 2]
      x <- x[g == 1]
    }
  }
  #####################
  
  outp <- list()
  outp$m <- m <- length(x)
  outp$n <- n <- length(y)
  outp$n.mc <- n.mc
  N <- outp$m + outp$n
  outp$ties <- (length(c(x, y)) != length(unique(c(x, y))))
  even <- (outp$m + outp$n + 1)%%2
  outp$stat.name <- "Lepage D"
  
  
  ##When the user doesn't give us any indication of which method to use, try to pick one.
  if(is.na(method)){
    if(choose(outp$m + outp$n, outp$n) <= 10000){
      method <- "Exact"
    }
    if(choose(outp$m + outp$n, outp$n) > 10000){
      method <- "Monte Carlo"
    }
  }
  #####################################################################
  outp$method <- method
  
  tmp.W <- rank(c(x, y))
  
  our.data <- rbind(c(x, y), c(rep(1, length(x)), rep(0, length(y))))
  sorted <- our.data[1, order(our.data[1, ]) ]
  x.labels <-our.data[2, order(our.data[1, ]) ]
  
  med <- ceiling(N / 2)
  if(even){no.ties <- c(1:med, med:1)}
  if(!even){no.ties <- c(1:med, (med - 1):1)}
  
  obs.group <- numeric(N)
  group.num <- 1
  for(i in 1:N){
    if(obs.group[i] == 0){
      obs.group[i] <- group.num
      for(j in i:N){
        if(sorted[i] == sorted[j]){
          obs.group[j] <- obs.group[i]
        }
      }
      group.num <- group.num + 1;
    }
  }
  
  group.ranks <- tapply(no.ties, obs.group, mean)
  
  tied.ranks <- numeric(N)
  for(i in 1:group.num){
    tied.ranks[which(obs.group == as.numeric(names(group.ranks)[i]))] <- group.ranks[i]
  }
  
  tmp.C <- c(tied.ranks[x.labels == 1], tied.ranks[x.labels == 0])
  
  ##Only needs to depend on y values
  D.calc <- function(C.vals, W.vals){
    
    if(even){
      exp.C <- n * (N + 2) / 4
      var.C <- m * n * (N + 2) * (N - 2) / (48 * (N - 1))
    }
    if(!even){
      exp.C <- n * (N + 1)^2 / (4 * N)
      var.C <- m * n * (N + 1) * (3 + N^2) / (48 * N^2)
    }
    W.obs <- sum(W.vals)
    W.star <- (W.obs - n * (N + 1) / 2) / sqrt(m * n * (N + 1) / 12)
    C.star <- (sum(C.vals) - exp.C) / sqrt(var.C)
    return(W.star^2 + C.star^2)
  }
  
  outp$obs.stat <- D.calc(tmp.C[(m + 1):N], tmp.W[(m + 1):N])
  
  if(outp$method == "Exact"){
    possible.orders <- gtools::combinations(outp$m + outp$n, outp$n)
    
    possible.C <- t(apply(possible.orders, 1, function(x) tmp.C[x]))
    possible.W <- t(apply(possible.orders, 1, function(x) tmp.W[x]))
    
    theor.dist <- numeric(nrow(possible.C))
    for(i in 1:nrow(possible.C)){
      theor.dist[i] <- D.calc(possible.C[i, ], possible.W[i, ])
    }
    
    outp$p.value <- mean(theor.dist >= outp$obs.stat)
  }
  
  if(outp$method == "Asymptotic"){
    outp$p.value <- (1 - pchisq(outp$obs.stat, 2))
  }
  
  if(outp$method == "Monte Carlo"){
    outp$p.value <- 0
    for(i in 1:n.mc){
      mc.sample <- sample(1:N, n)
      
      if(D.calc(tmp.C[mc.sample], tmp.W[mc.sample]) >= outp$obs.stat){
        outp$p.value = outp$p.value + 1 / n.mc
      }
    }
  }
  
  cat("\nLepage two-sample location-scale test\n")
  cat("\nNull hypothesis: The locations and scales of the two population distributions are equal.\n")
  cat("Alternative hypothesis: The locations and/or scales of the two population distributions differ.\n")
  cat(paste("\nD = ", round(outp$obs.stat, 3), ", p-value = ", round(outp$p.value, 4), "\n\n", sep=""))
  
  #class(outp)="NSM3Ch5p"
  return(outp)
}


combs <- combn(k, 2)

cluster_lep <- tibble(
  var = sort(unique(event_cdfs$var))
)

for(comb in seq(ncol(combs))){
  comb_lep <- event_cdfs %>%
    group_by(var) %>%
    summarize(L = lepage.test(value[cluster == combs[1,comb]], value[cluster == combs[2,comb]])$obs.stat,
              sig = lepage.test(value[cluster == combs[1,comb]], value[cluster == combs[2,comb]])$p.value)
  
  cluster_lep[,comb+1] <- comb_lep[,3]
  colnames(cluster_lep)[comb+1] <- paste0('sig_',combs[1,comb],'_',combs[2,comb])
}

#All significant pairwise differences
diff_vars_lep <- cluster_lep %>%
  pivot_longer(!var, names_to = 'combo', values_to = 'sig') %>%
  filter(sig <= 0.05)

vars_obs_lep <- diff_vars_lep %>%
  count(var) %>%
  arrange(desc(n)) %>%
  left_join(vars_obs, by = 'var', suffix = c('.lep', '.ks')) %>%
  left_join(vars_obs_cuc) %>%
  replace_na(list(n.ks = 0, n.lep = 0, n.cuc = 0))

method1 <- 'n.lep'
method2 <- 'n.cuc'
ggplot(data = vars_obs_lep) +
  geom_col(aes(x = fct_rev(fct_reorder(var, n.lep - n.cuc)), y = n.lep - n.cuc)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
