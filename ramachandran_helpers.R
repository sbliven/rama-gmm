require(ggplot2)
require(magrittr)
require(dplyr)
require(mixtools)  #GNN

# Transform phi/psi to ramachandran-friendly ranges
wrap.psi = function(psi) {
  psi + 360*(psi < -100)
}
wrap.phi = function(phi) {
  phi + 360*(phi < 0)
}

# Save a GMM model from mixtools to TSV
saveGMM = function(gmm, path) {
  #saveRDS(stripGMM(gmm), path, ascii=T, compress=F)
  sigma = gmm$sigma %>%
    sapply(as.numeric) %>%
    t %>% as.data.frame()
  names(sigma) = paste0("sigma",1:ncol(sigma))
  mu = gmm$mu %>%
    sapply(as.numeric) %>%
    t %>% as.data.frame()
  names(mu) = paste0("mu", 1:ncol(mu))
  
  model = cbind(data.frame(name=gmm$name,
                           lambda=gmm$lambda),
                mu,
                sigma)
  
  write_tsv(model, path)
}
# Laod a GMM model from TSV to mixtools-lite (no data or predictions)
loadGMM = function(path) {
  model = read_tsv(path, col_types=cols())
  
  name = model %>%
    select(grep("name", names(.)))
  if(!(nrow(name) > 0)) {
    name = NULL
  }
  lambda = model %>%
    select(grep("lambda", names(.)))
  mu = model %>%
    select(grep("mu", names(.))) %>%
    apply(1,list) %>%
    unlist(recursive = F)
  sigma = model %>%
    select(grep("sigma",names(.))) %>%
    apply(1,list) %>%
    lapply(function(l){matrix(l[[1]], 2,2)})
  
  list(name=name[[1]],
       lambda=lambda[[1]],
       mu=mu,
       sigma=sigma)
}

# Reduce a GMM model to just the parameters to save space
stripGMM = function(gmm) {
  list(lamda=gmm$lambda,
       mu=gmm$mu,
       sigma=gmm$sigma
  )
}

# Display a GMM nicely
showClasses = function(x, gmm, name=NULL) {
  # x: Data frame with 'phi' and 'psi' elements
  # gmm: Simple GMM-like object. If `predictions` are present and the same length as s, use them for coloring
  
  # ellipsis data from gmm
  ell <- cbind(data.frame(Class=factor(rep(gmm$name, each=250))), 
               do.call(rbind, mapply(ellipse, gmm$mu, gmm$sigma, alpha=.1,
                                     npoints=250, draw=FALSE, SIMPLIFY=FALSE)))
  # Ellipse centers
  centers = data.frame(cbind(do.call(rbind, gmm$mu), gmm$lambda)) %>%
    set_names(c("phi","psi","lambda")) %>%
    mutate(Name=gmm$name, Class=gmm$name)
  
  if( !is.null(gmm$posterior) && nrow(gmm$posterior) == nrow(x)) {
    posterior_class = factor(apply(gmm$posterior, 1, which.max),
                             levels=1:nlevels(ell$Class),
                             labels=levels(ell$Class))
  } else {
    posterior_class = factor(rep(NA, nrow(x)),
                             levels=1:nlevels(ell$Class),
                             labels=levels(ell$Class))
  }
  
  x %>% as.data.frame() %>%
    mutate(Class=posterior_class) %>%
    ggplot(aes(x=phi, y=psi, color=Class)) +
    xlim(-10,370) + ylim(-110,270) +
    scale_colour_discrete(labels = parse(text=gmm$name), breaks=gmm$name) +
    # Points
    geom_point(alpha=.5, size=.5) +
    # Learned ellipsoids
    geom_path(data=ell, aes(x=`1`, y=`2`, color=Class, group=Class), color="black", na.rm=TRUE) +
    # Learned mu
    geom_point(data=centers,
               aes(x=phi,y=psi, size=lambda),
               color="grey",
               alpha=.8, shape='+') +
    scale_size_continuous(range=c(0,20)) +
    # Mu text
    geom_text(data=centers,
              aes(label=Name),
              #hjust="center",
              hjust="left", nudge_x=5,
              color="black",
              parse=TRUE)
  
}

# Computes posterior probabilities for x from the given Multivariate Normal Mixture
mvnormalmix.posterior = function(x, gmm) {
  posterior <- mapply(function(n, l, mu, sigma, x=x) { l * dmvnorm(x, mu, sigma) },
         gmm$name, gmm$lambda, gmm$mu, gmm$sigma,
         MoreArgs = list(x=as.matrix(x)))
  posterior/rowSums(posterior)
}

# Probability density for each point
dmvnormmix = function(x, gmm) {
  posterior <- mapply(function(n, l, mu, sigma, x=x) { l * dmvnorm(x, mu, sigma) },
                      gmm$name, gmm$lambda, gmm$mu, gmm$sigma,
                      MoreArgs = list(x=as.matrix(x)))
  rowSums(posterior)
}

rmvnormmix.var = function(n, gmm) {
  k = length(gmm$lambda)
  p = length(gmm$mu[[1]])
  z = sample(k, n, replace=T, prob=gmm$lambda)
  x = matrix(nrow=n, ncol=p)
  for(i in 1:k) {
    x[z==i,] = rmvnorm(sum(z==i), mu=gmm$mu[[i]], sigma=gmm$sigma[[i]])
  }
  x
}
#rmvnormmix.var(10, gmm.full)


