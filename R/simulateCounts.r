#' Simulate RNA-seq counts via Tukey's G transformation and Negative Binomial distribution
#' 
#' @param mu base counts
#' @param groupSize vector containing each group size. Names of groupSize should be the group names
#' @param disp dispersion parameter for NB distribution
#' @param groupFC the fc for every gene in every group (use first group as reference)
#' @param chiDF degree of freedom for chisq gene-specific dispersion factor
#' @param techShapeBias a vector of length sum(groupSize) indicating g shape parameters
#' @param scaleBais a vector of length sum(groupSize) indicating scaling factors
#' @export
simulateCounts = function(mu, groupSize, disp=0.1, chiDF=40, groupFC=NULL, techShapeBias=NULL, scaleBias=NULL) {
  
  ngroups = length(groupSize)
  
  groupMu = matrix(rep(mu + 0.5, ngroups), ncol=ngroups)
  
  if (is.null(groupFC)) {
    
    groupFC = 1
    
  }else{
    
    if (!all(dim(groupFC) == dim(groupMu))) stop("K.Okrah: Check groupFC size.")
    
  }
  
  groupMu = groupMu * groupFC
  
  ugroups = names(groupSize)
  
  if (is.null(ugroups)) {
    
    ugroups = LETTERS[1:ngroups]
      
  }
  
  groups = rep(ugroups, groupSize)
  
  if (is.null(techShapeBias)) {
    
    techShapeBias = tapply(rep(0, sum(groupSize)), groups, identity)
    
  }else{
    
    if (!length(techShapeBias)==sum(groupSize)) stop("K.Okrah: Check length of techShapeBias.")
      
    techShapeBias = tapply(techShapeBias, groups, identity)
    
  }
  
  if (is.null(scaleBias)) {
    
    scaleBias = tapply(rep(1, sum(groupSize)), groups, identity)  
    
  }else{
    
    if (!length(scaleBias)==sum(groupSize)) stop("K.Okrah: Check length of scaleBias.")
    
    scaleBias = tapply(scaleBias, groups, identity) 
    
  }
  
  counts = c()
  
  for (k in 1:ncol(groupMu)) {
    
    for (j in 1:groupSize[k]) {
      
      z = log(groupMu[,k])
      
      mz = mean(z)
      
      sdz = sd(z)
      
      z = (z - mz) / sdz
      
      g = techShapeBias[[k]][j]
      
      if (g!=0) {
        z = (exp(z * g) - 1) / g
      }
      
      z = z * sdz + mz
      
      z = round(exp(z))
        
      z = z * scaleBias[[k]][j]
      
      DISP = dispf(z, a=sqrt(disp), chiDF=chiDF)
      
      counts = cbind(counts, rnbinom(length(z), mu=z, size = 1 / DISP))
      
    }
    
  }
  
  colnames(counts) = paste0(groups, ".", 1:sum(groupSize))
  
  rownames(counts) = paste0("gene.", 1:nrow(groupMu))
  
  list(counts=counts, groups=groups, groupMu=groupMu)
}
