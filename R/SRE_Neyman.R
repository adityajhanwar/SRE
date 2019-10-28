# @file: SRE_Neyman.R
# @description: Neymanian Point Estimation of Average Causal Effect for Stratified
#     and Post-Stratified Randomized Experiments

#' SRE and Post-SRE Neyman Analysis
#'
#' @param z Vector of treatment/control assignments
#' @param y Vector of reponses
#' @param x Vector of covariate to be blocked on
#' @return List containing Neyman point estimate of the Average Causal Effect
#'     and the respective conservative estimator
Neyman_SRE = function(z, y, x){

  stratums  = unique(x)
  pi.k      = c()
  tau.k     = c()
  var.tau.k = c()

  for(k in stratums){
    zk        = z[x == k]
    yk        = y[x == k]
    pi.k      = append(pi.k, length(zk)/length(z))
    tau.k     = append(tau.k, mean(yk[zk == 1]) - mean(yk[zk == 0]))
    var.tau.k = append(var.tau.k, var(yk[zk == 1])/sum(zk) + var(yk[zk == 0])/sum(1 - zk))
  }

  return(list('average.causal.effect' = sum(pi.k*tau.k),
              'variance.estimator'    = sum(pi.k^2*var.tau.k)))
}
