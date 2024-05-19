# For the discovery step, we need to find the discovery cutoff. Previously, we use grid search.
# Next, we consider the bisection search to improve the efficiency
MinMaxP.discov2 = function (P.max, P.min, method = c("FWER", "FDR"), alpha = 0.05)
{
  d = length(P.min)
  M = min(P.max)

  if (method == "FDR") {
    critFunc = function(x) {
      SAT:::Fmax_null(t = x, P.min = P.min)/(1/d + mean(P.max < x))
    }
  }else if(method == "FWER"){
    citeria_fun = function(x){ Cmax_null(t = x, P.min = P.min) * sum(P.min < x)  }
  }else{
    stop("Invalid method specified")
  }

  # we first narrow down the searching interval
  first_cutoff <- NA
  cutoff.choice.0 = sort( c( exp(log(M) * c(1:10) * 0.1), 1) )
  for(k in 1: length(cutoff.choice.0) ){
    criteria_value = critFunc( cutoff.choice.0[k] )
    if( criteria_value > alpha){
      first_cutoff <- cutoff.choice.0[k]
      break
    }
  }

  if (k == 1) {
    return(list(S = NULL, tau = M))
  }
  if (is.na(first_cutoff)) {
    print("No suitable cutoff found where criteria function exceeds alpha")
    return(list(S = NULL, tau = M))
  }

  # Bisection search to find the root for largest p such that c(p) <= alpha
  cutoff.upper = cutoff.choice.0[k]
  cutoff.lower = max(P.max[P.max <= cutoff.choice.0[k - 1]]) # we have to guarrantte that there exists a p-value achieving the lower
  cutoff_values <- sort(unique(P.max[P.max >= cutoff.lower & P.max <= cutoff.upper]))
  low <- 1; high <- length(cutoff_values)
  largest_cutoff <- NA  # Initialize largest_cutoff

  if(critFunc(cutoff_values[high]) <= alpha){
    largest_cutoff <- cutoff_values[high]
  }else{
    # Bisection search loop
    while (low <= high) {
      mid <- low + (high - low) %/% 2  # Find the middle index
      # Check the criteria function at the midpoint
      if (critFunc(cutoff_values[mid]) <= alpha) {
        largest_cutoff <- cutoff_values[mid]  # Update largest_cutoff if condition is met
        low <- mid + 1  # Move to the right half of the array
      } else{
        high <- mid - 1  # Move to the left half of the array
      }
    }
  }

  if(F){
    # due to the non-monotone of FDP, grid search for mid + 2, add one more BS
    low =  mid+2
    high = length(cutoff_values)
    if( low <= high & (critFunc(cutoff_values[low]) <= alpha) ){
      while (low <= high) {
        mid <- low + (high - low) %/% 2
        if (critFunc(cutoff_values[mid]) <= alpha){
          largest_cutoff <- cutoff_values[mid]  # Update largest_cutoff if condition is met
          low <- mid + 1  # Move to the right half of the array
        }else{
          high <- mid - 1  # Move to the left half of the array
        }
      }
    }
  }

  # sapply(cutoff_values,  critFunc)
  #sapply(cutoff.choice.0,  critFunc)

  tau = largest_cutoff
  S = which(P.max <= tau)
  return(list(S = S, tau = tau))
}
