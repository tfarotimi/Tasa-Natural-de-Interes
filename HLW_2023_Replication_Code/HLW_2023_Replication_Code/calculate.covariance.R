#------------------------------------------------------------------------------#
# File:        calculate.covariance.R
#
# Description: This function calculates the covariance matrix of the
#              initial state from the gradients of the likelihood function.
#------------------------------------------------------------------------------#
calculate.covariance <- function(initial.parameters, theta.lb, theta.ub,
                                 y.data, x.data, stage,
                                 lambda.g=NA,lambda.z=NA, xi.00,
                                 use.kappa=FALSE, kappa.inputs=NA, param.num){

  # Number of state variables
  n.state.vars <- length(xi.00)

  # Set covariance matrix equal to 0.2 times the identity matrix
  P.00 <- diag(0.2,n.state.vars,n.state.vars)

  # Get parameter estimates via maximum likelihood
  f <- function(theta) {return(-log.likelihood.wrapper(theta, y.data=y.data, x.data=x.data, stage=stage,
                                                       lambda.g=lambda.g, lambda.z=lambda.z,
                                                       xi.00=xi.00, P.00=P.00,
                                                       use.kappa=use.kappa, kappa.inputs=kappa.inputs,
                                                       param.num=param.num)$ll.cum)}
  
  xtol_values <- c(1.0e-08, 1.0e-07, 1.0e-06, 1.0e-05, 1.0e-04, 1.0e-03, 1.0e-02, 1.0e-01, 1.0)
  for (xtol in xtol_values) {
    nloptr.out <- nloptr(initial.parameters, f, eval_grad_f=function(x) {gradient(f, x)},
                         lb=theta.lb, ub=theta.ub, opts=list("algorithm"="NLOPT_LD_LBFGS", "xtol_rel"=xtol, "maxeval"=5000))
  if (is.na(nloptr.out$solution[param.num["a_r"]])) {
    if (all(abs(nloptr.out$solution) >= 0.0001)) {
      print("no excessively small parameters")
      break
    }
  } else {
    if (all(abs(nloptr.out$solution) >= 0.0001) && 
        nloptr.out$solution[param.num["a_r"]] < -0.0025) {
      print("no excessively small parameters")
      break
    }
  }
  }
  theta <- nloptr.out$solution

  print(theta[param.num["sigma_ystar"]])

  if (nloptr.out$status==-1 | nloptr.out$status==5) {
      print(paste0("Look at the termination conditions for nloptr in calculate.covariance, Stage ",as.character(stage)))
      stop(nloptr.out$message)
  } else {
    print(paste0("Stage ",as.character(stage),", calculate.covariance: The terminal conditions in nloptr are"))
    print(nloptr.out$message)
  }

  print("running kalman filter - Stage", stage)
  # Run Kalman filter with above covariance matrix and corresponding parameter
  # estimates
  states <- kalman.states.wrapper(parameters=theta, y.data=y.data, x.data=x.data, stage=stage,
                                  lambda.g=lambda.g, lambda.z=lambda.z, xi.00=xi.00, P.00=P.00,
                                  use.kappa=use.kappa, kappa.inputs=kappa.inputs, param.num=param.num)

  # Save initial covariance matrix
  P.00 <- states$filtered$P.ttm1[1:n.state.vars,]

  return(P.00)
}