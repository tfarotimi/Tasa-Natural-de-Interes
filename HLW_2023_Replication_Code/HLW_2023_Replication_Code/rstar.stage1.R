#------------------------------------------------------------------------------#
# File:        rstar.stage1.R
#
# Description: This file runs the model in the first stage of the HLW estimation.
#------------------------------------------------------------------------------#
rstar.stage1 <- function(log.output,
                         inflation,
                         covid.indicator,
                         b.y.constraint=NA,
                         sample.end,
                         use.kappa=FALSE,
                         kappa.inputs=NA,
                         fix.phi=NA,
                         xi.00.stage1=NA,
                         P.00.stage1=NA) {

  stage <- 1

  #----------------------------------------------------------------------------#
  # Obtain initial parameter values
  #----------------------------------------------------------------------------#

  # Data must start 4 quarters before the estimation period
  t.end <- length(log.output) - 4

  # Estimate log-linear trend in GDP via OLS
  x.og <- cbind(rep(1,t.end+4), 1:(t.end+4))
  y.og <- log.output
  output.gap <- (y.og - x.og %*% solve(t(x.og) %*% x.og, t(x.og) %*% y.og)) * 100


  # IS curve, estimated by nonlinear LS:
  y.is     <- (output.gap[5:(t.end+4)])
  y.is.l1  <- (output.gap[4:(t.end+3)])
  y.is.l2  <- (output.gap[3:(t.end+2)])
  d        <- covid.indicator[5:(t.end+4)]
  d.l1     <- covid.indicator[4:(t.end+3)]
  d.l2     <- covid.indicator[3:(t.end+2)]

  # If the start date is 2020:Q1 or later, run the model with COVID dummy included
  # Otherwise, initialize phi at zero
  if (is.na(fix.phi) & sample.end[1] >= 2020) {
    nls.is   <- nls(y.is ~ phi*d + a_1*(y.is.l1 - phi*d.l1) + a_2*(y.is.l2 - phi*d.l2),
                  start = list(phi=0,a_1=0,a_2=0))
    b.is     <- coef(nls.is)
  } else {
    nls.is   <- nls(y.is ~ a_1*(y.is.l1) + a_2*(y.is.l2),
                    start = list(a_1=0,a_2=0))
    b.is     <- coef(nls.is)
    b.is["phi"] <- 0
  }

  # Initialize phi at fix.phi, if not NA
  if (!is.na(fix.phi)) {
    b.is["phi"] <- fix.phi
  }

  r.is     <- y.is - predict(nls.is) # residuals
  s.is     <- sqrt(sum(r.is^2) / (length(r.is)-(length(b.is))))

  # Initialize kappa vector: 1 in all periods
  kappa.is.vec <- rep(1, length(r.is))


  # Phillips curve
  # Estimate by LS, modified to include the covid dummy:
  # pi(t) = B(L)pi(t-1) + b_y*(y(t-1)-q(t-1)-phi*d(t-1)) + eps_pi(t)
  y.ph <- inflation[5:(t.end+4)]
  x.ph <- cbind(inflation[4:(t.end+3)],
                (inflation[3:(t.end+2)]+inflation[2:(t.end+1)]+inflation[1:t.end])/3,
                y.is.l1 - b.is["phi"]*d.l1)

  b.ph <- solve(t(x.ph) %*% x.ph, t(x.ph) %*% y.ph)
  r.ph <- y.ph - x.ph %*% b.ph
  s.ph <- sqrt(sum(r.ph^2) / (length(r.ph)-(length(b.ph))))

  # Initialize kappa vector: 1 in all periods
  kappa.ph.vec <- rep(1, length(r.ph))


  # Initial parameter values
  initial.parameters <- c(b.is["a_1"], b.is["a_2"], b.ph[1], b.ph[3], 0.85, s.is, s.ph, 0.5, b.is["phi"])
  # Set parameter labels
  param.num <- c("a_y1"=1,"a_y2"=2,"b_pi"=3,"b_y"=4,"g"=5,"sigma_ygap"=6,"sigma_pi"=7,"sigma_ystar"=8,"phi"=9)
  n.params <- length(initial.parameters)


  if (any(!is.na(xi.00.stage1))) {
    print('Stage 1: Using xi.00 input')
    xi.00 <- xi.00.stage1
  } else {
    print('Stage 1: xi.00 from HP trend in log output')
    # Initialization of state vector for Kalman filter using HP trend of log output
    log.output.hp.trend <- hpfilter(log.output,freq=36000,type="lambda",drift=FALSE)$trend
    g.pot <- log.output.hp.trend[(g.pot.start.index):length(log.output.hp.trend)]
    xi.00 <- c(100*g.pot[3:1])
  }


  #----------------------------------------------------------------------------#
  # Build data matrices
  #----------------------------------------------------------------------------#
  y.data <- cbind(100 * log.output[5:(t.end+4)],
                  inflation[5:(t.end+4)])
  x.data <- cbind(100 * log.output[4:(t.end+3)],
                  100 * log.output[3:(t.end+2)],
                  inflation[4:(t.end+3)],
                  (inflation[3:(t.end+2)]+inflation[2:(t.end+1)]+inflation[1:t.end])/3,
                  covid.indicator[5:(t.end+4)],
                  covid.indicator[4:(t.end+3)],
                  covid.indicator[3:(t.end+2)])


  #----------------------------------------------------------------------------#
  # Estimate parameters using maximum likelihood
  #----------------------------------------------------------------------------#

  # Set an upper and lower bound on the parameter vectors:
  # The vector is unbounded unless values are otherwise specified
  theta.lb <- c(rep(-Inf,length(initial.parameters)))
  theta.ub <- c(rep(Inf,length(initial.parameters)))

  # Set a lower bound for the Phillips curve slope (b_y) of b.y.constraint, if not NA
  # In HLW, b.y.constraint = 0.025
  if (!is.na(b.y.constraint)) {
      if (initial.parameters[param.num["b_y"]] < b.y.constraint) {
          initial.parameters[param.num["b_y"]] <- b.y.constraint
      }
      theta.lb[param.num["b_y"]] <- b.y.constraint
  }

  # Set an upper and lower bound on phi parameter to fix phi, if fix.phi is not NA
  if (!is.na(fix.phi)) {
      print(paste0("Fixing phi at ",as.character(fix.phi)))
      theta.ub[param.num["phi"]] <- fix.phi
      theta.lb[param.num["phi"]] <- fix.phi
  }

  # Define Kappa
  if (use.kappa) {

    n.kappa <- dim(kappa.inputs)[1]

    for (k in 1:n.kappa) {
      theta.ind <- n.params + k
      kappa.inputs$theta.index[k] <- theta.ind # Store index position within theta in kappa.inputs
      param.num[kappa.inputs$name[k]] <- theta.ind # Store index position within theta in param.num

      initial.parameters[theta.ind] <- kappa.inputs$init[k]
      theta.lb[theta.ind] <- kappa.inputs$lower.bound[k] # default is 1 in kappa.inputs
      theta.ub[theta.ind] <- kappa.inputs$upper.bound[k] # default is Inf in kappa.inputs

      # Print statement for kappa value
      if (kappa.inputs$lower.bound[k]==kappa.inputs$upper.bound[k]) {
        print(paste0("Fixing ",kappa.inputs$name[k]," at ",as.character(kappa.inputs$lower.bound[k])))
      } else {
        print(paste0("Initializing ",kappa.inputs$name[k]," at ",as.character(kappa.inputs$init[k])))
      }

    }
  }


  # Set the initial covariance matrix (see HLW 2017, footnote 6)
  if (any(is.na(P.00.stage1))) {
    print('Stage 1: Initializing covariance matrix')
    P.00 <- calculate.covariance(initial.parameters=initial.parameters, theta.lb=theta.lb, theta.ub=theta.ub,
                                 y.data=y.data, x.data=x.data, stage=stage, lambda.g=NA, lambda.z=NA, xi.00=xi.00,
                                 use.kappa=use.kappa, kappa.inputs=kappa.inputs, param.num=param.num)
  } else {
    print('Stage 1: Using P.00 input')
    P.00 <- P.00.stage1
  }

  # Get parameter estimates via maximum likelihood
  f <- function(theta) {return(-log.likelihood.wrapper(parameters=theta, y.data=y.data, x.data=x.data, stage=stage,
                                                       lambda.g=NA, lambda.z=NA, xi.00=xi.00, P.00=P.00,
                                                       use.kappa=use.kappa, kappa.inputs=kappa.inputs,
                                                       param.num=param.num)$ll.cum)}
  nloptr.out <- nloptr(initial.parameters, f, eval_grad_f=function(x) {gradient(f, x)},
                       lb=theta.lb,ub=theta.ub,
                       opts=list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-8,"maxeval"=5000))
  theta <- nloptr.out$solution

  if (nloptr.out$status==-1 | nloptr.out$status==5) {
      print("Look at the termination conditions for nloptr in Stage 1")
      stop(nloptr.out$message)
  } else {
    print("Stage 1: The terminal conditions in nloptr are")
    print(nloptr.out$message)
  }

  # Kalman log likelihood value
  log.likelihood <- log.likelihood.wrapper(parameters=theta, y.data=y.data, x.data=x.data, stage=stage,
                                           lambda.g=NA, lambda.z=NA, xi.00=xi.00, P.00=P.00,
                                           use.kappa=use.kappa, kappa.inputs=kappa.inputs, param.num=param.num)$ll.cum

  #----------------------------------------------------------------------------#
  # Kalman filtering
  #----------------------------------------------------------------------------#

  # Get state vectors (xi.tt, xi.ttm1, xi.tT, P.tt, P.ttm1, P.tT) via Kalman filter
  states <- kalman.states.wrapper(parameters=theta, y.data=y.data, x.data=x.data, stage=stage,
                                  lambda.g=NA, lambda.z=NA, xi.00=xi.00, P.00=P.00,
                                  use.kappa=use.kappa, kappa.inputs=kappa.inputs, param.num=param.num)

  matrices <- unpack.parameters.stage1(parameters=theta, y.data=y.data, x.data=x.data,
                                       xi.00=xi.00, P.00=P.00,
                                       use.kappa=use.kappa, kappa.inputs=kappa.inputs, param.num=param.num)

  # One-sided (filtered) estimates
  potential.filtered  <- states$filtered$xi.tt[,1]/100 # y*_t ; not COVID-adjusted
  output.gap.filtered <- y.data[,1] - (potential.filtered * 100) - theta[param.num["phi"]]*covid.indicator[5:(t.end+4)] # COVID-adjusted output gap

  # Two-sided (smoothed) estimates
  potential.smoothed  <- as.vector(states$smoothed$xi.tT[,1])/100 # y*_t ; not COVID-adjusted
  output.gap.smoothed <- y.data[,1] - (potential.smoothed * 100) - theta[param.num["phi"]]*covid.indicator[5:(t.end+4)] # COVID-adjusted output gap

  # Save variables to return
  return.list                <- list()
  return.list$theta          <- theta
  return.list$param.num      <- param.num
  return.list$matrices       <- matrices
  return.list$log.likelihood <- log.likelihood
  return.list$states         <- states
  return.list$xi.00          <- xi.00
  return.list$P.00           <- P.00
  return.list$y.data              <- y.data
  return.list$x.data              <- x.data
  return.list$initial.parameters  <- initial.parameters
  return.list$kappa.inputs        <- kappa.inputs
  return.list$potential.filtered  <- potential.filtered
  return.list$output.gap.filtered <- output.gap.filtered
  return.list$potential.smoothed  <- potential.smoothed
  return.list$output.gap.smoothed <- output.gap.smoothed
  return(return.list)
}
