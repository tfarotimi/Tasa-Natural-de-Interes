#------------------------------------------------------------------------------
# File:        run.hlw.estimation.R
#
# Description: Runs the three stages of the HLW estimation using country-specific
#              inputs and returns output from each stage. Each country file
#              'run.hlw.XX.R' calls this file.
#------------------------------------------------------------------------------#
run.hlw.estimation <- function(log.output,
                               inflation,
                               real.interest.rate,
                               nominal.interest.rate,
                               covid.indicator,
                               a.r.constraint=NA,
                               b.y.constraint=NA,
                               g.pot.start.index,
                               use.kappa=FALSE,
                               kappa.inputs=NA,
                               fix.phi=NA,
                               xi.00.stage1=NA, P.00.stage1=NA,
                               xi.00.stage2=NA, P.00.stage2=NA,
                               xi.00.stage3=NA, P.00.stage3=NA,
                               run.se=TRUE,
                               sample.end) {


  # Check COVID-adjusted model flags
  if (!is.na(fix.phi) & !is.numeric(fix.phi)) {
      stop("fix.phi must be set as NA or a numeric value")
  }
  if ((use.kappa) & is.null(dim(kappa.inputs))) {
      stop("if use.kappa=TRUE, kappa.inputs must be specified")
  }


  # Running the stage 1 model
  out.stage1 <- rstar.stage1(log.output=log.output,
                             inflation=inflation,
                             covid.indicator=covid.indicator,
                             b.y.constraint=b.y.constraint,
                             sample.end=sample.end,
                             use.kappa=use.kappa,
                             kappa.inputs=kappa.inputs,
                             fix.phi=fix.phi,
                             xi.00.stage1=xi.00.stage1,
                             P.00.stage1=P.00.stage1)

  # Median unbiased estimate of lambda_g
  lambda.g <- median.unbiased.estimator.stage1(out.stage1$potential.smoothed)
  #assign out.stage1 initial parameters[param_num["a_r"] to lm_ar

  #   lambda.g <- 0.056
 xtol_values <- c(1.0e-08, 1.0e-07, 1.0e-06, 1.0e-05, 1.0e-04, 1.0e-03, 1.0e-02, 1.0e-01, 1.0)
 xtol_last <- c(1.0e-04,1.0e-03, 1.0e-02, 1.0e-01, 1.0)

 lambda.z <- 0
  # Running the stage 2 model
  while (lambda.z == 0 && length(xtol_values) > 0) {
    # Run the stage 2 model with the initial guess for lambda.z
    out.stage2 <- rstar.stage2(log.output=log.output,
                               inflation=inflation,
                               real.interest.rate=real.interest.rate,
                               covid.indicator=covid.indicator,
                               lambda.g=lambda.g,
                               a.r.constraint=a.r.constraint,
                               b.y.constraint=b.y.constraint,
                               g.pot.start.index=g.pot.start.index,
                               sample.end=sample.end,
                               use.kappa=use.kappa,
                               kappa.inputs=kappa.inputs,
                               fix.phi=fix.phi,
                               xi.00.stage2=xi.00.stage2,
                               P.00.stage2=P.00.stage2,
                               xtol_values=xtol_values)

    # Median unbiased estimate of lambda_z
    lambda.z <- median.unbiased.estimator.stage2(out.stage2$y, out.stage2$x, out.stage2$kappa.vec)

    xtol_values <- xtol_values[-1]
  }

  # If lambda.z piled at zero (Stock-Watson test couldn't reject for this sample),
  # use the last Stage 2 estimate and set a small positive floor to keep Stage 3 identified
  if (lambda.z == 0) {
    print("WARNING: lambda.z piled at zero — setting floor of 0.01 to keep Stage 3 identified")
    lambda.z <- 0.01
  }

  print(paste('xtol_values',xtol_values))
  # Running the stage 3 model
  out.stage3 <- rstar.stage3(log.output=log.output,
                            inflation=inflation,
                            real.interest.rate=real.interest.rate,
                            nominal.interest.rate=nominal.interest.rate,
                            covid.indicator=covid.indicator,
                            lambda.g=lambda.g,
                            lambda.z=lambda.z,
                            a.r.constraint=a.r.constraint,
                            b.y.constraint=b.y.constraint,
                            g.pot.start.index=g.pot.start.index,
                            run.se=run.se,
                            sample.end=sample.end,
                            use.kappa=use.kappa,
                            kappa.inputs=kappa.inputs,
                            fix.phi=fix.phi,
                            xi.00.stage3=xi.00.stage3,
                            P.00.stage3=P.00.stage3)

  # Return output from all three stages
  return(list(out.stage1=out.stage1,
              out.stage2=out.stage2,
              out.stage3=out.stage3,
              lambda.g=lambda.g,
              lambda.z=lambda.z))
}