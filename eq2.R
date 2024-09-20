#' @name eq2
#' @title Calculate Photosynthetic Rates Using a Nonlinear Model EQ2
#'
#' @description Uses the nonlinear least squares Michaelis-Menton model
#'     equation 2 from Lobo et. al (2013) to transform measured photosynthetic
#'     data into a smoothed function-valued trait with the following function:
#'       A~((Pgmax*PARi)/(PARi+I50))-Rd
#'     The function will return predicted values, calculated quantities,
#'     or both.
#'
#' @param pars A named vector of parameters.
#'     Default values are Pgmax = 19.5,I50 = 216.4, and Rd = 1.8.
#'     These serve as initial starting parameters for the function to rapidly
#'     assess your data through an iterative process. These values may be
#'     changed to fall within the minimum and maximum parameter values of your
#'     study system.
#' @param data A data frame containing the experimental data with at least two
#'     columns: 'PARi' for the incident light and 'A' for the photosynthetic
#'     rate.
#' @param PARi A numeric vector of incident light values. Defaults to a sequence
#'    from 0 to 2500.
#' @param return Character string indicating what the function should return.
#'     Options are "predict" for predicted values, "calc" for calculated
#'     quantities, and "all" for both. Defaults to "predict".
#'
#' @return Depending on the 'return' argument, the function returns:
#'   \itemize{
#'     \item \code{"predict"}: A numeric vector of predicted photosynthetic
#'                             rates.
#'     \item \code{"calc"}: A named vector of calculated quantities:
#'                          Pgmax,
#'                          Pmax,
#'                          Icomp,
#'                          phi_I0 (quantum yield calculated at I0),
#'                          phi_Icomp (quantum yield calculated at Icomp),
#'                          phi_I0_Icomp (quantum yield calculated by the range
#'                          of values between I0 and Icomp),
#'                          phi_Icomp_I200 (quantum yield calculated by the
#'                          range between Icomp and I200),
#'                          Rd (dark respiration),
#'                          Imax (Imax calculated),
#'                          Imax_obs (Imax observed),
#'                          P_Imax (assimilation value at maximum light),
#'                          Isat_x, x = .25, .50, .75, .85, .90, .95 (light
#'                          saturation at x percent of Pmax),
#'                          Ix, x = .25, .50, .75, .85, .90, .95 (light
#'                          intensity at x percent of Pmax)
#'     \item \code{"all"}: A list containing both the predicted values,
#'                         calculated quantities, and model fit statistics.
#'     }
#'
#' @usage
#'     # Example dataset
#'     example_data <- data.frame(
#'       PARi = c(0, 50, 100, 250, 500, 1000, 1500, 2000, 2500),
#'       A = c(1.8, 4.2, 7.5, 12.8, 16.2, 18.5, 19.3, 19.4, 19.5)
#'     )
#'
#'     # Predict photosynthetic rates given the parameters
#'     predicted_values <- eq2(pars = c(Pgmax = 20, I50 = 216.4, Rd = 2),
#'                        PARi = c(0, 100, 200, 400, 800), return = "predict")
#'     print(predicted_values)
#'
#'     # Use experimental data to predict photosynthetic rates and estimate
#'     # linear parameters
#'     result <- eq2(data = example_data, return = "all")
#'     print(result$calc)  # View calculated quantities
#'     print(result$fit)   # View fit statistics and optimized parameters
#'
#'     # Get calculated quantities directly
#'     calculated_quantities <- eq2(data = example_data, return = "calc")
#'     print(calculated_quantities)
#'
#' @details The function uses the provided data to estimate the parameters
#'     Pgmax, assimilation at half the maximum assimilation rate (I50), and Rd
#'     by minimizing the squared differences between observed and predicted
#'     photosynthetic rates. The model is then used to calculate a range of
#'     derived functional trait quantities such as the dark respiration rate
#'     (Rd), light compensation point (Icomp), maximum photosynthetic rate
#'     (Pmax), and curve derived parameters (Ix) among other calculated
#'     quantities.
#'
#' @references
#'     Lobo, F. de A., M. P. de Barros, H. J. Dalmagro,  .C. Dalmolin,
#'     W. E. Pereira, É.C. de Souza, G. L. Vourlitis and C. E. Rodriguez Ortiz
#'     2013 Fitting net photosynthetic light-response curves with Microsoft
#'     Excel – a critical look at the models. Photosynthetica 51 (3): 445-456.
#'
#'     Kaipiainen, E.L. 2009 Parameters of photosynthesis light curve in
#'     *Salix dasyclados* and their changes during the growth season.
#'     Russ. J. Plant Physiol. 56: 445-453
#'
#'     Davis, R.E., C. M. Mason, E. W. Goolsby 2024 Comparative evolution of
#'     photosynthetic light response curve: approaches and pitfalls in
#'     phylogenetic modeling of a function-valued trait. IJPS
#'
#' @export
eq2 <- function(pars =
                  c(Pgmax = 19.5,I50 = 216.4,Rd = 1.8),
                     data,
                     PARi = c(0, 50, 100, 250, 500, 1000, 1500, 2000, 2500),
                     return = c("predict","calc","all")[1])
{
  mod <- NULL
  if(!missing(data))
  {
    PARi <- data$PARi
    if(any(PARi==0)) if(any(!is.na(data$A[PARi==0]))) pars[["Rd"]] <-
        abs(data$A[which(!is.na(data$A[PARi==0]))][1])
    pars[["Pgmax"]] <- max(data$A,na.rm = TRUE) + pars[["Rd"]]

    pars <- optim(par = pars,function(pars) sum(((data$A) -
              (((pars[["Pgmax"]]*data$PARi)/(data$PARi+pars[["I50"]]))-
                 pars[["Rd"]]))^2),control=list(fnscale=1))$par
    try({
      pars <- optim(par = pars,function(pars) sum(((data$A) -
                (((pars[["Pgmax"]]*data$PARi)/(data$PARi+pars[["I50"]]))-
                   pars[["Rd"]]))^2),control=list(fnscale=1))$par
    },silent = TRUE)
    try({
      mod <- nls(control = nls.control(maxiter = 1000),
                 A~((Pgmax*PARi)/(PARi+I50))-Rd,data = d,start = pars)
      pars <- coef(mod)
    },silent = TRUE)
  }
  # parameters to optimize
  Pgmax <- pars[["Pgmax"]]
  I50 <- pars[["I50"]]
  Rd <- pars[["Rd"]]

  # predicted values
  pred <- ((Pgmax*PARi)/(PARi+I50))-Rd
  if(return == "predict") return(pred)

  phi_I0 <- (I50*Pgmax)/((0+I50)^2)

  # calculated quantities
  Icomp <- (I50*Rd)/(Pgmax-Rd)
  x <- c(.05,.1,.15,.25,.5,.75,.85,.9,.95)
  PARi_fine <- seq(0,6000,by = .1)
  pred_fine <- ((Pgmax*PARi_fine)/(PARi_fine+I50))-Rd
  Isat_x <- (I50*(x*Rd-x*Pgmax-Rd))/(Pgmax*(x-1)+Rd*(1-x))
  lo_inds <- which(PARi_fine <= (max(PARi_fine)-50))
  hi_inds <- lo_inds + (length(PARi_fine) - max(lo_inds))

  diffs <- (pred_fine[hi_inds] - pred_fine[lo_inds])
  if(any(diffs <= .067)) Imax <- PARi_fine[which(diffs <= .067)[1]]
  else Imax <- max(PARi_fine)

  Pmax <- max(pred_fine)
  P_Imax <- pred_fine[which.min(abs(Imax -PARi_fine))[1]]
  phi_Icomp <-(I50*Pgmax)/((Icomp+I50)^2)
  I0_Icomp <- if(Icomp <= 0) I0_Icomp <- 1:2 else which(PARi_fine <= Icomp)
  Icomp_I200 <- which(PARi_fine >= Icomp & PARi_fine <= 200)

  phi_I0_Icomp <- as.numeric(coef(lm(pred_fine[I0_Icomp] ~
                                       PARi_fine[I0_Icomp]))[2])
  phi_Icomp_I200 <- as.numeric(coef(lm(pred_fine[Icomp_I200] ~
                                         PARi_fine[Icomp_I200]))[2])

  P2500 <- pred_fine[PARi_fine == 2500]

  Ix <- sapply(x,function(x){
    ind_x <- which((pred_fine - (x*P2500))>0 & PARi_fine <
                     PARi_fine[which.max(pred_fine[PARi_fine <= 2500])])[1]
    if(is.na(ind_x)) return(NA)
    ind_x <- c((ind_x-1):(ind_x+1))
    Ix_mod <- lm(y~x,data.frame(x=PARi_fine[ind_x],y=pred_fine[ind_x]))
    (((x*P2500)-coef(Ix_mod)[[1]])/coef(Ix_mod)[[2]])
  })

  calc <- c(Pgmax = Pgmax,
            phi_I0 = phi_I0,
            Rd = Rd,
            Icomp = Icomp,
            Isat_25 = Isat_x[x==.25],Isat_50 = Isat_x[x==.5],
            Isat_75 = Isat_x[x==.75],Isat_85 = Isat_x[x==.85],
            Isat_90 = Isat_x[x==.9],
            Isat_95 = Isat_x[x==.95],
            I5 = Ix[x==.05],I10 = Ix[x==.1],I15 = Ix[x==.15],
            I25 = Ix[x==.25],I50 = Ix[x==.5],
            I75 = Ix[x==.75],I85 = Ix[x==.85],
            I90 = Ix[x==.9],
            I95 = Ix[x==.95],
            Imax = Imax,
            Pmax = Pmax,
            phi_Icomp = phi_Icomp,
            phi_I0_Icomp = phi_I0_Icomp,
            phi_Icomp_I200 = phi_Icomp_I200,
            P_Imax = P_Imax,
            Imax_obs = 2500,Pmax_obs = pred_fine[PARi_fine == 2500])

  if(return == "calc") return(calc)
  fit <- c(r2 = cor(data$A,pred,use = "pair")^2,
           mse = mean((data$A-pred)^2,na.rm = TRUE),list(pars = pars),
                    list(mod = mod))

  return(list(calc = calc,pred = pred,pred_fine = pred_fine,fit = fit))
}
