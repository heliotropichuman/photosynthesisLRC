#' @name eq11
#' @title Calculate Photosynthetic Rates Using a Nonlinear Model EQ11
#'
#' @description Uses the nonlinear least squares Ye model equation 11 from
#'     Lobo et. al (2013) to transform measured photosynthetic
#'     data into a smoothed function-valued trait with the following function:
#'       A~phi_I0_Icomp((1-beta(PARi))/(1+gamma(PARi)))(PARi-Icomp)
#'     The function will return predicted values, calculated quantities,
#'     or both.
#'
#' @param pars A named vector of parameters.
#'     Default values are phi_I0_Icomp = .0756, beta = .0000432, gamma = .0039,
#'     Icomp = 22.6.
#'     These serve as initial starting parameters for the function to rapidly
#'     assess your data through an iterative process. The empirical values
#'     beta and gamma range between 0 and 1, and are not explicitly described by
#'     Ye (2007) but are independent coefficients of I implemented to
#'     incorporate a more dynamic response to light. All of these values may be
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
#'     predicted_values <- eq11(pars = c(Pphi_I0_Icomp = .08,beta = .00004,
#'                        gamma = .0039,Icomp = 20.6),
#'                        PARi = c(0, 100, 200, 400, 800), return = "predict")
#'     print(predicted_values)
#'
#'     # Use experimental data to predict photosynthetic rates and estimate
#'     # linear parameters
#'     result <- eq11(data = example_data, return = "all")
#'     print(result$calc)  # View calculated quantities
#'     print(result$fit)   # View fit statistics and optimized parameters
#'
#'     # Get calculated quantities directly
#'     calculated_quantities <- eq11(data = example_data, return = "calc")
#'     print(calculated_quantities)
#'
#' @details The function uses the provided data to estimate the parameters
#'     phi_I0_Icomp, Icomp, and empirical parameters beta and gamma by
#'     minimizing the squared differences between observed and predicted
#'     photosynthetic rates. The model is then used to calculate a range of derived functional trait quantities such as the dark respiration rate (Rd), light compensation point (Icomp), maximum photosynthetic rate (Pmax), and curve derived parameters (Ix) among other calculated quantities.
#'
#' @references
#'     Lobo, F. de A., M. P. de Barros, H. J. Dalmagro,  .C. Dalmolin,
#'     W. E. Pereira, É.C. de Souza, G. L. Vourlitis and C. E. Rodriguez Ortiz
#'     2013 Fitting net photosynthetic light-response curves with Microsoft
#'     Excel – a critical look at the models. Photosynthetica 51 (3): 445-456.
#'
#'     Ye, Z.-P. 2007 A new model for relationship between irradiance and the
#'     rate of photosynthesis in *Oryza sativa*. Photosynthetica 45: 637-640.
#'
#'     Davis, R.E., C. M. Mason, E. W. Goolsby 2024 Comparative evolution of
#'     photosynthetic light response curve: approaches and pitfalls in
#'     phylogenetic modeling of a function-valued trait. IJPS
#'
#' @export
eq11 <- function(pars = c(phi_I0_Icomp = .0756,beta = .0000432,
                          gamma = .0039,Icomp = 22.6),
                 data,
                 PARi = c(0, 50, 100, 250, 500, 1000, 1500, 2000, 2500),
                 return = c("predict","calc","all")[1])
{
  mod <- NULL
  if(!missing(data))
  {
    try({
      mod <- nls(control = nls.control(maxiter = 1000),
              A ~ phi_I0_Icomp*((1-beta*PARi)/(1+gamma*PARi))*(PARi-Icomp),
              data = data,
              start = pars,
              lower = c(phi_I0_Icomp = .01,beta = 1e-5,gamma = 1e-5,Icomp = 5),
              upper = c(phi_I0_Icomp = .2,beta = 1,gamma = 1,Icomp = 150),
              algorithm = "port")
      pars <- coef(mod)
    },silent = TRUE)

    pars[["beta"]] <- log(pars[["beta"]])
    pars[["gamma"]] <- log(pars[["gamma"]])

    try({
      mod <- nls(control = nls.control(maxiter = 1000),
                A~phi_I0_Icomp*((1-exp(beta)*PARi)/
                                  (1+exp(gamma)*PARi))*(PARi-Icomp),
                data = data,start = pars)
      pars <- coef(mod)
    },silent = TRUE)

    PARi <- data$PARi
    if(!any(data$PARi==0))
    {
      I0 <- coef(lm(A~PARi,data = data[order(data$PARi)[1:2],]))[1]
      if(I0<0) pars[["Icomp"]] <- approx(x = c(I0,data$A),
                                         y = c(0,data$PARi),xout = 0)$y
    } else
    {
      I0 <- coef(lm(A~PARi,data = data[order(data$PARi)[1:2],]))[1]
      pars[["Icomp"]] <- approx(x = data$A,y = data$PARi,xout = 0)$y
    }

    pars[["phi_I0_Icomp"]] <- abs(I0) / pars[["Icomp"]]
    pars <- optim(par = pars,function(pars) sum(((data$A) -
            (pars[["phi_I0_Icomp"]]*((1-exp(pars[["beta"]])*data$PARi)/
            (1+exp(pars[["gamma"]])*data$PARi))*(data$PARi-pars[["Icomp"]])))^2),
            control=list(fnscale=1))$par
    try({ pars <- optim(par = pars,function(pars) sum(((data$A) -
          (pars[["phi_I0_Icomp"]]*((1-exp(pars[["beta"]])*data$PARi)/
            (1+exp(pars[["gamma"]])*data$PARi))*(data$PARi-pars[["Icomp"]])))^2),
              control=list(fnscale=1),method = "BFGS")$par },silent = TRUE)

    try({
      mod <- nls(control = nls.control(maxiter = 1000),
           A~phi_I0_Icomp*((1-exp(beta)*PARi)/(1+exp(gamma)*PARi))*(PARi-Icomp),
            data = data,start = pars)
      pars <- coef(mod)
    },silent = TRUE)

    pars[["beta"]] <- exp(pars[["beta"]])
    pars[["gamma"]] <- exp(pars[["gamma"]])

    try({
      mod <- nls(control = nls.control(maxiter = 1000),
                A~phi_I0_Icomp*((1-beta*PARi)/(1+gamma*PARi))*(PARi-Icomp),
                 data = data, start = pars, lower = c(phi_I0_Icomp = .01,
                                      beta = 1e-32,gamma = 1e-32,Icomp = 5),
                    upper = c(phi_I0_Icomp = .2,beta = 1,gamma = 1,Icomp = 150),
                     algorithm = "port")
      pars <- coef(mod)
    },silent = TRUE)
  }

  # parameters to optimize
  phi_I0_Icomp <- pars[["phi_I0_Icomp"]]
  Icomp <- pars[["Icomp"]]
  beta <- pars[["beta"]]
  gamma <- pars[["gamma"]]

  pred <- phi_I0_Icomp*((1-beta*PARi)/(1+gamma*PARi))*(PARi-Icomp)

  # predicted values
  if(return == "predict") return(pred)

  # calculated quantities
  x <- c(.05,.1,.15,.25,.5,.75,.85,.9,.95)
  PARi_fine <- seq(0,6000,by = .1)
  pred_fine <- phi_I0_Icomp*((1-beta*PARi_fine)/
                               (1+gamma*PARi_fine))*(PARi_fine-Icomp)
  Isat <- ((((beta+gamma)*(1+gamma*Icomp)/beta)^0.5)-1)/gamma
  Rd <- phi_I0_Icomp*Icomp
  Pgmax <- phi_I0_Icomp*((1-beta*Isat)/(1+gamma*Isat))*(Isat-Icomp)+Rd
  Isat_x <- (((phi_I0_Icomp*beta*Icomp)-(x*(Pgmax-Rd)*gamma)+phi_I0_Icomp)-
               (((x*(Pgmax-Rd)*gamma)-(phi_I0_Icomp*beta*Icomp)-phi_I0_Icomp)^2-
                  (4*phi_I0_Icomp*beta*(phi_I0_Icomp*Icomp+x*(Pgmax-Rd))))^0.5)/
                    (2*phi_I0_Icomp*beta)
  phi_Icomp <- phi_I0_Icomp*((1-2*beta*Icomp-beta*gamma*(Icomp^2)+
                                (gamma+beta)*Icomp)/(1+gamma*Icomp)^2)
  phi_I0 <- phi_I0_Icomp*((1-2*beta*0-beta*gamma*(0^2)+(gamma+beta)*0)/
                                  (1+gamma*0)^2)
  lo_inds <- which(PARi_fine <= (max(PARi_fine)-50))
  hi_inds <- lo_inds + (length(PARi_fine) - max(lo_inds))

  diffs <- (pred_fine[hi_inds] - pred_fine[lo_inds])
  if(any(diffs <= .067)) Imax <- PARi_fine[which(diffs <= .067)[1]]
  else Imax <- max(PARi_fine)

  Pmax <- max(pred_fine)
  P_Imax <- pred_fine[which.min(abs(Imax -PARi_fine))[1]]
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
            Imax_obs = 2500,Pmax_obs = pred_fine[PARi_fine == 2500],
            beta = beta,
            gamma = gamma)

  fit <- c(r2 = cor(data$A,pred,use = "pair")^2,
           mse = mean((data$A-pred)^2,na.rm = TRUE),list(pars = pars),
                  list(mod = mod))


  if(return == "calc") return(calc)

  return(list(calc = calc,pred = pred,pred_fine = pred_fine,fit = fit))
}
