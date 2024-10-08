% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eq9.R
\name{eq9}
\alias{eq9}
\title{Calculate Photosynthetic Rates Using a Nonlinear Model EQ9}
\usage{
# Example dataset
    example_data <- data.frame(
      PARi = c(0, 50, 100, 250, 500, 1000, 1500, 2000, 2500),
      A = c(1.8, 4.2, 7.5, 12.8, 16.2, 18.5, 19.3, 19.4, 19.5)
    )

    # Predict photosynthetic rates given the parameters
    predicted_values <- eq9(pars= c(Pgmax = 20,Icomp = 20,k = .0015,Rd = .2),
                       PARi = c(0, 100, 200, 400, 800), return = "predict")
    print(predicted_values)

    # Use experimental data to predict photosynthetic rates and estimate
    # linear parameters
    result <- eq9(data = example_data, return = "all")
    print(result$calc)  # View calculated quantities
    print(result$fit)   # View fit statistics and optimized parameters

    # Get calculated quantities directly
    calculated_quantities <- eq9(data = example_data, return = "calc")
    print(calculated_quantities)
}
\arguments{
\item{pars}{A named vector of parameters.
Default values are Pgmax = 22, Icomp = 10, k = .0015, and Rd = .1.
These serve as initial starting parameters for the function to rapidly
assess your data through an iterative process. The coefficient k
relates the photosynthetic capacity of a plant to the specific leaf mass,
and may increase with an increasing potential for photosynthetic activity
per unit leaf mass. All initial parameter values may be changed to fall
within the minimum and maximum parameter values of your study system.}

\item{data}{A data frame containing the experimental data with at least two
columns: 'PARi' for the incident light and 'A' for the photosynthetic
rate.}

\item{PARi}{A numeric vector of incident light values. Defaults to a sequence
from 0 to 2500.}

\item{return}{Character string indicating what the function should return.
Options are "predict" for predicted values, "calc" for calculated
quantities, and "all" for both. Defaults to "predict".}
}
\value{
Depending on the 'return' argument, the function returns:
  \itemize{
    \item \code{"predict"}: A numeric vector of predicted photosynthetic
                            rates.
    \item \code{"calc"}: A named vector of calculated quantities:
                         Pgmax,
                         Pmax,
                         Icomp,
                         phi_I0 (quantum yield calculated at I0),
                         phi_Icomp (quantum yield calculated at Icomp),
                         phi_I0_Icomp (quantum yield calculated by the range
                         of values between I0 and Icomp),
                         phi_Icomp_I200 (quantum yield calculated by the
                         range between Icomp and I200),
                         Rd (dark respiration),
                         Imax (Imax calculated),
                         Imax_obs (Imax observed),
                         P_Imax (assimilation value at maximum light),
                         Isat_x, x = .25, .50, .75, .85, .90, .95 (light
                         saturation at x percent of Pmax),
                         Ix, x = .25, .50, .75, .85, .90, .95 (light
                         intensity at x percent of Pmax)
    \item \code{"all"}: A list containing both the predicted values,
                        calculated quantities, and model fit statistics.
    }
}
\description{
Uses the nonlinear least squares exponential model
    equation 9 from Lobo et. al (2013) to transform measured photosynthetic
    data into a smoothed function-valued trait with the following function:
      A~Pgmax((1-exp(-k*(PARi-Icomp))))-Rd
    The function will return predicted values, calculated quantities,
    or both.
}
\details{
The function uses the provided data to estimate the parameters
    Pgmax, Icomp, Rd, and k by minimizing the squared differences between
    observed and predicted photosynthetic rates. The model is then used to
    calculate a range of derived functional trait quantities such as the dark
    respiration rate (Rd), light compensation point (Icomp), maximum
    photosynthetic rate (Pmax), and curve derived parameters (Ix) among other
    calculated quantities.
}
\references{
Lobo, F. de A., M. P. de Barros, H. J. Dalmagro,  .C. Dalmolin,
    W. E. Pereira, É.C. de Souza, G. L. Vourlitis and C. E. Rodriguez Ortiz
    2013 Fitting net photosynthetic light-response curves with Microsoft
    Excel – a critical look at the models. Photosynthetica 51 (3): 445-456.

    Prado, C. H. B. A., J. P. A. P. V. de Moraes 1997 Photosynthetic
    capacity and specific leaf mass in twenty woody species of cerrado
    vegetation under field conditions. Photosynthetica 33: 103-112.

    Davis, R.E., C. M. Mason, E. W. Goolsby 2024 Comparative evolution of
    photosynthetic light response curve: approaches and pitfalls in
    phylogenetic modeling of a function-valued trait. IJPS
}
