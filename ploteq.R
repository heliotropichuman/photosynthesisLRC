#' @name plot_eq
#'
#' @title Plot Photosynthetic Models
#'
#' @description This function plots the fit of a given photosynthetic light
#'     response equation for all or select species in a data set. This base R
#'     plot includes the SampleID and the equation number in the title.
#'
#' @usage
#'  ```R
#'     ploteq(
#'         eqX = eq1
#'         eq_name = "eq1"
#'         i,
#'         data,
#'         title = "_",
#'         species_subset = NULL)
#'  ```
#' @param eqX A function representing the photosynthetic light response equation
#'     (e.g., \code{eq1}, \code{eq2}).
#' @param eq_name A character string representing the name of the equation,
#'     to be included in the plot title.
#' @param i An integer specifying the index of the species in the \code{inds}
#'     vector.
#' @param data A data frame containing the experimental data with at least two
#'     columns: 'PARi' for the incident light and 'A' for the measured
#'     photosynthetic rate.
#' @param title An optional character string specifying the title of the plot
#'     (defaults to title in the format "i SampleID Equation X.
#' @param species_subset An optional vector of species names from \code{inds}
#'     to be plotted. If \code{NULL}, all species in \code{inds} will be used
#'     (default is \code{NULL}).
#'
#' @return A plot of the measured data points for the selected species (open
#'     points), with curve parameters from the fitted equation (black points),
#'     the NLS curve (red line), and the model fit (dashed blue line). It will
#'     also return the reconstructed model fit as a list.
#'
#' @details This function takes the equation of photosynthetic light response
#'     models and fits it to the data for a given species. It then plots the
#'     observed and predicted values, highlighting specific points on the curve
#'     (such as the model curve paramaters \code{I15}, \code{I25}, \code{I85},
#'     and \code{I95}), where the number (X) is the carbon assimilation rate at
#'     X percent of the maximum assimilation in the measured data. The equation
#'     name is included in the plot title, and an optional subset of
#'     species can be selected for plotting. The function also calculates
#'     various fit statistics and adds both the original and reconstructed
#'     predictions as curves to the plot.
#'
#' @examples
#' \dontrun{
#' ```R
#' # Example with eq1 and all species
#' for (i in 1:length(inds)) {
#'   ploteq(eq1, "eq1", i, data = my_obesrved_data, title = "Equation 1")
#'      }
#'
#' # Example of using the function for all equations with all species or
#'     a subset of species
#' highlight <- c("Species1", "Species2", "Species3")
#' par(mfrow = c(3, 3))
#' for (i in 1:length(highlight)) {
#'   # Add equation names to the function calls
#'   ploteq(eq1, "eq1", i, data = LRCdata, species_subset = highlight)
#'   ploteq(eq2, "eq2", i, data = LRCdata, species_subset = highlight)
#'   ploteq(eq3, "eq3", i, data = LRCdata, species_subset = highlight)
#'   ploteq(eq4, "eq4", i, data = LRCdata, species_subset = highlight)
#'   ploteq(eq5, "eq5", i, data = LRCdata, species_subset = highlight)
#'   ploteq(eq6, "eq6", i, data = LRCdata, species_subset = highlight)
#'   ploteq(eq8, "eq8", i, data = LRCdata, species_subset = highlight)
#'   ploteq(eq9, "eq9", i, data = LRCdata, species_subset = highlight)
#'   ploteq(eq11, "eq11", i, data = LRCdata, species_subset = highlight)
#'   dev.ofploteq()
#'   }
#'
#'}
#' ```
#'
#' @export

ploteq <- function(eqX, eq_name, i, data, title = "", species_subset = NULL) {
  # Optionally subset species from inds
  selected_inds <- if (!is.null(species_subset)) species_subset else inds

  # Fit measured data to the photosynthetic model equation
  fit <- eqX(data = data[data$SampleID == selected_inds[i],],
             return = "all")
  calc <- fit$calc

  # Placement of open points on the curve obtained as parameters for
  # I15, I25, I85, I95, and Pmax_obs
  x <- c(
    rep(calc[["I15"]], 2),
    rep(calc[["I25"]], 2),
    rep(calc[["I85"]], 2),
    rep(calc[["I95"]], 3),
    rep(2500, 4)
  )

  y <- c(
    rep(.15 * calc[["Pmax_obs"]], 2),
    rep(.25 * calc[["Pmax_obs"]], 2),
    rep(.85 * calc[["Pmax_obs"]], 2),
    rep(.95 * calc[["Pmax_obs"]], 3),
    rep(calc[["Pmax_obs"]], 4)
  )

  # Extract the equation number for a clean title
  eq_number <- sub("eq", "Equation ", eq_name)

  # Plot the data
  plot(data[data$SampleID == selected_inds[i], 2:3],
       main = paste(i, selected_inds[i], "\n", title, eq_number),
       cex.main = 0.9,
       xlim = c(0, 2500),
       ylim = c(-4, max(data[data$SampleID == selected_inds[i], 3],
                        na.rm = TRUE) + 2))

  points(x, y, pch = 19)

  if (any(is.na(x) | is.na(y))) {
    exclude <- which(is.na(x) | is.na(y))
    x <- x[-exclude]
    y <- y[-exclude]
  }

  eqX <- get(eq_name)  # Get the equation by its name
  recon <- eqX(dat = data.frame(PARi = x, A = y), return = "all")

  # Add lines to the plot
  points(PARi_fine, recon$pred_fine, type = 'l', col = 'red', lwd = 2)
  points(PARi_fine, fit$pred_fine, type = 'l', col = 'blue', lty = 2)

  return(recon)
}
