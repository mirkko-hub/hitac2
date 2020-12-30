

# Models ------------------------------------------------------------------


# Linear model ------------------------------------------------------------

#' group_model
#'
#' This functions fits a linear model to 'Activity' as function of 'Time'
#' @param df dataframe with minimum columns 'Activity' and 'Time'
#'
#' @export
#' @examples
#' group_model()

group_model <- function(df){
  lm(Activity ~ Time, data = df)
}


# Y=Vmax*S/(Km + S) - simple Michaelis menten model -----------------------

#' model.nls1
#'
#' This functions fits a Michaelis-Menten curve (non-linear least squares) to initial velocities
#' @param df
#'
#' @export
#' @examples
#' model.nls1()

model.nls1 <- function(df){
  nls(v ~ Vm * S/(K+S), data = df,
      start = list(K = max(df$v)/2, Vm = max(df$v)))
}

# Y=Vmax*S/(Km + S*(1+S/Ki)) - simple inhibition model --------------------

#' model.nls2
#'
#' This functions fits a simple inhibition model to initial velocities
#' @param df
#'
#' @export
#' @examples
#' model.nls2()

model.nls2 <- function(df){
  nlsLM(v ~ ((Vm * S) / (Km + S * (1+S/Ki))), data = df,
        start = list(Km = max(df$v)/2, Vm = max(df$v), Ki = mean(filter(df, S == max(S))$v)))
}

# changed nls1 model ------------------------------------------------------

#' model.nls1_test
#'
#' This functions fits a Michaelis-Menten curve (non-linear least squares) to initial velocities
#' @param df
#'
#' @export
#' @examples
#' model.nls1_test()

model.nls1_test <- function(df){
  nls(v ~ Vm * S/(Km+S), data = df,
      start = list(Km = max(df$S)/2, Vm = max(df$v)))
}



# Fit lm per concentration ------------------------------------------------


#' hitachi_extract_lm
#'
#' This function fits linear models for each concentration and extracts model coefficients and statistics (glance, tidy, augment)
#' @param data tibble containing minimum columns 'Time', 'Concentration', 'Activity'
#'
#' @export
#' @examples
#' hitachi_extract_lm()

hitachi_extract_lm <- function (data)
{
  data <- data %>%
    dplyr::select(ExpID, Concentration, Time, Activity) %>%
    dplyr::group_by(Concentration) %>%
    tidyr::nest()

  linear_models <- data %>%
    dplyr::mutate(
      mod = purrr::map(data, group_model),
      glance = purrr::map(mod, broom::glance),
      tidy = purrr::map(mod, broom::tidy),
      augment = purrr::map(mod, broom::augment),
      rsq = glance %>% purrr::map_dbl("r.squared")
    )
}


# extract initial velocities ----------------------------------------------

#' hitachi_extract_rate
#'
#' This function extracts slopes and their std.errors from linear models per concentration
#' @param data output from hitachi_extract_lm with one column tidy containing broom::tidy parameters
#'
#' @export
#' @examples
#' hitachi_extract_rate()

hitachi_extract_rate <- function (data) {
data <- data %>%
  tidyr::unnest(. , tidy) %>%
  dplyr::filter(. , term == "Time") %>%
  dplyr::select(S = Concentration, v = estimate, std.error)
}


# Build modelmatrix -------------------------------------------------------

#' hitachi_model
#'
#' This function fits models model.nls1 and model.nls2 to initial velocity in dependence of substrate concentration and provides a grid for fitting.
#' @param data tibble with column S (Substrate concentration), v (initial velocity) and std.error (of the slope of the lm)
#'
#' @export
#' @examples
#' hitachi_model()

hitachi_model <- function (data) {
  modelmatrix <- dplyr::tibble(S = seq(0, max(data$S), length.out = 100)) %>%
    modelr::gather_predictions(model.nls1(data), model.nls2(data), .pred = "Velocity_mod")
}
