#' hitachi_tibble
#'
#' takes left_join(hitchi_read / hitachi_read_xls, design, by = Hitachi_ID) as
#' input and does background subtraction and sa-calculation on basis of provided
#' parameters V_RuBP and C_RuBP. If provided, adds additional groupwise information by
#' left_join(., treatment, by = 'Group').
#' @param data left_join(hitchi_read / hitachi_read_xls, design, by = Hitachi_ID) as input. Minimum columns 'Group' and 'CPM'.
#' @param treatment additional file providing groupwise information, e.g substrate concentration. Joined with data by column 'Group'.
#' @param Join_bg default to column name 'Join_bg', provides ties for background counts (e.g. dbl values).
#' @param Join_sa default to column name 'Join_sa', provides ties for specific activity (sa) counts (e.g. dbl values).
#' @param V_RuBP Volume [L] of RuBP-solution added to sa-reactions - default to 0.000004 .
#' @param C_RuBP Concentration [M] of RuBP-solution added to sa-reactions - default to 0.0006 .
#' @param convert default to T - converts CPM in Activity_nmol fixed carbon after groupwise background subtraction (Group == "control"). Requires parameters V_RuBP and C_RuBP!
#' @export
#' @importFrom magrittr %>% %$%
#' @importFrom lubridate ymd_hms
#' @importFrom rlang enquo quo_name !!
#' @examples
#' hitachi_tibble()

hitachi_tibble <- function(data, treatment, Join_bg = Join_bg,
                           Join_sa = Join_sa, V_RuBP = 0.000004,
                           C_RuBP = 0.0006, convert = TRUE) {

  # stop if convert is not logical
  checkmate::assert_logical(convert, len = 1)

  # stop if treatment given, but no join by `Group` possible
  if (! missing(treatment)) {
    checkmate::assert_tibble(treatment)
    checkmate::assert_names(base::names(treatment), must.include = c("Group"))
    checkmate::assert_names(base::names(data), must.include = c("Group"))
  }

  if (convert == TRUE) {

    Join_bg <- rlang::enquo(Join_bg)
    Join_sa <- rlang::enquo(Join_sa)

    checkmate::assert_tibble(data)
    checkmate::assert_names(base::names(data), must.include = c("Group", "CPM"))
    checkmate::assert_names(data$Group %>% unique, must.include = c("sa", "bg"))
    checkmate::check_double(V_RuBP, lower = 0, len = 1)
    checkmate::check_double(C_RuBP, lower = 0, len = 1)

    checkmate::check_subset(c(Join_sa, Join_bg), setdiff(names(data),
                                    c("Set_ID", "Hitachi_ID", "Program",
                                      "Counting_rounds", "Counting_times",
                                      "Counts", "CPM", "Exp_date", "Exp_time")))

    control <- data %>%
      dplyr::filter(Group == "bg") %>%
      dplyr::group_by(!! Join_bg) %>%
      dplyr::summarise(Background_avg = base::mean(CPM), Background_sd = stats::sd(CPM))

    specific_activity <- data %>%
      dplyr::filter(Group == "sa") %>%
      dplyr::group_by(!! Join_sa) %>%
      dplyr::summarise(sa_avg = base::mean(CPM), sa_sd = stats::sd(CPM)) %>%
      dplyr::mutate(SpecificActivity = sa_avg / (C_RuBP * V_RuBP))

    table_activity <- data %>%
      dplyr::filter(Group != "bg" & Group != "sa") %>%
      dplyr::left_join(., control, by = rlang::quo_name(Join_bg)) %>%
      dplyr::left_join(., specific_activity, by = rlang::quo_name(Join_sa)) %>%
      dplyr::mutate(Activity_nmol = ((CPM - Background_avg) / SpecificActivity) * 10^9) %>%
      dplyr::mutate(., CO_2 = (CPM - Background_avg) / SpecificActivity)
  }

  else if (convert == FALSE) {
    table_activity <- data
  }

  # stop if `treatment` is given, but does not contain `Group`
  if (! missing(treatment)) {
    table_activity %>%
        dplyr::left_join(., treatment, by = "Group")
  }
}
