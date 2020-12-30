#' hitachi_read
#'
#' This function allows you to import data from Hitachi LSC-8000 (.txt files) in tidy format.
#' @param file Hitachi output, interpreted by and saved as .txt
#' @export
#' @importFrom magrittr %>% %$%
#' @importFrom lubridate ymd_hms
#' @examples
#' hitachi_read()

hitachi_read <- function(file) {
  readr::read_csv(file = file, skip = 2,
                  col_names = c("Program", "Hitachi_ID", "Counting_rounds", "NA1",
                                "NA2", "Counting_times", "Counts", "CPM",
                                paste0("NA", 3:13), "Exp_date", "Exp_time", "NA14",
                                "NA15"),
                  col_types = cols_only(Program = col_double(),
                                        Hitachi_ID = col_double(),
                                        Counting_rounds = col_double(),
                                        Counting_times = col_double(),
                                        Counts = col_double(),
                                        CPM = col_double(),
                                        Exp_date = col_date(format = ""),
                                        Exp_time = col_time(format = ""))) %>%
    dplyr::mutate(Set_ID = magrittr::extract(., 1,7:8) %>%
                    tidyr::unite(., col = "Set_ID", Exp_date, Exp_time) %>%
                    lubridate::ymd_hms())  %>%
    dplyr::select("Set_ID", "Hitachi_ID", "Program", "Counting_rounds", "Counting_times",
                  "Counts", "CPM", "Exp_date", "Exp_time")
}


