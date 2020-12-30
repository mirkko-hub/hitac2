library(tidyverse)
library(modelr)
library(broom)
library(magrittr)
library(jasco2)
library(checkmate)


# import data -------------------------------------------------------------
ha_data <- hitachi_read("/media/mirkko/ToshibaWorkCrypt/Rfolder/NewPackages/hitachi_data/20180109_raw.csv")
design <- read_csv("/media/mirkko/ToshibaWorkCrypt/Rfolder/NewPackages/hitachi_data/20180109_design.csv")
treatment <- read_csv("/media/mirkko/ToshibaWorkCrypt/Rfolder/NewPackages/hitachi_data/20180109_treatment.csv")

# process -----------------------------------------------------------------
data <- left_join(ha_data, design, by = "Hitachi_ID")
table_activity <- hitachi_tibble(data = data, treatment = treatment, V_RuBP = 8e-06, C_RuBP = 5e-04, convert = TRUE)


# plot --------------------------------------------------------------------
# controls
ggplot(data = data %>% filter(Group %in% c("sa", "bg")),
       mapping = aes(x = Group, y = CPM)) +
  geom_boxplot() +
  geom_point() +
  facet_wrap("Group", scales = "free")

ggplot(data = table_activity, mapping = aes(x = Time_s, y = Activity_nmol, colour = factor(Concentration_RuBP_mM))) +
  geom_point()


# experimental plotting and modelling -------------------------------------



# ordered -----------------------------------------------------------------


lm_wrapper <- function(df, response, predictor, ...){
  response <- ensym(response)
  predictor <- ensym(predictor)
  lm(expr(!!response ~ !!predictor), df)
}

hitachi_lm <- function(data, response, predictor, min, max, ...){
  response <- enquo(response)
  predictor <- enquo(predictor)

  data <- data %>%
    filter(., !!predictor < max & !!predictor > min) %>%
    select(-c(Set_ID, Hitachi_ID, Program, Counting_rounds, Counting_times,
              Counts, Exp_date, Exp_time, Join_bg, Join_sa, Notes, FLAG)) %>%
    nest(c(Exp_ID, !!predictor, !!response, Activity_nmol, CO_2, CPM, Time_s))

  data %>%
    mutate(
      mod = map(data, lm_wrapper, predictor = !!predictor, response = !!response),
      glance = map(mod, glance),
      tidy = map(mod, tidy),
      augment = map(mod, augment),
      rsq = glance %>% map_dbl("r.squared")
    )
}

hitachi_lm(data = table_activity, response = CO_2, predictor = Time_s, min = 20, max = 120) %>% View()


# improve with regard to Time_s and proper nesting ------------------------

# implement Time_s requirement for design? / hitachi
# implement Timecourse requirement for hitachi_lm

design <- design %>% mutate(., Timecourse = c(rep("-", 12), rep(1:30, each = 6)))
data <- left_join(ha_data, design, by = "Hitachi_ID")
table_activity <- hitachi_tibble(data = data, treatment = treatment, V_RuBP = 8e-06, C_RuBP = 5e-04, convert = TRUE)

df <- hitachi_lm(data = table_activity, response = CO_2, predictor = Time_s, min = 20, max = 120)
# df %>% ggplot(., aes(y = rsq, x = log(Concentration_RuBP_mM))) + geom_point()

table_activity %>% ggplot(., aes(x = Time_s, y = Activity_nmol, colour = Concentration_RuBP_mM)) + geom_point()

# extract params
hitachi_get_turnover <- function(.rate, .protein_pmol, .unit) {
  unit <- enquo(.unit)

  time <- tibble(unit = c("h", "min", "sec", "msec"),
                 calc = c(1/60, 1, 60, 3600 ))
  calc <- time %>% filter(unit == !!unit) %$% calc
  (((.rate / (.protein_pmol * 1e-12)) * calc))
}


hitachi_extract_params <- function(.data, .bgcorr, .protein_pmol, .unit, ...){

  assert(
    check_names(names(.data), must.include = "tidy"),
    check_names(.data %$% tidy %>% extract2(1) %>% names(), must.include = "estimate"),
    check_logical(.bgcorr, len = 1, any.missing = FALSE),
    check_double(.protein_pmol, lower = 0.001, upper = 1000, any.missing = FALSE, len = 1),
    check_choice(.unit, c("h", "min", "sec", "msec")),
    check_character(.unit, any.missing = FALSE, len = 1),
    combine = "and"
  )

  group_vars <- enquos(..., .named = TRUE)

  df <- .data %>%
    unnest(., tidy) %>%
    filter(., term != "(Intercept)") %>%
    mutate(rate = estimate)

  if (.bgcorr == TRUE) {
    blank_id <- df[str_which(df$Group, pattern = regex("blank", ignore_case = T)), ] %$% Exp_ID
    blank_rate <- df %>% filter(., Exp_ID %in% blank_id) %>% summarise(., rate = mean(rate)) %$% rate
    df <- df %>% filter(! Exp_ID %in% blank_id) %>% mutate(., rate = rate - blank_rate)
  }

  df <- df %>% mutate(turnover = rate %>% map_dbl(., hitachi_get_turnover, .protein_pmol, .unit))

  # Modify the names of the list of quoted dots
  names(group_vars) <- paste0("groups_", names(group_vars))

  df_summary <- df %>%
    group_by(!!!group_vars) %>%  # Unquote-splice as usual
    summarise(.,
              rate_avg = mean(rate),
              rate_sd = sd(rate),
              turnover_avg = mean(turnover),
              turnover_sd = sd(turnover))

  list(data = df,
       summary = df_summary)

}

df_ex <- df %>% hitachi_extract_params(., .bgcorr = FALSE, .protein_pmol = 30, .unit = "sec", Concentration_RuBP_mM)
df_ex$data %>% ggplot(., aes(x = Concentration_RuBP_mM, y = turnover)) + geom_point()

df_ex$summary %>% ggplot() + geom_point(aes(x = groups_Concentration_RuBP_mM, y = turnover_avg)) +
  geom_errorbar(aes(x = groups_Concentration_RuBP_mM, ymin = turnover_avg - turnover_sd, ymax = turnover_avg + turnover_sd)) +
  geom_line(data = mm1, aes(x = S, y = v), colour = "red", linetype = 3, alpha = 0.5)

# models

model_in <- df_ex$summary %>% select(S = groups_Concentration_RuBP_mM, v = turnover_avg, turnover_sd)

model.nls1 <- nls(v ~ Vm * S/(K+S), data = model_in,
                  start = list(K = max(model_in$S)*0.1, Vm = max(model_in$v)))
summary(model.nls1)

mm1 <- tibble(S = seq(0, max(model_in$S), length.out = 100))
mm1$v <- predict(model.nls1, newdata = mm1)











MM_data <- filter(models_parameters, term == "Time" ) %>%
  select(S = Concentration_RuBP, v = estimate, std.error)

model.nls1 <- nls(v ~ Vm * S/(K+S), data = MM_data,
                  start = list(K = max(MM_data$v)/2, Vm = max(MM_data$v)))
summary(model.nls1)
mm1 <- data.frame(S = seq(0, max(MM_data$S), length.out = 100))
mm1$v <- predict(model.nls1, newdata = mm1)

model.nls2 <- nls(v ~ Vm * S/(K+S), data = filter(MM_data, S != max(S) ),
                  start = list(K = max(MM_data$v)/2, Vm = max(MM_data$v)))
summary(model.nls2)
mm2 <- data.frame(S = seq(0, max(MM_data$S), length.out = 100))
mm2$v <- predict(model.nls2, newdata = mm2)



ggplot(MM_data, aes(x = S, y = v)) +
  theme_bw() +
  labs(y = expression(paste("initial velocity ( ", nmol / sec, " )")), x = expression(paste("[RuBP] ( ", mM, " )"))) +
  geom_point(size = 0.5) +
  geom_errorbar(data = MM_data, aes(x = S, ymin = v - std.error, ymax = v + std.error)) +
  geom_line(data = mm1, aes(x = S, y = v), colour = "red", linetype = 3, alpha = 0.2) +
  geom_line(data = mm2, aes(x = S, y = v), colour = "blue", linetype = 3, alpha = 0.2) +
  Theme_1
# ggsave(filename = "/home/mirkkoadmin/Desktop/Rfolder/Project4/Output/20180109_Vmax_temporary_3.tiff", width = 11, height = 7, units = "cm")








at_mod <- lm(Activity_nmol ~ Time_s * Group, data = table_activity)

table_activity <- table_activity %>%
  group_by(Group) %>%
  add_predictions(at_mod) %>%
  add_residuals(at_mod)



