################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Combine results from simulation study
################################################################################

merge_results = function(filename){
  load(filename)
  map_dfr(results, bind_rows)
}

scenarios = list.files(here::here("results", "20260415"), full.names = TRUE)
sim_df = map_dfr(scenarios, merge_results)

sim_df1 <- sim_df %>%
  pivot_longer(cols = matches("group"),
               names_to = "group") %>%
  mutate(true_mean = case_match(group,
                                "group_1" ~ 0,
                                "group_2" ~ 5,
                                "group_3" ~ 10,
                                "group_4" ~ 20))

# calculating bias of group means
bias_df = sim_df1 %>%
  filter(statistic == "mean") %>%
  group_by(method, n, group) %>%
  summarize(nsim = n(),
            bias = mean(value, na.rm = TRUE) - true_mean,
            var_bias = sum((value - mean(value, na.rm = TRUE))^2)/(nsim*(nsim-1)),
            se_bias = sqrt(var_bias)) %>%
  ungroup() %>%
  distinct()

# calculate coverage of beta_hat for all three methods
coverage_df = sim_df1 %>%
  filter(statistic %in% c("2.5%", "97.5%")) %>%
  mutate(statistic = if_else(statistic == "2.5%", "lower", "upper")) %>%
  arrange(method, n, seed, group, statistic) %>%
  mutate(has_coverage = case_when(statistic == "lower" ~ NA_real_,
                                  (true_mean <= value &
                                     true_mean >= lag(value, 1)) ~ 1,
                                   T ~ 0)) %>%
  group_by(method, n, group, statistic) %>%
  summarize(coverage = mean(has_coverage),
            nsim = sum(!is.na(has_coverage))) %>%
  ungroup() %>%
  mutate(var_coverage = coverage * (1-coverage)/n(),
         se_coverage = sqrt(var_coverage)) %>%
  filter(statistic == "upper") %>%
  select(-statistic)
