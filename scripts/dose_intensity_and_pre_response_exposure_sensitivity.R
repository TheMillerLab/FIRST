# ==============================================================================
# DOSE INTENSITY + PRE-RESPONSE EXPOSURE SENSITIVITY ANALYSIS
# ==============================================================================
#
# Purpose:
#   1. Build a regimen-standardized dose-intensity covariate reflecting schedule
#      adherence across the full delivered treatment course.
#   2. Build a response-truncated exposure sensitivity dataset in which only
#      doses administered on or before sys_resp_clin_date are counted.
#
# Study cohort:
#   The study dataframe (df) defines the analytic cohort.
#   The operational infusion table is used only to inform schedule
#   classification when available.
#
# Key outputs:
#   - study_anchor
#   - missing_mrns
#   - unmatched_after_filter
#   - med_episode_check
#   - pembro_med_summary
#   - expected_interval_crosswalk
#   - dose_intensity_patient
#   - dose_count_check
#   - dose_level_timing
#   - dose_level_summary
#   - dose_intensity_pre_response
#   - pre_response_summary
#   - df_model
# ==============================================================================

suppressPackageStartupMessages({
  source(here::here("scripts/load_packages.R"))
})

# ==============================================================================
# 0. FILE PATHS AND DATA LOAD
# ==============================================================================

dir <- "/Users/davidmiller/Partners HealthCare Dropbox/David Miller/mLab/Projects/FIRST and Neoadjuvant Outcomes in CSCC at MGH-MEEI"

dir_2 <- "/Users/davidmiller/Partners HealthCare Dropbox/David Miller/David Partners DropBox/Skin Cancer Program Operational Tool"

df <- open_recent_file(
  directory = file.path(dir, "files/Post Patient Characteristic Processing Data Frame")
)

df_medications_final <- open_recent_file(
  directory = file.path(dir_2, "files/Infusion Data/cleaned")
)

# ==============================================================================
# 1. DEFINE STUDY COHORT (SOURCE OF TRUTH)
# ==============================================================================

study_anchor <- df %>%
  transmute(
    MRN = str_trim(as.character(MRN)),
    patient,
    study_sys_med = str_to_lower(str_trim(sys_med)),
    study_regimen = str_to_lower(str_trim(regimen)),
    sys_start_date = as.Date(sys_start_date),
    dose_1_date = as.Date(dose_1_date),
    sys_resp_clin_date = as.Date(sys_resp_clin_date),
    num_doses
  ) %>%
  distinct()

# ==============================================================================
# 2. CLEAN OPERATIONAL MEDICATION TABLE
# ==============================================================================

med_clean_all <- df_medications_final %>%
  mutate(
    MRN = str_trim(as.character(MRN)),
    Patient = str_trim(as.character(Patient)),
    medication_clean = str_to_lower(str_trim(medication)),
    regimen_clean = str_to_lower(str_trim(regimen)),
    dose_clean = str_to_upper(str_trim(dose)),
    Date = as.Date(Date)
  )

mrn_check <- study_anchor %>%
  distinct(MRN) %>%
  left_join(
    med_clean_all %>%
      distinct(MRN) %>%
      mutate(in_med_table = TRUE),
    by = "MRN"
  ) %>%
  mutate(in_med_table = replace_na(in_med_table, FALSE))

missing_mrns <- mrn_check %>%
  filter(!in_med_table)
missing_mrns
med_clean <- med_clean_all %>%
  semi_join(study_anchor, by = "MRN")

# ==============================================================================
# 3. KEEP RELEVANT IMMUNOTHERAPY ROWS
# ==============================================================================

med_io <- med_clean %>%
  filter(
    medication_clean %in% c(
      "pembrolizumab",
      "cemiplimab",
      "nivolumab",
      "ipilimumab",
      "nivolumab + ipilimumab"
    )
  )

# ==============================================================================
# 4. MATCH OPERATIONAL MEDICATION ROWS TO STUDY DRUG
# ==============================================================================

med_matched <- med_io %>%
  inner_join(study_anchor, by = "MRN") %>%
  filter(
    case_when(
      str_detect(study_sys_med, "pembrolizumab") ~ medication_clean == "pembrolizumab",
      str_detect(study_sys_med, "cemiplimab") ~ medication_clean == "cemiplimab",
      str_detect(study_sys_med, "nivolumab") & str_detect(study_sys_med, "ipilimumab") ~
        medication_clean %in% c("nivolumab", "ipilimumab", "nivolumab + ipilimumab"),
      str_detect(study_sys_med, "nivolumab") ~ medication_clean == "nivolumab",
      str_detect(study_sys_med, "ipilimumab") ~ medication_clean == "ipilimumab",
      TRUE ~ FALSE
    )
  )

unmatched_after_filter <- study_anchor %>%
  distinct(MRN, patient, study_sys_med, sys_start_date, dose_1_date) %>%
  anti_join(
    med_matched %>% distinct(MRN),
    by = "MRN"
  )

# ==============================================================================
# 5. ALIGN OPERATIONAL MEDICATION RECORDS TO STUDY TREATMENT EPISODE
# ==============================================================================

med_episode_check <- med_matched %>%
  group_by(MRN, patient, study_sys_med, sys_start_date, dose_1_date) %>%
  arrange(Date, .by_group = TRUE) %>%
  summarise(
    first_med_date = first(Date),
    days_from_sys_start = as.numeric(first_med_date - sys_start_date),
    days_from_dose1 = as.numeric(first_med_date - dose_1_date),
    episode_match_7d = abs(days_from_sys_start) <= 7,
    episode_match_14d = abs(days_from_sys_start) <= 14,
    .groups = "drop"
  )

med_episode_aligned <- med_matched %>%
  inner_join(
    med_episode_check %>%
      filter(episode_match_7d) %>%
      select(MRN, first_med_date),
    by = "MRN"
  ) %>%
  filter(Date >= first_med_date)

# ==============================================================================
# 6. PEMBROLIZUMAB SCHEDULE CLASSIFICATION USING OPERATIONAL TABLE
# ==============================================================================

pembro_interval_screen_df <- med_episode_aligned %>%
  filter(medication_clean == "pembrolizumab") %>%
  arrange(MRN, Date) %>%
  group_by(MRN, patient) %>%
  mutate(
    interval_days = as.numeric(Date - lag(Date))
  ) %>%
  ungroup()

pembro_interval_screen_df %>%
  count(dose_clean, sort = TRUE)

pembro_med_summary <- pembro_interval_screen_df %>%
  group_by(MRN, patient) %>%
  summarise(
    n_rows = n(),
    any_200mg = any(dose_clean == "200 MG", na.rm = TRUE),
    any_400mg = any(dose_clean == "400 MG", na.rm = TRUE),
    median_interval_days = median(interval_days, na.rm = TRUE),
    n_q3_like = sum(between(interval_days, 14, 28), na.rm = TRUE),
    n_q6_like = sum(between(interval_days, 35, 49), na.rm = TRUE),
    likely_schedule = case_when(
      any_200mg & any_400mg ~ "mixed",
      any_400mg ~ "q6",
      any_200mg ~ "q3",
      n_q6_like >= 2 & n_q6_like > n_q3_like ~ "q6",
      n_q3_like >= 1 & n_q3_like >= n_q6_like ~ "q3",
      TRUE ~ "unknown"
    ),
    .groups = "drop"
  )

other_io_summary <- med_episode_aligned %>%
  filter(medication_clean %in% c("cemiplimab", "nivolumab", "ipilimumab", "nivolumab + ipilimumab")) %>%
  group_by(MRN, patient, medication_clean) %>%
  summarise(
    any_350mg = any(dose_clean == "350 MG", na.rm = TRUE),
    any_240mg = any(dose_clean == "240 MG", na.rm = TRUE),
    any_480mg = any(dose_clean == "480 MG", na.rm = TRUE),
    any_1mgkg = any(dose_clean == "1 MG/KG", na.rm = TRUE),
    any_3mgkg = any(dose_clean == "3 MG/KG", na.rm = TRUE),
    .groups = "drop"
  )

# ==============================================================================
# 6B. PEMBROLIZUMAB INTERVAL-LEVEL EXPECTED TIMING (HYBRID METHOD)
# ==============================================================================

pembro_interval_screen_df_clean <- pembro_interval_screen_df %>%
  distinct(MRN, patient, Date, dose_clean, .keep_all = TRUE)

pembro_interval_screen_df_clean |> group_by(MRN) |> slice_head() |> nrow()

pembro_interval_expected <- pembro_interval_screen_df_clean %>%
  arrange(MRN, Date) %>%
  group_by(MRN, patient) %>%
  mutate(
    next_date = lead(Date),
    next_dose_clean = lead(dose_clean),
    observed_interval_days = as.numeric(next_date - Date),
    
    # observed-interval pattern for fallback classification
    observed_interval_pattern = case_when(
      between(observed_interval_days, 14, 28) ~ "q3_like",
      between(observed_interval_days, 35, 49) ~ "q6_like",
      TRUE ~ "other"
    ),
    
    # hybrid expected interval:
    # prioritize recorded dose; if missing, fall back to observed interval
    expected_interval_days_interval = case_when(
      dose_clean == "200 MG" ~ 21,
      dose_clean == "400 MG" ~ 42,
      is.na(dose_clean) & observed_interval_pattern == "q3_like" ~ 21,
      is.na(dose_clean) & observed_interval_pattern == "q6_like" ~ 42,
      TRUE ~ NA_real_
    ),
    
    expected_interval_source = case_when(
      dose_clean == "200 MG" ~ "recorded_200mg",
      dose_clean == "400 MG" ~ "recorded_400mg",
      is.na(dose_clean) & observed_interval_pattern == "q3_like" ~ "inferred_from_interval_q3",
      is.na(dose_clean) & observed_interval_pattern == "q6_like" ~ "inferred_from_interval_q6",
      TRUE ~ "unclassified"
    )
  ) %>%
  ungroup() %>%
  filter(!is.na(next_date))

pembro_interval_expected |> group_by(MRN) |> slice_head() |> nrow()

pembro_interval_screen_df_clean %>%
  count(MRN, patient) %>%
  filter(n == 1)

pembro_interval_patient_core <- pembro_interval_expected %>%
  group_by(MRN, patient) %>%
  summarise(
    n_intervals_observed = n(),
    n_intervals_with_expected = sum(!is.na(expected_interval_days_interval)),
    n_recorded_200mg = sum(expected_interval_source == "recorded_200mg", na.rm = TRUE),
    n_recorded_400mg = sum(expected_interval_source == "recorded_400mg", na.rm = TRUE),
    n_inferred_q3 = sum(expected_interval_source == "inferred_from_interval_q3", na.rm = TRUE),
    n_inferred_q6 = sum(expected_interval_source == "inferred_from_interval_q6", na.rm = TRUE),
    n_unclassified = sum(expected_interval_source == "unclassified", na.rm = TRUE),
    expected_days_sum = sum(expected_interval_days_interval, na.rm = TRUE),
    observed_days_sum = sum(observed_interval_days, na.rm = TRUE),
    dose_intensity_pembro_interval = case_when(
      observed_days_sum <= 0 ~ NA_real_,
      expected_days_sum <= 0 ~ NA_real_,
      TRUE ~ pmin(1, expected_days_sum / observed_days_sum)
    ),
    dose_intensity_pembro_interval_uncapped = case_when(
      observed_days_sum <= 0 ~ NA_real_,
      expected_days_sum <= 0 ~ NA_real_,
      TRUE ~ expected_days_sum / observed_days_sum
    ),
    dose_intensity_pembro_interval_uncapped_wins = case_when(
      is.na(dose_intensity_pembro_interval_uncapped) ~ NA_real_,
      TRUE ~ pmin(dose_intensity_pembro_interval_uncapped, 1.5)
    ),
    dose_intensity_pembro_interval_uncapped_centered =
      dose_intensity_pembro_interval_uncapped -
      mean(dose_intensity_pembro_interval_uncapped, na.rm = TRUE
           ),
    .groups = "drop"
  )

pembro_interval_patient <- pembro_interval_screen_df_clean %>%
  distinct(MRN, patient) %>%
  left_join(
    pembro_interval_patient_core,
    by = c("MRN", "patient")
  ) %>%
  mutate(
    n_intervals_observed = replace_na(n_intervals_observed, 0),
    n_intervals_with_expected = replace_na(n_intervals_with_expected, 0),
    n_recorded_200mg = replace_na(n_recorded_200mg, 0),
    n_recorded_400mg = replace_na(n_recorded_400mg, 0),
    n_inferred_q3 = replace_na(n_inferred_q3, 0),
    n_inferred_q6 = replace_na(n_inferred_q6, 0),
    n_unclassified = replace_na(n_unclassified, 0),
    expected_days_sum = replace_na(expected_days_sum, 0),
    observed_days_sum = replace_na(observed_days_sum, 0),
    dose_intensity_pembro_interval = case_when(
      n_intervals_observed == 0 ~ 1,
      TRUE ~ dose_intensity_pembro_interval
    )
  )

pembro_interval_patient |> distinct(MRN, patient) |> nrow()

pembro_interval_expected %>%
  count(expected_interval_source, sort = TRUE)

# ==============================================================================
# 7. ASSIGN EXPECTED INTERVAL DAYS
# ==============================================================================

expected_interval_crosswalk <- study_anchor %>%
  select(MRN, patient, study_sys_med) %>%
  left_join(
    pembro_med_summary %>%
      select(MRN, likely_schedule, any_200mg, any_400mg),
    by = "MRN"
  ) %>%
  mutate(
    likely_schedule = replace_na(likely_schedule, "unknown"),
    any_200mg = replace_na(any_200mg, FALSE),
    any_400mg = replace_na(any_400mg, FALSE),
    expected_interval_days = case_when(
      str_detect(study_sys_med, "cemiplimab") ~ 21,
      str_detect(study_sys_med, "pembrolizumab") & likely_schedule == "q6" ~ 42,
      str_detect(study_sys_med, "pembrolizumab") & likely_schedule == "q3" ~ 21,
      str_detect(study_sys_med, "pembrolizumab") & likely_schedule == "mixed" ~ NA_real_,
      str_detect(study_sys_med, "pembrolizumab") & likely_schedule == "unknown" ~ 21,
      str_detect(study_sys_med, "nivolumab") ~ 14,
      TRUE ~ NA_real_
    ),
    interval_source = case_when(
      str_detect(study_sys_med, "pembrolizumab") & likely_schedule %in% c("q3", "q6") ~ "medication_table_patient_level",
      str_detect(study_sys_med, "pembrolizumab") & likely_schedule == "mixed" ~ "medication_table_interval_level_required",
      str_detect(study_sys_med, "pembrolizumab") & likely_schedule == "unknown" ~ "default_pembro_q3",
      str_detect(study_sys_med, "cemiplimab") ~ "default_cemiplimab_q3",
      str_detect(study_sys_med, "nivolumab") ~ "default_nivolumab_q2",
      TRUE ~ "unclassified"
    )
  )

# ==============================================================================
# 8. BUILD ANALYSIS DATAFRAME WITH EXPECTED INTERVAL
# ==============================================================================

df_dose_intensity <- df %>%
  mutate(MRN = str_trim(as.character(MRN))) %>%
  left_join(
    expected_interval_crosswalk %>%
      select(
        MRN, expected_interval_days, interval_source,
        likely_schedule, any_200mg, any_400mg
      ),
    by = "MRN"
  ) %>%
  mutate(
    has_med_match = MRN %in% med_matched$MRN,
    med_episode_match_7d = MRN %in% (med_episode_check %>% filter(episode_match_7d) %>% pull(MRN))
  )

# ==============================================================================
# 9. PIVOT STUDY DOSE DATES LONG (FULL COURSE)
# ==============================================================================

dose_long <- df_dose_intensity %>%
  select(
    MRN, patient, sys_med, regimen, num_doses,
    expected_interval_days, interval_source, likely_schedule,
    sys_resp_clin, sys_resp_clin_date, first_resp_date, death_date,
    starts_with("dose_") & ends_with("_date")
  ) %>%
  pivot_longer(
    cols = starts_with("dose_") & ends_with("_date"),
    names_to = "dose_number",
    values_to = "dose_date"
  ) %>%
  mutate(
    dose_number = str_extract(dose_number, "\\d+") %>% as.integer(),
    dose_date = as.Date(dose_date),
    sys_resp_clin_date = as.Date(sys_resp_clin_date),
    first_resp_date = as.Date(first_resp_date),
    death_date = as.Date(death_date),
    response_cutoff_date = case_when(
      !is.na(sys_resp_clin_date) ~ sys_resp_clin_date,
      is.na(sys_resp_clin_date) & !is.na(first_resp_date) ~ first_resp_date,
      is.na(sys_resp_clin_date) & is.na(first_resp_date) &
        str_detect(str_to_lower(sys_resp_clin), "progress") & !is.na(death_date) ~ death_date,
      TRUE ~ as.Date(NA)
    ),
    dose_before_response = dose_date <= response_cutoff_date
  ) %>%
  filter(!is.na(dose_date))
# ==============================================================================
# 10. PATIENT-LEVEL DOSE INTENSITY (FULL COURSE)
# ==============================================================================

dose_intensity_patient <- dose_long %>%
  group_by(MRN, patient) %>%
  arrange(dose_number, .by_group = TRUE) %>%
  summarise(
    sys_med = first(sys_med),
    regimen = first(regimen),
    num_doses_recorded = n(),
    num_doses_df = first(num_doses),
    expected_interval_days = first(expected_interval_days),
    interval_source = first(interval_source),
    likely_schedule = first(likely_schedule),
    dose_1_date = first(dose_date),
    final_dose_date = last(dose_date),
    observed_elapsed_days = as.numeric(final_dose_date - dose_1_date),
    expected_elapsed_days = (num_doses_recorded - 1) * expected_interval_days,
    dose_intensity_final = case_when(
      num_doses_recorded <= 1 ~ 1,
      is.na(expected_interval_days) ~ NA_real_,
      observed_elapsed_days <= 0 ~ NA_real_,
      TRUE ~ pmin(1, expected_elapsed_days / observed_elapsed_days)
    ),
    dose_intensity_final_uncapped = case_when(
      num_doses_recorded <= 1 ~ 1,
      is.na(expected_interval_days) ~ NA_real_,
      observed_elapsed_days <= 0 ~ NA_real_,
      TRUE ~ expected_elapsed_days / observed_elapsed_days
    ),
    dose_intensity_final_uncapped_wins = case_when(
      is.na(dose_intensity_final_uncapped) ~ NA_real_,
      TRUE ~ pmin(dose_intensity_final_uncapped, 1.5)
    ),
    dose_intensity_final_uncapped_centered =
      dose_intensity_final_uncapped -
      mean(dose_intensity_final_uncapped, na.rm = TRUE
           ),
    mean_interval_days = if_else(
      num_doses_recorded >= 2,
      observed_elapsed_days / (num_doses_recorded - 1),
      NA_real_
    ),
    interval_deviation_days = case_when(
      num_doses_recorded <= 1 ~ 0,
      is.na(expected_interval_days) ~ NA_real_,
      TRUE ~ mean_interval_days - expected_interval_days
    ),
    .groups = "drop"
  )

dose_count_check <- dose_intensity_patient %>%
  mutate(dose_count_match = num_doses_recorded == num_doses_df)

# ==============================================================================
# 11. DOSE-LEVEL TIMING TABLE (FULL COURSE)
# ==============================================================================

dose_level_timing <- dose_long %>%
  group_by(MRN, patient) %>%
  arrange(dose_number, .by_group = TRUE) %>%
  mutate(
    first_dose_date = first(dose_date),
    observed_day_from_dose1 = as.numeric(dose_date - first_dose_date)
  ) %>%
  ungroup() %>%
  mutate(
    expected_day_from_dose1 = (dose_number - 1) * expected_interval_days,
    deviation_days = observed_day_from_dose1 - expected_day_from_dose1,
    within_3_days = abs(deviation_days) <= 3,
    within_7_days = abs(deviation_days) <= 7,
    within_10_days = abs(deviation_days) <= 10
  )

dose_level_summary <- dose_level_timing %>%
  filter(dose_number >= 2) %>%
  group_by(dose_number) %>%
  summarise(
    n_patients = n(),
    median_deviation_days = median(deviation_days, na.rm = TRUE),
    iqr_deviation_days = IQR(deviation_days, na.rm = TRUE),
    prop_within_3_days = mean(within_3_days, na.rm = TRUE),
    prop_within_7_days = mean(within_7_days, na.rm = TRUE),
    prop_within_10_days = mean(within_10_days, na.rm = TRUE),
    .groups = "drop"
  )

# ==============================================================================
# 12. PRE-RESPONSE EXPOSURE SENSITIVITY
# ==============================================================================

dose_intensity_pre_response <- dose_long %>%
  filter(dose_before_response) %>%
  group_by(MRN, patient) %>%
  arrange(dose_number, .by_group = TRUE) %>%
  summarise(
    sys_med = first(sys_med),
    regimen = first(regimen),
    num_doses_pre_response = n(),
    expected_interval_days = first(expected_interval_days),
    first_dose_pre_response = first(dose_date),
    last_dose_pre_response = last(dose_date),
    observed_elapsed_days_pre_response =
      as.numeric(last_dose_pre_response - first_dose_pre_response),
    expected_elapsed_days_pre_response =
      (num_doses_pre_response - 1) * expected_interval_days,
    dose_intensity_pre_response = case_when(
      num_doses_pre_response <= 1 ~ 1,
      is.na(expected_interval_days) ~ NA_real_,
      observed_elapsed_days_pre_response <= 0 ~ NA_real_,
      TRUE ~ pmin(1, expected_elapsed_days_pre_response /
                    observed_elapsed_days_pre_response)
    ),
    dose_intensity_pre_response_uncapped = case_when(
      num_doses_pre_response <= 1 ~ 1,
      is.na(expected_interval_days) ~ NA_real_,
      observed_elapsed_days_pre_response <= 0 ~ NA_real_,
      TRUE ~ expected_elapsed_days_pre_response /
        observed_elapsed_days_pre_response
    ),
    dose_intensity_pre_response_uncapped_wins = case_when(
      is.na(dose_intensity_pre_response_uncapped) ~ NA_real_,
      TRUE ~ pmin(dose_intensity_pre_response_uncapped, 1.5)
    ),
    dose_intensity_pre_response_uncapped_centered =
      dose_intensity_pre_response_uncapped -
      mean(dose_intensity_pre_response_uncapped, na.rm = TRUE
           ),
    .groups = "drop"
  )

# ==============================================================================
# 13. MODEL-READY DATAFRAME
# ==============================================================================

df_model <- df_dose_intensity %>%
  mutate(response_cutoff_date = case_when(
    !is.na(sys_resp_clin_date) ~ sys_resp_clin_date,
    is.na(sys_resp_clin_date) & !is.na(first_resp_date) ~ first_resp_date,
    is.na(sys_resp_clin_date) & is.na(first_resp_date) &
      str_detect(str_to_lower(sys_resp_clin), "progress") & !is.na(death_date) ~ death_date,
    TRUE ~ as.Date(NA)
  )) %>%
  left_join(
    dose_intensity_patient %>%
      select(
        MRN, patient,
        num_doses_recorded,
        dose_intensity_final,
        mean_interval_days,
        interval_deviation_days
      ),
    by = c("MRN", "patient")
  ) %>%
  left_join(
    dose_intensity_pre_response %>%
      select(
        MRN, patient,
        num_doses_pre_response,
        dose_intensity_pre_response,
        dose_intensity_pre_response_uncapped,
        dose_intensity_pre_response_uncapped_wins,
        dose_intensity_pre_response_uncapped_centered
      ),
    by = c("MRN", "patient")
  ) %>%
  left_join(
    pembro_interval_patient,
    by = c("MRN", "patient")
  ) %>%
  mutate(
    dose_intensity_final = case_when(
      str_detect(str_to_lower(sys_med), "pembrolizumab") &
        likely_schedule == "mixed" ~ dose_intensity_pembro_interval,
      TRUE ~ dose_intensity_final
    ),
    dose_intensity_pre_response = case_when(
      str_detect(str_to_lower(sys_med), "pembrolizumab") &
        likely_schedule == "mixed" ~ dose_intensity_pembro_interval,
      TRUE ~ dose_intensity_pre_response
    ),
    dose_intensity_pre_response_uncapped = case_when(
      str_detect(str_to_lower(sys_med), "pembrolizumab") &
        likely_schedule == "mixed" ~ dose_intensity_pembro_interval_uncapped,
      TRUE ~ dose_intensity_pre_response_uncapped
    ),
    dose_intensity_pre_response_uncapped_wins = case_when(
      str_detect(str_to_lower(sys_med), "pembrolizumab") &
        likely_schedule == "mixed" ~ dose_intensity_pembro_interval_uncapped_wins,
      TRUE ~ dose_intensity_pre_response_uncapped_wins
    ),
    dose_intensity_pre_response_uncapped_centered =
      dose_intensity_pre_response_uncapped -
      mean(dose_intensity_pre_response_uncapped, na.rm = TRUE),
    extra_doses_after_response = num_doses_recorded - num_doses_pre_response,
    dose_intensity_source = case_when(
      str_detect(str_to_lower(sys_med), "pembrolizumab") &
        likely_schedule == "mixed" ~ "pembro_interval_level",
      TRUE ~ "patient_level_expected_interval"
    ),
    agent_schedule = case_when(
      str_detect(sys_med, "Pembrolizumab") & likely_schedule == "q3" ~ "pembro_q3",
      str_detect(sys_med, "Pembrolizumab") & likely_schedule == "q6" ~ "pembro_q6",
      str_detect(sys_med, "Pembrolizumab") & likely_schedule == "mixed" ~ "pembro_mixed",
      str_detect(sys_med, "Cemiplimab") ~ "cemiplimab",
      str_detect(sys_med, "Ipilimumab \\+ Nivolumab") ~ "Ipi-nivo",
      TRUE ~ "other"
    )
  )
# For analyses evaluating post-response exposure, treatment doses were truncated at the date of clinical response adjudication (sys_resp_clin_date). When a final adjudication date was unavailable, the initial response assessment date was used as the truncation point.
# ==============================================================================
# 14. PRE-RESPONSE SUMMARY OBJECTS
# ==============================================================================

pre_response_summary <- df_model %>%
  summarise(
    n_total = n(),
    n_any_extra_after_response = sum(num_doses_recorded > num_doses_pre_response, na.rm = TRUE),
    prop_any_extra_after_response = mean(num_doses_recorded > num_doses_pre_response, na.rm = TRUE)
  )

extra_doses_distribution <- df_model %>%
  count(extra_doses_after_response, sort = TRUE)
extra_doses_distribution

extra_doses_summary_positive <- df_model %>%
  filter(extra_doses_after_response > 0) %>%
  summarise(
    median_extra = median(extra_doses_after_response),
    iqr_extra = IQR(extra_doses_after_response),
    max_extra = max(extra_doses_after_response)
  )

# ==============================================================================
# 15. OPTIONAL QC OBJECTS
# ==============================================================================

cohort_qc_summary <- tibble(
  n_study_mrns = n_distinct(study_anchor$MRN),
  n_missing_from_med_table = nrow(missing_mrns),
  n_unmatched_after_drug_filter = nrow(unmatched_after_filter),
  n_episode_match_7d = med_episode_check %>% filter(episode_match_7d) %>% nrow(),
  n_pembro_classified_q3 = pembro_med_summary %>% filter(likely_schedule == "q3") %>% nrow(),
  n_pembro_classified_q6 = pembro_med_summary %>% filter(likely_schedule == "q6") %>% nrow(),
  n_pembro_classified_mixed = pembro_med_summary %>% filter(likely_schedule == "mixed") %>% nrow(),
  n_pembro_classified_unknown = pembro_med_summary %>% filter(likely_schedule == "unknown") %>% nrow()
)
cohort_qc_summary

med_episode_check %>%
  filter(episode_match_7d) %>%
  distinct(MRN) %>%
  summarise(n_episode_matched_mrns = n())

df_model %>%
  summarise(
    n_with_response_cutoff = sum(!is.na(response_cutoff_date)),
    n_missing_response_cutoff = sum(is.na(response_cutoff_date))
  )

# ==============================================================================
# 16B. ONE-PAGE PIPELINE QC SUMMARY
# ==============================================================================

pipeline_qc_summary <- tibble(
  metric = c(
    "Study cohort MRNs",
    "Study MRNs missing from medication table",
    "Study MRNs present in medication table",
    "Study MRNs unmatched after drug filter",
    "Study MRNs with episode match within 7 days",
    "Infusion rows matched within 7 days",
    "Pembrolizumab patients expected clinically",
    "Pembrolizumab patients in medication-based pipeline",
    "Pembrolizumab patients classified q3",
    "Pembrolizumab patients classified q6",
    "Pembrolizumab patients classified mixed",
    "Pembrolizumab patients classified unknown",
    "Patients with response cutoff date",
    "Patients missing response cutoff date",
    "Patients with any extra doses after response",
    "Patients with no extra doses after response"
  ),
  value = c(
    n_distinct(study_anchor$MRN),
    nrow(missing_mrns),
    n_distinct(med_clean$MRN),
    nrow(unmatched_after_filter),
    med_episode_check %>% filter(episode_match_7d) %>% distinct(MRN) %>% nrow(),
    med_episode_check %>% filter(episode_match_7d) %>% nrow(),
    59,  # replace if you want this derived from df instead of hard-coded
    nrow(pembro_med_summary),
    pembro_med_summary %>% filter(likely_schedule == "q3") %>% nrow(),
    pembro_med_summary %>% filter(likely_schedule == "q6") %>% nrow(),
    pembro_med_summary %>% filter(likely_schedule == "mixed") %>% nrow(),
    pembro_med_summary %>% filter(likely_schedule == "unknown") %>% nrow(),
    df_model %>% summarise(n = sum(!is.na(response_cutoff_date))) %>% pull(n),
    df_model %>% summarise(n = sum(is.na(response_cutoff_date))) %>% pull(n),
    df_model %>% summarise(n = sum(extra_doses_after_response > 0, na.rm = TRUE)) %>% pull(n),
    df_model %>% summarise(n = sum(extra_doses_after_response == 0, na.rm = TRUE)) %>% pull(n)
  )
)

pipeline_qc_summary

# ==============================================================================
# 16. OPTIONAL FIGURE
# ==============================================================================

plot_extra_doses_after_response <- df_model %>%
  ggplot(aes(extra_doses_after_response)) +
  geom_histogram(binwidth = 1) +
  labs(
    x = "Additional doses after response adjudication",
    y = "Number of patients",
    title = "Additional Immunotherapy Doses After Response Adjudication"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
plot_extra_doses_after_response
# ==============================================================================
# 17. SAVE KEY OUTPUTS
# ==============================================================================

save_files(
  save_object = df_model,
  filename = "Dose Intensity Sensitivity Model Data",
  subD = "files/Dose Intensity Sensitivity"
)

save_files(
  save_object = dose_level_summary,
  filename = "Dose Level Timing Summary",
  subD = "files/Dose Intensity Sensitivity"
)

save_files(
  save_object = pembro_med_summary,
  filename = "Pembrolizumab Schedule Classification Summary",
  subD = "files/Dose Intensity Sensitivity"
)

save_files(
  save_object = cohort_qc_summary,
  filename = "Dose Intensity Cohort QC Summary",
  subD = "files/Dose Intensity Sensitivity"
)

save_files(
  save_object = pre_response_summary,
  filename = "Pre Response Exposure Summary",
  subD = "files/Dose Intensity Sensitivity"
)

save_files(
  save_object = extra_doses_distribution,
  filename = "Extra Doses After Response Distribution",
  subD = "files/Dose Intensity Sensitivity"
)

save_files(
  save_object = extra_doses_summary_positive,
  filename = "Extra Doses After Response Positive Summary",
  subD = "files/Dose Intensity Sensitivity"
)

save_files(
  save_object = plot_extra_doses_after_response,
  filename = "Extra Doses After Response Histogram",
  subD = "files/Dose Intensity Sensitivity"
)

# ==============================================================================
# 18. OPTIONAL INSPECTION
# ==============================================================================

cohort_qc_summary
pre_response_summary
extra_doses_distribution
extra_doses_summary_positive
dose_count_check %>% count(dose_count_match)
dose_level_summary

df_model %>%
  select(
    MRN, patient, sys_med, num_doses,
    num_doses_recorded, num_doses_pre_response,
    expected_interval_days, interval_source, likely_schedule,
    any_200mg, any_400mg,
    dose_intensity_final, dose_intensity_pre_response,
    extra_doses_after_response,
    has_med_match, med_episode_match_7d
  ) %>%
  arrange(num_doses, patient)

pembro_med_summary %>%
  count(likely_schedule, sort = TRUE)

response_qc_df <- df_model %>%
  select(
    MRN,
    patient,
    sys_med,
    regimen,
    num_doses_recorded,
    num_doses_pre_response,
    extra_doses_after_response,
    first_resp_date,
    sys_resp_clin_date,
    response_cutoff_date
  ) %>%
  arrange(desc(extra_doses_after_response))

dose_intensity_qc_df <- df_model %>%
  select(
    MRN,
    patient,
    sys_med,
    regimen,
    sys_resp_clin,
    event_after_mgmt,
    num_doses_recorded,
    expected_interval_days,
    mean_interval_days,
    interval_deviation_days,
    dose_intensity_final
  ) %>%
  arrange(desc(interval_deviation_days))

pembro_schedule_qc_df <- df_model %>%
  filter(str_detect(str_to_lower(sys_med), "pembro")) %>%
  select(
    MRN,
    patient,
    sys_resp_clin,
    event_after_mgmt,
    likely_schedule,
    any_200mg,
    any_400mg,
    n_recorded_200mg,
    n_recorded_400mg,
    n_inferred_q3,
    n_inferred_q6,
    dose_intensity_pembro_interval
  ) %>%
  arrange(desc(n_recorded_400mg))

clinical_spotcheck_df <- df_model %>%
  select(
    MRN,
    patient,
    sys_med,
    likely_schedule,
    sys_resp_clin,
    event_after_mgmt,
    num_doses_recorded,
    num_doses_pre_response,
    extra_doses_after_response,
    dose_intensity_pre_response,
    dose_intensity_final
  ) %>%
  arrange(desc(num_doses_recorded))


# ==============================================================================
# EXPOSURE VARIABLES FOR DOWNSTREAM MODELING
# ==============================================================================

# The dose-intensity pipeline generates several related exposure variables.
# These capture treatment adherence and exposure before clinical response.

# 1. dose_intensity_final
# --------------------------------------------------
# Regimen-standardized dose intensity calculated across the full treatment
# course using expected vs observed interval timing.
#
# Interpretation:
#   1.0 = perfectly on-schedule treatment
#   <1.0 = treatment delays or interruptions
#
# Note:
#   For mixed pembrolizumab schedules (200 mg Q3 + 400 mg Q6),
#   interval-level expected timing is used rather than a single
#   patient-level expected interval.

# 2. dose_intensity_pre_response
# --------------------------------------------------
# Dose intensity calculated only using doses administered on or before
# the response cutoff date (response_cutoff_date).
#
# This variable is preferred for causal analyses because it avoids
# post-response exposure influencing the exposure metric.

# 3. num_doses_pre_response
# --------------------------------------------------
# Number of doses delivered prior to the response cutoff date.
#
# This provides a simpler exposure measure that is robust to schedule
# differences and may be used in sensitivity analyses.

# Recommended modeling strategy:
# --------------------------------------------------
# Primary exposure:
#   dose_intensity_pre_response
#
# Sensitivity analyses:
#   dose_intensity_final
#   num_doses_pre_response

# 4. agent_schedule
# --------------------------------------------------
# Treatment agent including schedule-specific pembrolizumab categories.
#
# Levels:
#   pembro_q3
#   pembro_q6
#   pembro_mixed
#   cemiplimab
#   nivolumab
#   ipilimumab
#
# This allows modeling of schedule differences that may influence
# dose intensity and response.

# ==============================================================================