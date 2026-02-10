suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(purrr)
})

# Use the same logical pattern object you already have
neoadj_pattern <- regex(
  "\\b(neoadjuvant|neo-?adjuvant|pre\\s*operative|pre-?op(?:erative)?)\\b",
  ignore_case = TRUE
)

# For highlighting in base gsub (supports \\0 backref)
neoadj_pattern_str <- "(?i)\\b(neoadjuvant|neo-?adjuvant|pre\\s*operative|pre-?op(?:erative)?)\\b"

# ---- Helpers ----

# Build one or more snippets that each include ±pre/post chars around the match.
# Overlapping windows are merged. Terms inside snippets can be highlighted.
extract_neoadj_snippets <- function(text,
                                    pre = 200,
                                    post = 200,
                                    highlight = TRUE,
                                    mark_start = "«",
                                    mark_end   = "»") {
  if (is.null(text) || is.na(text) || !nzchar(text)) return(NA_character_)
  locs <- stringr::str_locate_all(text, neoadj_pattern)[[1]]
  if (is.null(dim(locs)) || nrow(locs) == 0) return(NA_character_)
  
  n <- nchar(text)
  starts <- pmax(locs[, 1] - pre, 1L)
  ends   <- pmin(locs[, 2] + post, n)
  
  # Merge overlapping intervals
  o <- order(starts)
  starts <- starts[o]; ends <- ends[o]
  merged_s <- integer(0); merged_e <- integer(0)
  cur_s <- starts[1]; cur_e <- ends[1]
  if (length(starts) > 1) {
    for (i in 2:length(starts)) {
      if (starts[i] <= cur_e + 1L) {
        cur_e <- max(cur_e, ends[i])
      } else {
        merged_s <- c(merged_s, cur_s); merged_e <- c(merged_e, cur_e)
        cur_s <- starts[i]; cur_e <- ends[i]
      }
    }
  }
  merged_s <- c(merged_s, cur_s); merged_e <- c(merged_e, cur_e)
  
  snippets <- map2_chr(merged_s, merged_e, function(s, e) {
    frag <- stringr::str_sub(text, s, e)
    frag <- stringr::str_replace_all(frag, "[\\r\\n]+", " ")
    frag <- stringr::str_squish(frag)
    if (highlight) {
      # Wrap each matched term inside the snippet
      frag <- stringr::str_replace_all(
        frag, neoadj_pattern,
        paste0(mark_start, "\\1", mark_end)  # group 1 = the matched term
      )
    }
    paste0(if (s > 1L) "…" else "", frag, if (e < n) "…" else "")
  })
  
  paste(snippets, collapse = " || ")
}

extract_neoadj_snippets_safely <- purrr::possibly(extract_neoadj_snippets, otherwise = NA_character_)

# Exact matched terms as they appear in the note (case-preserving, unique)
extract_neoadj_terms_exact <- function(text) {
  if (is.null(text) || is.na(text) || !nzchar(text)) return(NA_character_)
  m <- gregexpr(neoadj_pattern_str, text, perl = TRUE)
  hits <- regmatches(text, m)[[1]]
  if (length(hits) == 0) return(NA_character_)
  paste(unique(hits), collapse = "; ")
}
extract_neoadj_terms_exact_safely <- purrr::possibly(extract_neoadj_terms_exact, otherwise = NA_character_)




           # =============================================================================
# Neoadjuvant module — functions
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(rlang); library(yaml)
})

# ---- tiny utils --------------------------------------------------------------
`%||%` <- function(x, y) if (is.null(x)) y else x
`%NA%` <- function(x, y) ifelse(is.na(x), y, x)

safe_min_date <- function(x) {
  if (length(x) == 0 || all(is.na(x))) as.Date(NA) else suppressWarnings(as.Date(min(x, na.rm = TRUE)))
}

# ---- Config loaders ----------------------------------------------------------
tr_params <- function(yml = here::here("config/tr_response_params.yml")) {
  if (!file.exists(yml)) return(
    list(gap_days=70, pre_grace_days=45, post_grace_days=120,
         keep_only_overlapping=TRUE, count_only_labeled=FALSE,
         max_steps_in_trajectory=12)
  )
  yaml::read_yaml(yml)
}

proj_params <- function(
    yml = Sys.getenv("PROJ_PARAMS", here::here("config/project_params.yml"))
) {
  if (!file.exists(yml)) {
    return(list(
      force_neoadj_empis = character(),
      exclusions = list(),
      files = list(force_neoadj_empis_csv = NULL, exclusions_csv = NULL)
    ))
  }
  p <- yaml::read_yaml(yml)
  
  # normalize IDs
  p$force_neoadj_empis <- unique(as.character(p$force_neoadj_empis %||% character()))
  
  # exclusions -> tibble
  excl <- p$exclusions %||% list()
  p$exclusions_df <- if (length(excl)) {
    tibble::tibble(
      mrn        = vapply(excl, \(x) as.character(x$mrn %||% NA), character(1)),
      medication = vapply(excl, \(x) tolower(as.character(x$medication %||% NA)), character(1))
    ) |> dplyr::filter(!is.na(mrn), !is.na(medication)) |> dplyr::distinct()
  } else tibble::tibble(mrn=character(), medication=character())
  
  # optional CSV merges
  if (!is.null(p$files$force_neoadj_empis_csv) && file.exists(p$files$force_neoadj_empis_csv)) {
    ids <- readr::read_csv(p$files$force_neoadj_empis_csv, show_col_types = FALSE)
    if ("mrn" %in% names(ids)) p$force_neoadj_empis <- unique(c(p$force_neoadj_empis, as.character(ids$mrn)))
  }
  if (!is.null(p$files$exclusions_csv) && file.exists(p$files$exclusions_csv)) {
    ex <- readr::read_csv(p$files$exclusions_csv, show_col_types = FALSE)
    need <- c("mrn","medication")
    if (all(need %in% names(ex))) {
      ex <- ex |> transmute(mrn = as.character(mrn), medication = tolower(as.character(medication)))
      p$exclusions_df <- dplyr::bind_rows(p$exclusions_df, ex) |> dplyr::distinct()
    }
  }
  p
}

# ---- Date parsing ------------------------------------------------------------
parse_date_safely <- function(x) {
  if (inherits(x, "Date"))   return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  if (is.list(x))            x <- unlist(x, use.names = FALSE)
  if (is.numeric(x)) {
    if (is.finite(stats::median(x, na.rm = TRUE)) && stats::median(x, na.rm = TRUE) > 10000)
      return(as.Date(x, origin = "1899-12-30"))  # Excel
    return(as.Date(x, origin = "1970-01-01"))    # Unix days
  }
  if (is.character(x)) {
    xx <- trimws(x)
    suppressWarnings({
      d <- as.Date(xx, "%Y-%m-%d"); if (any(!is.na(d))) return(d)
      d <- as.Date(xx, "%m/%d/%Y"); if (any(!is.na(d))) return(d)
      d <- as.Date(xx, "%d.%m.%Y"); if (any(!is.na(d))) return(d)
    })
    return(as.Date(xx)) # base fallback
  }
  as.Date(as.character(x))
}

pick_and_parse_date <- function(df, candidates = c(
  "Date","DATE","date","ServiceDate","InfusionDate","Infusion_Date",
  "Appt Date","Appt_Date","DOS","Encounter_Date","encounter_date",
  "OpDate","SurgeryDate","OperativeDate"
)) {
  found <- intersect(candidates, names(df))
  if (!length(found)) {
    # fallback to any col containing "date"
    date_like <- grep("date", names(df), ignore.case = TRUE, value = TRUE)
    if (length(date_like)) found <- date_like[1]
  }
  if (!length(found)) stop("No date-like column found in df. Available: ", paste(names(df), collapse = ", "))
  src <- found[1]
  df %>% dplyr::mutate(Date = parse_date_safely(.data[[src]]), .date_src = src)
}

# ---- Medication parsing / reconciliation ------------------------------------
.med_token_dict <- c(
  # PD-1 / CTLA-4 / LAG3
  "nivo"="nivolumab","nivolumab"="nivolumab","opdivo"="nivolumab",
  "pembro"="pembrolizumab","pembrolizumab"="pembrolizumab","keytruda"="pembrolizumab",
  "cemi"="cemiplimab","cemiplimab"="cemiplimab","libtayo"="cemiplimab",
  "ipi"="ipilimumab","ipilimumab"="ipilimumab","yervoy"="ipilimumab",
  "rela"="relatlimab","relatlimab"="relatlimab",
  # chemo / targeted / misc
  "pacli"="paclitaxel","taxol"="paclitaxel","paclitaxel"="paclitaxel",
  "carbo"="carboplatin","carboplatin"="carboplatin",
  "cisplatin"="cisplatin","cddp"="cisplatin",
  "cetux"="cetuximab","erbitux"="cetuximab","cetuximab"="cetuximab",
  "t-vec"="talimogene","tvec"="talimogene","talimogene"="talimogene",
  "dabraf"="dabrafenib","dabrafenib"="dabrafenib",
  "tram"="trametinib","trametinib"="trametinib",
  "encorafenib"="encorafenib","binimetinib"="binimetinib",
  "etoposide"="etoposide","vp16"="etoposide",
  "docetaxel"="docetaxel",
  "cobimetinib"="cobimetinib",
  "cyclophosphamide"="cyclophosphamide","fludarabine"="fludarabine",
  "lifileucel"="lifileucel","il-2"="il-2","il2"="il-2",
  "5fu"="5fu","mito"="mitomycin","mitomycin"="mitomycin",
  "avelumab"="avelumab","atezolizumab"="atezolizumab","aldesleukin"="aldesleukin",
  "xtx101"="xtx101","il-12-abp"="il-12-abp","tebentafusp"="tebentafusp"
)

.noise_patterns <- c(
  "(?i)infusion\\s+nos","(?i)symptom\\s+(assessment|management)",
  "(?i)virtual\\s+appt","(?i)treatment\\s+held","(?i)same\\s*day\\s*add\\s*on",
  "(?i)follow\\s*up","(?i)hsr\\b","(?i)investigational\\s+product\\b"
)

.tokenize <- function(x) {
  if (is.na(x) || !nzchar(x)) return(character(0))
  s <- tolower(trimws(as.character(x)))
  s <- gsub("\\s*(\\+|/|,|&|;|\\-|plus)\\s*", "/", s)
  toks <- unlist(strsplit(s, "/", fixed = TRUE), use.names = FALSE)
  toks <- trimws(toks)
  toks[nzchar(toks)]
}

.map_tokens <- function(tokens) {
  m <- .med_token_dict[tokens]
  known <- !is.na(m)
  unique(m[known])
}

canon_regimen_key <- function(x, keep_unknown = FALSE) {
  s <- tolower(paste(x, collapse = " "))
  if (any(vapply(.noise_patterns, \(p) grepl(p, s, perl = TRUE), logical(1L)))) return(NA_character_)
  toks  <- .tokenize(x)
  gens  <- .map_tokens(toks)
  gens  <- gsub("^cisplatin_etopiside$", "cisplatin/etoposide", gens)
  gens  <- unlist(strsplit(gens, "/", fixed = TRUE), use.names = FALSE)
  gens  <- unique(sort(gens))
  if (!length(gens)) return(NA_character_)
  paste(gens, collapse = " + ")
}

add_regimen_key <- function(df, col) {
  col <- enquo(col)
  df |> mutate(regimen_key = vapply(as.character(!!col), canon_regimen_key, character(1)))
}

make_daily_combo_names <- function(df,
                                   mrn_col  = "MRN",
                                   empi_col = "EMPI",
                                   date_col = "DATE",
                                   med_col  = "medication",
                                   dose_col = "dose") {
  stopifnot(all(c(mrn_col, empi_col, date_col, med_col, dose_col) %in% names(df)))
  
  df %>%
    mutate(
      !!date_col := as.Date(.data[[date_col]]),
      !!med_col  := as.character(.data[[med_col]]),
      !!dose_col := as.character(.data[[dose_col]])
    ) %>%
    group_by(.data[[mrn_col]], .data[[empi_col]], .data[[date_col]]) %>%
    summarise(
      {
        g     <- cur_data_all()
        meds  <- g[[med_col]]
        comps <- unlist(lapply(meds, .tokenize), use.names = FALSE)
        toks  <- sort(unique(.map_tokens(comps)))
        regimen_display <- if (length(toks)) paste(tools::toTitleCase(toks), collapse = " + ") else NA_character_
        tibble::tibble(
          MRN             = g[[mrn_col]][1],
          EMPI            = g[[empi_col]][1],
          tokens          = list(toks),
          regimen_display = regimen_display
        )
      },
      .groups = "drop"
    )
}


# ---- tolerant episodes ↔ responses (nearest start within tol) ---------------
tolerant_left_join_responses <- function(episodes, regsum, tol_days = 3, quiet_many_to_many = TRUE) {
  stopifnot(all(c("EMPI","regimen_key","start_date") %in% names(episodes)))
  stopifnot(all(c("EMPI","regimen_key","start_date") %in% names(regsum)))
  
  epi <- episodes %>%
    mutate(
      start_date = as.Date(start_date),
      end_date   = if ("end_date" %in% names(.)) as.Date(.data$end_date) else as.Date(start_date)
    ) %>%
    arrange(EMPI, regimen_key, start_date) %>%
    mutate(.epi_id = row_number()) %>%
    rename(epi_start_date = start_date)
  
  rs <- regsum %>%
    mutate(start_date = as.Date(start_date)) %>%
    rename(trt_start_date = start_date)
  
  cand <- (if (quiet_many_to_many) suppressWarnings else identity)(
    left_join(epi, rs, by = c("EMPI","regimen_key"))
  ) %>%
    mutate(
      epi_start_date = as.Date(epi_start_date),
      trt_start_date = as.Date(trt_start_date),
      delta_days_num = if_else(is.na(trt_start_date), NA_integer_,
                               abs(as.integer(trt_start_date - epi_start_date)))
    )
  
  best <- cand %>%
    group_by(.epi_id) %>%
    arrange(coalesce(delta_days_num, 2147483647L), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(within_tol = !is.na(delta_days_num) & delta_days_num <= tol_days)
  
  prefer_epi_cols <- c("EMPI","regimen_key",".epi_id","epi_start_date","end_date","n_admins","regimen_display")
  possible_resp_cols <- c(
    "BOR","BOR_date","BOR_evidence","ever_respond","best_response","best_response_date",
    "response_trajectory","response_trajectory_dated",
    "n_notes_total","n_notes_labeled","n_note_dates_total","n_note_dates_labeled",
    "n_PD","n_CR","n_PR","n_SD","bor_diag"
  )
  resp_cols <- intersect(names(best), possible_resp_cols)
  
  matched <- best %>%
    filter(within_tol) %>%
    transmute(
      across(any_of(prefer_epi_cols)),
      across(any_of(resp_cols))
    ) %>%
    rename(start_date = epi_start_date)
  
  unmatched <- epi %>%
    anti_join(select(matched, .epi_id), by = ".epi_id") %>%
    transmute(
      EMPI, regimen_key,
      start_date = epi_start_date,
      across(any_of(c("end_date","n_admins","regimen_display"))),
      !!!setNames(rep(list(NA_character_), sum(resp_cols %in% c(
        "BOR","BOR_evidence","best_response","response_trajectory","response_trajectory_dated","bor_diag"))),
        resp_cols[resp_cols %in% c("BOR","BOR_evidence","best_response","response_trajectory","response_trajectory_dated","bor_diag")]),
      !!!setNames(rep(list(as.Date(NA)), sum(resp_cols %in% c("BOR_date","best_response_date"))),
                  resp_cols[resp_cols %in% c("BOR_date","best_response_date")]),
      !!!setNames(rep(list(NA_integer_), sum(resp_cols %in% c(
        "n_notes_total","n_notes_labeled","n_note_dates_total","n_note_dates_labeled","n_PD","n_CR","n_PR","n_SD"))),
        resp_cols[resp_cols %in% c("n_notes_total","n_notes_labeled","n_note_dates_total","n_note_dates_labeled","n_PD","n_CR","n_PR","n_SD")]),
      !!!setNames(rep(list(NA), sum(resp_cols %in% c("ever_respond"))),
                  resp_cols[resp_cols %in% c("ever_respond")])
    )
  
  bind_rows(matched, unmatched) %>% arrange(EMPI, regimen_key, start_date)
}

# ---- BOR normalization for plan-level joins ---------------------------------
normalize_BOR <- function(x) {
  x0 <- tolower(trimws(as.character(x)))
  x0 <- gsub("\\s+", " ", x0)
  dplyr::case_when(
    grepl("^cr$|^complete( response)?$", x0) ~ "CR",
    grepl("^pr$|^partial( response)?$", x0)  ~ "PR",
    grepl("^sd$|^stable( disease)?$", x0)    ~ "SD",
    grepl("^pd$|^progress(ion|ive)( disease)?$", x0) ~ "PD",
    grepl("^mixed", x0) ~ "Mixed",
    grepl("non.?evaluable|not evaluable|ne", x0) ~ "NE",
    TRUE ~ NA_character_
  )
}

# ---- Linker: doses ↔ surgery + notes (neoadj) --------------------------------
link_immuno_to_surgery <- function(dose_df,
                                   op_df,
                                   notes_df = NULL,
                                   windows = c(45, 90, 180, 365),
                                   neoadj_window_days = 30,
                                   surgery_window_days = 90) {
  suppressPackageStartupMessages({ require(dplyr); require(tidyr); require(fuzzyjoin) })
  
  # ---------------------------------------------------------------------------
  # 1) Dose-side: ensure Date + IDs
  # ---------------------------------------------------------------------------
  doses <- dose_df %>%
    pick_and_parse_date() %>%
    mutate(
      medication = as.character(medication),
      MRN        = as.character(MRN),
      EMPI       = as.character(EMPI)
    ) %>%
    filter(!is.na(Date)) %>%
    select(MRN, EMPI, medication, Date) %>%
    distinct()
  
  all_patients  <- doses %>% distinct(MRN, EMPI)
  all_dose_keys <- doses %>% distinct(MRN, EMPI, medication, Date)
  
  # ---------------------------------------------------------------------------
  # 2) Surgery-side: flexible date picker
  # ---------------------------------------------------------------------------
  surgeries <- op_df %>%
    pick_and_parse_date() %>%
    select(EMPI, Date) %>%
    filter(!is.na(Date)) %>%
    distinct() %>%
    rename(encounter_date = Date)
  
  # ---------------------------------------------------------------------------
  # 3) Non-equi join (dose <= surgery date), pick nearest surgery after dose
  # ---------------------------------------------------------------------------
  surgeries2 <- surgeries %>% rename(EMPI_surg = EMPI)
  dose_to_surg <- fuzzyjoin::fuzzy_left_join(
    doses, surgeries2,
    by = c("EMPI" = "EMPI_surg", "Date" = "encounter_date"),
    match_fun = list(`==`, `<=`)
  ) %>%
    mutate(days_to_surgery = as.integer(encounter_date - Date)) %>%
    group_by(MRN, EMPI, medication, Date) %>%
    summarise(
      next_surgery_date = safe_min_date(encounter_date),
      days_to_surgery   = suppressWarnings(min(days_to_surgery, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      has_surgery_after_this_dose =
        !is.na(next_surgery_date) & !is.na(days_to_surgery) & days_to_surgery >= 0
    )
  
  # Re-expand to include all doses
  dose_to_surg <- all_dose_keys %>%
    left_join(dose_to_surg, by = c("MRN","EMPI","medication","Date"))
  
  # For patients with no EMPI → leave adjudication as NA
  dose_to_surg <- dose_to_surg %>%
    mutate(
      has_surgery_after_this_dose =
        if_else(is.na(EMPI), NA, has_surgery_after_this_dose)
    )
  
  # ---------------------------------------------------------------------------
  # 4) Window flags per dose
  # ---------------------------------------------------------------------------
  if (length(windows)) {
    for (w in windows) {
      nm <- paste0("surgery_within_", w, "d")
      dose_to_surg[[nm]] <- with(
        dose_to_surg,
        has_surgery_after_this_dose & !is.na(days_to_surgery) & days_to_surgery <= w
      )
      dose_to_surg[[nm]] <- if_else(is.na(dose_to_surg$EMPI), NA, dose_to_surg[[nm]])
    }
  }
  
  # ---------------------------------------------------------------------------
  # 5) Patient-level rollup
  # ---------------------------------------------------------------------------
  patient_summary <- dose_to_surg %>%
    group_by(MRN) %>%
    summarise(
      EMPI = first(EMPI[!is.na(EMPI)]),  # keep non-NA EMPI if present
      any_surgery_after_any_dose = any(has_surgery_after_this_dose %in% TRUE, na.rm = TRUE),
      first_surgery_date_after_any_dose = safe_min_date(next_surgery_date),
      earliest_days_to_surgery = {
        v <- suppressWarnings(min(days_to_surgery, na.rm = TRUE))
        if (is.infinite(v)) NA_real_ else v
      },
      .groups = "drop"
    )
  
  # Re-expand to include all patients
  patient_summary <- all_patients %>%
    left_join(patient_summary, by = "MRN") %>%
    mutate(
      EMPI = coalesce(EMPI.x, EMPI.y),
      any_surgery_after_any_dose =
        if_else(is.na(EMPI), NA, coalesce(any_surgery_after_any_dose, FALSE)),
      .keep = "unused"
    )
  
  # Last dose before first surgery
  if (nrow(patient_summary)) {
    patient_summary <- patient_summary %>%
      left_join(
        doses %>% group_by(MRN) %>% summarise(all_doses = list(Date), .groups = "drop"),
        by = "MRN"
      ) %>%
      rowwise() %>%
      mutate(
        last_dose_before_first_surgery = {
          fs <- first_surgery_date_after_any_dose
          if (is.na(fs) || length(all_doses) == 0) NA_Date_ else {
            prior <- all_doses[all_doses <= fs]
            if (length(prior)) max(prior) else NA_Date_
          }
        }
      ) %>%
      ungroup() %>%
      select(-all_doses)
  }
  
  # Surgery window rollups
  if (nrow(dose_to_surg) && any(startsWith(names(dose_to_surg), "surgery_within_"))) {
    win_rollups <- dose_to_surg %>%
      select(MRN, EMPI, starts_with("surgery_within_")) %>%
      group_by(MRN) %>%
      summarise(across(starts_with("surgery_within_"),
                       ~ any(.x %in% TRUE, na.rm = TRUE)), .groups = "drop")
    patient_summary <- patient_summary %>% left_join(win_rollups, by = "MRN")
    # if no EMPI → reset to NA
    patient_summary <- patient_summary %>%
      mutate(across(starts_with("surgery_within_"),
                    ~ if_else(is.na(EMPI), NA, .x)))
  }
  
  # ---------------------------------------------------------------------------
  # 6) Neoadj from notes (optional)
  # ---------------------------------------------------------------------------
  if (!is.null(notes_df)) {
    neo_col <- if ("neoadj" %in% names(notes_df)) "neoadj" else
      if ("neoadjuvant" %in% names(notes_df)) "neoadjuvant" else
        stop("`notes_df` must have logical column `neoadj` or `neoadjuvant`.")
    
    neo_notes <- notes_df %>%
      pick_and_parse_date(candidates = c("encounter_date","Date","DATE")) %>%
      transmute(EMPI, encounter_date = Date, neoadj = as.logical(.data[[neo_col]])) %>%
      filter(!is.na(EMPI), !is.na(encounter_date)) %>%
      distinct()
    
    neo_notes2 <- neo_notes %>% rename(EMPI_note = EMPI)
    dose_neoadj <- fuzzyjoin::fuzzy_left_join(
      doses, neo_notes2,
      by = c("EMPI" = "EMPI_note", "Date" = "encounter_date"),
      match_fun = list(`==`, `>=`)
    ) %>%
      mutate(days_from_note_to_dose = as.integer(Date - encounter_date)) %>%
      filter(!is.na(days_from_note_to_dose),
             days_from_note_to_dose >= 0,
             days_from_note_to_dose <= neoadj_window_days) %>%
      group_by(MRN, EMPI, medication, Date) %>%
      summarise(
        neoadj_near_dose = any(neoadj, na.rm = TRUE),
        nearest_note_days = suppressWarnings(min(days_from_note_to_dose, na.rm = TRUE)),
        .groups = "drop"
      )
    
    dose_to_surg <- dose_to_surg %>%
      left_join(dose_neoadj, by = c("MRN","EMPI","medication","Date")) %>%
      mutate(
        neoadj_near_dose = if_else(is.na(EMPI), NA, coalesce(neoadj_near_dose, FALSE))
      )
    
    likely_flag_name <- paste0(
      "likely_neoadjuvant__note_", neoadj_window_days,
      "d_before__surg_", surgery_window_days, "d_after"
    )
    dose_to_surg[[likely_flag_name]] <- with(
      dose_to_surg,
      neoadj_near_dose %in% TRUE & has_surgery_after_this_dose %in% TRUE &
        !is.na(days_to_surgery) & days_to_surgery <= surgery_window_days
    )
    dose_to_surg[[likely_flag_name]] <- if_else(is.na(dose_to_surg$EMPI),
                                                NA, dose_to_surg[[likely_flag_name]])
    
    # Patient-level neoadj rollups
    patient_summary <- patient_summary %>%
      left_join(
        dose_to_surg %>%
          group_by(MRN) %>%
          summarise(
            neoadj_any_near_dose = any(neoadj_near_dose %in% TRUE, na.rm = TRUE),
            !!likely_flag_name   := any(.data[[likely_flag_name]] %in% TRUE, na.rm = TRUE),
            .groups = "drop"
          ),
        by = "MRN"
      ) %>%
      mutate(
        neoadj_any_near_dose = if_else(is.na(EMPI), NA, neoadj_any_near_dose),
        !!likely_flag_name   := if_else(is.na(EMPI), NA, .data[[likely_flag_name]])
      )
  }
  
  # ---------------------------------------------------------------------------
  # Return
  # ---------------------------------------------------------------------------
  list(
    dose_level    = dose_to_surg %>% arrange(MRN, EMPI, Date, medication),
    patient_level = patient_summary %>% arrange(MRN)
  )
}



# ---- Plan ↔ response nearest-date join (plan_start to episode start) --------
nearest_join_plan_response <- function(plans, episodes, tol_days = 7L, keep_unmatched_plans = TRUE) {
  stopifnot(all(c("EMPI","regimen_key","plan_start") %in% names(plans)))
  stopifnot(all(c("EMPI","regimen_key","start_date") %in% names(episodes)))
  
  cand <- dplyr::left_join(
    plans %>% dplyr::mutate(.plan_id = dplyr::row_number(),
                            plan_start = as.Date(plan_start)),
    episodes %>% dplyr::mutate(.epi_id = dplyr::row_number(),
                               start_date = as.Date(start_date)),
    by = c("EMPI","regimen_key")
  ) %>%
    dplyr::mutate(delta_days = abs(as.integer(start_date - plan_start)))
  
  best <- cand %>%
    dplyr::group_by(.plan_id) %>%
    dplyr::arrange(dplyr::coalesce(delta_days, 2147483647L), .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(within_tol = !is.na(delta_days) & delta_days <= tol_days)
  
  if (keep_unmatched_plans) best else dplyr::filter(best, within_tol)
}

# ---- Output glossary ---------------------------------------------------------
make_out_glossary <- function(out) {
  suppressPackageStartupMessages({ library(dplyr); library(purrr); library(stringr); library(tidyr) })
  build <- function(df, dataset_name) {
    tibble(variable = names(df)) %>%
      mutate(
        type = map_chr(variable, ~ paste(class(df[[.x]]), collapse = ", ")),
        description = dplyr::case_when(
          variable == "EMPI"  ~ "Patient identifier (EMPI).",
          variable == "medication" ~ "Systemic therapy administered (lowercase).",
          variable == "Date"  ~ "Date of the immunotherapy dose.",
          variable == "next_surgery_date" ~ "Nearest surgery date on/after the dose date.",
          variable == "days_to_surgery"    ~ "Days from dose to nearest subsequent surgery.",
          variable == "has_surgery_after_this_dose" ~ "TRUE if a surgery occurred on/after the dose date.",
          variable == "neoadj_near_dose" ~ "TRUE if a neoadjuvant note within the pre-dose window.",
          variable == "nearest_note_days" ~ "Days between qualifying note and dose.",
          variable == "any_surgery_after_any_dose" ~ "Patient-level: any dose followed by surgery.",
          variable == "first_surgery_date_after_any_dose" ~ "Earliest surgery date after any dose.",
          variable == "earliest_days_to_surgery" ~ "Min days from dose to subsequent surgery.",
          variable == "last_dose_before_first_surgery" ~ "Latest dose on/before first post-dose surgery.",
          TRUE ~ NA_character_
        ),
        description = if_else(
          is.na(description) & str_starts(variable, "surgery_within_"),
          paste0("TRUE if a surgery occurred within ",
                 stringr::str_match(variable, "surgery_within_(\\d+)d$")[,2],
                 " days after this dose."),
          description
        ),
        description = if_else(
          is.na(description) & str_detect(variable, "^likely_neoadjuvant__note_\\d+d_before__surg_\\d+d_after$"),
          "TRUE if note-within-window AND surgery-within-window.",
          description
        ),
        description = tidyr::replace_na(description, "Derived variable (see pipeline)."),
        dataset = dataset_name,
        .before = 1
      )
  }
  dplyr::bind_rows(build(out$dose_level, "dose_level"), build(out$patient_level, "patient_level")) %>%
    arrange(dataset, variable)
}

suppressPackageStartupMessages({
  library(shiny)
  library(dplyr)
  library(readr)
  library(here)
})

# Path to override directory
override_dir <- here::here("files/manual_overrides")
dir.create(override_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# Load all overrides (coalesce multiple RDS files, enforce schema)
# -------------------------------------------------------------------
load_all_overrides <- function(directory = here::here("files/manual_overrides")) {
  files <- list.files(directory, pattern = "med_overrides.*\\.rds$", full.names = TRUE)
  
  if (length(files) == 0) {
    return(tibble(
      MRN                   = character(),
      Patient               = character(),
      medication            = character(),
      neoadj_intent_manual  = character(),
      borderline_manual     = character(),
      surgery_manual        = character(),
      surgery_date_manual   = as.Date(character()),
      trt_simple_med_manual = character(),
      trt_med_type5_manual  = character(),
      manual_override_flag  = logical(),
      timestamp             = as.POSIXct(character())
    ))
  }
  
  overrides_list <- lapply(files, readRDS)
  overrides <- bind_rows(overrides_list)
  
  # --- Enforce schema ---
  expected_cols <- c(
    "MRN","Patient","medication",
    "neoadj_intent_manual","borderline_manual","surgery_manual",
    "surgery_date_manual","trt_simple_med_manual","trt_med_type5_manual",
    "manual_override_flag","timestamp"
  )
  for (col in expected_cols) {
    if (!col %in% names(overrides)) overrides[[col]] <- NA
  }
  
  # --- Coerce types ---
  overrides <- overrides %>%
    mutate(
      MRN                   = as.character(MRN),
      Patient               = as.character(Patient),
      medication            = as.character(medication),
      neoadj_intent_manual  = as.character(neoadj_intent_manual),
      borderline_manual     = as.character(borderline_manual),
      surgery_manual        = as.character(surgery_manual),
      surgery_date_manual   = as.Date(surgery_date_manual),
      trt_simple_med_manual = as.character(trt_simple_med_manual),
      trt_med_type5_manual  = as.character(trt_med_type5_manual),
      manual_override_flag  = as.logical(manual_override_flag),
      timestamp             = as.POSIXct(timestamp)
    )
  overrides <- overrides %>%
    mutate(
      timestamp = as.POSIXct(timestamp, origin = "1970-01-01", tz = "UTC")
    )
  
  
  # --- Deduplicate (keep most recent per MRN+medication) ---
  overrides %>%
    arrange(MRN, medication, desc(timestamp)) %>%
    distinct(MRN, medication, .keep_all = TRUE)
}



# -------------------------------------------------------------------
# Save overrides (always enforce schema)
# -------------------------------------------------------------------
save_overrides <- function(df, directory = here::here("files/manual_overrides")) {
  fs::dir_create(directory)
  file_date <- format(Sys.time(), "%Y-%m-%d-%H-%M-%S")
  file_name_save <- file.path(directory, paste0("med_overrides-", file_date, ".rds"))
  
  df <- df %>%
    mutate(
      MRN                   = as.character(MRN),
      Patient               = as.character(Patient),
      medication            = as.character(medication),
      neoadj_intent_manual  = as.character(neoadj_intent_manual),
      borderline_manual     = as.character(borderline_manual),
      surgery_manual        = as.character(surgery_manual),
      surgery_date_manual   = as.Date(surgery_date_manual),
      trt_simple_med_manual = as.character(trt_simple_med_manual),
      trt_med_type5_manual  = as.character(trt_med_type5_manual),
      manual_override_flag  = as.logical(manual_override_flag),
      timestamp             = as.POSIXct(timestamp, origin = "1970-01-01", tz = "UTC")
    )
  
  saveRDS(df, file_name_save)
  message("Overrides saved to: ", file_name_save)
  invisible(file_name_save)
}


           # =============================================================================
# CSCC Systemic Therapy — Neoadjuvant pipeline EXECUTE
# (analysis build + calls to 80_neo_visuals.R for plots)
# =============================================================================

suppressPackageStartupMessages({
  source(here::here("scripts/Functions and Packages/load_packages.R"))
  source(here::here("scripts/Treatment Response Assessment Scripts/10_mcode_tr_response_functions.R"))
  source(here::here("scripts/Infusion Related Scripts/neoadjuvant/10_neo_functions.R"))
  source(here::here("scripts/Infusion Related Scripts/neoadjuvant/80_neo_visuals.R"))
  source(here::here("scripts/Infusion Related Scripts/neoadjuvant/app_neoadj_intent_overrides.R"))
  library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(ggplot2); library(ggpattern); library(scales)
})

# ------------------------------ Parameters -----------------------------------
tr_p   <- tr_params()
proj_p <- proj_params()
gap_days <- tr_p$gap_days %||% 70L

# --------------------------- Helper: safe save -------------------------------
save_png <- function(plot, name, width = 12, height = 7, dpi = 300, units = "in") {
  dir.create(here::here("files","data visualizations"), showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(
    filename = here::here("files","data visualizations", paste0(name, ".png")),
    plot     = plot,
    width    = width, height = height, dpi = dpi, units = units
  )
}

# -----------------------------------------------------------------------------
# 1) LOAD, NORMALIZE (dose rows → daily combos)
# -----------------------------------------------------------------------------
message("• Loading infusion dataset...")
systemic_tx_load <- open_recent_file(directory = "files/Infusion Data/cleaned")
if (!inherits(systemic_tx_load$Date, "Date")) systemic_tx_load$Date <- as.Date(systemic_tx_load$Date)

# QC dmm's patients
dmm_systemic_tx_load <- systemic_tx_load |> filter(`Primary Oncologist` =="David Miller")

systemic_TX <- systemic_tx_load %>%
  transmute(
    Patient,
    EMPI = canon_empi(EMPI),
    MRN  = as.character(MRN),
    DATE = as.Date(Date),
    medication,
    dose
  ) %>%
  filter(!is.na(MRN), nzchar(EMPI))

start_date <- suppressWarnings(min(systemic_TX$DATE, na.rm = TRUE))
last_date  <- suppressWarnings(max(systemic_TX$DATE, na.rm = TRUE))

message("• Building daily regimen display labels...")
daily_regimens <- make_daily_combo_names(systemic_TX)

systemic_TX_with_combo <- 
  systemic_TX %>% 
  select(-EMPI) |> 
  left_join(daily_regimens, by = c("MRN","DATE")) %>%
  mutate(medication_std = coalesce(regimen_display, medication))

# collapse to distinct dose rows
systemic_tx_pre <- systemic_TX_with_combo %>%
  group_by(MRN, DATE, medication_std) %>%
  distinct() %>%
  ungroup()

# carry labels from the loader for later slicing
systemic_tx <- systemic_tx_pre %>%
  left_join(
    systemic_tx_load %>%
      group_by(MRN, Date) %>% slice_head() %>% ungroup() %>%
      select(Patient, MRN, EMPI, `Encounter Provider`, `Primary Oncologist`,
             `Primary Diagnosis`, Date, `Appt Time`),
    by = join_by(MRN, DATE == Date, Patient, EMPI)
  )

# -----------------------------------------------------------------------------
# 2) EPISODES (gap rule) + RECONCILE to regimen_summary (tolerant nearest-date)
# -----------------------------------------------------------------------------
message("• Building regimen keys and gap-based episodes...")
systemic_tx_keys_daily <- systemic_tx %>%
  mutate(DATE = as.Date(DATE)) %>%
  add_regimen_key(medication_std) %>%
  distinct(MRN, DATE, regimen_key, .keep_all = TRUE)

systemic_keys_episodes <- 
  systemic_tx_keys_daily %>%
  filter(!is.na(regimen_key)) %>%
  arrange(MRN, regimen_key, DATE) %>%
  group_by(MRN, regimen_key) %>%
  mutate(
    gap     = as.integer(DATE - lag(DATE)),
    new_ep  = if_else(is.na(gap) | gap > gap_days, 1L, 0L),
    episode_id = cumsum(replace_na(new_ep, 1L))
  ) %>%
  group_by(Patient, MRN, regimen_key, episode_id, EMPI) %>%
  summarise(
    start_date      = min(DATE),
    end_date        = max(DATE),
    n_admins        = n_distinct(DATE),
    regimen_display = first(na.omit(medication_std)),
    .groups = "drop"
  )

message("• Loading regimen_summary and reconciling to episodes...")
regimen_summary <- open_recent_file(directory = here::here("files/treatment response assessments"))

regimen_summary_keys <- 
  regimen_summary %>%
  mutate(start_date = as.Date(start_date)) %>%
  add_regimen_key(regimen)

systemic_with_responses <- 
  tolerant_left_join_responses(
  episodes = systemic_keys_episodes,
  regsum   = regimen_summary_keys,
  tol_days = 3L
)

# -----------------------------------------------------------------------------
# 3) CSCC subset, exclusions, normalized med labels
# -----------------------------------------------------------------------------
message("• Creating CSCC subset and applying exclusions...")
excluded_meds <- c(
  "investigational product","infusion nos","sonidegib","vismodegib",
  "symptom assessment","treatment held - symptom management",
  "same day add on - symptom management","follow up","virtual appt no charge",
  "hsr - mild","hsr - moderate"
)

.valid_meds <- c(
  "XTX101","IL-12-ABP","Paclitaxel","Ipilimumab","Nivolumab","Talimogene",
  "Cisplatin (CDDP)","Cetuximab","Pembrolizumab","Cemiplimab","Other"
)

norm_med_display <- function(x) {
  y <- stringr::str_squish(x)
  y <- stringr::str_replace_all(y, "\\s*-\\s*", "-")
  rules <- list(
    "^xtx\\s*101$"="XTX101","^il\\s*-?\\s*12\\s*-?\\s*abp$"="IL-12-ABP",
    "^paclitaxel$"="Paclitaxel","^ipilimumab$"="Ipilimumab",
    "^nivolumab(\\s*\\(iv/?sc\\)\\s*240\\s*mg)?$"="Nivolumab",
    "^talimogene$"="Talimogene","^cetuximab$"="Cetuximab",
    "^pembrolizumab$"="Pembrolizumab","^cemiplimab$"="Cemiplimab",
    "^cisplatin(\\s*\\(cddp\\))?$"="Cisplatin (CDDP)"
  )
  for (rx in names(rules)) y <- stringr::str_replace(y, stringr::regex(rx, TRUE), rules[[rx]])
  y[!y %in% .valid_meds] <- "Other"
  y
}

cscc <- systemic_tx %>% filter(`Primary Diagnosis` == "CSCC")

cscc_filter <- cscc %>%
  # 1. Remove rows where the medication name (case-insensitive)
  #    is in the predefined `excluded_meds` list (e.g., "investigational product",
  #    "symptom assessment", etc.).
  filter(!(tolower(medication) %in% tolower(excluded_meds))) %>%
  
  # 2. Remove any records that appear in the project-level
  #    exclusions table (`proj_p$exclusions_df`), matching by MRN and medication.
  #    (anti_join = keep everything that does NOT match).
  anti_join(proj_p$exclusions_df,
            by = c("MRN" = "mrn", "medication" = "medication")) %>%
  
  # 3. Normalize the medication display names to a small set of
  #    canonical labels (Nivolumab, Cemiplimab, Cetuximab, etc.).
  #    Anything not recognized is mapped to "Other".
  mutate(medication = norm_med_display(medication))


# -----------------------------------------------------------------------------
# 4) PROGRESS NOTES & OP NOTES → dose→surgery linking + plan-level intent
# -----------------------------------------------------------------------------
message("• Linking doses to surgeries and computing neoadjuvant windows...")
source(here::here("scripts/Infusion Related Scripts/neoadjuvant/Neoadjuvant Strings from Progress Notes.R"))
progress_note_rpdr <- open_recent_file(
  directory = "~/Partners HealthCare Dropbox/David Miller/mLab/Projects/Text Mining Skin Cancer Patients from RPDR/files/progress notes"
)
skin.cancer.providers <- readRDS(
  file = here::here("files/providers/Skin Cancer Provider Code Look Up Table.rds")
)

cscc_filter_patients <- cscc_filter %>% group_by(MRN) %>% slice_head() %>% ungroup()
cscc_pn <- progress_note_rpdr %>% filter(EMPI %in% cscc_filter_patients$EMPI)
rm(progress_note_rpdr); gc()

cscc_pn_skin <- cscc_pn %>% filter(provider %in% skin.cancer.providers$`Provider/Resource`)
rm(cscc_pn); gc()

neoadj_flag_notes <- cscc_pn_skin %>%
  mutate(
    neoadj = str_detect(text, regex("neoadjuvant|preoperative|pre-operative", TRUE)),
    neoadj_200char_context = purrr::map_chr(text, ~ extract_neoadj_snippets_safely(.x, pre = 100, post = 100, highlight = TRUE)),
    neoadj_terms_exact     = purrr::map_chr(text,   extract_neoadj_terms_exact_safely)
  ) %>%
  select(-text)

op_notes_rpdr <- open_recent_file(
  directory = "/Users/davidmiller/Partners HealthCare Dropbox/David Miller/mLab/Projects/Text Mining Skin Cancer Patients from RPDR/files/op notes"
)

out_link <- link_immuno_to_surgery(
  dose_df = cscc_filter,
  op_df   = op_notes_rpdr,
  notes_df= neoadj_flag_notes,
  windows = c(45, 90, 180, 365),
  neoadj_window_days = 30,
  surgery_window_days = 90
)

dose_level_pre    <- out_link$dose_level
patient_level_pre <- out_link$patient_level
glossary      <- make_out_glossary(out_link)

# Ensure all CSCC doses are present in dose_level (even if EMPI is NA)
dose_level <- cscc_filter %>%
  transmute(Patient, MRN, EMPI, medication, Date = DATE) %>% 
  distinct() %>%
  left_join(dose_level_pre, by = c("MRN","EMPI","medication","Date")) %>%
  mutate(
    neoadj_near_dose            = coalesce(neoadj_near_dose, FALSE),
    has_surgery_after_this_dose = coalesce(has_surgery_after_this_dose, FALSE),
    across(starts_with("surgery_within_"), ~ coalesce(.x, FALSE))
  )

# Ensure all CSCC patients are present in patient_level
patient_level <- cscc_filter %>%
  distinct(MRN, EMPI) %>%
  left_join(patient_level_pre, by = c("MRN","EMPI")) %>%
  mutate(
    any_surgery_after_any_dose = coalesce(any_surgery_after_any_dose, FALSE),
    earliest_days_to_surgery = {
      v <- suppressWarnings( if (is.null(earliest_days_to_surgery)) NA_real_ else earliest_days_to_surgery )
      ifelse(is.infinite(v), NA_real_, v)
    },
    across(starts_with("surgery_within_"), ~ coalesce(.x, FALSE)),
    across(starts_with("likely_neoadjuvant__"), ~ coalesce(.x, FALSE)),
    neoadj_any_near_dose = coalesce(neoadj_any_near_dose, FALSE)
  )


# Plan-level rollup (EMPI × medication)
safe_min_date <- function(x) if (length(x) == 0 || all(is.na(x))) as.Date(NA) else suppressWarnings(min(x, na.rm = TRUE))
likely_col <- grep("^likely_neoadjuvant__", names(dose_level), value = TRUE)[1]
has_likely <- length(likely_col) == 1 && !is.na(likely_col)

med_intent <- dose_level %>%
  group_by(Patient, MRN, medication) %>%
  summarise(
    EMPI = first(EMPI[!is.na(EMPI)]),   # keep one EMPI if present
    doses_total   = n(),
    first_dose    = safe_min_date(Date),
    doses_neoadj  = sum(neoadj_near_dose %in% TRUE, na.rm = TRUE),
    any_neoadj    = any(neoadj_near_dose %in% TRUE, na.rm = TRUE),
    surg_after_neoadj = any((neoadj_near_dose %in% TRUE) & (has_surgery_after_this_dose %in% TRUE), na.rm = TRUE),
    likely_neoadj = if (has_likely) any(.data[[likely_col]] %in% TRUE, na.rm = TRUE) else
      any((neoadj_near_dose %in% TRUE) & (has_surgery_after_this_dose %in% TRUE) &
            !is.na(days_to_surgery) & days_to_surgery <= 45, na.rm = TRUE),
    first_neoadj_dose = safe_min_date(Date[neoadj_near_dose %in% TRUE]),
    first_surg_after  = safe_min_date(next_surgery_date[neoadj_near_dose %in% TRUE]),
    .groups = "drop"
  ) %>%
  mutate(
    neoadj_intent_med           = any_neoadj,
    neoadj_intent_plus_surg_med = likely_neoadj | surg_after_neoadj,
    neoadj_intent_no_surg_med   = neoadj_intent_med & !neoadj_intent_plus_surg_med
  )

# overrides from project params
forced <- tibble::tibble(EMPI = proj_p$force_neoadj_empis, medication = NA_character_)
med_intent <- med_intent %>%
  mutate(is_forced_neoadj =
           EMPI %in% forced$EMPI |
           paste(EMPI, medication) %in%
           paste(forced$EMPI[!is.na(forced$medication)],
                 forced$medication[!is.na(forced$medication)])) %>%
  mutate(neoadj_intent_med = neoadj_intent_med | is_forced_neoadj)

# final plan table for plotting
# -----------------------------------------------------------------------------
# Merge manual overrides into med_episodes
# -----------------------------------------------------------------------------

# 1. Create med_episodes from pipeline logic
med_episodes <- med_intent %>%
  mutate(
    trt_simple_med = if_else(neoadj_intent_med, "Neoadjuvant", "Definitive"),
    neoadj_intent_mod_med = if_else(
      is_forced_neoadj, TRUE,
      if_else(doses_total > 4, FALSE, neoadj_intent_med)
    ),
    trt_med_type5 = case_when(
      neoadj_intent_mod_med  &  neoadj_intent_plus_surg_med ~ "Neoadjuvant Intent - Surgery",
      neoadj_intent_mod_med  & !neoadj_intent_plus_surg_med ~ "Neoadjuvant Intent - No Surgery",
      !neoadj_intent_mod_med &  neoadj_intent_med &  neoadj_intent_plus_surg_med ~ "Borderline Neoadjuvant Intent - Surgery",
      !neoadj_intent_mod_med &  neoadj_intent_med & !neoadj_intent_plus_surg_med ~ "Borderline Neoadjuvant Intent - No Surgery",
      TRUE ~ "Definitive Treatment"
    )
  ) %>%
  arrange(desc(doses_total), MRN, medication) %>%
  mutate(Plan_rank = row_number())

# 2. Load + coalesce all override files (keeps most recent row per MRN+medication)
# -----------------------------------------------------------------------------
# Merge manual overrides into med_episodes (character-based schema)
# -----------------------------------------------------------------------------

# Load + coalesce all override files (keeps most recent per MRN+medication)
overrides <- load_all_overrides(here::here("files/manual_overrides"))

# --- enforce schema ---
expected_cols <- c(
  "MRN","Patient","medication",
  "neoadj_intent_manual","borderline_manual","surgery_manual",
  "surgery_date_manual","trt_simple_med_manual","trt_med_type5_manual",
  "manual_override_flag","timestamp"
)

for (col in expected_cols) {
  if (!col %in% names(overrides)) overrides[[col]] <- NA
}

# --- coerce types (all character flags stay as Yes/No) ---
overrides <- overrides %>%
  mutate(
    MRN                   = as.character(MRN),
    Patient               = as.character(Patient),
    medication            = as.character(medication),
    neoadj_intent_manual  = as.character(neoadj_intent_manual),
    borderline_manual     = as.character(borderline_manual),
    surgery_manual        = as.character(surgery_manual),
    surgery_date_manual   = as.Date(surgery_date_manual),
    trt_simple_med_manual = as.character(trt_simple_med_manual),
    trt_med_type5_manual  = as.character(trt_med_type5_manual),
    manual_override_flag  = as.logical(manual_override_flag),
    timestamp             = as.POSIXct(timestamp)
  )

# --- merge into med_episodes ---
med_episodes <- med_episodes %>%
  left_join(overrides, by = c("Patient","MRN","medication")) %>%
  mutate(
    # Manual always takes precedence
    neoadj_intent_med    = case_when(
      neoadj_intent_manual == "Yes" ~ TRUE,
      neoadj_intent_manual == "No"  ~ FALSE,
      TRUE ~ neoadj_intent_med
    ),
    surg_after_neoadj    = case_when(
      surgery_manual == "Yes" ~ TRUE,
      surgery_manual == "No"  ~ FALSE,
      TRUE ~ surg_after_neoadj
    ),
    first_surg_after     = coalesce(surgery_date_manual, first_surg_after),
    trt_simple_med       = coalesce(trt_simple_med_manual, trt_simple_med),
    trt_med_type5        = coalesce(trt_med_type5_manual, trt_med_type5),
    manual_override_flag = coalesce(manual_override_flag, FALSE)
  )


# -----------------------------------------------------------------------------
# 4. Launch the manual override app (optional, interactive step)
# -----------------------------------------------------------------------------

neoadj_intent_overrides(med_episodes)

# 2. Load + coalesce all override files (keeps most recent row per MRN+medication)
overrides <- load_all_overrides(here::here("files/manual_overrides"))
# --- enforce schema ---
expected_cols <- c(
  "MRN","Patient","medication",
  "neoadj_intent_manual","borderline_manual","surgery_manual",
  "surgery_date_manual","trt_simple_med_manual","trt_med_type5_manual",
  "manual_override_flag","timestamp"
)

for (col in expected_cols) {
  if (!col %in% names(overrides)) overrides[[col]] <- NA
}
# --- coerce types (all character flags stay as Yes/No) ---
overrides <- overrides %>%
  mutate(
    MRN                   = as.character(MRN),
    Patient               = as.character(Patient),
    medication            = as.character(medication),
    neoadj_intent_manual  = as.character(neoadj_intent_manual),
    borderline_manual     = as.character(borderline_manual),
    surgery_manual        = as.character(surgery_manual),
    surgery_date_manual   = as.Date(surgery_date_manual),
    trt_simple_med_manual = as.character(trt_simple_med_manual),
    trt_med_type5_manual  = as.character(trt_med_type5_manual),
    manual_override_flag  = as.logical(manual_override_flag),
    timestamp             = as.POSIXct(timestamp)
  )

# --- merge into med_episodes ---
override_only_cols <- c(
  "neoadj_intent_manual","borderline_manual","surgery_manual",
  "surgery_date_manual","trt_simple_med_manual","trt_med_type5_manual",
  "manual_override_flag","timestamp"
)

med_episodes <- med_episodes %>%
  select(-any_of(override_only_cols)) %>%   # keep Patient/MRN/medication intact
  left_join(overrides, by = c("Patient","MRN","medication")) %>%
  mutate(
    neoadj_intent_med = case_when(
      neoadj_intent_manual == "Yes" ~ TRUE,
      neoadj_intent_manual == "No"  ~ FALSE,
      TRUE ~ neoadj_intent_med
    ),
    surg_after_neoadj = case_when(
      surgery_manual == "Yes" ~ TRUE,
      surgery_manual == "No"  ~ FALSE,
      TRUE ~ surg_after_neoadj
    ),
    first_surg_after     = coalesce(surgery_date_manual, first_surg_after),
    trt_simple_med       = coalesce(trt_simple_med_manual, trt_simple_med),
    trt_med_type5        = coalesce(trt_med_type5_manual, trt_med_type5),
    manual_override_flag = coalesce(manual_override_flag, FALSE)
  )
# Clean trt_simple_med
med_episodes <- med_episodes %>%
  mutate(
    trt_simple_med = case_when(
      str_detect(trt_med_type5, "Neoadjuvant Intent") ~ "Neoadjuvant",
      str_detect(trt_med_type5, "Borderline Neoadjuvant Intent") ~ "Neoadjuvant",
      trt_med_type5 == "Definitive Treatment" ~ "Definitive",
      TRUE ~ trt_simple_med   # fallback if new categories appear
    )
  )
# if there are duplicate rows and one has an EMPI, keep that one
med_episodes <- med_episodes %>%
  group_by(MRN) %>%
  # if any row has a non-NA EMPI, keep only those
  filter(!(is.na(EMPI) & any(!is.na(EMPI)))) %>%
  ungroup()

# -----------------------------------------------------------------------------
# 5) PLAN ↔ RESPONSE pairing (for responder-pattern waterfall)
# -----------------------------------------------------------------------------
message("• Nearest plan↔response pairing...")
plan_intent <- med_episodes %>%
  mutate(
    regimen_key = vapply(medication, canon_regimen_key, character(1)),
    plan_start  = as.Date(first_dose)
  ) %>%
  select(Patient, EMPI, MRN, medication, doses_total, trt_simple_med, trt_med_type5, regimen_key, plan_start) #%>%
  #filter(!is.na(regimen_key))
# Save this data frame (to be used to screen for neoadjuvant treated patients)
save_files(
  save_object = plan_intent,
  subD = "files/Infusion Data/neoadjuvant/cscc/plans with intent"
)

plan_with_resp <- nearest_join_plan_response(
  plans    = plan_intent,
  episodes = systemic_with_responses,
  tol_days = 7L,
  keep_unmatched_plans = TRUE
) %>%
  mutate(
    BOR_std = normalize_BOR(BOR),
    responder_bor = BOR_std %in% c("CR","PR"),
    responder_any = dplyr::coalesce(as.logical(ever_respond), responder_bor)
  )

# Neo-only plans + response binning for pattern overlay
neo_levels <- c("Neoadjuvant Intent - Surgery", "Neoadjuvant Intent - No Surgery")

plan_with_resp_bin <- plan_with_resp %>%
  mutate(
    resp_bin = case_when(
      BOR_std %in% c("CR","PR")           ~ "Responder (CR/PR)",
      BOR_std %in% c("SD","PD","Mixed")   ~ "Non-responder (SD/PD/Mixed)",
      TRUE                                ~ "Unknown"
    ),
    resp_bin = factor(resp_bin,
                      levels = c("Responder (CR/PR)", "Non-responder (SD/PD/Mixed)", "Unknown"))
  )

neo_plans_bin <- med_episodes %>%
  filter(trt_med_type5 %in% neo_levels) %>%
  mutate(
    regimen_key = vapply(medication, canon_regimen_key, character(1)),
    plan_start  = as.Date(first_dose)
  ) %>%
  left_join(
    plan_with_resp_bin %>% select(EMPI, regimen_key, plan_start, resp_bin),
    by = c("EMPI","regimen_key","plan_start")
  ) %>%
  arrange(desc(doses_total), EMPI, medication) %>%
  mutate(Neo_rank = row_number())

# Histogram data (counts & proportions per doses_total bin)
neo_hist <- med_episodes %>%
  filter(trt_med_type5 %in% neo_levels) %>%
  transmute(
    EMPI, medication,
    doses_total = as.integer(doses_total),
    intent = trt_med_type5
  ) %>%
  count(doses_total, intent, name = "plans") %>%
  group_by(doses_total) %>%
  mutate(
    total_in_bin = sum(plans),
    pct          = plans / total_in_bin
  ) %>%
  ungroup()


# ---- Summarize doses and patients by medication -----------------------------
tx_summary <- cscc_filter %>%
  dplyr::group_by(medication) %>%
  dplyr::summarise(
    total_doses     = dplyr::n(),
    unique_patients = dplyr::n_distinct(MRN),
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(total_doses)) %>%
  dplyr::mutate(
    percent_doses    = round(total_doses / sum(total_doses) * 100, 1),
    percent_patients = round(unique_patients / sum(unique_patients) * 100, 1),
    label = paste0(scales::comma(total_doses), " (", percent_doses, "%)")
  )

# -----------------------------------------------------------------------------
# 6) PLOTS (call the utilities in 80_neo_visuals.R)
# -----------------------------------------------------------------------------
message("• Rendering plots via 80_neo_visuals.R ...")
#----Tx Plan---------------------
tx_plot        <- plot_tx_by_med(tx_summary, start_date, last_date)
tx_plot
tx_nested_plot <- plot_tx_nested(tx_summary, start_date, last_date)
tx_nested_plot
#----Time Series-----------------
# compute one shared y bound (once) and pass to both
ts_ymax <- max(med_intent %>% arrange(as.Date(first_dose)) %>% mutate(n = row_number()) %>% pull(n))

p_ts_all  <- plot_plan_timeseries_all(med_intent,  y_max = ts_ymax)
p_ts_all

p_ts_both <- plot_plan_timeseries_both(med_intent, y_max = ts_ymax)
p_ts_both

# ---------------- STRICT NEODJ (exclude borderlines) ----------------
neo_levels_strict <- c("Neoadjuvant Intent - Surgery",
                       "Neoadjuvant Intent - No Surgery")

neo_strict <- med_episodes %>%
  dplyr::filter(trt_med_type5 %in% neo_levels_strict) %>%
  dplyr::transmute(
    EMPI,
    first_dose = as.Date(first_dose)
  ) %>%
  dplyr::arrange(first_dose) %>%
  dplyr::mutate(Neo_Seq = dplyr::row_number())

# --- optional single-curve viz just for strict neoadjuvant (no borderlines) ---
plot_plan_timeseries_neo_strict <- function(neo_df, y_max = NULL, pad = 0.05) {
  stopifnot(all(c("first_dose","Neo_Seq") %in% names(neo_df)))
  
  max_val <- max(neo_df$Neo_Seq, na.rm = TRUE)
  y_max   <- y_max %||% max_val * (1 + pad)   # add a little headroom
  
  ggplot2::ggplot(neo_df, ggplot2::aes(first_dose, Neo_Seq)) +
    ggplot2::geom_point(size = 2, color = "firebrick") +
    ggplot2::scale_y_continuous(
      limits = c(0, y_max),
      breaks = seq(0, y_max, by = 10),   # adjust step size here
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::labs(
      title = "Neoadjuvant Systemic Therapy Plans Over Time",
      subtitle = paste0("(", "n = ", nrow(neo_df), ")"),
      x = "Date", y = "Cumulative Count"
    ) +
    neo_base_theme() +
    ggplot2::theme(legend.position = "none")
}

# Render ONLY if you want the strict-neo curve as its own figure
p_ts_neo_strict <- plot_plan_timeseries_neo_strict(neo_strict)
p_ts_neo_strict
# Window: inclusive June 1, 2024 through June 30, 2025
date_min <- as.Date("2024-06-01")
date_max <- as.Date("2025-05-31")

strict_neo_counts <- neo_strict %>%
  dplyr::filter(first_dose >= date_min, first_dose <= date_max) %>%
  dplyr::summarise(
    patients  = dplyr::n_distinct(EMPI),
    plans     = dplyr::n(),
    .groups = "drop"
  )

strict_neo_counts

# --- Waterfalls:  -------------------------------
p_waterfall_mono <- plot_waterfall_mono(med_episodes) + 
  scale_x_continuous(breaks = seq(0, 310, by = 50))
p_waterfall_mono 

p_waterfall_binary <- plot_waterfall_binary(med_episodes) + 
  scale_x_continuous(breaks = seq(0, 310, by = 50))
p_waterfall_binary
p_waterfall_5   <- plot_waterfall_5(med_episodes) + 
  scale_x_continuous(breaks = seq(0, 310, by = 50))
p_waterfall_5
p_waterfall_medications <- plot_waterfall_medication_pattern(med_episodes) + 
  scale_x_continuous(breaks = seq(0, 310, by = 50))
p_waterfall_medications
p_neo_resp      <- plot_neo_waterfall_with_response(neo_plans_bin)
p_neo_resp
p_neo_hist_ct   <- plot_neo_hist_counts(neo_hist)
p_neo_hist_ct
p_neo_hist_ct_side <- plot_neo_hist_counts_dodged(neo_hist, dose_range = 1:14)
p_neo_hist_ct_side
# -----------------------------------------------------------------------------
# 7) SAVE FIGURES
# -----------------------------------------------------------------------------
message("• Saving figures...")
save_png(tx_plot, "systemimc_tx_for_cscc", height = 7)
save_png(tx_nested_plot, "systemimc_tx_for_cscc_nested", height = 7)
save_png(p_ts_all, "time_series_plot", height = 7)
save_png(p_ts_both, "time_series_plot_both", height = 7)
save_png(p_waterfall_mono, "plans_B_monochrome", height = 7)
save_png(p_waterfall_binary, "plans_B_binary", height = 7)
save_png(p_waterfall_5,   "plans_C_detailed_intent", height = 7)
save_png(p_waterfall_medications,   "plans_C_detailed_intent_by_meds", height = 7)
save_png(p_neo_resp,      "plans_neoadjuvant_waterfall_binned", height = 7)
save_png(p_neo_hist_ct,   "plans_neoadjuvant_hist_counts", height = 7)
save_png(p_neo_hist_ct_side,   "plans_neoadjuvant_hist_side_by_side", height = 7)

# ---- Save figures (dual output) ---------------------------------------------

# tx_* plots: 12 × 8
out_tx_by_med         <- save_plot_dual(tx_plot,        "systemimc_tx_for_cscc",        height = 8)
out_tx_nested         <- save_plot_dual(tx_nested_plot, "systemimc_tx_for_cscc_nested", height = 8)
out_ts                <- save_plot_dual(p_ts_all,"time_series_plot",                    height = 8)
out_ts_both           <- save_plot_dual(p_ts_both,"time_series_plot_both",              height = 8)
out_ts_neo_strict     <- save_plot_dual(p_ts_neo_strict, "time_series_neo_only",        height = 8)
# plan-level waterfalls: 12 × 7
out_pA           <- save_plot_dual(p_waterfall_mono,       "plans_B_monochrome")
out_pB           <- save_plot_dual(p_waterfall_binary,     "plans_B_binary")
out_pC           <- save_plot_dual(p_waterfall_5,          "plans_C_detailed_intent")
out_pC_med       <- save_plot_dual(p_waterfall_medications,"plans_C_detailed_intent_by_meds")

# neoadjuvant-specific plots: 12 × 7
out_plan_neo_resp <- save_plot_dual(p_neo_resp,    "plans_neoadjuvant_waterfall_binned")
out_plan_neo_ct   <- save_plot_dual(p_neo_hist_ct, "plans_neoadjuvant_hist_counts")
out_plan_neo_side_by_side   <- save_plot_dual(p_neo_hist_ct_side, "plans_neoadjuvant_hist_side_by_side")

message("✓ 90_neo_execute.R complete.")
# =============================================================================
# End execute
# =============================================================================
