library(tidyverse)
library(survival)
library(adjrct)
### 03-Using-R.Rmd #############################################################
data_url <-
  paste0("https://github.com/jbetz-jhu/CovariateAdjustmentTutorial",
         "/raw/main/Simulated_MISTIE_III_v1.2.csv")

sim_miii_full <- read.csv(file = url(data_url))

# Read in data: Recast categorical variables as factors
sim_miii_full <-
  sim_miii_full %>%
  dplyr::tibble() %>%
  dplyr::mutate(
    # Convert variables from binary indicators to labeled categorical variables
    male =
      factor(
        x = male,
        levels = 0:1,
        labels = c("0. Female", "1. Male")
      ),
    across(
      .cols =
        all_of(
          x = c("hx_cvd", "hx_hyperlipidemia",
                "on_anticoagulants", "on_antiplatelets")
        ),
      .fns = function(x) factor(x, levels = 0:1, labels = c("0. No", "1. Yes"))
    ),
    # Convert GCS and MRS variables from character data to categorical variables
    across(
      .cols = starts_with("gcs") | starts_with("mrs"),
      .fns = factor
    ),
    ich_location =
      factor(
        x = ich_location,
        levels = c("Deep", "Lobar")
      ),
    arm =
      factor(
        x = arm,
        levels = c("medical", "surgical")
      ),
    tx = 1*(arm == "surgical")
  )


# Take the first 500 rows
sim_miii <-
  sim_miii_full %>%
  dplyr::slice(1:500)


#-------------------------------------------------------------------------------


wilcox_to_auc <-
  function(data, indices = NULL, formula){
    # Input data must be a data.frame
    if(!all(class(data) == "data.frame")){
      stop("`data` must be a data.frame: use `as.data.frame()` for a tibble.")
    }

    # If bootstrap indices not supplied, use entire dataset
    if(is.null(indices)) indices <- 1:nrow(data)

    # Extract Outcome/Treatment from Formula
    outcome <- all.vars(update(formula, . ~ 0))
    treatment <- all.vars(update(formula, 0 ~ .))
    stopifnot(length(treatment) == 1)

    # Convert outcome to numeric using levels: Assumes levels are ordered
    if(!is.numeric(data[, outcome])){
      data[, outcome] <- as.numeric(data[, outcome])
    }

    # Run Wilcoxon Rank Sum on data using the bootstrap indices
    wrst_result <-
      wilcox.test(
        formula = formula,
        data = data[indices,]
      )

    # Compute AUC statistic
    return(wrst_result$statistic/prod(table(data[indices, treatment])))
  }

n_boot_samples <- 10000

mrs_365d_auc_boot <-
  boot::boot(
    data = as.data.frame(sim_miii),
    statistic = wilcox_to_auc,
    R = n_boot_samples,
    formula = mrs_365d ~ arm
  )


#-------------------------------------------------------------------------------


mrs_365d_auc_boot_ci <-
  boot::boot.ci(
    boot.out = mrs_365d_auc_boot,
    conf = 0.95,
    type = "bca",
    index = 1
  )




### 04-Standardization.Rmd #####################################################
data_url <-
  paste0("https://github.com/jbetz-jhu/CovariateAdjustmentTutorial/",
         "raw/main/SIMULATED_CTN03_220506.Rdata")

load(file = url(data_url))


#-------------------------------------------------------------------------------


# Write a function to produce the ANCOVA estimate
margins_fun <-
  function(data, indices = NULL, formula, family, term){
    # Input data must be a data.frame
    if(!all(class(data) == "data.frame")){
      stop("`data` must be a data.frame: use `as.data.frame()` for a tibble.")
    }

    # If bootstrap indices not supplied, use entire dataset
    if(is.null(indices)) indices <- 1:nrow(data)

    data <- data[indices,]

    glm_fit <-
      glm(
        formula = formula,
        family = family,
        data = data
      )

    tx_levels <- levels(data[, term])

    e_y_1 <-
      predict(
        object = glm_fit,
        newdata =
          within(
            data,
            expr = assign(x = term, value = tx_levels[2])
          ),

        type = "response"
      )

    e_y_0 <-
      predict(
        object = glm_fit,
        newdata =
          within(
            data,
            expr = assign(x = term, value = tx_levels[1])
          ),

        type = "response"
      )

    return(mean(e_y_1) - mean(e_y_0))
  }

vas_ancova_boot <-
  boot::boot(
    data = ctn03_sim_mar,
    statistic = margins_fun,
    R = 10000,
    formula = vas_crave_opiates_eot ~ arm + vas_crave_opiates_bl,
    family = gaussian(link = "identity"),
    term = "arm"
  )




### 06-Time-To-Event-Outcomes.Rmd ##############################################
surv_metadata_unadj <-
  adjrct::survrct(
    outcome.formula =
      Surv(days_on_study, died_on_study) ~ tx,
    trt.formula = tx ~ 1,
    data = sim_miii,
  )

surv_prob_unadj <-
  adjrct::survprob(
    metadata = surv_metadata_unadj,
    horizon = 90
  )

#-------------------------------------------------------------------------------
rmst_unadj <-
  adjrct::rmst(
    metadata = surv_metadata_unadj,
    horizon = 90
  )
#-------------------------------------------------------------------------------
surv_metadata_adj <-
  adjrct::survrct(
    outcome.formula =
      Surv(days_on_study, died_on_study) ~
      tx + age + male + hx_cvd + hx_hyperlipidemia +
      on_anticoagulants + on_antiplatelets + ich_location +
      ich_s_volume + ivh_s_volume + gcs_category,
    trt.formula =
      tx ~
      age + male + hx_cvd + hx_hyperlipidemia +
      on_anticoagulants + on_antiplatelets + ich_location +
      ich_s_volume + ivh_s_volume + gcs_category,
    data = sim_miii
  )
#-------------------------------------------------------------------------------
surv_prob_adj <-
  adjrct::survprob(
    metadata = surv_metadata_adj,
    horizon = 90
  )
#-------------------------------------------------------------------------------
rmst_adj <-
  adjrct::rmst(
    metadata = surv_metadata_adj,
    horizon = 90
  )

save(
  list =
    c("sim_miii_full",
      "ctn03_sim_mar",
      "mrs_365d_auc_boot",
      "mrs_365d_auc_boot_ci",
      "vas_ancova_boot",
      "surv_metadata_unadj",
      "surv_prob_unadj",
      "rmst_unadj",
      "surv_metadata_adj",
      "surv_prob_adj",
      "rmst_adj"),
  file = cached_results_path
)
