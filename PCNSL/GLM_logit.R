
Sys.setenv(LANG = "en") # change language

library(survival)
library(survminer)


### (1)read file
load_df <- function(filename) {
  df <- read.csv(file = filename, row.names = 1)
  df <- as.data.frame(t(df))
  return(df)
}


### (2) select feature
feature_select <- function(df, feature, time = "PFS", label = "Label") {
  lm_df <- df[, c(time, feature, label)]

  lm_df[[label]] <- as.numeric(factor(lm_df[[label]],
                    labels = c("Good", "Poor"), ordered = TRUE))
  
  lm_df[[label]] <- lm_df[[label]] - 1
  lm_df[, c(time, feature)] <- as.data.frame(
                                lapply(lm_df[, c(time, feature)], as.numeric))

  return(lm_df)
}



###(3) Surival Anaylsis(Multi-feature)
surival_anaylsis <- function(df, feature, model, threshold = 0.06,
                             label = "Label", time = "PFS") {

  ###calculate Z-score
  glm1_coef <- model$coefficients
  glm1_x <- as.matrix(df[, feature])
  glm1_z <- glm1_x %*% glm1_coef

  ### Fit
  pfs <- df[[time]]
  sur_df <- data.frame(event = df[[label]], Zscore = glm1_z, PFS = pfs)
  sur_df[["scoreType"]] <- as.numeric(sur_df[["Zscore"]] > threshold)
  fit <- survfit(Surv(PFS, event) ~ scoreType,  # Create survival object
               data = sur_df) # dataset

  ###plot survival analysis curve
  mysurplot <-  ggsurvplot(fit, data = sur_df,
                surv.median.line = "hv",
                pval = TRUE,
                risk.table = TRUE,
                palette = "npg",
                conf.int = TRUE)

  return(mysurplot)
}

###(4) Surival Anaylsis(one-feature)
surival_anaylsis_single <- function(df, feature, model, threshold = 0.06,
                             label = "Label", time = "PFS") {

  ###calculate Z-score
  glm1_coef <- model$coefficients
  glm1_x <- as.matrix(df[, feature])
  glm1_z <- glm1_x %*% glm1_coef

   ### Fit
  pfs <- df[[time]]
  sur_df <- data.frame(event = df[[label]], Zscore = glm1_z, PFS = pfs)
  sur_df[["scoreType"]] <- as.numeric(sur_df[["Zscore"]] > threshold)
  fit <- survfit(Surv(PFS, event) ~ scoreType,  # Create survival object
               data = sur_df) # dataset

  ###plot survival analysis curve
  mysurplot <-  ggsurvplot(fit, data = sur_df,
                surv.median.line = "hv",
                pval = TRUE,
                risk.table = TRUE,
                palette = "npg",
                conf.int = TRUE)

  return(mysurplot)
}

### (5) total procedure
procedure <- function(feature) {
  #(1) read files & set parameters
  rc_df <- load_df("./data/batch/RC_rm_batch.csv")
  pc_df <- load_df("./data/batch/PC_rm_batch.csv")
  time <- c("PFS")
  label <- c("Label")

  ### set name of output files

  outfile <- "0"
  if (length(feature) > 1) {
    outfile <- "mixed_feature"
  }else {
    outfile <- feature
  }

  rc_lm_df <- feature_select(rc_df, feature, time = time, label = label)

  ###Logistics Regression
  formula1 <- "Label ~ . - PFS - 1"
  glm1 <- glm(formula = formula1, family = binomial("logit"), data = rc_lm_df)

  ###survival analysis in retrospective cohort
  rc_survivalplot <- surival_anaylsis(rc_lm_df, feature, glm1, threshold = 0.06,
                                      time = time, label = label)
  pdf(file = paste0("./result/rc_", outfile, ".pdf"))
  print(rc_survivalplot)
  dev.off()

  ### survival analysis in prospective cohort
  pc_lm_df <- feature_select(pc_df, feature)
  pc_survivalplot <- surival_anaylsis(pc_lm_df, feature, glm1, threshold = 0.06,
                                      time = time, label = label)
  pdf(file = paste0("./result/pc_", outfile, ".pdf"))
  print(pc_survivalplot)
  dev.off()
}


features <- c("4-Pyridoxic acid", "p-hydroxybenzoate", "cholesterol",
             "2-Isopropylmalic acid", "AMP", "dehydroascorbic acid")
procedure(features)

for (feature in features){
  procedure(feature)
}
