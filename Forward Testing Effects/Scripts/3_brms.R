# -----------------------------------------------------------------------------
# brms analyses of Kemp et al. aging + interpolated testing
# written by Vanessa Loaiza & Paige Kemp Feb 2023
# -----------------------------------------------------------------------------

# packages
packages <- c("tidyverse", "brms", "emmeans", "tidybayes")
purrr::map(packages, library, character.only = TRUE)

# directories
myDir <- dirname(rstudioapi::getSourceEditorContext()$path) 
setwd(myDir) # set to local directory only once in the script

# settings
which_exp <- 2     # flexibly read in, wrangle, and analyze either exp (1 or 2)
seed_repr <- 1     # set a seed for reproducible results
run_model <- TRUE # run the models (TRUE) or more quickly load results (FALSE)
save_data <- TRUE # save the output (TRUE) or skip if already saved (FALSE)
options(cli.num_colors=1) # this is to make printing the tibbles in the output work right

# -----------------------------------------------------------------------------

# read in data
data_inter <- read_csv(paste0('../Data/Exp',which_exp,'_Interpolated_data.csv'))
data_test  <- read_csv(paste0('../Data/Exp',which_exp,'_Test_data.csv'))

# checks
table(data_inter$ItemType, data_inter$Subject)
table(data_inter$ItemType, data_inter$AgeGroup)
table(data_test$ItemType, data_test$Subject)
table(data_test$ItemType, data_test$AgeGroup)

# bind the data, rejig a bit, get rid of stuff that we don't need
data_all <- left_join(data_test, data_inter, 
                      by = c("Subject", "AgeGroup", "ItemType", "Cue")) %>%
  select(-contains("blanks")) %>% # don't need this as far as I can tell
  rename(Target_L2 = Target.x, Target_L1 = Target.y) %>%
  mutate_at(vars(Subject, AgeGroup, ItemType, Cue), factor)
rm(data_inter, data_test)
glimpse(data_all)

# -----------------------------------------------------------------------------

# I hadn't used brms in a while and had trouble to get anything running
# seems like it was an R update issue as the advice in this link word for me FYI:
# https://stackoverflow.com/questions/72775716/cannot-get-install-of-brms-or-rstan-to-compile 

# general stuff for analysis -- priors following recs of Bartsch & Oberauer (2022)
pr <- c(prior(cauchy(0,.353), class = b), 
        prior(gamma(1,0.04), class = sd))
it <- 50000
ch <- 4
wu <- 1000

# check on our coding of the factors -- default is treatment coding
contrasts(data_all$ItemType)
contrasts(data_all$AgeGroup)

# function to make next analyses a little tidier since all the models use the
# same model setup 
fit_model <- function(formula, name, prior, data){
  m_fit <- brm(formula = formula,
               data = data,
               family = bernoulli(link = "logit"),
               prior = prior, iter = it, chains = ch, warmup = wu, cores = 4,
               control = list(max_treedepth = 15, adapt_delta = 0.99),
               save_pars = save_pars(all = TRUE), seed = seed_repr)
  print(paste0('../Output/Exp',which_exp,'_brms_',name,'.rds'))
  saveRDS(m_fit, file = paste0('../Output/Exp',which_exp,'_brms_',name,'.rds'))
  return(m_fit)
}

# ----------> Analysis 1: Interpolated test (list 1) performance

# formulas of models to test with predictors of interest
# eventually we will compare them to get BFs 
m1_int  <- bf(L1_Inter_Test_Acc ~ 1 + (1 | Subject) + (1 | Cue))
m1_age  <- bf(L1_Inter_Test_Acc ~ AgeGroup + (1 | Subject) + (1 | Cue))
# E1 has two AB_AD test conditions (2x2 design), E2 only has one (t-test situation), so adjust accordingly
if (which_exp == 1){
  m1_item <- bf(L1_Inter_Test_Acc ~ ItemType + (1 | Subject) + (1 | Cue))
  m1_both <- bf(L1_Inter_Test_Acc ~ AgeGroup + ItemType + (1 | Subject) + (1 | Cue))
  m1_full <- bf(L1_Inter_Test_Acc ~ AgeGroup * ItemType + (1 | Subject) + (1 | Cue))
  
  m1_formulas <- list(m1_int, m1_age, m1_item, m1_both, m1_full)
  m1_names    <- c("m1_int", "m1_age", "m1_item", "m1_both", "m1_full")
  m1_priors   <- list(NULL, pr, pr, pr, pr) # no prior for intercept-only model
} else if (which_exp == 2){
  m1_formulas <- list(m1_int, m1_age)
  m1_names    <- c("m1_int", "m1_age")
  m1_priors   <- list(NULL, pr) # no prior for intercept-only model
}

# fit the models
if (run_model){
  # pmap loops the fit_model function over the models with the same data
  # just keeps things tidy and avoids a lot of the copy/pastes 
  m1_fits <- pmap(list(m1_formulas, m1_names, m1_priors), fit_model,
                  # data: only tested L1 items (A_B_A_D_Test conditions)
                  data=data_all %>% filter(grepl("AB_AD_Test", ItemType)) %>% droplevels)
  # fit_model saves the individual models, but let's save them all here too
  saveRDS(m1_fits, file = paste0('../Output/Exp',which_exp,'_brms_m1_models.rds'))
} else {
  m1_fits <- readRDS(file = paste0('../Output/Exp',which_exp,'_brms_m1_models.rds'))
}

# # summarize results, inspect model fit and convergence visually
# # note that it takes a while to print the plots etc. so just skip if done already
# map(m1_fits, summary)
# map(m1_fits, pp_check)
# map(m1_fits, plot)

if (run_model){
  # compare best model (null) to the fixed effects models
  bf_m1_10 <- bayes_factor(m1_fits[[2]], m1_fits[[1]])
  if (which_exp == 1){
    bf_m1_20 <- bayes_factor(m1_fits[[3]], m1_fits[[1]])
    bf_m1_30 <- bayes_factor(m1_fits[[4]], m1_fits[[1]])
    bf_m1_40 <- bayes_factor(m1_fits[[5]], m1_fits[[1]])
    # save the data since those BFs take forever
    bf_m1 <- list(bf_m1_10, bf_m1_20, bf_m1_30, bf_m1_40)
  } else if (which_exp == 2) {
    bf_m1 <- bf_m1_10
  }
  
  saveRDS(bf_m1, file = paste0('../Output/Exp',which_exp,'_brms_m1_BFs.rds'))
} else {
  bf_m1 <- readRDS(file = paste0('../Output/Exp',which_exp,'_brms_m1_BFs.rds'))
}

# print the output
if (save_data){
  res_file <- paste0('../Output/Exp',which_exp,'_brms_m1_results.txt')
  cat("\n***************************************\n", file = res_file, append = TRUE)
  cat("\nAnalysis 1: Interpolated test (list 1) performance\n", file = res_file, append = TRUE)
  cat("\n***************************************\n", file = res_file, append = TRUE)
  cat("\nOutput from the models\n\n", file = res_file, append = TRUE)
  capture.output(m1_fits, file = res_file, append = TRUE)
  cat("\n***************************************\n", file = res_file, append = TRUE)
  cat("\nBayes factors comparing models\n\n", file = res_file, append = TRUE)
  cat("\n[[1]]=null, [[2]]=age, [[3]]=item, [[4]=both, [[5]]=full\n\n", file = res_file, append = TRUE)
  capture.output(bf_m1, file = res_file, append = TRUE)
}

# ----------> Analysis 2: Correct list 2 recall at final test

# formulas of models to test with predictors of interest
# eventually we will compare them to get BFs 
m2_int  <- bf(L2_Test_Acc ~ 1 + (1 | Subject) + (1 | Cue))
m2_age  <- bf(L2_Test_Acc ~ AgeGroup + (1 | Subject) + (1 | Cue))
m2_item <- bf(L2_Test_Acc ~ ItemType + (1 | Subject) + (1 | Cue))
m2_both <- bf(L2_Test_Acc ~ AgeGroup + ItemType + (1 | Subject) + (1 | Cue))
m2_full <- bf(L2_Test_Acc ~ AgeGroup * ItemType + (1 | Subject) + (1 | Cue))

m2_formulas <- list(m2_int, m2_age, m2_item, m2_both, m2_full)
m2_names    <- c("m2_int", "m2_age", "m2_item", "m2_both", "m2_full")
m2_priors   <- list(NULL, pr, pr, pr, pr) # no prior for intercept-only model

# fit the models
if (run_model){
  # pmap loops the fit_model function over the models with the same data
  # just keeps things tidy and avoids a lot of copy/pastes 
  m2_fits <- pmap(list(m2_formulas, m2_names, m2_priors), fit_model,
                  # data: all item types except for fillers (AB_AB items)
                  data=data_all %>% filter(ItemType != "AB_AB") %>% droplevels)
  # fit_model saves the individual models, but let's save them all here too
  saveRDS(m2_fits, file = paste0('../Output/Exp',which_exp,'_brms_m2_models.rds'))
} else {
  m2_fits <- readRDS(file = paste0('../Output/Exp',which_exp,'_brms_m2_models.rds'))
}

# # summarize results, inspect model fit and convergence visually
# # note that it takes a while to print the plots etc. so just skip if done already
# map(m2_fits, summary)
# map(m2_fits, pp_check)
# map(m2_fits, plot)

if (run_model){
  # compare each fixed effect model to the null model
  bf_m2_10 <- bayes_factor(m2_fits[[2]], m2_fits[[1]])
  bf_m2_20 <- bayes_factor(m2_fits[[3]], m2_fits[[1]])
  bf_m2_30 <- bayes_factor(m2_fits[[4]], m2_fits[[1]])
  bf_m2_40 <- bayes_factor(m2_fits[[5]], m2_fits[[1]])
  
  # save the data since those BFs take forever
  bf_m2 <- list(bf_m2_10, bf_m2_20, bf_m2_30, bf_m2_40)
  saveRDS(bf_m2, file = paste0('../Output/Exp',which_exp,'_brms_m2_BFs.rds'))
} else {
  bf_m2 <- readRDS(file = paste0('../Output/Exp',which_exp,'_brms_m2_BFs.rds'))
}

# pairwise comparisons -- all the simple effects
m2_pw <- m2_fits[[5]] %>%
  emmeans( ~ AgeGroup * ItemType) %>%
  contrast(method = "pairwise") %>%
  gather_emmeans_draws %>%
  tidybayes::mean_qi() %>% select(contrast, .value, .lower, .upper) %>%
  mutate(credible = ifelse(.lower < 0 & .upper < 0 | .lower > 0 & .upper > 0, 1, 0))
m2_pw

# age effect for each relevant item type contrast
m2_ageeff <- m2_fits[[5]] %>% 
  emmeans( ~ AgeGroup * ItemType) %>%
  contrast(interaction = "pairwise") %>%
  gather_emmeans_draws %>%
  tidybayes::mean_qi() %>% select(AgeGroup_pairwise, ItemType_pairwise, .value, .lower, .upper) %>%
  mutate(credible = ifelse(.lower < 0 & .upper < 0 | .lower > 0 & .upper > 0, 1, 0))
m2_ageeff

# print the output
if (save_data){
  res_file <- paste0('../Output/Exp',which_exp,'_brms_m2_results.txt')
  cat("\n***************************************\n", file = res_file, append = TRUE)
  cat("\nAnalysis 2: Correct list 2 recall\n", file = res_file, append = TRUE)
  cat("\n***************************************\n", file = res_file, append = TRUE)
  cat("\nOutput from the models\n\n", file = res_file, append = TRUE)
  capture.output(m2_fits, file = res_file, append = TRUE)
  cat("\n***************************************\n", file = res_file, append = TRUE)
  cat("\nBayes factors comparing models\n\n", file = res_file, append = TRUE)
  cat("\n[[1]]=null, [[2]]=age, [[3]]=item, [[4]=both, [[5]]=full\n\n", file = res_file, append = TRUE)
  capture.output(bf_m2, file = res_file, append = TRUE)
  cat("\n***************************************\n", file = res_file, append = TRUE)
  cat("\nPairwise comparisons\n\n", file = res_file, append = TRUE)
  capture.output(print(as_tibble(m2_pw), n=28), file = res_file, append = TRUE)
  cat("\n***************************************\n", file = res_file, append = TRUE)
  cat("\nage effect for each relevant item type contrast\n\n", file = res_file, append = TRUE)
  capture.output(print(as_tibble(m2_ageeff), n=28), file = res_file, append = TRUE)
}

# ----------> Analysis 3: Intrusions from list 1 during final test of list 2

# formulas of models to test with predictors of interest
# eventually we will compare them to get BFs 
m3_int  <- bf(L1_Intrusion ~ 1 + (1 | Subject) + (1 | Cue))
m3_age  <- bf(L1_Intrusion ~ AgeGroup + (1 | Subject) + (1 | Cue))
m3_item <- bf(L1_Intrusion ~ ItemType + (1 | Subject) + (1 | Cue))
m3_both <- bf(L1_Intrusion ~ AgeGroup + ItemType + (1 | Subject) + (1 | Cue))
m3_full <- bf(L1_Intrusion ~ AgeGroup * ItemType + (1 | Subject) + (1 | Cue))

m3_formulas <- list(m3_int, m3_age, m3_item, m3_both, m3_full)
m3_names    <- c("m3_int", "m3_age", "m3_item", "m3_both", "m3_full")
m3_priors   <- list(NULL, pr, pr, pr, pr) # no prior for intercept-only model

# fit the models
if (run_model){
  # pmap loops the fit_model function over the models with the same data
  # just keeps things tidy and avoids a lot of copy/pastes 
  m3_fits <- pmap(list(m3_formulas, m3_names, m3_priors), fit_model,
                  # data: all AB_AD item types (excluding AB_AB and AB_CD)
                  data=data_all %>% filter(grepl("AB_AD", ItemType)) %>% droplevels)
  # fit_model saves the individual models, but let's save them all here too
  saveRDS(m3_fits, file = paste0('../Output/Exp',which_exp,'_brms_m3_models.rds'))
} else {
  m3_fits <- readRDS(file = paste0('../Output/Exp',which_exp,'_brms_m3_models.rds'))
}

# # summarize results, inspect model fit and convergence visually
# # note that it takes a while to print the plots etc. so just skip if done already
# map(m3_fits, summary)
# map(m3_fits, pp_check)
# map(m3_fits, plot)

if (run_model){
  # compare each fixed effect model to the null model
  bf_m3_10 <- bayes_factor(m3_fits[[2]], m3_fits[[1]])
  bf_m3_20 <- bayes_factor(m3_fits[[3]], m3_fits[[1]])
  bf_m3_30 <- bayes_factor(m3_fits[[4]], m3_fits[[1]])
  bf_m3_40 <- bayes_factor(m3_fits[[5]], m3_fits[[1]])
  
  # save the data since those BFs take forever
  bf_m3 <- list(bf_m3_10, bf_m3_20, bf_m3_30, bf_m3_40) 
  saveRDS(bf_m3, file = paste0('../Output/Exp',which_exp,'_brms_m3_BFs.rds'))
} else {
  bf_m3 <- readRDS(file = paste0('../Output/Exp',which_exp,'_brms_m3_BFs.rds'))
}

# pairwise comparisons 
m3_pw <- m3_fits[[5]] %>%
  emmeans( ~ AgeGroup * ItemType) %>%
  contrast(method = "pairwise") %>%
  gather_emmeans_draws %>%
  tidybayes::mean_qi() %>% select(contrast, .value, .lower, .upper) %>%
  mutate(credible = ifelse(.lower < 0 & .upper < 0 | .lower > 0 & .upper > 0, 1, 0))
m3_pw

# age effect for each relevant item type contrast
m3_ageeff <- m3_fits[[5]] %>% 
  emmeans( ~ AgeGroup * ItemType) %>%
  contrast(interaction = "pairwise") %>%
  gather_emmeans_draws %>%
  tidybayes::mean_qi() %>% select(AgeGroup_pairwise, ItemType_pairwise, .value, .lower, .upper) %>%
  mutate(credible = ifelse(.lower < 0 & .upper < 0 | .lower > 0 & .upper > 0, 1, 0))
m3_ageeff

# print the output
if (save_data){
  res_file <- paste0('../Output/Exp',which_exp,'_brms_m3_results.txt')
  cat("\n***************************************\n", file = res_file, append = TRUE)
  cat("\nAnalysis 3: Intrusions from List 1\n", file = res_file, append = TRUE)
  cat("\n***************************************\n", file = res_file, append = TRUE)
  cat("\nOutput from the models\n\n", file = res_file, append = TRUE)
  capture.output(m3_fits, file = res_file, append = TRUE)
  cat("\n***************************************\n", file = res_file, append = TRUE)
  cat("\nBayes factors comparing models\n\n", file = res_file, append = TRUE)
  cat("\n[[1]]=null, [[2]]=age, [[3]]=item, [[4]=both, [[5]]=full\n\n", file = res_file, append = TRUE)
  capture.output(bf_m3, file = res_file, append = TRUE)
  cat("\n***************************************\n", file = res_file, append = TRUE)
  cat("\nPairwise comparisons\n\n", file = res_file, append = TRUE)
  capture.output(print(as_tibble(m3_pw), n=28), file = res_file, append = TRUE)
  cat("\n***************************************\n", file = res_file, append = TRUE)
  cat("\nage effect for each relevant item type contrast\n\n", file = res_file, append = TRUE)
  capture.output(print(as_tibble(m3_ageeff), n=28), file = res_file, append = TRUE)
}

# ----------> Analysis 4: Change recollection of previous list 1 items

# formulas of models to test with predictors of interest
# eventually we will compare them to get BFs 
m4_int  <- bf(ChangeRecollected ~ 1 + (1 | Subject) + (1 | Cue))
m4_age  <- bf(ChangeRecollected ~ AgeGroup + (1 | Subject) + (1 | Cue))
m4_item <- bf(ChangeRecollected ~ ItemType + (1 | Subject) + (1 | Cue))
m4_both <- bf(ChangeRecollected ~ AgeGroup + ItemType + (1 | Subject) + (1 | Cue))
m4_full <- bf(ChangeRecollected ~ AgeGroup * ItemType + (1 | Subject) + (1 | Cue))

m4_formulas <- list(m4_int, m4_age, m4_item, m4_both, m4_full)
m4_names    <- c("m4_int", "m4_age", "m4_item", "m4_both", "m4_full")
m4_priors   <- list(NULL, pr, pr, pr, pr) # no prior for intercept-only model

# fit the models
if (run_model){
  # pmap loops the fit_model function over the models with the same data
  # just keeps things tidy and avoids a lot of copy/pastes 
  m4_fits <- pmap(list(m4_formulas, m4_names, m4_priors), fit_model,
                  # data: all AB_AD item types (where there was a change)
                  data=data_all %>% filter(grepl("AB_AD", ItemType)) %>% droplevels)
  # fit_model saves the individual models, but let's save them all here too
  saveRDS(m4_fits, file = paste0('../Output/Exp',which_exp,'_brms_m4_models.rds'))
} else {
  m4_fits <- readRDS(file = paste0('../Output/Exp',which_exp,'_brms_m4_models.rds'))
}

# # summarize results, inspect model fit and convergence visually
# # note that it takes a while to print the plots etc. so just skip if done already
# map(m4_fits, summary)
# map(m4_fits, pp_check)
# map(m4_fits, plot)

if (run_model){
  # compare each fixed effect model to the null model
  bf_m4_10 <- bayes_factor(m4_fits[[2]], m4_fits[[1]])
  bf_m4_20 <- bayes_factor(m4_fits[[3]], m4_fits[[1]])
  bf_m4_30 <- bayes_factor(m4_fits[[4]], m4_fits[[1]])
  bf_m4_40 <- bayes_factor(m4_fits[[5]], m4_fits[[1]])
  
  # save the data since those BFs take forever
  bf_m4 <- list(bf_m4_10, bf_m4_20, bf_m4_30, bf_m4_40)
  saveRDS(bf_m4, file = paste0('../Output/Exp',which_exp,'_brms_m4_BFs.rds'))
} else {
  bf_m4 <- readRDS(file = paste0('../Output/Exp',which_exp,'_brms_m4_BFs.rds'))
}

# pairwise comparisons 
m4_pw <- m4_fits[[5]] %>%
  emmeans( ~ AgeGroup * ItemType) %>%
  contrast(method = "pairwise") %>%
  gather_emmeans_draws %>%
  tidybayes::mean_qi() %>% select(contrast, .value, .lower, .upper) %>%
  mutate(credible = ifelse(.lower < 0 & .upper < 0 | .lower > 0 & .upper > 0, 1, 0))
m4_pw

# age effect for each relevant item type contrast
m4_ageeff <- m4_fits[[5]] %>% 
  emmeans( ~ AgeGroup * ItemType) %>%
  contrast(interaction = "pairwise") %>%
  gather_emmeans_draws %>%
  tidybayes::mean_qi() %>% select(AgeGroup_pairwise, ItemType_pairwise, .value, .lower, .upper) %>%
  mutate(credible = ifelse(.lower < 0 & .upper < 0 | .lower > 0 & .upper > 0, 1, 0))
m4_ageeff

# print the output
if (save_data){
  res_file <- paste0('../Output/Exp',which_exp,'_brms_m4_results.txt')
  cat("\n***************************************\n", file = res_file, append = TRUE)
  cat("\nAnalysis 4: Change recollection\n", file = res_file, append = TRUE)
  cat("\n***************************************\n", file = res_file, append = TRUE)
  cat("\nOutput from the models\n\n", file = res_file, append = TRUE)
  capture.output(m4_fits, file = res_file, append = TRUE)
  cat("\n***************************************\n", file = res_file, append = TRUE)
  cat("\nBayes factors comparing models\n\n", file = res_file, append = TRUE)
  cat("\n[[1]]=null, [[2]]=age, [[3]]=item, [[4]=both, [[5]]=full\n\n", file = res_file, append = TRUE)
  capture.output(bf_m4, file = res_file, append = TRUE)
  cat("\n***************************************\n", file = res_file, append = TRUE)
  cat("\nPairwise comparisons\n\n", file = res_file, append = TRUE)
  capture.output(print(as_tibble(m4_pw), n=15), file = res_file, append = TRUE)
  cat("\n***************************************\n", file = res_file, append = TRUE)
  cat("\nage effect for each relevant item type contrast\n\n", file = res_file, append = TRUE)
  capture.output(print(as_tibble(m4_ageeff), n=28), file = res_file, append = TRUE)
}

# ----------> Analysis 5a: Correct list 2 recall cond. on change recollection

# formulas of models to test with predictors of interest
# eventually we will compare them to get BFs 
m5a_int  <- bf(L2_Test_Acc ~ 1 + (1 | Subject) + (1 | Cue))
m5a_age  <- bf(L2_Test_Acc ~ AgeGroup + (1 | Subject) + (1 | Cue))
m5a_item <- bf(L2_Test_Acc ~ ItemType + (1 | Subject) + (1 | Cue))
m5a_both <- bf(L2_Test_Acc ~ AgeGroup + ItemType + (1 | Subject) + (1 | Cue))
m5a_full <- bf(L2_Test_Acc ~ AgeGroup * ItemType + (1 | Subject) + (1 | Cue))

m5a_formulas <- list(m5a_int, m5a_age, m5a_item, m5a_both, m5a_full)
m5a_names    <- c("m5a_int", "m5a_age", "m5a_item", "m5a_both", "m5a_full")
m5a_priors   <- list(NULL, pr, pr, pr, pr) # no prior for intercept-only model

# fit the models
if (run_model){
  # pmap loops the fit_model function over the models with the same data
  # just keeps things tidy and avoids a lot of copy/pastes 
  m5a_fits <- pmap(list(m5a_formulas, m5a_names, m5a_priors), fit_model,
                  # data: changes that were recollected of AB_AD item types
                  data=data_all %>% filter(grepl("Change Recollected", Change_Class),
                                           grepl("AB_AD", ItemType)) %>% droplevels)
  # fit_model saves the individual models, but let's save them all here too
  saveRDS(m5a_fits, file = paste0('../Output/Exp',which_exp,'_brms_m5a_models.rds'))
} else {
  m5a_fits <- readRDS(file = paste0('../Output/Exp',which_exp,'_brms_m5a_models.rds'))
}

# # summarize results, inspect model fit and convergence visually
# # note that it takes a while to print the plots etc. so just skip if done already
# map(m5a_fits, summary)
# map(m5a_fits, pp_check)
# map(m5a_fits, plot)

if (run_model){
  # compare each fixed effect model to the null model
  bf_m5a_10 <- bayes_factor(m5a_fits[[2]], m5a_fits[[1]])
  bf_m5a_20 <- bayes_factor(m5a_fits[[3]], m5a_fits[[1]])
  bf_m5a_30 <- bayes_factor(m5a_fits[[4]], m5a_fits[[1]])
  bf_m5a_40 <- bayes_factor(m5a_fits[[5]], m5a_fits[[1]])
  
  # save the data since those BFs take forever
  bf_m5a <- list(bf_m5a_10, bf_m5a_20, bf_m5a_30, bf_m5a_40)
  saveRDS(bf_m5a, file = paste0('../Output/Exp',which_exp,'_brms_m5a_BFs.rds'))
} else {
  bf_m5a <- readRDS(file = paste0('../Output/Exp',which_exp,'_brms_m5a_BFs.rds'))
}

# pairwise comparisons
m5a_pw <- m5a_fits[[5]] %>%
  emmeans( ~ AgeGroup * ItemType) %>%
  contrast(method = "pairwise") %>%
  gather_emmeans_draws %>%
  tidybayes::mean_qi() %>% select(contrast, .value, .lower, .upper) %>%
  mutate(credible = ifelse(.lower < 0 & .upper < 0 | .lower > 0 & .upper > 0, 1, 0))
m5a_pw

# age effect for each relevant item type contrast
m5a_ageeff <- m5a_fits[[5]] %>% 
  emmeans( ~ AgeGroup * ItemType) %>%
  contrast(interaction = "pairwise") %>%
  gather_emmeans_draws %>%
  tidybayes::mean_qi() %>% select(AgeGroup_pairwise, ItemType_pairwise, .value, .lower, .upper) %>%
  mutate(credible = ifelse(.lower < 0 & .upper < 0 | .lower > 0 & .upper > 0, 1, 0))
m5a_ageeff

# print the output
if (save_data){
  res_file <- paste0('../Output/Exp',which_exp,'_brms_m5a_results.txt')
  cat("\n***************************************\n", file = res_file, append = TRUE)
  cat("\nAnalysis 5a: List 2 recall conditionalized on list 1 change recollection\n", file = res_file, append = TRUE)
  cat("\n***************************************\n", file = res_file, append = TRUE)
  cat("\nOutput from the models\n\n", file = res_file, append = TRUE)
  capture.output(m5a_fits, file = res_file, append = TRUE)
  cat("\n***************************************\n", file = res_file, append = TRUE)
  cat("\nBayes factors comparing models\n\n", file = res_file, append = TRUE)
  cat("\n[[1]]=null, [[2]]=age, [[3]]=item, [[4]=both, [[5]]=full\n\n", file = res_file, append = TRUE)
  capture.output(bf_m5a, file = res_file, append = TRUE)
  cat("\n***************************************\n", file = res_file, append = TRUE)
  cat("\nPairwise comparisons\n\n", file = res_file, append = TRUE)
  capture.output(print(as_tibble(m5a_pw), n=15), file = res_file, append = TRUE)
  cat("\n***************************************\n", file = res_file, append = TRUE)
  cat("\nage effect for each relevant item type contrast\n\n", file = res_file, append = TRUE)
  capture.output(print(as_tibble(m5a_ageeff), n=28), file = res_file, append = TRUE)
}

# ----------> Analysis 5b: Correct list 2 recall cond. on change remembered

# formulas of models to test with predictors of interest
# eventually we will compare them to get BFs 
m5b_int  <- bf(L2_Test_Acc ~ 1 + (1 | Subject) + (1 | Cue))
m5b_age  <- bf(L2_Test_Acc ~ AgeGroup + (1 | Subject) + (1 | Cue))
m5b_item <- bf(L2_Test_Acc ~ ItemType + (1 | Subject) + (1 | Cue))
m5b_both <- bf(L2_Test_Acc ~ AgeGroup + ItemType + (1 | Subject) + (1 | Cue))
m5b_full <- bf(L2_Test_Acc ~ AgeGroup * ItemType + (1 | Subject) + (1 | Cue))

m5b_formulas <- list(m5b_int, m5b_age, m5b_item, m5b_both, m5b_full)
m5b_names    <- c("m5b_int", "m5b_age", "m5b_item", "m5b_both", "m5b_full")
m5b_priors   <- list(NULL, pr, pr, pr, pr) # no prior for intercept-only model

# fit the models
if (run_model){
  # pmap loops the fit_model function over the models with the same data
  # just keeps things tidy and avoids a lot of copy/pastes 
  m5b_fits <- pmap(list(m5b_formulas, m5b_names, m5b_priors), fit_model,
                   # data: changes that were remembered of AB_AD item types
                   data=data_all %>% filter(grepl("Change Remembered", Change_Class),
                                            grepl("AB_AD", ItemType)) %>% droplevels)
  # fit_model saves the individual models, but let's save them all here too
  saveRDS(m5b_fits, file = paste0('../Output/Exp',which_exp,'_brms_m5b_models.rds'))
} else {
  m5b_fits <- readRDS(file = paste0('../Output/Exp',which_exp,'_brms_m5b_models.rds'))
}

# # summarize results, inspect model fit and convergence visually
# # note that it takes a while to print the plots etc. so just skip if done already
# map(m5b_fits, summary)
# map(m5b_fits, pp_check)
# map(m5b_fits, plot)

if (run_model){
  # compare each fixed effect model to the null model
  bf_m5b_10 <- bayes_factor(m5b_fits[[2]], m5b_fits[[1]])
  bf_m5b_20 <- bayes_factor(m5b_fits[[3]], m5b_fits[[1]])
  bf_m5b_30 <- bayes_factor(m5b_fits[[4]], m5b_fits[[1]])
  bf_m5b_40 <- bayes_factor(m5b_fits[[5]], m5b_fits[[1]])
  
  # save the data since those BFs take forever
  bf_m5b <- list(bf_m5b_10, bf_m5b_20, bf_m5b_30, bf_m5b_40)
  saveRDS(bf_m5b, file = paste0('../Output/Exp',which_exp,'_brms_m5b_BFs.rds'))
} else {
  bf_m5b <- readRDS(file = paste0('../Output/Exp',which_exp,'_brms_m5b_BFs.rds'))
}

# pairwise comparisons
m5b_pw <- m5b_fits[[5]] %>%
  emmeans( ~ AgeGroup * ItemType) %>%
  contrast(method = "pairwise") %>%
  gather_emmeans_draws %>%
  tidybayes::mean_qi() %>% select(contrast, .value, .lower, .upper) %>%
  mutate(credible = ifelse(.lower < 0 & .upper < 0 | .lower > 0 & .upper > 0, 1, 0))
m5b_pw

# age effect for each relevant item type contrast
m5b_ageeff <- m5b_fits[[5]] %>% 
  emmeans( ~ AgeGroup * ItemType) %>%
  contrast(interaction = "pairwise") %>%
  gather_emmeans_draws %>%
  tidybayes::mean_qi() %>% select(AgeGroup_pairwise, ItemType_pairwise, .value, .lower, .upper) %>%
  mutate(credible = ifelse(.lower < 0 & .upper < 0 | .lower > 0 & .upper > 0, 1, 0))
m5b_ageeff

# print the output
if (save_data){
  res_file <- paste0('../Output/Exp',which_exp,'_brms_m5b_results.txt')
  cat("\n***************************************\n", file = res_file, append = TRUE)
  cat("\nAnalysis 5b: List 2 recall conditionalized on list 1 change recollection\n", file = res_file, append = TRUE)
  cat("\n***************************************\n", file = res_file, append = TRUE)
  cat("\nOutput from the models\n\n", file = res_file, append = TRUE)
  capture.output(m5b_fits, file = res_file, append = TRUE)
  cat("\n***************************************\n", file = res_file, append = TRUE)
  cat("\nBayes factors comparing models\n\n", file = res_file, append = TRUE)
  cat("\n[[1]]=null, [[2]]=age, [[3]]=item, [[4]=both, [[5]]=full\n\n", file = res_file, append = TRUE)
  capture.output(bf_m5b, file = res_file, append = TRUE)
  cat("\n***************************************\n", file = res_file, append = TRUE)
  cat("\nPairwise comparisons\n\n", file = res_file, append = TRUE)
  capture.output(print(as_tibble(m5b_pw), n=15), file = res_file, append = TRUE)
  cat("\n***************************************\n", file = res_file, append = TRUE)
  cat("\nage effect for each relevant item type contrast\n\n", file = res_file, append = TRUE)
  capture.output(print(as_tibble(m5b_ageeff), n=28), file = res_file, append = TRUE)
}

# ----------> Analysis 5c: Correct list 2 recall cond. on change not remembered

# formulas of models to test with predictors of interest
# eventually we will compare them to get BFs 
m5c_int  <- bf(L2_Test_Acc ~ 1 + (1 | Subject) + (1 | Cue))
m5c_age  <- bf(L2_Test_Acc ~ AgeGroup + (1 | Subject) + (1 | Cue))
m5c_item <- bf(L2_Test_Acc ~ ItemType + (1 | Subject) + (1 | Cue))
m5c_both <- bf(L2_Test_Acc ~ AgeGroup + ItemType + (1 | Subject) + (1 | Cue))
m5c_full <- bf(L2_Test_Acc ~ AgeGroup * ItemType + (1 | Subject) + (1 | Cue))

m5c_formulas <- list(m5c_int, m5c_age, m5c_item, m5c_both, m5c_full)
m5c_names    <- c("m5c_int", "m5c_age", "m5c_item", "m5c_both", "m5c_full")
m5c_priors   <- list(NULL, pr, pr, pr, pr) # no prior for intercept-only model

# fit the models
if (run_model){
  # pmap loops the fit_model function over the models with the same data
  # just keeps things tidy and avoids a lot of copy/pastes 
  m5c_fits <- pmap(list(m5c_formulas, m5c_names, m5c_priors), fit_model,
                   # data: changes that were not remembered of AB_AD item types
                   data=data_all %>% filter(grepl("Change Not Remembered", Change_Class),
                                            grepl("AB_AD", ItemType)) %>% droplevels)
  # fit_model saves the individual models, but let's save them all here too
  saveRDS(m5c_fits, file = paste0('../Output/Exp',which_exp,'_brms_m5c_models.rds'))
} else {
  m5c_fits <- readRDS(file = paste0('../Output/Exp',which_exp,'_brms_m5c_models.rds'))
}

# # summarize results, inspect model fit and convergence visually
# # note that it takes a while to print the plots etc. so just skip if done already
# map(m5c_fits, summary)
# map(m5c_fits, pp_check)
# map(m5c_fits, plot)

if (run_model){
  # compare each fixed effect model to the null model
  bf_m5c_10 <- bayes_factor(m5c_fits[[2]], m5c_fits[[1]])
  bf_m5c_20 <- bayes_factor(m5c_fits[[3]], m5c_fits[[1]])
  bf_m5c_30 <- bayes_factor(m5c_fits[[4]], m5c_fits[[1]])
  bf_m5c_40 <- bayes_factor(m5c_fits[[5]], m5c_fits[[1]])
  
  # compare the best model to the next best model 
  bf_m5c_23 <- bayes_factor(m5c_fits[[3]], m5c_fits[[4]])
  bf_m5c_24 <- bayes_factor(m5c_fits[[3]], m5c_fits[[5]])
  
  # save the data since those BFs take forever
  bf_m5c <- list(bf_m5c_10, bf_m5c_20, bf_m5c_30, bf_m5c_40, bf_m5c_23, bf_m5c_24)
  saveRDS(bf_m5c, file = paste0('../Output/Exp',which_exp,'_brms_m5c_BFs.rds'))
} else {
  bf_m5c <- readRDS(file = paste0('../Output/Exp',which_exp,'_brms_m5c_BFs.rds'))
}

# pairwise comparisons
m5c_pw <- m5c_fits[[5]] %>%
  emmeans( ~ AgeGroup * ItemType) %>%
  contrast(method = "pairwise") %>%
  gather_emmeans_draws %>%
  tidybayes::mean_qi() %>% select(contrast, .value, .lower, .upper) %>%
  mutate(credible = ifelse(.lower < 0 & .upper < 0 | .lower > 0 & .upper > 0, 1, 0))
m5c_pw

# age effect for each relevant item type contrast
m5c_ageeff <- m5c_fits[[5]] %>% 
  emmeans( ~ AgeGroup * ItemType) %>%
  contrast(interaction = "pairwise") %>%
  gather_emmeans_draws %>%
  tidybayes::mean_qi() %>% select(AgeGroup_pairwise, ItemType_pairwise, .value, .lower, .upper) %>%
  mutate(credible = ifelse(.lower < 0 & .upper < 0 | .lower > 0 & .upper > 0, 1, 0))
m5c_ageeff

# print the output
if (save_data){
  res_file <- paste0('../Output/Exp',which_exp,'_brms_m5c_results.txt')
  cat("\n***************************************\n", file = res_file, append = TRUE)
  cat("\nAnalysis 5c: List 2 recall conditionalized on list 1 change recollection\n", file = res_file, append = TRUE)
  cat("\n***************************************\n", file = res_file, append = TRUE)
  cat("\nOutput from the models\n\n", file = res_file, append = TRUE)
  capture.output(m5c_fits, file = res_file, append = TRUE)
  cat("\n***************************************\n", file = res_file, append = TRUE)
  cat("\nBayes factors comparing models\n\n", file = res_file, append = TRUE)
  cat("\n[[1]]=null, [[2]]=age, [[3]]=item, [[4]=both, [[5]]=full\n\n", file = res_file, append = TRUE)
  capture.output(bf_m5c, file = res_file, append = TRUE)
  cat("\n***************************************\n", file = res_file, append = TRUE)
  cat("\nPairwise comparisons\n\n", file = res_file, append = TRUE)
  capture.output(print(as_tibble(m5c_pw), n=15), file = res_file, append = TRUE)
  cat("\n***************************************\n", file = res_file, append = TRUE)
  cat("\nage effect for each relevant item type contrast\n\n", file = res_file, append = TRUE)
  capture.output(print(as_tibble(m5c_ageeff), n=28), file = res_file, append = TRUE)
}

