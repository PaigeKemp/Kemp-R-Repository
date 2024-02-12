#-------------------------------------------------------------------------------------------------------------
### 2004 Analyses - MPT
#-------------------------------------------------------------------------------------------------------------
# this script is to fit the processed data to an MPT model
# programmed by Vanessa Loaiza & Paige Kemp August 2023

# packages
pacman::p_load("tidyverse", "TreeBUGS", "ggplot2")

# global parameters
run_models <- FALSE # set to true if you want to run the models, or FALSE to read in the results
one_fill   <- "#F85A3E"
three_fill <- "#9D96B8"

# set to local working directory and then never again
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# load data -- keep only what is needed to keep environment tidy
load("../Data/2004_All_Data_df.RData") 
rm(exp1_combined_df_n, exp1_familiarity_df, exp1_test_data_df_belief_2, exp1_test_data_df_correction,
   exp2_combined_df_n, exp2_familiarity_df, exp2_test_data_df_belief_2, exp2_test_data_df_correction)

#-------------------------------------------------------------------------------------------------------------
### data wrangling

# wrangle the data -- join them up to keep it easy
exp1_test_data_df <- exp1_test_data_df %>% 
  select(subject, age_group, statementtype, recog_1_correct:recog_1_other) %>%
  # these next renames are going to make wrangling/modeling easier later
  rename(targ = recog_1_correct, lure = recog_1_misinfo, new = recog_1_other) %>%
  add_column(exp = 1)
exp2_test_data_df %>%
  # get rid of blank recall responses (they should not count as "new")
  mutate(testresponse = na_if(testresponse, "")) %>%
  drop_na(testresponse) %>%
  select(subject, age_group, statementtype, recall_1_correct:recall_1_other) %>%
  # these next renames are going to make wrangling/modeling easier later
  rename(targ = recall_1_correct, lure = recall_1_misinfo, new = recall_1_other) %>%
  add_column(exp = 2) %>% 
  # bind to E1 to make it easy for later
  bind_rows(exp1_test_data_df) %>% arrange(exp) -> exp_test_data
rm(exp1_test_data_df, exp2_test_data_df)

# some quick sanity checks
exp_test_data %>% group_by(exp, age_group, statementtype) %>% summarize_at(vars(targ:new), sum)

# carry on wrangling -- getting ready for analysis
exp_MPT_data <- exp_test_data %>%
  mutate(statementtype = factor(statementtype, labels = c("misinfo1", "misinfo3", "repeated"))) %>% 
  # repeated items don't make sense for the MPT model because the lure/other are equally unfamiliar
  filter(statementtype != "repeated") %>% 
  # aggregate frequencies for each participant
  group_by(exp, subject, age_group, statementtype) %>%
  summarize_all(sum) %>% ungroup %>%
  pivot_wider(names_from = statementtype, values_from = targ:new) %>%
  # there is one YA in E2 who had no observations in most of the cells, so we exclude 
  drop_na()

#-------------------------------------------------------------------------------------------------------------
### MPT model prep and fitting

# define the model equations directly within R
model_eqn <- 
  "# Fake News Independence Model
   misinfo1   targ_misinfo1    Pr_misinfo1
   misinfo1   targ_misinfo1    (1-Pr_misinfo1) * Pf_misinfo1 * Gc
   misinfo1   lure_misinfo1    (1-Pr_misinfo1) * Pf_misinfo1 * (1-Gc)
   misinfo1   targ_misinfo1    (1-Pr_misinfo1) * (1-Pf_misinfo1) * Gr * Gc
   misinfo1   lure_misinfo1    (1-Pr_misinfo1) * (1-Pf_misinfo1) * Gr * (1-Gc)
   misinfo1   new_misinfo1     (1-Pr_misinfo1) * (1-Pf_misinfo1) * (1-Gr)
   misinfo3   targ_misinfo3    Pr_misinfo3
   misinfo3   targ_misinfo3    (1-Pr_misinfo3) * Pf_misinfo3 * Gc
   misinfo3   lure_misinfo3    (1-Pr_misinfo3) * Pf_misinfo3 * (1-Gc)
   misinfo3   targ_misinfo3    (1-Pr_misinfo3) * (1-Pf_misinfo3) * Gr * Gc
   misinfo3   lure_misinfo3    (1-Pr_misinfo3) * (1-Pf_misinfo3) * Gr * (1-Gc)
   misinfo3   new_misinfo3     (1-Pr_misinfo3) * (1-Pf_misinfo3) * (1-Gr)"

# inspect the model
readEQN(model_eqn, paramOrder = TRUE)

# set restrictions
readEQN(model_eqn, restrictions = list("Gc = 0.5", "Gr = 0.33333"), paramOrder = TRUE)

# function to check the heterogeneity of each age group in each experiment
check_het <- function(d){
  print(testHetChi(freq = d, 
             tree = c("T_misinfo1", "T_misinfo1", "T_misinfo1",
                      "T_misinfo3", "T_misinfo3", "T_misinfo3")))
  plotFreq(x = d, eqn = model_eqn)
}

check_het(d = exp_MPT_data %>% 
            filter(exp == 1 & age_group == "younger") %>%
            select(targ_misinfo1:new_misinfo3))
check_het(d = exp_MPT_data %>% 
            filter(exp == 1 & age_group == "older") %>%
            select(targ_misinfo1:new_misinfo3))
check_het(d = exp_MPT_data %>% 
            filter(exp == 2 & age_group == "younger") %>%
            select(targ_misinfo1:new_misinfo3))
check_het(d = exp_MPT_data %>% 
            filter(exp == 2 & age_group == "older") %>%
            select(targ_misinfo1:new_misinfo3))

# function to fit the models for each age group in each experiment
fit_model <- function(data, whichModel){
  m <- traitMPT(eqnfile = model_eqn,
                data = data, 
                modelfilename = paste0("MPT/",whichModel,".jags"),
                parEstFile = paste0("MPT/",whichModel,"_results.txt"),
                restrictions = list("Gc = 0.5", "Gr = 0.33333"),
                transformedParameters = list("Pr = Pr_misinfo3 - Pr_misinfo1",
                                             "Pf = Pf_misinfo3 - Pf_misinfo1"),
                n.chain = 4, n.iter = 100000, n.adapt = 20000,
                n.burnin = 2000, n.thin = 5, ppp = 5000, dic = TRUE)
  saveRDS(m, file = paste0("MPT/",whichModel,".rds"))
  return(m)
}

# fit the models, setting a seed before each run for reproducible results
if (run_models){
  set.seed(123)
  m1y <- fit_model(data = exp_MPT_data %>%
                     filter(exp == 1 & age_group == "younger") %>%
                     select(targ_misinfo1:new_misinfo3),
                   whichModel = "E1_model_YA")
  set.seed(123)
  m1o <- fit_model(data = exp_MPT_data %>%
                     filter(exp == 1 & age_group == "older") %>%
                     select(targ_misinfo1:new_misinfo3),
                   whichModel = "E1_model_OA")
  set.seed(123)
  m2y <- fit_model(data = exp_MPT_data %>%
                     filter(exp == 2 & age_group == "younger") %>%
                     select(targ_misinfo1:new_misinfo3),
                   whichModel = "E2_model_YA")
  set.seed(123)
  m2o <- fit_model(data = exp_MPT_data %>%
                     filter(exp == 2 & age_group == "older") %>%
                     select(targ_misinfo1:new_misinfo3),
                   whichModel = "E2_model_OA")
} else {
  m1y <- readRDS("MPT/E1_model_YA.rds")
  m1o <- readRDS("MPT/E1_model_OA.rds")
  m2y <- readRDS("MPT/E2_model_YA.rds")
  m2o <- readRDS("MPT/E2_model_OA.rds")
}


#-------------------------------------------------------------------------------------------------------------
### MPT model output

# model summaries to check model convergence and fit
summary(m1y)
summary(m1o)
summary(m2y)
summary(m2o)

# plot to inspect convergence
plot(m1y)
plot(m1o)
plot(m2y)
plot(m2o)

# plot to inspect fit
plotFit(m1y)
plotFit(m1o)
plotFit(m2y)
plotFit(m2o)
plotFit(m1y, stat = "cov")
plotFit(m1o, stat = "cov")
plotFit(m2y, stat = "cov") # some misfit in the covariances
plotFit(m2o, stat = "cov") # some misfit in the covariances

# compare age groups within experiments for each parameter
m1_Pr_misinfo1 <- betweenSubjectMPT(m1y, m1o, par1 = "Pr_misinfo1")
m1_Pr_misinfo3 <- betweenSubjectMPT(m1y, m1o, par1 = "Pr_misinfo3")
m1_Pf_misinfo1 <- betweenSubjectMPT(m1y, m1o, par1 = "Pf_misinfo1")
m1_Pf_misinfo3 <- betweenSubjectMPT(m1y, m1o, par1 = "Pf_misinfo3")
m2_Pr_misinfo1 <- betweenSubjectMPT(m2y, m2o, par1 = "Pr_misinfo1")
m2_Pr_misinfo3 <- betweenSubjectMPT(m2y, m2o, par1 = "Pr_misinfo3")
m2_Pf_misinfo1 <- betweenSubjectMPT(m2y, m2o, par1 = "Pf_misinfo1")
m2_Pf_misinfo3 <- betweenSubjectMPT(m2y, m2o, par1 = "Pf_misinfo3")
m1_Pr_misinfo1
m1_Pr_misinfo3
m1_Pf_misinfo1
m1_Pf_misinfo3
m2_Pr_misinfo1
m2_Pr_misinfo3
m2_Pf_misinfo1
m2_Pf_misinfo3

# get individual estimates for the plots
m_means_indiv <- as_tibble(getParam(m1y, parameter = "theta", stat = "mean")) %>%
  rownames_to_column("ID") %>% add_column(exp = 1, age_group = "younger") %>% 
  relocate(exp:age_group, .before = "ID")
as_tibble(getParam(m1o, parameter = "theta", stat = "mean")) %>%
  rownames_to_column("ID") %>% add_column(exp = 1, age_group = "older") %>% 
  relocate(exp:age_group, .before = "ID") %>% 
  bind_rows(m_means_indiv) -> m_means_indiv
as_tibble(getParam(m2y, parameter = "theta", stat = "mean")) %>%
  rownames_to_column("ID") %>% add_column(exp = 2, age_group = "younger") %>% 
  relocate(exp:age_group, .before = "ID") %>% 
  bind_rows(m_means_indiv) -> m_means_indiv
as_tibble(getParam(m2o, parameter = "theta", stat = "mean")) %>%
  rownames_to_column("ID") %>% add_column(exp = 2, age_group = "older") %>% 
  relocate(exp:age_group, .before = "ID") %>% 
  bind_rows(m_means_indiv) %>% arrange(exp) -> m_means_indiv

m_means_indiv <- m_means_indiv %>%
  pivot_longer(!exp:ID, names_to = "parm", values_to = "value") %>%
  separate(col = parm, into = c("parm", "headline"), sep = "_") %>%
  mutate(parm = ifelse(parm == "Pr", "Recollection", "Familiarity"),
         parm = factor(parm, levels = c("Recollection", "Familiarity")),
         headline = factor(headline, labels = c("One Fake News Exposure",
                                                "Three Fake News Exposures")),
         age_group = factor(age_group, levels = c("younger", "older")),
         age_group = factor(age_group, labels = c("Younger", "Older")),
         exp = factor(exp, labels = c("Experiment 1", "Experiment 2")))
m_means_indiv

# get group mean estimates for the plots
m_means_group <- as_tibble(getParam(m1y, parameter = "mean", stat = "summary")) %>%
  add_column(exp = 1, age_group = "younger", 
             parm = c("Pf_misinfo1", "Pf_misinfo3", "Pr_misinfo1", "Pr_misinfo3")) %>%
  relocate(exp:parm, .before = "Mean")
as_tibble(getParam(m1o, parameter = "mean", stat = "summary")) %>%
  add_column(exp = 1, age_group = "older",
             parm = c("Pf_misinfo1", "Pf_misinfo3", "Pr_misinfo1", "Pr_misinfo3")) %>%
  relocate(exp:parm, .before = "Mean") %>% 
  bind_rows(m_means_group) -> m_means_group
as_tibble(getParam(m2y, parameter = "mean", stat = "summary")) %>%
  add_column(exp = 2, age_group = "younger",
             parm = c("Pf_misinfo1", "Pf_misinfo3", "Pr_misinfo1", "Pr_misinfo3")) %>%
  relocate(exp:parm, .before = "Mean") %>% 
  bind_rows(m_means_group) -> m_means_group
as_tibble(getParam(m2o, parameter = "mean", stat = "summary")) %>%
  add_column(exp = 2, age_group = "older",
             parm = c("Pf_misinfo1", "Pf_misinfo3", "Pr_misinfo1", "Pr_misinfo3")) %>%
  relocate(exp:parm, .before = "Mean") %>% 
  bind_rows(m_means_group) %>% arrange(exp) -> m_means_group

m_means_group <- m_means_group %>%
  separate(col = parm, into = c("parm", "headline"), sep = "_") %>%
  mutate(parm = ifelse(parm == "Pr", "Recollection", "Familiarity"),
         parm = factor(parm, levels = c("Recollection", "Familiarity")),
         headline = factor(headline, labels = c("One Fake News Exposure",
                                                "Three Fake News Exposures")),
         age_group = factor(age_group, levels = c("younger", "older")),
         age_group = factor(age_group, labels = c("Younger", "Older")),
         exp = factor(exp, labels = c("Experiment 1", "Experiment 2"))) %>%
  rename(lowerCI = `2.5%`, upperCI = `97.5%`, value = Mean)
m_means_group

m_plot <- ggplot(data = m_means_group, aes(x = age_group, y = value, fill = headline, color = headline)) +
  facet_grid(exp ~ parm) +
  geom_point(position = position_dodge(width = 0.5), stroke = 0.3, shape = 21, size = 4) +
  geom_point(data = m_means_indiv, aes(fill = headline), position = position_jitterdodge(dodge.width = 0.5), alpha = 0.08, size = 2) +
  geom_errorbar(aes(ymin = lowerCI, ymax = upperCI), position = position_dodge(width = 0.5), linewidth = 0.3, width = 0) +
  scale_fill_manual(values = c("One Fake News Exposure" = one_fill, "Three Fake News Exposures" = three_fill)) +
  scale_color_manual(values = c("One Fake News Exposure" = one_fill, "Three Fake News Exposures" = three_fill)) +
  scale_y_continuous(name = "Posterior Parameter Estimate\n", limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  xlab("\nAge Group") +
  theme(legend.position = c(.75,.58),
        legend.key.size = unit(.75, 'cm'),
        legend.title=element_blank(),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.margin = margin(0,0,0,0),
        strip.background=element_blank(),
        strip.text.x=element_text(size = 12, color = "black", margin = margin(0, 0, 8, 0)),
        strip.text.y=element_text(size = 12, color = "black", margin = margin(0, 0, 0, 8)),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_line(size = .3,color="black"),
        axis.text.x=element_text(size = 10,color="black"),
        axis.text.y=element_text(size = 10,color="black"),
        axis.title.x=element_text(size = 12, margin = margin(8, 0, 0, 0), hjust = .5),
        axis.title.y=element_text(size = 12, margin = margin(0, 8, 0, 0)),
        plot.margin = unit(c(.1, .1, .1, .1), "in"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.spacing = unit(-.0005, "in"),
        panel.border = element_rect(linewidth = .3, fill = NA, color = "black"))
m_plot
ggsave("MPT/plot.jpg", plot = m_plot, height = 14, width = 18, units = "cm", dpi = 300)

