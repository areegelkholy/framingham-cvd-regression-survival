library("RwR") 
library("knitr")
library("markdown")
library("ggplot2")
library("gridExtra")
library("broom")
library("splines")
library("skimr")
library("tibble")
library(GGally)
library(readr)
library(dplyr)
library(hexbin)
library(riskCommunicator)
library(tidyr)
library(effects)

#Importing data
data("framingham", package = "riskCommunicator")
skim(framingham)
?framingham

#Making the relevant columns into factors
factor_columns <- c("SEX","CURSMOKE","DIABETES","BPMEDS","PREVCHD","PREVAP","PREVMI",
                    "PREVSTRK","PREVHYP","DEATH","ANGINA","HOSPMI","MI_FCHD","ANYCHD",
                    "STROKE","CVD","HYPERTEN","educ","PERIOD")
framingham[factor_columns] <- lapply(framingham[factor_columns], factor)
skim(framingham)

#Removing columns corresponding to responses we do not care about / do
#not have information about at time 0.
#This includes: DEATH, ANGINA, HOSPMI, MI_FCHD, ANYCHD, STROKE, HYPERTEN, 
#TIMEAP, TIMEMI, TIMEMIFC, TIMECHD, TIMESTRK, TIMEDTH, TIMEHYP.
drop <- c("DEATH","ANGINA","HOSPMI","MI_FCHD","ANYCHD","STROKE","HYPERTEN",
          "TIMEAP","TIMEMI","TIMEMIFC","TIMECHD","TIMESTRK","TIMEDTH","TIMEHYP","PREVMI")

framingham <- framingham[ , !(names(framingham) %in% drop)]
skim(framingham)

idx <- framingham$PERIOD %in% c(1, 2)

framingham$HDLC[idx][is.na(framingham$HDLC[idx])] <- 0
framingham$LDLC[idx][is.na(framingham$LDLC[idx])] <- 0

#Removing missing Values
framingham_clean <- na.omit(framingham)
skim(framingham_clean)

#Spearman correlation
framingham_clean_Spearman <- framingham_clean
framingham_clean_Spearman$HDLC <- NULL
framingham_clean_Spearman$LDLC <- NULL
framingham_clean_Spearman$RANDID <- NULL
#Code taken from Chap 3.
cp_full <- cor(
  data.matrix(framingham_clean_Spearman), 
  use = "complete.obs", 
  method = "spearman"
)

corrplot::corrplot(
  cp_full, 
  diag = FALSE, 
  order = "hclust", 
  addrect = 4, 
  tl.srt = 45, 
  tl.col = "black", 
  tl.cex = 0.8
)

#-------------------------------------------------------------------------------
#Regression models
framingham_clean_period1 <- framingham_clean[framingham_clean$PERIOD == 1, ]
form_period1 <- CVD ~ PREVHYP + PREVAP + SEX + DIABP + DIABETES +
        BPMEDS + educ + AGE + BMI + GLUCOSE + CIGPDAY + TOTCHOL
model_period1 <- glm(form_period1, data = framingham_clean_period1, family = binomial())
summary(model_period1)

#Lets look at some residuals, code from chap 10.
model1_aug <- augment(model_period1, type.residuals = "pearson") |>
  mutate(.fitted_group = cut(.fitted, quantile(.fitted, probs = seq(5, 95, 1) / 100)))


p1 <- model1_aug |>
  ggplot(aes(.fitted, .resid)) + 
  geom_point() + 
  geom_smooth() +
  xlab("fitted values") + 
  ylab("Pearson residuals")

p4 <- model1_aug |>
  group_by(.fitted_group) |>
  summarize(.fitted_local = mean(.fitted), .var_local = var(.resid)) |>
  na.omit() |> 
  ggplot(aes(.fitted_local, .var_local)) + 
  geom_point() + 
  geom_hline(yintercept = 1, linetype = 2) +
  geom_smooth() + 
  xlab("fitted values") +
  ylab("variance")

gridExtra::grid.arrange(p1, p4, ncol = 2)

#From Chap 3:
drop1(model_period1, test = "Chisq") |> 
  tidy() |> 
  arrange(p.value) |>
  filter(term != "<none>") |>
  select(-AIC) |>
  knitr::kable()
#We see that GLUCOSE, educ and BPMEDS should not be included.

#Lets determine if this model is significantly better:
form_period1_new <- CVD ~ AGE + SEX + TOTCHOL + PREVHYP + PREVAP + DIABP + DIABETES + CIGPDAY + BMI
model_period1_new <- glm(form_period1_new, data = framingham_clean_period1, family = binomial())

anova(model_period1, model_period1_new, test = "Chisq")
#The anova test shows that the original model is more accurate. However this
#improvement is fairly small compared to the total deviance, only around 0.29%.

#Nonlinear Considerations:
model_period1_new_diag <- augment(model_period1_new)
p1 <- ggplot(model_period1_new_diag, aes(AGE, .std.resid)) +
  stat_binhex(bins = 20) + geom_smooth(linewidth =  1, fill = "blue") +
  xlab("AGE") + ylab("")
p2 <- ggplot(model_period1_new_diag, aes(SEX, .std.resid)) +
  stat_binhex(bins = 20) + geom_smooth(linewidth =  1, fill = "blue") +
  xlab("SEX") + ylab("")
p3 <- ggplot(model_period1_new_diag, aes(TOTCHOL, .std.resid)) +
  stat_binhex(bins = 20) + geom_smooth(linewidth =  1, fill = "blue") +
  xlab("TOTCHOL") + ylab("")
p4 <- ggplot(model_period1_new_diag, aes(PREVHYP, .std.resid)) +
  stat_binhex(bins = 20) + geom_smooth(linewidth =  1, fill = "blue") +
  xlab("PREVHYP") + ylab("")
p5 <- ggplot(model_period1_new_diag, aes(PREVAP, .std.resid)) +
  stat_binhex(bins = 20) + geom_smooth(linewidth =  1, fill = "blue") +
  xlab("PREVAP") + ylab("")
p6 <- ggplot(model_period1_new_diag, aes(DIABP, .std.resid)) +
  stat_binhex(bins = 20) + geom_smooth(linewidth =  1, fill = "blue") +
  xlab("DIABP") + ylab("")
p7 <- ggplot(model_period1_new_diag, aes(DIABETES, .std.resid)) +
  stat_binhex(bins = 20) + geom_smooth(linewidth =  1, fill = "blue") +
  xlab("DIABETES") + ylab("")
p8 <- ggplot(model_period1_new_diag, aes(CIGPDAY, .std.resid)) +
  stat_binhex(bins = 20) + geom_smooth(linewidth =  1, fill = "blue") +
  xlab("CIGPDAY") + ylab("")
p9 <- ggplot(model_period1_new_diag, aes(BMI, .std.resid)) +
  stat_binhex(bins = 20) + geom_smooth(linewidth =  1, fill = "blue") +
  xlab("BMI") + ylab("")
gridExtra::grid.arrange(p1, p3, p6, p8, p9, ncol = 2)
#We see that there is no significant nonlinear contributions. Maybe note that 
#this is different from when we attempted to include Period 2 and 3.

#Lets look at how much interactions would improve the model:
form_period1_new_int <- CVD ~ (AGE + SEX + TOTCHOL + PREVHYP + PREVAP + DIABP + DIABETES + CIGPDAY + BMI)^2
model_period1_new_int <- glm(form_period1_new_int, data = framingham_clean_period1, family = binomial())

drop1(model_period1_new_int, test = "Chisq") |> 
  tidy() |> 
  arrange(p.value) |>
  filter(term != "<none>") |>
  select(-AIC) |>
  knitr::kable()

#Graphs of interactions if needed.
plot(Effect(c("SEX", "CIGPDAY"), model_period1_new))

anova(model_period1_new, model_period1_new_int, test = "Chisq")
#We see a residual deviance of 1.40% difference between the models.
#This is probably not worth losing interpretability over.


#-------------------------------------------------------------------------------
#Period 3 models
#Importing data
data("framingham", package = "riskCommunicator")
skim(framingham)
?framingham

#Making the relevant columns into factors
factor_columns <- c("SEX","CURSMOKE","DIABETES","BPMEDS","PREVCHD","PREVAP","PREVMI",
                    "PREVSTRK","PREVHYP","DEATH","ANGINA","HOSPMI","MI_FCHD","ANYCHD",
                    "STROKE","CVD","HYPERTEN","educ","PERIOD")
framingham[factor_columns] <- lapply(framingham[factor_columns], factor)
skim(framingham)

#Removing columns corresponding to responses we do not care about / do
#not have information about at time 0.
#This includes: DEATH, ANGINA, HOSPMI, MI_FCHD, ANYCHD, STROKE, HYPERTEN, 
#TIMEAP, TIMEMI, TIMEMIFC, TIMECHD, TIMESTRK, TIMEDTH, TIMEHYP.
drop <- c("DEATH","ANGINA","HOSPMI","MI_FCHD","ANYCHD","STROKE","HYPERTEN",
          "TIMEAP","TIMEMI","TIMEMIFC","TIMECHD","TIMESTRK","TIMEDTH","TIMEHYP","PREVMI")

framingham <- framingham[ , !(names(framingham) %in% drop)]

#Removing missing Values
framingham_clean <- na.omit(framingham)

framingham_clean_period3 <- framingham_clean[framingham_clean$PERIOD == 3, ]
head(framingham_clean_period3)

form_period3_TOTCHOL <- CVD ~ PREVHYP + PREVAP + SEX + DIABP + DIABETES +
  BPMEDS + educ + AGE + BMI + GLUCOSE + CIGPDAY + TOTCHOL
model_period3_TOTCHOL <- glm(form_period3_TOTCHOL, data = framingham_clean_period3, family = binomial())

#TOTCHOL = LDLC + HDLC

#We exclude LDLC
form_period3_HDLC <- CVD ~ PREVHYP + PREVAP + SEX + DIABP + DIABETES +
  BPMEDS + educ + AGE + BMI + GLUCOSE + CIGPDAY + HDLC
model_period3_HDLC <- glm(form_period3_HDLC, data = framingham_clean_period3, family = binomial())

#Lets look at some residuals, code from chap 10.
model3_TOTCHOL_aug <- augment(model_period3_TOTCHOL, type.residuals = "pearson") |>
  mutate(.fitted_group = cut(.fitted, quantile(.fitted, probs = seq(5, 95, 1) / 100)))

p1_period3_TOTCHOL <- model3_TOTCHOL_aug |>
  ggplot(aes(.fitted, .resid)) + 
  geom_point() + 
  geom_smooth() +
  xlab("fitted values") + 
  ylab("Pearson residuals")

p4_period3_TOTCHOL <- model3_TOTCHOL_aug |>
  group_by(.fitted_group) |>
  summarize(.fitted_local = mean(.fitted), .var_local = var(.resid)) |>
  na.omit() |> 
  ggplot(aes(.fitted_local, .var_local)) + 
  geom_point() + 
  geom_hline(yintercept = 1, linetype = 2) +
  geom_smooth() + 
  xlab("fitted values") +
  ylab("variance")

gridExtra::grid.arrange(p1_period3_TOTCHOL, p4_period3_TOTCHOL, ncol = 2)

model3_HDLC_aug <- augment(model_period3_HDLC, type.residuals = "pearson") |>
  mutate(.fitted_group = cut(.fitted, quantile(.fitted, probs = seq(5, 95, 1) / 100)))

p1_period3_HDLC <- model3_HDLC_aug |>
  ggplot(aes(.fitted, .resid)) + 
  geom_point() + 
  geom_smooth() +
  xlab("fitted values") + 
  ylab("Pearson residuals")

p4_period3_HDLC <- model3_HDLC_aug |>
  group_by(.fitted_group) |>
  summarize(.fitted_local = mean(.fitted), .var_local = var(.resid)) |>
  na.omit() |> 
  ggplot(aes(.fitted_local, .var_local)) + 
  geom_point() + 
  geom_hline(yintercept = 1, linetype = 2) +
  geom_smooth() + 
  xlab("fitted values") +
  ylab("variance")

gridExtra::grid.arrange(p1_period3_HDLC, p4_period3_HDLC, ncol = 2)

#From Chap 3:
drop1(model_period3_TOTCHOL, test = "Chisq") |> 
  tidy() |> 
  arrange(p.value) |>
  filter(term != "<none>") |>
  select(-AIC) |>
  knitr::kable()
#We see that BMI, GLUCOSE, DIABP ought not to be included.

drop1(model_period3_HDLC, test = "Chisq") |> 
  tidy() |> 
  arrange(p.value) |>
  filter(term != "<none>") |>
  select(-AIC) |>
  knitr::kable()
#We see that DIABP, GLUCOSE, BMI ought not to be included.

form_period3_TOTCHOL_new <- CVD ~ PREVHYP + PREVAP + SEX + DIABETES +
  BPMEDS + educ + AGE + + CIGPDAY + TOTCHOL
model_period3_TOTCHOL_new <- glm(form_period3_TOTCHOL_new, data = framingham_clean_period3, family = binomial())

form_period3_HDLC_new <- CVD ~ PREVHYP + PREVAP + SEX + DIABETES +
  BPMEDS + educ + AGE + CIGPDAY + HDLC
model_period3_HDLC_new <- glm(form_period3_HDLC_new, data = framingham_clean_period3, family = binomial())

anova(model_period3_TOTCHOL, model_period3_TOTCHOL_new, test = "Chisq") #Not significant
anova(model_period3_HDLC, model_period3_HDLC_new, test = "Chisq") #Not significant

#Nonlinear considerations
model_period3_TOTCHOL_new_diag <- augment(model_period3_TOTCHOL_new)
p1_3_TOTCHOL <- ggplot(model_period3_TOTCHOL_new_diag, aes(AGE, .std.resid)) +
  stat_binhex(bins = 20) + geom_smooth(linewidth =  1, fill = "blue") +
  xlab("AGE") + ylab("")
p2_3_TOTCHOL <- ggplot(model_period3_TOTCHOL_new_diag, aes(SEX, .std.resid)) +
  stat_binhex(bins = 20) + geom_smooth(linewidth =  1, fill = "blue") +
  xlab("SEX") + ylab("")
p3_3_TOTCHOL <- ggplot(model_period3_TOTCHOL_new_diag, aes(TOTCHOL, .std.resid)) +
  stat_binhex(bins = 20) + geom_smooth(linewidth =  1, fill = "blue") +
  xlab("TOTCHOL") + ylab("")
p4_3_TOTCHOL <- ggplot(model_period3_TOTCHOL_new_diag, aes(PREVHYP, .std.resid)) +
  stat_binhex(bins = 20) + geom_smooth(linewidth =  1, fill = "blue") +
  xlab("PREVHYP") + ylab("")
p5_3_TOTCHOL <- ggplot(model_period3_TOTCHOL_new_diag, aes(PREVAP, .std.resid)) +
  stat_binhex(bins = 20) + geom_smooth(linewidth =  1, fill = "blue") +
  xlab("PREVAP") + ylab("")
p6_3_TOTCHOL <- ggplot(model_period3_TOTCHOL_new_diag, aes(BPMEDS, .std.resid)) +
  stat_binhex(bins = 20) + geom_smooth(linewidth =  1, fill = "blue") +
  xlab("BPMEDS") + ylab("")
p7_3_TOTCHOL <- ggplot(model_period3_TOTCHOL_new_diag, aes(DIABETES, .std.resid)) +
  stat_binhex(bins = 20) + geom_smooth(linewidth =  1, fill = "blue") +
  xlab("DIABETES") + ylab("")
p8_3_TOTCHOL <- ggplot(model_period3_TOTCHOL_new_diag, aes(CIGPDAY, .std.resid)) +
  stat_binhex(bins = 20) + geom_smooth(linewidth =  1, fill = "blue") +
  xlab("CIGPDAY") + ylab("")
p9_3_TOTCHOL <- ggplot(model_period3_TOTCHOL_new_diag, aes(educ, .std.resid)) +
  stat_binhex(bins = 20) + geom_smooth(linewidth =  1, fill = "blue") +
  xlab("educ") + ylab("")
gridExtra::grid.arrange(p1_3_TOTCHOL, p3_3_TOTCHOL, p8_3_TOTCHOL, ncol = 2)

#Nonlinear considerations
model_period3_HDLC_new_diag <- augment(model_period3_HDLC_new)
p1_3_HDLC <- ggplot(model_period3_HDLC_new_diag, aes(AGE, .std.resid)) +
  stat_binhex(bins = 20) + geom_smooth(linewidth =  1, fill = "blue") +
  xlab("AGE") + ylab("")
p2_3_HDLC <- ggplot(model_period3_HDLC_new_diag, aes(SEX, .std.resid)) +
  stat_binhex(bins = 20) + geom_smooth(linewidth =  1, fill = "blue") +
  xlab("SEX") + ylab("")
p3_3_HDLC <- ggplot(model_period3_HDLC_new_diag, aes(HDLC, .std.resid)) +
  stat_binhex(bins = 20) + geom_smooth(linewidth =  1, fill = "blue") +
  xlab("HDLC") + ylab("")
p4_3_HDLC <- ggplot(model_period3_HDLC_new_diag, aes(PREVHYP, .std.resid)) +
  stat_binhex(bins = 20) + geom_smooth(linewidth =  1, fill = "blue") +
  xlab("PREVHYP") + ylab("")
p5_3_HDLC <- ggplot(model_period3_HDLCL_new_diag, aes(PREVAP, .std.resid)) +
  stat_binhex(bins = 20) + geom_smooth(linewidth =  1, fill = "blue") +
  xlab("PREVAP") + ylab("")
p6_3_HDLC <- ggplot(model_period3_HDLC_new_diag, aes(BPMEDS, .std.resid)) +
  stat_binhex(bins = 20) + geom_smooth(linewidth =  1, fill = "blue") +
  xlab("BPMEDS") + ylab("")
p7_3_HDLC <- ggplot(model_period3_HDLC_new_diag, aes(DIABETES, .std.resid)) +
  stat_binhex(bins = 20) + geom_smooth(linewidth =  1, fill = "blue") +
  xlab("DIABETES") + ylab("")
p8_3_HDLC <- ggplot(model_period3_HDLC_new_diag, aes(CIGPDAY, .std.resid)) +
  stat_binhex(bins = 20) + geom_smooth(linewidth =  1, fill = "blue") +
  xlab("CIGPDAY") + ylab("")
p9_3_HDLC <- ggplot(model_period3_HDLC_new_diag, aes(educ, .std.resid)) +
  stat_binhex(bins = 20) + geom_smooth(linewidth =  1, fill = "blue") +
  xlab("educ") + ylab("")
gridExtra::grid.arrange(p1_3_HDLC, p3_3_HDLC, p8_3_HDLC, ncol = 2)
#We see HDLC has some nonlinear effects, which we will include via splines.

form_period3_HDLC_new_nonlin <- CVD ~ PREVHYP + PREVAP + SEX + DIABETES +
  BPMEDS + educ + AGE + CIGPDAY + ns(HDLC, 3)
model_period3_HDLC_new_nonlin <- glm(form_period3_HDLC_new_nonlin, data = framingham_clean_period3, family = binomial())

anova(model_period3_HDLC_new, model_period3_HDLC_new_nonlin, test = "Chisq") #Significant

#Lets look at how much interactions would improve the model:
form_period3_HDLC_new_nonlin_int <- CVD ~ (PREVHYP + PREVAP + SEX + DIABETES +
                                  BPMEDS + educ + AGE + CIGPDAY + ns(HDLC, 3))^2
model_period3_HDLC_new_nonlin_int <- glm(form_period3_HDLC_new_nonlin_int, data = framingham_clean_period3, family = binomial())

drop1(model_period3_HDLC_new_nonlin_int, test = "Chisq") |> 
  tidy() |> 
  arrange(p.value) |>
  filter(term != "<none>") |>
  select(-AIC) |>
  knitr::kable()

#Graphs of interactions if needed.
plot(Effect(c("SEX", "CIGPDAY"), model_period1_new))

anova(model_period3_HDLC_new_nonlin, model_period3_HDLC_new_nonlin_int, test = "Chisq")
#Not significant.

#-------------------------------------------------------------------------------
#Calibration
#For this part we make CVD numeric
framingham_clean_period1$CVD <- as.numeric(as.character(framingham_clean_period1$CVD))
#Code from 10.3.3
group_p <- function(y, p, n_probs = 100) {
  qt <- quantile(p, probs = seq(0, n_probs, 1) / n_probs)
  tibble(
    y = y,
    p = p,
    p_group = cut(p, qt, include.lowest = TRUE)
  ) |>
    group_by(p_group) |>
    summarize(
      p_local = mean(p), 
      obs_local = mean(y), 
      se_local = sqrt(p_local * (1 - p_local) / n()),
    )
}

cal_plot <- ggplot(mapping = aes(p_local, obs_local)) + 
  geom_ribbon(aes(
    ymin = p_local - 1.96 * se_local, 
    ymax = p_local + 1.96 * se_local
  ), fill = "orange", alpha = 0.1) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "orange") + 
  geom_smooth() + 
  xlab("predicted probability") +
  ylab("observed probability")

gridExtra::grid.arrange(
  cal_plot %+% group_p(framingham_clean_period1$CVD, fitted(model_period1_new), n_probs = 50) + 
    coord_cartesian(xlim = c(0, 0.35), ylim = c(0, 0.35)),
  cal_plot %+% group_p(framingham_clean_period1$CVD, fitted(model_period1_new), n_probs = 100) + 
    coord_cartesian(xlim = c(0, 0.35), ylim = c(0, 0.35)),
  ncol = 2
)     
cal_plot %+% group_p(framingham_clean_period1$CVD, fitted(model_period1_new), n_probs = 100) + 
  coord_cartesian(xlim = c(0, 0.35), ylim = c(0, 0.5))

#-------------------------------------------------------------------------------
#Predictive performance
#Code from 10.3.4
auc <- function(y, eta) {
  eta1 <- eta[y == 1]
  eta0 <- eta[y == 0]
  wilcox.test(eta1, eta0)$statistic / (length(eta1) * length(eta0))
}
err <- tibble(
  model = c("Additive (train)", "Interaction (train)"),
  Pearson = c(
    residuals(model_period1_new, type = "pearson")^2 |> mean(),
    residuals(model_period1_new_int, type = "pearson")^2 |> mean()
  ),
  Deviance = c(
    residuals(model_period1_new)^2 |> mean(),
    residuals(model_period1_new_int)^2 |> mean()
  ),
  AUC = c(
    auc(framingham_clean_period1$CVD, predict(model_period1_new)), 
    auc(framingham_clean_period1$CVD, predict(model_period1_new_int))
  )
)
err
#We see that the interaction model is slightly better. Both models have an AUC
#of around 0.75. Indicating some usefulness. Which you prefer depends on wether
#you prefer the slightly higher model accuracy of the interaction version. Or if
#you prefer the ease of interpretability of the additive model.

summary(model_period1_new)
#Which variables has the highest effect in relation to CVD?
#DIABETES: 1.23.
#SEX: -0.865.
#PREVAP: 0.795.
#PREVHYP: 0.484.
