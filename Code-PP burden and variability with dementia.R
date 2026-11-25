

#1 Definition of PP burden and variability 
#example
CHARLS_PP_burden <- CHARLS_data_6[
  !is.na(PP) & !is.na(Age)][
    order(ID, Age), {
      pp_excess <- pmax(PP - 60, 0)
      yr <- Age
      high <- PP >= 60
      PP_mean <- mean(PP, na.rm = TRUE)
      follow_time <- Age - min(Age)
      PP_SD <- if (.N >= 2) {
        sd(PP, na.rm = TRUE)
      } else {
        NA_real_
      }
      PP_CV <- if (.N >= 2 && PP_mean != 0) {
        PP_SD / PP_mean * 100} else { NA_real_ }
      # ARV
      PP_ARV <- if (.N >= 2) {
        mean(abs(diff(PP)), na.rm = TRUE)} else {NA_real_}
      #slope
      PP_slope <- if (.N >= 2 && length(unique(follow_time)) >= 2) {
        coef(lm(PP ~ follow_time))[["follow_time"]]
      } else {
        NA_real_
      }
      
      if (.N == 1) {
        cPP_auc <- NA_real_
        cPP_excess_auc <- NA_real_
        fu <- NA_real_} 
      else {
        cPP_auc <- sum(((PP[-.N] + PP[-1]) / 2) * diff(yr))
        cPP_excess_auc <- sum(((pp_excess[-.N] + pp_excess[-1]) / 2) * diff(yr))
        fu <- max(yr) - min(yr) }
      .(cPP = cPP_auc,
        cPP_load = cPP_excess_auc / cPP_auc,  #PP >60 mmHg 的超额暴露面积 / 总 PP 暴露面积
        cPP_time = fu,
        PP_mean = PP_mean,
        PP_SD = PP_SD,
        PP_CV = PP_CV,
        PP_ARV = PP_ARV,
        PP_slope = PP_slope)
    },by = ID]
summary(CHARLS_PP_burden)
CHARLS_VIM_model <- lm(
  log(PP_SD) ~ log(PP_mean),
  data = CHARLS_PP_burden[
    !is.na(PP_SD) &
      !is.na(PP_mean) &
      PP_SD > 0 &
      PP_mean > 0 ])
CHARLS_beta <- coef(CHARLS_VIM_model)[2]
# Calculate VIM for each participant
CHARLS_PP_burden[, PP_VIM :=  100 * PP_SD / (PP_mean ^CHARLS_beta)]


#2\ create categorical variables

Pool_burden_data_wide[, cPP_load_tertile := {
  median_nonzero <- median(cPP_load[cPP_load > 0], na.rm = TRUE)
  factor( fifelse(cPP_load == 0, "PP_T1",
                  fifelse(cPP_load < median_nonzero, "PP_T2", "PP_T3")),
          levels = c("PP_T1", "PP_T2", "PP_T3"))}, by = Study]
summary(Pool_burden_data_wide$cPP_load_tertile)
Pool_burden_data_wide[, PP_VIM_tertile := {
  q <- quantile(PP_VIM, probs = c(1/3, 2/3), na.rm = TRUE)
  factor(
    fifelse(PP_VIM <= q[1], "VIM_T1",
            fifelse(PP_VIM <= q[2], "VIM_T2", "VIM_T3")),
    levels = c("VIM_T1", "VIM_T2", "VIM_T3")
  )}, by = Study]
summary(Pool_burden_data_wide$PP_VIM_tertile)

#3*3
Pool_burden_data_wide[, PP_burden_VIM_9cat := factor(
  paste0(cPP_load_tertile, "_", PP_VIM_tertile),
  levels = as.vector(
    outer(
      c("PP_T1", "PP_T2", "PP_T3"),
      c("VIM_T1", "VIM_T2", "VIM_T3"),
      FUN = function(x, y) paste0(x, "_", y)
    )
  )
)]




#3\ Survival Analysis - Cox Model

#two-step IPD
coxph(Surv(follow_time,IncidentDementia_num) ~cPP_load+Age +SEX+Covariates,
      data =Study)

frame_1 <-  metagen(TE = B,seTE = std,data = frame_1_CIND,sm = "HR",studlab = Study,common =F,subgroup = `Subgroups`)

meta::forest(frame_1,leftcols="studlab",test.subgroup.random=F,overall=F,col.subgroup="black",
             print.subgroup.name=F,overall.hetstat=F,fontsize=10,header.line="below",digits = 3)

#one-step IPD

coxph(Surv(follow_time,IncidentDementia_num) ~cPP_load+Age +SEX+strata(Study),
      data = Pool)




