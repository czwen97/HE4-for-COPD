library(TwoSampleMR)
#library(gwasglue)
#library(gwasvcf)
#library(ieugwasr)
library(dplyr)
library(MRPRESSO)
set.seed(1234)

setwd("D://work//2025//2025_05//HE4_COPD//Data//")

expd1 <- TwoSampleMR::extract_instruments(outcomes = "prot-a-3224", p1 =1e-05)
r = get_r_from_pn(p = expd1$pval.exposure, n = 3301)
summary(r)
F_s =  (3301 - nrow(expd1) - 1) /  (nrow(expd1) * ( r / (1-r)))
summary(sqrt(F_s))

outd1 <- TwoSampleMR::extract_outcome_data(expd1$SNP, "finn-b-J10_COPD", proxies=FALSE)
dat1 <- TwoSampleMR::harmonise_data(expd1, outd1)

res <- TwoSampleMR::mr(dat = dat1, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_simple_mode"))
res
need_res <- res[,c(5, 7:9)]
setdiff(expd1$SNP, dat1$SNP)

presso_result <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure" , SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                           OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat1, NbDistribution = 1000,  SignifThreshold = 0.05)
need_presso_result <- presso_result$`Main MR results`[1, c(1, 3,4,6)]
need_presso_result[1,1] <- "MR-PRESSO"
colnames(need_presso_result) <- colnames(need_res)
need_res <- rbind(need_res, need_presso_result)

OR = exp(need_res$b)
upper_CI <- exp(need_res$b + 1.96 * need_res$se)
lower_CI <- exp(need_res$b - 1.96 * need_res$se)

OR_str = paste0(round(OR, digits = 4), "[", round(lower_CI,  digits = 4), "-", round(upper_CI,  digits = 4),"]")

need_info <- cbind(need_res[,1:3], OR_str, need_res[,4])
write.table(need_info, "13.txt", col.names = T, row.names = F, sep = "\t", quote = F)

#write.table(dat1, "SNP.txt", col.names = T, row.names = F, sep = "\t", quote = F)

###
mr_heterogeneity(dat1)


###

mr_pleiotropy_test(dat1)


##

mr_leaveoneout_plot(mr_leaveoneout(dat1))

###
mr_scatter_plot(res, dat1)

mr_forest_plot(mr_singlesnp(dat1))
###
library(forestplot)
data.frame(coef = round(OR, digits = 4),
           Beta = round(need_info$b, digits = 4),
           se = round(need_info$se, digits = 4),
           p =  round(need_info$`need_res[, 4]`, digits = 4),
           low = lower_CI,
           high = upper_CI,
           boxsize = rep(0.2, 5),
           OR_str,
           variables = c("MR Egger", "IVW", "Weighted median", "Simple mode", "MR-PRESSO")) |>
  forestplot(labeltext = c(variables, Beta, se, coef, OR_str, p),
             mean = coef,
             lower = low,
             upper = high,
             boxsize = boxsize,
             clip = c(0, 2),
             zero = 1,
             xlog = F) |>
  fp_set_style(lines = "#08403E", box = "#B57114") |>
  fp_add_header(coef = "MR info" |> fp_txt_plain() |> fp_align_center(),
                variables = "Methods")

###
singlesnp_results <- mr_singlesnp(
  dat1,
  parameters = default_parameters(),
  single_method = "mr_wald_ratio",
  all_method = c("mr_ivw", "mr_egger_regression")
)
mr_funnel_plot(singlesnp_results)
