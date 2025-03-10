
hhmodel <- housepi::HouseholdModel$new("chain-binomial")

hhmodel$load_param_template()
hhmodel$loglikelihood(housepi::survey)

combined <- hhmodel$.__enclos_env__$private$.model$.__enclos_env__$private$.rcpp$export_house_size_report()


combined$pct <- sapply(combined$pct, function(x) sprintf("%.1f%%", x))

ft <- flextable::flextable(combined)
ft <- flextable::separate_header(ft)
ft <- flextable::theme_vanilla(ft)
ft <- flextable::align(ft, align = "left", part = "all")
ft <- flextable::bold(ft, j = ~ comment + value, bold = TRUE)
ft <- flextable::autofit(ft, add_w = 0, add_h=0)
ft <- flextable::labelizor(
  ft, part = "all",
  labels = c(comment = "", value = "n.", pct = "[%]")
)
ft <- flextable::colformat_num(ft, big.mark = "")
ft

