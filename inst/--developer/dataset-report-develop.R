
hhmodel <- housepi::HouseholdModel$new("chain-binomial")
hhmodel$load_param_template()
hhmodel$loglikelihood(housepi::survey)

combined <- hhmodel$.__enclos_env__$private$.model$.__enclos_env__$private$.rcpp$export_data_report()

hhmodel$.__enclos_env__$private$.model$.__enclos_env__$private$.rcpp$export_house_size_report()

all(hhmodel$dataset == housepi::survey)

df <- hhmodel$simulate()

combined$pct <- sapply(combined$pct, function(x) sprintf("%.1f%%", x))

ft <- flextable::flextable(combined)
ft <- flextable::separate_header(ft)
ft <- flextable::theme_vanilla(ft)
ft <- flextable::align(ft, align = "left", part = "all")
ft <- flextable::merge_v(ft, j = ~ type)
ft <- flextable::bold(ft, j = ~ type, bold = TRUE)
ft <- flextable::autofit(ft, add_w = 0, add_h=0)
ft <- flextable::labelizor(
  ft, part = "all",
  labels = c(type = "", comment = "", pct = "[%]")
)
ft <- flextable::colformat_num(ft, big.mark = "")
ft

