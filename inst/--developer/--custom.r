library(housepi)

hhmodel_1 <- new_housepi_model("inst/examples/custom")
expose_function("growth", hhmodel = hhmodel_1)$func(seq(3))

hhmodel_2 <- new_housepi_model("blank")
expose_function("growth", hhmodel = hhmodel_2)$func(seq(3))

hhmodel_3 <- new_housepi_model("inst/examples/custom")
expose_function("growth", hhmodel = hhmodel_3)$func(seq(3))
