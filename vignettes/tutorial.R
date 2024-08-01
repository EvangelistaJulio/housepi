## ----echo=FALSE, out.width = "100%"-------------------------------------------
knitr::include_graphics("images/cover-light.png")

## -----------------------------------------------------------------------------
library(housepi)

## -----------------------------------------------------------------------------
hemi <- generate_model("chain_binom")

## ----eval=FALSE---------------------------------------------------------------
#  template_model(which = "chain_binom", where = "path/to/model")

## ----eval=FALSE---------------------------------------------------------------
#  hemi <- generate_model("path/to/model")

## -----------------------------------------------------------------------------
build_model(hemi,
            data = "data/database.csv",
            params = "data/params.csv",
            setup = "data/setup.csv")

# show the parameter definitions
hemi$show_params()

# show the setup definitions
hemi$show_setup()

