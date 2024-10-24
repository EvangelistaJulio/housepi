
library(cli)
# cli theme
model_descr_cli_theme = list(
  #ul = list("list-style-type" = "."),
  .gro = list("color" = "cyan"),
  .par = list("font-style" = "bold", "color" = "orange"),
  .aug = list("font-style" = "bold", "color" = "violet")
)
cli_div(theme = model_descr_cli_theme)

#----- Title -----#
cli_h1("Augmented Onsets - Single Introduction {.emph (housepi template)}")

#----- Description -----#
cli_text(
  paste(
    col_red("[MISSING DESCRIPTION]")
  )
)
cli_rule("Likelihood")
ulid <- cli_ul()
cli_li("Force of infection: FOI_i(t) = {.gro Susceptibility}_{{i}} {.par beta} n^(-{.par gamma}) ∑ {.gro Contact}_{{k{symbol$arrow_right}i}} {.gro Infectivity}_{{k}} f_k(t)")
cli_li("Sum over all possible infectors: ∑ = ∑_{{k : t_i > t_k}}")
cli_li("Time of infection: t_i = {.aug onset_date} - {.aug incub_period}")
cli_li("Infectivity profile: f_k(t)")
cli_li("{.gro Susceptibility}_{{i}}")
grolid <- cli_ol()
cli_li("{.param sus_vacc}: vaccinated")
cli_li("[REFERENCE]: others")
cli_end(grolid)
cli_li("{.gro Infectivity}_{{k}}")
grolid <- cli_ol()
cli_li("{.param inf_vacc}: vaccinated")
cli_li("[REFERENCE]: others")
cli_end(grolid)
cli_li("{.gro Contact}_{{k{symbol$arrow_right}i}}")
grolid <- cli_ol()
cli_li("{.param cont_adult_child}: adult {symbol$arrow_right} child")
cli_li("{.param cont_child_adult}: child {symbol$arrow_right} adult")
cli_li("{.param cont_child_child}: child {symbol$arrow_right} child")
cli_li("[REFERENCE]: others")
cli_end(grolid)
cli_end(ulid)
cli_alert_danger("Missing data and asymptomatic states are not modeled.")

#----- Augmented description -----#
cli_rule("Augmented Variables")
cli_dl(
  c(
    "{.aug onset_date}" = "time of symptom onset of infectious individual;",
    "{.aug incub_period}" = "incubation period of infectious individual."
  )
)

#----- Parameter description -----#
cli_rule("Parameters")
cli_dl(
  c(
    "{.emph {.par beta}}" = "person to person transmission rate;",
    "{.emph {.par gamma}}" = "contact dilution parameter;",
    "{.emph {.par sus_vacc}}" = "susceptibility of vaccinated individuals;",
    "{.emph {.par inf_vacc}}" = "infectivity of vaccinated individuals;",
    "{.emph {.par cont_adult_child}} (c)" = "contact rate from adult to child;",
    "{.emph {.par cont_child_adult}} (c)" = "contact rate from child to adult;",
    "{.emph {.par cont_child_child}} (c)" = "contact rate from adult to child;",
    "{.emph {.par study_period}} (c)" = "study period after inclusion;",
    "{.emph {.par child_age}} (c)" = "children/adult age threshold."
  )
)

#----- References -----#
cli_rule("References")
cli_text(
  paste(
    # authors
    "{.emph Cauchemez S, Carrat F, Viboud C, Valleron AJ, Boelle PY}.",
    # title
    "A Bayesian MCMC approach to study transmission of influenza: application to household longitudinal data.",
    # journal, year, volume/pages
    "Statist Med. 2004;23: 3469–3487.",
    # url
    "{.url https://doi.org/10.1002/sim.1912}."
  )
)

cli_end()
