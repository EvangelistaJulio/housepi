# Contributing to housepi

*housepi* is developed by the Mathematical Modelling for Infectious Diseases unit at Institut Pasteur.
Report bugs and request features on the
[GitLab issue tracker](https://gitlab.pasteur.fr/mmmi-pasteur/housepi/-/work_items).
We close merge requests from outside the team without review, so open an issue instead.

## Bug reports

Open the issue with the [Bug template](https://gitlab.pasteur.fr/mmmi-pasteur/housepi/-/issues/new?issuable_template=Bug)
and fill in every section:

- the *housepi* version, from `packageVersion("housepi")`;
- a minimal example that reproduces the bug from `library(housepi)`, with `set.seed()` if the result depends on random draws;
- the full error or warning, or the output you got and the output you expected;
- the output of `sessionInfo()`.

The example should run on its own.
If you can't share your data, try to reproduce it on households from `run_simulation()` or on one of the bundled datasets.
If the bug only appears on your data, describe its columns and the number of households and members.

Try the [development version](https://gitlab.pasteur.fr/mmmi-pasteur/housepi#installation) first, since it may already fix the bug.

## Feature requests

Describe what you want to model or compute, the data you have, and what you do today instead.
