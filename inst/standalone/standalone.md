# Housepi, standalone C++ version

To run your from terminal or to debug your code, you may want to use the standalone version of the `housepi` package.
Here a quick tutorial of the basic commands as if the where run within the folder "inst/standalone".

## Compile a model

```
MODEL=../models/chain_binom
make MODEL=$MODEL
make clean MODEL=$MODEL
```
If you want to compile with debug flags (same used by R package compiler), use `make debug` instead.

## Compute Loglikelihood

```
PARAMS=../../vignettes/data/params.csv
SETUP=../../vignettes/data/setup.csv
DATABASE=../../vignettes/data/database.csv
${MODEL}/model.exe -s ${SETUP} -p ${PARAMS} -d ${DATABASE} -v
```

## Run Simulation 

```
${MODEL}/model.exe -s ${SETUP} -p ${PARAMS} --simulate -d ${DATABASE} -o simulated.tmp -v
```

## Run Inference 

```
${MODEL}/model.exe -s ${SETUP} -p ${PARAMS} --infer -d simulated.tmp -o inferstate.tmp -v
```
