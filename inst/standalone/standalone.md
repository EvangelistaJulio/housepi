# How to compile a model from C++

```
MODEL=../inst/models/chain_binom
make MODEL=$MODEL
make clean MODEL=$MODEL
```
If you want to compile with debug flags (same used by R package compiler), use `make debug` instead.

## Run Simulation 

```
PARAMS=../../vignettes/data/params.csv
SETUP=../../vignettes/data/setup.csv
DATABASE=../../vignettes/data/database.csv
${MODEL}/model.exe -simulate -v -p ${PARAMS} -d ${DATABASE} -s ${SETUP} -o simulated.csv
```

## Run Inference 

```
PARAMS=../../vignettes/data/params.csv
SETUP=$../../vignettes/data/setup.csv
${MODEL}/model.exe -infer -v -p ${PARAMS} -d simulated.csv -s ${SETUP} -o delete-me.csv
```

