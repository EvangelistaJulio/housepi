# Housepi, standalone C++ version

To run your from terminal or to debug your code, you may want to use the standalone version of the `housepi` package.
Here a quick tutorial of the basic commands as if the where run within the folder "inst/standalone".

## Compile a model

```
MODEL=../models/chain_binom
make clean MODEL=$MODEL
make debug MODEL=$MODEL
```

If you don't want tocompile with debug flags (same used by R package compiler), use `make` instead. 
To load an example you can use the following files:

```
PARAMS=../../vignettes/data/params.csv
SETUP=../../vignettes/data/setup.csv
DATABASE=../../vignettes/data/data.csv
```

## Run Simulation 

```
${MODEL}/model.exe -s ${SETUP} -p ${PARAMS} --simulate -d ${DATABASE} -o simulated.tmp -v
```

## Compute Loglikelihood

```
${MODEL}/model.exe -s ${SETUP} -p ${PARAMS} --evaluate -d simulated.tmp -v
```

## Run Inference 

```
${MODEL}/model.exe -s ${SETUP} -p ${PARAMS} --infer -d simulated.tmp -o inferstate.tmp -v
```

## Marginal likelihood

```
${MODEL}/model.exe -s ${SETUP} -p ${PARAMS} --marginal 25 -d simulated.tmp        
```
