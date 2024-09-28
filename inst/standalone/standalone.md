# Housepi, standalone C++ version

To run your from terminal or to debug your code, you may want to use the standalone version of the `housepi` package.
Here a quick tutorial of the basic commands as if the where run within the folder "inst/standalone".

## Finite Size Model

### Compile a model
```
MODEL=../models/augm_onset_single_intr
make clean MODEL=$MODEL; make debug MODEL=$MODEL
```
If you don't want to compile with debug flags (same used by R package compiler), use `make` instead. 
Load your configuration from files:
```
PARAMS=params.tmp
AUGM=augm.tmp
SETUP=setup.tmp
DATA=../../vignettes/survey.csv
```

### Run Simulation 
```
${MODEL}/model.exe -p ${PARAMS} -a ${AUGM} -s ${SETUP} -d ${DATA} --simulate -o simulated.tmp -v
```
### Compute Loglikelihood
```
${MODEL}/model.exe -p ${PARAMS} -a ${AUGM} -s ${SETUP} -d simulated.tmp --evaluate -v
```
### Run Inference 
```
${MODEL}/model.exe -p ${PARAMS} -a ${AUGM} -s ${SETUP} -d simulated.tmp --infer -v
```
### Marginal likelihood
```
${MODEL}/model.exe -p ${PARAMS} -a ${AUGM} -s ${SETUP} -d simulated.tmp --marginal 25 -v
```