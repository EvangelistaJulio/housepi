# Housepi, standalone C++ version
To run your code from terminal or to debug, you may want to use the standalone version of the `housepi` package.
Here a quick tutorial of the basic commands as if they where run within the folder "inst/standalone".

## Compile a model
```
MODEL=../models/augm_onset_single_intr
make -s MODEL=$MODEL
```

This will compile with debug flags (same used by R package compiler), otherwise use
```
make fast -s MODEL=$MODEL
```
To clean use
```
make clean
make clean -s MODEL=$MODEL
```

## Run a command
Load your configuration from files, e. g. create temporary files:
```
PARAM=param.tmp
AUGM=augm.tmp
SETUP=setup.tmp
DATASET=dataset.tmp
DATALIST=datalist.tmp
GROUP=group_rates.tmp
INFEC=infec.tmp
CONFIG="-s ${SETUP} -p ${PARAM} -d ${DATASET} -l ${datalist} -a ${AUGM} -g ${GROUP} -i ${INFEC}"
```
and store your output in
```
OUT=out.tmp
```
### Compute Loglikelihood
```
$MODEL/model.exe ${CONFIG}
```
### Run Simulation 
```
$MODEL/model.exe ${CONFIG} --simulate -o ${OUT}
```
### Run Inference 
```
$MODEL/model.exe ${CONFIG} --infer -o ${OUT}
```
