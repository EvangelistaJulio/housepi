# How to compile a model from C++

```
MODEL=../models/chain_binom
EXE=${MODEL}/model.exe
make MODEL=$MODEL
make clean MODEL=$MODEL
```
If you want to compile with debug flags (same used by R package compiler), use `make debug` instead.

## Run Simulation 

```
PARAMS=../example/params.csv
DATABASE=../example/database.csv
SETUP=../example/setup.csv
${EXE} -simulate -v -p ${PARAMS} -d ${DATABASE} -s ${SETUP} -o simulated.csv
```

## Run Inference 

```
PARAMS=../example/params.csv
SETUP=../example/setup.csv
${EXE} -infer -v -p ${PARAMS} -d simulated.csv -s ${SETUP} -o delete-me.csv
```

