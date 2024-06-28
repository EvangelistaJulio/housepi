# How to compile a model from C++

```
MODEL=../models/chain_binomial-heter/model.cpp
EXE=${MODEL%.cpp}.exe
make MODEL=$MODEL
make clean MODEL=$MODEL
```
If you want to compile with debug flags (same used by R package compiler), use `make debug` instead.

```
PARAMS=../example/params.csv
DATABASE=../example/database.csv
SETUP=../example/params.csv
${EXE} -infer -v -p ${PARAMS} -d ${DATABASE} -s ${SETUP} -o delete-me.csv
```