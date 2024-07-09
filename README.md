# housepi
Epidemiological modeling and inference of household surveys.

[[INCOMPLETE]]

## Standalone C++ Compilation

In order to compile, move to the standalone folder and run the following:

```
MY_MODEL="<path/to/model>"
make MODEL=${MY_MODEL}
make clean MODEL=${MY_MODEL}
```

where `<path/to/model>` is the path to the implementation **.cpp** file, or header to the header **.h** file. Both need to be in the same folder. 
The command will produce an executable file `<model>.exe` in the same folder of the model.

## Definition files

In order to run your MCMC job, you while need to provide the following files:

* **Database** (e. g.):
```{csv}
index,status,house_index,age,sex
1,0,1,2,55,1
2,1,1,0,60,1
...
```
where the only mandatory columns are **index**, **status** and **house_index**. The order of the columns is irrelevant. 

* **Parameters** (e. g.):
```{csv}
alpha,0.3,TRUE,uniform,low/0/up/1,lognorm,sdlog/2,low/0/up/1
incub_shape,3.5,FALSE,,,,,
incub_scale,1.1,FALSE,,,,,
...
```

Optionally a setup file

* **Setup** (e. g.):
```{csv}
seed,None
delta_t,0.01
num_iter,40000
num_adapt_iter,5000
thinning,25
```

## Run

To launch the **inference**, supposing the executable to be in the same folder as the control files, use the following command from terminal:
```
<model>.exe -infer -s <setup.csv> -d <database.csv> -p <parameters.csv> -o <infer.csv> 
```
In order to print a summary of the execution use the verbose flag `-v`, for other available options use the flag `--help`.

To launch a **simulation** use the flag `simulate`:
```
<model>.exe -simulate -d <database.csv> -p <parameters.csv. -o <simulation.csv>
```