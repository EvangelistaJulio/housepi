# Epidemic Modeling of Longitudinal Surveys
Epidemics mathematical modeling of houselhold surveys.

## Compile

In order to compile, move to the program folder and run the following:
```
make MODEL="<path/to/model>"
make clean
```
where `<path/to/model>` is the path to the C++ file (with extension **.cpp**).
The command will produce an executable file `<model>.exe` in the same folder of the model.

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

## Control files

In order to run your MCMC job, you while need certain control files.

* **Database** (e. g.):
```{csv}
index,status,house_idx,house_size,onset,followup,vaccine,age,sex
1,0,1,3,1000,247,2,55,1,0
2,1,1,3,228,247,0,60,1,1
...
```
where the only mandatory columns are **index** and **status**, while the names of the others are not.
The order of columns is irrelevant. 

* **Parameters** (e. g.):
```{csv}
alpha,0.3,TRUE,uniform,low/0/up/1,lognorm,sdlog/2,low/0/up/1
incub_shape,3.5,FALSE,,,,,
incub_scale,1.1,FALSE,,,,,
...
```

* **Setup** (e. g.):
```{csv}
seed,None
delta_t,0.01
num_iter,40000
num_adapt_iter,5000
thinning,25
```
All entries of the setup are mandatory.

## Analyze

To return an analysis of the inference file you can run the following command
```
Rscript R/analyzer.R -i <infer.csv> -s <setup.csv> -p <parameters.csv> --save_fig
```
For the other available options use the flag `--help`.

## Compute Bayes Factor

To be updated. 
