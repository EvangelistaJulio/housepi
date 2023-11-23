# Epidemic Modeling of Longitudinal Surveys
Epidemics mathematical modeling of houselhold surveys.

## Compile

In order to compile, move to the program folder and run the following:
```
make MODEL="<path/to/model>"
make clean
```
where `<path/to/model>` is the path to the C++ file (with extension *.cpp*).
The command will produce an executable file `<model>.exe` in the same folder of the model.

## Run

To launch the inference from the terminal run the executable file followed by the flag `-infer`:
```
<model>.exe -infer -s <path/to/setup.csv> -d <path/to/database.csv> -p <path/to/parameters.csv> -o <path/to/infer.csv> 
```
In order to print a summary of the execution use the verbose flag `-v`, for other available options use the flag `--help`.

To launch a simulation from the terminal use the flag `simulate`:
```
<model>.exe -simulate -d <path/to/database.csv> -p <path/to/parameters.csv. -o <path/to/simulation.csv>
```

## Analyze

To return an analysis of the inference file you can run the following command
```
Rscript R/analyzer.R -i <path/to/infer.csv> -s <path/to/setup.csv> -p <path/to/parameters.csv> --save_fig
```
For the other available options use the flag `--help`.

## Compute Bayes Factor

To be updated.  