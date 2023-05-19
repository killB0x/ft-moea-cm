# SymLearn: Data-driven inference of fault tree models for symmetric systems

The SymLearn toolchain infers fault tree models from given failure data sets.

## Preparation
SymLearn requires the following Python libraries.
- pandas
- numpy
- sympy
- pyeda
- scipy
- sortedcontainers

Install the missing dependencies with:
```
python setup.py develop
```

## Usage
An example invocation of the toolchain is:
```
python learn_ft.py ../input_data/simple_case_OR.mat
```

The following commandline arguments are useful:
- `--learn-approach ftmoea` to define the back-end (either `ftmoea`,`sympy` or `espresso`)
- `--disable-symmetries` to disable the usage of symmetries
- `--disable-modules` to disable the usage of modules
- `--disable-recursion` to disable any recursive calls

More arguments are listed when executing the help function:
```
python learn_ft.py -h
```
