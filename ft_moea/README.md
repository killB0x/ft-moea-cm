# FT-MOEA: Automatic inference of fault tree models via multi-objective evolutionary algorithms

FT-MOEA [1] is an algorithm in Python library that infers Fault Tree models using multi-objective evolutionary algorithms.

## Repository overview

This repository counts with two main folders. In the folder _scripts_ find the implementation in Python of the FT-MOEA together a couple of examples. In the folder _figures_ find the scripts in Matlab and data to reproduce the figures in [1].


# Introduction

We present a novel approach to automatically infer Fault Tree (FT) models from a failure data set via multi-objective evolutionary algorithms. This with the aim to overcome one of the main drawbacks of FTs related to their construction, traditionally carried out in conjunction with domain expertise and in a hand-crafted manner, resulting in a tedious and time-consuming task. In the case of complex industrial systems, manual development of these models can lead to incompleteness, inconsistencies, and even errors.

The first attempt using evolutionary algorithms with this purpose was carried out by Linard et al. [2], here the authors created an algorithm ([FT-GA](https://gitlab.science.ru.nl/alinard/learning-ft)) to generate FT models from a labeled binary failure data set using a uni-dimensional cost function based on the accuracy, which is the proportion of correctly predicted top events by a given FT. A drawback of this approach is that it focuses solely on achieving high accuracy without considering the FT structure. The latter has negative implications such as: running the algorithm twice for the same failure data set may yield considerably different FT structures; it is prone to complexity explosion, leading to massive FTs that are difficult to handle; long computational time, and bad convergence of the algorithm.

We tackle this challenge via multi-objective evolutionary algorithms (MOEAs), where the main contributions are that (i) we show it is possible to achieve more consistent and efficient FT structures via MOEAs, which simultaneously minimize the fault tree size, the error computed based on the failure data set, and the error based on Minimal Cut Sets (ii) we propose a metric to compare FT structures via Minimal Cut Sets using the RV-coefficient, (iii) we carry out a parametric analysis that explains the performance of the algorithm under different assumptions.

# Downloading sources

You can clone the repository via the following command:

```
$ git clone https://gitlab.utwente.nl/fmt-ut/ft-moea.git
```

## Dependencies
* Python 3

## Additional packages
* [NumPy](https://numpy.org/install/): ` pip install numpy `
* [regex](https://pypi.org/project/regex/): ` pip install regex `
* [sympy](https://www.sympy.org/en/index.html): ` pip install sympy `
* [scipy](https://www.scipy.org/install.html): ` pip install scipy `
* [itertools](https://pypi.org/project/more-itertools/): ` pip install more-itertools `

# Multi-objective evolutionary algorithm
Below are enlisted the input parameters and output of the FT-MOEA.

## Input
We need to provide the algorithm with (i) failure data set, and (ii) initial parameters

### Failure data set
You can generate your own failure data set by proposing a fault tree and the Monte Carlo method. In the file _scripts/genFailureDataSet_MC.py_ is provided a self-explanatory example on how to do this (see also example _scripts/example_ft_moea_from_ft_string.py_). 

```
python3 example_ft_moea_from_ft_string.py
```

If you have a failure data set with the required format, you can explicitly indicate the path and the file name as part of the FT-MOEA input (see example _scripts/example_ft_moea_covid_19.py_). Examples of these failure data sets can be found into the folder scripts/failure_data_sets. 

```
python3 example_ft_moea_covid_19.py
```

### Initial parameters

* Population_size (ps): the minimal population size maintained throughout generations after applying the genetic operations
* Minimal Cut Sets (MCSs): this is optional, only applicable on noise-free and complete failure data sets. Based on this metric the error based on Minimal Cut Sets (\phi_c) is computed.
* Parent Fault Tree (parentFT): this is optional, entry a Fault Tree string to be used as parent Fault Tree. If empty, it is going to be initialized with two-parent Fault Trees, one with all BEs connected to an AND gate, and the other with all BEs connected to an OR gate. 
* Max. number of generations (ng): Terminates the optimization process if the number of generations exceeds $ng$ and none of the other convergence criteria is met. 
* Max. generations with unchanged best candidate (ug): if after $uc$ number of generations the best individual (i.e., the FT with the smallest size, and smallest error(s) within the best Pareto set) remains unchanged, then we assume the process has converged and is therefore terminated.
* multi_objective_function: It is a list with three terms that switches on (1) and off (0) the metrics to be optimized.
    - The first term corresponds to the error based on Minimal Cut Sets (\phi_c)
    - The second term corresponds to the size of the Fault Tree (\phi_s)
    - The third term corresponds to the error based on the failure data set (\phi_d)
* Failure data set (dataset): corresponds to the failure data set.
* Basic events (be_s): associated with the basic events considered in the inference process.
* Path to save results externally (path_save_results): when the path of a folder is provided, the results per generation are saved in a .mat file.

## Output

* Inferred Fault Trees (fts): this list contains all the objects of inferred Fault Trees of the last generation, sorted from the worst to the best Fault Tree.
* Convergence time per generation (t): computational time per generation.
* Metrics: is an array of three columns, corresponding respectively to \phi_c, \phi_s, and \phi_d. And the rows are respectively associated with the inferred Fault Trees contained in (fts).


# Other useful files

* The file _scripts/visualization_ft_graph.m_ is a Matlab script that takes as input a string and returns an schematic representation of a Fault Tree. Example:

```
visualization_ft_graph('AND(BE4,BE8,OR(BE2,AND(BE8,OR(BE7,BE5,BE9,BE6,BE3),BE4,BE1)))')
```

# References
[1] L.A. Jimenez-Roa, T. Heskes, T. Tinga, and M. Stoelinga, "Automatic inference of fault tree models via multi-objective evolutionary algorithms", 2021 (submitted and under review).

[[2] A. Linard, D. Bucur, and M. Stoelinga, “Fault trees from data: Efficient learning with an evolutionary algorithm,” in International Symposium on Dependable Software Engineering: Theories, Tools, and Applications. Springer, 2019, pp. 19–37.](https://link.springer.com/chapter/10.1007/978-3-030-35540-1_2)

# Other sources used in this work

* From _pymoo - Multi-Objective Optimization Framework_ the [NSGA-II](https://github.com/anyoptimization/pymoo/blob/20abef1ade71915352217400c11ece4c2f35163e/pymoo/algorithms/nsga2.py) algorithm

# Additional notes

This research has been partially funded by Dutch Research Council (NWO) under the grant PrimaVera ([https://primavera-project.com](https://primavera-project.com)) number NWA.1160.18.238.


# Contact

[Lisandro Jimenez](https://people.utwente.nl/l.jimenezroa) <br>
l.jimenezroa@utwente.nl <br>
PhD candidate in Computer Science <br>
Formal Methods & Tools (FMT) group <br>
Faculty of Electrical Engineering, Mathematics and Computer Science (EEMCS) <br>
University of Twente, The Netherlands <br>
+31 53 489 15 60‬


