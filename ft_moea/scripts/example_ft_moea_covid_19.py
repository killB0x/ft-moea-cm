#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 27

EXAMPLE USING FT-MOEA FOR USING THE CASE STUDY COVID-19

Here is presented an example to infer a Fault Tree model providing a failure data set.
Here it is used the one associated to one of our case studies, the COVID-19 infection 
risk fault tree by Bakeli et al. [1].

Input and output parameters of FT-MOEA are as follows:

Input parameters:
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
Output:
    * Inferred Fault Trees (fts): this list contains all the objects of inferred Fault Trees of the last generation, sorted from the worst to the best Fault Tree.
    * Convergence time per generation (t): computational time per generation.
    * Metrics: is an array of three columns, corresponding respectively to \phi_c, \phi_s, and \phi_d. And the rows are respectively associated with the inferred Fault Trees contained in (fts).


References:
[1] T. Bakeli, A. A. Hafidi et al., “Covid-19 infection risk management during 
construction activities: An approach based on fault tree analysis (fta),” Journal 
of Emergency Management, vol. 18, no. 7, pp. 161–176, 2020.

@author: Lisandro A. Jimenez Roa (l.jimenezroa@utwente.nl)
"""
from ft_moea import read_dataset_counts, FTMOEA, getMCSs
from scipy.io import loadmat

#%% Load the failure data set
# Generate the failure data set:
data = loadmat('failure_data_sets/covid-19_noise_0_0sfv.mat')['data_before_noise']
#%% FT-MOEA input parameters
ps = 250
ft_as_input = ''
max_gen = 100
ug = 10
obj_functions = [1,1,1] # {\phi_c,\phi_s, \phi_d}
#%% Running FT-MOEA
be_ss ,failure_data_set  = read_dataset_counts(data)
fts,t,metrics = FTMOEA(population_size=ps,
MCSs = getMCSs(data[data[:,-2]==1,0:-2]),
parentFT = '',
ng=max_gen,
convergence_criterion=ug,
multi_objective_function = obj_functions, 
dataset = failure_data_set,
be_s = be_ss,
path_save_results = 'Example',
)

# Inferred FT
print(fts[-1])
print(metrics[-1])
