This folder contains the different failure data sets generated through the Monte Carlo simulation. Here we used the ground truth Fault Trees and used p=0.5 for all the basic events failure probability. We draw N=250.000 samples plus Minimal Cut Sets (to ensure completeness in the failure data set).

The way to read the name of the files is as follows:

[name-of-the-case-study]_noise_X_Ysfv.mat

Where:
- [name-of-the-case-study]: we have 6 case studies. In the paper can be found the names and references.
- X is the level of noise in the failure data set. Being zero for all the cases.
- Y is the number of superfluous variables in the failure data set.

From the content of the .mat, the variable "data_before_noise" is the one needed to input in the ft-moea algorithm. The string_ft corresponds to the ground truth in its string form.


