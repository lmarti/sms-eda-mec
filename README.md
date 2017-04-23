# `sms-eda-mec`
#### Matlab reference implementation for the *S*-Metric Selection Estimation of Distribution Algorithm based on Multivariate Extension of Copulas (SMS-EDA-MEC)


[![DOI](https://zenodo.org/badge/DOI/10.1109/CEC.2016.7744261.svg)](https://doi.org/10.1109/CEC.2016.7744261)
[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg)](http://creativecommons.org/licenses/by-nc-sa/4.0/)

## Summary

It can be argued that  in order to produce a substantial improvement in multi-objective estimation of distribution algorithms it is necessary to focus on a particular group of issues, in particular, on the weaknesses derived from multi-objective fitness assignment and selection methods, the incorrect treatment of relevant but isolated (precursor) individuals; the loss of population diversity, and the use of 'general purpose' modeling algorithms without taking note of the particular requirements of the task.

In this work we introduce the S-Metric Selection Estimation of Distribution Algorithm based on Multivariate Extension of Copulas (SMS-EDA-MEC). SMS-EDA-MEC was devised with the intention of dealing with those issues in mind. It builds the population model relying on the comprehensive Clayton's copula and incorporates methods for automatic population restarting and for priming precursor individuals. The experimental studies presented show that SMS-EDA-MEC yields better results than current and 'traditional' approaches.

For further information and results see: [http://lmarti.com/smsedamec-cec2016](http://lmarti.com/smsedamec-cec2016).

## Citing SMS-EDA-MEC

The following paper contains a description of the algorithm.

* Luis Marti, Harold D. de Mello Jr., Nayat Sanchez-Pi and Marley Vellasco (2016) SMS-EDA-MEC: Extending Copula-based EDAs to Multi-Objective Optimization, 2016 IEEE Conference on Evolutionary Computation (CEC'2016), part of 2016 IEEE World Congress on Computational Intelligence (WCCI'2016), Vancouver, Canada, pp. 3726--3733. doi: [10.1109/CEC.2016.7744261](http://dx.doi.org/10.1109/CEC.2016.7744261)
```bibtex
@inproceedings{marti-2016:cec,
  author = {Luis Mart\'i and Dias de Mello Jr, Harold and Nayat Sanchez-Pi and Rebuzzi Vellasco, Marley~M.~B.},
  title = {{SMS-EDA-MEC}: {E}xtending Copula-based {EDAs} to Multi-Objective Optimization},
  booktitle = {Proceedings of the 2016 IEEE Congress on Evolutionary Computation (IEEE CEC 2016) part of the 2016 IEEE World Congress on Computational Intelligence (IEEE WCCI 2016)},
  year = {2016},
  publisher = {IEEE Press},
  month = {7},
  location = {Vancouver, Canada},
  pages={3726--3733},
  doi={10.1109/CEC.2016.7744261},
}
```

## Requirements

You will need to compile `hv.cpp` and `paretofront.c` (distributed along) within Matlab.

For compiling the hypervolume (`hv.cpp`) run:
```matlab
>> mex -I. hv.cpp Hypervolume.cpp
```

For the `paretofront.c` function run:
```matlab
>> mex paretofront.c
```

## Usage

The function `sms_eda_mec.m` is the main entry point. Try running in Matlab
```matlab
>> help sms_eda_mec
```
to get the parameters details.

For example, running
```matlab
>> sms_eda_mec('WFG1')
```
will run SMS-EDA-MEC on WFG1 with default parameters.

Similarly, running
```matlab
>> sms_eda_mec('WFG1', struct('show_plots', true))
```
will show the populations as evolution takes place if dealing with a 2-objective problem.

Problems currently implemented are WFG1-9, Dent, GSP and OKA2.

## Credits

Parts of SMS-EDA-MEC code is based on the Matlab implementation of SMS-EMOA and OCD by Fabian Kretzschmar and Tobias Wagner. See [https://ls11-www.cs.uni-dortmund.de/rudolph/hypervolume/start](https://ls11-www.cs.uni-dortmund.de/rudolph/hypervolume/start) for further details.

* WFG problems implementation (`wfg.m`) by Robin Purshouse (University of Sheffield).
* Hypervolume implementation by Thomas Vo&szlig; (Ruhr-Universit&auml;t Bochum) as part of the [Shark](http://image.diku.dk/shark/) library.
* `paretofront.m` by Yi Cao (Cranfield University).
* For other cases, check the function documentation for author acknowledgements.
