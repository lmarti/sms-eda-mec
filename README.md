# `sms-eda-mec`
### Matlab source code for the SMS-EDA-MEC multi-objective estimation of distribution algorithm

## Summary

It can be argued that  in order to produce a substantial improvement in multi-objective estimation of distribution algorithms it is necessary to focus on a particular group of issues, in particular, on the weaknesses derived from multi-objective fitness assignment and selection methods, the incorrect treatment of relevant but isolated (precursor) individuals; the loss of population diversity, and the use of 'general purpose' modeling algorithms without taking note of the particular requirements of the task. 

In this work we introduce the S-Metric Selection Estimation of Distribution Algorithm based on Multivariate Extension of Copulas (SMS-EDA-MEC). SMS-EDA-MEC was devised with the intention of dealing with those issues in mind. It builds the population model relying on the comprehensive Clayton's copula and incorporates methods for automatic population restarting and for priming precursor individuals. The experimental studies presented show that SMS-EDA-MEC yields better results than current and 'traditional' approaches. 

For further information and results see: [http://lmarti.com/sms-eda-mec](http://lmarti.com/sms-eda-mec).

## Papers

* Luis Martí, Harold D. de Mello Jr., Nayat Sanchez-Pi and Marley Vellasco (2016) SMS-EDA-MEC: Extending Copula-based EDAs to Multi-Objective Optimization, 2016 IEEE Conference on Evolutionary Computation (CEC'2016), part of 2016 IEEE World Congress on Computational Intelligence (WCCI'2016), Vancouver, Canada.

## Usage

The function `multi_mec_eda_clayton.m` is the main entry point. Try running in matlab
```
>> help multi_mec_eda_clayton
```
to get the parameters details.

Problems currently implemented are WFG1-9, Dent , GSP and OKA2.

**Note:** You will need to compile `hv.cpp` and `paretofront.c` within matlab. Refer to the `matlab\hv.readme` and `matlab\paretofront.readme` files for instructions.

## Credits

Parts of the code is based on the Matlab implementations of SMS-EMOA and OCD by Fabian Kretzschmar and Tobias Wagner. See [https://ls11-www.cs.uni-dortmund.de/rudolph/hypervolume/start](https://ls11-www.cs.uni-dortmund.de/rudolph/hypervolume/start) for further details.
