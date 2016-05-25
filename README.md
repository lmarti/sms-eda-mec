# `sms-eda-mec`
#### Matlab reference implementation for the S-Metric Selection Estimation of Distribution Algorithm based on Multivariate Extension of Copulas (SMS-EDA-MEC)

## Summary

It can be argued that  in order to produce a substantial improvement in multi-objective estimation of distribution algorithms it is necessary to focus on a particular group of issues, in particular, on the weaknesses derived from multi-objective fitness assignment and selection methods, the incorrect treatment of relevant but isolated (precursor) individuals; the loss of population diversity, and the use of 'general purpose' modeling algorithms without taking note of the particular requirements of the task. 

In this work we introduce the S-Metric Selection Estimation of Distribution Algorithm based on Multivariate Extension of Copulas (SMS-EDA-MEC). SMS-EDA-MEC was devised with the intention of dealing with those issues in mind. It builds the population model relying on the comprehensive Clayton's copula and incorporates methods for automatic population restarting and for priming precursor individuals. The experimental studies presented show that SMS-EDA-MEC yields better results than current and 'traditional' approaches. 

For further information and results see: [http://lmarti.com/smsedamec-cec2016](http://lmarti.com/smsedamec-cec2016).

## Papers

The following papers have to do with the algorithm: 

* Luis Martí, Harold D. de Mello Jr., Nayat Sanchez-Pi and Marley Vellasco (2016) SMS-EDA-MEC: Extending Copula-based EDAs to Multi-Objective Optimization, 2016 IEEE Conference on Evolutionary Computation (CEC'2016), part of 2016 IEEE World Congress on Computational Intelligence (WCCI'2016), Vancouver, Canada, *in press*.

## Usage

The function `multi_mec_eda_clayton.m` is the main entry point. Try running in matlab
```
>> help multi_mec_eda_clayton
```
to get the parameters details.

For example, runnning 
```
>> multi_mec_eda_clayton('WFG1')
```
will run SMS-EDA-MEC on WFG1 with default parameters.

Problems currently implemented are WFG1-9, Dent, GSP and OKA2.

**Note:** You will need to compile `hv.cpp` and `paretofront.c` within Matlab.

For compiling the hypervolume run:
```bash
>> mex -I. hv.cpp Hypervolume.cpp
``` 

For the `paretofront` function run:
```bash
>> mex paretofront.c
```


## Credits

Parts of SMS-EDA-MEC code is based on the Matlab implementation of SMS-EMOA and OCD by Fabian Kretzschmar and Tobias Wagner. See [https://ls11-www.cs.uni-dortmund.de/rudolph/hypervolume/start](https://ls11-www.cs.uni-dortmund.de/rudolph/hypervolume/start) for further details.

* WFG problems implementation (`wfg.m`) by Robin Purshouse (University of Sheffield).
* Hypervolume implementation by Thomas Vo&szlig; (Ruhr-Universit&auml;t Bochum) as part of the [Shark](http://image.diku.dk/shark/) library.
* `paretofront.m` by Yi Cao (Cranfield University).
* For other cases, check the function documentation for author acknowledgements.
