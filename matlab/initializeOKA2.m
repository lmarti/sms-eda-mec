function [nVar, rngMin, rngMax, isInt, nObj, algoCall] = initializeOKA2()
nVar = 3;
rngMin = [-pi -5 -5];
rngMax = [pi 5 5];
isInt = zeros(1,nVar);
nObj = 2;
algoCall = 'OKA2';