function [nVar rngMin rngMax isInt nObj problem_call] = initializeDent()
nVar = 2;
rngMin = ones(1,nVar).*(-1.5);
rngMax = ones(1,nVar).*(1.5);
isInt = zeros(1,nVar);
nObj = 2;
problem_call = @(x) Dent(x);