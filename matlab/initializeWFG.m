function [nVar rngMin rngMax isInt nObj problem_call] = initializeWFG(prob_num)
dist_params = 20;
pos_params = 10;

nVar = dist_params + pos_params;
rngMin = zeros(1,nVar);
rngMax = 2 * (1:nVar);
isInt = zeros(1,nVar);
nObj = 2;

problem_call = @(x) wfg(x, nObj, pos_params, dist_params, str2num(prob_num));