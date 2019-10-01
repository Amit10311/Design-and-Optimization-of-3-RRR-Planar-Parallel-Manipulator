close all
clear all
clc

% Options including stopping criteria
options = optimset('Display','iter',...
                'Tolfun',1e-6,...
                'Tolx',1e-6,...
                'MaxFunEval',200,...
                'MaxIter',100);

% Starting point
L0 = [30,30];     %Given by the expert, physical problem

% minimum and maximum feasible values of X
lb = [40,60];
ub = [200,200];

% objfun minimization
L = fmincon('objfun',L0,[],[],[],[],lb,ub,'constraints',options);

L          % Final Solution
