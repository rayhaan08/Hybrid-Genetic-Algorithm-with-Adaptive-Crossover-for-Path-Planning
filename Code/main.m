%% This code calls the GA_Berlin for r runs and m salesman.
% ------------------------------------------------------------------------------------------------------------------
% INPUTS
% m: No of salesman
% r: No of runs to be performed
% ------------------------------------------------------------------------------------------------------------------
% OUTPUTS
% Opt: Optimum longest tour
% Opt_design: Optimum route in two-part chromosome form
%
% Copyright - Rayhaan Iqbal (2020)
% ADAMS Lab, UB
%% Enter parameters
m = 5;
r = 1;
Opt, Opt_design = GA_Berlin(m,r);
