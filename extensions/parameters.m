% this file contains the parameters for the reaction system 
% defined in "Stability and Stabilization of Constrained Runs"

% volume-- can adjust this to get more/less refinement in the stocastic
% portion of the simulation (larger V = more refinement)
% V = 10000;
V = 100000;

% rate constants
b1 = 5 / V;
d1 = 0.0009 / V;
e1 = 0.1 / V;
f1 = 0.1 / V;
b_1 = 10.6;
% d_1 = 0.005;
d_1 = 0.05;
% e_1 = 0.05;
e_1 = 0.5;
f_1 = 0.01;
b2 = 0.4;
% d2 = 0.085;
d2 = 0.85;
e2 = 0.05;
f2 = 2;

% total concentrations
ST = 1 * V;
ET = 1 * V;
DT = 1 * V;
FT = 0.02 * V;
