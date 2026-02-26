
clear all; close all; clc;

addpath("module_geometry_preprocess\");
addpath("Q3D\");

    % b, c_root, c_tip, twist_tip, v_inf, alpha
x = [4,1,0.7,5,20,3];

[LD,stress_crit,L] = run_solver(x);