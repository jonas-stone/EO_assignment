
clear all; close all; clc;

addpath("module_geometry_preprocess\");
addpath("Q3D\");

% b2, c_root, c_tip, twist_tip, v_inf, alpha
x0 = [7.5, 0.9, 0.35, -2, 25, 3];

[LD,stress_crit,L] = run_model(x0);