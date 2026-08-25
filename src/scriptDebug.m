clear all
close all

%% Page 1 Left side old mat spring damped, and perturb in [dx,dy] of eigenvector
simEigenmode("./OLD.mat", scalModeFreq = .3, loadPolyData=true)

%% Page 1 Right side old mat spring damped, and perturb in [dx,-dy] of eigenvector
simEigenmode("./OLD.mat", scalModeFreq = .3, loadPolyData=true, negDy=true)

%% Page 2 left side old mat spring damped, and perturbe in [dx,-dy] of eigenvector
simEigenmode("./NEW.mat", scalModeFreq = .3, loadPolyData=true)
% compare this to page 1 right hand side - matches

%% produces plot from page 3 left hand side
simEigenmode("./results_2D_iso_N100_P0.1_Seed5_gamma_3.00e-03.mat", scalModeFreq=.3)

%% Produces plot from page 3 right hand side
simEigenmode("./results_2D_iso_N100_P0.1_Seed5_gamma_1.00e-02.mat", scalModeFreq=.3, negDy=true)





