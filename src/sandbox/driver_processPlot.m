% driver_processPlot.m  (run from the repo root)
% Runs processData on the new files, then plots attenuation vs omega.
addpath('./src');

inDir  = './data/simMD/3d/hooke_local';
outMat = './data/simMD/3d/hooke/combined_simData_new.mat';

options = struct('threeD', true);

fprintf('=== processData ===\n');
combined = processData(inDir, outMat, options);

% Quick summary of what we got.
omegaVals   = sort(unique(combined.omega(:)));
gammaVals   = sort(unique(combined.gamma(:)));
pressVals   = sort(unique(combined.pressure_actual(:)));
fprintf('Files processed: %d\n', numel(combined.omega));
fprintf('Unique omega: ');  disp(omegaVals');
fprintf('Unique gamma: ');  disp(gammaVals');
fprintf('Unique P:     ');  disp(pressVals');

fprintf('\n=== plotAttenuation ===\n');
o        = struct();
o.lightMode    = true;   % white background
o.pressureArray = [];    % all pressures

plotAttenuation(outMat, o);

% Export the figure.
fig = gcf;
exportgraphics(fig, './attenuation_vs_omega_highres.png');
fprintf('Saved ./attenuation_vs_omega_highres.png\n');
