% tests processEigenModesDampedPara.m
% should take all the packings (*.mat  files) in a directory (the first argument)
% and then calculate the eigenmode data (eigenvectors, eigenvalues) and output file or files as specificed by the parameters

clearvars

% Base directory = location of this test file
thisFile = mfilename('fullpath');
[thisDir, ~, ~] = fileparts(thisFile);

packingDir = fullfile(thisDir, 'packing');
outputDir  = fullfile(thisDir, 'testOutput');

processEigenModesDampedPara( ...
    packingDir, ...
    outputDir, ...
    [0.01, 0.001], ...
    "periodic", true, ...
    "serial", false, ...
    "singleFiles", true);

matFile = fullfile(outputDir, "2D_damped_eigenstuff_N100_10by10_K100_M1.mat");

% 0) output file exists
assert(exist(matFile, 'file') == 2, ...
    'Expected eigenstuff MAT-file not found');

load(matFile, 'outData');

% 1) outData exists in the workspace
assert(exist('outData','var') == 1, ...
    'outData not found in workspace after load');

% 2) outData has an eigenValues field
assert(isfield(outData, 'eigenValues'), ...
    'outData.eigenValues field is missing');

% 3) eigenValues is a cell array with at least one entry
assert(iscell(outData.eigenValues) && ~isempty(outData.eigenValues), ...
    'outData.eigenValues must be a non-empty cell array');

% 4) the first entry is non-empty numeric eigenvalues
assert(~isempty(outData.eigenValues{1}), ...
    'First eigenvalue set (eigenValues{1}) is empty');
assert(isnumeric(outData.eigenValues{1}), ...
    'First eigenvalue set (eigenValues{1}) is not numeric');

ev1 = outData.eigenValues{1};
assert(isvector(ev1) && numel(ev1) == 400, ...
    'eigenValues{1} should contain 400 eigenvalues');

disp('TESTmeetingAug2026: ALL PASSED');

