% Focused regression for 2D frictional packing control.
% Run from repository root with:
%   octave --quiet --eval "run('tests/testPack/testPackFrictionFocused.m')"

clear;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(scriptDir));
addpath(fullfile(repoRoot, 'src'));

outDir = [tempname() filesep];
mkdir(outDir);

N = 100;
K = 100;
D = 1;
G = 1.4;
M = 1;
P_target = 1e-3;
seed = 1;

opts = struct();
opts.flagFrictionOn = true;
opts.scalFricCoef = 0.5;
opts.scalTangentialK = 1/3;
opts.scalGammaNormal = 0;
opts.scalGammaTangential = 0;
opts.saveFrictionalState = false;
opts.saveFullState = false;

pack(N, K, D, G, M, P_target, seed, false, 1, 1, 0, false, [outDir filesep], opts);

matFile = fullfile(outDir, '2D_N100_P0.001_Width10_Seed1.mat');
assert(isfile(matFile), sprintf('Expected output file missing: %s', matFile));

foo = load(matFile);
assert(isfield(foo, 'scalPressure') && foo.scalPressure > 0, 'Expected positive saved pressure.');
assert(abs(foo.scalPressure - P_target) / P_target < 0.20, ...
    sprintf('Frictional pressure %.4e not within 20%% of target %.4e', foo.scalPressure, P_target));
assert(foo.scalMeanCoordNum > 2.5 && foo.scalMeanCoordNum < 4.5, ...
    sprintf('Frictional mean coordination number %.4f out of expected range', foo.scalMeanCoordNum));
assert(foo.scalPackingFraction > 0.70 && foo.scalPackingFraction < 0.90, ...
    sprintf('Frictional packing fraction %.4f out of expected range', foo.scalPackingFraction));

rmdir(outDir, 's');

disp('testPackFrictionFocused.m: PASSED');
