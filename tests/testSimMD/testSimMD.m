addpath("../../src/") % so we can load simMD.m

clear all


% OLD data file is  GranE/tests/testSimMD/2D_N100_P0.1_Width10_Seed1.mat
% NEW data file is  GranE/tests/testSimMD/2D_N100_P0.001_Width10_Seed1.mat

scalSpringConst = 100;
scalMass = 1;
scalDamping = .02;
scalFreqDrive = 1;
scalNumPart = 100;
scalPressure = 0.001;
scalWidth = 10;
scalSeed = 1;
stringInPath = "~/repos/GranE/tests/testSimMD/data/";
stringOutPath = "~/repos/GranE/tests/testSimMD/data/";

% this one should NOT detect attenuation
% assert that it doesn't detect attenuation? How?
simMD(scalSpringConst, scalMass, scalDamping, scalFreqDrive, scalNumPart, scalPressure, scalWidth, scalSeed, stringInPath, stringOutPath, struct('cleanRats', true, 'shear', false, 'maxAmpTracking', true))




scalNumPart = 100;
simMD(scalSpringConst, scalMass, scalDamping, scalFreqDrive, scalNumPart, scalPressure, scalWidth, scalSeed, stringInPath, stringOutPath, struct('cleanRats', true, 'shear', false, 'maxAmpTracking', true))

% ── Argument-validation smoke tests for the positional options struct ──
% (mirrors src/pack.m's calling convention; see simMD.m arguments block)

% 10-argument call: options should fall back to all defaults.
simMD(scalSpringConst, scalMass, scalDamping, scalFreqDrive, scalNumPart, scalPressure, scalWidth, scalSeed, stringInPath, stringOutPath)

% 11-argument call with an empty struct(): every field must be backfilled.
simMD(scalSpringConst, scalMass, scalDamping, scalFreqDrive, scalNumPart, scalPressure, scalWidth, scalSeed, stringInPath, stringOutPath, struct())

% 11-argument call with a partial struct: missing fields must be backfilled.
simMD(scalSpringConst, scalMass, scalDamping, scalFreqDrive, scalNumPart, scalPressure, scalWidth, scalSeed, stringInPath, stringOutPath, struct('shear', true))


