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
simMD(scalSpringConst, scalMass, scalDamping, scalFreqDrive, scalNumPart, scalPressure, scalWidth, scalSeed, stringInPath, stringOutPath, 'cleanRats', true, 'shear', false, 'maxAmpTracking', true)




scalNumPart = 100;
simMD(scalSpringConst, scalMass, scalDamping, scalFreqDrive, scalNumPart, scalPressure, scalWidth, scalSeed, stringInPath, stringOutPath, 'cleanRats', true, 'shear', false, 'maxAmpTracking', true)


