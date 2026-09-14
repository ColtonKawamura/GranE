% this tests the pack.m function
% matlab should be cd into ~/repos/ in order to run with:
% run("GranE/tests/testPack/testPack.m")

clear all

% ---- test findNeighbors3D.m first

run("testFindNeighbors3D.m") % testFindNeighbors3D.m shoudl sit next to this test


% ----------------- 2D frictionless test -----------------
scalNumParts = 100;
scalSpringConstant = 100;
scalDiamSmall = 1;
scalDiamBig = 1.4;
scalMass = 1;
scalPressTarg = .001;
scalSeed = 1;
scalXMult = 1;
scalYMult = 1;
scalZMult = 0;
boolCalcEig = false;

%% without friction
pack(scalNumParts, scalSpringConstant, scalDiamSmall, scalDiamBig, scalMass, scalPressTarg, scalSeed, false, scalXMult, scalYMult, scalZMult, boolCalcEig, 'data/')

foo = load("data/2D_N100_P0.001_Width10_Seed1.mat");
delete("data/2D_N100_P0.001_Width10_Seed1.mat");

assert(foo.scalPackingFraction > 0.6 && foo.scalPackingFraction < 0.8, ...
    sprintf('2D scalPackingFraction = %.4f, expected between 0.6 and 0.8', ...
            foo.scalPackingFraction));

% assert(foo.scalMeanCoordNum > 3.8 && foo.scalMeanCoordNum < 4.4, ...
%     sprintf('2D mean coordination number = %.4f, expected between 3.8 and 4.4', ...
%             foo.scalMeanCoordNum));

scal2DPackFrac = foo.scalPackingFraction;
clear foo
% delete("GranE/tests/testPack/data/*.png");

% ----------------- 3D frictionless test -----------------
scalNumParts = 6^3;
scalZMult = 1;
pack(scalNumParts, scalSpringConstant, scalDiamSmall, scalDiamBig, scalMass, scalPressTarg, scalSeed, false, scalXMult, scalYMult, scalZMult, boolCalcEig, 'data/')


foo = load("data/3D_N216_P0.001_Width6_Seed1.mat");
delete("data/3D_N216_P0.001_Width6_Seed1.mat");

assert(foo.scalPackingFraction > 0.55 && foo.scalPackingFraction < 0.65, ...
    sprintf('3D scalPackingFraction = %.4f, expected between 0.55 and 0.65', ...
            foo.scalPackingFraction));

% assert(foo.scalMeanCoordNum > 5.5 && foo.scalMeanCoordNum < 6.5, ...
%     sprintf('3D mean coordination number = %.4f, expected between 5.5 and 6.5', ...
%             foo.scalMeanCoordNum));

% delete("GranE/tests/testPack/data/*.png");
clear foo

%% with friction

opts.flagFrictionOn        = true;   % master switch (2D only)
opts.scalFricCoef         = 0.5;    % mu, Coulomb coeff (default 0.50)
opts.scalTangentialK      = 1/3;    % Kt/K ratio, Cundall-Strack default
opts.scalGammaNormal      = 0;      % normal dashpot (0 keeps original)
opts.scalGammaTangential  = 0;      % optional tangential dashpot
opts.saveFrictionalState  = true;   % export fricState sidecar for K^fric

% ----------------- 2D with friction test -----------------
scalNumParts = 10^2;
scalZMult = 0;

pack(scalNumParts, scalSpringConstant, scalDiamSmall, scalDiamBig, scalMass, scalPressTarg, scalSeed, false, scalXMult, scalYMult, scalZMult, boolCalcEig, 'data/', opts)

foo = load("data/2D_N100_P0.001_Width10_Seed1.mat");
delete("data/2D_N100_P0.001_Width10_Seed1.mat");

assert(foo.scalPackingFraction < scal2DPackFrac, ...
    sprintf('Expected frictional scalPackingFraction < frictionless (%.4f), got %.4f', ...
            scal2DPackFrac, foo.scalPackingFraction));

% assert(foo.scalMeanCoordNum > 2 && foo.scalMeanCoordNum < 4, ...
%     sprintf('Frictional mean coordination number = %.4f, expected between 5.5 and 6.5', ...
%             foo.scalMeanCoordNum));

clear foo
% delete("GranE/tests/testPack/data/*.png");

disp('pack.m: ALL PASSED');
