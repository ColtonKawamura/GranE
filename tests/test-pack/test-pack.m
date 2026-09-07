% this tests the pack.m function
% matlab should be cd into ~/repos/GranE/ in order to run with:
% run("GranE/tests/test-pack/test-pack.m")

clear all


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

% ----------------- 2D frictionless test -----------------
pack(scalNumParts, scalSpringConstant, scalDiamSmall, scalDiamBig, scalMass, scalPressTarg, scalSeed, false, scalXMult, scalYMult, scalZMult, boolCalcEig, 'GranE/tests/test-pack/data/')

foo = load("GranE/tests/test-pack/data/2D_N100_P0.001_Width10_Seed1.mat");

assert(foo.scalPackingFraction > 0.6 && foo.scalPackingFraction < 0.8, ...
    sprintf('2D scalPackingFraction = %.4f, expected between 0.6 and 0.8', ...
            foo.scalPackingFraction));

% assert(coordNumber > 3.8 && coordNumber < 4.4, ...
%     sprintf('2D mean coordination number = %.4f, expected between 3.8 and 4.4', ...
%             coordNumber));

scal2DPackFrac = foo.scalPackingFraction;
clear foo
delete("GranE/tests/test-pack/data/2D_N100_P0.001_Width10_Seed1.mat");

% ----------------- 3D frictionless test -----------------
scalNumParts = 10^3;
scalZMult = 1;
pack(scalNumParts, scalSpringConstant, scalDiamSmall, scalDiamBig, scalMass, scalPressTarg, scalSeed, false, scalXMult, scalYMult, scalZMult, boolCalcEig, 'GranE/tests/test-pack/data/')


foo = load("GranE/tests/test-pack/data/3D_N1000_P0.001_Width10_Seed1.mat");

assert(foo.scalPackingFraction > 0.55 && foo.scalPackingFraction < 0.65, ...
    sprintf('3D scalPackingFraction = %.4f, expected between 0.55 and 0.65', ...
            foo.scalPackingFraction));

% assert(coordNumber > 5.5 && coordNumber < 6.5, ...
%     sprintf('3D mean coordination number = %.4f, expected between 5.5 and 6.5', ...
%             coordNumber));

delete("GranE/tests/test-pack/data/3D_N1000_P0.001_Width10_Seed1.mat");
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

pack(scalNumParts, scalSpringConstant, scalDiamSmall, scalDiamBig, scalMass, scalPressTarg, scalSeed, false, scalXMult, scalYMult, scalZMult, boolCalcEig, 'GranE/tests/test-pack/data/', opts)

foo = load("GranE/tests/test-pack/data/2D_N100_P0.001_Width10_Seed1.mat");

assert(foo.scalPackingFraction < scal2DPackFrac, ...
    sprintf('Expected frictional scalPackingFraction < frictionless (%.4f), got %.4f', ...
            scal2DPackFrac, foo.scalPackingFraction));


clear foo
delete("GranE/tests/test-pack/data/2D_N100_P0.001_Width10_Seed1.mat");

disp('pack.m: ALL PASSED');
