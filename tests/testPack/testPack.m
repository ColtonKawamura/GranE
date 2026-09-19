% this tests the pack.m function
% matlab should be cd into ~/repos/ in order to run with:
% run("GranE/tests/testPack/testPack.m")

clear all

% ---- test findNeighbors3D.m first

run("testFindNeighbors3D.m") % testFindNeighbors3D.m shoudl sit next to this test
delete("~/repos/GranE/tests/testPack/data/*.mat");


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

assert(foo.scalPackingFraction > 0.6 && foo.scalPackingFraction < 0.8, ...
    sprintf('2D scalPackingFraction = %.4f, expected between 0.6 and 0.8', ...
            foo.scalPackingFraction));

assert(foo.scalMeanCoordNum > 3.6 && foo.scalMeanCoordNum < 4.4, ...
    sprintf('2D mean coordination number = %.4f, expected between 3.6 and 4.4', ...
            foo.scalMeanCoordNum));

scal2DPackFrac = foo.scalPackingFraction;
scalPressFrictionless = foo.scalPressure;
clear foo
% delete("GranE/tests/testPack/data/*.png");

% Test packRepeatTile
scalXMult = 2;
scalYMult = 2;
scalWidth = 10;
boolCalcEig = false;
stringInPath = "~/repos/GranE/tests/testPack/data/";
stringOutPath = "~/repos/GranE/tests/testPack/data/";

packRepeatTile(scalNumParts, scalSpringConstant, scalPressTarg, scalWidth, scalSeed, scalXMult, scalYMult, boolCalcEig, stringInPath , stringOutPath);
delete("data/2D_N100_P0.001_Width10_Seed1.mat");

foo = load("~/repos/GranE/tests/testPack/data/2D_N368_P0.001_Width20_Seed1.mat");

assert(foo.scalPackingFraction > 0.6 && foo.scalPackingFraction < 0.8, ...
    sprintf('2D scalPackingFraction = %.4f, expected between 0.6 and 0.8', ...
            foo.scalPackingFraction));

assert(foo.scalMeanCoordNum > 3.6 && foo.scalMeanCoordNum < 4.4, ...
    sprintf('2D mean coordination number = %.4f, expected between 3.6 and 4.4', ...
            foo.scalMeanCoordNum));

% ----------------- 3D frictionless test -----------------
scalNumParts = 6^3;
scalZMult = 1;
pack(scalNumParts, scalSpringConstant, scalDiamSmall, scalDiamBig, scalMass, scalPressTarg, scalSeed, false, scalXMult, scalYMult, scalZMult, boolCalcEig, 'data/')


foo = load("data/3D_N216_P0.001_Width6_Seed1.mat");
delete("data/3D_N216_P0.001_Width6_Seed1.mat");

assert(foo.scalPackingFraction > 0.55 && foo.scalPackingFraction < 0.65, ...
    sprintf('3D scalPackingFraction = %.4f, expected between 0.55 and 0.65', ...
            foo.scalPackingFraction));

assert(foo.scalMeanCoordNum > 5.5 && foo.scalMeanCoordNum < 6.5, ...
    sprintf('3D mean coordination number = %.4f, expected between 5.5 and 6.5', ...
            foo.scalMeanCoordNum));

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


scalPressFrictional = foo.scalPressure;

pressRelDiff = abs(scalPressFrictional - scalPressFrictionless) / scalPressFrictionless;

if pressRelDiff < 0.05
    % Pressures are close: compare φ directly
    assert(foo.scalPackingFraction <= scal2DPackFrac + 1e-3, ...
        sprintf(['At matched pressure: expected frictional scalPackingFraction <= ', ...
                 'frictionless (%.4f), got %.4f'], ...
                scal2DPackFrac, foo.scalPackingFraction));
else
    % Pressures differ: compare to literature-style φ range near jamming
    warning('Skipping φ comparison vs frictionless: frictional p = %.4g, frictionless p = %.4g (rel diff = %.3f)', ...
            scalPressFrictional, scalPressFrictionless, pressRelDiff);

    % 2D frictional random disk packings near jamming: φ typically ~0.75–0.82
    assert(foo.scalPackingFraction > 0.7 && foo.scalPackingFraction < 0.85, ...
        sprintf('Frictional scalPackingFraction = %.4f, expected between 0.7 and 0.85 near jamming', ...
                foo.scalPackingFraction));
end

assert(foo.scalMeanCoordNum > 2 && foo.scalMeanCoordNum < 4.2, ...
    sprintf('Frictional mean coordination number = %.4f, expected between 2 and 4.2', ...
            foo.scalMeanCoordNum));

clear foo
% delete("GranE/tests/testPack/data/*.png");


disp('pack.m: ALL PASSED');


