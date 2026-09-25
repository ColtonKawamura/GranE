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
% Reset the 2D-tile multipliers: the 2D-tile phase left scalXMult/scalYMult = 2,
% and a 3D BASE packing must tile 1x1x1, not inherit 2D multipliers (which
% would otherwise trigger the 3D repeat-tile path below).
scalXMult = 1;
scalYMult = 1;
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

% ----------------- 3D repeat-tile invariance test -----------------
% A 10x10x10 (1000-particle) base tile is tiled 9x in z -> 9000 particles.
% If the tiling lines up periodically across tile boundaries, the packing
% fraction and the coordination number are INVARIANT: tiling scales particle
% volume and box volume by the same factor (PF unchanged) and replicates the
% contact network under full PBC (Zn unchanged). This is the requested smoke
% test on a small packing (<= 9000 particles) that verifies the tiling is
% physically consistent, not merely error-free.
scalNumParts3D = 10^3;               % 1000 = 10x10x10 base tile
optsTile = struct('saveFullState', true);  % 3D repeat-tile needs the pre-cleanRats state
% z_mult = 9 -> 3D mode (z_mult ~= 0) AND 9 copies in z -> 9000 particles.
pack(scalNumParts3D, scalSpringConstant, scalDiamSmall, scalDiamBig, scalMass, ...
    scalPressTarg, scalSeed, false, 1, 1, 9, false, 'data/', optsTile)

% Load the base FULL (pre-cleanRats) state: the tile source (1000 particles).
strBaseFull = "data/3D_N1000_P0.001_Width10_Seed1_Full.mat";
assert(isfile(strBaseFull), sprintf('base full-state tile missing: %s', strBaseFull));
base = load(strBaseFull);
basePF = computePackingFraction(base.vecDiameter, base.scalBoxWidthX, base.scalBoxHeightY, base.scalBoxDepthZ);
baseZn = computeMeanCoordNum(base.vecPosX, base.vecPosY, base.vecDiameter, base.scalBoxWidthX, base.scalBoxHeightY, base.vecPosZ, base.scalBoxDepthZ);

% Load the TILED (9000-particle superlattice) and recompute PF + Zn from its
% stored positions (independent of the stored summary numbers).
strTiled = "data/3D_N9000_P0.001_Width10_Seed1_TiledX1Y1Z9.mat";
assert(isfile(strTiled), sprintf('tiled 9000-particle file missing: %s', strTiled));
til = load(strTiled);
assert(til.N == 9000, sprintf('tiled N = %d, expected 9000', til.N));
tilPF = computePackingFraction(til.vecDiameter, til.scalBoxWidthX, til.scalBoxHeightY, til.scalBoxDepthZ);
tilZn = computeMeanCoordNum(til.vecPosX, til.vecPosY, til.vecDiameter, til.scalBoxWidthX, til.scalBoxHeightY, til.vecPosZ, til.scalBoxDepthZ);

fprintf('3D tile invariance: base PF=%.6f Zn=%.4f | tiled PF=%.6f Zn=%.4f\n', basePF, baseZn, tilPF, tilZn);

% The stored summary metrics should agree with the geometry recomputation
% (guards against a save() that recorded the wrong state).
assert(abs(til.scalPackingFraction - tilPF) < 1e-6, ...
    sprintf('tiled stored PF=%.6f differs from recomputed %.6f', til.scalPackingFraction, tilPF));
assert(abs(til.scalMeanCoordNum - tilZn) < 1e-3, ...
    sprintf('tiled stored Zn=%.4f differs from recomputed %.4f', til.scalMeanCoordNum, tilZn));

% The invariance the user asked to verify: PF and Zn unchanged by tiling.
assert(abs(tilPF - basePF) < 1e-9, ...
    sprintf('PF changed under tiling: base=%.6f tiled=%.6f', basePF, tilPF));
assert(abs(tilZn - baseZn) < 1e-6, ...
    sprintf('Zn changed under tiling: base=%.4f tiled=%.4f', baseZn, tilZn));

% Clean up so the suite can re-run without stale-file skips.
delete(strBaseFull);
delete(strTiled);
delete("data/3D_N1000_P0.001_Width10_Seed1.mat");
clear base til

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


