function runFrictionalDOS(strContactLaw, options)
% runFrictionalDOS  Build 2D frictional packings and save their vibrational spectra.
%
%   runFrictionalDOS("hooke")
%   runFrictionalDOS("hertz")
%   runFrictionalDOS("hooke", "vecMu", [0.1 1], "vecPressure", [1e-3 1e-4], "vecSeed", 1)
%
%   Reproduces the numerical experiment of Henkes, van Hecke & van Saarloos,
%   EPL 90, 14003 (2010) with GranE's pack.m and matSpringDampMass:
%
%   1. pack.m makes a frictional packing for every (mu, P, seed) of the grid
%      (Hookean + Cundall-Strack contacts for "hooke"; Hertzian normal force +
%      Mindlin tangential stiffness for "hertz"). Existing packings are reused
%      (pack.m skips a packing whose .mat exists), so an interrupted run
%      resumes where it stopped.
%   2. Rattlers are removed (particles with fewer than 2 contacts, repeated
%      until none are left) and the pack.m contact forces (fricState) are
%      re-indexed to the backbone.
%   3. matSpringDampMass builds the stiffness matrix of Henkes et al. eq. (2)
%      (pre-stress term on) in four versions:
%         stick        every contact keeps k_t ("non-slipping", Somfai et al.)
%         slip         fully mobilized contacts (m = 1) get k_t = 0 (Henkes)
%         allSlip      every contact gets k_t = 0 (the k_t -> 0 limit, Fig. 3d)
%         frictionless the same backbone without friction (2 DOF per particle)
%   4. calcEigenFrequencies gives omega for each version; the spectra and the
%      packing statistics (z, n_m, pressure, ...) are appended to
%      <strDataDir>/<strContactLaw>/spectra.mat after every packing.
%
%   plotFrictionalDOS makes the figures from spectra.mat.
%
%   Options:
%       vecMu          friction coefficients (default [0.01 0.1 0.3 1 10])
%       vecPressure    pack.m target pressures (default [1e-2 1e-3 1e-4])
%       vecSeed        random seeds (default [1 2])
%       matExtraRuns   extra [mu, P, seed] rows (default: mu = 0.3 at
%                      P = 3e-5, seeds 1 and 2, for the pressure series)
%       scalN          particles per packing (default 500)
%       scalK          pack.m spring constant K (default 100)
%       scalKtOverK    tangential/normal stiffness ratio (default 1/3)
%       scalMobTol     m >= 1 - scalMobTol is fully mobilized (default 1e-6)
%       strDataDir     output directory (default <GranE>/data/frictionalDOS)

    arguments
        strContactLaw (1, 1) string {mustBeMember(strContactLaw, ["hooke", "hertz"])}
        options.vecMu        (1, :) double = [0.01 0.1 0.3 1 10]
        options.vecPressure  (1, :) double = [1e-2 1e-3 1e-4]
        options.vecSeed      (1, :) double = [1 2]
        options.matExtraRuns (:, 3) double = [0.3 3e-5 1; 0.3 3e-5 2]
        options.scalN        (1, 1) double = 500
        options.scalK        (1, 1) double = 100
        options.scalKtOverK  (1, 1) double = 1/3
        options.scalMobTol   (1, 1) double = 1e-6
        options.strDataDir   (1, 1) string = ""
    end

    strRepoRoot = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(strRepoRoot, 'src'));
    strDataDir = options.strDataDir;
    if strDataDir == ""
        strDataDir = fullfile(strRepoRoot, 'data', 'frictionalDOS');
    end
    strLawDir = fullfile(strDataDir, strContactLaw);
    if ~isfolder(strLawDir)
        mkdir(strLawDir);
    end
    strSpectraFile = fullfile(strLawDir, 'spectra.mat');

    % The run list: full grid, then the extra rows.
    [matMuGrid, matPGrid, matSeedGrid] = ndgrid(options.vecMu, options.vecPressure, options.vecSeed);
    matRuns = [matMuGrid(:), matPGrid(:), matSeedGrid(:); options.matExtraRuns];
    matRuns = unique(matRuns, 'rows', 'stable');

    % Resume: keep the spectra already computed.
    sRuns = struct([]);
    if isfile(strSpectraFile)
        sLoaded = load(strSpectraFile, 'sRuns');
        sRuns = sLoaded.sRuns;
    end

    for idxRun = 1:size(matRuns, 1)
        scalMu = matRuns(idxRun, 1);
        scalP  = matRuns(idxRun, 2);
        scalSeed = matRuns(idxRun, 3);
        if ~isempty(sRuns) && any([sRuns.scalMu] == scalMu & [sRuns.scalPressureTarget] == scalP & [sRuns.scalSeed] == scalSeed)
            fprintf('[runFrictionalDOS] %s mu = %g, P = %g, seed %d: already done\n', strContactLaw, scalMu, scalP, scalSeed);
            continue;
        end
        fprintf('[runFrictionalDOS] %s mu = %g, P = %g, seed %d (%d of %d)\n', ...
            strContactLaw, scalMu, scalP, scalSeed, idxRun, size(matRuns, 1));
        tic;

        % 1) the packing
        strPackDir = fullfile(strLawDir, 'packings', sprintf('mu%g', scalMu));
        if ~isfolder(strPackDir)
            mkdir(strPackDir);
        end
        sPackOptions = struct('flagFrictionOn', true, 'scalFricCoef', scalMu, ...
            'scalTangentialK', options.scalKtOverK, 'scalGammaNormal', 0, 'scalGammaTangential', 0, ...
            'saveFrictionalState', true, 'saveFullState', true, ...
            'hertzian', strContactLaw == "hertz", 'hertzMindlin', strContactLaw == "hertz");
        pack(options.scalN, options.scalK, 1, 1.4, 1, scalP, scalSeed, false, 1, 1, 0, false, ...
            string([char(strPackDir) filesep]), sPackOptions);
        close all;
        strStem = fullfile(strPackDir, sprintf('2D_N%d_P%s_Width%d_Seed%d', options.scalN, num2str(scalP), ...
            round(sqrt(options.scalN)), scalSeed));
        if strContactLaw == "hertz"
            sFull = load([char(strStem) '_Full_HertzMindlin.mat']);
            sFric = load([char(strStem) '_HertzMindlin_FricState.mat']);
        else
            sFull = load([char(strStem) '_Full.mat']);
            sFric = load([char(strStem) '_FricState.mat']);
        end

        % 2) backbone and re-indexed contact forces
        [sPacking, sFricBackbone, sStats] = frictionalBackbone(sFull, sFric.fricState, scalMu, options.scalMobTol);

        % 3) the four stiffness matrices (pre-stress on, Γ = 0)
        cellCommon = {"periodic", true, "uniformMass", true, "friction", true, "stressed", true, ...
            "fricState", sFricBackbone, "contactLaw", strContactLaw, ...
            "tangentialSpringConstant", options.scalKtOverK*options.scalK, "mobilizationTol", options.scalMobTol};
        [matKStick, ~, matMass, sInfo] = matSpringDampMass(sPacking, 0, options.scalK, 1, cellCommon{:});
        matKSlip    = matSpringDampMass(sPacking, 0, options.scalK, 1, cellCommon{:}, "slipMobilized", true);
        matKAllSlip = matSpringDampMass(sPacking, 0, options.scalK, 1, cellCommon{1:8}, "contactLaw", strContactLaw, ...
            "fricState", sFricBackbone, "slipMobilized", true, "mobilizationTol", 1);
        [matKFrictionless, ~, matMassFrictionless] = matSpringDampMass(sPacking, 0, options.scalK, 1, ...
            "periodic", true, "uniformMass", true, "stressed", true, "contactLaw", strContactLaw);
        assert(sInfo.numContacts == sStats.scalNc && sInfo.numMobilized == sStats.scalNmContacts && sInfo.numNotInFricState == 0, ...
            'runFrictionalDOS:ContactMismatch', 'matSpringDampMass found a different contact network than fricState.');

        % 4) spectra and statistics
        sRun = sStats;
        sRun.strContactLaw = char(strContactLaw);
        sRun.scalMu = scalMu;
        sRun.scalPressureTarget = scalP;
        sRun.scalSeed = scalSeed;
        sRun.sStick = calcEigenFrequencies(matKStick, matMass);
        sRun.sSlip = calcEigenFrequencies(matKSlip, matMass);
        sRun.sAllSlip = calcEigenFrequencies(matKAllSlip, matMass);
        sRun.sFrictionless = calcEigenFrequencies(matKFrictionless, matMassFrictionless);
        sRun.scalElapsed = toc;
        if isempty(sRuns)
            sRuns = sRun;
        else
            sRuns(end + 1) = sRun; %#ok<AGROW>
        end
        save(strSpectraFile, 'sRuns');
        fprintf(['[runFrictionalDOS]   Ni = %d, z = %.3f, n_m = %.3f, z - (3 + n_m) = %.3f, p = %.3e, <k_n> = %.3e | ' ...
            'zero/negative modes: stick %d/%d, slip %d/%d | %.0f s\n'], ...
            sRun.scalNi, sRun.scalZ, sRun.scalNm, sRun.scalZ - 3 - sRun.scalNm, sRun.scalPressureVirial, sRun.scalMeanKn, ...
            sRun.sStick.scalNumZero, sRun.sStick.scalNumNegative, sRun.sSlip.scalNumZero, sRun.sSlip.scalNumNegative, sRun.scalElapsed);
    end
    fprintf('[runFrictionalDOS] %s done: %d packings in %s\n', strContactLaw, numel(sRuns), strSpectraFile);
end

function [sPacking, sFricBackbone, sStats] = frictionalBackbone(sFull, sFricState, scalMu, scalMobTol)
% Remove rattlers (fewer than 2 contacts, repeated) from the pre-cleanRats
% state and re-index pack.m's per-contact forces to the surviving particles.
    scalN = numel(sFull.vecPosX);
    vecI = sFricState.vecContactNN(:);
    vecJ = sFricState.vecContactMM(:);
    vecRadius = sFull.vecDiameter(:)/2;
    vecBoolKeep = true(scalN, 1);
    while true
        vecBoolLive = vecBoolKeep(vecI) & vecBoolKeep(vecJ);
        vecZi = accumarray([vecI(vecBoolLive); vecJ(vecBoolLive)], 1, [scalN 1]);
        vecBoolDrop = vecBoolKeep & vecZi < 2;
        if ~any(vecBoolDrop)
            break;
        end
        vecBoolKeep(vecBoolDrop) = false;
    end
    vecBoolLive = vecBoolKeep(vecI) & vecBoolKeep(vecJ);
    vecNewIndex = zeros(scalN, 1);
    vecNewIndex(vecBoolKeep) = 1:sum(vecBoolKeep);

    sPacking = struct('vecPosX', sFull.vecPosX(vecBoolKeep), 'vecPosY', sFull.vecPosY(vecBoolKeep), ...
        'scalBoxWidthX', sFull.scalBoxWidthX, 'scalBoxHeightY', sFull.scalBoxHeightY, ...
        'vecDiameter', sFull.vecDiameter(vecBoolKeep));
    sFricBackbone = sFricState;
    sFricBackbone.N = sum(vecBoolKeep);
    sFricBackbone.vecRadius = vecRadius(vecBoolKeep);
    sFricBackbone.vecContactNN = vecNewIndex(vecI(vecBoolLive));
    sFricBackbone.vecContactMM = vecNewIndex(vecJ(vecBoolLive));
    sFricBackbone.vecFn = sFricState.vecFn(vecBoolLive);
    sFricBackbone.vecFt = sFricState.vecFt(vecBoolLive);
    sFricBackbone.vecOverlap = sFricState.vecOverlap(vecBoolLive);
    sFricBackbone.vecUnitTanX = sFricState.vecUnitTanX(vecBoolLive);
    sFricBackbone.vecUnitTanY = sFricState.vecUnitTanY(vecBoolLive);
    if isfield(sFricBackbone, 'matDispTan')
        sFricBackbone = rmfield(sFricBackbone, 'matDispTan');   % N x N, not needed downstream
    end

    % Statistics of the backbone (2D).
    vecFn = abs(sFricBackbone.vecFn);
    vecOverlap = sFricBackbone.vecOverlap;
    vecRi = vecRadius(vecI(vecBoolLive));
    vecRj = vecRadius(vecJ(vecBoolLive));
    vecDist = vecRi + vecRj - vecOverlap;
    vecMob = abs(sFricBackbone.vecFt)./(scalMu*vecFn);
    sStats = struct();
    sStats.scalN = scalN;
    sStats.scalNi = sum(vecBoolKeep);
    sStats.scalNc = numel(vecFn);
    sStats.scalNmContacts = sum(vecMob >= 1 - scalMobTol);
    sStats.scalZ = 2*sStats.scalNc/sStats.scalNi;
    sStats.scalNm = sStats.scalNmContacts/sStats.scalNi;        % n_m, mobilized contacts per particle
    sStats.scalPhi = sum(pi*vecRadius.^2)/(sFull.scalBoxWidthX*sFull.scalBoxHeightY);   % rattlers included
    % Virial pressure p = (1/(d A)) sum_c f_n r_c, d = 2 (normal forces only:
    % r_c is along n, so the tangential forces do not contribute to the trace).
    sStats.scalPressureVirial = sum(vecFn.*vecDist)/(2*sFull.scalBoxWidthX*sFull.scalBoxHeightY);
    % Mean normal contact stiffness: K (Hooke) or 2 K sqrt(R_eff delta) (Hertz).
    if isfield(sFricState, 'strContactLaw') && ~strcmp(sFricState.strContactLaw, 'hooke')
        vecRadiusEff = vecRi.*vecRj./(vecRi + vecRj);
        sStats.scalMeanKn = mean(2*sFricState.K*sqrt(vecRadiusEff.*vecOverlap));
    else
        sStats.scalMeanKn = sFricState.K;
    end
    sStats.vecMobilization = vecMob;
    % Particles whose contacts all slip spin freely when slipping is on.
    vecZiB = accumarray([sFricBackbone.vecContactNN; sFricBackbone.vecContactMM], 1, [sStats.scalNi 1]);
    vecBoolSlip = vecMob >= 1 - scalMobTol;
    vecMiB = accumarray([sFricBackbone.vecContactNN(vecBoolSlip); sFricBackbone.vecContactMM(vecBoolSlip)], 1, [sStats.scalNi 1]);
    sStats.scalNumFreeSpin = sum(vecMiB == vecZiB);
end
