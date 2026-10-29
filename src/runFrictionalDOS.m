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
%      Mindlin tangential stiffness for "hertz") and saves ONE file,
%      2D_N<N>_P<P>_Width<W>_Seed<s>_Friction.mat ("hertz": ..._Hertz_Friction.mat).
%      Existing packings are reused (pack.m skips a packing whose .mat
%      exists), so an interrupted run resumes where it stopped.
%      A packing that does not reach force and torque balance within
%      pack.m's step cap is not saved (pack:FrictionNotConverged). It is
%      replaced by the packing of seed + 1000 (then + 2000, ...), up to
%      scalSeedRetries times; sRuns records the requested seed (scalSeed,
%      the grid label) and the seed actually packed (scalSeedUsed).
%   2. That file already holds the backbone (rattlers, particles with fewer
%      than 2 contacts, removed repeatedly by cleanRats), the pack.m contact
%      forces re-indexed to the backbone (fricState) and the full
%      pre-cleanRats state (sFullState, for the rattlers-included numbers).
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
%       scalSeedRetries replacement seeds tried when a packing does not
%                      converge (seed + 1000*k, k = 1..scalSeedRetries; default 2)
%       strDataDir     output directory (default <GranE>/data/frictionalDOS)
%       scalForceTol   pack.m's frictional force/torque balance tolerance, in
%                      units of <|f_n|> (default 1e-4). The file names do not
%                      record it: use a separate strDataDir for another value.
%       scalMaxSteps   pack.m's frictional step cap (default 3e6)

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
        options.scalSeedRetries (1, 1) double {mustBeInteger, mustBeNonnegative} = 2
        options.strDataDir   (1, 1) string = ""
        options.scalForceTol (1, 1) double {mustBePositive} = 1e-4
        options.scalMaxSteps (1, 1) double {mustBePositive} = 3e6
    end

    strRepoRoot = fileparts(fileparts(mfilename('fullpath'))); % climb twice up to add src too path
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
    [matMuGrid, matPGrid, matSeedGrid] = ndgrid(options.vecMu, options.vecPressure, options.vecSeed); % makes every combination of three vectors
    matRuns = [matMuGrid(:), matPGrid(:), matSeedGrid(:); options.matExtraRuns];
    matRuns = unique(matRuns, 'rows', 'stable');

    % Resume: keep the spectra already computed.
    % check to see if the spectra.mat file already exsis, and if so,
    % just append the new runs to it
    sRuns = struct([]);
    if isfile(strSpectraFile)
        sLoaded = load(strSpectraFile, 'sRuns');
        sRuns = sLoaded.sRuns;
        if ~isempty(sRuns) && ~isfield(sRuns, 'scalSeedUsed')
            % spectra.mat written before replacement seeds existed: every
            % packing used its requested seed
            [sRuns.scalSeedUsed] = sRuns.scalSeed;
        end
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
            'hertzian', strContactLaw == "hertz", 'hertzMindlin', strContactLaw == "hertz", ...
            'scalFrictionForceTol', options.scalForceTol, 'scalFrictionMaxSteps', options.scalMaxSteps);
        % pack.m refuses to save a packing that hits its step cap
        % (pack:FrictionNotConverged). Replace it with the packing of another
        % seed rather than raising the cap: seed + 1000*k cannot collide with
        % the small grid seeds, and the file name records the seed used.
        for idxRetry = 0:options.scalSeedRetries
            scalSeedUsed = scalSeed + 1000*idxRetry;
            try
                pack(options.scalN, options.scalK, 1, 1.4, 1, scalP, scalSeedUsed, false, 1, 1, 0, false, ...
                    string([char(strPackDir) filesep]), sPackOptions);
                break;
            catch sPackME
                close all;
                if ~strcmp(sPackME.identifier, 'pack:FrictionNotConverged') || idxRetry == options.scalSeedRetries
                    rethrow(sPackME);
                end
                fprintf('[runFrictionalDOS]   seed %d did not converge; replacing it with seed %d\n', ...
                    scalSeedUsed, scalSeed + 1000*(idxRetry + 1));
            end
        end
        close all;
        strFileSuffix = '_Friction';
        if strContactLaw == "hertz"
            strFileSuffix = '_Hertz_Friction';
        end
        strPackingFile = fullfile(strPackDir, sprintf('2D_N%d_P%s_Width%d_Seed%d%s.mat', options.scalN, num2str(scalP), ...
            round(sqrt(options.scalN)), scalSeedUsed, strFileSuffix));
        sPackingFile = load(strPackingFile);

        % 2) backbone, its re-indexed contact forces and the statistics
        [sPacking, sFricBackbone, sStats] = frictionalBackbone(sPackingFile, scalMu, options.scalMobTol);

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
        sRun.scalSeed = scalSeed;            % requested seed (grid label)
        sRun.scalSeedUsed = scalSeedUsed;    % seed actually packed (differs after a replacement)
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

function [sPacking, sFricBackbone, sStats] = frictionalBackbone(sPackingFile, scalMu, scalMobTol)
% Backbone geometry and contact forces from pack.m's one frictional file
% (rattlers already removed by cleanRats and fricState already re-indexed
% onto the backbone by pack.m), plus the packing statistics. The full
% pre-cleanRats state (sFullState) gives the rattlers-included numbers.
    sFull = sPackingFile.sFullState;           % pre-cleanRats state (all particles)
    scalN = numel(sFull.vecPosX);
    vecRadius = sFull.vecDiameter(:)/2;        % [N x 1] all particles, rattlers included

    sPacking = struct('vecPosX', sPackingFile.vecPosX(:), 'vecPosY', sPackingFile.vecPosY(:), ...
        'scalBoxWidthX', sPackingFile.scalBoxWidthX, 'scalBoxHeightY', sPackingFile.scalBoxHeightY, ...
        'vecDiameter', sPackingFile.vecDiameter(:));
    sFricBackbone = sPackingFile.fricState;    % backbone numbering, no matDispTan
    sFricState = sFricBackbone;                % material constants (K, strContactLaw) for the stats below

    % Statistics of the backbone (2D).
    vecFn = abs(sFricBackbone.vecFn);
    vecOverlap = sFricBackbone.vecOverlap;
    vecRadiusBackbone = sPacking.vecDiameter/2;
    vecRi = vecRadiusBackbone(sFricBackbone.vecContactNN);
    vecRj = vecRadiusBackbone(sFricBackbone.vecContactMM);
    vecDist = vecRi + vecRj - vecOverlap;
    vecMob = abs(sFricBackbone.vecFt)./(scalMu*vecFn);
    sStats = struct();
    sStats.scalN = scalN;
    sStats.scalNi = numel(sPacking.vecPosX);
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
