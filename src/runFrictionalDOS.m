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
%   5. With boolUnresolvedSpectra (default true) the four spectra are
%      computed a second time without the unresolved particles: particles
%      whose backbone contacts all carry f_n < scalForceFloor*<f>, i.e.
%      forces inside pack.m's equilibration noise (Agnolin & Roux, PRE 76,
%      061302 (2007)), followed by cleanRats' z_i < 2 cascade (see
%      unresolvedParticles). They are stored next to the originals as
%      sStickCut, sSlipCut, sAllSlipCut, sFrictionlessCut with the reduced
%      network's statistics (scalNiCut, scalZCut, ...). Both versions are
%      always reported, because for Hertz contacts these few particles bias
%      the low-frequency DOS and omega*. Runs already in spectra.mat that do
%      not have the second set (or have it for another scalForceFloor) get
%      it from their saved packing, without re-packing. The packings
%      themselves are unchanged: cleanRats' own rule stays opt-in.
%
%   plotFrictionalDOS makes the figures from spectra.mat; plotUnresolvedDOS
%   (with runUnresolvedDOS) makes the with/without comparison figures.
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
%       boolUnresolvedSpectra  also store the spectra without the unresolved
%                      particles (default true; step 5)
%       scalForceFloor force floor of the unresolved-particle rule, in units
%                      of <f> = mean |f_n| of the full state (default 1e-4,
%                      pack.m's balance tolerance)

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
        options.boolUnresolvedSpectra (1, 1) logical = true
        options.scalForceFloor (1, 1) double {mustBePositive} = 1e-4
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

    % Resume, part 2: runs saved without the spectra minus the unresolved
    % particles (or with another force floor) get them from their packing.
    if options.boolUnresolvedSpectra
        for idxOld = 1:numel(sRuns)
            if isfield(sRuns, 'scalNiCut') && ~isempty(sRuns(idxOld).scalNiCut) ...
                    && isequal(sRuns(idxOld).scalForceFloor, options.scalForceFloor)
                continue;
            end
            sOld = sRuns(idxOld);
            sPackingFile = load(packingFileName(strLawDir, strContactLaw, sOld.scalN, sOld.scalMu, ...
                sOld.scalPressureTarget, sOld.scalSeedUsed));
            [sPacking, sFricBackbone] = frictionalBackbone(sPackingFile, sOld.scalMu, options.scalMobTol);
            sCut = unresolvedSpectra(sPackingFile, sPacking, sFricBackbone, sOld, strContactLaw, options);
            for strField = string(fieldnames(sCut)).'
                sRuns(idxOld).(strField) = sCut.(strField);
            end
            save(strSpectraFile, 'sRuns');
            fprintf('[runFrictionalDOS] %s mu = %g, P = %g, seed %d: added the spectra without %d unresolved particles\n', ...
                strContactLaw, sOld.scalMu, sOld.scalPressureTarget, sOld.scalSeed, sCut.scalNumUnresolvedRemoved);
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
        sPackingFile = load(packingFileName(strLawDir, strContactLaw, options.scalN, scalMu, scalP, scalSeedUsed));

        % 2) backbone, its re-indexed contact forces and the statistics
        [sPacking, sFricBackbone, sStats] = frictionalBackbone(sPackingFile, scalMu, options.scalMobTol);

        % 3) the four stiffness matrices (pre-stress on, Γ = 0) and their spectra
        [sSpectra, sInfo] = fourSpectra(sPacking, sFricBackbone, strContactLaw, options);
        assert(sInfo.numContacts == sStats.scalNc && sInfo.numMobilized == sStats.scalNmContacts && sInfo.numNotInFricState == 0, ...
            'runFrictionalDOS:ContactMismatch', 'matSpringDampMass found a different contact network than fricState.');

        % 4) spectra and statistics
        sRun = sStats;
        sRun.strContactLaw = char(strContactLaw);
        sRun.scalMu = scalMu;
        sRun.scalPressureTarget = scalP;
        sRun.scalSeed = scalSeed;            % requested seed (grid label)
        sRun.scalSeedUsed = scalSeedUsed;    % seed actually packed (differs after a replacement)
        sRun.sStick = sSpectra.sStick;
        sRun.sSlip = sSpectra.sSlip;
        sRun.sAllSlip = sSpectra.sAllSlip;
        sRun.sFrictionless = sSpectra.sFrictionless;
        sRun.scalElapsed = toc;

        % 5) the same spectra without the unresolved particles
        strCutNote = '';
        if options.boolUnresolvedSpectra
            sCut = unresolvedSpectra(sPackingFile, sPacking, sFricBackbone, sRun, strContactLaw, options);
            for strField = string(fieldnames(sCut)).'
                sRun.(strField) = sCut.(strField);
            end
            strCutNote = sprintf(' | unresolved removed %d', sCut.scalNumUnresolvedRemoved);
        end
        sRuns = appendRun(sRuns, sRun);
        save(strSpectraFile, 'sRuns');
        fprintf(['[runFrictionalDOS]   Ni = %d, z = %.3f, n_m = %.3f, z - (3 + n_m) = %.3f, p = %.3e, <k_n> = %.3e | ' ...
            'zero/negative modes: stick %d/%d, slip %d/%d%s | %.0f s\n'], ...
            sRun.scalNi, sRun.scalZ, sRun.scalNm, sRun.scalZ - 3 - sRun.scalNm, sRun.scalPressureVirial, sRun.scalMeanKn, ...
            sRun.sStick.scalNumZero, sRun.sStick.scalNumNegative, sRun.sSlip.scalNumZero, sRun.sSlip.scalNumNegative, ...
            strCutNote, sRun.scalElapsed);
    end
    fprintf('[runFrictionalDOS] %s done: %d packings in %s\n', strContactLaw, numel(sRuns), strSpectraFile);
end

function strPackingFile = packingFileName(strLawDir, strContactLaw, scalN, scalMu, scalP, scalSeedUsed)
% The one file pack.m writes for a frictional packing.
    strFileSuffix = '_Friction';
    if strContactLaw == "hertz"
        strFileSuffix = '_Hertz_Friction';
    end
    strPackingFile = fullfile(strLawDir, 'packings', sprintf('mu%g', scalMu), ...
        sprintf('2D_N%d_P%s_Width%d_Seed%d%s.mat', scalN, num2str(scalP), round(sqrt(scalN)), scalSeedUsed, strFileSuffix));
end

function [sSpectra, sInfo] = fourSpectra(sPacking, sFric, strContactLaw, options)
% Frictional (stick), Partial Slip (slip), Total Slip (allSlip) and
% post-pack frictionless spectra of one backbone (pre-stress on, Γ = 0).
    cellCommon = {"periodic", true, "uniformMass", true, "friction", true, "stressed", true, ...
        "fricState", sFric, "contactLaw", strContactLaw, ...
        "tangentialSpringConstant", options.scalKtOverK*options.scalK, "mobilizationTol", options.scalMobTol};
    [matKStick, ~, matMass, sInfo] = matSpringDampMass(sPacking, 0, options.scalK, 1, cellCommon{:});
    matKSlip    = matSpringDampMass(sPacking, 0, options.scalK, 1, cellCommon{:}, "slipMobilized", true);
    matKAllSlip = matSpringDampMass(sPacking, 0, options.scalK, 1, cellCommon{1:8}, "contactLaw", strContactLaw, ...
        "fricState", sFric, "slipMobilized", true, "mobilizationTol", 1);
    [matKFrictionless, ~, matMassFrictionless] = matSpringDampMass(sPacking, 0, options.scalK, 1, ...
        "periodic", true, "uniformMass", true, "stressed", true, "contactLaw", strContactLaw);
    sSpectra = struct();
    sSpectra.sStick = calcEigenFrequencies(matKStick, matMass);
    sSpectra.sSlip = calcEigenFrequencies(matKSlip, matMass);
    sSpectra.sAllSlip = calcEigenFrequencies(matKAllSlip, matMass);
    sSpectra.sFrictionless = calcEigenFrequencies(matKFrictionless, matMassFrictionless);
end

function sCut = unresolvedSpectra(sPackingFile, sPacking, sFric, sRun, strContactLaw, options)
% The four spectra and the network statistics without the unresolved
% particles (unresolvedParticles' "unresolved" rule + cascade). sRun holds
% the spectra with every particle; they are reused when nothing is removed.
    scalNi = numel(sPacking.vecPosX);
    scalForceMean = mean(abs(sPackingFile.sFullState.fricState.vecFn));     % <f>: all contacts, full state
    [vecBoolKeep, vecBoolHit] = unresolvedParticles(sFric, scalNi, scalForceMean, "unresolved", options.scalForceFloor);
    sCut = struct();
    sCut.scalForceFloor = options.scalForceFloor;
    sCut.scalNumUnresolvedHit = nnz(vecBoolHit);             % removed by the rule itself
    sCut.scalNumUnresolvedRemoved = nnz(~vecBoolKeep);       % rule + cascade
    if all(vecBoolKeep)
        sNet = networkStats(sPacking, sFric, sRun.scalMu, options.scalMobTol);
        sSpectra = struct('sStick', sRun.sStick, 'sSlip', sRun.sSlip, 'sAllSlip', sRun.sAllSlip, ...
            'sFrictionless', sRun.sFrictionless);
    elseif nnz(vecBoolKeep) < 2
        error('runFrictionalDOS:UnresolvedCollapse', ...
            'Removing the unresolved particles leaves fewer than 2 particles (mu = %g, P = %g, seed %d).', ...
            sRun.scalMu, sRun.scalPressureTarget, sRun.scalSeed);
    else
        [sPackCut, sFricCut] = subsetBackbone(sPacking, sFric, vecBoolKeep);
        sNet = networkStats(sPackCut, sFricCut, sRun.scalMu, options.scalMobTol);
        [sSpectra, sInfo] = fourSpectra(sPackCut, sFricCut, strContactLaw, options);
        assert(sInfo.numContacts == sNet.scalNc && sInfo.numNotInFricState == 0, ...
            'runFrictionalDOS:ContactMismatch', 'matSpringDampMass found a different contact network than the reduced fricState.');
    end
    sCut.scalNiCut = sNet.scalNi;
    sCut.scalNcCut = sNet.scalNc;
    sCut.scalZCut = sNet.scalZ;
    sCut.scalNmCut = sNet.scalNm;
    sCut.scalMeanKnCut = sNet.scalMeanKn;
    sCut.scalNumFreeSpinCut = sNet.scalNumFreeSpin;
    sCut.sStickCut = sSpectra.sStick;
    sCut.sSlipCut = sSpectra.sSlip;
    sCut.sAllSlipCut = sSpectra.sAllSlip;
    sCut.sFrictionlessCut = sSpectra.sFrictionless;
end

function sRuns = appendRun(sRuns, sRun)
% Append one run to the struct array; a field only one side has is added to
% the other as [] (e.g. a run made with boolUnresolvedSpectra = false).
    if isempty(sRuns)
        sRuns = sRun;
        return;
    end
    for strField = string(setdiff(fieldnames(sRuns), fieldnames(sRun))).'
        sRun.(strField) = [];
    end
    for strField = string(setdiff(fieldnames(sRun), fieldnames(sRuns))).'
        [sRuns.(strField)] = deal([]);
    end
    sRuns(end + 1) = orderfields(sRun, sRuns(1));
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

    % Statistics of the backbone (2D).
    sNet = networkStats(sPacking, sFricBackbone, scalMu, scalMobTol);
    vecFn = abs(sFricBackbone.vecFn);
    vecRadiusBackbone = sPacking.vecDiameter/2;
    vecDist = vecRadiusBackbone(sFricBackbone.vecContactNN) + vecRadiusBackbone(sFricBackbone.vecContactMM) ...
        - sFricBackbone.vecOverlap;
    sStats = struct();
    sStats.scalN = scalN;
    sStats.scalNi = sNet.scalNi;
    sStats.scalNc = sNet.scalNc;
    sStats.scalNmContacts = sNet.scalNmContacts;
    sStats.scalZ = sNet.scalZ;
    sStats.scalNm = sNet.scalNm;                                % n_m, mobilized contacts per particle
    sStats.scalPhi = sum(pi*vecRadius.^2)/(sFull.scalBoxWidthX*sFull.scalBoxHeightY);   % rattlers included
    % Virial pressure p = (1/(d A)) sum_c f_n r_c, d = 2 (normal forces only:
    % r_c is along n, so the tangential forces do not contribute to the trace).
    sStats.scalPressureVirial = sum(vecFn.*vecDist)/(2*sFull.scalBoxWidthX*sFull.scalBoxHeightY);
    sStats.scalMeanKn = sNet.scalMeanKn;
    sStats.vecMobilization = sNet.vecMobilization;
    sStats.scalNumFreeSpin = sNet.scalNumFreeSpin;
end

function sNet = networkStats(sPacking, sFric, scalMu, scalMobTol)
% Contact-network statistics of a backbone and its re-indexed fricState:
% particle and contact counts, z, n_m, <k_n> and the free spinners.
    vecFn = abs(sFric.vecFn);
    vecOverlap = sFric.vecOverlap;
    vecRadiusBackbone = sPacking.vecDiameter/2;
    vecRi = vecRadiusBackbone(sFric.vecContactNN);
    vecRj = vecRadiusBackbone(sFric.vecContactMM);
    vecMob = abs(sFric.vecFt)./(scalMu*vecFn);
    sNet = struct();
    sNet.scalNi = numel(sPacking.vecPosX);
    sNet.scalNc = numel(vecFn);
    sNet.scalNmContacts = sum(vecMob >= 1 - scalMobTol);
    sNet.scalZ = 2*sNet.scalNc/sNet.scalNi;
    sNet.scalNm = sNet.scalNmContacts/sNet.scalNi;
    % Mean normal contact stiffness: K (Hooke) or 2 K sqrt(R_eff delta) (Hertz).
    if isfield(sFric, 'strContactLaw') && ~strcmp(sFric.strContactLaw, 'hooke')
        vecRadiusEff = vecRi.*vecRj./(vecRi + vecRj);
        sNet.scalMeanKn = mean(2*sFric.K*sqrt(vecRadiusEff.*vecOverlap));
    else
        sNet.scalMeanKn = sFric.K;
    end
    sNet.vecMobilization = vecMob;
    % Particles whose contacts all slip spin freely when slipping is on.
    vecZiB = accumarray([sFric.vecContactNN; sFric.vecContactMM], 1, [sNet.scalNi 1]);
    vecBoolSlip = vecMob >= 1 - scalMobTol;
    vecMiB = accumarray([sFric.vecContactNN(vecBoolSlip); sFric.vecContactMM(vecBoolSlip)], 1, [sNet.scalNi 1]);
    sNet.scalNumFreeSpin = sum(vecMiB == vecZiB);
end
