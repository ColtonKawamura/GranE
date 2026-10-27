function sResults = runUnresolvedDOS(strContactLaw, options)
% runUnresolvedDOS  Vibrational spectra with and without unresolved particles.
%
%   sResults = runUnresolvedDOS("hertz")
%   sResults = runUnresolvedDOS("hooke", "strRule", "twoContact")
%   sResults = runUnresolvedDOS("hertz", "strDataDir", "data/frictionalDOS_N2000")
%
%   For every packing of a runFrictionalDOS data set, removes a set of
%   particles from the backbone, recomputes the Frictional (sStick, every
%   contact keeps k_t) and Partial Slip (sSlip, fully mobilized contacts get
%   k_t = 0) spectra with runFrictionalDOS' own matSpringDampMass and
%   calcEigenFrequencies calls, and saves both versions side by side.
%
%   Two removal rules (strRule):
%     "unresolved"  Agnolin & Roux, PRE 76, 061302 (2007): a particle whose
%                   backbone contacts ALL carry f_n < scalForceFloor*<f> is
%                   held only by forces inside pack.m's equilibration noise
%                   (its force balance tolerance is 1e-4 <|f_n|>). It is
%                   removed, then cleanRats' z_i < 2 rule is re-applied
%                   until nothing changes (the cascade). This is the rule
%                   cleanRats applies with sUncertOptions.boolRemoveUnresolved.
%     "twoContact"  Shimada et al., arXiv:2608.30151 (2026) and Ishima et al.,
%                   arXiv:2207.06632 (2023): every particle with z_i <= 2 is
%                   removed, iterated. Loaded two-contact particles go too.
%                   The network can dismantle completely; such packings are
%                   flagged boolCollapsed and get no spectrum.
%   <f> is the mean |f_n| over all contacts of the full pre-cleanRats state,
%   the normalization of cleanRats and plotUncertContacts.
%
%   Also stored, for plotUnresolvedDOS:
%     - the modes of the ORIGINAL Frictional dynamical matrix: their
%       frequency, participation ratio and the share of each mode on the
%       removed particles (boolModes, dense eig: about 1 min at N = 2000),
%     - every backbone contact's f_n/<f> and k_n/<k_n> (the mechanism: Hertz
%       k_n ~ f_n^(1/3) makes noise-level contacts soft, Hooke k_n = K does not),
%     - tan(beta) of every two-contact particle, kept or removed (the Coulomb
%       cone test of Gilabert et al. 2007: a two-contact disk can carry load
%       only if tan(beta) <= mu, beta = (pi - theta)/2, theta the angle between
%       its two contact normals).
%
%   Frequencies are omega_hat = omega/sqrt(<k_n>/m) with the <k_n> of the
%   network they come from (as in plotFrictionalDOS).
%
%   Writes <strDataDir>/<strContactLaw>/unresolved_<strRule>.mat (sResults,
%   sOptions). spectra.mat and the packings are only read.
%
%   Options:
%       strDataDir      runFrictionalDOS output (default <GranE>/data/frictionalDOS)
%       strRule         "unresolved" (default) or "twoContact"
%       scalForceFloor  noise floor in units of <f> (default 1e-4)
%       scalMobTol      m >= 1 - scalMobTol is fully mobilized (default 1e-6,
%                       as runFrictionalDOS)
%       boolModes       store the mode localization (default true)
%       scalCollapseFraction  a packing that keeps fewer than this fraction of
%                       its backbone is flagged collapsed (default 0.5)

    arguments
        strContactLaw (1, 1) string {mustBeMember(strContactLaw, ["hooke", "hertz"])}
        options.strDataDir     (1, 1) string = ""
        options.strRule        (1, 1) string {mustBeMember(options.strRule, ["unresolved", "twoContact"])} = "unresolved"
        options.scalForceFloor (1, 1) double {mustBePositive} = 1e-4
        options.scalMobTol     (1, 1) double {mustBePositive} = 1e-6
        options.boolModes      (1, 1) logical = true
        options.scalCollapseFraction (1, 1) double = 0.5
    end

    strRepoRoot = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(strRepoRoot, 'src'));
    strDataDir = options.strDataDir;
    if strDataDir == ""
        strDataDir = fullfile(strRepoRoot, 'data', 'frictionalDOS');
    end
    strLawDir = fullfile(strDataDir, strContactLaw);
    sSpectra = load(fullfile(strLawDir, 'spectra.mat'), 'sRuns');
    sRuns = sSpectra.sRuns;
    if ~isfield(sRuns, 'scalSeedUsed')
        [sRuns.scalSeedUsed] = sRuns.scalSeed;
    end
    strSuffix = '_Friction';
    if strContactLaw == "hertz"
        strSuffix = '_Hertz_Friction';
    end
    strOutFile = fullfile(strLawDir, sprintf('unresolved_%s.mat', options.strRule));
    sOptions = options;

    sResults = struct([]);
    for idxRun = 1:numel(sRuns)
        sRun = sRuns(idxRun);
        strFile = fullfile(strLawDir, 'packings', sprintf('mu%g', sRun.scalMu), ...
            sprintf('2D_N%d_P%s_Width%d_Seed%d%s.mat', sRun.scalN, num2str(sRun.scalPressureTarget), ...
            round(sqrt(sRun.scalN)), sRun.scalSeedUsed, strSuffix));
        sFile = load(strFile);
        sPack = struct('vecPosX', sFile.vecPosX(:), 'vecPosY', sFile.vecPosY(:), ...
            'scalBoxWidthX', sFile.scalBoxWidthX, 'scalBoxHeightY', sFile.scalBoxHeightY, ...
            'vecDiameter', sFile.vecDiameter(:));
        sFric = sFile.fricState;
        scalNi = numel(sPack.vecPosX);
        vecI = sFric.vecContactNN(:);
        vecJ = sFric.vecContactMM(:);
        scalFMean = mean(abs(sFile.sFullState.fricState.vecFn));      % <f>: all contacts, full state
        vecForceRel = abs(sFric.vecFn(:))/scalFMean;                  % [Nc x 1]
        vecKn = contactStiffness(sPack, sFric, strContactLaw);        % [Nc x 1]
        vecZ = accumarray([vecI; vecJ], 1, [scalNi 1]);

        % 1) the particles the rule removes, then the cascade
        if options.strRule == "unresolved"
            vecBoolWeak = vecForceRel < options.scalForceFloor;
            vecZWeak = accumarray([vecI(vecBoolWeak); vecJ(vecBoolWeak)], 1, [scalNi 1]);
            vecBoolHit = vecZ > 0 & vecZWeak == vecZ;
            scalZMin = 2;                                             % cleanRats' rattler rule
        else
            vecBoolHit = vecZ <= 2;
            scalZMin = 3;                                             % z_i <= 2 removed, iterated
        end
        vecBoolKeep = cascade(~vecBoolHit, vecI, vecJ, scalNi, scalZMin);

        s = struct();
        s.strContactLaw = char(strContactLaw);
        s.strRule = char(options.strRule);
        s.scalMu = sRun.scalMu;
        s.scalPressureTarget = sRun.scalPressureTarget;
        s.scalSeed = sRun.scalSeed;
        s.scalSeedUsed = sRun.scalSeedUsed;
        s.scalN = sRun.scalN;
        s.scalNi = scalNi;
        s.scalNumHit = nnz(vecBoolHit);
        s.scalNumRemoved = nnz(~vecBoolKeep);
        s.scalNumCascade = nnz(~vecBoolKeep & ~vecBoolHit);
        s.boolCollapsed = nnz(vecBoolKeep) < options.scalCollapseFraction*scalNi;
        s.scalZ = sRun.scalZ;
        s.scalNm = sRun.scalNm;
        s.scalMeanKn = sRun.scalMeanKn;
        s.sStick = sRun.sStick;
        s.sSlip = sRun.sSlip;

        % 2) contact-level data for the mechanism figure
        vecBoolOnRemoved = ~vecBoolKeep(vecI) | ~vecBoolKeep(vecJ);
        s.vecForceRel = vecForceRel;
        s.vecKnRel = vecKn/mean(vecKn);
        s.vecBoolContactRemoved = vecBoolOnRemoved;
        vecLoad = accumarray([vecI; vecJ], [vecForceRel; vecForceRel], [scalNi 1]);
        s.vecLoadRemoved = vecLoad(~vecBoolKeep);                     % sum f_n/<f> on each removed particle
        [vecTanBeta, vecIdxTwo] = twoContactCone(sPack, sFric, vecZ);
        s.vecTanBetaKept = vecTanBeta(vecBoolKeep(vecIdxTwo));
        s.vecTanBetaRemoved = vecTanBeta(~vecBoolKeep(vecIdxTwo));

        % 3) the modes of the original Frictional matrix, localized or not
        cellCommon = {"periodic", true, "uniformMass", true, "friction", true, "stressed", true, ...
            "contactLaw", strContactLaw, "tangentialSpringConstant", sFric.K/3, "mobilizationTol", options.scalMobTol};
        if options.boolModes
            [matK, ~, matMass] = matSpringDampMass(sPack, 0, sFric.K, 1, cellCommon{:}, "fricState", sFric);
            [s.vecModeOmegaHat, s.vecModeShareRemoved, s.vecModePR] = modeLocalization(matK, matMass, ~vecBoolKeep, s.scalMeanKn);
        end

        % 4) the spectra without the removed particles
        s.scalNumKept = nnz(vecBoolKeep);
        if s.boolCollapsed
            s.scalZCut = NaN; s.scalNmCut = NaN; s.scalMeanKnCut = NaN;
            s.sStickCut = []; s.sSlipCut = [];
        elseif all(vecBoolKeep)
            s.scalZCut = s.scalZ; s.scalNmCut = s.scalNm; s.scalMeanKnCut = s.scalMeanKn;
            s.sStickCut = s.sStick; s.sSlipCut = s.sSlip;
        else
            [sPackCut, sFricCut] = subsetBackbone(sPack, sFric, vecBoolKeep);
            [matKStick, ~, matMassCut, sInfo] = matSpringDampMass(sPackCut, 0, sFric.K, 1, cellCommon{:}, "fricState", sFricCut);
            matKSlip = matSpringDampMass(sPackCut, 0, sFric.K, 1, cellCommon{:}, "fricState", sFricCut, "slipMobilized", true);
            assert(sInfo.numNotInFricState == 0, 'runUnresolvedDOS:ContactMismatch', ...
                'matSpringDampMass found contacts that are not in the reduced fricState.');
            s.sStickCut = calcEigenFrequencies(matKStick, matMassCut);
            s.sSlipCut = calcEigenFrequencies(matKSlip, matMassCut);
            vecFnCut = abs(sFricCut.vecFn(:));
            s.scalZCut = 2*numel(vecFnCut)/numel(sPackCut.vecPosX);
            s.scalNmCut = nnz(abs(sFricCut.vecFt(:))./(sFricCut.mu*vecFnCut) >= 1 - options.scalMobTol)/numel(sPackCut.vecPosX);
            s.scalMeanKnCut = mean(contactStiffness(sPackCut, sFricCut, strContactLaw));
        end

        if isempty(sResults)
            sResults = s;
        else
            sResults(end + 1) = s; %#ok<AGROW>
        end
        save(strOutFile, 'sResults', 'sOptions');
        fprintf(['[runUnresolvedDOS] %s %s mu = %g, P = %g, seed %d: removed %d (rule %d + cascade %d) of %d', ...
            ' | z %.3f -> %.3f%s\n'], strContactLaw, options.strRule, s.scalMu, s.scalPressureTarget, s.scalSeed, ...
            s.scalNumRemoved, s.scalNumHit, s.scalNumCascade, scalNi, s.scalZ, s.scalZCut, repmat(' COLLAPSED', 1, s.boolCollapsed));
    end
    fprintf('[runUnresolvedDOS] %s done: %d packings -> %s\n', strContactLaw, numel(sResults), strOutFile);
end

function vecBoolKeep = cascade(vecBoolKeep, vecI, vecJ, scalNi, scalZMin)
% Remove kept particles with fewer than scalZMin contacts among kept
% particles, repeatedly, until nothing changes.
    while true
        vecBoolLive = vecBoolKeep(vecI) & vecBoolKeep(vecJ);
        vecZ = zeros(scalNi, 1);
        if any(vecBoolLive)
            vecZ = accumarray([vecI(vecBoolLive); vecJ(vecBoolLive)], 1, [scalNi 1]);
        end
        vecBoolDrop = vecBoolKeep & vecZ < scalZMin;
        if ~any(vecBoolDrop)
            return;
        end
        vecBoolKeep(vecBoolDrop) = false;
    end
end

function vecKn = contactStiffness(sPack, sFric, strContactLaw)
% Normal contact stiffness: K (Hooke) or 2 K sqrt(R_eff delta) (Hertz), as
% matSpringDampMass and runFrictionalDOS' <k_n>.
    if strContactLaw == "hertz"
        vecR = sPack.vecDiameter/2;
        vecRi = vecR(sFric.vecContactNN(:));
        vecRj = vecR(sFric.vecContactMM(:));
        vecKn = 2*sFric.K*sqrt(vecRi.*vecRj./(vecRi + vecRj).*sFric.vecOverlap(:));
    else
        vecKn = sFric.K*ones(numel(sFric.vecFn), 1);
    end
end

function [vecTanBeta, vecIdxTwo] = twoContactCone(sPack, sFric, vecZ)
% tan(beta) of every two-contact particle. Its two forces must act along the
% chord joining its contact points; that chord makes the angle
% beta = (pi - theta)/2 with each inward normal (theta: angle between the two
% outward contact normals), so the forces fit inside the Coulomb cones only
% if tan(beta) <= mu (Gilabert, Roux & Castellanos 2007).
    vecIdxTwo = find(vecZ == 2);
    vecTanBeta = NaN(numel(vecIdxTwo), 1);
    vecI = sFric.vecContactNN(:);
    vecJ = sFric.vecContactMM(:);
    for idxK = 1:numel(vecIdxTwo)
        idxP = vecIdxTwo(idxK);
        vecC = find(vecI == idxP | vecJ == idxP);
        matN = zeros(2, 2);
        for idxC = 1:2
            idxOther = vecJ(vecC(idxC));
            if idxOther == idxP
                idxOther = vecI(vecC(idxC));
            end
            vecSep = [sPack.vecPosX(idxOther) - sPack.vecPosX(idxP), sPack.vecPosY(idxOther) - sPack.vecPosY(idxP)];
            vecSep = vecSep - [sPack.scalBoxWidthX, sPack.scalBoxHeightY].*round(vecSep./[sPack.scalBoxWidthX, sPack.scalBoxHeightY]);
            matN(idxC, :) = vecSep/norm(vecSep);
        end
        scalTheta = acos(max(-1, min(1, dot(matN(1, :), matN(2, :)))));
        vecTanBeta(idxK) = tan((pi - scalTheta)/2);
    end
end

function [vecOmegaHat, vecShare, vecPR] = modeLocalization(matK, matMass, vecBoolRemoved, scalMeanKn)
% Modes of M^(-1/2) K M^(-1/2): frequency, participation ratio over
% particles, and the share of each (unit) mode on the removed particles.
    vecInvSqrtM = 1./sqrt(full(diag(matMass)));
    matD = full(matK).*(vecInvSqrtM*vecInvSqrtM.');
    matD = (matD + matD.')/2;
    [matV, matL] = eig(matD);
    vecOmegaSq = diag(matL);
    scalNi = numel(vecBoolRemoved);
    scalDof = size(matD, 1)/scalNi;                                   % 3 (x, y, theta)
    matWeight = reshape(sum(reshape(matV.^2, scalDof, scalNi, []), 1), scalNi, []);   % [Ni x modes]
    vecPR = (sum(matWeight, 1).^2./(scalNi*sum(matWeight.^2, 1))).';
    vecShare = sum(matWeight(vecBoolRemoved, :), 1).';
    vecOmegaHat = sign(vecOmegaSq).*sqrt(abs(vecOmegaSq))/sqrt(scalMeanKn);
end

function [sPackCut, sFricCut] = subsetBackbone(sPack, sFric, vecBoolKeep)
% The backbone without the removed particles; fricState re-indexed onto it
% (as pack.m's reindexFricState).
    vecIdxNew = zeros(numel(vecBoolKeep), 1);
    vecIdxNew(vecBoolKeep) = 1:nnz(vecBoolKeep);
    sPackCut = sPack;
    sPackCut.vecPosX = sPack.vecPosX(vecBoolKeep);
    sPackCut.vecPosY = sPack.vecPosY(vecBoolKeep);
    sPackCut.vecDiameter = sPack.vecDiameter(vecBoolKeep);
    vecBoolLive = vecBoolKeep(sFric.vecContactNN) & vecBoolKeep(sFric.vecContactMM);
    sFricCut = sFric;
    for strField = ["vecOverlap", "vecFn", "vecFt", "vecUnitTanX", "vecUnitTanY"]
        sFricCut.(strField) = sFric.(strField)(vecBoolLive);
    end
    sFricCut.vecContactNN = vecIdxNew(sFric.vecContactNN(vecBoolLive));
    sFricCut.vecContactMM = vecIdxNew(sFric.vecContactMM(vecBoolLive));
    sFricCut.vecRadius = sFric.vecRadius(vecBoolKeep);
    sFricCut.vecInertia = sFric.vecInertia(vecBoolKeep);
    sFricCut.N = nnz(vecBoolKeep);
end
