% testMatSpringDampMassFrictional.m
%
% Tests the frictional-eigenmode options of matSpringDampMass.m ("stressed",
% "fricState", "slipMobilized") on small frictional packings made by pack.m,
% and checks those packings against known results.
%
% Run from the GranE repo root:
%     matlab -batch "run('tests/testMatSpringDampMass/testMatSpringDampMassFrictional.m')"
%
% The packings are written to tests/testMatSpringDampMass/data/mu<mu>/ and are
% kept after the run. Every run deletes and rebuilds them first, because pack.m
% skips a packing whose .mat already exists (it would replay old numbers).
%
% References
%   [S07] Shundyak, van Hecke & van Saarloos, PRE 75, 010301(R) (2007)
%   [H10] Henkes, van Hecke & van Saarloos, EPL 90, 14003 (2010)
%   [W05] Wyart, Silbert, Nagel & Witten, PRE 72, 051306 (2005)

clearvars

strTestDir  = fileparts(mfilename('fullpath'));
strRepoRoot = fileparts(fileparts(strTestDir));
addpath(fullfile(strRepoRoot, 'src'));
strDataDir  = fullfile(strTestDir, 'data');

% ---- packing parameters (pack.m units: D_small = 1, mass = 1) ----
scalN    = 100;            % particles
scalK    = 100;            % normal spring constant k_n
scalP    = 1e-3;           % target pressure
scalSeed = 1;
vecMu    = [0.1 0.5 1.0];  % friction coefficients, low to high
scalTolMob = 1e-6;         % m >= 1 - tol counts as fully mobilized

vecPhi = zeros(size(vecMu));   % packing fraction (rattlers included)
vecZ   = zeros(size(vecMu));   % backbone contact number z = 2 Nc / Ni
vecNm  = zeros(size(vecMu));   % fully mobilized contacts per particle

for idxMu = 1:numel(vecMu)
    scalMu = vecMu(idxMu);
    fprintf('\n===== mu = %g =====\n', scalMu);

    % ------------------------------------------------------------------
    % 1) Build a fresh frictional packing with pack.m
    % ------------------------------------------------------------------
    strSub = fullfile(strDataDir, sprintf('mu%g', scalMu));
    if ~isfolder(strSub)
        mkdir(strSub);
    end
    delete(fullfile(strSub, '*.mat'));
    sPackOpts = struct('flagFrictionOn', true, 'scalFricCoef', scalMu, 'scalTangentialK', 1/3, ...
        'scalGammaNormal', 0, 'scalGammaTangential', 0, ...
        'saveFrictionalState', true, 'saveFullState', true);
    pack(scalN, scalK, 1, 1.4, 1, scalP, scalSeed, false, 1, 1, 0, false, [strSub filesep], sPackOpts);
    close all;

    strStem = fullfile(strSub, sprintf('2D_N%d_P%s_Width%d_Seed%d', scalN, num2str(scalP), round(sqrt(scalN)), scalSeed));
    sPackFile = load([strStem '_Friction.mat']);     % pack.m's one output file
    sFull = sPackFile.sFullState;                    % all particles, before cleanRats
    fs = sFull.fricState;                            % per-contact forces, same indexing
    vecR = sFull.vecDiameter/2;
    vecBox = [sFull.scalBoxWidthX, sFull.scalBoxHeightY];
    matPos = [sFull.vecPosX, sFull.vecPosY];

    % ------------------------------------------------------------------
    % 2) Audit pack.m: is the saved state a frictional mechanical equilibrium?
    % ------------------------------------------------------------------
    vecI = fs.vecContactNN(:);
    vecJ = fs.vecContactMM(:);

    % 2a) fricState contacts are exactly the geometric overlaps (minimum image)
    [vecIg, vecJg, matNg, vecRg] = overlapList(matPos, vecR, vecBox);
    assert(isequal(sortrows([vecI vecJ]), sortrows([vecIg vecJg])), 'fricState contacts differ from geometric overlaps');
    matN = [-fs.vecUnitTanY(:), fs.vecUnitTanX(:)];   % fricState stores t = (n_y, -n_x)
    [~, idxG] = ismember([vecI vecJ], [vecIg vecJg], 'rows');
    assert(max(vecnorm(matN - matNg(idxG, :), 2, 2)) < 1e-12, 'fricState normals differ from geometry');
    assert(max(abs(fs.vecOverlap(:) - (vecR(vecI) + vecR(vecJ) - vecRg(idxG)))) < 1e-12, 'fricState overlaps differ from geometry');

    % 2b) Coulomb: |f_t| <= mu f_n on every contact
    vecMob = abs(fs.vecFt(:))./(scalMu*abs(fs.vecFn(:)));
    assert(max(vecMob) <= 1 + 1e-12, 'Coulomb bound violated');

    % 2c) force and 2d) torque balance on every particle. pack.m only checks
    % force balance while it runs, so torque balance is an independent check.
    matT = [fs.vecUnitTanX(:), fs.vecUnitTanY(:)];
    matFc = fs.vecFn(:).*matN + fs.vecFt(:).*matT;    % force on particle i of contact (i, j)
    matFnet = [accumarray(vecI, matFc(:,1), [scalN 1]) - accumarray(vecJ, matFc(:,1), [scalN 1]), ...
               accumarray(vecI, matFc(:,2), [scalN 1]) - accumarray(vecJ, matFc(:,2), [scalN 1])];
    vecTnet = -(accumarray(vecI, vecR(vecI).*fs.vecFt(:), [scalN 1]) + accumarray(vecJ, vecR(vecJ).*fs.vecFt(:), [scalN 1]));
    scalFbar = mean(abs(fs.vecFn));
    scalForceRatio  = max(vecnorm(matFnet, 2, 2))/scalFbar;
    scalTorqueRatio = max(abs(vecTnet))/(scalFbar*mean(vecR));
    fprintf('pack.m audit: max|F_net|/<f_n> = %.2e, max|tau_net|/(<f_n><R>) = %.2e, max m = %.15f\n', ...
        scalForceRatio, scalTorqueRatio, max(vecMob));
    assert(scalForceRatio < 1e-2, 'packing is not in force balance');
    assert(scalTorqueRatio < 1e-2, 'packing is not in torque balance');

    % 2e) pressure inside pack.m's +-15% dead band
    assert(abs(sFull.scalPressure/sFull.P_target - 1) < 0.15, 'pressure outside the dead band');

    % ------------------------------------------------------------------
    % 3) Backbone: drop particles with < 2 contacts until none are left
    %    (the cleanRats rule for friction), and reindex fricState to match.
    % ------------------------------------------------------------------
    boolKeep = true(scalN, 1);
    while true
        boolLive = boolKeep(vecI) & boolKeep(vecJ);
        vecZi = accumarray([vecI(boolLive); vecJ(boolLive)], 1, [scalN 1]);
        boolDrop = boolKeep & vecZi < 2;
        if ~any(boolDrop)
            break;
        end
        boolKeep(boolDrop) = false;
    end
    boolLive = boolKeep(vecI) & boolKeep(vecJ);
    vecNew = zeros(scalN, 1);
    vecNew(boolKeep) = 1:sum(boolKeep);
    scalNi = sum(boolKeep);
    sPk = struct('vecPosX', sFull.vecPosX(boolKeep), 'vecPosY', sFull.vecPosY(boolKeep), ...
        'scalBoxWidthX', sFull.scalBoxWidthX, 'scalBoxHeightY', sFull.scalBoxHeightY, ...
        'vecDiameter', sFull.vecDiameter(boolKeep));
    fsB = fs;
    fsB.N = scalNi;
    fsB.vecRadius = vecR(boolKeep);
    fsB.vecContactNN = vecNew(vecI(boolLive));
    fsB.vecContactMM = vecNew(vecJ(boolLive));
    fsB.vecFn = fs.vecFn(boolLive);
    fsB.vecFt = fs.vecFt(boolLive);
    vecRB = fsB.vecRadius;
    vecIB = min(fsB.vecContactNN, fsB.vecContactMM);
    vecJB = max(fsB.vecContactNN, fsB.vecContactMM);
    boolSlip = abs(fsB.vecFt)./(scalMu*abs(fsB.vecFn)) >= 1 - scalTolMob;
    scalNc = numel(vecIB);
    scalNm = sum(boolSlip);

    % ------------------------------------------------------------------
    % 4) Known results for 2D frictional jamming [S07]
    % ------------------------------------------------------------------
    vecPhi(idxMu) = sum(pi*vecR.^2)/prod(vecBox);
    vecZ(idxMu)   = 2*scalNc/scalNi;
    vecNm(idxMu)  = scalNm/scalNi;
    fprintf('phi = %.4f, backbone Ni = %d, Nc = %d, Nm = %d, z = %.3f, n_m = %.3f, z - (3 + n_m) = %.3f\n', ...
        vecPhi(idxMu), scalNi, scalNc, scalNm, vecZ(idxMu), vecNm(idxMu), vecZ(idxMu) - 3 - vecNm(idxMu));
    % [S07] Figs. 3-4: at p -> 0, z runs from 3.98 (mu -> 0) to 3.00 (mu -> inf)
    % and phi (rattlers included) from 0.8395 to 0.758; finite p raises both a little.
    assert(vecZ(idxMu) >= 3 && vecZ(idxMu) <= 4.25, 'backbone z outside the 2D frictional range [3, 4]');
    assert(vecPhi(idxMu) >= 0.74 && vecPhi(idxMu) <= 0.86, 'packing fraction outside the 2D frictional range');
    % [S07] eq. (1): a stable packing needs z >= z_iso^m = 3 + n_m (2D)
    assert(vecZ(idxMu) >= 3 + vecNm(idxMu), 'packing is left of the generalized isostaticity line');

    % ------------------------------------------------------------------
    % 5) matSpringDampMass on the backbone
    % ------------------------------------------------------------------
    cellCommon = {"periodic", true, "uniformMass", true, "friction", true, "algorithm", "particle", "fricState", fsB};
    [matKplain, ~, matM]    = matSpringDampMass(sPk, 0, scalK, 1, "periodic", true, "uniformMass", true, "friction", true, "algorithm", "particle");
    [matKns, ~, ~, sInfo]   = matSpringDampMass(sPk, 0, scalK, 1, cellCommon{:});                                       % no slip, unstressed
    matKnsS = matSpringDampMass(sPk, 0, scalK, 1, cellCommon{:}, "stressed", true);                                      % no slip, stressed
    matKsl  = matSpringDampMass(sPk, 0, scalK, 1, cellCommon{:}, "slipMobilized", true);                                 % slip, unstressed
    matKslS = matSpringDampMass(sPk, 0, scalK, 1, cellCommon{:}, "slipMobilized", true, "stressed", true);               % slip, stressed [H10]

    % 5a) a fricState alone changes nothing; contactInfo counts are right
    assert(isequal(matKns, matKplain), 'passing fricState without slipMobilized changed K');
    assert(sInfo.numContacts == scalNc && sInfo.numMobilized == scalNm && sInfo.numNotInFricState == 0, 'contactInfo counts are wrong');
    assert(isequal(sInfo.vecNumContacts, accumarray([vecIB; vecJB], 1, [scalNi 1])), 'contactInfo.vecNumContacts is wrong');

    % 5b) cell and particle algorithms agree with every option on
    matKcell = matSpringDampMass(sPk, 0, scalK, 1, cellCommon{1:6}, "algorithm", "cell", "fricState", fsB, ...
        "slipMobilized", true, "stressed", true);
    assert(norm(matKcell - matKslS, 'fro') < 1e-12*norm(matKslS, 'fro'), 'cell and particle K disagree');

    % 5c) K reproduces eq. (2) of [H10] contact by contact:
    %     r'Kr = sum_c [ k_n (dr.n)^2 - (f_n/r)(dr.t)^2 + k_t,c (dt)^2 ],
    %     with k_t,c = 0 on slipping contacts; Gamma keeps its tangential dashpot.
    scalKt = scalK/3;
    scalGn = 0.2;
    scalGt = 0.2;
    vecRand = randn(3*scalNi, 1);
    vecVel  = randn(3*scalNi, 1);
    [matKd, matGd, matMd] = matSpringDampMass(sPk, scalGn, scalK, 1, cellCommon{:}, "slipMobilized", true, ...
        "stressed", true, "tangentialDampingConstant", scalGt);
    [scalEK, scalEG] = henkesQuadraticForm(vecRand, vecVel, [sPk.vecPosX sPk.vecPosY], vecRB, vecBox, ...
        vecIB, vecJB, boolSlip, scalK, scalKt, scalGn, scalGt);
    assert(norm(matKd - matKd.', 'fro') < 1e-12*norm(matKd, 'fro'), 'K is not symmetric');
    assert(abs(vecRand.'*matKd*vecRand - scalEK) < 1e-9*abs(scalEK), 'K does not match Henkes eq. (2)');
    assert(abs(vecVel.'*matGd*vecVel - scalEG) < 1e-9*abs(scalEG), 'Gamma changed on slipping contacts');

    % 5d) the pre-stress term makes K the exact Hessian of the Hookean energy
    %     of the contact network: compare u'Ku with a central second difference
    %     of E along random directions u. E sums over the contacts present now
    %     (the network K linearizes), so it stays smooth across the stencil.
    matKfl  = matSpringDampMass(sPk, 0, scalK, 1, "periodic", true, "algorithm", "particle", "stressed", true);
    matKfl0 = matSpringDampMass(sPk, 0, scalK, 1, "periodic", true, "algorithm", "particle");
    vecX0 = reshape([sPk.vecPosX sPk.vecPosY].', [], 1);
    fnEnergy = @(x) networkEnergy(reshape(x, 2, []).', vecRB, vecBox, vecIB, vecJB, scalK);
    scalH = 1e-4;
    for idxDir = 1:3
        vecU = randn(2*scalNi, 1);
        vecU = vecU/norm(vecU);
        scalD2 = (fnEnergy(vecX0 + scalH*vecU) - 2*fnEnergy(vecX0) + fnEnergy(vecX0 - scalH*vecU))/scalH^2;
        scalErrS = abs(vecU.'*matKfl*vecU - scalD2)/abs(scalD2);
        scalErrU = abs(vecU.'*matKfl0*vecU - scalD2)/abs(scalD2);
        fprintf('  d2E/du2 check: stressed rel err = %.1e (unstressed %.1e)\n', scalErrS, scalErrU);
        assert(scalErrS < 1e-6, 'stressed K is not the Hessian of the Hookean energy');
    end

    % ------------------------------------------------------------------
    % 6) Undamped spectra (Gamma = 0): omega^2 = eig(M^-1/2 K M^-1/2)
    % ------------------------------------------------------------------
    matMh = diag(1./sqrt(diag(matM)));
    fnW2 = @(matKx) sort(eig(matMh*((matKx + matKx.')/2)*matMh));
    vecW2ns  = fnW2(matKns);
    vecW2nsS = fnW2(matKnsS);
    vecW2sl  = fnW2(matKsl);
    scalW2max = max(vecW2ns);
    scalTolZero = 1e-9*scalW2max;
    fnZeros = @(w2) sum(abs(w2) < scalTolZero);

    % 6a) slipping only removes stiffness (K_noslip - K_slip = sum k_t t t' >= 0),
    %     so every sorted omega^2 can only go down (Courant-Fischer) [H10 Fig. 3a vs 3b]
    assert(max(vecW2sl - vecW2ns) <= scalTolZero, 'slipping raised a frequency');

    % 6b) no slip, unstressed: positive semi-definite, and the zero modes are
    %     exactly what constraint counting gives (2 rigid translations when hyperstatic)
    assert(min(vecW2ns) > -scalTolZero, 'unstressed K has a negative eigenvalue');
    assert(fnZeros(vecW2ns) == max(2, 3*scalNi - 2*scalNc), 'no-slip zero-mode count differs from constraint counting');

    % 6c) slip, unstressed: counting gives a lower bound on the floppy modes,
    %     globally (3 Ni DOF vs 2 Nc - Nm constraints) and locally (a particle
    %     with z_i contacts, m_i of them slipping, has 2 z_i - m_i constraints)
    vecZiB = accumarray([vecIB; vecJB], 1, [scalNi 1]);
    vecMiB = accumarray([vecIB(boolSlip); vecJB(boolSlip)], 1, [scalNi 1]);
    scalLocal = sum(max(0, 3 - (2*vecZiB - vecMiB)));
    % a particle whose contacts all slip has no rotational stiffness: free spin
    scalFreeSpin = sum(vecMiB == vecZiB);
    fprintf('  zero modes: no slip %d, slip %d (global count >= %d, local count >= %d, free spins %d)\n', ...
        fnZeros(vecW2ns), fnZeros(vecW2sl), max(2, 3*scalNi - (2*scalNc - scalNm)), 2 + scalLocal, scalFreeSpin);
    assert(fnZeros(vecW2sl) >= max(2, 3*scalNi - (2*scalNc - scalNm)), 'slip zero modes below the global count');
    assert(fnZeros(vecW2sl) >= 2 + scalLocal, 'slip zero modes below the local count');
    assert(fnZeros(vecW2sl) >= 2 + scalFreeSpin, 'slip zero modes below the free-spin count');

    % 6d) all contacts slipping (k_t -> 0, [H10] Fig. 3d): rotations decouple, so
    %     the spectrum is the frictionless one plus Ni zero-frequency rotations
    matKall = matSpringDampMass(sPk, 0, scalK, 1, cellCommon{:}, "slipMobilized", true, "mobilizationTol", 1);
    matKfr  = matSpringDampMass(sPk, 0, scalK, 1, "periodic", true, "uniformMass", true, "algorithm", "particle");
    vecW2all = fnW2(matKall);
    vecW2fr  = sort(eig((matKfr + matKfr.')/2));    % uniform mass 1
    vecNzAll = vecW2all(abs(vecW2all) >= scalTolZero);
    vecNzFr  = vecW2fr(abs(vecW2fr) >= scalTolZero);
    assert(fnZeros(vecW2all) == scalNi + fnZeros(vecW2fr), 'all-slip zero modes != Ni + frictionless zero modes');
    assert(numel(vecNzAll) == numel(vecNzFr) && norm(vecNzAll - vecNzFr) < 1e-10*norm(vecNzFr), ...
        'all-slip spectrum differs from the frictionless spectrum');

    % 6e) no slip, stressed: the packing relaxed with sticking contacts, so it
    %     is mechanically stable and the stressed K stays positive semi-definite
    assert(min(vecW2nsS) > -1e-6*scalW2max, 'stressed no-slip K is unstable');

    % ------------------------------------------------------------------
    % 7) Damped QEP: the number of growing modes (Re lambda > 0) equals the
    %    number of negative omega^2 of K (Kelvin-Tait-Chetaev theorem)
    % ------------------------------------------------------------------
    for boolSlipOn = [false true]
        [matKq, matGq, matMq] = matSpringDampMass(sPk, scalGn, scalK, 1, cellCommon{:}, ...
            "slipMobilized", boolSlipOn, "stressed", true, "tangentialDampingConstant", scalGt);
        vecLambda = polyeig(matKq, matGq, matMq);
        scalGrowing = sum(real(vecLambda) > 1e-6*max(abs(vecLambda)));
        scalNegative = sum(fnW2(matKq) < -scalTolZero);
        fprintf('  damped QEP (slip = %d, stressed): growing modes = %d, negative omega^2 = %d\n', ...
            boolSlipOn, scalGrowing, scalNegative);
        assert(numel(vecLambda) == 6*scalNi, 'QEP should have 6 Ni eigenvalues');
        assert(scalGrowing == scalNegative, 'growing damped modes != negative omega^2');
    end
end

% ----------------------------------------------------------------------
% 8) Trends with friction [S07 Figs. 1, 3, 4]: phi, z and n_m all fall as mu rises
% ----------------------------------------------------------------------
fprintf('\nmu  = %s\nphi = %s\nz   = %s\nn_m = %s\n', mat2str(vecMu), mat2str(vecPhi, 4), mat2str(vecZ, 4), mat2str(vecNm, 3));
assert(vecPhi(1) > vecPhi(end) && vecZ(1) > vecZ(end) && vecNm(1) > vecNm(end), ...
    'phi, z and n_m should fall from the lowest to the highest friction');

% ----------------------------------------------------------------------
% 9) misuse is caught
% ----------------------------------------------------------------------
assertErrorId(@() matSpringDampMass(sPk, 0, scalK, 1, "periodic", true, "friction", true, "slipMobilized", true), ...
    'matSpringDampMass:SlipNeedsFricState');
assertErrorId(@() matSpringDampMass(sPk, 0, scalK, 1, "periodic", true, "slipMobilized", true, "fricState", fsB), ...
    'matSpringDampMass:SlipNeedsFricState');
fsBad = fsB;
fsBad.N = scalNi - 1;
assertErrorId(@() matSpringDampMass(sPk, 0, scalK, 1, "periodic", true, "friction", true, "fricState", fsBad), ...
    'matSpringDampMass:FricStateMismatch');

fprintf('testMatSpringDampMassFrictional passed.\n');

% ======================================================================
function [vecI, vecJ, matN, vecDist] = overlapList(matPos, vecR, vecBox)
% All overlapping pairs i < j (periodic minimum image), unit normal from i to j.
    scalN = size(matPos, 1);
    vecI = []; vecJ = []; matN = zeros(0, 2); vecDist = [];
    for i = 1:scalN-1
        matD = matPos(i+1:scalN, :) - matPos(i, :);
        matD = matD - round(matD./vecBox).*vecBox;
        vecDst = sqrt(sum(matD.^2, 2));
        vecHit = find(vecDst < vecR(i) + vecR(i+1:scalN));
        vecHit = vecHit(:);    % find() on a 1x1 comparison returns a 1x0 row
        vecI = [vecI; i*ones(numel(vecHit), 1)]; %#ok<AGROW>
        vecJ = [vecJ; i + vecHit]; %#ok<AGROW>
        matN = [matN; matD(vecHit, :)./vecDst(vecHit)]; %#ok<AGROW>
        vecDist = [vecDist; vecDst(vecHit)]; %#ok<AGROW>
    end
end

function scalE = networkEnergy(matPos, vecR, vecBox, vecI, vecJ, scalK)
% E = sum over the given contacts (i, j) of 1/2 k (R_i + R_j - r)^2 (periodic).
    matD = matPos(vecJ, :) - matPos(vecI, :);
    matD = matD - round(matD./vecBox).*vecBox;
    vecOv = vecR(vecI) + vecR(vecJ) - sqrt(sum(matD.^2, 2));
    scalE = 0.5*scalK*sum(vecOv.^2);
end

function [scalEK, scalEG] = henkesQuadraticForm(vecR3, vecV3, matPos, vecR, vecBox, vecI, vecJ, boolSlip, kn, kt, gn, gt)
% 2V of [H10] eq. (2) with k_t = 0 on slipping contacts, and the matching
% dashpot form (normal + tangential, kept on every contact).
    scalEK = 0;
    scalEG = 0;
    for c = 1:numel(vecI)
        i = vecI(c);
        j = vecJ(c);
        vecD = (matPos(j, :) - matPos(i, :)).';
        vecD = vecD - round(vecD./vecBox.').*vecBox.';
        scalR = norm(vecD);
        vecN = vecD/scalR;
        vecT = [-vecN(2); vecN(1)];
        scalFn = kn*(vecR(i) + vecR(j) - scalR);
        vecDu = vecR3(3*j-2:3*j-1) - vecR3(3*i-2:3*i-1);
        scalDt = vecT.'*vecDu - (vecR(i)*vecR3(3*i) + vecR(j)*vecR3(3*j));
        scalEK = scalEK + kn*(vecN.'*vecDu)^2 - (scalFn/scalR)*(vecT.'*vecDu)^2 + kt*(~boolSlip(c))*scalDt^2;
        vecDv = vecV3(3*j-2:3*j-1) - vecV3(3*i-2:3*i-1);
        scalDtv = vecT.'*vecDv - (vecR(i)*vecV3(3*i) + vecR(j)*vecV3(3*j));
        scalEG = scalEG + gn*(vecN.'*vecDv)^2 + gt*scalDtv^2;
    end
end

function assertErrorId(fnCall, strId)
    try
        fnCall();
    catch err
        assert(strcmp(err.identifier, strId), 'expected error %s, got %s', strId, err.identifier);
        return;
    end
    error('expected error %s, but none was thrown', strId);
end
