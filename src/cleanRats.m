function [matPositions, vecRadii, vecBoolKeep, sUncert] = cleanRats(matPositions, vecRadii, scalBoxHeightY, scalBoxWidthX, scalBoxDepthZ, boolFullyPeriodic, scalMu, sUncertOptions)
% cleanRats(pos, radii, Ly, Lx)                   — 2D, PBC in y only, frictionless
% cleanRats(pos, radii, Ly, Lx, Lz)              — 3D, PBC in y and z, frictionless
% cleanRats(pos, radii, Ly, Lx, Lz, boolFullyP)  — + full PBC flag (Lz = [] for 2D)
% cleanRats(..., mu)                              — mu > 0: frictional variant.
% cleanRats(..., mu, sUncertOptions)                — + statistics of the weak contacts
% [pos, radii, keep, sUncert] = cleanRats(...)      — (force below a noise floor), see below.
%
% Removes rattlers, repeatedly, until every remaining particle has enough
% contacts. Returns the surviving positions and radii (input order kept) and
% vecBoolKeep [N_original x 1], true for the particles that survive, so
% callers can re-index per-particle or per-contact data (pack.m uses it to
% re-index the frictional contact forces onto the backbone).
%
% A contact is a real overlap, r_i + r_j - |d_ij| > 0, with the minimum-image
% separation along every periodic axis. This is the same definition as
% computeMeanCoordNum and pack.m's buildFricState, so for a frictional
% packing the backbone is exactly the one obtained by dropping particles with
% fewer than 2 entries in fricState's contact list. (Before, near-touching
% pairs with a gap up to 1e-3 also counted as contacts, so particles that
% touch nothing were kept.)
%
% Periodicity: y (and z in 3D) are always periodic. x is periodic only when
% boolFullyPeriodic is true; otherwise x has walls and particles touching a
% wall get 2 contacts credited. pack.m's packings are fully periodic, so it
% calls with boolFullyPeriodic = true.
%
% Frictional variant (2D): when mu > 0 the removal threshold is Zn > 1,
% not the frictionless Zn > 2.
%
% PHYSICS — isostatic coordination number (Maxwell counting, 2D disks):
%   Frictionless: DOF = 2N, constraints = Nz/2  =>  z_iso = 4.
%   Frictional:   DOF = 3N, constraints = Nz    =>  z_iso = 3.
% In this gravity-free, isotropically compressed packing, a particle with
% Zn <= 1 cannot satisfy force–torque balance: with one contact,
% sum(F) = F_1 and sum(tau) = r x F_1, both nonzero unless F_1 = 0. Such
% particles are rattlers even with friction, so for mu > 0 we keep Zn >= 2
% grains and remove Zn <= 1 (floaters + one-contact rattlers) as
% frictional rattlers.
%
% mu is the 7th positional argument (callers that omit it get mu = 0, the
% frictionless rule); scalBoxDepthZ = [] signals a 2D call.
%
% UNCERTAIN (SUB-NOISE) CONTACTS — optional, off unless asked for. A contact whose
% normal force is below scalForceFloor * <f> (default 1e-4, pack.m's force-
% balance tolerance) is below the equilibration noise: whether it exists at
% all is not resolved, yet it is counted like any other contact. With an
% 8th argument sUncertOptions (struct; all fields optional), or when the 4th
% output is requested, cleanRats also returns sUncert, the statistics of
% those contacts on THIS packing, and with sUncertOptions.boolPlot = true it
% draws them with plotUncertContacts. The forces come from the overlaps alone:
% f ~ delta (Hooke) or sqrt(R_eff) delta^(3/2) (Hertz), so f/<f> needs no
% spring constant. <f> is the mean over ALL contacts of the input packing.
%   sUncertOptions.boolPlot       draw plotUncertContacts for this packing (false)
%   sUncertOptions.strContactLaw  'hooke' or 'hertz' ('hooke')
%   sUncertOptions.scalPressure   pressure label, used for the colour (NaN)
%   sUncertOptions.scalForceFloor noise floor in units of <f> (1e-4)
%   sUncertOptions.boolLightMode  light (true) or dark (false) figure
%   sUncertOptions.strSavePath    PNG to export the figure to ("": none)
% The removal itself (positions, radii, vecBoolKeep) is unchanged by these
% options.
    if nargin < 6
        boolFullyPeriodic = false;
    end
    if nargin < 7
        scalMu = 0;
    end
    boolThreeD = nargin >= 5 && ~isempty(scalBoxDepthZ);
    if nargin < 8
        sUncertOptions = struct();
    end
    sUncertOptions = fillUncertOptions(sUncertOptions);
    boolUncert = nargout >= 4 || sUncertOptions.boolPlot;
    if boolUncert
        % contact list of the INPUT packing (before any removal)
        matPositionsIn = matPositions;
        vecRadiiIn     = vecRadii;
    end

    boolFriction = (scalMu > 0);
    if boolFriction
        fprintf('[cleanRats] mu=%.4f: frictional variant — keeping Zn>=2 backbone, removing Zn<=1 (floaters + one-contact rattlers)\n', scalMu);
    end

    boolChanged        = true;
    scalPassNumber     = 0;
    scalTotalRattlers  = 0;
    scalNumOriginal    = size(matPositions, 1);
    vecIdxSurvivor     = (1:scalNumOriginal)';   % [N x 1] original index of each remaining particle

    while boolChanged
        scalPassNumber = scalPassNumber + 1;
        fprintf('[cleanRats] Pass %d\n', scalPassNumber);
        scalNumParticles = size(matPositions, 1);

        % Recompute wall lists on current positions (only used when NOT fully periodic)
        fprintf('[cleanRats] Recomputing wall lists...\n');
        vecBoolLeftWall  = matPositions(:,1) < vecRadii;
        vecBoolRightWall = matPositions(:,1) > scalBoxWidthX - vecRadii;

        vecZn = zeros(scalNumParticles, 1);
        if ~boolFullyPeriodic
            vecZn(vecBoolLeftWall | vecBoolRightWall) = 2;
        end

        fprintf('[cleanRats] Computing contact lists...\n');
        for idxI = 1:scalNumParticles
            for idxJ = idxI+1:scalNumParticles
                scalSepX = matPositions(idxI,1) - matPositions(idxJ,1);
                if boolFullyPeriodic
                    scalSepX = scalSepX - round(scalSepX / scalBoxWidthX) * scalBoxWidthX;
                end
                scalSepY = matPositions(idxI,2) - matPositions(idxJ,2);
                scalSepY = scalSepY - round(scalSepY / scalBoxHeightY) * scalBoxHeightY;

                if boolThreeD
                    scalSepZ = matPositions(idxI,3) - matPositions(idxJ,3);
                    scalSepZ = scalSepZ - round(scalSepZ / scalBoxDepthZ) * scalBoxDepthZ;
                    scalDist = sqrt(scalSepX^2 + scalSepY^2 + scalSepZ^2);
                else
                    scalDist = sqrt(scalSepX^2 + scalSepY^2);
                end

                if vecRadii(idxI) + vecRadii(idxJ) - scalDist > 0   % real overlap only
                    vecZn(idxI) = vecZn(idxI) + 1;
                    vecZn(idxJ) = vecZn(idxJ) + 1;
                end
            end
        end

        fprintf('[cleanRats] Zn distribution: min=%d, max=%d, mean=%.2f\n', min(vecZn), max(vecZn), mean(vecZn));

        if boolFriction
             % Frictional 2D: remove Zn<=1 (floaters + one-contact rattlers).
             % In this gravity-free model a 1-contact disk cannot satisfy
             % force–torque balance, so it is a rattler even with friction.
             % Keep only Zn >= 2 grains (frictional 2D isostatic z_iso = 3).
            vecBoolToKeep = vecZn > 1;
        elseif boolThreeD
            vecBoolToKeep = vecZn > 3;
        else
            vecBoolToKeep = vecZn > 2;
        end

        boolChanged       = ~all(vecBoolToKeep);
        matPositions      = matPositions(vecBoolToKeep, :);
        vecRadii          = vecRadii(vecBoolToKeep);
        vecIdxSurvivor    = vecIdxSurvivor(vecBoolToKeep);
        scalTotalRattlers = scalTotalRattlers + sum(~vecBoolToKeep);
        fprintf('[cleanRats] Pass %d: %d rattlers removed.\n', scalPassNumber, sum(~vecBoolToKeep));
    end

    vecBoolKeep = false(scalNumOriginal, 1);
    vecBoolKeep(vecIdxSurvivor) = true;

    sUncert = struct([]);
    if boolUncert
        sUncert = uncertContactStats(matPositionsIn, vecRadiiIn, scalBoxHeightY, scalBoxWidthX, scalBoxDepthZ, ...
            boolThreeD, boolFullyPeriodic, scalMu, vecBoolKeep, sUncertOptions);
        if sUncertOptions.boolPlot
            plotUncertContacts(sUncert, 'boolLightMode', sUncertOptions.boolLightMode, 'strSavePath', sUncertOptions.strSavePath);
        end
    end

    fprintf('[cleanRats] %d total rattlers removed.\n', scalTotalRattlers);
    fprintf('[cleanRats] Total percentage of rattlers removed: %.2f%%\n', scalTotalRattlers / scalNumOriginal * 100);
end


function sOpt = fillUncertOptions(sOpt)
% Defaults for the optional weak-contact statistics.
    sDefault = struct('boolPlot', false, 'strContactLaw', 'hooke', 'scalPressure', NaN, ...
        'scalForceFloor', 1e-4, 'boolLightMode', true, 'strSavePath', "");
    cellNames = fieldnames(sDefault);
    for idxName = 1:numel(cellNames)
        if ~isfield(sOpt, cellNames{idxName})
            sOpt.(cellNames{idxName}) = sDefault.(cellNames{idxName});
        end
    end
    assert(any(strcmpi(sOpt.strContactLaw, {'hooke', 'hertz'})), 'cleanRats:ContactLaw', ...
        'sUncertOptions.strContactLaw must be ''hooke'' or ''hertz''.');
end

function sUncert = uncertContactStats(matPositions, vecRadii, scalBoxHeightY, scalBoxWidthX, scalBoxDepthZ, ...
        boolThreeD, boolFullyPeriodic, scalMu, vecBoolKeep, sOpt)
% Statistics of the contacts whose force is below sOpt.scalForceFloor * <f>,
% and what removing them would do to the backbone. Uses the same contact
% definition (real overlap, minimum image) and removal rule as cleanRats.
    scalN = size(matPositions, 1);
    vecL = [scalBoxWidthX, scalBoxHeightY];
    if boolThreeD
        vecL = [vecL, scalBoxDepthZ];
    end
    vecBoolWrap = true(size(vecL));
    vecBoolWrap(1) = boolFullyPeriodic;

    % contact list, i < j (vectorized over j for each i: O(N) memory)
    cellI = cell(scalN, 1); cellJ = cell(scalN, 1); cellDelta = cell(scalN, 1);
    for idxI = 1:scalN - 1
        vecJ = (idxI + 1:scalN)';
        matSep = matPositions(vecJ, :) - matPositions(idxI, :);
        for idxDim = find(vecBoolWrap)
            matSep(:, idxDim) = matSep(:, idxDim) - vecL(idxDim) * round(matSep(:, idxDim) / vecL(idxDim));
        end
        vecDelta = vecRadii(idxI) + vecRadii(vecJ) - sqrt(sum(matSep.^2, 2));
        vecBoolC = vecDelta > 0;
        cellI{idxI} = idxI * ones(nnz(vecBoolC), 1);
        cellJ{idxI} = vecJ(vecBoolC);
        cellDelta{idxI} = vecDelta(vecBoolC);
    end
    vecI = vertcat(cellI{:}); vecJ = vertcat(cellJ{:}); vecDelta = vertcat(cellDelta{:});

    % relative force and stiffness of every contact (no spring constant needed)
    vecREff = vecRadii(vecI) .* vecRadii(vecJ) ./ (vecRadii(vecI) + vecRadii(vecJ));
    if strcmpi(sOpt.strContactLaw, 'hertz')
        vecForce = sqrt(vecREff) .* vecDelta.^(3/2);   % f / ((4/3) K)
        vecStiff = sqrt(vecREff .* vecDelta);          % k_n / (2 K)
    else
        vecForce = vecDelta;                           % f / K
        vecStiff = ones(size(vecDelta));               % k_n / K
    end
    vecForceRel = vecForce / mean(vecForce);
    vecBoolWeak = vecForceRel < sOpt.scalForceFloor;
    vecBoolBackbone = vecBoolKeep(vecI) & vecBoolKeep(vecJ);

    % the rattler rule of cleanRats, on a contact list (walls credited as in cleanRats)
    vecZWall = zeros(scalN, 1);
    if ~boolFullyPeriodic
        vecZWall(matPositions(:, 1) < vecRadii | matPositions(:, 1) > scalBoxWidthX - vecRadii) = 2;
    end
    if scalMu > 0
        scalZMin = 2;
    elseif boolThreeD
        scalZMin = 4;
    else
        scalZMin = 3;
    end
    vecBoolKeepCut = removeRattlers(vecI(~vecBoolWeak), vecJ(~vecBoolWeak), vecZWall, scalZMin);

    % coordination of the backbone with and without the weak contacts
    scalNi = nnz(vecBoolKeep);
    scalNcB = nnz(vecBoolBackbone);
    scalNiCut = nnz(vecBoolKeepCut);
    scalNcCut = nnz(~vecBoolWeak & vecBoolKeepCut(vecI) & vecBoolKeepCut(vecJ));
    scalD = 2 + boolThreeD;
    if scalMu > 0
        scalZIso = scalD + 1;    % frictional: d + 1 (no mobilized contacts counted here)
    else
        scalZIso = 2 * scalD;    % frictionless: 2d
    end
    scalZ = 2 * scalNcB / max(scalNi, 1);
    scalZCut = 2 * scalNcCut / max(scalNiCut, 1);

    sUncert = struct();
    sUncert.strContactLaw = lower(char(sOpt.strContactLaw));
    sUncert.scalMu = scalMu;
    sUncert.scalPressure = sOpt.scalPressure;
    sUncert.scalForceFloor = sOpt.scalForceFloor;
    sUncert.vecForceRel = vecForceRel;              % [Nc x 1] f/<f>, every contact of the input packing
    sUncert.vecStiffRel = vecStiff / mean(vecStiff(vecBoolBackbone));   % [Nc x 1] k_n/<k_n>_backbone
    sUncert.vecBoolBackbone = vecBoolBackbone;      % [Nc x 1] both disks survive cleanRats
    sUncert.vecBoolWeak = vecBoolWeak;              % [Nc x 1] f < floor * <f>
    sUncert.scalNumContacts = numel(vecI);
    sUncert.scalNumContactsBackbone = scalNcB;
    sUncert.scalNumWeak = nnz(vecBoolWeak);
    sUncert.scalNumWeakBackbone = nnz(vecBoolWeak & vecBoolBackbone);
    sUncert.scalStiffShareWeak = sum(vecStiff(vecBoolWeak & vecBoolBackbone)) / sum(vecStiff(vecBoolBackbone));
    sUncert.scalNi = scalNi;
    sUncert.scalNiCut = scalNiCut;                  % backbone disks left if the weak contacts were cut
    sUncert.scalZ = scalZ;
    sUncert.scalZCut = scalZCut;
    sUncert.scalZIso = scalZIso;
    sUncert.scalDeltaZRel = (scalZ - scalZCut) / (scalZ - scalZIso);   % share of z - z_iso the cut removes
end

function vecBoolKeep = removeRattlers(vecI, vecJ, vecZWall, scalZMin)
% cleanRats' repeated removal on a fixed contact list.
    scalN = numel(vecZWall);
    vecBoolKeep = true(scalN, 1);
    while true
        vecBoolLive = vecBoolKeep(vecI) & vecBoolKeep(vecJ);
        vecZ = vecZWall;
        if any(vecBoolLive)
            vecZ = vecZ + accumarray([vecI(vecBoolLive); vecJ(vecBoolLive)], 1, [scalN 1]);
        end
        vecBoolDrop = vecBoolKeep & vecZ < scalZMin;
        if ~any(vecBoolDrop)
            break;
        end
        vecBoolKeep(vecBoolDrop) = false;
    end
end
