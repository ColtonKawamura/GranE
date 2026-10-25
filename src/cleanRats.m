function [matPositions, vecRadii, vecBoolKeep] = cleanRats(matPositions, vecRadii, scalBoxHeightY, scalBoxWidthX, scalBoxDepthZ, boolFullyPeriodic, scalMu)
% cleanRats(pos, radii, Ly, Lx)                   — 2D, PBC in y only, frictionless
% cleanRats(pos, radii, Ly, Lx, Lz)              — 3D, PBC in y and z, frictionless
% cleanRats(pos, radii, Ly, Lx, Lz, boolFullyP)  — + full PBC flag (Lz = [] for 2D)
% cleanRats(..., mu)                              — mu > 0: frictional variant.
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
    if nargin < 6
        boolFullyPeriodic = false;
    end
    if nargin < 7
        scalMu = 0;
    end
    boolThreeD = nargin >= 5 && ~isempty(scalBoxDepthZ);

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

    fprintf('[cleanRats] %d total rattlers removed.\n', scalTotalRattlers);
    fprintf('[cleanRats] Total percentage of rattlers removed: %.2f%%\n', scalTotalRattlers / scalNumOriginal * 100);
end
