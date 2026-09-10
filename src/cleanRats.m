function [positions, radii] = cleanRats(positions, radii, Ly, Lx, Lz, boolFullyPeriodic, mu)
% cleanRats(positions, radii, Ly, Lx)                   — 2D, PBC in y only, frictionless
% cleanRats(positions, radii, Ly, Lx, Lz)              — 3D, PBC in y and z, frictionless
% cleanRats(positions, radii, Ly, Lx, Lz, boolFullyP)  — 3D + full PBC flag, frictionless
% cleanRats(...)                                        — mu>0: frictional variant.
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
% All existing callers that omit mu use mu <= 0, so the frictionless
% behavior (byte-for-byte identical). mu is the 7th positional argument;
% Lz= []  signals a 2D call (is3D = false).
    if nargin < 6
        boolFullyPeriodic = false;
    end
    if nargin < 7
        mu = 0;
    end
    is3D = nargin >= 5 && ~isempty(Lz);

    boolFriction = (mu > 0);
    if boolFriction
        fprintf('[cleanRats] mu=%.4f: frictional variant — keeping Zn>=2 backbone, removing Zn<=1 (floaters + one-contact rattlers)\n', mu);
    end

    changed = true;
    passNumber = 0;
    totalRattlers = 0;
    N_original = size(positions, 1);

    while changed
        passNumber    = passNumber + 1;
        fprintf('[cleanRats] Pass %d\n', passNumber);
        N             = size(positions, 1);

         % Recompute wall lists on current positions (only when NOT fully periodic)
        fprintf('[cleanRats] Recomputing wall lists...\n');
        left_wall_list   = positions(:,1) < radii;
        right_wall_list  = positions(:,1) > Lx - radii;

        Zn = zeros(N, 1);
        if ~boolFullyPeriodic
            Zn(left_wall_list | right_wall_list) = 2;
        end

        fprintf('[cleanRats] Computing contact lists...\n');
        for i = 1:N
            for j = i+1:N
                dx = positions(i,1) - positions(j,1);
                if boolFullyPeriodic
                    dx = dx - round(dx / Lx) * Lx;
                end
                dy = positions(i,2) - positions(j,2);
                dy = dy - round(dy / Ly) * Ly;

                if is3D
                    dz = positions(i,3) - positions(j,3);
                    dz = dz - round(dz / Lz) * Lz;
                    r   = sqrt(dx^2 + dy^2 + dz^2);
                else
                    r   = sqrt(dx^2 + dy^2);
                end

                if radii(i) + radii(j) - r > -1e-3
                    Zn(i) = Zn(i) + 1;
                    Zn(j) = Zn(j) + 1;
                end
            end
        end

        fprintf('[cleanRats] Zn distribution: min=%d, max=%d, mean=%.2f\n', min(Zn), max(Zn), mean(Zn));

        if boolFriction
             % Frictional 2D: remove Zn<=1 (floaters + one-contact rattlers).
             % In this gravity-free model a 1-contact disk cannot satisfy
             % force–torque balance, so it is a rattler even with friction.
             % Keep only Zn >= 2 grains (frictional 2D isostatic z_iso = 3).
            to_keep = Zn > 1;
        elseif is3D
            to_keep = Zn > 3;
        else
            to_keep = Zn > 2;
        end

        changed      = ~all(to_keep);
        positions     = positions(to_keep, :);
        radii         = radii(to_keep);
        totalRattlers = totalRattlers + sum(~to_keep);
        fprintf('[cleanRats] Pass %d: %d rattlers removed.\n', passNumber, sum(~to_keep));
    end

    fprintf('[cleanRats] %d total rattlers removed.\n', totalRattlers);
    fprintf('[cleanRats] Total percentage of rattlers removed: %.2f%%\n', totalRattlers / N_original * 100);
end
