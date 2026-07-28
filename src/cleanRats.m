function [positions, radii] = cleanRats(positions, radii, Ly, Lx, Lz)
% cleanRats(positions, radii, Ly, Lx)       — 2D: PBC in y only
% cleanRats(positions, radii, Ly, Lx, Lz)   — 3D: PBC in y and z

    is3D = nargin == 5;

    changed = true;
    passNumber = 0;
    totalRattlers = 0;
    N_original = size(positions, 1);

    while changed
        passNumber = passNumber + 1;
        fprintf('[cleanRats] Pass %d\n', passNumber);
        N = size(positions, 1);

        % Recompute wall lists on current positions
        fprintf('[cleanRats] Recomputing wall lists...\n');
        left_wall_list  = positions(:,1) < radii;
        right_wall_list = positions(:,1) > Lx - radii;

        Zn = zeros(N, 1);
        Zn(left_wall_list | right_wall_list) = 2; % wall particles get 2 contacts so they aren't removed

        fprintf('[cleanRats] Computing contact lists...\n');
        for i = 1:N
            for j = i+1:N
                dx = positions(i,1) - positions(j,1);
                dy = positions(i,2) - positions(j,2);
                dy = dy - round(dy / Ly) * Ly;
                if is3D
                    dz = positions(i,3) - positions(j,3);
                    dz = dz - round(dz / Lz) * Lz;
                    r  = sqrt(dx^2 + dy^2 + dz^2);
                else
                    r  = sqrt(dx^2 + dy^2);
                end
                if radii(i) + radii(j) - r > 0
                    Zn(i) = Zn(i) + 1;
                    Zn(j) = Zn(j) + 1;
                end
            end
        end

        to_keep = Zn > 2;
        changed = ~all(to_keep);
        positions = positions(to_keep, :);
        radii     = radii(to_keep);
        totalRattlers = totalRattlers + sum(~to_keep);
        fprintf('[cleanRats] Pass %d: %d rattlers removed.\n', passNumber, sum(~to_keep));
    end

    fprintf('[cleanRats] %d total rattlers removed.\n', totalRattlers);
    fprintf('[cleanRats] Total percentage of rattlers removed: %.2f%%\n', totalRattlers / N_original * 100);
end



