function [positions, radii] = cleanRats(positions, radii, Ly, Lx, Lz, boolFullyPeriodic)
% cleanRats(positions, radii, Ly, Lx)       — 2D: PBC in y only
% cleanRats(positions, radii, Ly, Lx, Lz)   — 3D: PBC in y and z
    if nargin < 6
        boolFullyPeriodic = false;
    end
    is3D = nargin >= 5 && ~isempty(Lz);


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

                %----------------- debug
                % if i <= 3 && j <= 4
                %     fprintf('  i=%d j=%d dx=%.4f dy=%.4f\n', i, j, dx, dy);
                % end
                %----------------- debug

                if is3D
                    dz = positions(i,3) - positions(j,3);
                    dz = dz - round(dz / Lz) * Lz;
                    r  = sqrt(dx^2 + dy^2 + dz^2);
                else
                    r  = sqrt(dx^2 + dy^2);
                end

                %----------------- debug
                % if i <= 3 && j <= 4
                %     fprintf('  i=%d j=%d r=%.6f ri+rj=%.6f overlap=%.6e\n', i, j, r, radii(i)+radii(j), radii(i)+radii(j)-r);
                % end
                %----------------- debug

                if radii(i) + radii(j) - r > -1e-3
                    Zn(i) = Zn(i) + 1;
                    Zn(j) = Zn(j) + 1;
                end
            end
        end

        fprintf('[cleanRats] Zn distribution: min=%d, max=%d, mean=%.2f\n', min(Zn), max(Zn), mean(Zn));
        % fprintf('[cleanRats] Zn values: ');
        % disp(Zn');
        if is3D
            to_keep = Zn > 3;
        else
            to_keep = Zn > 2;
        end

        changed = ~all(to_keep);
        positions = positions(to_keep, :);
        radii     = radii(to_keep);
        totalRattlers = totalRattlers + sum(~to_keep);
        fprintf('[cleanRats] Pass %d: %d rattlers removed.\n', passNumber, sum(~to_keep));
    end

    fprintf('[cleanRats] %d total rattlers removed.\n', totalRattlers);
    fprintf('[cleanRats] Total percentage of rattlers removed: %.2f%%\n', totalRattlers / N_original * 100);
end
