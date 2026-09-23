function simMD(K, M, Bv, w_D, N, P, W, seed, in_path, out_path, options)
%   simMD(100, 1, 1, 1.28, 80000, 0.1, 40, 1, "in/2d_poly_20by20/", "out/junk_yard")

    arguments
        K (1,1) double
        M (1,1) double
        Bv (1,1) double
        w_D (1,1) double
        N (1,1) double
        P (1,1) double
        W (1,1) double
        seed (1,1) double
        in_path (1,1) string
        out_path (1,1) string
        options.visSim          (1,1) logical = false
        options.shear           (1,1) logical = false
        options.cleanRats       (1,1) logical = false
        options.plotProbes      (1,1) logical = false
        options.maxAmpTracking  (1,1) logical = false
        options.fullSpectrum    (1,1) logical = false
        options.plotSpectrum    (1,1) logical = false
    end

%% Logging
    fprintf('[simMD] START  K=%g M=%g Bv=%g w_D=%g N=%d P=%g W=%d seed=%d\n', K, M, Bv, w_D, N, P, W, seed);
    fprintf('[simMD]   in_path  = %s\n', in_path);
    fprintf('[simMD]   out_path = %s\n', out_path);
    simMD_tStart = tic;

try

%% GPU setup
    useGPU = canUseGPU();
    if useGPU
        gdev = gpuDevice();
        fprintf('[simMD] GPU detected: %s (%.1f GB)\n', gdev.Name, gdev.TotalMemory/1024^3);
    else
        fprintf('[simMD] No GPU found — running on CPU.\n');
    end

%% Check for previous output
    packing_name    = string(sprintf("2D_N%d_P%s_Width%d_Seed%d", N, num2str(P), W, seed));
    filename_output = string(sprintf("%s_K%d_Bv%d_wD%.2f_M%d.mat", packing_name, K, Bv, w_D, M));
    save_path       = fullfile(out_path, filename_output);

    if exist(char(out_path + filename_output), 'file')
        display("output already exists: " + out_path + filename_output);
        return
    end
    input_pressure = P;
    filename = in_path + packing_name + ".mat";

    fprintf('[simMD] Loading packing: %s\n', filename);
    load(filename);

    %% Normalize input to the NEW data-file format
    %   Two on-disk formats exist:
    %     OLD: x, y, Dn (row vectors), Lx, Ly, P
    %     NEW: vecPosX, vecPosY, vecDiameter (column vectors),
    %          scalBoxWidthX, scalBoxHeightY, scalPressure
    %   The NEW format is preferred. If the file is OLD (detected by the
    %   presence of 'x'), convert its variables into the NEW names/shapes so
    %   the rest of the function works uniformly on the NEW format.
    if exist('x', 'var')
        fprintf('[simMD] OLD input format (x/y/Dn) detected; converting to NEW format names.\n');
        vecPosX        = x(:);
        vecPosY        = y(:);
        vecDiameter    = Dn(:);
        scalBoxWidthX  = Lx;
        scalBoxHeightY = Ly;
        if ~exist('scalPressure', 'var')
            scalPressure = P;
        end
    elseif exist('vecPosX', 'var')
        fprintf('[simMD] NEW input format (vecPosX) detected.\n');
    else
        error('simMD:badInput', ...
            'Input file %s contains neither OLD (x) nor NEW (vecPosX) format variables.', ...
            char(filename));
    end

    % The particle arrays are the source of truth for the count.
    N = length(vecPosX);

    fprintf('[simMD] Packing loaded (%d particles).\n', N);

%% Clean rattlers
    if options.cleanRats
        positions = [vecPosX, vecPosY];
        radii = vecDiameter/2;
        [positions, radii] = cleanRats(positions, radii, scalBoxHeightY, scalBoxWidthX);
        vecPosX = positions(:,1);
        vecPosY = positions(:,2);
        vecDiameter = 2*radii;
        N  = length(vecPosX);
        display("cleaning rattlers");
    end

    % Enforce column vectors throughout — critical for GPU indexing
    vecPosX     = vecPosX(:);
    vecPosY     = vecPosY(:);
    vecDiameter = vecDiameter(:);

    % update masses based on diameter
    mass = (pi/4) .* vecDiameter.^2;
    inv_mass = 1 ./ mass;
    mass_particle_average = mean(mass);

%% Simulation parameters
    B   = 0;
    A   = P_target/100;
    % dt  = pi*sqrt(M/K)*0.005;
    m_min = min(mass); 
    dt = pi*sqrt(m_min/K)*0.005;
    c_0 = min(vecDiameter)*sqrt(K/m_min);
    Nt  = round(.8*(scalBoxWidthX/c_0)/dt);

    fprintf('[simMD] Nt = %d, dt = %g\n', Nt, dt);

    % Pre-compute scalar constants (avoids repeated division inside loop)
    % inv_M    = 1 / M;
    dt2_half = dt^2 / 2;
    dt_half  = dt / 2;

%% Memory estimate
    traj_bytes = double(N) * double(Nt) * 8;
    fprintf('[simMD] Avoided trajectory arrays: %.2f GB each, %.2f GB total.\n', ...
        traj_bytes/1024^3, 2*traj_bytes/1024^3);

%% Wall / bulk masks
    if seed == 0
        left_wall_list  = (vecPosX < vecDiameter/1);
    else
        left_wall_list  = (vecPosX < vecDiameter/2);
    end
    right_wall_list = (vecPosX > scalBoxWidthX - vecDiameter/2);
    bulk_list       = ~(left_wall_list | right_wall_list);  %#ok<NASGU>

    % Index vectors (faster than logical masks for GPU indexed assignment)
    left_wall_idx  = find(left_wall_list);
    right_wall_idx = find(right_wall_list);

%% Build neighbor lists (CPU, done once)
    fprintf('[simMD] Building neighbor lists ...\n');
    skin = 0;
    Zn_list           = zeros(N, 1);
    neighbor_list_all = cell(1, N);
    spring_list_all   = cell(1, N);

    for nn = 1:N
        neighbor_list_nn = [];
        spring_list_nn   = [];
        for mm = [1:nn-1, nn+1:N]
            dvecY  = vecPosY(mm) - vecPosY(nn);
            dvecY  = dvecY - round(dvecY/scalBoxHeightY)*scalBoxHeightY;
            Dnm = (1+skin)*(vecDiameter(nn)+vecDiameter(mm))/2;
            if abs(dvecY) <= Dnm
                dvecX  = vecPosX(mm) - vecPosX(nn);
                dnm = dvecX^2 + dvecY^2;
                if dnm < Dnm^2
                    neighbor_list_nn = [neighbor_list_nn, mm];         %#ok<AGROW>
                    spring_list_nn   = [spring_list_nn,   sqrt(dnm)];  %#ok<AGROW>
                end
            end
        end
        neighbor_list_all{nn} = neighbor_list_nn;
        spring_list_all{nn}   = spring_list_nn;
        Zn_list(nn)           = length(spring_list_nn);
    end
    fprintf('[simMD] Neighbor lists built.\n');

%% Convert cell lists → flat edge vectors (enables full GPU vectorization)
    % -----------------------------------------------------------------
    % src_flat(e), dst_flat(e), D0_flat(e) describe directed edge e.
    % Each undirected pair (i,j) appears TWICE: as (i→j) and (j→i).
    % This means accumarray scatter-add works without any special
    % symmetry handling — each particle accumulates its own forces.
    % Ep is divided by 2 at the end to correct for double-counting.
    % -----------------------------------------------------------------
    max_nb  = max(Zn_list);
    nbr_idx = zeros(N, max_nb, 'int32');
    spr_len = zeros(N, max_nb);

    for nn = 1:N
        nl = neighbor_list_all{nn};
        sl = spring_list_all{nn};
        nb = length(nl);
        nbr_idx(nn, 1:nb) = int32(nl);
        spr_len(nn, 1:nb) = sl;
    end
    clear neighbor_list_all spring_list_all;  % free cell memory

    valid_mask = (nbr_idx > 0);
    src_rows   = repmat((1:N)', 1, max_nb);
    src_flat   = src_rows(valid_mask);          % E×1
    dst_flat   = double(nbr_idx(valid_mask));   % E×1
    D0_flat    = spr_len(valid_mask);           % E×1
    clear nbr_idx spr_len valid_mask src_rows;

    fprintf('[simMD] Edge list: %d directed edges (avg %.1f neighbors/particle).\n', ...
        length(src_flat), length(src_flat)/N);

% ---- CONTINUES IN CHUNK 2 ----
% ---- CONTINUES FROM CHUNK 1 ----

%% Verlet state arrays
    vecPosX0 = vecPosX;   vecPosY0 = vecPosY;
    ax_old = zeros(N, 1);
    ay_old = zeros(N, 1);
    vx     = zeros(N, 1);
    vy     = zeros(N, 1);
    Ek     = zeros(Nt, 1);
    Ep     = zeros(Nt, 1);
    g      = 0;

%% Full-spectrum trajectory storage
    if options.fullSpectrum
        traj_bytes = double(N) * double(Nt) * 8; % double = 64 bits = 64/8 bytes
        fprintf('[simMD] fullSpectrum enabled — allocating %.2f GB per trajectory array.\n', ...
            traj_bytes/1024^3); % 1024 GB/ bytes
        traj_disp_x = zeros(N, Nt);
        traj_disp_y = zeros(N, Nt);
    end

%% DFT accumulators
    fprintf('[simMD] Integer-cycle movement-triggered on-the-fly DFT at exact driving freq %g Hz.\n', w_D/(2*pi));

    dft_x           = complex(zeros(N, 1));
    dft_y           = complex(zeros(N, 1));
    n_dft_samples   = zeros(N, 1);
    running_sum_x   = zeros(N, 1);
    running_sum_y   = zeros(N, 1);
    running_sumsq_x = zeros(N, 1);
    running_sumsq_y = zeros(N, 1);

%% DFT trigger parameters
    T_period_int             = round(2*pi / (w_D * dt));
    movementTriggerThreshold = 1e-6;
    dft_delay_cycles         = 3;
    arrival_safety_factor    = 1.2;

    dft_active   = false(N, 1);
    nt_dft_start = zeros(N, 1);
    nt_dft_on    = zeros(N, 1);
    nt_dft_end   = zeros(N, 1);

    % Wall particles: driven sinusoid from t=0, no onset delay needed
    dft_active(left_wall_list)   = true;
    nt_dft_start(left_wall_list) = 1;
    nt_dft_on(left_wall_list)    = 1;
    n_complete_wall              = floor(Nt / T_period_int);
    nt_dft_end(left_wall_list)   = n_complete_wall * T_period_int;

    expected_arrival_nt = floor((vecPosX0 - min(vecPosX0(left_wall_list))) / c_0 / dt * arrival_safety_factor);
    expected_arrival_nt = max(expected_arrival_nt, 1);

    fprintf('[simMD] T_period_int = %d samples, c_0 = %g.\n', T_period_int, c_0);
    fprintf('[simMD] arrival_safety_factor = %.2f, movementTriggerThreshold = %g, dft_delay_cycles = %d.\n', ...
        arrival_safety_factor, movementTriggerThreshold, dft_delay_cycles);

%% Max-amplitude tracking init
    if options.maxAmpTracking
        maxAmpXAfter2 = zeros(N, 1);
        maxAmpYAfter2 = zeros(N, 1);
        twoCycles     = 2 * T_period_int;
    end

%% Probe particles
    [~, idx] = sort(vecPosX0);
    if options.plotProbes
        probe_targets = [20, 100, 300, 400];
        probe_idx     = zeros(1, 4);
        for pp = 1:4
            [~, probe_idx(pp)] = min(abs(vecPosX0 - probe_targets(pp)));
        end
        probe_traj_x = zeros(4, Nt);
        probe_traj_y = zeros(4, Nt);
        fprintf('[simMD] Probe particles at x0 = [%.2f, %.2f, %.2f, %.2f].\n', ...
            vecPosX0(probe_idx(1)), vecPosX0(probe_idx(2)), vecPosX0(probe_idx(3)), vecPosX0(probe_idx(4)));
    else
        fprintf('[simMD] No probe particles being recorded.\n');
    end

%% GPU transfer (everything hot goes to GPU once, stays there until after the loop)
    if useGPU
        % Particle state
        vecPosX     = gpuArray(vecPosX);
        vecPosY     = gpuArray(vecPosY);
        vecPosX0_g  = gpuArray(vecPosX0);
        vecPosY0_g  = gpuArray(vecPosY0);
        vx     = gpuArray(vx);
        vy     = gpuArray(vy);
        ax_old = gpuArray(ax_old);
        ay_old = gpuArray(ay_old);
        Ek     = gpuArray(Ek);
        Ep     = gpuArray(Ep);
        inv_mass = gpuArray(inv_mass);
        mass_g = gpuArray(mass);

        % Edge list — transferred once, read every timestep
        src_flat = gpuArray(src_flat);
        dst_flat = gpuArray(dst_flat);
        D0_flat  = gpuArray(D0_flat);

        % DFT accumulators
        dft_x           = gpuArray(dft_x);
        dft_y           = gpuArray(dft_y);
        n_dft_samples   = gpuArray(n_dft_samples);
        running_sum_x   = gpuArray(running_sum_x);
        running_sum_y   = gpuArray(running_sum_y);
        running_sumsq_x = gpuArray(running_sumsq_x);
        running_sumsq_y = gpuArray(running_sumsq_y);

        % DFT trigger arrays
        dft_active          = gpuArray(dft_active);
        nt_dft_start        = gpuArray(nt_dft_start);
        nt_dft_on           = gpuArray(nt_dft_on);
        nt_dft_end          = gpuArray(nt_dft_end);
        expected_arrival_nt = gpuArray(expected_arrival_nt);

        % Wall index lists
        left_wall_idx_g  = gpuArray(left_wall_idx);
        right_wall_idx_g = gpuArray(right_wall_idx);

        % Max-amp tracking
        if options.maxAmpTracking
            maxAmpXAfter2 = gpuArray(maxAmpXAfter2);
            maxAmpYAfter2 = gpuArray(maxAmpYAfter2);
        end

        fprintf('[simMD] Arrays transferred to GPU.\n');
    else
        % CPU path: just alias names so the loop code is identical
        mass_g = mass;
        vecPosX0_g     = vecPosX0;
        vecPosY0_g     = vecPosY0;
        left_wall_idx_g  = left_wall_idx;
        right_wall_idx_g = right_wall_idx;
    end


%% Main Verlet loop
    fprintf('[simMD] Starting main loop (%d timesteps) ...\n', Nt);
    simMD_loop_tic     = tic;
    simMD_log_interval = max(1, floor(Nt/10));

    for nt = 1:Nt

        if mod(nt, simMD_log_interval) == 0
            fprintf('[simMD]   timestep %d / %d  (%.1f%%, elapsed %.1f s)\n', ...
                nt, Nt, 100*nt/Nt, toc(simMD_loop_tic));
        end

        if options.visSim
            % gather only what visualizer needs — minimal CPU↔GPU transfer
            visualizeSim(gather(vecPosX), vecPosX0, gather(vecPosY), vecPosY0, idx, A, Bv, P, w_D, K, true);
        end

        % ── Verlet step 1: update positions ──────────────────────────────
        vecPosX = vecPosX + vx.*dt + ax_old.*dt2_half;
        vecPosY = vecPosY + vy.*dt + ay_old.*dt2_half;

        % ── Forced wall displacements ─────────────────────────────────────
        t_now = nt * dt;
        if options.shear
            vecPosX(left_wall_idx_g) = vecPosX0_g(left_wall_idx_g);
            vecPosY(left_wall_idx_g) = vecPosY0_g(left_wall_idx_g) + A*sin(w_D*t_now);
        else
            vecPosX(left_wall_idx_g) = vecPosX0_g(left_wall_idx_g) + A*sin(w_D*t_now);
            vecPosY(left_wall_idx_g) = vecPosY0_g(left_wall_idx_g);
        end
        vecPosX(right_wall_idx_g) = vecPosX0_g(right_wall_idx_g);
        vecPosY(right_wall_idx_g) = vecPosY0_g(right_wall_idx_g);

        % ── Vectorized GPU force kernel ───────────────────────────────────
        %
        %  Replaces the entire nested for nn / for mm_counter loop.
        %
        %  For each directed edge e = (src, dst) with equilibrium length D0:
        %    dx, dy   = position vector src → dst  (with PBC in y)
        %    dnm      = current distance
        %    Fmag     = -K*(D0/dnm - 1)            harmonic spring
        %    Fx_edge  = Fmag*dx - Bv*(vx_src-vx_dst)
        %    Fy_edge  = Fmag*dy - Bv*(vy_src-vy_dst)
        %
        %  accumarray scatter-adds all edge contributions onto each particle.
        %  All E edges computed in parallel on GPU — zero serial loops.
        % -----------------------------------------------------------------

        % Gather src and dst positions/velocities (single indexed read each)
        xs  = vecPosX(src_flat);    ys  = vecPosY(src_flat);
        xd  = vecPosX(dst_flat);    yd  = vecPosY(dst_flat);
        vxs = vx(src_flat);   vys = vy(src_flat);
        vxd = vx(dst_flat);   vyd = vy(dst_flat);

        % Displacement vector with periodic boundary in y
        dxv = xd - xs;
        dyv = yd - ys;
        dyv = dyv - round(dyv./scalBoxHeightY).*scalBoxHeightY;

        % Distance and spring force magnitude
        dnm  = sqrt(dxv.^2 + dyv.^2);
        Fmag = -K .* (D0_flat./dnm - 1);

        % Damping
        dvxv = vxs - vxd;
        dvyv = vys - vyd;

        % Per-edge force and potential energy
        Fx_edge = Fmag.*dxv - Bv.*dvxv;
        Fy_edge = Fmag.*dyv - Bv.*dvyv;
        Ep_edge = 0.5*K.*(D0_flat - dnm).^2;

        % Scatter-add: sum edge contributions per source particle
        Fx = accumarray(src_flat, Fx_edge, [N, 1]);
        Fy = accumarray(src_flat, Fy_edge, [N, 1]);

        % Global background damping (B=0 by default, kept for generality)
        Fx = Fx - B.*vx;
        Fy = Fy - B.*vy;

        % Wall particles are kinematically driven — zero their forces
        Fx(left_wall_idx_g)  = 0;
        Fx(right_wall_idx_g) = 0;

        % Energy (Ep /2 because each undirected pair appears twice in edge list)
        Ep(nt) = sum(Ep_edge) / (2*N);
        % Ek(nt) = (0.5*M * sum(vx.^2 + vy.^2)) / N;
        Ek(nt) = (0.5 * sum(mass_g .* (vx.^2 + vy.^2))) / N;

        % ── Accelerations ────────────────────────────────────────────────
        % ax = Fx .* inv_M;
        % ay = Fy .* inv_M - g;
        ax = Fx .* inv_mass;
        ay = Fy .* inv_mass - g;

        % ── Verlet step 2: update velocities ─────────────────────────────
        vx = vx + (ax_old + ax) .* dt_half;
        vy = vy + (ay_old + ay) .* dt_half;

        ax_old = ax;
        ay_old = ay;

        % ── DFT accumulation ─────────────────────────────────────────────
        % if GPU is used,  vecPosX and vecPosX0_g were transfered to the GPU
        % so disp_x is created on the GPU
        disp_x = vecPosX - vecPosX0_g; % [N x 1]
        disp_y = vecPosY - vecPosY0_g;

        % ── Full-spectrum trajectory recording  ─────────────────
          %       traj_disp_x  =  N × Nt  matrix
          %
          %                nt=1  nt=2  nt=3 ...
          % particle 1  [  u1    u1    u1  ... ]
          % particle 2  [  u2    u2    u2  ... ]
          % ...
          % particle N  [  uN    uN    uN  ... ]
        if options.fullSpectrum
            if useGPU
                traj_disp_x(:, nt) = gather(disp_x); % [N x nt]
                traj_disp_y(:, nt) = gather(disp_y);
            else
                traj_disp_x(:, nt) = disp_x; 
                traj_disp_y(:, nt) = disp_y;
            end
        end

        % Movement trigger (time-gated to block soliton)
        if options.shear
            newly_active = (~dft_active) ...
                & (nt > expected_arrival_nt) ...
                & (abs(disp_y) > movementTriggerThreshold);
        else
            newly_active = (~dft_active) ...
                & (nt > expected_arrival_nt) ...
                & (abs(disp_x) > movementTriggerThreshold);
        end

        if any(newly_active)
            nt_dft_start(newly_active) = nt;
            delay_nt = nt + dft_delay_cycles * T_period_int;
            nt_dft_on(newly_active)    = delay_nt;
            n_complete_new             = max(0, floor((Nt - delay_nt + 1) / T_period_int));
            nt_dft_end(newly_active)   = delay_nt + n_complete_new * T_period_int - 1;
            dft_active(newly_active)   = true;
        end

        active_mask = double(dft_active & (nt >= nt_dft_on) & (nt <= nt_dft_end));
        twiddle     = exp(-1i * w_D * (nt-1) * dt);   % scalar, no GPU overhead

        dft_x           = dft_x           + disp_x .* active_mask .* twiddle;
        dft_y           = dft_y           + disp_y .* active_mask .* twiddle;
        n_dft_samples   = n_dft_samples   + active_mask;
        running_sum_x   = running_sum_x   + vecPosX .* active_mask;
        running_sum_y   = running_sum_y   + vecPosY .* active_mask;
        running_sumsq_x = running_sumsq_x + vecPosX.^2 .* active_mask;
        running_sumsq_y = running_sumsq_y + vecPosY.^2 .* active_mask;

        % ── Max-amplitude tracking (optional) ────────────────────────────
        if options.maxAmpTracking
            afterTwoCycles = dft_active & (nt >= nt_dft_start + twoCycles);
            absDispX = abs(disp_x);
            absDispY = abs(disp_y);
		replaceX = afterTwoCycles & (absDispX > maxAmpXAfter2);
            replaceY = afterTwoCycles & (absDispY > maxAmpYAfter2);
            replaceX_idx = find(replaceX);
            if ~isempty(replaceX_idx)
                maxAmpXAfter2(replaceX_idx) = absDispX(replaceX_idx);
            end
            replaceY_idx = find(replaceY);
            if ~isempty(replaceY_idx)
                maxAmpYAfter2(replaceY_idx) = absDispY(replaceY_idx);
            end
        end

        % ── Probe recording (optional, minimal gather) ────────────────────
        if options.plotProbes
            probe_traj_x(:, nt) = gather(vecPosX(probe_idx(:)));
            probe_traj_y(:, nt) = gather(vecPosY(probe_idx(:)));
        end

    end % nt loop

    fprintf('[simMD] Main loop finished (%.1f s).\n', toc(simMD_loop_tic));

%% Gather GPU arrays back to CPU
    if useGPU
        vecPosX         = gather(vecPosX);
        vecPosY         = gather(vecPosY);
        vecPosX0        = gather(vecPosX0_g);
        vecPosY0        = gather(vecPosY0_g);
        vx              = gather(vx);
        vy              = gather(vy);
        dft_x           = gather(dft_x);
        dft_y           = gather(dft_y);
        n_dft_samples   = gather(n_dft_samples);
        running_sum_x   = gather(running_sum_x);
        running_sum_y   = gather(running_sum_y);
        running_sumsq_x = gather(running_sumsq_x);
        running_sumsq_y = gather(running_sumsq_y);
        dft_active      = gather(dft_active);
        nt_dft_start    = gather(nt_dft_start);
        Ek              = gather(Ek);
        Ep              = gather(Ep);
        if options.maxAmpTracking
            maxAmpXAfter2 = gather(maxAmpXAfter2);
            maxAmpYAfter2 = gather(maxAmpYAfter2);
        end
        fprintf('[simMD] GPU arrays gathered to CPU.\n');
    end


%% Probe plotting
    if options.plotProbes
        time_vector_full = (1:Nt) * dt;
        plotProbes(time_vector_full, probe_traj_x, probe_traj_y, vecPosX0(probe_idx(:)));
    end


%% DFT normalization
    safe_n = max(n_dft_samples, 1); % clamp to 1 to avoid dividing by 0 for non-active particles
    dft_x  = dft_x ./ safe_n; % we zero-them out anyways in dft_x(~valid) = 0
    dft_y  = dft_y ./ safe_n;
    var_x  = running_sumsq_x ./ safe_n - (running_sum_x ./ safe_n).^2;
    var_y  = running_sumsq_y ./ safe_n - (running_sum_y ./ safe_n).^2;
    clear running_sum_x running_sum_y running_sumsq_x running_sumsq_y;

    min_cycles = 5;
    % T_period_int = timesteps / period
    valid      = n_dft_samples >= min_cycles * T_period_int;
    dft_x(~valid) = 0; 
    dft_y(~valid) = 0;

    fprintf('[simMD] DFT normalized; active: %d / %d, passing min_cycles=%d: %d / %d.\n', ...
        sum(dft_active), N, min_cycles, sum(valid), N);
    clear n_dft_samples safe_n valid;

%% Full frequency spectrum 
    if options.fullSpectrum
        fprintf('[simMD] Computing full FFT spectrum for all particles ...\n');

        Nfft = Nt;

        % get all possible frequencies from 0 to nyquist
        % this is just integer multipes  * 1/(N * deltaX)
        % which is all possible frequncies that fit into Nt*deltaX
        % create vector from largest wavenumber (lowest freq) possible in Nt
        % lowest freq is one oscillation in (Nt*dt)
        % higst is Nyquist = 1/2 * 1/dt
        freq_axis = (0:floor(Nfft/2)) / (Nfft * dt);   % 1/time, one-side

        % FFT along time dimension, normalizd by Nfft
        % "2" argument applies fft() to columns (time) for a single paticle
        % since traj_disp_x is [N x Nt]
        F_x = fft(traj_disp_x, Nfft, 2) / Nfft; % [N x Nfft] double sided!!
        F_y = fft(traj_disp_y, Nfft, 2) / Nfft; % so [particle_n x it's fft pectrum]

        % Single-sided amplitude: 2*|X_k|
        % keep all rows, but only 1 to half of frequency spectrum
        % take abs and dbouel to give the ampltiude back
        spec_x = 2 * abs(F_x(:, 1:floor(Nfft/2)+1)); % +1 because matlab is 1-based
        spec_y = 2 * abs(F_y(:, 1:floor(Nfft/2)+1));

        % Correct DC (no factor of 2)
        % because that's where the fold happens so they dont'
        % have a mirror "peak" on double sided spectrum
        spec_x(:, 1)   = abs(F_x(:, 1));
        spec_y(:, 1)   = abs(F_y(:, 1));

        % if Nt = Nfft is even, there is a bin at exactly Nyquist
        % If so, then there will be a fold there and no mirror
        % if odd, then there will be no fold and don't need to corret
        if mod(Nfft, 2) == 0
            spec_x(:, end) = abs(F_x(:, floor(Nfft/2)+1));
            spec_y(:, end) = abs(F_y(:, floor(Nfft/2)+1));
        end

        fprintf('[simMD] Full spectrum computed: %d particles x %d freq bins (df = %.6f Hz).\n', ...
            N, numel(freq_axis), freq_axis(2));

        % DEBUG
        fprintf('max traj_disp_x = %e\n', max(traj_disp_x(:)));
        fprintf('max traj_disp_y = %e\n', max(traj_disp_y(:)))

        clear traj_disp_x traj_disp_y F_x F_y;
    end

%% Distance–frequency spectrum plot (optional)
    if options.plotSpectrum && options.fullSpectrum
        fprintf('[simMD] Plotting distance-frequency spectrum ...\n');

        initial_distance = vecPosX0(:);

        [dist_sorted, sort_idx] = sort(initial_distance);

        % sort the spectrum rows by initial distances
        spec_x_sorted = spec_x(sort_idx, :);
        spec_y_sorted = spec_y(sort_idx, :);

        % only plot 10 times driving frequency
        omega_axis = 2*pi * freq_axis;
        omega_mask = omega_axis <= 10 * w_D;

        % X displacement spectrum
        figure;
        imagesc(dist_sorted,...
            omega_axis(omega_mask), ...
            log10(spec_x_sorted(:, omega_mask).' + 1e-20));
        axis xy;
        clim([-8, -3]); % only plot from 10^-8 to 10^-3
        xlabel('Initial distance from oscillating wall',...
            'Interpreter', 'latex');
        ylabel('$\omega_x$',...
            'Interpreter', 'latex');

        title_str = sprintf('$p=%.2f$, $\\sqrt{p}=%.2f$, $w_d=%.2f$, $\\gamma=%.4f$, $L_x=%.2f$, $L_y=%.2f$', ...
            P, ...
            sqrt(P), ...
            w_D,...
            Bv,... 
            scalBoxWidthX,...
            scalBoxHeightY);

        title(title_str,...
            'Interpreter', 'latex');

        cmap = [linspace(1,0,256)',...
            zeros(256,1),...
            linspace(0,1,256)'];
        colormap(flipud(cmap));
        cb = colorbar;
        ylabel(cb,...
            '$\log_{10}$ amplitude',...
            'Interpreter', 'latex');
        set(gca,...
            'TickLabelInterpreter', 'latex');
        exportgraphics(gcf, fullfile(out_path, sprintf('%s_spectrum_x.png', packing_name)), 'Resolution', 300);

        % Y displacement spectrum
        figure;
        imagesc(dist_sorted,...
            omega_axis(omega_mask), ...
            log10(spec_y_sorted(:, omega_mask).' + 1e-20));
        axis xy;
        clim([-8, -3]);
        xlabel('Initial distance from oscillating wall',...
            'Interpreter', 'latex');
        ylabel('$\omega_y$',...
            'Interpreter', 'latex');
        title(title_str,...
            'Interpreter', 'latex');
        colormap(flipud(cmap));
        cb = colorbar;
        ylabel(cb,...
            '$\log_{10}$ amplitude',...
            'Interpreter', 'latex');
        set(gca, ...
            'TickLabelInterpreter', 'latex');
        exportgraphics(gcf, fullfile(out_path, sprintf('%s_spectrum_y.png', packing_name)), 'Resolution', 300);
    end

%% Post-processing
    fprintf('[simMD] Starting post-processing ...\n');
    addpath('./src/matlab_functions');

    % ── X direction ──────────────────────────────────────────────────────
    fprintf('[simMD] DFT processing X direction ...\n');
    [~, index_particles]              = sort(vecPosX0);
    index_oscillating_wall            = left_wall_list;
    driving_frequency                 = w_D / (2*pi);
    driving_amplitude                 = A;
    initial_distance_from_oscillation = vecPosX0;

    [fitted_attenuation, wavenumber, attenuation_fit_line, ...
     initial_distance_from_oscillation_output, amplitude_vector, ...
     unwrapped_phase_vector, cleaned_particle_index, ...
     x_fft_initial_y, x_fft_initial_z] = ...
        processDFT(dft_x, var_x, driving_amplitude, index_particles, ...
                       index_oscillating_wall, initial_distance_from_oscillation, vecPosY0, vecPosY0);

    if isempty(cleaned_particle_index) && ~options.shear
        fprintf('Simulation P=%d, Omega=%d, Gamma=%d, Seed=%d did not detect attenuation\n', P, w_D, Bv, seed);
        if options.fullSpectrum
            save(save_path, 'freq_axis', 'spec_x', 'spec_y', '-append', '-v7.3');
            fprintf('[simMD] Saved spectrum before early exit.\n');
        end
        return
    end

    attenuation_x        = fitted_attenuation;
    attenuation_fit_line_x = attenuation_fit_line;
    wavenumber_x         = wavenumber;
    unwrapped_phase_vector_x = unwrapped_phase_vector;
    wavespeed_x          = driving_frequency*2*pi*sqrt(mass_particle_average/K) / wavenumber;
    initial_distance_from_oscillation_output_x_fft = initial_distance_from_oscillation_output;
    amplitude_vector_x   = amplitude_vector;
    cleaned_particle_index_x = cleaned_particle_index;

    % ── Y direction ──────────────────────────────────────────────────────
    fprintf('[simMD] DFT processing Y direction ...\n');
    [~, index_particles]              = sort(vecPosY0);
    initial_distance_from_oscillation = vecPosX0;

    [fitted_attenuation, wavenumber, attenuation_fit_line, ...
     initial_distance_from_oscillation_output, amplitude_vector, ...
     unwrapped_phase_vector, cleaned_particle_index, ...
     y_fft_initial_y, y_fft_initial_z] = ...
        processDFT(dft_y, var_y, driving_amplitude, index_particles, ...
                       index_oscillating_wall, initial_distance_from_oscillation, vecPosY0, vecPosY0);

    if isempty(cleaned_particle_index)
        fprintf('Simulation P=%d, Omega=%d, Gamma=%d, Seed=%d did not detect attenuation\n for y-direction', P, w_D, Bv, seed);
        if options.fullSpectrum
            save(save_path, 'freq_axis', 'spec_x', 'spec_y');
            fprintf('[simMD] Saved spectrum before early exit.\n');
        end
        return
    end

    attenuation_y        = fitted_attenuation;
    attenuation_fit_line_y = attenuation_fit_line;
    wavenumber_y         = wavenumber;
    unwrapped_phase_vector_y = unwrapped_phase_vector;
    wavespeed_y          = driving_frequency*2*pi*sqrt(mass_particle_average/K) / wavenumber;
    initial_distance_from_oscillation_output_y_fft = initial_distance_from_oscillation_output;
    amplitude_vector_y   = amplitude_vector;
    cleaned_particle_index_y = cleaned_particle_index;

%% Max-amplitude plots (optional)
    if options.maxAmpTracking
        vecPosX0_col = vecPosX0(:);

        % X direction fit
        valid_x = maxAmpXAfter2(:) > 0;
        if any(valid_x)
            x0_max_x   = max(vecPosX0_col(valid_x));
            fit_mask_x = valid_x & (vecPosX0_col <= 0.75*x0_max_x);
            if sum(fit_mask_x) >= 2
                p_x = polyfit(vecPosX0_col(fit_mask_x), log(maxAmpXAfter2(fit_mask_x)), 1);
                attenuationMaxAmpX = p_x(1);
                fit_line_maxX      = exp(polyval(p_x, vecPosX0_col(fit_mask_x)));
            else
                attenuationMaxAmpX = NaN;  fit_mask_x = false(size(vecPosX0_col));  fit_line_maxX = [];
            end
        else
            attenuationMaxAmpX = NaN;  fit_mask_x = false(size(vecPosX0_col));  fit_line_maxX = [];
        end

        % Y direction fit
        valid_y = maxAmpYAfter2(:) > 0;
        if any(valid_y)
            x0_max_y   = max(vecPosX0_col(valid_y));
            fit_mask_y = valid_y & (vecPosX0_col <= 0.75*x0_max_y);
            if sum(fit_mask_y) >= 2
                p_y = polyfit(vecPosX0_col(fit_mask_y), log(maxAmpYAfter2(fit_mask_y)), 1);
                attenuationMaxAmpY = p_y(1);
                fit_line_maxY      = exp(polyval(p_y, vecPosX0_col(fit_mask_y)));
            else
                attenuationMaxAmpY = NaN;  fit_mask_y = false(size(vecPosX0_col));  fit_line_maxY = [];
            end
        else
            attenuationMaxAmpY = NaN;  fit_mask_y = false(size(vecPosX0_col));  fit_line_maxY = [];
        end

        title_str_x = sprintf(['$N=%d,\\ P=%g,\\ \\omega_D=%g,\\ K=%g,\\ B_v=%g,\\ \\mathrm{seed}=%d$\n' ...
            '$\\alpha_x^{\\max}=%.4f,\\ \\alpha_x^{\\mathrm{DFT}}=%.4f$'], ...
            N, P, w_D, K, Bv, seed, -attenuationMaxAmpX, -attenuation_x);
        title_str_y = sprintf(['$N=%d,\\ P=%g,\\ \\omega_D=%g,\\ K=%g,\\ B_v=%g,\\ \\mathrm{seed}=%d$\n' ...
            '$\\alpha_y^{\\max}=%.4f,\\ \\alpha_y^{\\mathrm{DFT}}=%.4f$'], ...
            N, P, w_D, K, Bv, seed, -attenuationMaxAmpY, -attenuation_y);

        figure;
        if any(valid_x)
            semilogy(vecPosX0_col(valid_x), maxAmpXAfter2(valid_x), 'o', ...
                'DisplayName', '$|\Delta x|_{\max}$');  hold on;
            if any(fit_mask_x)
                semilogy(vecPosX0_col(fit_mask_x), fit_line_maxX, '-', ...
                    'DisplayName', sprintf('Fit: $\\alpha_x^{\\max}=%.4f$', -attenuationMaxAmpX));
            end
        else; hold on; end
        semilogy(initial_distance_from_oscillation_output_x_fft, amplitude_vector_x, '*', ...
            'DisplayName', '$|\Delta x|$ (DFT)');
        hold off;
        xlabel('$x_0$','Interpreter','latex');  ylabel('Amplitude','Interpreter','latex');
        title(title_str_x,'Interpreter','latex');
        legend('Interpreter','latex','Location','best');
        set(gca,'TickLabelInterpreter','latex');

        figure;
        if any(valid_y)
            semilogy(vecPosX0_col(valid_y), maxAmpYAfter2(valid_y), 'o', ...
                'DisplayName', '$|u_y|_{\max}$');  hold on;
            if any(fit_mask_y)
                semilogy(vecPosX0_col(fit_mask_y), fit_line_maxY, '-', ...
                    'DisplayName', sprintf('Fit: $\\alpha_y^{\\max}=%.4f$', -attenuationMaxAmpY));
            end
        else; hold on; end
        semilogy(initial_distance_from_oscillation_output_y_fft, amplitude_vector_y, '*', ...
            'DisplayName', '$|u_y|$ (DFT)');
        hold off;
        xlabel('$x_0$','Interpreter','latex');  ylabel('Amplitude','Interpreter','latex');
        title(title_str_y,'Interpreter','latex');
        legend('Interpreter','latex','Location','best');
        set(gca,'TickLabelInterpreter','latex');
    end

%% Dimensionless quantities
    diameter_average                        = mean(vecDiameter);
    attenuation_x_dimensionless             = attenuation_x * diameter_average;
    attenuation_y_dimensionless             = attenuation_y * diameter_average;
    wavenumber_x_dimensionless              = wavenumber_x  * diameter_average;
    wavenumber_y_dimensionless              = wavenumber_y  * diameter_average;
    driving_angular_frequency_dimensionless = w_D * sqrt(mass_particle_average/K);
    gamma_dimensionless                     = Bv / sqrt(K * mass_particle_average);
    pressure_dimensionless                  = P;

%% Save
    fprintf('[simMD] Saving output to: %s\n', save_path);

    save(save_path, ...
        'gamma_dimensionless', 'index_particles', ...
        'attenuation_x_dimensionless', 'attenuation_y_dimensionless', ...
        'wavenumber_x_dimensionless',  'wavenumber_y_dimensionless', ...
        'wavespeed_x', 'wavespeed_y', ...
        'driving_angular_frequency_dimensionless', ...
        'attenuation_fit_line_x', ...
        'initial_distance_from_oscillation_output_x_fft', ...
        'initial_distance_from_oscillation_output_y_fft', ...
        'amplitude_vector_x', 'amplitude_vector_y', ...
        'pressure_dimensionless', 'seed', 'input_pressure', ...
        'unwrapped_phase_vector_x', 'unwrapped_phase_vector_y', ...
        'x_fft_initial_y', 'x_fft_initial_z', ...
        'y_fft_initial_y', 'y_fft_initial_z');

    if options.fullSpectrum
        save(save_path, 'freq_axis', 'spec_x', 'spec_y', '-append');
        fprintf('[simMD] Full spectrum saved.\n');
    end

    if options.maxAmpTracking
        attenuationMaxAmpX_dimensionless = attenuationMaxAmpX * diameter_average;
        attenuationMaxAmpY_dimensionless = attenuationMaxAmpY * diameter_average;
        save(save_path, 'vecPosX0', 'maxAmpXAfter2', 'maxAmpYAfter2', ...
            'attenuationMaxAmpX_dimensionless', 'attenuationMaxAmpY_dimensionless', '-append');
    end
    fprintf('[simMD] Output saved.\n');

    try
        memTest('memlog.txt');
    catch
        % ignore if memTest unavailable
    end

    fprintf('[simMD] DONE (total %.1f s)\n', toc(simMD_tStart));

%% Error handler
catch simMD_ME
    fprintf(2, '[simMD] ERROR: %s\n', simMD_ME.message);
    for simMD_k = 1:length(simMD_ME.stack)
        fprintf(2, '[simMD]   in %s (line %d)\n', ...
            simMD_ME.stack(simMD_k).name, simMD_ME.stack(simMD_k).line);
    end
    rethrow(simMD_ME);
end

end  % function






