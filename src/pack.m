function pack(N, K, D, G, M, P_target, seed, plotit, x_mult, y_mult, z_mult, calc_eig, save_path, options)


    arguments
        N        (1,1) double {mustBeInteger, mustBePositive} = 100
        K        (1,1) double {mustBePositive} = 100
        D        (1,1) double {mustBePositive} = 1
        G        (1,1) double {mustBePositive} = 1.4
        M        (1,1) double {mustBePositive} = 1
        P_target (1,1) double {mustBePositive} = 0.0001
        seed     (1,1) double {mustBeInteger, mustBePositive} = 1
        plotit   (1,1) logical = false
        x_mult   (1,1) double = 1
        y_mult   (1,1) double = 1
        z_mult   (1,1) double = 0
        calc_eig (1,1) logical = false
        save_path (1,1) string = "./junkyard"
        % options is passed as a plain positional struct (the frictional test
        % calls pack(...,save_path,opts)). A fully-defaulted struct means the
        % frictionless 13-arg calls get all fields, and a caller may pass a
        % partial struct; missing fields are backfilled below.
        options (1,1) struct = struct('hertzian', false, ...
            'flagFrictionOn', false, ...
            'scalFricCoef', 0.50, ...
            'scalTangentialK', 1/3, ...
            'scalGammaNormal', 0, ...
            'scalGammaTangential', 0, ...
            'saveFrictionalState', false, ...
            'saveFullState', false)
    end

     % Backfill any option fields a caller omitted so both the frictionless
     % (13-arg) and frictional (14-arg with partial opts) call styles work.
    if ~isfield(options, 'hertzian')
        options.hertzian = false;
    end
    if ~isfield(options, 'flagFrictionOn')
        options.flagFrictionOn = false;
    end
    if ~isfield(options, 'scalFricCoef')
        options.scalFricCoef = 0.50;
    end
    if ~isfield(options, 'scalTangentialK')
        options.scalTangentialK = 1/3;
    end
    if ~isfield(options, 'scalGammaNormal')
        options.scalGammaNormal = 0;
    end
    if ~isfield(options, 'scalGammaTangential')
        options.scalGammaTangential = 0;
    end
    if ~isfield(options, 'saveFrictionalState')
        options.saveFrictionalState = false;
    end
    if ~isfield(options, 'saveFullState')
        options.saveFullState = false;
    end

    % check to see if 3d path is needed
    boolThreeD = (z_mult ~= 0);

    %% Guard: Cundall-Strack friction is implemented for 2D packings only.
    %% 3D friction (rotation about 3 axes) requires a different model.
    if boolThreeD && options.flagFrictionOn
        error('pack:Friction3DNotSupported', 'flagFrictionOn = true is 2D only.');
    end


    % Check if packing already exists — skip if so
    if boolThreeD
        scalRoundedWidth = round(N^(1/3));
        if options.hertzian
            strFilename = sprintf('%s3D_N%d_P%s_Width%d_Seed%d_Hertz.mat', ...
                save_path, N, num2str(P_target), scalRoundedWidth, seed);
        else
            strFilename = sprintf('%s3D_N%d_P%s_Width%d_Seed%d.mat', ...
                save_path, N, num2str(P_target), scalRoundedWidth, seed);
        end
        % Full-state (pre-cleanRats) filename for the repeat-tile source file.
        % A '_Full' tag keeps it distinct from the backbone .mat so re-running
        % pack on the same parameters does not clobber the tile source.
        if options.hertzian
            strFullFilename = sprintf('%s3D_N%d_P%s_Width%d_Seed%d_Full_Hertz.mat', ...
                save_path, N, num2str(P_target), scalRoundedWidth, seed);
        else
            strFullFilename = sprintf('%s3D_N%d_P%s_Width%d_Seed%d_Full.mat', ...
                save_path, N, num2str(P_target), scalRoundedWidth, seed);
        end
    else
        scalRoundedWidth = round(sqrt(N));
        if options.hertzian
            strFilename = sprintf('%s2D_N%d_P%s_Width%d_Seed%d_Hertz.mat', ...
                save_path, N, num2str(P_target), scalRoundedWidth, seed);
        else
            strFilename = sprintf('%s2D_N%d_P%s_Width%d_Seed%d.mat', ...
                save_path, N, num2str(P_target), scalRoundedWidth, seed);
        end
    end
    if isfile(strFilename)
        fprintf('Packing already exists, skipping: %s\n', strFilename);
        return;
    end

    tic

    rng(seed)

%% Box and particle setup
    if boolThreeD
        scalBoxWidthX  = 2*N^(1/3)*D;   % Lx
        scalBoxHeightY = 2*N^(1/3)*D;   % Ly
        scalBoxDepthZ  = 2*N^(1/3)*D;   % Lz
    else
        scalBoxWidthX  = 2*sqrt(N)*D;   % Lx
        scalBoxHeightY = 2*sqrt(N)*D;   % Ly
    end

    scalDissipationVelocity = 0.1;  % Bv: velocity-dependent dissipation prefactor
    scalDissipationAbsolute = 0.5;  % B:  absolute (drag) dissipation
    scalTemperature = 1;    % T:  initial velocity scale

    %% Equal number of small and large particles
    scalNumSmall = N/2;

    % Assign diameters: smallest half get D, largest half get D*G
    [~, vecSortIdx] = sort(rand(N, 1)); % [N x 1] randomize the particle indices
    vecDiameter = D * G * ones(N, 1); % [N x 1] default all to large
    vecDiameter(vecSortIdx(1:scalNumSmall)) = D; % overwrite bottom half with small

    if options.hertzian
        vecRadii = vecDiameter / 2; % [N x 1]
    end

    % Precompute contact-distance matrix: matContactDist(i,j) = (r_i + r_j)
    % Avoids recomputing inside the force loop every timestep
    matContactDist = (vecDiameter + vecDiameter') / 2;  % [N x N]

    %% Physical parameters
    scalGravity = 0;
    scalPressure = 0;
    scalPressureFastGrow = P_target / 50;
    scalCompressionRate = P_target;
    %% Frictional (Cundall-Strack) compression protocol — 2D only.
    %%
    %%
    %%  Protocol (literature: OverDamp.cpp's Acc_max < Fthresh, Vinutha
    %% & Sastry DEM relaxation "<|F_tot|> < threshold", Silbert et al. pressure-
    %% controlled compression): drive the box only while P is outside a dead-band
    %% around P_target (compress when loose, expand when dense), HOLD the box
    %% while P is in-band, and accept only when the grains are in FORCE BALANCE —
    %% max_i |F_net,i| < tol * mean_contact_force — sustained for several hundred
    %% consecutive in-band steps, with a percolating contact network
    %% (mean Zn >= scalFrictionZmin). Force balance is the physically correct
    %% "settled" signal: a balanced packing has ~zero accelerations, so KE
    %% decays and P becomes a smooth function of box size (no limit cycle).
    scalFrictionRate         = 0.005;      % box change per step while P is outside the dead-band
    scalFrictionDeadBand     = 0.15;        % P in-band: |P - P_target|/P_target < 15%
    scalFrictionForceTol     = 0.01;        % accept when max|F_net|/mean|F_contact| < 1%
    scalFrictionBalCount      = 0;          % consecutive in-band steps satisfying force balance
    scalFrictionBalWindow     = 300;        % sustained force-balance steps to accept
    scalFrictionZmin          = 2.5;        % percolation guard: mean Zn above this to accept
                                             % (frictional 2D isostatic z_iso = 3)
    scalFrictionMaxSteps      = 3e6;        % hard cap for the frictional phase (safety only)
    scalFrictionVelDecay      = 0.10;       % per-step velocity/omega decay for the frictional
                                             % relaxation. Damping rate = decay/dt ~ 16 per
                                             % time unit >> contact spring frequency sqrt(K/M)=10,
                                             % so the Cundall-Strack springs and rotational DOF
                                             % are OVERDAMPED: energy drains out instead of
                                             % feeding the P limit-cycle. The renormalized drag
                                             % above is a no-op at this dt (Bn*dt ~ 3e-3).
                                             % Gated on boolFrictionOn: the frictionless path
                                             % is untouched.
    scalSlowConvSteps           = 20000;        % for the FRICTIONLESS slow phase:
                                              % consecutive steps the box must be
                                              % frozen (|Lx-Lx_prev| < 1e-8) to call
                                              % it converged. The energy-< 1e-20 gate
                                              % is physically unreachable, so without
                                              % this the 3D phase sits on its flat
                                              % pressure plateau (P/P_target ~1.53)
                                              % and would otherwise run to scalMaxSteps.
    scalLxFrozenSlow             = 0;             % running count for the frictionless
                                               % frozen-box convergence
    scalLxPrevSlow               = 0;             % previous Lx for the frictionless
                                               % frozen-box convergence check
    scalCompressionRateFast = 0.01;

    boolConverged = false; % only update plot after each compression step
    boolCellUpdateNeeded = true; % make sure to update cell list on first step
    boolFastCompressPhase = true;
    scalMeanCoordNum = NaN;  % initiated here so it exisits before the loop for echoing to screen
                             % mean particle-particle coordination number,
                             % tracked on EVERY pathway (frictionless /
                             % frictional, 2D / 3D) for the force-balance
                             % convergence logic; the value SAVED with the
                             % packing is recomputed after cleanRats via
                             % the shared computeMeanCoordNum (backbone
                             % only, matching the stored particles)
             %% Cundall-Strack tangential friction parameters
            %%   boolFrictionOn       = master switch; when false all friction
            %%          paths below are skipped, relaxation identical to original.
            %%   scalMu             = Coulomb coefficient; caps |F_t| <= mu*|F_n|.
            %%   scalKtOverK        = tangential/normal stiffness ratio K_t/K
            %%                        (Cundall-Strack reference: 1/3).
            %%   scalGammaNormal    = normal dashpot prefactor; 0 keeps original.
            %%   scalGammaTangential = tangential dashpot prefactor (optional).
            %%   boolSaveFricState   = export fricState sidecar .mat before cleanRats.
    boolFrictionOn          = options.flagFrictionOn;
    scalMu                  = options.scalFricCoef;
    scalKtOverK             = options.scalTangentialK;
    scalKt                  = K * scalKtOverK;
    scalGammaNormal          = options.scalGammaNormal;
    scalGammaTangential     = options.scalGammaTangential;
    boolSaveFricState       = options.saveFrictionalState;



%% Display / simulation parameters
    boolPlotKE = false;
    scalPlotSkip = 1000;   % timesteps between plot updates
    scalCellUpdateInterval = 1;

    % time step should be 1/100 of a particle oscillation period
    % For Hertzian, k_eff ~ 2K*sqrt(Reff*delta) is unknown at init,
    % so use a conservative estimate based on largest R and largest delta
    % largest Reff ~ G*D/4 (large-large contact)
    % and largest delta ~ 1e-1 (max target pressure)
    if options.hertzian
        scalTimestep = 2*pi * sqrt(M / (2*K*sqrt(G*D/4 .* 1E-1))) * 0.01;
    else
        scalTimestep = 2*pi * sqrt(M/K) * 0.01;
    end
    scalMaxSteps = 5e6; % sane hard cap: every phase now converges on a frozen
                        % box (scalSlowConvSteps / frictional controller) well
                        % before this. The old 1e8 cap let a non-converging
                        % phase run ~24h and exhaust memory (crashed the machine).

%% Initial conditions — place particles on a grid then shuffle
    %  Keep them D/2 from the walls
    % Frist point at D/2 from left wall
    % Last point at (Lx - D/2) from right wall
    % Dividing by G*D gives the number of particles in the row/column
    % The spacing between particles is G*D
    % Result is x matrix of [N x N] positions where the row is the x-position
    % and the column is the y-position
    % meshgrid() is more inuitive (row is y-position,etc),
    % bu ndgrid() is faster
    if boolThreeD
    [vecPosX, vecPosY, vecPosZ] = ndgrid(D/2 : G*D : scalBoxWidthX-D/2, ...
                                          D/2 : G*D : scalBoxHeightY-D/2, ...
                                          D/2 : G*D : scalBoxDepthZ-D/2);
    else
        [vecPosX, vecPosY] = ndgrid(D/2 : G*D : scalBoxWidthX-D/2, ...
                                    D/2 : G*D : scalBoxHeightY-D/2);
    end

    % shuffle the particles to avoid crystallization
    % shuffle the indices of the particles
    [~, vecShuffleIdx] = sort(rand(numel(vecPosX), 1));  % [numel x 1]
    vecPosX = vecPosX(vecShuffleIdx(1:N));   % [N x 1] of particle x-positions
    vecPosY = vecPosY(vecShuffleIdx(1:N));   % [N x 1] etc
    if boolThreeD
        vecPosZ = vecPosZ(vecShuffleIdx(1:N));
    end

    % assign random initial velocities
    vecVelX = sqrt(scalTemperature) * randn(N, 1);  % [N x 1]
    vecVelX = vecVelX - mean(vecVelX);
    vecVelY = sqrt(scalTemperature) * randn(N, 1);  % [N x 1]
    vecVelY = vecVelY - mean(vecVelY);
    if boolThreeD
        vecVelZ = sqrt(scalTemperature) * randn(N, 1);
        vecVelZ = vecVelZ - mean(vecVelZ);  % [N x 1]
    end

    % start with zero accelerations
    vecAccelXPrev = zeros(N, 1);  % [N x 1]
    vecAccelYPrev = zeros(N, 1);  % [N x 1]
    if boolThreeD
        vecAccelZPrev = zeros(N, 1);  % [N x 1]
    end

    %% Rotational DOFs for 2D disks (out-of-plane z rotation)
    %%   vecOmega        = angular velocity omega_i [N x 1]
    %%   vecAlphaPrev    = previous angular acceleration, Verlet half-step.
    %%   Solid disk: I_i = (M_i * r_i^2) / 2,  alpha = torque / I
    if boolFrictionOn
        vecOmega       = zeros(N, 1);
        vecAlphaPrev   = zeros(N, 1);
    end



    vecKineticEnergyHistory   = zeros(scalMaxSteps, 1);  % [scalMaxSteps x 1]
    vecPotentialEnergyHistory = zeros(scalMaxSteps, 1);  % [scalMaxSteps x 1]


             %% Cundall-Strack tangential spring state
            %%   matDispTan(i,j)   = tangential spring displacement for pair (i,j)
            %%                   stored at linear index i + N*(j-1)
            %%   matDispTanStuck    = pair has overlapped at least once
            %%   N x N matrices store the displacement for all N*(N-1)/2 pairs,
            %%   matching Energy_Disk_VL.cpp.
    if boolFrictionOn
        matDispTan       = zeros(N, N);
        matDispTanStuck  = false(N, N);
    end

%% Verlet cell list setup
    % Determine cell size rounded to be at least 1*G*D
    % to avoid missing interactions, I tried this out 
    % many times and 3 works best
    scalRawCellWidth = 3 * G * D; % Changing this will mess up findNeighbors2D. Update findNeighbors2D to be like 3D version (so it doesn't double count neighbors) before changing this.

    % Divide the box into integer number of cells
    % so that the cell width is a multiple of scalRawCellWidth
    scalNumCellsX  = round(scalBoxWidthX  / scalRawCellWidth);
    scalCellWidthX = scalBoxWidthX  / scalNumCellsX;
    scalNumCellsY  = round(scalBoxHeightY / scalRawCellWidth);
    scalCellWidthY = scalBoxHeightY / scalNumCellsY;
    if boolThreeD
        scalNumCellsZ  = round(scalBoxDepthZ  / scalRawCellWidth);
        scalCellWidthZ = scalBoxDepthZ  / scalNumCellsZ;
    end

    % --- Warn about coarse cell grids (especially 2D double-count edge case) ---
    if ~boolThreeD
        % 2D: old findNeighbors2D logic will double-count neighbors when an axis has 2 cells
        if scalNumCellsX == 2 || scalNumCellsY == 2
            warning('GranE:CellGridTwoCells2D', ...
                ['Cell grid has only 2 cells along at least one axis: ', ...
                 'scalNumCellsX = %d, scalNumCellsY = %d. ', ...
                 'Current findNeighbors2D will double-count neighbor pairs in this regime. ', ...
                 'Consider updating findNeighbors2D to the 3D-style stencil before using this setup.'], ...
                 scalNumCellsX, scalNumCellsY);
        end
    else
        % 3D: this is now handled correctly by findNeighbors3D + unique(...,'rows'),
        % but you might still want a "coarse grid" warning if any axis has only 1-2 cells.
        if scalNumCellsX <= 2 || scalNumCellsY <= 2 || scalNumCellsZ <= 2
            warning('GranE:CellGridCoarse3D', ...
                ['Cell grid is very coarse in 3D: scalNumCells = [%d %d %d]. ', ...
                 'findNeighbors3D will still be correct, but performance/contact statistics ', ...
                 'may be affected. Consider using more cells per axis if feasible.'], ...
                 scalNumCellsX, scalNumCellsY, scalNumCellsZ);
        end
    end
    % Wrap positions into [0, L) before rebuilding cell list
    % mod(x, L) == x - L*floor(x/L), handles both positive and negative overshoot
    vecPosX = mod(vecPosX, scalBoxWidthX);   % [N x 1]
    vecPosY = mod(vecPosY, scalBoxHeightY);  % [N x 1]
    if boolThreeD
        vecPosZ = mod(vecPosZ, scalBoxDepthZ);
    end

    % Map positions to cell indices in [1, scalNumCells]
    % ceil gives 1-based index; clamp handles the x=0 edge case where ceil returns 0
    % vecCellIdxX = min(max(ceil(vecPosX / scalCellWidthX), 1), scalNumCellsX);  % [N x 1]
    % vecCellIdxY = min(max(ceil(vecPosY / scalCellWidthY), 1), scalNumCellsY);  % [N x 1]

    % Convert 2D cell index to linear index, then bucket particles in one pass
    % Grid (3x3 example):
    %
    % (1,3) (2,3) (3,3)      7  8  9
    % (1,2) (2,2) (3,2)  ->  4  5  6
    % (1,1) (2,1) (3,1)      1  2  3
    %
    % formula: idx = x + numCellsX * (y - 1)
    % e.g. cell (2,3): 2 + 3*(3-1) = 2 + 6 = 8
    % vecCellLinearIdx = vecCellIdxX + scalNumCellsX * (vecCellIdxY - 1);  % [N x 1]

    % Particle:    1    2    3    4    5
    % Cell index:  3    1    3    2    1
    %
    % accumarray groups them:
    % cell 1 -> [2, 5]
    % cell 2 -> [4]
    % cell 3 -> [1, 3]
    % accumarray groups particle indices by cell — O(N) instead of O(N * numCells)
    % cellParticleList = accumarray(vecCellLinearIdx, (1:N)', [scalNumCellsX*scalNumCellsY 1], @(x){x});

    % reshape just converts that flat list back into a 2D grid so you can look up neighbors naturally by (ix, iy) index.
    % cellParticleList = reshape(cellParticleList, scalNumCellsX, scalNumCellsY);  % [scalNumCellsX x scalNumCellsY]

    %% Setup plotting
    if plotit
        figure(1), clf;
        hPlotHandles = gobjects(N, 1);  % [N x 1]
        for np = 1:N
            hPlotHandles(np) = rectangle( ...
                'Position',  [vecPosX(np) - 0.5*vecDiameter(np), ...
                               vecPosY(np) - 0.5*vecDiameter(np), ...
                               vecDiameter(np), vecDiameter(np)], ...
                'Curvature', [1 1], 'EdgeColor', 'b');
        end
        axis equal; axis([0 scalBoxWidthX 0 scalBoxHeightY]);
        figure(2), clf;
    end

    %% Main time-integration loop
    scalLastCompressStep = 0;

    % Pre-allocate pair buffers — N*12 is safe upper bound for 2D jamming
    % because for 2D jamming, the maximum number of contacts is 6N.
    % so twice that is a safe upper bound for 2D jamming
    % we're goingt to store pairs like this:
    % vecPairIdxSource = [1, 1, 2, 3, ...]   <- first particle of each pair
    % vecPairIdxDest = [2, 3, 3, 4, ...]   <- second particle of each pair
    vecPairIdxSource = zeros(N*12, 1);  % [N*12 x 1]
    vecPairIdxDest = zeros(N*12, 1);  % [N*12 x 1]
    scalMaxPairs = N*12;

    fprintf('Starting main integration loop (max %d steps)...\n', scalMaxSteps);
    scalLogInterval = round(0.05 * scalMaxSteps);  % 5% of max steps
    for nt = 1:scalMaxSteps

        % Progress logging
        if mod(nt, 5000) == 0
            fprintf('  step %d | P=%.4e | P_target=%.4e | P/P_target=%.4f\n', ...
                nt, scalPressure, P_target, scalPressure/P_target);
        end
        boolConverged = false;

        %% Plotting
        if plotit && mod(nt, scalPlotSkip) == 0
            if boolConverged
                figure(1);
                for np = 1:N
                    set(hPlotHandles(np), 'Position', ...
                        [vecPosX(np) - 0.5*vecDiameter(np), ...
                         vecPosY(np) - 0.5*vecDiameter(np), ...
                         vecDiameter(np), vecDiameter(np)]);
                end
                ylim([0 scalBoxHeightY]); xlim([0 scalBoxWidthX]);
                title(num2str(scalBoxHeightY));
            end
            figure(2);
            semilogy(nt, vecKineticEnergyHistory(nt-1),   'ro'); hold on;
            semilogy(nt, vecPotentialEnergyHistory(nt-1), 'bs');
            plot(nt, scalPressure, 'kx');
            drawnow;
        elseif boolPlotKE && mod(nt, scalPlotSkip) == 0
            figure(1), plot(vecPosX, vecPosY, 'k.'); drawnow;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% First step in Verlet integration %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        vecPosX = vecPosX + vecVelX*scalTimestep + vecAccelXPrev.*(scalTimestep^2/2);
        vecPosY = vecPosY + vecVelY*scalTimestep + vecAccelYPrev.*(scalTimestep^2/2);
        if boolThreeD
            vecPosZ = vecPosZ + vecVelZ*scalTimestep + vecAccelZPrev.*(scalTimestep^2/2);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Re-assign particles to cells %%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if boolCellUpdateNeeded || mod(nt, scalCellUpdateInterval) == 0
            if boolThreeD
                [vecPosX, vecPosY, vecPosZ, cellParticleList, scalNumCellsX, scalNumCellsY, scalNumCellsZ] = ...
                    rebuildCellList3D(vecPosX, vecPosY, vecPosZ, scalBoxWidthX, scalBoxHeightY, scalBoxDepthZ, scalRawCellWidth, scalTimestep, N);
            else
                [vecPosX, vecPosY, cellParticleList, scalNumCellsX, scalNumCellsY] = ...
                    rebuildCellList(vecPosX, vecPosY, scalBoxWidthX, scalBoxHeightY, scalRawCellWidth, scalTimestep, N);
            end
            boolCellUpdateNeeded = false;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Build candidate pair list %%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Collect unique pairs (idxNN < idxMM) from neighboring cells.
        % Enforcing idxNN < idxMM means each pair is visited once,
        % halving the number of force evaluations.

        if boolThreeD
            [vecPairIdxSource, vecPairIdxDest, scalNumPairs, scalMaxPairs] = findNeighbors3D( ...
                cellParticleList, scalNumCellsX, scalNumCellsY, scalNumCellsZ, ...
                vecPairIdxSource, vecPairIdxDest, scalMaxPairs);
        else
            [vecPairIdxSource, vecPairIdxDest, scalNumPairs, scalMaxPairs] = findNeighbors2D( ...
                cellParticleList, scalNumCellsX, scalNumCellsY, ...
                vecPairIdxSource, vecPairIdxDest, scalMaxPairs);
        end

        vecActivePairSource = vecPairIdxSource(1:scalNumPairs);  % [scalNumPairs x 1]
        vecActivePairDest = vecPairIdxDest(1:scalNumPairs);  % [scalNumPairs x 1]

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Vectorized force evaluation %%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % All arithmetic operates on [scalNumPairs x 1] vectors — no inner loops.
        % no mod() here because we need signed distances
        vecSepX = vecPosX(vecActivePairDest) - vecPosX(vecActivePairSource);
        vecSepX = vecSepX - scalBoxWidthX  * round(vecSepX / scalBoxWidthX);
        vecSepY = vecPosY(vecActivePairDest) - vecPosY(vecActivePairSource);
        vecSepY = vecSepY - scalBoxHeightY * round(vecSepY / scalBoxHeightY);
        if boolThreeD
            vecSepZ = vecPosZ(vecActivePairDest) - vecPosZ(vecActivePairSource);
            vecSepZ = vecSepZ - scalBoxDepthZ * round(vecSepZ / scalBoxDepthZ);
        end

        vecContactDist = matContactDist(vecActivePairSource + N*(vecActivePairDest-1));

        if boolThreeD
            vecSepDistSq = vecSepX.^2 + vecSepY.^2 + vecSepZ.^2;
        else
            vecSepDistSq = vecSepX.^2 + vecSepY.^2;
        end
        boolContact = vecSepDistSq < vecContactDist.^2;      % [scalNumPairs x 1] overlapping pairs only
        scalNumContacts = sum(boolContact);      % [1 x 1] number of active contact pairs
        scalTangentialPE = 0.0;   % tangential spring energy this step, excluded from pressure

         % Trim all arrays to only pairs in contact
        vecSepX = vecSepX(boolContact); % [scalNumContacts x 1]
        vecSepY = vecSepY(boolContact); % [scalNumContacts x 1]
        if boolThreeD
            vecSepZ = vecSepZ(boolContact);
        end
        vecSepDistSq = vecSepDistSq(boolContact); % [scalNumContacts x 1]
        vecContactDist = vecContactDist(boolContact); % [scalNumContacts x 1]
        vecContactNN = vecActivePairSource(boolContact);% [scalNumContacts x 1]
        vecContactMM = vecActivePairDest(boolContact);% [scalNumContacts x 1]

        if options.hertzian
            vecRadiiNN = vecRadii(vecContactNN); % [scalNumContacts x 1]
            vecRadiiMM = vecRadii(vecContactMM); % [scalNumContacts x 1]
            % this falls out of the math for two parabaloids https://en.wikipedia.org/wiki/Contact_mechanics
            vecRadiiEff = (vecRadiiNN .* vecRadiiMM) ./ (vecRadiiNN + vecRadiiMM); % [scalNumContacts x 1] 
        end

        vecSepDist = sqrt(vecSepDistSq); % [scalNumContacts x 1]
        vecOverlap = vecContactDist - vecSepDist; % [scalNumContacts x 1] positive when overlapping

        % Force magnitude and potential energy per contact
        if options.hertzian
            vecForceMag = -(4/3) .* K .* sqrt(vecRadiiEff) .* vecOverlap.^(3/2);
            vecPotentialContact = (4/3) .* (2/5) * K .* sqrt(vecRadiiEff) .* vecOverlap.^(5/2);
        else
            vecForceMag = -K .* vecOverlap; % [scalNumContacts x 1]
            vecPotentialContact = 0.5 * K .* vecOverlap.^2; % [scalNumContacts x 1]
        end

        % Unit vectors along contact normal
        vecNormalX = vecSepX ./ vecSepDist; % [scalNumContacts x 1]
        vecNormalY = vecSepY ./ vecSepDist; % [scalNumContacts x 1]
        if boolThreeD
            vecNormalZ = vecSepZ ./ vecSepDist;
        end

        % Velocity-dependent dissipation projected onto contact normal
        scalReducedMass = M / 2; % may need to change this for non-uniform mass if the future
        vecRelVelDotNormal = (vecVelX(vecContactNN) - vecVelX(vecContactMM)) .* vecNormalX ...
                           + (vecVelY(vecContactNN) - vecVelY(vecContactMM)) .* vecNormalY;
        if boolThreeD
            vecRelVelDotNormal = vecRelVelDotNormal ...
                               + (vecVelZ(vecContactNN) - vecVelZ(vecContactMM)) .* vecNormalZ;
        end
        vecForceDissipative = scalDissipationVelocity * scalReducedMass .* vecRelVelDotNormal;

        % Net contact force components
        vecForceContactX = (vecForceMag - vecForceDissipative) .* vecNormalX;
        vecForceContactY = (vecForceMag - vecForceDissipative) .* vecNormalY;
        if boolThreeD
            vecForceContactZ = (vecForceMag - vecForceDissipative) .* vecNormalZ;
        end

        % Distribute force via Newton's 3rd law
        vecForceX = accumarray(vecContactNN, vecForceContactX, [N 1]) ...
                  - accumarray(vecContactMM, vecForceContactX, [N 1]);
        vecForceY = accumarray(vecContactNN, vecForceContactY, [N 1]) ...
                  - accumarray(vecContactMM, vecForceContactY, [N 1]);
        if boolThreeD
            vecForceZ = accumarray(vecContactNN, vecForceContactZ, [N 1]) ...
                      - accumarray(vecContactMM, vecForceContactZ, [N 1]);
        end

        
              %% =============================
              %%  Cundall-Strack Tangential Friction
              %%  Ref: CundallStrack_2D/Energy_Disk_VL.cpp
              %%
              %%  Tangential spring coord U_ij evolves:  dU/dt = v_t^total
              %%  Capped by Coulomb:  |F_t| = K_t*|disp|  <=  mu*|F_n|
              %%  Tangential viscous damping (optional):
              %%    F_t^visc = -gamma_t * m_red * v_t^total
              %%  Torque:  tau_i = r_i * F_t
              %%  t_hat = (n_y, -n_x)  normal rotated -90 degrees
              %% =============================
        if boolFrictionOn && ~boolThreeD
            vecTorque  = zeros(N, 1);
        end

        if boolFrictionOn && ~boolThreeD && scalNumContacts > 0

            % Per-contact state lookup from the [N x N] tangential-displacement matrix
            vecLinIdx = vecContactNN + N * (vecContactMM - 1);   % [scalNumContacts x 1] linear index into matDispTan
            vecDispTan = matDispTan(vecLinIdx);                 % [scalNumContacts x 1] per-contact tangential displacement
            vecStuck   = matDispTanStuck(vecLinIdx);            % [scalNumContacts x 1] per-contact "has-contacted" flag
            vecDispTan(~vecStuck) = 0;

            % Unit tangent: normal rotated -90 degrees in 2D,  t_hat = (n_y, -n_x)
            vecUnitTanX =  vecNormalY;   % [scalNumContacts x 1] unit tangent x  ( = n_y)
            vecUnitTanY = -vecNormalX;   % [scalNumContacts x 1] unit tangent y  ( = -n_x)

            % Total tangential slip rate (velocity): (v_i - v_j) . t_hat
            vecVelTan = ...                   % [scalNumContacts x 1] total tangential slip rate
                (vecVelX(vecContactNN) - vecVelX(vecContactMM)) .* vecUnitTanX + ...
                (vecVelY(vecContactNN) - vecVelY(vecContactMM)) .* vecUnitTanY;

            % Rotational slip rate: -(omega_i*r_i + omega_j*r_j)
            vecRi = vecDiameter(vecContactNN) / 2;   % [scalNumContacts x 1] radius of grain i
            vecRj = vecDiameter(vecContactMM) / 2;   % [scalNumContacts x 1] radius of grain j
            vecVelTan = vecVelTan - (...     % subtract rotational contribution
                vecOmega(vecContactNN) .* vecRi + vecOmega(vecContactMM) .* vecRj);

            % Advance tangential displacement:  disp(t+dt) = disp(t) + vel(t)*dt
            vecDispTan = vecDispTan + vecVelTan * scalTimestep;

            % Coulomb cap: |F_t| = K_t*|disp| <= mu*|F_n|  ->  clamp the displacement
            vecFnAbs = abs(vecForceMag);                         % [scalNumContacts x 1] |F_n|
            vecDispTanCap = scalMu .* vecFnAbs ./ scalKt;        % [scalNumContacts x 1] Coulomb cap = mu*|F_n|/K_t
            vecDispTan = sign(vecDispTan) .* min(vecDispTanCap, abs(vecDispTan));
            matDispTanStuck(vecLinIdx) = true;

            % Tangential spring force: F_t = -K_t * disp
            vecFtMag = -scalKt .* vecDispTan;      % [scalNumContacts x 1] tangential force magnitude

            % Tangential contact energy: 0.5 * K_t * disp^2
            vecPotentialContact = vecPotentialContact + 0.5 * scalKt .* (vecDispTan .^ 2);
            scalTangentialPE = sum(0.5 * scalKt .* (vecDispTan .^ 2));   % [1 x 1] tangential PE this step

            % Optional tangential viscous damping (dashpot on the slip rate)
            if scalGammaTangential > 0
                vecFtMag = vecFtMag - (scalGammaTangential * (M / 2)) .* vecVelTan;
            end

            % Distribute tangential force via Newton's 3rd law
            vecFtX = vecFtMag .* vecUnitTanX;   % [scalNumContacts x 1] tangential force x-component
            vecFtY = vecFtMag .* vecUnitTanY;   % [scalNumContacts x 1] tangential force y-component
            vecForceX = vecForceX + ...
                accumarray(vecContactNN, vecFtX, [N 1]) - ...
                accumarray(vecContactMM, vecFtX, [N 1]);
            vecForceY = vecForceY + ...
                accumarray(vecContactNN, vecFtY, [N 1]) - ...
                accumarray(vecContactMM, vecFtY, [N 1]);

            % Contact torques: tau_i = r_i * F_t,  tau_j = r_j * F_t
            vecTorque = accumarray(vecContactNN,  vecRi .* vecFtMag, [N 1]) + ...
                           accumarray(vecContactMM, vecRj .* vecFtMag, [N 1]);

            % Persist updated tangential displacement
            matDispTan(vecLinIdx) = vecDispTan;   % write per-contact displacement back to the [N x N] matrix
        end

        % Cundall-Strack: reset the tangential spring state for candidate
        % pairs that have SEPARATED this step. A stale displacement would
        % otherwise re-engage at full strength on re-contact and inject
        % energy — a driver of the P limit-cycle. Only cell-list candidate
        % pairs can become contacts on the next step, so resetting just those
        % (O(pairs)) is sufficient.
        if boolFrictionOn && ~boolThreeD
            vecSeparating = vecActivePairSource(~boolContact) + N * (vecActivePairDest(~boolContact) - 1);
            if ~isempty(vecSeparating)
                matDispTan(vecSeparating) = 0;
                matDispTanStuck(vecSeparating) = false;
            end
        end

            % ============ Rotational velocity-Verlet half-step ============
            if boolFrictionOn && ~boolThreeD
            vecInertiaC = 0.5 * M .* (vecDiameter / 2) .^ 2;
             % Rotational damping: ref OverDamp.cpp L104-106:
               %   W_n+1 = (T_n - Bt*W_n)/Bt_denorm,  Bt_denorm = 1 + Bt*dt/2
               % This is the over-damped form that makes the rotational DOF
               % unconditionally stable. Without it omega oscillates forever.
            scalBt = scalDissipationAbsolute / sqrt(3);
            if scalGammaTangential > 0
                scalBt = scalGammaTangential;    % user override
            end
            Bt_den = 1 + scalBt * scalTimestep / 2;
            vecTorque = (vecTorque - scalBt .* vecOmega) / Bt_den;
            vecAlpha       = vecTorque ./ vecInertiaC;
            vecOmega       = vecOmega + (vecAlphaPrev + vecAlpha) .* (scalTimestep / 2);
            vecAlphaPrev   = vecAlpha;
            end

% Contact count per particle (coordination number)
        vecCoordNum = accumarray(vecContactNN, 1, [N 1]) ...
                    + accumarray(vecContactMM, 1, [N 1]);

        vecPotentialEnergyHistory(nt) = sum(vecPotentialContact) / N;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Drag, boundaries, energy %%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % OverDamped integrator (ref: OverDamp.cpp L97-L100):
          %   F_n+1 = (F_n - Bn*Vel_n) / (1 + Bn*dt/2)
          %   W_n+1 = (W_n - Bt*W_n) / (1 + Bt*dt/2)
          % For frictional contacts the renormalized form prevents the
          % contact-damped oscillation that causes P/P_target to cycle.
          % The frictionless additive form is identical (Bn_den -> 1).
        if boolFrictionOn
            Bn_den = 1 + scalDissipationAbsolute * scalTimestep / 2;
            vecForceX = (vecForceX - scalDissipationAbsolute .* vecVelX) / Bn_den;
            vecForceY = (vecForceY - scalDissipationAbsolute .* vecVelY) / Bn_den;
            if boolThreeD
                vecForceZ = (vecForceZ - scalDissipationAbsolute .* vecVelZ) / Bn_den;
            end
        else
            vecForceX = vecForceX - scalDissipationAbsolute .* vecVelX;  % [N x 1]
            vecForceY = vecForceY - scalDissipationAbsolute .* vecVelY;  % [N x 1]
            if boolThreeD
                vecForceZ = vecForceZ - scalDissipationAbsolute .* vecVelZ;
            end
        end

        % TODO: get rid of these since this is isotropic
        boolLeftWallContact  = vecPosX < vecDiameter/2;
        boolRightWallContact = vecPosX > scalBoxWidthX - vecDiameter/2;

        vecPosX = vecPosX - scalBoxWidthX  .* floor(vecPosX / scalBoxWidthX);
        vecPosY = vecPosY - scalBoxHeightY .* floor(vecPosY / scalBoxHeightY);
        if boolThreeD
            vecPosZ = vecPosZ - scalBoxDepthZ .* floor(vecPosZ / scalBoxDepthZ);
        end

        if boolThreeD
            vecKineticEnergyHistory(nt) = 0.5 * M * sum(vecVelX.^2 + vecVelY.^2 + vecVelZ.^2) / N;
        else
            vecKineticEnergyHistory(nt) = 0.5 * M * sum(vecVelX.^2 + vecVelY.^2) / N;
        end

        % Rotational kinetic energy (friction) must be included in the
        % convergence check: omega carries energy that the translational KE
        % misses; without it scalEk reads ~0 while omega still oscillates.
        if boolFrictionOn && ~boolThreeD
            vecInertiaC = 0.5 * M .* (vecDiameter / 2) .^ 2;   % solid-disk moment of inertia
            vecKineticEnergyHistory(nt) = vecKineticEnergyHistory(nt) ...
                + 0.5 * sum(vecInertiaC .* vecOmega.^2) / N;
        end

        vecAccelX = vecForceX ./ M;
        vecAccelY = vecForceY ./ M - scalGravity;
        if boolThreeD
            vecAccelZ = vecForceZ ./ M;  % no gravity in Z
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Second step in Verlet integration %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        vecVelX = vecVelX + (vecAccelXPrev + vecAccelX) .* (scalTimestep/2);  % [N x 1]
        vecVelY = vecVelY + (vecAccelYPrev + vecAccelY) .* (scalTimestep/2);
        if boolThreeD
                vecVelZ = vecVelZ + (vecAccelZPrev + vecAccelZ) .* (scalTimestep/2);  % ← correct
        end

        % Over-damped frictional relaxation: drain energy out of the
        % translational and rotational DOF every step so the Cundall-Strack
        % springs settle into force balance instead of oscillating. Damping
        % rate scalFrictionVelDecay/dt >> contact spring frequency, so the
        % system is over-damped. Gated on boolFrictionOn so the frictionless
        % path is untouched.
        if boolFrictionOn
            vecVelX = vecVelX * (1 - scalFrictionVelDecay);
            vecVelY = vecVelY * (1 - scalFrictionVelDecay);
            if boolThreeD
                vecVelZ = vecVelZ * (1 - scalFrictionVelDecay);
            end
            vecOmega = vecOmega * (1 - scalFrictionVelDecay);
        end

        % Zero out rattlers (no contacts)
        boolRattler = (vecCoordNum == 0);          % [N x 1]
        vecVelX(boolRattler) = 0;
        vecVelY(boolRattler) = 0;
        vecAccelX(boolRattler) = 0;
        vecAccelY(boolRattler) = 0;
        if boolThreeD
            vecVelZ(boolRattler) = 0;
            vecAccelZ(boolRattler) = 0;
        end             %% Zero rotational rates for rattlers (no contacts => no torque)
        if boolFrictionOn && ~boolThreeD
            vecOmega(boolRattler)       = 0;
            vecAlphaPrev(boolRattler)   = 0;
        end



        vecAccelXPrev = vecAccelX;  % [N x 1]
        vecAccelYPrev = vecAccelY;  % [N x 1]
        if boolThreeD
            vecAccelZPrev = vecAccelZ;
        end

        scalTotalContacts = sum(vecCoordNum) / 2;
        scalWallContacts = sum(boolLeftWallContact) + sum(boolRightWallContact);
        scalNumRattlers = sum(boolRattler);
        scalExcessContacts = scalTotalContacts + scalWallContacts - 2*(N - scalNumRattlers);

        % Mean particle-particle coordination number (full packing, rattlers
        % included). Tracked on every pathway — frictionless and frictional,
        % 2D and 3D — so the converged value is saved with the packing for
        % comparison against the jamming literature (2D frictionless z_iso ~
        % 4, 2D frictional z_iso = 3, 3D frictionless z_iso ~ 6).
        scalMeanCoordNum = mean(vecCoordNum);

        % Pressure estimate from mean potential energy
        scalEp = vecPotentialEnergyHistory(nt);
        % Under friction, the tangential spring energy is a constraint DOF,
        % not a compressive load: exclude it from the box-control pressure so
        % the compression target P_target is reached on the NORMAL contacts.
        if boolFrictionOn
           scalEp = scalEp - scalTangentialPE / N;
        end
        if options.hertzian
            scalPressure = (scalEp * (5/2) / K)^(2/5); % this has an implied d= 1 in the denominator
        else
            scalPressure = sqrt(2 * scalEp / K); % this has an implied d= 1 in the denominator
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% COMPRESSION DECISIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        scalEk = vecKineticEnergyHistory(nt);

        % ===== ENV-GUARDED DIAGNOSTIC (removable; no effect unless GRAN_DIAG set) =====
        if ~isempty(getenv('GRAN_DIAG')) && mod(nt, 5000) == 0
            scZnm = mean(vecCoordNum);  % mean particle-particle coordination
            fprintf('DIAG step=%d Ek=%.4e P=%.4e P/Ptrgt=%.4f Lx=%.4f Zn=%.2f\n', nt, scalEk, scalPressure, scalPressure/P_target, scalBoxWidthX, scZnm);
        end
        % ===== end diagnostic =====

        % Frictional case: drive the box purely by relative-pressure control
        % in the slow phase from the start. The fast-compress two-stage
        % transition to slow requires scalEk < 1e-10, but a frictional packing
        % carries a permanent rotational/tangential KE floor (OverDamp.cpp
        % converges on force, Acc_max < Fthresh, not on energy), so that gate
        % is unreachable and the slow-phase expansion branch -- the only thing
        % that can relieve an over-compressed box -- never fires. Gating on
        % ~boolFrictionOn leaves the frictionless path byte-for-byte identical.
        if boolFrictionOn
              % Frictional box controller + FORCE-BALANCE convergence.
              % The box is driven toward P_target with a relative-pressure
              % dead-band (compress when loose, expand when dense, HOLD when
              % in-band). Convergence is NOT "box happened to stop moving" or
              % "P in a band" (both are transients of the compression); it is
              % per-particle FORCE BALANCE: max_i |F_net,i| / mean |F_contact|
              % below a small tolerance, sustained for several hundred
              % consecutive in-band steps, with a percolating contact network.
              % This matches OverDamp.cpp's Acc_max < Fthresh and the DEM
              % relaxation stopping criteria in the jamming literature.
              boolBoxMoved = false;
              if abs(scalPressure - P_target) / P_target < scalFrictionDeadBand
                  % in-band: hold the box fixed; let the grains relax to balance
                  boolBoxMoved = false;
              elseif scalPressure < P_target * (1 - scalFrictionDeadBand)
                  rate = -scalFrictionRate;   % compress (box too loose)
                  scalBoxWidthX = scalBoxWidthX  *(1 + rate);
                  scalBoxHeightY= scalBoxHeightY*(1 + rate);
                  vecPosX = vecPosX *(1 + rate);
                  vecPosY = vecPosY *(1 + rate);
                  if boolThreeD
                     scalBoxDepthZ = scalBoxDepthZ*(1 + rate);
                     vecPosZ = vecPosZ *(1 + rate);
                  end
                  boolCellUpdateNeeded = true;
                  boolBoxMoved = true;
              else
                  rate =  scalFrictionRate;   % expand (box too dense)
                  scalBoxWidthX = scalBoxWidthX  *(1 + rate);
                  scalBoxHeightY= scalBoxHeightY*(1 + rate);
                  vecPosX = vecPosX *(1 + rate);
                  vecPosY = vecPosY *(1 + rate);
                  if boolThreeD
                     scalBoxDepthZ = scalBoxDepthZ*(1 + rate);
                     vecPosZ = vecPosZ *(1 + rate);
                  end
                  boolCellUpdateNeeded = true;
                  boolBoxMoved = true;
              end

              % Force-balance measure: net contact force (normal + tangential)
              % on each grain. A jammed, settled packing has ~zero net force on
              % every grain (mechanical equilibrium under PBC). Normalized by
              % the mean contact force so the threshold is scale-free.
              if scalNumContacts > 0
                  vecNetFx = accumarray(vecContactNN, vecForceContactX, [N 1]) ...
                           - accumarray(vecContactMM, vecForceContactX, [N 1]) ...
                           + accumarray(vecContactNN, vecFtX, [N 1]) ...
                           - accumarray(vecContactMM, vecFtX, [N 1]);
                  vecNetFy = accumarray(vecContactNN, vecForceContactY, [N 1]) ...
                           - accumarray(vecContactMM, vecForceContactY, [N 1]) ...
                           + accumarray(vecContactNN, vecFtY, [N 1]) ...
                           - accumarray(vecContactMM, vecFtY, [N 1]);
                  vecFnetMag     = sqrt(vecNetFx.^2 + vecNetFy.^2);
                  scalMaxFnet    = max(vecFnetMag);
                  scalMeanFc     = mean([abs(vecForceMag); abs(vecFtMag)]);
                  scalForceRatio = scalMaxFnet / max(scalMeanFc, eps);
              else
                  scalForceRatio = inf;   % no contacts: not balanced, not percolating
              end

              % Acceptance: P in-band, contact network percolates, box held, and
              % sustained force balance. The percolation guard (mean Zn) plus the
              % "box held" guard prevent accepting a loose, unjammed state whose
              % net force is trivially small simply because it carries no load.
              boolInBand   = abs(scalPressure - P_target) / P_target < scalFrictionDeadBand;
              boolPercol   = (scalMeanCoordNum >= scalFrictionZmin);
              boolBalanced = (scalForceRatio < scalFrictionForceTol);
              if boolInBand && boolPercol && boolBalanced && ~boolBoxMoved
                  scalFrictionBalCount = scalFrictionBalCount + 1;
              else
                  scalFrictionBalCount = 0;
              end

              if mod(nt, 5000) == 0
                  fprintf('  [fric] step %d | P/Pt=%.3f | Lx=%.4f | maxFnet/meanFc=%.3e | meanCoordNum=%.2f | balCount=%d\n', ...
                      nt, scalPressure / P_target, scalBoxWidthX, scalForceRatio, scalMeanCoordNum, scalFrictionBalCount);
              end
              if scalFrictionBalCount >= scalFrictionBalWindow
                  fprintf('Frictional convergence (FORCE BALANCE) at step %d | P=%.4e (P/Pt=%.3f) Lx=%.4f maxFnet/meanFc=%.3e meanCoordNum=%.2f\n', ...
                      nt, scalPressure, scalPressure / P_target, scalBoxWidthX, scalForceRatio, scalMeanCoordNum);
                  break;
              elseif nt >= scalFrictionMaxSteps
                  fprintf('Frictional MAX-STEP cap at step %d | P=%.4e (P/Pt=%.3f) maxFnet/meanFc=%.3e meanCoordNum=%.2f\n', ...
                      nt, scalPressure, scalPressure / P_target, scalForceRatio, scalMeanCoordNum);
                  break;
              end
        elseif boolFastCompressPhase
            % ===== Fast-compress phase (frictionless two-stage) =====
            if scalPressure < P_target/50
                scalBoxWidthX= scalBoxWidthX * (1-scalCompressionRateFast);
                scalBoxHeightY = scalBoxHeightY * (1-scalCompressionRateFast);
                vecPosX = vecPosX * (1-scalCompressionRateFast);
                vecPosY = vecPosY * (1-scalCompressionRateFast);
                if boolThreeD
                    scalBoxDepthZ = scalBoxDepthZ * (1-scalCompressionRateFast);
                    vecPosZ = vecPosZ * (1-scalCompressionRateFast);
                end
                boolCellUpdateNeeded = true;
                scalLastCompressStep = nt;
            elseif scalPressure < P_target && scalEk < 1e-8
                scalBoxWidthX   = scalBoxWidthX * (1-scalCompressionRateFast);
                scalBoxHeightY = scalBoxHeightY * (1-scalCompressionRateFast);
                vecPosX = vecPosX * (1-scalCompressionRateFast);
                vecPosY = vecPosY * (1-scalCompressionRateFast);
                if boolThreeD
                    scalBoxDepthZ = scalBoxDepthZ * (1-scalCompressionRateFast);
                    vecPosZ = vecPosZ * (1-scalCompressionRateFast);
                end
                boolCellUpdateNeeded = true;
                scalLastCompressStep = nt;
            elseif scalPressure > P_target && scalEk < 1e-10 && nt > (scalLastCompressStep+100)
                scalBoxWidthX   = scalBoxWidthX * (1+scalCompressionRateFast);
                scalBoxHeightY = scalBoxHeightY * (1+scalCompressionRateFast);
                vecPosX = vecPosX * (1+scalCompressionRateFast);
                vecPosY = vecPosY * (1+scalCompressionRateFast);
                if boolThreeD
                    scalBoxDepthZ = scalBoxDepthZ * (1+scalCompressionRateFast);
                    vecPosZ = vecPosZ * (1+scalCompressionRateFast);
                end
                boolCellUpdateNeeded = true;
                scalLastCompressStep = nt;
                boolFastCompressPhase = false;
            end
        else
             % ===== Slow-phase compression / expansion (frictionless) =====
             if scalPressure < scalPressureFastGrow
                scalBoxWidthX    = scalBoxWidthX    * (1 - scalCompressionRate);
                scalBoxHeightY   = scalBoxHeightY    * (1 - scalCompressionRate);
                vecPosX = vecPosX * (1 - scalCompressionRate);
                vecPosY = vecPosY * (1 - scalCompressionRate);
                if boolThreeD
                    scalBoxDepthZ = scalBoxDepthZ * (1 - scalCompressionRate);
                    vecPosZ = vecPosZ * (1 - scalCompressionRate);
                end
                boolCellUpdateNeeded = true;
                scalLastCompressStep = nt;
             elseif scalPressure < P_target && scalEk < 1e-8
                scalBoxWidthX    = scalBoxWidthX     * (1 - scalCompressionRate);
                scalBoxHeightY    = scalBoxHeightY    * (1 - scalCompressionRate);
                vecPosX = vecPosX * (1 - scalCompressionRate);
                vecPosY = vecPosY * (1 - scalCompressionRate);
                if boolThreeD
                    scalBoxDepthZ = scalBoxDepthZ * (1 - scalCompressionRate);
                    vecPosZ = vecPosZ * (1 - scalCompressionRate);
                end
                boolCellUpdateNeeded = true;
                scalLastCompressStep = nt;
             elseif scalPressure > P_target
                 % Frictionless slow phase: a flat pressure plateau is the sign
                 % of a settled packing (energy < 1e-20 is unreachable). Track a
                 % frozen box — |Lx - Lx_prev| staying near zero for
                 % scalSlowConvSteps is the physically correct convergence.
                 if abs(scalBoxWidthX - scalLxPrevSlow) < 1e-8
                     scalLxFrozenSlow = scalLxFrozenSlow + 1;
                 else
                     scalLxFrozenSlow = 0;
                 end
                 scalLxPrevSlow = scalBoxWidthX;
                 if scalLxFrozenSlow >= scalSlowConvSteps
                     fprintf('Frictionless converged at step %d | P=%.4e P/P_target=%.3f Lx=%.4f\n', ...
                        nt, scalPressure, scalPressure / P_target, scalBoxWidthX);
                     break;
                 end
             end
         end
    end
               %% Build and save frictional state on ORIGINAL indices
            %%   (BEFORE cleanRats renumbers)
            %%  Downstream frictional linear-response pipeline uses
            %%  fricState to assemble the full frictional Hessian.
    if boolFrictionOn && boolSaveFricState
        fricState = buildFricState( ...
            vecPosX, vecPosY, vecDiameter, ...
            scalBoxWidthX, scalBoxHeightY, ...
            matDispTan, K, scalKt, scalMu, ...
            scalGammaNormal, scalGammaTangential, M);
        strFricFilename = [strFilename(1:end-4) '_FricState.mat'];
        save(strFricFilename, 'fricState');
        fprintf('fricState saved to: %s\n', strFricFilename);
    end

    fprintf('Loop finished at step %d.\n', nt);
    N_original = numel(vecPosX);  % particle count before cleanRats (for plot titles)

    %% Snapshot the FULL jammed state (all N particles, rattlers included)
    %    BEFORE cleanRats, for tile-based repetition (see the tiling block at
    %    the end of this file). The per-step box rescalings already wrap
    %    positions into [0, L), so every coordinate is in-box; the 3D tiler
    %    reproduces this state exactly across tile boundaries. Gated on
    %    options.saveFullState so default callers keep byte-identical output.
    if options.saveFullState
        if boolThreeD
            save(strFullFilename, 'vecPosX', 'vecPosY', 'vecPosZ', 'vecDiameter', ...
                'scalBoxWidthX', 'scalBoxHeightY', 'scalBoxDepthZ', ...
                'K', 'P_target', 'scalPressure', 'N_original', 'seed', ...
                'scalRoundedWidth');
        else
            save(strFullFilename, 'vecPosX', 'vecPosY', 'vecDiameter', ...
                'scalBoxWidthX', 'scalBoxHeightY', ...
                'K', 'P_target', 'scalPressure', 'N_original', 'seed', ...
                'scalRoundedWidth');
        end
        % Keep an in-memory copy for the 3D repeat-tile block at the end of
        % this file: after cleanRats the local position/diameter variables are
        % rebound to the BACKBONE, so the full (rattler-inclusive) state must
        % be captured here, before cleanRats runs.
        vecTileSrcX = vecPosX;
        vecTileSrcY = vecPosY;
        vecTileSrcD = vecDiameter;
        NTileSrc    = N;   % full per-tile particle count (rattlers included)
        if boolThreeD
            vecTileSrcZ = vecPosZ;
        end
        fprintf('Full-state tile saved to: %s\n', strFullFilename);
    end

    %% Remove rattlers before saving
    % Shared metric functions (computePackingFraction / computeMeanCoordNum)
    % are the single source of truth for these definitions — pack.m and
    % packRepeatTile.m both call them, so updating one definition updates
    % both algorithms.
    % Packing fraction of the FULL jammed state (all N particles, before
    % rattler removal). This is the number comparable to the literature
    % (e.g. Silbert 2010, 2D bidisperse: 0.843 frictionless -> 0.767 at
    % mu=10); the post-cleanRats scalPackingFraction is the eigen-analysis
    % backbone fraction and is NOT a jamming-state property.
    if boolThreeD
        scalPackingFractionFull = computePackingFraction(vecDiameter, scalBoxWidthX, scalBoxHeightY, scalBoxDepthZ);
    else
        scalPackingFractionFull = computePackingFraction(vecDiameter, scalBoxWidthX, scalBoxHeightY);
    end   % [1 x 1] PF before cleanRats
    fprintf('Packing fraction before cleanRats: %.4f\n', scalPackingFractionFull);

    % The SAVED coordination number is recomputed AFTER cleanRats (see
    % below) so it describes exactly the particles stored in the file —
    % the same convention packRepeatTile.m uses when it recomputes on
    % the particles it loads. In-loop tracking of scalMeanCoordNum
    % (full state, rattlers included) is kept: the force-balance
    % convergence logic reads it every step.

    fprintf('Running cleanRats...\n');
    if boolThreeD
        fprintf('Box dims: Lx=%.4f, Ly=%.4f, Lz=%.4f\n', scalBoxWidthX, scalBoxHeightY, scalBoxDepthZ);
    else
        fprintf('Box dims: Lx=%.4f, Ly=%.4f\n', scalBoxWidthX, scalBoxHeightY);
    end
    fprintf('Radii range: min=%.4f, max=%.4f\n', min(vecDiameter/2), max(vecDiameter/2));

    %% Plot packing BEFORE cleanRats
    % Full jammed state on the ORIGINAL particle indices (all N particles,
    % rattlers included). Same visual style as the post-cleanRats plot
    % below: blue main particles + red periodic ghost tiles.
    figure;
    hold on;
    if boolThreeD
        [sx, sy, sz] = sphere(16);

        % Main particles (full pre-cleanRats packing)
        for np = 1:N
            r = vecDiameter(np)/2;
            surf(r*sx + vecPosX(np), r*sy + vecPosY(np), r*sz + vecPosZ(np), ...
                'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.6);
        end

        % Ghost particles on +x, +y, +z faces
        vecOffsets = [scalBoxWidthX, 0, 0; ...
                      0, scalBoxHeightY, 0; ...
                      0, 0, scalBoxDepthZ];  % [3 x 3] one offset per face

        for iface = 1:3
            ox = vecOffsets(iface, 1);
            oy = vecOffsets(iface, 2);
            oz = vecOffsets(iface, 3);
            for np = 1:N
                r = vecDiameter(np)/2;
                surf(r*sx + vecPosX(np) + ox, ...
                     r*sy + vecPosY(np) + oy, ...
                     r*sz + vecPosZ(np) + oz, ...
                    'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.15);
            end
        end

        axis equal;
        axis([-scalBoxWidthX*0.0 2*scalBoxWidthX ...
              -scalBoxHeightY*0.0 2*scalBoxHeightY ...
              -scalBoxDepthZ*0.0  2*scalBoxDepthZ]);
        axis manual;
        xlabel('x'); ylabel('y'); zlabel('z');
        lighting gouraud;
        camlight;
        rotate3d on;
        title(sprintf('Before cleanRats: N=%d, phi=%.4f', N, scalPackingFractionFull));

    else
        % Main particles (full pre-cleanRats packing)
        for np = 1:N
            rectangle('Position', [vecPosX(np) - vecDiameter(np)/2, ...
                                    vecPosY(np) - vecDiameter(np)/2, ...
                                    vecDiameter(np), vecDiameter(np)], ...
                'Curvature', [1 1], 'FaceColor', 'b', 'EdgeColor', 'none');
        end

        % Ghost particles: all 8 surrounding tiles
        vecOffsets2D = [scalBoxWidthX,  0; ...
                       -scalBoxWidthX,  0; ...
                        0,  scalBoxHeightY; ...
                        0, -scalBoxHeightY; ...
                        scalBoxWidthX,  scalBoxHeightY; ...
                       -scalBoxWidthX,  scalBoxHeightY; ...
                        scalBoxWidthX, -scalBoxHeightY; ...
                       -scalBoxWidthX, -scalBoxHeightY];

        for iface = 1:8
            ox = vecOffsets2D(iface, 1);
            oy = vecOffsets2D(iface, 2);
            for np = 1:N
                rectangle('Position', [vecPosX(np) + ox - vecDiameter(np)/2, ...
                                        vecPosY(np) + oy - vecDiameter(np)/2, ...
                                        vecDiameter(np), vecDiameter(np)], ...
                    'Curvature', [1 1], 'FaceColor', 'r', 'EdgeColor', 'none', ...
                    'FaceAlpha', 0.15);
            end
        end

        axis equal;
        axis([-scalBoxWidthX 2*scalBoxWidthX -scalBoxHeightY 2*scalBoxHeightY]);
        title(sprintf('Before cleanRats: N=%d, phi=%.4f', N, scalPackingFractionFull));
    end
    drawnow;
    hold off;

    % Export the before-cleanRats packing photo. strFilename already includes
    % save_path, so just swap the extension; the export mirrors the .mat's
    % own path resolution, and a failed export never aborts the run.
    % Frictional runs get a '_Fric' tag so their photos don't clobber the
    % frictionless ones (the .mat names for the two 2D pathways collide).
    try
        strPngTag = '';
        if boolFrictionOn
            strPngTag = '_Fric';
        end
        strBeforePng = [strFilename(1:end-4) strPngTag '_BeforeCleanRats.png'];
        print(gcf, strBeforePng, '-dpng', '-r120');
        fprintf('Packing photo saved to: %s\n', strBeforePng);
    catch ME
        warning('pack:PNGExportFailed', 'Could not export before-cleanRats plot: %s', ME.message);
    end

    % TODO: need to decide logic if I want to use PBC or not
    % for now, I'll assume fully periodic since weh're mostly done with DEM
    if boolThreeD
        matPositions = [vecPosX, vecPosY, vecPosZ];
        vecRadii = vecDiameter ./ 2;
        [matPositions, vecRadii] = cleanRats(matPositions, vecRadii, scalBoxHeightY, scalBoxWidthX, scalBoxDepthZ, true);
        vecPosX = matPositions(:,1);
        vecPosY = matPositions(:,2);
        vecPosZ = matPositions(:,3);
    else
        matPositions = [vecPosX, vecPosY];
        vecRadii = vecDiameter ./ 2;
         mu_cleanRats = scalMu * boolFrictionOn;
         [matPositions, vecRadii] = cleanRats(matPositions, vecRadii, scalBoxHeightY, scalBoxWidthX, [], false, mu_cleanRats);
        vecPosX = matPositions(:,1);
        vecPosY = matPositions(:,2);
    end
    vecDiameter = vecRadii .* 2;
    N_clean = size(matPositions, 1);
    fprintf('cleanRats complete. %d particles remaining (of %d original).\n', N_clean, N);
    if N_clean == 0
        warning('All particles removed by cleanRats — packing did not jam. Skipping save.');
        return;
    end

    % Coordination number of the SAVED (backbone) packing, via the shared
    % function — the same convention packRepeatTile.m uses on the particles
    % it loads, so the two files now agree for the same packing.
    if boolThreeD
        scalMeanCoordNum = computeMeanCoordNum(vecPosX, vecPosY, vecDiameter, scalBoxWidthX, scalBoxHeightY, vecPosZ, scalBoxDepthZ);
    else
        scalMeanCoordNum = computeMeanCoordNum(vecPosX, vecPosY, vecDiameter, scalBoxWidthX, scalBoxHeightY);
    end
    fprintf('Mean coordination number after cleanRats: %.4f\n', scalMeanCoordNum);

%% Compute linearized Hertzian contact stiffnesses at jammed state
    if options.hertzian

        % update data because cleanRats may have removed particles
        vecRadii_final  = vecDiameter ./ 2;
        matContactDist_final = (vecDiameter + vecDiameter') / 2;
        scalNumFinal = numel(vecPosX);

        vecHertzNN = zeros(scalNumFinal * 12, 1);
        vecHertzMM = zeros(scalNumFinal * 12, 1);
        vecHertzKeff  = zeros(scalNumFinal * 12, 1);
        scalNumHertzContacts = 0;

        for ii = 1:scalNumFinal
            for jj = ii+1:scalNumFinal
                dx = vecPosX(jj) - vecPosX(ii);
                dy = vecPosY(jj) - vecPosY(ii);
                dx = dx - scalBoxWidthX * round(dx / scalBoxWidthX);
                dy = dy - scalBoxHeightY * round(dy / scalBoxHeightY);
                if boolThreeD
                    dz = vecPosZ(jj) - vecPosZ(ii);
                    dz = dz - scalBoxDepthZ * round(dz / scalBoxDepthZ);
                    scalDist = sqrt(dx^2 + dy^2 + dz^2);
                else
                    scalDist = sqrt(dx^2 + dy^2);
                end

                scalSumRadii = matContactDist_final(ii, jj); % grab the minimum distanced needed for contact
                scalDelta = scalSumRadii - scalDist; % if negative, no contact

                % Go through each contact and assign k = dF/d(delta) for F_hertzian
                if scalDelta > 0
                    scalReff = (vecRadii_final(ii) * vecRadii_final(jj)) / scalSumRadii; % effective radius from curvature
                    scalKeff = 2 * K * sqrt(scalReff * scalDelta); % k_eff = dF/d(delta)
                    scalNumHertzContacts = scalNumHertzContacts + 1; % this is indexes, so +1 because matlab isbase 1 

                    % assign row (scalNumHertzContacts) a particle, the particle it's in contact with, and an k_eff
                    vecHertzNN(scalNumHertzContacts) = ii; 
                    vecHertzMM(scalNumHertzContacts) = jj;
                    vecHertzKeff(scalNumHertzContacts) = scalKeff;
                end

            end
        end

        % these vectors where initiliazed with zeros
        % trim them down so they only contact the contacts to save space in output file
        vecHertzNN = vecHertzNN(1:scalNumHertzContacts);
        vecHertzMM  = vecHertzMM(1:scalNumHertzContacts);
        vecHertzKeff = vecHertzKeff(1:scalNumHertzContacts);

        fprintf('Hertzian contacts: %d | mean k_eff=%.4f | min=%.4f | max=%.4f\n', ...
            scalNumHertzContacts, mean(vecHertzKeff), min(vecHertzKeff), max(vecHertzKeff));
    end

%% Final plot (AFTER cleanRats — backbone only, rattlers removed)
    figure;
    hold on;
    if boolThreeD
        [sx, sy, sz] = sphere(16);

        % Main particles
        for np = 1:N_clean
            r = vecDiameter(np)/2;
            surf(r*sx + vecPosX(np), r*sy + vecPosY(np), r*sz + vecPosZ(np), ...
                'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.6);
        end

        % Ghost particles on +x, +y, +z faces
        vecOffsets = [scalBoxWidthX, 0, 0; ...
                      0, scalBoxHeightY, 0; ...
                      0, 0, scalBoxDepthZ];  % [3 x 3] one offset per face

        for iface = 1:3
            ox = vecOffsets(iface, 1);
            oy = vecOffsets(iface, 2);
            oz = vecOffsets(iface, 3);
            for np = 1:N_clean
                r = vecDiameter(np)/2;
                surf(r*sx + vecPosX(np) + ox, ...
                     r*sy + vecPosY(np) + oy, ...
                     r*sz + vecPosZ(np) + oz, ...
                    'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.15);
            end
        end

        axis equal;
        axis([-scalBoxWidthX*0.0 2*scalBoxWidthX ...
              -scalBoxHeightY*0.0 2*scalBoxHeightY ...
              -scalBoxDepthZ*0.0  2*scalBoxDepthZ]);
        axis manual;
        xlabel('x'); ylabel('y'); zlabel('z');
        lighting gouraud;
        camlight;
        rotate3d on;
        title(sprintf('After cleanRats: N=%d (of %d original), mean coord num=%.2f', N_clean, N_original, scalMeanCoordNum));

    else
        % Main particles
        for np = 1:N_clean
            rectangle('Position', [vecPosX(np) - vecDiameter(np)/2, ...
                                    vecPosY(np) - vecDiameter(np)/2, ...
                                    vecDiameter(np), vecDiameter(np)], ...
                'Curvature', [1 1], 'FaceColor', 'b', 'EdgeColor', 'none');
        end

        % Ghost particles: all 8 surrounding tiles
        vecOffsets2D = [scalBoxWidthX,  0; ...
                       -scalBoxWidthX,  0; ...
                        0,  scalBoxHeightY; ...
                        0, -scalBoxHeightY; ...
                        scalBoxWidthX,  scalBoxHeightY; ...
                       -scalBoxWidthX,  scalBoxHeightY; ...
                        scalBoxWidthX, -scalBoxHeightY; ...
                       -scalBoxWidthX, -scalBoxHeightY];

        for iface = 1:8
            ox = vecOffsets2D(iface, 1);
            oy = vecOffsets2D(iface, 2);
            for np = 1:N_clean
                rectangle('Position', [vecPosX(np) + ox - vecDiameter(np)/2, ...
                                        vecPosY(np) + oy - vecDiameter(np)/2, ...
                                        vecDiameter(np), vecDiameter(np)], ...
                    'Curvature', [1 1], 'FaceColor', 'r', 'EdgeColor', 'none', ...
                    'FaceAlpha', 0.15);
            end
        end

        axis equal;
        axis([-scalBoxWidthX 2*scalBoxWidthX -scalBoxHeightY 2*scalBoxHeightY]);
        title(sprintf('After cleanRats: N=%d (of %d original), mean coord num=%.2f', N_clean, N_original, scalMeanCoordNum));
    end
    drawnow;
    hold off;

    % Export the after-cleanRats packing photo (backbone only)
    try
        strPngTag = '';
        if boolFrictionOn
            strPngTag = '_Fric';
        end
        strAfterPng = [strFilename(1:end-4) strPngTag '_AfterCleanRats.png'];
        print(gcf, strAfterPng, '-dpng', '-r120');
        fprintf('Packing photo saved to: %s\n', strAfterPng);
    catch ME
        warning('pack:PNGExportFailed', 'Could not export after-cleanRats plot: %s', ME.message);
    end

%% Save results

    % if boolThreeD
    %     scalRoundedWidth = round(N^(1/3));
    %     strFilename = sprintf('%s3D_N%d_P%s_Width%d_Seed%d.mat', ...
    %         save_path, N, num2str(P_target), scalRoundedWidth, seed);
    % else
    %     scalRoundedWidth = round(sqrt(N));
    %     strFilename = sprintf('%s2D_N%d_P%s_Width%d_Seed%d.mat', ...
    %         save_path, N, num2str(P_target), scalRoundedWidth, seed);
    % end

    N_original = N;
    N_clean = size(matPositions, 1);
    N = N_clean;

    % Compute packing fraction after cleanRats
    if boolThreeD
        scalPackingFraction = computePackingFraction(vecDiameter, scalBoxWidthX, scalBoxHeightY, scalBoxDepthZ);
    else
        scalPackingFraction = computePackingFraction(vecDiameter, scalBoxWidthX, scalBoxHeightY);
    end
    fprintf('Packing fraction after cleanRats: %.4f\n', scalPackingFraction);

    if calc_eig
        % 3D hessian not yet implemented — skip eigenmodes
        if boolThreeD
            warning('calc_eig not supported for 3D yet — saving positions only.');
            if options.hertzian
                save(strFilename, 'vecPosX', 'vecPosY', 'vecPosZ', 'vecDiameter', ...
                    'scalBoxWidthX', 'scalBoxHeightY', 'scalBoxDepthZ', ...
                    'K', 'P_target', 'scalPressure', 'N', 'N_original', 'scalPackingFraction', 'scalPackingFractionFull', 'scalMeanCoordNum', ...
                    'vecHertzNN', 'vecHertzMM', 'vecHertzKeff', 'boolFrictionOn', 'scalMu', 'scalKt');
            else
                save(strFilename, 'vecPosX', 'vecPosY', 'vecPosZ', 'vecDiameter', ...
                    'scalBoxWidthX', 'scalBoxHeightY', 'scalBoxDepthZ', ...
                    'K', 'P_target', 'scalPressure', 'N', 'N_original', 'scalPackingFraction', 'scalPackingFractionFull', 'scalMeanCoordNum', 'boolFrictionOn', 'scalMu', 'scalKt');
            end
        else
            matPositions = [vecPosX, vecPosY];
            vecRadii = vecDiameter ./ 2;
            [matPositions, vecRadii] = cleanRats(matPositions, vecRadii, K, scalBoxHeightY, scalBoxWidthX);
            matHessian = hess2d(matPositions, vecRadii, K, scalBoxHeightY, scalBoxWidthX);
            [matEigenVectors, matEigenValues] = eig(matHessian);
            if options.hertzian
                save(strFilename, 'vecPosX', 'vecPosY', 'vecDiameter', ...
                    'scalBoxWidthX', 'scalBoxHeightY', 'K', 'P_target', 'scalPressure', 'N', 'N_original', ...
                    'scalPackingFraction', 'scalPackingFractionFull', 'scalMeanCoordNum', 'matEigenVectors', 'matEigenValues', ...
                    'vecHertzNN', 'vecHertzMM', 'vecHertzKeff', 'boolFrictionOn', 'scalMu', 'scalKt');
            else
                save(strFilename, 'vecPosX', 'vecPosY', 'vecDiameter', ...
                    'scalBoxWidthX', 'scalBoxHeightY', 'K', 'P_target', 'scalPressure', 'N', 'N_original', ...
                    'scalPackingFraction', 'scalPackingFractionFull', 'scalMeanCoordNum', 'matEigenVectors', 'matEigenValues', 'boolFrictionOn', 'scalMu', 'scalKt');
            end
        end
    else
        if boolThreeD
            if options.hertzian
                save(strFilename, 'vecPosX', 'vecPosY', 'vecPosZ', 'vecDiameter', ...
                    'scalBoxWidthX', 'scalBoxHeightY', 'scalBoxDepthZ', ...
                    'K', 'P_target', 'scalPressure', 'N', 'N_original', 'scalPackingFraction', 'scalPackingFractionFull', 'scalMeanCoordNum', ...
                    'vecHertzNN', 'vecHertzMM', 'vecHertzKeff', 'boolFrictionOn', 'scalMu', 'scalKt');
            else
                save(strFilename, 'vecPosX', 'vecPosY', 'vecPosZ', 'vecDiameter', ...
                    'scalBoxWidthX', 'scalBoxHeightY', 'scalBoxDepthZ', ...
                    'K', 'P_target', 'scalPressure', 'N', 'N_original', 'scalPackingFraction', 'scalPackingFractionFull', 'scalMeanCoordNum', 'boolFrictionOn', 'scalMu', 'scalKt');
            end
        else
            if options.hertzian
                save(strFilename, 'vecPosX', 'vecPosY', 'vecDiameter', ...
                    'scalBoxWidthX', 'scalBoxHeightY', 'K', 'P_target', 'scalPressure', ...
                    'N', 'N_original', 'scalPackingFraction', 'scalPackingFractionFull', 'scalMeanCoordNum', ...
                    'vecHertzNN', 'vecHertzMM', 'vecHertzKeff', 'boolFrictionOn', 'scalMu', 'scalKt');
            else
                save(strFilename, 'vecPosX', 'vecPosY', 'vecDiameter', ...
                    'scalBoxWidthX', 'scalBoxHeightY', 'K', 'P_target', 'scalPressure', ...
                    'N', 'N_original', 'scalPackingFraction', 'scalPackingFractionFull', 'scalMeanCoordNum', 'boolFrictionOn', 'scalMu', 'scalKt');
            end
        end
    end

    fprintf('File saved to: %s\n', strFilename);

    %% Repeat-tile the packing in x, y, and/or z for superlattice packings.
    %   The tile SOURCE is the FULL jammed state captured before cleanRats
    %   (rattlers included), NOT the backbone saved to strFilename: the
    %   user-visible superlattice repeats the exact simulated state, and
    %   rattler positions must line up across the periodic boundary.
    if boolThreeD
        if x_mult ~= 1 || y_mult ~= 1 || z_mult ~= 1
            if ~options.saveFullState
                error('pack:NeedFullStateFor3DTile', ...
                    ['3D repeat-tile requires the full (pre-cleanRats) state. ', ...
                     'Call pack(..., options) with options.saveFullState = true.']);
            end
            [vecPosXFinal, vecPosYFinal, vecPosZFinal, vecDiameterFinal, ...
                scalBoxWidthXTiled, scalBoxHeightYFinal, scalBoxDepthZFinal, ...
                scalNTiled] = tile3D( ...
                vecTileSrcX, vecTileSrcY, vecTileSrcZ, vecTileSrcD, ...
                scalBoxWidthX, scalBoxHeightY, scalBoxDepthZ, ...
                x_mult, y_mult, z_mult, NTileSrc);

            % Metrics on the tiled packing (rattlers included, matching the
            % stored particles). Packing fraction is identical to the base
            % tile by construction; coordination number under full PBC equals
            % the base tile's full-state Zn (tiling replicates the contact
            % network) — recomputed anyway via the shared function so the
            % file is self-consistent with what it stores.
            scalPackingFractionTiled  = computePackingFraction(vecDiameterFinal, scalBoxWidthXTiled, scalBoxHeightYFinal, scalBoxDepthZFinal);
            scalPackingFractionFullTiled = scalPackingFractionTiled;
            scalMeanCoordNumTiled = computeMeanCoordNum(vecPosXFinal, vecPosYFinal, vecDiameterFinal, ...
                scalBoxWidthXTiled, scalBoxHeightYFinal, vecPosZFinal, scalBoxDepthZFinal);
            fprintf('3D tiled packing: N=%d (of %d per tile), PF=%.4f, mean coordination number=%.4f\n', ...
                scalNTiled, NTileSrc, scalPackingFractionTiled, scalMeanCoordNumTiled);

            % Tiled output filename: 3D_N%d_P%s_Width%d_Seed%d_TiledX<xm>Y<ym>Z<zm>.mat
            % (backbone file naming is 3D_N%d_P%s_Width%d_Seed%d[_Hertz].mat,
            %  so the _TiledX..Y..Z.. tag keeps the two distinct and records
            %  every multiplier — a pure-z tiling is ..._TiledX1Y1Z9.)
            if options.hertzian
                strTiledFilename = sprintf('%s3D_N%d_P%s_Width%d_Seed%d_TiledX%dY%dZ%d_Hertz.mat', ...
                    save_path, scalNTiled, num2str(P_target), scalRoundedWidth, seed, x_mult, y_mult, z_mult);
            else
                strTiledFilename = sprintf('%s3D_N%d_P%s_Width%d_Seed%d_TiledX%dY%dZ%d.mat', ...
                    save_path, scalNTiled, num2str(P_target), scalRoundedWidth, seed, x_mult, y_mult, z_mult);
            end

            % Alias the canonical variable names to the TILED state so the
            % save() below records the superlattice, not the backbone.
            vecPosX        = vecPosXFinal;
            vecPosY        = vecPosYFinal;
            vecPosZ        = vecPosZFinal;
            vecDiameter    = vecDiameterFinal;
            scalBoxWidthX  = scalBoxWidthXTiled;
            scalBoxHeightY = scalBoxHeightYFinal;
            scalBoxDepthZ  = scalBoxDepthZFinal;
            N              = scalNTiled;
            N_original     = scalNTiled;   % full (pre-cleanRats) tiled state
            scalPackingFraction    = scalPackingFractionTiled;
            scalPackingFractionFull = scalPackingFractionFullTiled;
            scalMeanCoordNum       = scalMeanCoordNumTiled;

            % Save the TILED (superlattice) packing — the full pre-cleanRats
            % state repeated across the tile grid — not the backbone. The
            % variable names match the backbone .mat convention (vecPosX, N,
            % scalPackingFraction, ...) so downstream loaders are unchanged.
            save(strTiledFilename, ...
                'vecPosX', 'vecPosY', 'vecPosZ', 'vecDiameter', ...
                'scalBoxWidthX', 'scalBoxHeightY', 'scalBoxDepthZ', ...
                'K', 'P_target', 'scalPressure', 'N', 'N_original', ...
                'scalPackingFraction', 'scalPackingFractionFull', 'scalMeanCoordNum', ...
                'boolFrictionOn', 'scalMu', 'scalKt', ...
                'x_mult', 'y_mult', 'z_mult', ...
                'vecPosXFinal', 'vecPosYFinal', 'vecPosZFinal', 'vecDiameterFinal', ...
                'scalBoxWidthXTiled', 'scalBoxHeightYFinal', 'scalBoxDepthZFinal', ...
                'scalNTiled', 'NTileSrc', ...
                'scalPackingFractionTiled', 'scalPackingFractionFullTiled', 'scalMeanCoordNumTiled');
            fprintf('3D tiled packing saved to: %s\n', strTiledFilename);
        end
    else
        % 2D path: unchanged — reuse the existing packRepeatTile on the
        % backbone (its documented behavior).
        if x_mult ~= 1 || y_mult ~= 1
            packRepeatTile(N_original, K, P_target, scalRoundedWidth, seed, x_mult, y_mult, ...
                calc_eig, save_path, save_path);
            disp("Tile saved to: " + strFilename);
        end
    end

    toc
end

function [vecPosX, vecPosY, cellParticleList, scalNumCellsX, scalNumCellsY] = ...
        rebuildCellList(vecPosX, vecPosY, scalBoxWidthX, scalBoxHeightY, scalRawCellWidth, scalTimestep, N)

    scalNumCellsX  = round(scalBoxWidthX  / scalRawCellWidth);
    scalCellWidthX = scalBoxWidthX  / scalNumCellsX;
    scalNumCellsY  = round(scalBoxHeightY / scalRawCellWidth);
    scalCellWidthY = scalBoxHeightY / scalNumCellsY;

    % Sanity check: no particle should move more than one box length per step
    vecFloorX = floor(vecPosX / scalBoxWidthX);
    vecFloorY = floor(vecPosY / scalBoxHeightY);
    if any(abs(vecFloorX) > 1) || any(abs(vecFloorY) > 1)
        error(['Particle moved more than one box length in a single timestep.\n' ...
               'Max x overshoot: %.2f box lengths\n' ...
               'Max y overshoot: %.2f box lengths\n' ...
               'Reduce scalTimestep (currently %.4e).'], ...
               max(abs(vecFloorX)), max(abs(vecFloorY)), scalTimestep);
    end

    % Wrap positions into [0, L) before rebuilding
    vecPosX = mod(vecPosX, scalBoxWidthX);   % [N x 1]
    vecPosY = mod(vecPosY, scalBoxHeightY);  % [N x 1]

    % Map to cell indices, clamped to [1, scalNumCells]
    vecCellIdxX = min(max(ceil(vecPosX / scalCellWidthX), 1), scalNumCellsX);  % [N x 1]
    vecCellIdxY = min(max(ceil(vecPosY / scalCellWidthY), 1), scalNumCellsY);  % [N x 1]

    % O(N) bucketing via accumarray
    vecCellLinearIdx = vecCellIdxX + scalNumCellsX * (vecCellIdxY - 1);        % [N x 1]
    cellParticleList = accumarray(vecCellLinearIdx, (1:N)', [scalNumCellsX*scalNumCellsY 1], @(x){x});
    cellParticleList = reshape(cellParticleList, scalNumCellsX, scalNumCellsY); % [scalNumCellsX x scalNumCellsY]

end


function [vecPairIdxSource, vecPairIdxDest, scalNumPairs, scalMaxPairs] = findNeighbors2D( ...
        cellParticleList, scalNumCellsX, scalNumCellsY, ...
        vecPairIdxSource, vecPairIdxDest, scalMaxPairs)

    scalNumPairs = 0;

    for idxCellX = 1:scalNumCellsX
        for idxCellY = 1:scalNumCellsY

            scalCellLeft  = mod(idxCellX-2, scalNumCellsX)+1;
            scalCellRight = mod(idxCellX,   scalNumCellsX)+1;
            scalCellDown  = mod(idxCellY-2, scalNumCellsY)+1;
            scalCellUp    = mod(idxCellY,   scalNumCellsY)+1;

            vecCurrentCellPartIdx = cellParticleList{idxCellX, idxCellY};

            vecNeighborList = [ ...
                cellParticleList{scalCellLeft,  scalCellDown}; ...
                cellParticleList{scalCellLeft,  idxCellY};     ...
                cellParticleList{scalCellLeft,  scalCellUp};   ...
                cellParticleList{idxCellX,      scalCellDown}; ...
                vecCurrentCellPartIdx;                                ...
                cellParticleList{idxCellX,      scalCellUp};   ...
                cellParticleList{scalCellRight, scalCellDown}; ...
                cellParticleList{scalCellRight, idxCellY};     ...
                cellParticleList{scalCellRight, scalCellUp}];

            for idxNN = vecCurrentCellPartIdx'
                vecContactCandidatePartIdx = vecNeighborList(vecNeighborList > idxNN);
                scalNumCandidates = numel(vecContactCandidatePartIdx);
                if scalNumCandidates == 0; continue; end
                if scalNumPairs + scalNumCandidates > scalMaxPairs
                    scalMaxPairs = 2 * scalMaxPairs;
                    vecPairIdxSource(scalMaxPairs) = 0;
                    vecPairIdxDest(scalMaxPairs) = 0;
                end
                vecPairIdxSource(scalNumPairs+1 : scalNumPairs+scalNumCandidates) = idxNN;
                vecPairIdxDest(scalNumPairs+1 : scalNumPairs+scalNumCandidates) = vecContactCandidatePartIdx;
                scalNumPairs = scalNumPairs + scalNumCandidates;
            end
        end
    end
end

function fricState = buildFricState( ...
        vecPosX, vecPosY, vecDiameter, ...
        scalBoxWidthX, scalBoxHeightY, ...
        matDispTan, K, Kt, mu, ...
        gammaNormal, gammaTang, M)
% buildFricState -- Rebuild the final contact graph at ORIGINAL
% particle indices (before cleanRats) and export a fricState struct.
%
% O(N^2) pairwise loop. For large N, replace with a cell-list.
%
% Output struct fields (all on original N particles):
%   vecContactNN / vecContactMM   [Nc x 1]  pair indices (i < j)
%   vecOverlap                    [Nc x 1]  normal overlaps
%   vecFn / vecFt                 [Nc x 1]  normal / tangential forces
%   vecUnitTanX / vecUnitTanY     [Nc x 1]  unit tangent [n_y, -n_x]
%   K, Kt, mu, gammaNormal, gammaTang, M
%   vecRadius / vecInertia        [N x 1]
%   N   original particle count
%   boxLx / boxLy
%   matDispTan  [N x N]  tangential spring displacement history

    N           = size(vecPosX, 1);
    vecRadius       = vecDiameter / 2;
    vecInertia      = 0.5 * M .* (vecRadius .^ 2);

    fricState = struct();
    fricState.vecContactNN = zeros(0,1);
    fricState.vecContactMM = zeros(0,1);
    fricState.vecOverlap     = zeros(0,1);
    fricState.vecFn          = zeros(0,1);
    fricState.vecFt          = zeros(0,1);
    fricState.vecUnitTanX    = zeros(0,1);
    fricState.vecUnitTanY    = zeros(0,1);

    for ii = 1:(N-1)
        for jj = ii+1:N
            dx = vecPosX(jj) - vecPosX(ii);
            dy = vecPosY(jj) - vecPosY(ii);
            dx = dx - scalBoxWidthX * round(dx / scalBoxWidthX);
            dy = dy - scalBoxHeightY * round(dy / scalBoxHeightY);
            dist = sqrt(dx^2 + dy^2);
            if dist < 1e-12, continue; end
            % Contact distance is the SUM OF RADII = (d_i + d_j)/2, NOT the sum
            % of diameters. (vecDiameter(ii)+vecDiameter(jj)) would double-count
            % and report ~2x overlaps, so buildFricState's forces/overlaps are
            % only correct if this is (r_i + r_j).
            D_ij   = (vecDiameter(ii) + vecDiameter(jj)) / 2;   % = r_i + r_j
            delta  = D_ij - dist;
            if delta <= 0, continue; end
            nx = dx / dist;
            ny = dy / dist;
            tx =  ny;
            ty = -nx;
            Fn      = -K * delta;
            DispTan = matDispTan(ii + N * (jj - 1));   % [1 x 1] pair tangential displacement
            Ft      = -Kt * DispTan;
            fricState.vecContactNN = [fricState.vecContactNN; ii];
            fricState.vecContactMM = [fricState.vecContactMM; jj];
            fricState.vecOverlap     = [fricState.vecOverlap;     delta];
            fricState.vecFn          = [fricState.vecFn;          Fn];
            fricState.vecFt          = [fricState.vecFt;          Ft];
            fricState.vecUnitTanX    = [fricState.vecUnitTanX;    tx];
            fricState.vecUnitTanY    = [fricState.vecUnitTanY;    ty];
        end
    end

    fricState.K       = K;
    fricState.Kt      = Kt;
    fricState.mu      = mu;
    fricState.gammaNormal = gammaNormal;
    fricState.gammaTang   = gammaTang;
    fricState.M       = M;
    fricState.vecRadius = vecRadius;
    fricState.vecInertia = vecInertia;
    fricState.N     = N;
    fricState.boxLx = scalBoxWidthX;
    fricState.boxLy = scalBoxHeightY;
    fricState.matDispTan = matDispTan;
end
