function pack(scalNumParticles, scalSpringConstant, scalDiameterSmall, scalDiameterRatio, scalMass, ...
        scalPressureTarget, scalSeed, boolPlotIt, scalXMult, scalYMult, scalZMult, boolCalcEig, strSavePath, sOptions)


    arguments
        scalNumParticles   (1,1) double {mustBeInteger, mustBePositive} = 100     % N, number of particles
        scalSpringConstant (1,1) double {mustBePositive} = 100                    % K, contact spring constant
        scalDiameterSmall  (1,1) double {mustBePositive} = 1                      % D, small-particle diameter
        scalDiameterRatio  (1,1) double {mustBePositive} = 1.4                    % G, large/small diameter ratio
        scalMass           (1,1) double {mustBePositive} = 1                      % M, particle mass
        scalPressureTarget (1,1) double {mustBePositive} = 0.0001                 % P_target, target pressure
        scalSeed           (1,1) double {mustBeInteger, mustBePositive} = 1       % random seed
        boolPlotIt         (1,1) logical = false                                  % draw the packing as it compresses
        scalXMult          (1,1) double = 1                                       % repeat-tile multiplier in x
        scalYMult          (1,1) double = 1                                       % repeat-tile multiplier in y
        scalZMult          (1,1) double = 0                                       % repeat-tile multiplier in z; ~= 0 selects 3D
        boolCalcEig        (1,1) logical = false                                  % also save Hessian eigenmodes (2D only)
        strSavePath        (1,1) string = "./junkyard"                            % prefix of every saved file, e.g. "data/"
        % sOptions is passed as a plain positional struct (the frictional test
        % calls pack(...,strSavePath,opts)). A fully-defaulted struct means the
        % frictionless 13-arg calls get all fields, and a caller may pass a
        % partial struct; missing fields are backfilled below.
        sOptions (1,1) struct = struct('hertzian', false, ...
            'flagFrictionOn', false, ...
            'scalFricCoef', 0.50, ...
            'scalTangentialK', 1/3, ...
            'scalGammaNormal', 0, ...
            'scalGammaTangential', 0, ...
            'saveFrictionalState', false, ...
            'saveFullState', true)
    end

     % Backfill any option fields a caller omitted so both the frictionless
     % (13-arg) and frictional (14-arg with partial opts) call styles work.
    if ~isfield(sOptions, 'hertzian')
        sOptions.hertzian = false;
    end
    if ~isfield(sOptions, 'flagFrictionOn')
        sOptions.flagFrictionOn = false;
    end
    if ~isfield(sOptions, 'scalFricCoef')
        sOptions.scalFricCoef = 0.50;
    end
    if ~isfield(sOptions, 'scalTangentialK')
        sOptions.scalTangentialK = 1/3;
    end
    if ~isfield(sOptions, 'scalGammaNormal')
        sOptions.scalGammaNormal = 0;
    end
    if ~isfield(sOptions, 'scalGammaTangential')
        sOptions.scalGammaTangential = 0;
    end
    if ~isfield(sOptions, 'saveFrictionalState')
        sOptions.saveFrictionalState = false;
    end
    if ~isfield(sOptions, 'saveFullState')
        sOptions.saveFullState = true;
    end
    boolHertzian      = sOptions.hertzian;       % Hertzian (vs. linear) contact law
    boolSaveFullState = sOptions.saveFullState;  % also save the pre-cleanRats tile source

    % check to see if 3d path is needed
    boolThreeD = (scalZMult ~= 0);

    %% Guard: Cundall-Strack friction is implemented for 2D packings only.
    %% 3D friction (rotation about 3 axes) requires a different model.
    if boolThreeD && sOptions.flagFrictionOn
        error('pack:Friction3DNotSupported', 'flagFrictionOn = true is 2D only.');
    end


    % Check if packing already exists — skip if so
    if boolThreeD
        scalRoundedWidth = round(scalNumParticles^(1/3));
        if boolHertzian
            strFilename = sprintf('%s3D_N%d_P%s_Width%d_Seed%d_Hertz.mat', ...
                strSavePath, scalNumParticles, num2str(scalPressureTarget), scalRoundedWidth, scalSeed);
        else
            strFilename = sprintf('%s3D_N%d_P%s_Width%d_Seed%d.mat', ...
                strSavePath, scalNumParticles, num2str(scalPressureTarget), scalRoundedWidth, scalSeed);
        end
        % Full-state (pre-cleanRats) filename for the repeat-tile source file.
        % A '_Full' tag keeps it distinct from the backbone .mat so re-running
        % pack on the same parameters does not clobber the tile source.
        if boolHertzian
            strFullFilename = sprintf('%s3D_N%d_P%s_Width%d_Seed%d_Full_Hertz.mat', ...
                strSavePath, scalNumParticles, num2str(scalPressureTarget), scalRoundedWidth, scalSeed);
        else
            strFullFilename = sprintf('%s3D_N%d_P%s_Width%d_Seed%d_Full.mat', ...
                strSavePath, scalNumParticles, num2str(scalPressureTarget), scalRoundedWidth, scalSeed);
        end
    else
        scalRoundedWidth = round(sqrt(scalNumParticles));
        if boolHertzian
            strFilename = sprintf('%s2D_N%d_P%s_Width%d_Seed%d_Hertz.mat', ...
                strSavePath, scalNumParticles, num2str(scalPressureTarget), scalRoundedWidth, scalSeed);
        else
            strFilename = sprintf('%s2D_N%d_P%s_Width%d_Seed%d.mat', ...
                strSavePath, scalNumParticles, num2str(scalPressureTarget), scalRoundedWidth, scalSeed);
        end
        if boolHertzian
            strFullFilename = sprintf('%s2D_N%d_P%s_Width%d_Seed%d_Full_Hertz.mat', ...
                strSavePath, scalNumParticles, num2str(scalPressureTarget), scalRoundedWidth, scalSeed);
        else
            strFullFilename = sprintf('%s2D_N%d_P%s_Width%d_Seed%d_Full.mat', ...
                strSavePath, scalNumParticles, num2str(scalPressureTarget), scalRoundedWidth, scalSeed);
        end
    end
    if isfile(strFilename)
        fprintf('Packing already exists, skipping: %s\n', strFilename);
        return;
    end

    tic

    rng(scalSeed)

%% Box and particle setup
    if boolThreeD
        scalBoxWidthX  = 2*scalNumParticles^(1/3)*scalDiameterSmall;   % Lx
        scalBoxHeightY = 2*scalNumParticles^(1/3)*scalDiameterSmall;   % Ly
        scalBoxDepthZ  = 2*scalNumParticles^(1/3)*scalDiameterSmall;   % Lz
    else
        scalBoxWidthX  = 2*sqrt(scalNumParticles)*scalDiameterSmall;   % Lx
        scalBoxHeightY = 2*sqrt(scalNumParticles)*scalDiameterSmall;   % Ly
    end

    scalDissipationVelocity = 0.1;  % Bv: velocity-dependent dissipation prefactor
    scalDissipationAbsolute = 0.5;  % B:  absolute (drag) dissipation
    scalTemperature = 1;    % T:  initial velocity scale

    %% Equal number of small and large particles
    scalNumSmall = scalNumParticles/2;

    % Assign diameters: smallest half get D, largest half get D*G
    [~, vecSortIdx] = sort(rand(scalNumParticles, 1)); % [N x 1] randomize the particle indices
    vecDiameter = scalDiameterSmall * scalDiameterRatio * ones(scalNumParticles, 1); % [N x 1] default all to large
    vecDiameter(vecSortIdx(1:scalNumSmall)) = scalDiameterSmall; % overwrite bottom half with small

    if boolHertzian
        vecRadii = vecDiameter / 2; % [N x 1]
    end

    % Precompute contact-distance matrix: matContactDist(i,j) = (r_i + r_j)
    % Avoids recomputing inside the force loop every timestep
    matContactDist = (vecDiameter + vecDiameter') / 2;  % [N x N]

    %% Physical parameters
    scalGravity = 0;
    scalPressure = 0;
    scalPressureFastGrow = scalPressureTarget / 50;
    scalCompressionRate = scalPressureTarget;
    %% Frictional (Cundall-Strack) compression protocol — 2D only.
    %%
    %%  Friction is ON for the whole compression (no frictionless
    %% pre-compression). The box is resized only once the grains have RELAXED
    %% after the previous resize, so every decision reads a settled pressure
    %% (the frictionless path likewise only compresses when Ek < 1e-8):
    %% compress when P is below a dead-band around P_target, expand when P is
    %% above it, HOLD while P is in-band. The strain step starts at the
    %% frictionless fast rate and is HALVED every time the controller reverses
    %% direction (it overshot P_target), i.e. a bisection on the box size, so
    %% P lands in the band instead of bouncing across it. The old controller
    %% resized by 0.5% every step while P was out of band; one such step moves
    %% P by several P_target, which is the limit cycle of issue #16.
    %%  Convergence (unchanged) requires FORCE BALANCE — max_i |F_net,i| < tol *
    %% mean contact force — sustained for several hundred consecutive in-band
    %% steps with the box held and a percolating contact network
    %% (mean Zn >= scalFrictionZmin), as OverDamp.cpp's Acc_max < Fthresh.
    scalFrictionDeadBand     = 0.15;        % P in-band: |P - P_target|/P_target < 15%
    scalFrictionForceTol     = 0.01;        % accept when max|F_net|/mean|F_contact| < 1%
    scalFrictionBalCount      = 0;          % consecutive in-band steps satisfying force balance
    scalFrictionBalWindow     = 300;        % sustained force-balance steps to accept
    scalFrictionZmin          = 2.5;        % percolation guard: mean Zn above this to accept
                                             % (frictional 2D isostatic z_iso = 3)
    scalFrictionMaxSteps      = 3e6;        % hard cap for the frictional phase (safety only)
    scalFrictionStrain        = 0.01;       % current box strain per resize (starts at the
                                             % frictionless fast rate, halved on each reversal)
    scalFrictionStrainMin     = 1e-2 * scalPressureTarget;  % floor for the bisected strain step
    scalFrictionRelaxSteps    = 100;        % minimum steps between resizes (one contact
                                             % oscillation period: dt = period/100)
    scalFrictionRelaxKE       = 1e-2;       % relaxed once kinetic energy per grain (incl.
                                             % rotation) < this * elastic energy per grain at
                                             % P_target, or once in force balance
    scalFrictionLastDir       = 0;          % last resize: -1 compress, +1 expand, 0 none yet
    scalFrictionLastResize    = 0;          % step of the last resize
    % elastic energy per grain at which the pressure estimate below reads P_target
    if boolHertzian
        scalFrictionEpTarget = (2/5) * scalSpringConstant * scalPressureTarget^(5/2);
    else
        scalFrictionEpTarget = 0.5 * scalSpringConstant * scalPressureTarget^2;
    end
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
    boolFrictionOn          = sOptions.flagFrictionOn;
    scalMu                  = sOptions.scalFricCoef;
    scalKtOverK             = sOptions.scalTangentialK;
    scalKt                  = scalSpringConstant * scalKtOverK;
    scalGammaNormal          = sOptions.scalGammaNormal;
    scalGammaTangential     = sOptions.scalGammaTangential;
    boolSaveFricState       = sOptions.saveFrictionalState;



%% Display / simulation parameters
    boolPlotKE = false;
    scalPlotSkip = 200;   % timesteps between plot updates
    scalCellUpdateInterval = 1;

    % time step should be 1/100 of a particle oscillation period
    % For Hertzian, k_eff ~ 2K*sqrt(Reff*delta) is unknown at init,
    % so use a conservative estimate based on largest R and largest delta
    % largest Reff ~ G*D/4 (large-large contact)
    % and largest delta ~ 1e-1 (max target pressure)
    if boolHertzian
        scalTimestep = 2*pi * sqrt(scalMass / (2*scalSpringConstant*sqrt(scalDiameterRatio*scalDiameterSmall/4 .* 1E-1))) * 0.01;
    else
        scalTimestep = 2*pi * sqrt(scalMass/scalSpringConstant) * 0.01;
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
    [vecPosX, vecPosY, vecPosZ] = ndgrid(scalDiameterSmall/2 : scalDiameterRatio*scalDiameterSmall : scalBoxWidthX-scalDiameterSmall/2, ...
                                          scalDiameterSmall/2 : scalDiameterRatio*scalDiameterSmall : scalBoxHeightY-scalDiameterSmall/2, ...
                                          scalDiameterSmall/2 : scalDiameterRatio*scalDiameterSmall : scalBoxDepthZ-scalDiameterSmall/2);
    else
        [vecPosX, vecPosY] = ndgrid(scalDiameterSmall/2 : scalDiameterRatio*scalDiameterSmall : scalBoxWidthX-scalDiameterSmall/2, ...
                                    scalDiameterSmall/2 : scalDiameterRatio*scalDiameterSmall : scalBoxHeightY-scalDiameterSmall/2);
    end

    % shuffle the particles to avoid crystallization
    % shuffle the indices of the particles
    [~, vecShuffleIdx] = sort(rand(numel(vecPosX), 1));  % [numel x 1]
    vecPosX = vecPosX(vecShuffleIdx(1:scalNumParticles));   % [N x 1] of particle x-positions
    vecPosY = vecPosY(vecShuffleIdx(1:scalNumParticles));   % [N x 1] etc
    if boolThreeD
        vecPosZ = vecPosZ(vecShuffleIdx(1:scalNumParticles));
    end

    % assign random initial velocities
    vecVelX = sqrt(scalTemperature) * randn(scalNumParticles, 1);  % [N x 1]
    vecVelX = vecVelX - mean(vecVelX);
    vecVelY = sqrt(scalTemperature) * randn(scalNumParticles, 1);  % [N x 1]
    vecVelY = vecVelY - mean(vecVelY);
    if boolThreeD
        vecVelZ = sqrt(scalTemperature) * randn(scalNumParticles, 1);
        vecVelZ = vecVelZ - mean(vecVelZ);  % [N x 1]
    end

    % start with zero accelerations
    vecAccelXPrev = zeros(scalNumParticles, 1);  % [N x 1]
    vecAccelYPrev = zeros(scalNumParticles, 1);  % [N x 1]
    if boolThreeD
        vecAccelZPrev = zeros(scalNumParticles, 1);  % [N x 1]
    end

    %% Rotational DOFs for 2D disks (out-of-plane z rotation)
    %%   vecOmega        = angular velocity omega_i [N x 1]
    %%   vecAlphaPrev    = previous angular acceleration, Verlet half-step.
    %%   Solid disk: I_i = (M_i * r_i^2) / 2,  alpha = torque / I
    if boolFrictionOn
        vecOmega       = zeros(scalNumParticles, 1);
        vecAlphaPrev   = zeros(scalNumParticles, 1);
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
        matDispTan       = zeros(scalNumParticles, scalNumParticles);
        matDispTanStuck  = false(scalNumParticles, scalNumParticles);
    end

%% Verlet cell list setup
    % Determine cell size rounded to be at least 1*G*D
    % to avoid missing interactions, I tried this out 
    % many times and 3 works best
    scalRawCellWidth = 3 * scalDiameterRatio * scalDiameterSmall; % Changing this will mess up findNeighbors2D. Update findNeighbors2D to be like 3D version (so it doesn't double count neighbors) before changing this.

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
    % cellParticleList = accumarray(vecCellLinearIdx, (1:scalNumParticles)', [scalNumCellsX*scalNumCellsY 1], @(vecCellMembers){vecCellMembers});

    % reshape just converts that flat list back into a 2D grid so you can look up neighbors naturally by (ix, iy) index.
    % cellParticleList = reshape(cellParticleList, scalNumCellsX, scalNumCellsY);  % [scalNumCellsX x scalNumCellsY]

    %% Setup plotting
    % boolPlotIt shows only the packing, redrawn every scalPlotSkip steps: in
    % figure 1, or for 2D frictional runs in the compression-GIF figure below
    % (disks with a diameter line that shows rotation).
    boolFricGif     = boolPlotIt && boolFrictionOn && ~boolThreeD;
    boolPlotPacking = boolPlotIt && ~boolFricGif;
    if boolPlotPacking
        figure(1), clf;
        hPlotHandles = gobjects(scalNumParticles, 1);  % [N x 1]
        for idxParticle = 1:scalNumParticles
            hPlotHandles(idxParticle) = rectangle( ...
                'Position',  [vecPosX(idxParticle) - 0.5*vecDiameter(idxParticle), ...
                               vecPosY(idxParticle) - 0.5*vecDiameter(idxParticle), ...
                               vecDiameter(idxParticle), vecDiameter(idxParticle)], ...
                'Curvature', [1 1], 'EdgeColor', 'b');
        end
        axis equal; axis([0 scalBoxWidthX 0 scalBoxHeightY]);
        hAxPacking = gca;
    end

    %% Frictional compression movie (2D, friction on, boolPlotIt = true)
    %  Every scalPlotSkip steps the packing is drawn with a line across each
    %  disk's diameter at its accumulated rotation angle vecTheta, and the
    %  frame is appended to an animated GIF next to the .mat, so grain
    %  rotation and the approach to the converged packing can be checked by
    %  eye. vecTheta is only integrated while this movie is being recorded.
    if boolFricGif
        vecTheta = zeros(scalNumParticles, 1);          % [N x 1] rotation angle (rad, counter-clockwise +)
        strGifFilename = [strFilename(1:end-4) '_Fric_Compression.gif'];
        hFigGif = figure('Color', 'w', 'Name', 'Frictional compression', 'Position', [100 100 600 640]);
        vecGifFrameSize = [];            % [rows cols] of the first frame; later frames match it
        [vecGifFrameSize, boolFricGif] = writeFricGifFrame(hFigGif, strGifFilename, vecGifFrameSize, ...
            vecPosX, vecPosY, vecDiameter, vecTheta, scalBoxWidthX, scalBoxHeightY, ...
            fricGifTitle('Frictional compression', 0, scalPressure, scalPressureTarget, vecDiameter, scalBoxWidthX, scalBoxHeightY, 0, scalMu), 0.1);
    end

    %% Main time-integration loop
    scalLastCompressStep = 0;

    % Pre-allocate pair buffers — N*12 is safe upper bound for 2D jamming
    % because for 2D jamming, the maximum number of contacts is 6N.
    % so twice that is a safe upper bound for 2D jamming
    % we're goingt to store pairs like this:
    % vecPairIdxSource = [1, 1, 2, 3, ...]   <- first particle of each pair
    % vecPairIdxDest = [2, 3, 3, 4, ...]   <- second particle of each pair
    vecPairIdxSource = zeros(scalNumParticles*12, 1);  % [N*12 x 1]
    vecPairIdxDest = zeros(scalNumParticles*12, 1);  % [N*12 x 1]
    scalMaxPairs = scalNumParticles*12;

    fprintf('Starting main integration loop (max %d steps)...\n', scalMaxSteps);
    scalLogInterval = round(0.05 * scalMaxSteps);  % 5% of max steps
    for idxStep = 1:scalMaxSteps

        % Progress logging
        if mod(idxStep, 5000) == 0
            fprintf('  step %d | P=%.4e | P_target=%.4e | P/P_target=%.4f\n', ...
                idxStep, scalPressure, scalPressureTarget, scalPressure/scalPressureTarget);
        end

        %% Plotting: redraw the packing
        if boolPlotPacking && mod(idxStep, scalPlotSkip) == 0
            % if the figure was closed, stop redrawing but keep packing
            boolPlotPacking = all(ishandle(hPlotHandles));
            if boolPlotPacking
                for idxParticle = 1:scalNumParticles
                    set(hPlotHandles(idxParticle), 'Position', ...
                        [vecPosX(idxParticle) - 0.5*vecDiameter(idxParticle), ...
                         vecPosY(idxParticle) - 0.5*vecDiameter(idxParticle), ...
                         vecDiameter(idxParticle), vecDiameter(idxParticle)]);
                end
                axis(hAxPacking, [0 scalBoxWidthX 0 scalBoxHeightY]);
                title(hAxPacking, sprintf('step %d, P/P_{target} = %.3f, L_y = %.4f', ...
                    idxStep, scalPressure / scalPressureTarget, scalBoxHeightY));
                drawnow;
            end
        elseif boolPlotKE && mod(idxStep, scalPlotSkip) == 0
            figure(1), plot(vecPosX, vecPosY, 'k.'); drawnow;
        end
        if boolFricGif && mod(idxStep, scalPlotSkip) == 0
            [vecGifFrameSize, boolFricGif] = writeFricGifFrame(hFigGif, strGifFilename, vecGifFrameSize, ...
                vecPosX, vecPosY, vecDiameter, vecTheta, scalBoxWidthX, scalBoxHeightY, ...
                fricGifTitle('Frictional compression', idxStep, scalPressure, scalPressureTarget, vecDiameter, scalBoxWidthX, scalBoxHeightY, scalMeanCoordNum, scalMu), 0.1);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% First step in Verlet integration %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        vecPosX = vecPosX + vecVelX*scalTimestep + vecAccelXPrev.*(scalTimestep^2/2);
        vecPosY = vecPosY + vecVelY*scalTimestep + vecAccelYPrev.*(scalTimestep^2/2);
        if boolThreeD
            vecPosZ = vecPosZ + vecVelZ*scalTimestep + vecAccelZPrev.*(scalTimestep^2/2);
        end
        if boolFricGif
            % rotation angle, same Verlet position step as x and y
            vecTheta = vecTheta + vecOmega*scalTimestep + vecAlphaPrev.*(scalTimestep^2/2);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Re-assign particles to cells %%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if boolCellUpdateNeeded || mod(idxStep, scalCellUpdateInterval) == 0
            if boolThreeD
                [vecPosX, vecPosY, vecPosZ, cellParticleList, scalNumCellsX, scalNumCellsY, scalNumCellsZ] = ...
                    rebuildCellList3D(vecPosX, vecPosY, vecPosZ, scalBoxWidthX, scalBoxHeightY, scalBoxDepthZ, scalRawCellWidth, scalTimestep, scalNumParticles);
            else
                [vecPosX, vecPosY, cellParticleList, scalNumCellsX, scalNumCellsY] = ...
                    rebuildCellList(vecPosX, vecPosY, scalBoxWidthX, scalBoxHeightY, scalRawCellWidth, scalTimestep, scalNumParticles);
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

        vecContactDist = matContactDist(vecActivePairSource + scalNumParticles*(vecActivePairDest-1));

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

        if boolHertzian
            vecRadiiNN = vecRadii(vecContactNN); % [scalNumContacts x 1]
            vecRadiiMM = vecRadii(vecContactMM); % [scalNumContacts x 1]
            % this falls out of the math for two parabaloids https://en.wikipedia.org/wiki/Contact_mechanics
            vecRadiiEff = (vecRadiiNN .* vecRadiiMM) ./ (vecRadiiNN + vecRadiiMM); % [scalNumContacts x 1] 
        end

        vecSepDist = sqrt(vecSepDistSq); % [scalNumContacts x 1]
        vecOverlap = vecContactDist - vecSepDist; % [scalNumContacts x 1] positive when overlapping

        % Force magnitude and potential energy per contact
        if boolHertzian
            vecForceMag = -(4/3) .* scalSpringConstant .* sqrt(vecRadiiEff) .* vecOverlap.^(3/2);
            vecPotentialContact = (4/3) .* (2/5) * scalSpringConstant .* sqrt(vecRadiiEff) .* vecOverlap.^(5/2);
        else
            vecForceMag = -scalSpringConstant .* vecOverlap; % [scalNumContacts x 1]
            vecPotentialContact = 0.5 * scalSpringConstant .* vecOverlap.^2; % [scalNumContacts x 1]
        end

        % Unit vectors along contact normal
        vecNormalX = vecSepX ./ vecSepDist; % [scalNumContacts x 1]
        vecNormalY = vecSepY ./ vecSepDist; % [scalNumContacts x 1]
        if boolThreeD
            vecNormalZ = vecSepZ ./ vecSepDist;
        end

        % Velocity-dependent dissipation projected onto contact normal
        scalReducedMass = scalMass / 2; % may need to change this for non-uniform mass if the future
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
        vecForceX = accumarray(vecContactNN, vecForceContactX, [scalNumParticles 1]) ...
                  - accumarray(vecContactMM, vecForceContactX, [scalNumParticles 1]);
        vecForceY = accumarray(vecContactNN, vecForceContactY, [scalNumParticles 1]) ...
                  - accumarray(vecContactMM, vecForceContactY, [scalNumParticles 1]);
        if boolThreeD
            vecForceZ = accumarray(vecContactNN, vecForceContactZ, [scalNumParticles 1]) ...
                      - accumarray(vecContactMM, vecForceContactZ, [scalNumParticles 1]);
        end

        
              %% =============================
              %%  Cundall-Strack Tangential Friction
              %%  Ref: CundallStrack_2D/Energy_Disk_VL.cpp
              %%
              %%  Tangential spring coord U_ij evolves:  dU/dt = v_t^total
              %%  Capped by Coulomb:  |F_t| = K_t*|disp|  <=  mu*|F_n|
              %%  Tangential viscous damping (optional):
              %%    F_t^visc = -gamma_t * m_red * v_t^total
              %%  Torque:  tau_i = -r_i * F_t,  tau_j = -r_j * F_t
              %%  t_hat = (n_y, -n_x)  normal rotated -90 degrees
              %% =============================
        if boolFrictionOn && ~boolThreeD
            vecTorque  = zeros(scalNumParticles, 1);
        end

        if boolFrictionOn && ~boolThreeD && scalNumContacts > 0

            % Per-contact state lookup from the [N x N] tangential-displacement matrix
            vecLinIdx = vecContactNN + scalNumParticles * (vecContactMM - 1);   % [scalNumContacts x 1] linear index into matDispTan
            vecDispTan = matDispTan(vecLinIdx);                 % [scalNumContacts x 1] per-contact tangential displacement
            vecStuck   = matDispTanStuck(vecLinIdx);            % [scalNumContacts x 1] per-contact "has-contacted" flag
            vecDispTan(~vecStuck) = 0;

            % Unit tangent: normal rotated -90 degrees in 2D,  t_hat = (n_y, -n_x)
            vecUnitTanX =  vecNormalY;   % [scalNumContacts x 1] unit tangent x  ( = n_y)
            vecUnitTanY = -vecNormalX;   % [scalNumContacts x 1] unit tangent y  ( = -n_x)

            % Slip is evaluated with the HALF-STEP velocities v + a*dt/2 and
            % omega + alpha*dt/2, i.e. the same increment the Verlet position
            % update just applied (x(t+dt) - x(t) = (v + a*dt/2)*dt), so the
            % spring follows the actual relative motion of the contact points.
            % Using v(t) alone lags the spring by dt/2, which acts as NEGATIVE
            % damping (~ K_t*dt/2 per contact): grains with 4-5+ contacts then
            % vibrate and creep forever and the packing never reaches force
            % balance.
            vecVelXHalf  = vecVelX  + vecAccelXPrev * (scalTimestep/2);   % [N x 1]
            vecVelYHalf  = vecVelY  + vecAccelYPrev * (scalTimestep/2);   % [N x 1]
            vecOmegaHalf = vecOmega + vecAlphaPrev  * (scalTimestep/2);   % [N x 1]

            % Total tangential slip rate (velocity): (v_i - v_j) . t_hat
            vecVelTan = ...                   % [scalNumContacts x 1] total tangential slip rate
                (vecVelXHalf(vecContactNN) - vecVelXHalf(vecContactMM)) .* vecUnitTanX + ...
                (vecVelYHalf(vecContactNN) - vecVelYHalf(vecContactMM)) .* vecUnitTanY;

            % Rotational slip rate: -(omega_i*r_i + omega_j*r_j)
            vecRi = vecDiameter(vecContactNN) / 2;   % [scalNumContacts x 1] radius of grain i
            vecRj = vecDiameter(vecContactMM) / 2;   % [scalNumContacts x 1] radius of grain j
            vecVelTan = vecVelTan - (...     % subtract rotational contribution
                vecOmegaHalf(vecContactNN) .* vecRi + vecOmegaHalf(vecContactMM) .* vecRj);

            % Advance tangential displacement:  disp(t+dt) = disp(t) + vel(t+dt/2)*dt
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
                vecFtMag = vecFtMag - (scalGammaTangential * (scalMass / 2)) .* vecVelTan;
            end

            % Distribute tangential force via Newton's 3rd law
            vecFtX = vecFtMag .* vecUnitTanX;   % [scalNumContacts x 1] tangential force x-component
            vecFtY = vecFtMag .* vecUnitTanY;   % [scalNumContacts x 1] tangential force y-component
            vecForceX = vecForceX + ...
                accumarray(vecContactNN, vecFtX, [scalNumParticles 1]) - ...
                accumarray(vecContactMM, vecFtX, [scalNumParticles 1]);
            vecForceY = vecForceY + ...
                accumarray(vecContactNN, vecFtY, [scalNumParticles 1]) - ...
                accumarray(vecContactMM, vecFtY, [scalNumParticles 1]);

            % Contact torques: tau_i = -r_i * F_t,  tau_j = -r_j * F_t.
            % The force on i, F_t*t_hat, acts at r_i*n_hat from its centre and
            % (n_hat x t_hat)_z = -1; grain j gets -F_t*t_hat at -r_j*n_hat.
            % This is the sign consistent with the slip rate above
            % (-(omega_i*r_i + omega_j*r_j)): the tangential spring then
            % stores/returns energy exactly, while a +r*F_t torque does work
            % 2*F_t*(r_i*omega_i + r_j*omega_j) and spins the grains up.
            vecTorque = -(accumarray(vecContactNN, vecRi .* vecFtMag, [scalNumParticles 1]) + ...
                          accumarray(vecContactMM, vecRj .* vecFtMag, [scalNumParticles 1]));

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
            vecSeparating = vecActivePairSource(~boolContact) + scalNumParticles * (vecActivePairDest(~boolContact) - 1);
            if ~isempty(vecSeparating)
                matDispTan(vecSeparating) = 0;
                matDispTanStuck(vecSeparating) = false;
            end
        end

            % ============ Rotational velocity-Verlet half-step ============
            if boolFrictionOn && ~boolThreeD
            vecInertiaC = 0.5 * scalMass .* (vecDiameter / 2) .^ 2;
             % Rotational damping: ref OverDamp.cpp L104-106:
               %   W_n+1 = (T_n - Bt*W_n)/Bt_denorm,  Bt_denorm = 1 + Bt*dt/2
               % This is the over-damped form that makes the rotational DOF
               % unconditionally stable. Without it omega oscillates forever.
            scalBt = scalDissipationAbsolute / sqrt(3);
            if scalGammaTangential > 0
                scalBt = scalGammaTangential;    % user override
            end
            scalBtDenominator = 1 + scalBt * scalTimestep / 2;
            vecTorque = (vecTorque - scalBt .* vecOmega) / scalBtDenominator;
            vecAlpha       = vecTorque ./ vecInertiaC;
            vecOmega       = vecOmega + (vecAlphaPrev + vecAlpha) .* (scalTimestep / 2);
            vecAlphaPrev   = vecAlpha;
            end

% Contact count per particle (coordination number)
        vecCoordNum = accumarray(vecContactNN, 1, [scalNumParticles 1]) ...
                    + accumarray(vecContactMM, 1, [scalNumParticles 1]);

        vecPotentialEnergyHistory(idxStep) = sum(vecPotentialContact) / scalNumParticles;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Drag, boundaries, energy %%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % OverDamped integrator (ref: OverDamp.cpp L97-L100):
          %   F_n+1 = (F_n - Bn*Vel_n) / (1 + Bn*dt/2)
          %   W_n+1 = (W_n - Bt*W_n) / (1 + Bt*dt/2)
          % For frictional contacts the renormalized form prevents the
          % contact-damped oscillation that causes P/P_target to cycle.
          % The frictionless additive form is identical (scalBnDenominator -> 1).
        if boolFrictionOn
            scalBnDenominator = 1 + scalDissipationAbsolute * scalTimestep / 2;
            vecForceX = (vecForceX - scalDissipationAbsolute .* vecVelX) / scalBnDenominator;
            vecForceY = (vecForceY - scalDissipationAbsolute .* vecVelY) / scalBnDenominator;
            if boolThreeD
                vecForceZ = (vecForceZ - scalDissipationAbsolute .* vecVelZ) / scalBnDenominator;
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
            vecKineticEnergyHistory(idxStep) = 0.5 * scalMass * sum(vecVelX.^2 + vecVelY.^2 + vecVelZ.^2) / scalNumParticles;
        else
            vecKineticEnergyHistory(idxStep) = 0.5 * scalMass * sum(vecVelX.^2 + vecVelY.^2) / scalNumParticles;
        end

        % Rotational kinetic energy (friction) must be included in the
        % convergence check: omega carries energy that the translational KE
        % misses; without it scalEk reads ~0 while omega still oscillates.
        if boolFrictionOn && ~boolThreeD
            vecInertiaC = 0.5 * scalMass .* (vecDiameter / 2) .^ 2;   % solid-disk moment of inertia
            vecKineticEnergyHistory(idxStep) = vecKineticEnergyHistory(idxStep) ...
                + 0.5 * sum(vecInertiaC .* vecOmega.^2) / scalNumParticles;
        end

        vecAccelX = vecForceX ./ scalMass;
        vecAccelY = vecForceY ./ scalMass - scalGravity;
        if boolThreeD
            vecAccelZ = vecForceZ ./ scalMass;  % no gravity in Z
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Second step in Verlet integration %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        vecVelX = vecVelX + (vecAccelXPrev + vecAccelX) .* (scalTimestep/2);  % [N x 1]
        vecVelY = vecVelY + (vecAccelYPrev + vecAccelY) .* (scalTimestep/2);
        if boolThreeD
                vecVelZ = vecVelZ + (vecAccelZPrev + vecAccelZ) .* (scalTimestep/2);  % ← correct
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
        scalExcessContacts = scalTotalContacts + scalWallContacts - 2*(scalNumParticles - scalNumRattlers);

        % Mean particle-particle coordination number (full packing, rattlers
        % included). Tracked on every pathway — frictionless and frictional,
        % 2D and 3D — so the converged value is saved with the packing for
        % comparison against the jamming literature (2D frictionless z_iso ~
        % 4, 2D frictional z_iso = 3, 3D frictionless z_iso ~ 6).
        scalMeanCoordNum = mean(vecCoordNum);

        % Pressure estimate from mean potential energy
        scalEp = vecPotentialEnergyHistory(idxStep);
        % Under friction, the tangential spring energy is a constraint DOF,
        % not a compressive load: exclude it from the box-control pressure so
        % the compression target P_target is reached on the NORMAL contacts.
        if boolFrictionOn
           scalEp = scalEp - scalTangentialPE / scalNumParticles;
        end
        if boolHertzian
            scalPressure = (scalEp * (5/2) / scalSpringConstant)^(2/5); % this has an implied d= 1 in the denominator
        else
            scalPressure = sqrt(2 * scalEp / scalSpringConstant); % this has an implied d= 1 in the denominator
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% COMPRESSION DECISIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        scalEk = vecKineticEnergyHistory(idxStep);

        % ===== ENV-GUARDED DIAGNOSTIC (removable; no effect unless GRAN_DIAG set) =====
        if ~isempty(getenv('GRAN_DIAG')) && mod(idxStep, 5000) == 0
            scalDiagMeanCoordNum = mean(vecCoordNum);  % mean particle-particle coordination
            fprintf('DIAG step=%d Ek=%.4e P=%.4e P/Ptrgt=%.4f Lx=%.4f Zn=%.2f\n', idxStep, scalEk, scalPressure, scalPressure/scalPressureTarget, scalBoxWidthX, scalDiagMeanCoordNum);
        end
        % ===== end diagnostic =====

        % Frictional case: friction is on from the first step and the box is
        % driven by the relax-then-resize controller below instead of the
        % frictionless fast/slow phases. Gating on ~boolFrictionOn leaves the
        % frictionless path byte-for-byte identical.
        if boolFrictionOn
              % Frictional box controller + FORCE-BALANCE convergence.
              % Convergence is NOT "box happened to stop moving" or "P in a
              % band" (both are transients of the compression); it is
              % per-particle FORCE BALANCE: max_i |F_net,i| / mean |F_contact|
              % below a small tolerance, sustained for several hundred
              % consecutive in-band steps, with a percolating contact network.
              % This matches OverDamp.cpp's Acc_max < Fthresh and the DEM
              % relaxation stopping criteria in the jamming literature.

              % Force-balance measure: net contact force (normal + tangential)
              % on each grain. A jammed, settled packing has ~zero net force on
              % every grain (mechanical equilibrium under PBC). Normalized by
              % the mean contact force so the threshold is scale-free.
              if scalNumContacts > 0
                  vecNetFx = accumarray(vecContactNN, vecForceContactX, [scalNumParticles 1]) ...
                           - accumarray(vecContactMM, vecForceContactX, [scalNumParticles 1]) ...
                           + accumarray(vecContactNN, vecFtX, [scalNumParticles 1]) ...
                           - accumarray(vecContactMM, vecFtX, [scalNumParticles 1]);
                  vecNetFy = accumarray(vecContactNN, vecForceContactY, [scalNumParticles 1]) ...
                           - accumarray(vecContactMM, vecForceContactY, [scalNumParticles 1]) ...
                           + accumarray(vecContactNN, vecFtY, [scalNumParticles 1]) ...
                           - accumarray(vecContactMM, vecFtY, [scalNumParticles 1]);
                  vecFnetMag     = sqrt(vecNetFx.^2 + vecNetFy.^2);
                  scalMaxFnet    = max(vecFnetMag);
                  scalMeanFc     = mean([abs(vecForceMag); abs(vecFtMag)]);
                  scalForceRatio = scalMaxFnet / max(scalMeanFc, eps);
              else
                  scalForceRatio = inf;   % no contacts: not balanced, not percolating
              end

              % Relax-then-resize box controller. While P < P_target/50 there
              % is no load-bearing network yet, so keep compressing every step
              % (as the frictionless fast phase does). Otherwise resize only
              % once the grains have settled since the last resize — force
              % balance, or kinetic energy small against the elastic energy
              % at P_target — so each decision reads a relaxed P:
              % expand above the band; compress below it, or while the contact
              % network does not percolate; hold in-band. Reversing direction
              % means P_target was overshot, so the strain step is halved
              % (bisection on the box size) down to scalFrictionStrainMin.
              boolRelaxed = (idxStep - scalFrictionLastResize >= scalFrictionRelaxSteps) && ...
                  (scalForceRatio < scalFrictionForceTol || ...
                   scalEk < scalFrictionRelaxKE * scalFrictionEpTarget);
              scalFrictionDir = 0;    % -1 compress, +1 expand, 0 hold
              if scalPressure < scalPressureFastGrow
                  scalFrictionDir = -1;
              elseif boolRelaxed
                  if scalPressure > scalPressureTarget * (1 + scalFrictionDeadBand)
                      scalFrictionDir = 1;
                  elseif scalPressure < scalPressureTarget * (1 - scalFrictionDeadBand) || ...
                          scalMeanCoordNum < scalFrictionZmin
                      scalFrictionDir = -1;
                  end
              end
              boolBoxMoved = (scalFrictionDir ~= 0);
              if boolBoxMoved
                  if scalFrictionDir == -scalFrictionLastDir
                      scalFrictionStrain = max(scalFrictionStrain / 2, scalFrictionStrainMin);
                  end
                  scalStrainFactor = 1 + scalFrictionDir * scalFrictionStrain;
                  scalBoxWidthX  = scalBoxWidthX  * scalStrainFactor;
                  scalBoxHeightY = scalBoxHeightY * scalStrainFactor;
                  vecPosX = vecPosX * scalStrainFactor;
                  vecPosY = vecPosY * scalStrainFactor;
                  boolCellUpdateNeeded = true;
                  scalFrictionLastDir    = scalFrictionDir;
                  scalFrictionLastResize = idxStep;
              end

              % Acceptance: P in-band, contact network percolates, box held, and
              % sustained force balance. The percolation guard (mean Zn) plus the
              % "box held" guard prevent accepting a loose, unjammed state whose
              % net force is trivially small simply because it carries no load.
              boolInBand   = abs(scalPressure - scalPressureTarget) / scalPressureTarget < scalFrictionDeadBand;
              boolPercol   = (scalMeanCoordNum >= scalFrictionZmin);
              boolBalanced = (scalForceRatio < scalFrictionForceTol);
              if boolInBand && boolPercol && boolBalanced && ~boolBoxMoved
                  scalFrictionBalCount = scalFrictionBalCount + 1;
              else
                  scalFrictionBalCount = 0;
              end

              if mod(idxStep, 5000) == 0
                  fprintf('  [fric] step %d | P/Pt=%.3f | Lx=%.4f | maxFnet/meanFc=%.3e | meanCoordNum=%.2f | balCount=%d\n', ...
                      idxStep, scalPressure / scalPressureTarget, scalBoxWidthX, scalForceRatio, scalMeanCoordNum, scalFrictionBalCount);
              end
              if scalFrictionBalCount >= scalFrictionBalWindow
                  fprintf('Frictional convergence (FORCE BALANCE) at step %d | P=%.4e (P/Pt=%.3f) Lx=%.4f maxFnet/meanFc=%.3e meanCoordNum=%.2f\n', ...
                      idxStep, scalPressure, scalPressure / scalPressureTarget, scalBoxWidthX, scalForceRatio, scalMeanCoordNum);
                  break;
              elseif idxStep >= scalFrictionMaxSteps
                  fprintf('Frictional MAX-STEP cap at step %d | P=%.4e (P/Pt=%.3f) maxFnet/meanFc=%.3e meanCoordNum=%.2f\n', ...
                      idxStep, scalPressure, scalPressure / scalPressureTarget, scalForceRatio, scalMeanCoordNum);
                  break;
              end
        elseif boolFastCompressPhase
            % ===== Fast-compress phase (frictionless two-stage) =====
            if scalPressure < scalPressureTarget/50
                scalBoxWidthX= scalBoxWidthX * (1-scalCompressionRateFast);
                scalBoxHeightY = scalBoxHeightY * (1-scalCompressionRateFast);
                vecPosX = vecPosX * (1-scalCompressionRateFast);
                vecPosY = vecPosY * (1-scalCompressionRateFast);
                if boolThreeD
                    scalBoxDepthZ = scalBoxDepthZ * (1-scalCompressionRateFast);
                    vecPosZ = vecPosZ * (1-scalCompressionRateFast);
                end
                boolCellUpdateNeeded = true;
                scalLastCompressStep = idxStep;
            elseif scalPressure < scalPressureTarget && scalEk < 1e-8
                scalBoxWidthX   = scalBoxWidthX * (1-scalCompressionRateFast);
                scalBoxHeightY = scalBoxHeightY * (1-scalCompressionRateFast);
                vecPosX = vecPosX * (1-scalCompressionRateFast);
                vecPosY = vecPosY * (1-scalCompressionRateFast);
                if boolThreeD
                    scalBoxDepthZ = scalBoxDepthZ * (1-scalCompressionRateFast);
                    vecPosZ = vecPosZ * (1-scalCompressionRateFast);
                end
                boolCellUpdateNeeded = true;
                scalLastCompressStep = idxStep;
            elseif scalPressure > scalPressureTarget && scalEk < 1e-10 && idxStep > (scalLastCompressStep+100)
                scalBoxWidthX   = scalBoxWidthX * (1+scalCompressionRateFast);
                scalBoxHeightY = scalBoxHeightY * (1+scalCompressionRateFast);
                vecPosX = vecPosX * (1+scalCompressionRateFast);
                vecPosY = vecPosY * (1+scalCompressionRateFast);
                if boolThreeD
                    scalBoxDepthZ = scalBoxDepthZ * (1+scalCompressionRateFast);
                    vecPosZ = vecPosZ * (1+scalCompressionRateFast);
                end
                boolCellUpdateNeeded = true;
                scalLastCompressStep = idxStep;
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
                scalLastCompressStep = idxStep;
             elseif scalPressure < scalPressureTarget && scalEk < 1e-8
                scalBoxWidthX    = scalBoxWidthX     * (1 - scalCompressionRate);
                scalBoxHeightY    = scalBoxHeightY    * (1 - scalCompressionRate);
                vecPosX = vecPosX * (1 - scalCompressionRate);
                vecPosY = vecPosY * (1 - scalCompressionRate);
                if boolThreeD
                    scalBoxDepthZ = scalBoxDepthZ * (1 - scalCompressionRate);
                    vecPosZ = vecPosZ * (1 - scalCompressionRate);
                end
                boolCellUpdateNeeded = true;
                scalLastCompressStep = idxStep;
             elseif scalPressure > scalPressureTarget
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
                        idxStep, scalPressure, scalPressure / scalPressureTarget, scalBoxWidthX);
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
        sFricState = buildFricState( ...
            vecPosX, vecPosY, vecDiameter, ...
            scalBoxWidthX, scalBoxHeightY, ...
            matDispTan, scalSpringConstant, scalKt, scalMu, ...
            scalGammaNormal, scalGammaTangential, scalMass);
        strFricFilename = [strFilename(1:end-4) '_FricState.mat'];
        sFricStateFile.fricState = sFricState;   % stored in the .mat as 'fricState'
        save(strFricFilename, '-struct', 'sFricStateFile');
        fprintf('fricState saved to: %s\n', strFricFilename);
    end

    fprintf('Loop finished at step %d.\n', idxStep);
    if boolFricGif
        % final (converged) frame, held on screen longer than the others
        [~, boolFricGif] = writeFricGifFrame(hFigGif, strGifFilename, vecGifFrameSize, ...
            vecPosX, vecPosY, vecDiameter, vecTheta, scalBoxWidthX, scalBoxHeightY, ...
            fricGifTitle('Final frictional packing', idxStep, scalPressure, scalPressureTarget, vecDiameter, scalBoxWidthX, scalBoxHeightY, scalMeanCoordNum, scalMu), 2);
        if boolFricGif
            fprintf('Frictional compression GIF saved to: %s\n', strGifFilename);
        end
    end
    scalNumParticlesOriginal = numel(vecPosX);  % particle count before cleanRats (for plot titles)

    %% Snapshot the FULL jammed state (all N particles, rattlers included)
    %    BEFORE cleanRats, for tile-based repetition (see the tiling block at
    %    the end of this file). The per-step box rescalings already wrap
    %    positions into [0, L), so every coordinate is in-box; the 3D tiler
    %    reproduces this state exactly across tile boundaries. Gated on
    %    sOptions.saveFullState so default callers keep byte-identical output.
    if boolSaveFullState
        % The variable names stored in the .mat (K, P_target, N_original,
        % seed, ...) are the file format that packRepeatTile and other
        % loaders read, so they are kept as the field names of the struct
        % that save -struct writes out (one variable per field; the field
        % list is passed explicitly so Octave keeps this order too).
        sFullStateFile = struct();
        sFullStateFile.vecPosX = vecPosX;
        sFullStateFile.vecPosY = vecPosY;
        if boolThreeD
            sFullStateFile.vecPosZ = vecPosZ;
        end
        sFullStateFile.vecDiameter    = vecDiameter;
        sFullStateFile.scalBoxWidthX  = scalBoxWidthX;
        sFullStateFile.scalBoxHeightY = scalBoxHeightY;
        if boolThreeD
            sFullStateFile.scalBoxDepthZ = scalBoxDepthZ;
        end
        sFullStateFile.K                = scalSpringConstant;
        sFullStateFile.P_target         = scalPressureTarget;
        sFullStateFile.scalPressure     = scalPressure;
        sFullStateFile.N_original       = scalNumParticlesOriginal;
        sFullStateFile.seed             = scalSeed;
        sFullStateFile.scalRoundedWidth = scalRoundedWidth;
        cellFullStateFields = fieldnames(sFullStateFile);
        save(strFullFilename, '-struct', 'sFullStateFile', cellFullStateFields{:});
        % Keep an in-memory copy for the 3D repeat-tile block at the end of
        % this file: after cleanRats the local position/diameter variables are
        % rebound to the BACKBONE, so the full (rattler-inclusive) state must
        % be captured here, before cleanRats runs.
        vecTileSrcX = vecPosX;
        vecTileSrcY = vecPosY;
        vecTileSrcD = vecDiameter;
        scalNumParticlesTileSrc = scalNumParticles;   % full per-tile particle count (rattlers included)
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
        [matSphereX, matSphereY, matSphereZ] = sphere(16);

        % Main particles (full pre-cleanRats packing)
        for idxParticle = 1:scalNumParticles
            scalRadius = vecDiameter(idxParticle)/2;
            surf(scalRadius*matSphereX + vecPosX(idxParticle), ...
                 scalRadius*matSphereY + vecPosY(idxParticle), ...
                 scalRadius*matSphereZ + vecPosZ(idxParticle), ...
                'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.6);
        end

        % Ghost particles on +x, +y, +z faces
        matOffsets3D = [scalBoxWidthX, 0, 0; ...
                        0, scalBoxHeightY, 0; ...
                        0, 0, scalBoxDepthZ];  % [3 x 3] one offset per face

        for idxFace = 1:3
            scalOffsetX = matOffsets3D(idxFace, 1);
            scalOffsetY = matOffsets3D(idxFace, 2);
            scalOffsetZ = matOffsets3D(idxFace, 3);
            for idxParticle = 1:scalNumParticles
                scalRadius = vecDiameter(idxParticle)/2;
                surf(scalRadius*matSphereX + vecPosX(idxParticle) + scalOffsetX, ...
                     scalRadius*matSphereY + vecPosY(idxParticle) + scalOffsetY, ...
                     scalRadius*matSphereZ + vecPosZ(idxParticle) + scalOffsetZ, ...
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
        title(sprintf('Before cleanRats: N=%d, phi=%.4f', scalNumParticles, scalPackingFractionFull));

    else
        % Main particles (full pre-cleanRats packing)
        for idxParticle = 1:scalNumParticles
            rectangle('Position', [vecPosX(idxParticle) - vecDiameter(idxParticle)/2, ...
                                    vecPosY(idxParticle) - vecDiameter(idxParticle)/2, ...
                                    vecDiameter(idxParticle), vecDiameter(idxParticle)], ...
                'Curvature', [1 1], 'FaceColor', 'b', 'EdgeColor', 'none');
        end

        % Ghost particles: all 8 surrounding tiles
        matOffsets2D = [scalBoxWidthX,  0; ...
                       -scalBoxWidthX,  0; ...
                        0,  scalBoxHeightY; ...
                        0, -scalBoxHeightY; ...
                        scalBoxWidthX,  scalBoxHeightY; ...
                       -scalBoxWidthX,  scalBoxHeightY; ...
                        scalBoxWidthX, -scalBoxHeightY; ...
                       -scalBoxWidthX, -scalBoxHeightY];

        for idxFace = 1:8
            scalOffsetX = matOffsets2D(idxFace, 1);
            scalOffsetY = matOffsets2D(idxFace, 2);
            for idxParticle = 1:scalNumParticles
                rectangle('Position', [vecPosX(idxParticle) + scalOffsetX - vecDiameter(idxParticle)/2, ...
                                        vecPosY(idxParticle) + scalOffsetY - vecDiameter(idxParticle)/2, ...
                                        vecDiameter(idxParticle), vecDiameter(idxParticle)], ...
                    'Curvature', [1 1], 'FaceColor', 'r', 'EdgeColor', 'none', ...
                    'FaceAlpha', 0.15);
            end
        end

        axis equal;
        axis([-scalBoxWidthX 2*scalBoxWidthX -scalBoxHeightY 2*scalBoxHeightY]);
        title(sprintf('Before cleanRats: N=%d, phi=%.4f', scalNumParticles, scalPackingFractionFull));
    end
    drawnow;
    hold off;

    % Export the before-cleanRats packing photo. strFilename already includes
    % strSavePath, so just swap the extension; the export mirrors the .mat's
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
    catch sBeforePngME
        warning('pack:PNGExportFailed', 'Could not export before-cleanRats plot: %s', sBeforePngME.message);
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
         scalMuCleanRats = scalMu * boolFrictionOn;
         [matPositions, vecRadii] = cleanRats(matPositions, vecRadii, scalBoxHeightY, scalBoxWidthX, [], false, scalMuCleanRats);
        vecPosX = matPositions(:,1);
        vecPosY = matPositions(:,2);
    end
    vecDiameter = vecRadii .* 2;
    scalNumParticlesClean = size(matPositions, 1);
    fprintf('cleanRats complete. %d particles remaining (of %d original).\n', scalNumParticlesClean, scalNumParticles);
    if scalNumParticlesClean == 0
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
    if boolHertzian

        % update data because cleanRats may have removed particles
        vecRadiiFinal       = vecDiameter ./ 2;
        matContactDistFinal = (vecDiameter + vecDiameter') / 2;
        scalNumFinal = numel(vecPosX);

        vecHertzNN = zeros(scalNumFinal * 12, 1);
        vecHertzMM = zeros(scalNumFinal * 12, 1);
        vecHertzKeff  = zeros(scalNumFinal * 12, 1);
        scalNumHertzContacts = 0;

        for idxNN = 1:scalNumFinal
            for idxMM = idxNN+1:scalNumFinal
                scalSepX = vecPosX(idxMM) - vecPosX(idxNN);
                scalSepY = vecPosY(idxMM) - vecPosY(idxNN);
                scalSepX = scalSepX - scalBoxWidthX * round(scalSepX / scalBoxWidthX);
                scalSepY = scalSepY - scalBoxHeightY * round(scalSepY / scalBoxHeightY);
                if boolThreeD
                    scalSepZ = vecPosZ(idxMM) - vecPosZ(idxNN);
                    scalSepZ = scalSepZ - scalBoxDepthZ * round(scalSepZ / scalBoxDepthZ);
                    scalDist = sqrt(scalSepX^2 + scalSepY^2 + scalSepZ^2);
                else
                    scalDist = sqrt(scalSepX^2 + scalSepY^2);
                end

                scalSumRadii = matContactDistFinal(idxNN, idxMM); % grab the minimum distanced needed for contact
                scalDelta = scalSumRadii - scalDist; % if negative, no contact

                % Go through each contact and assign k = dF/d(delta) for F_hertzian
                if scalDelta > 0
                    scalReff = (vecRadiiFinal(idxNN) * vecRadiiFinal(idxMM)) / scalSumRadii; % effective radius from curvature
                    scalKeff = 2 * scalSpringConstant * sqrt(scalReff * scalDelta); % k_eff = dF/d(delta)
                    scalNumHertzContacts = scalNumHertzContacts + 1; % this is indexes, so +1 because matlab isbase 1 

                    % assign row (scalNumHertzContacts) a particle, the particle it's in contact with, and an k_eff
                    vecHertzNN(scalNumHertzContacts) = idxNN; 
                    vecHertzMM(scalNumHertzContacts) = idxMM;
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
        [matSphereX, matSphereY, matSphereZ] = sphere(16);

        % Main particles
        for idxParticle = 1:scalNumParticlesClean
            scalRadius = vecDiameter(idxParticle)/2;
            surf(scalRadius*matSphereX + vecPosX(idxParticle), ...
                 scalRadius*matSphereY + vecPosY(idxParticle), ...
                 scalRadius*matSphereZ + vecPosZ(idxParticle), ...
                'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.6);
        end

        % Ghost particles on +x, +y, +z faces
        matOffsets3D = [scalBoxWidthX, 0, 0; ...
                        0, scalBoxHeightY, 0; ...
                        0, 0, scalBoxDepthZ];  % [3 x 3] one offset per face

        for idxFace = 1:3
            scalOffsetX = matOffsets3D(idxFace, 1);
            scalOffsetY = matOffsets3D(idxFace, 2);
            scalOffsetZ = matOffsets3D(idxFace, 3);
            for idxParticle = 1:scalNumParticlesClean
                scalRadius = vecDiameter(idxParticle)/2;
                surf(scalRadius*matSphereX + vecPosX(idxParticle) + scalOffsetX, ...
                     scalRadius*matSphereY + vecPosY(idxParticle) + scalOffsetY, ...
                     scalRadius*matSphereZ + vecPosZ(idxParticle) + scalOffsetZ, ...
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
        title(sprintf('After cleanRats: N=%d (of %d original), mean coord num=%.2f', scalNumParticlesClean, scalNumParticlesOriginal, scalMeanCoordNum));

    else
        % Main particles
        for idxParticle = 1:scalNumParticlesClean
            rectangle('Position', [vecPosX(idxParticle) - vecDiameter(idxParticle)/2, ...
                                    vecPosY(idxParticle) - vecDiameter(idxParticle)/2, ...
                                    vecDiameter(idxParticle), vecDiameter(idxParticle)], ...
                'Curvature', [1 1], 'FaceColor', 'b', 'EdgeColor', 'none');
        end

        % Ghost particles: all 8 surrounding tiles
        matOffsets2D = [scalBoxWidthX,  0; ...
                       -scalBoxWidthX,  0; ...
                        0,  scalBoxHeightY; ...
                        0, -scalBoxHeightY; ...
                        scalBoxWidthX,  scalBoxHeightY; ...
                       -scalBoxWidthX,  scalBoxHeightY; ...
                        scalBoxWidthX, -scalBoxHeightY; ...
                       -scalBoxWidthX, -scalBoxHeightY];

        for idxFace = 1:8
            scalOffsetX = matOffsets2D(idxFace, 1);
            scalOffsetY = matOffsets2D(idxFace, 2);
            for idxParticle = 1:scalNumParticlesClean
                rectangle('Position', [vecPosX(idxParticle) + scalOffsetX - vecDiameter(idxParticle)/2, ...
                                        vecPosY(idxParticle) + scalOffsetY - vecDiameter(idxParticle)/2, ...
                                        vecDiameter(idxParticle), vecDiameter(idxParticle)], ...
                    'Curvature', [1 1], 'FaceColor', 'r', 'EdgeColor', 'none', ...
                    'FaceAlpha', 0.15);
            end
        end

        axis equal;
        axis([-scalBoxWidthX 2*scalBoxWidthX -scalBoxHeightY 2*scalBoxHeightY]);
        title(sprintf('After cleanRats: N=%d (of %d original), mean coord num=%.2f', scalNumParticlesClean, scalNumParticlesOriginal, scalMeanCoordNum));
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
    catch sAfterPngME
        warning('pack:PNGExportFailed', 'Could not export after-cleanRats plot: %s', sAfterPngME.message);
    end

%% Save results

    % if boolThreeD
    %     scalRoundedWidth = round(scalNumParticles^(1/3));
    %     strFilename = sprintf('%s3D_N%d_P%s_Width%d_Seed%d.mat', ...
    %         strSavePath, scalNumParticles, num2str(scalPressureTarget), scalRoundedWidth, scalSeed);
    % else
    %     scalRoundedWidth = round(sqrt(scalNumParticles));
    %     strFilename = sprintf('%s2D_N%d_P%s_Width%d_Seed%d.mat', ...
    %         strSavePath, scalNumParticles, num2str(scalPressureTarget), scalRoundedWidth, scalSeed);
    % end

    scalNumParticlesOriginal = scalNumParticles;
    scalNumParticlesClean = size(matPositions, 1);

    % Compute packing fraction after cleanRats
    if boolThreeD
        scalPackingFraction = computePackingFraction(vecDiameter, scalBoxWidthX, scalBoxHeightY, scalBoxDepthZ);
    else
        scalPackingFraction = computePackingFraction(vecDiameter, scalBoxWidthX, scalBoxHeightY);
    end
    fprintf('Packing fraction after cleanRats: %.4f\n', scalPackingFraction);

    if boolCalcEig
        % 3D hessian not yet implemented — skip eigenmodes
        if boolThreeD
            warning('calc_eig not supported for 3D yet — saving positions only.');
        else
            matPositions = [vecPosX, vecPosY];
            vecRadii = vecDiameter ./ 2;
            [matPositions, vecRadii] = cleanRats(matPositions, vecRadii, scalSpringConstant, scalBoxHeightY, scalBoxWidthX);
            matHessian = hess2d(matPositions, vecRadii, scalSpringConstant, scalBoxHeightY, scalBoxWidthX);
            [matEigenVectors, matEigenValues] = eig(matHessian);
        end
    end

    % Backbone .mat: the stored variable names (K, P_target, N, N_original,
    % ...) are the file format that simMD, packRepeatTile and the tests
    % load, so they are kept as the field names of the struct that
    % save -struct writes out (one variable per field; the field list is
    % passed explicitly so Octave keeps this order too).
    sPackingFile = struct();
    sPackingFile.vecPosX = vecPosX;
    sPackingFile.vecPosY = vecPosY;
    if boolThreeD
        sPackingFile.vecPosZ = vecPosZ;
    end
    sPackingFile.vecDiameter    = vecDiameter;
    sPackingFile.scalBoxWidthX  = scalBoxWidthX;
    sPackingFile.scalBoxHeightY = scalBoxHeightY;
    if boolThreeD
        sPackingFile.scalBoxDepthZ = scalBoxDepthZ;
    end
    sPackingFile.K                       = scalSpringConstant;
    sPackingFile.P_target                = scalPressureTarget;
    sPackingFile.scalPressure            = scalPressure;
    sPackingFile.N                       = scalNumParticlesClean;      % backbone count
    sPackingFile.N_original              = scalNumParticlesOriginal;   % count before cleanRats
    sPackingFile.scalPackingFraction     = scalPackingFraction;
    sPackingFile.scalPackingFractionFull = scalPackingFractionFull;
    sPackingFile.scalMeanCoordNum        = scalMeanCoordNum;
    if boolCalcEig && ~boolThreeD
        sPackingFile.matEigenVectors = matEigenVectors;
        sPackingFile.matEigenValues  = matEigenValues;
    end
    if boolHertzian
        sPackingFile.vecHertzNN   = vecHertzNN;
        sPackingFile.vecHertzMM   = vecHertzMM;
        sPackingFile.vecHertzKeff = vecHertzKeff;
    end
    sPackingFile.boolFrictionOn = boolFrictionOn;
    sPackingFile.scalMu         = scalMu;
    sPackingFile.scalKt         = scalKt;
    cellPackingFields = fieldnames(sPackingFile);
    save(strFilename, '-struct', 'sPackingFile', cellPackingFields{:});

    fprintf('File saved to: %s\n', strFilename);

    %% Repeat-tile the packing in x, y, and/or z for superlattice packings.
    %   The tile SOURCE is the FULL jammed state captured before cleanRats
    %   (rattlers included), NOT the backbone saved to strFilename: the
    %   user-visible superlattice repeats the exact simulated state, and
    %   rattler positions must line up across the periodic boundary.
    if boolThreeD
        if scalXMult ~= 1 || scalYMult ~= 1 || scalZMult ~= 1
            if ~boolSaveFullState
                error('pack:NeedFullStateFor3DTile', ...
                    ['3D repeat-tile requires the full (pre-cleanRats) state. ', ...
                     'Call pack(..., options) with options.saveFullState = true.']);
            end
            [vecPosXFinal, vecPosYFinal, vecPosZFinal, vecDiameterFinal, ...
                scalBoxWidthXTiled, scalBoxHeightYFinal, scalBoxDepthZFinal, ...
                scalNTiled] = tile3D( ...
                vecTileSrcX, vecTileSrcY, vecTileSrcZ, vecTileSrcD, ...
                scalBoxWidthX, scalBoxHeightY, scalBoxDepthZ, ...
                scalXMult, scalYMult, scalZMult, scalNumParticlesTileSrc);

            % Metrics on the tiled packing (rattlers included, matching the
            % stored particles). Packing fraction is identical to the base
            % tile by construction; coordination number under full PBC equals
            % the base tile's full-state Zn (tiling replicates the contact
            % network).
            scalPackingFractionTiled  = computePackingFraction(vecDiameterFinal, scalBoxWidthXTiled, scalBoxHeightYFinal, scalBoxDepthZFinal);
            scalPackingFractionFullTiled = scalPackingFractionTiled;
            % NOTE: the direct recompute on the tiled packing is ONLY for small
            % packings used by the unit tests (tests/testPack/testPack.m).
            % computeMeanCoordNum is O(N^2) in memory (N x N ndgrid), so for
            % production-sized tilings (N > 1200; e.g. 1000 x 160 = 160000
            % particles needs ~400 GB) it is skipped and the base tile's
            % full-state Zn is used instead — exact under full PBC.
            scalMaxNumForDirectCoordNum = 1200;
            if scalNTiled <= scalMaxNumForDirectCoordNum
                scalMeanCoordNumTiled = computeMeanCoordNum(vecPosXFinal, vecPosYFinal, vecDiameterFinal, ...
                    scalBoxWidthXTiled, scalBoxHeightYFinal, vecPosZFinal, scalBoxDepthZFinal);
            else
                scalMeanCoordNumTiled = computeMeanCoordNum(vecTileSrcX, vecTileSrcY, vecTileSrcD, ...
                    scalBoxWidthX, scalBoxHeightY, vecTileSrcZ, scalBoxDepthZ);
            end
            fprintf('3D tiled packing: N=%d (of %d per tile), PF=%.4f, mean coordination number=%.4f\n', ...
                scalNTiled, scalNumParticlesTileSrc, scalPackingFractionTiled, scalMeanCoordNumTiled);

            % Tiled output filename: the normal packing name,
            % 3D_N%d_P%s_Width%d_Seed%d[_Hertz].mat, with N = the TILED particle
            % count, so simMD finds it like any other packing. Width and Seed
            % are the base tile's; the multipliers are stored inside the file
            % (x_mult, y_mult, z_mult). N differs from the base tile's, so this
            % never overwrites the backbone file.
            if boolHertzian
                strTiledFilename = sprintf('%s3D_N%d_P%s_Width%d_Seed%d_Hertz.mat', ...
                    strSavePath, scalNTiled, num2str(scalPressureTarget), scalRoundedWidth, scalSeed);
            else
                strTiledFilename = sprintf('%s3D_N%d_P%s_Width%d_Seed%d.mat', ...
                    strSavePath, scalNTiled, num2str(scalPressureTarget), scalRoundedWidth, scalSeed);
            end

            % Save the TILED (superlattice) packing — the full pre-cleanRats
            % state repeated across the tile grid — not the backbone. The
            % stored variable names match the backbone .mat convention
            % (vecPosX, N, scalPackingFraction, ...) so downstream loaders are
            % unchanged; save -struct writes one variable per field, in the
            % order of the explicit field list.
            sTiledFile = struct();
            sTiledFile.vecPosX                 = vecPosXFinal;
            sTiledFile.vecPosY                 = vecPosYFinal;
            sTiledFile.vecPosZ                 = vecPosZFinal;
            sTiledFile.vecDiameter             = vecDiameterFinal;
            sTiledFile.scalBoxWidthX           = scalBoxWidthXTiled;
            sTiledFile.scalBoxHeightY          = scalBoxHeightYFinal;
            sTiledFile.scalBoxDepthZ           = scalBoxDepthZFinal;
            sTiledFile.K                       = scalSpringConstant;
            sTiledFile.P_target                = scalPressureTarget;
            sTiledFile.scalPressure            = scalPressure;
            sTiledFile.N                       = scalNTiled;
            sTiledFile.N_original              = scalNTiled;   % full (pre-cleanRats) tiled state
            sTiledFile.scalPackingFraction     = scalPackingFractionTiled;
            sTiledFile.scalPackingFractionFull = scalPackingFractionFullTiled;
            sTiledFile.scalMeanCoordNum        = scalMeanCoordNumTiled;
            sTiledFile.boolFrictionOn          = boolFrictionOn;
            sTiledFile.scalMu                  = scalMu;
            sTiledFile.scalKt                  = scalKt;
            sTiledFile.x_mult                  = scalXMult;
            sTiledFile.y_mult                  = scalYMult;
            sTiledFile.z_mult                  = scalZMult;
            sTiledFile.vecPosXFinal            = vecPosXFinal;
            sTiledFile.vecPosYFinal            = vecPosYFinal;
            sTiledFile.vecPosZFinal            = vecPosZFinal;
            sTiledFile.vecDiameterFinal        = vecDiameterFinal;
            sTiledFile.scalBoxWidthXTiled      = scalBoxWidthXTiled;
            sTiledFile.scalBoxHeightYFinal     = scalBoxHeightYFinal;
            sTiledFile.scalBoxDepthZFinal      = scalBoxDepthZFinal;
            sTiledFile.scalNTiled              = scalNTiled;
            sTiledFile.NTileSrc                = scalNumParticlesTileSrc;
            sTiledFile.scalPackingFractionTiled     = scalPackingFractionTiled;
            sTiledFile.scalPackingFractionFullTiled = scalPackingFractionFullTiled;
            sTiledFile.scalMeanCoordNumTiled        = scalMeanCoordNumTiled;
            cellTiledFields = fieldnames(sTiledFile);
            save(strTiledFilename, '-struct', 'sTiledFile', cellTiledFields{:});
            fprintf('3D tiled packing saved to: %s\n', strTiledFilename);
        end
    else
        % 2D path: unchanged — reuse the existing packRepeatTile on the
        % backbone (its documented behavior).
        if scalXMult ~= 1 || scalYMult ~= 1
            packRepeatTile(scalNumParticlesOriginal, scalSpringConstant, scalPressureTarget, scalRoundedWidth, scalSeed, scalXMult, scalYMult, ...
                boolCalcEig, strSavePath, strSavePath);
            disp("Tile saved to: " + strFilename);
        end
    end

    toc
end

function [vecPosX, vecPosY, cellParticleList, scalNumCellsX, scalNumCellsY] = ...
        rebuildCellList(vecPosX, vecPosY, scalBoxWidthX, scalBoxHeightY, scalRawCellWidth, scalTimestep, scalNumParticles)

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
    cellParticleList = accumarray(vecCellLinearIdx, (1:scalNumParticles)', [scalNumCellsX*scalNumCellsY 1], @(vecCellMembers){vecCellMembers});
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

function sFricState = buildFricState( ...
        vecPosX, vecPosY, vecDiameter, ...
        scalBoxWidthX, scalBoxHeightY, ...
        matDispTan, scalSpringConstant, scalKt, scalMu, ...
        scalGammaNormal, scalGammaTangential, scalMass)
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

    scalNumParticles = size(vecPosX, 1);
    vecRadius        = vecDiameter / 2;
    vecInertia       = 0.5 * scalMass .* (vecRadius .^ 2);

    sFricState = struct();
    sFricState.vecContactNN = zeros(0,1);
    sFricState.vecContactMM = zeros(0,1);
    sFricState.vecOverlap     = zeros(0,1);
    sFricState.vecFn          = zeros(0,1);
    sFricState.vecFt          = zeros(0,1);
    sFricState.vecUnitTanX    = zeros(0,1);
    sFricState.vecUnitTanY    = zeros(0,1);

    for idxNN = 1:(scalNumParticles-1)
        for idxMM = idxNN+1:scalNumParticles
            scalSepX = vecPosX(idxMM) - vecPosX(idxNN);
            scalSepY = vecPosY(idxMM) - vecPosY(idxNN);
            scalSepX = scalSepX - scalBoxWidthX * round(scalSepX / scalBoxWidthX);
            scalSepY = scalSepY - scalBoxHeightY * round(scalSepY / scalBoxHeightY);
            scalSepDist = sqrt(scalSepX^2 + scalSepY^2);
            if scalSepDist < 1e-12, continue; end
            % Contact distance is the SUM OF RADII = (d_i + d_j)/2, NOT the sum
            % of diameters. (vecDiameter(idxNN)+vecDiameter(idxMM)) would double-count
            % and report ~2x overlaps, so buildFricState's forces/overlaps are
            % only correct if this is (r_i + r_j).
            scalContactDist = (vecDiameter(idxNN) + vecDiameter(idxMM)) / 2;   % = r_i + r_j
            scalOverlap     = scalContactDist - scalSepDist;
            if scalOverlap <= 0, continue; end
            scalNormalX  = scalSepX / scalSepDist;
            scalNormalY  = scalSepY / scalSepDist;
            scalUnitTanX =  scalNormalY;
            scalUnitTanY = -scalNormalX;
            scalFn       = -scalSpringConstant * scalOverlap;
            scalDispTan  = matDispTan(idxNN + scalNumParticles * (idxMM - 1));   % [1 x 1] pair tangential displacement
            scalFt       = -scalKt * scalDispTan;
            sFricState.vecContactNN = [sFricState.vecContactNN; idxNN];
            sFricState.vecContactMM = [sFricState.vecContactMM; idxMM];
            sFricState.vecOverlap     = [sFricState.vecOverlap;     scalOverlap];
            sFricState.vecFn          = [sFricState.vecFn;          scalFn];
            sFricState.vecFt          = [sFricState.vecFt;          scalFt];
            sFricState.vecUnitTanX    = [sFricState.vecUnitTanX;    scalUnitTanX];
            sFricState.vecUnitTanY    = [sFricState.vecUnitTanY;    scalUnitTanY];
        end
    end

    sFricState.K           = scalSpringConstant;
    sFricState.Kt          = scalKt;
    sFricState.mu          = scalMu;
    sFricState.gammaNormal = scalGammaNormal;
    sFricState.gammaTang   = scalGammaTangential;
    sFricState.M           = scalMass;
    sFricState.vecRadius = vecRadius;
    sFricState.vecInertia = vecInertia;
    sFricState.N     = scalNumParticles;
    sFricState.boxLx = scalBoxWidthX;
    sFricState.boxLy = scalBoxHeightY;
    sFricState.matDispTan = matDispTan;
end

function [vecFrameSize, boolOk] = writeFricGifFrame(hFig, strGifFilename, vecFrameSize, ...
        vecPosX, vecPosY, vecDiameter, vecTheta, scalBoxWidthX, scalBoxHeightY, cellTitle, scalDelay)
% writeFricGifFrame -- Draw the 2D frictional packing and append it to an
% animated GIF. Each disk gets a line across its diameter at its rotation
% angle vecTheta, so grain rotation is visible from frame to frame; disks
% crossing the periodic boundary are drawn with their images. The first
% call (empty vecFrameSize) creates the file, later calls append frames
% cropped/padded to the first frame's size. Frames are quantized to a fixed
% 6x6x6 RGB palette (works in MATLAB and Octave). Any failure (no display,
% figure closed, write error) warns and returns boolOk = false so the caller
% stops recording; it never aborts the packing run.
    boolOk = true;
    boolFirstFrame = isempty(vecFrameSize);
    try
        vecRadius = vecDiameter / 2;
        boolLarge = vecDiameter > min(vecDiameter) * (1 + 1e-9);

        % Periodic images: keep every copy whose disk overlaps the box
        vecX = []; vecY = []; vecR = []; vecT = []; vecL = false(0, 1);
        for scalOffsetX = [-1 0 1] * scalBoxWidthX
            for scalOffsetY = [-1 0 1] * scalBoxHeightY
                boolIn = (vecPosX + scalOffsetX + vecRadius > 0) & (vecPosX + scalOffsetX - vecRadius < scalBoxWidthX) & ...
                         (vecPosY + scalOffsetY + vecRadius > 0) & (vecPosY + scalOffsetY - vecRadius < scalBoxHeightY);
                vecX = [vecX; vecPosX(boolIn) + scalOffsetX];   %#ok<AGROW>
                vecY = [vecY; vecPosY(boolIn) + scalOffsetY];   %#ok<AGROW>
                vecR = [vecR; vecRadius(boolIn)];      %#ok<AGROW>
                vecT = [vecT; vecTheta(boolIn)];       %#ok<AGROW>
                vecL = [vecL; boolLarge(boolIn)];      %#ok<AGROW>
            end
        end

        clf(hFig);
        hAx = axes('Parent', hFig);
        hold(hAx, 'on');
        vecCircle = linspace(0, 2*pi, 33)';
        vecCircle(end) = [];
        % disks (one patch per size class), then one NaN-separated line object
        % for all the diameters
        if any(~vecL)
            patch(hAx, cos(vecCircle) * vecR(~vecL)' + vecX(~vecL)', ...
                       sin(vecCircle) * vecR(~vecL)' + vecY(~vecL)', ...
                  [0.6 0.8 1.0], 'EdgeColor', 'k');
        end
        if any(vecL)
            patch(hAx, cos(vecCircle) * vecR(vecL)' + vecX(vecL)', ...
                       sin(vecCircle) * vecR(vecL)' + vecY(vecL)', ...
                  [1.0 0.8 0.6], 'EdgeColor', 'k');
        end
        vecDx = vecR .* cos(vecT);
        vecDy = vecR .* sin(vecT);
        matLineX = [vecX - vecDx, vecX + vecDx, nan(size(vecX))]';
        matLineY = [vecY - vecDy, vecY + vecDy, nan(size(vecY))]';
        line(hAx, matLineX(:), matLineY(:), 'Color', [0.8 0 0], 'LineWidth', 1.5);
        axis(hAx, 'equal');
        axis(hAx, [0 scalBoxWidthX 0 scalBoxHeightY]);
        box(hAx, 'on');
        hTitle = title(hAx, cellTitle, 'Interpreter', 'latex');
        set(hTitle, 'Color', [0 0 0], 'FontWeight', 'bold', 'FontSize', 14);
        drawnow;
        drawnow;

        imgFrame = getframe(hFig);
        imgFrame = imgFrame.cdata;
        if boolFirstFrame
            vecFrameSize = [size(imgFrame, 1) size(imgFrame, 2)];
        end
        % the figure may be resized mid-run: crop/pad onto a white canvas
        imgCanvas = 255 * ones([vecFrameSize 3], 'uint8');
        scalRows = min(vecFrameSize(1), size(imgFrame, 1));
        scalCols = min(vecFrameSize(2), size(imgFrame, 2));
        imgCanvas(1:scalRows, 1:scalCols, :) = imgFrame(1:scalRows, 1:scalCols, :);

        [imgIndexed, matColorMap] = rgb2ind(imgCanvas, 256);
        if boolFirstFrame
            imwrite(imgIndexed, matColorMap, strGifFilename, 'gif', ...
                'LoopCount', Inf, 'DelayTime', scalDelay);
            fprintf('Frictional compression GIF path: %s\n', strGifFilename);
        else
            imwrite(imgIndexed, matColorMap, strGifFilename, 'gif', ...
                'WriteMode', 'append', 'DelayTime', scalDelay);
        end
    catch sGifFrameME
        warning('pack:GIFExportFailed', ...
            'Could not write frictional compression GIF frame (recording stopped): %s', sGifFrameME.message);
        boolOk = false;
    end
end

function cellTitle = fricGifTitle(strLabel, idxStep, scalPressure, scalPressureTarget, vecDiameter, ...
        scalBoxWidthX, scalBoxHeightY, scalMeanCoordNum, scalMu)

    scalPhi = sum(pi * vecDiameter.^2 / 4) / (scalBoxWidthX * scalBoxHeightY);

    cellTitle = { ...
        sprintf('%s ($\\mu = %.2f$), step %d', ...
                strLabel, scalMu, idxStep), ...
        sprintf('$P/P_{\\mathrm{target}} = %.3f \\quad \\phi = %.4f \\quad Z = %.2f$', ...
                scalPressure / scalPressureTarget, scalPhi, scalMeanCoordNum) ...
    };
end

