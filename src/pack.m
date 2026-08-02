function pack(N, K, D, G, M, P_target, seed, plotit, x_mult, y_mult, z_mult, calc_eig, save_path, options)

    %Function to create 2D packing with following input parameters:
    % N,        Number of particles
    % K,        Spring constant (Hookean) or Hertz prefactor
    % D,        Average diameter
    % G,        Ratio of large to small particle diameter (typically 1.4)
    % M,        Mass of particles
    % P_target, Targeted threshold pressure
    % seed,     RNG seed
    % plotit,   Boolean: plot during simulation
    % x_mult,   Tile repeat factor in x
    % y_mult,   Tile repeat factor in y
    % calc_eig, Boolean: compute and save eigenmodes
    % save_path, Output directory
    % options.hertzian, Boolean: use Hertzian (nonlinear) contact law
    %                Force = -4/3 K sqrt(R_eff) delta^(3/2)
    %                Potential = 4/3 * 2/5 K sqrt(R_eff) delta^(5/2)
    %                
    %
    % Example:
    % pack(400,100,1,1.4,1,.0001,1,false,1,1,false,'in/tiles/')

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
        options.hertzian (1,1) logical = false
    end

    % check to see if 3d path is needed
    boolThreeD = (z_mult ~= 0);

    % Check if packing already exists — skip if so
    if boolThreeD
        scalRoundedWidth = round(N^(1/3));
        strFilename = sprintf('%s3D_N%d_P%s_Width%d_Seed%d.mat', ...
            save_path, N, num2str(P_target), scalRoundedWidth, seed);
    else
        scalRoundedWidth = round(sqrt(N));
        strFilename = sprintf('%s2D_N%d_P%s_Width%d_Seed%d.mat', ...
            save_path, N, num2str(P_target), scalRoundedWidth, seed);
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
    scalCompressionRateFast = 0.01;

    boolConverged = false; % only update plot after each compression step
    boolCellUpdateNeeded = true; % make sure to update cell list on first step
    boolFastCompressPhase = true;

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
    scalMaxSteps = 1e8; % enough to ensure convergence

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

    vecKineticEnergyHistory   = zeros(scalMaxSteps, 1);  % [scalMaxSteps x 1]
    vecPotentialEnergyHistory = zeros(scalMaxSteps, 1);  % [scalMaxSteps x 1]

%% Verlet cell list setup
    % Determine cell size rounded to be at least 1*G*D
    % to avoid missing interactions, I tried this out 
    % many times and 3 works best
    scalRawCellWidth = 3 * G * D;

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
    % vecPairNN = [1, 1, 2, 3, ...]   <- first particle of each pair
    % vecPairMM = [2, 3, 3, 4, ...]   <- second particle of each pair
    vecPairNN = zeros(N*12, 1);  % [N*12 x 1]
    vecPairMM = zeros(N*12, 1);  % [N*12 x 1]
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
            [vecPairNN, vecPairMM, scalNumPairs, scalMaxPairs] = findNeighbors3D( ...
                cellParticleList, scalNumCellsX, scalNumCellsY, scalNumCellsZ, ...
                vecPairNN, vecPairMM, scalMaxPairs);
        else
            [vecPairNN, vecPairMM, scalNumPairs, scalMaxPairs] = findNeighbors2D( ...
                cellParticleList, scalNumCellsX, scalNumCellsY, ...
                vecPairNN, vecPairMM, scalMaxPairs);
        end

        vecActivePairNN = vecPairNN(1:scalNumPairs);  % [scalNumPairs x 1]
        vecActivePairMM = vecPairMM(1:scalNumPairs);  % [scalNumPairs x 1]

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Vectorized force evaluation %%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % All arithmetic operates on [scalNumPairs x 1] vectors — no inner loops.
        % no mod() here because we need signed distances
        vecSepX = vecPosX(vecActivePairMM) - vecPosX(vecActivePairNN);
        vecSepX = vecSepX - scalBoxWidthX  * round(vecSepX / scalBoxWidthX);
        vecSepY = vecPosY(vecActivePairMM) - vecPosY(vecActivePairNN);
        vecSepY = vecSepY - scalBoxHeightY * round(vecSepY / scalBoxHeightY);
        if boolThreeD
            vecSepZ = vecPosZ(vecActivePairMM) - vecPosZ(vecActivePairNN);
            vecSepZ = vecSepZ - scalBoxDepthZ * round(vecSepZ / scalBoxDepthZ);
        end

        vecContactDist = matContactDist(vecActivePairNN + N*(vecActivePairMM-1));

        if boolThreeD
            vecSepDistSq = vecSepX.^2 + vecSepY.^2 + vecSepZ.^2;
        else
            vecSepDistSq = vecSepX.^2 + vecSepY.^2;
        end
        boolContact = vecSepDistSq < vecContactDist.^2;  % [scalNumPairs x 1] overlapping pairs only

        % Trim all arrays to only pairs in contac
        vecSepX = vecSepX(boolContact); % [scalNumContacts x 1]
        vecSepY = vecSepY(boolContact); % [scalNumContacts x 1]
        if boolThreeD
            vecSepZ = vecSepZ(boolContact);
        end
        vecSepDistSq = vecSepDistSq(boolContact); % [scalNumContacts x 1]
        vecContactDist = vecContactDist(boolContact); % [scalNumContacts x 1]
        vecContactNN = vecActivePairNN(boolContact);% [scalNumContacts x 1]
        vecContactMM = vecActivePairMM(boolContact);% [scalNumContacts x 1]

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

        % Contact count per particle (coordination number)
        vecCoordNum = accumarray(vecContactNN, 1, [N 1]) ...
                    + accumarray(vecContactMM, 1, [N 1]);

        vecPotentialEnergyHistory(nt) = sum(vecPotentialContact) / N;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Drag, boundaries, energy %%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        vecForceX = vecForceX - scalDissipationAbsolute .* vecVelX;  % [N x 1]
        vecForceY = vecForceY - scalDissipationAbsolute .* vecVelY;  % [N x 1]
        if boolThreeD
            vecForceZ = vecForceZ - scalDissipationAbsolute .* vecVelZ;
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

        vecAccelX = vecForceX ./ M;
        vecAccelY = vecForceY ./ M - scalGravity;
        if boolThreeD
            vecAccelZ = vecForceZ ./ M;  % no gravity in Z
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Second step in Verlet integration %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        vecVelX = vecVelX + (vecAccelXPrev + vecAccelX) .* (scalTimestep/2);  % [N x 1]
        vecVelY = vecVelY + (vecAccelYPrev + vecAccelY) .* (scalTimestep/2);  % [N x 1]
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

        % Pressure estimate from mean potential energy
        scalEp = vecPotentialEnergyHistory(nt);
        if options.hertzian
            scalPressure = (scalEp * (5/2) / K)^(2/5); % this has an implied d= 1 in the denominator
        else
            scalPressure = sqrt(2 * scalEp / K); % this has an implied d= 1 in the denominator
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% COMPRESSION DECISIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        scalEk = vecKineticEnergyHistory(nt);

        if boolFastCompressPhase
            if scalPressure < P_target/50
                scalBoxWidthX= scalBoxWidthX * (1-scalCompressionRateFast);
                scalBoxHeightY = scalBoxHeightY * (1-scalCompressionRateFast);
                vecPosX = vecPosX * (1-scalCompressionRateFast);  % [N x 1]
                vecPosY = vecPosY * (1-scalCompressionRateFast);  % [N x 1]
                if boolThreeD
                    scalBoxDepthZ = scalBoxDepthZ * (1-scalCompressionRateFast);
                    vecPosZ = vecPosZ * (1-scalCompressionRateFast);
                end
                boolCellUpdateNeeded = true;
                scalLastCompressStep = nt;
            elseif scalPressure < P_target && scalEk < 1e-8
                scalBoxWidthX  = scalBoxWidthX  * (1-scalCompressionRateFast);
                scalBoxHeightY = scalBoxHeightY * (1-scalCompressionRateFast);
                vecPosX = vecPosX * (1-scalCompressionRateFast);  % [N x 1]
                vecPosY = vecPosY * (1-scalCompressionRateFast);  % [N x 1]
                if boolThreeD
                    scalBoxDepthZ = scalBoxDepthZ * (1-scalCompressionRateFast);
                    vecPosZ = vecPosZ * (1-scalCompressionRateFast);
                end
                boolCellUpdateNeeded = true;
                scalLastCompressStep = nt;
            elseif scalPressure > P_target && scalEk < 1e-10 && nt > (scalLastCompressStep+100)
                scalBoxWidthX  = scalBoxWidthX  * (1+scalCompressionRateFast);
                scalBoxHeightY = scalBoxHeightY * (1+scalCompressionRateFast);
                vecPosX = vecPosX * (1+scalCompressionRateFast);  % [N x 1]
                vecPosY = vecPosY * (1+scalCompressionRateFast);  % [N x 1]
                if boolThreeD
                    scalBoxDepthZ = scalBoxDepthZ * (1+scalCompressionRateFast);
                    vecPosZ = vecPosZ * (1+scalCompressionRateFast);
                end
                boolCellUpdateNeeded = true;
                scalLastCompressStep = nt;
                boolFastCompressPhase = false;
            end
        else
            if scalPressure < scalPressureFastGrow
                scalBoxWidthX  = scalBoxWidthX  * (1-scalCompressionRate);
                scalBoxHeightY = scalBoxHeightY * (1-scalCompressionRate);
                vecPosX = vecPosX * (1-scalCompressionRate);  % [N x 1]
                vecPosY = vecPosY * (1-scalCompressionRate);  % [N x 1]
                if boolThreeD
                    scalBoxDepthZ = scalBoxDepthZ * (1-scalCompressionRate);
                    vecPosZ = vecPosZ * (1-scalCompressionRate);
                end
                boolCellUpdateNeeded = true;
                scalLastCompressStep = nt;
            elseif scalPressure < P_target && scalEk < 1e-8
                scalBoxWidthX  = scalBoxWidthX * (1-scalCompressionRate);
                scalBoxHeightY = scalBoxHeightY * (1-scalCompressionRate);
                vecPosX = vecPosX * (1-scalCompressionRate);  % [N x 1]
                vecPosY = vecPosY * (1-scalCompressionRate);  % [N x 1]
                if boolThreeD
                    scalBoxDepthZ = scalBoxDepthZ * (1-scalCompressionRate);
                    vecPosZ = vecPosZ * (1-scalCompressionRate);
                end
                boolCellUpdateNeeded = true;
                scalLastCompressStep = nt;
            elseif scalPressure > P_target && scalEk < 1e-20
                fprintf('Converged at step %d | P=%.4e\n', nt, scalPressure);
                break;
            end
        end
    end
    fprintf('Loop finished at step %d.\n', nt);

    %% Remove rattlers before saving
    if boolThreeD
        scalVolumeSpheres = sum((4/3)*pi*(vecDiameter/2).^3);
        scalVolumeBox = scalBoxWidthX * scalBoxHeightY * scalBoxDepthZ;
    else
        scalVolumeSpheres = sum(pi*(vecDiameter/2).^2);
        scalVolumeBox = scalBoxWidthX * scalBoxHeightY;
    end
    fprintf('Packing fraction before cleanRats: %.4f\n', scalVolumeSpheres/scalVolumeBox);

    fprintf('Running cleanRats...\n');
    if boolThreeD
        fprintf('Box dims: Lx=%.4f, Ly=%.4f, Lz=%.4f\n', scalBoxWidthX, scalBoxHeightY, scalBoxDepthZ);
    else
        fprintf('Box dims: Lx=%.4f, Ly=%.4f\n', scalBoxWidthX, scalBoxHeightY);
    end
    fprintf('Radii range: min=%.4f, max=%.4f\n', min(vecDiameter/2), max(vecDiameter/2));
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
        [matPositions, vecRadii] = cleanRats(matPositions, vecRadii, scalBoxHeightY, scalBoxWidthX);
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

%% Final plot
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
    end
    drawnow;
    hold off;

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
        scalVolumeSpheres_clean = sum((4/3) * pi * (vecDiameter/2).^3);
        scalVolumeBox_clean = scalBoxWidthX * scalBoxHeightY * scalBoxDepthZ;
    else
        scalVolumeSpheres_clean = sum(pi * (vecDiameter/2).^2);
        scalVolumeBox_clean = scalBoxWidthX * scalBoxHeightY;
    end
    packingFraction = scalVolumeSpheres_clean / scalVolumeBox_clean;
    fprintf('Packing fraction after cleanRats: %.4f\n', packingFraction);

    if calc_eig
        % 3D hessian not yet implemented — skip eigenmodes
        if boolThreeD
            warning('calc_eig not supported for 3D yet — saving positions only.');
            if options.hertzian
                save(strFilename, 'vecPosX', 'vecPosY', 'vecPosZ', 'vecDiameter', ...
                    'scalBoxWidthX', 'scalBoxHeightY', 'scalBoxDepthZ', ...
                    'K', 'P_target', 'scalPressure', 'N', 'N_original', 'packingFraction', ...
                    'vecHertzNN', 'vecHertzMM', 'vecHertzKeff');
            else
                save(strFilename, 'vecPosX', 'vecPosY', 'vecPosZ', 'vecDiameter', ...
                    'scalBoxWidthX', 'scalBoxHeightY', 'scalBoxDepthZ', ...
                    'K', 'P_target', 'scalPressure', 'N', 'N_original', 'packingFraction');
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
                    'packingFraction', 'matEigenVectors', 'matEigenValues', ...
                    'vecHertzNN', 'vecHertzMM', 'vecHertzKeff');
            else
                save(strFilename, 'vecPosX', 'vecPosY', 'vecDiameter', ...
                    'scalBoxWidthX', 'scalBoxHeightY', 'K', 'P_target', 'scalPressure', 'N', 'N_original', ...
                    'packingFraction', 'matEigenVectors', 'matEigenValues');
            end
        end
    else
        if boolThreeD
            if options.hertzian
                save(strFilename, 'vecPosX', 'vecPosY', 'vecPosZ', 'vecDiameter', ...
                    'scalBoxWidthX', 'scalBoxHeightY', 'scalBoxDepthZ', ...
                    'K', 'P_target', 'scalPressure', 'N', 'N_original', 'packingFraction', ...
                    'vecHertzNN', 'vecHertzMM', 'vecHertzKeff');
            else
                save(strFilename, 'vecPosX', 'vecPosY', 'vecPosZ', 'vecDiameter', ...
                    'scalBoxWidthX', 'scalBoxHeightY', 'scalBoxDepthZ', ...
                    'K', 'P_target', 'scalPressure', 'N', 'N_original', 'packingFraction');
            end
        else
            if options.hertzian
                save(strFilename, 'vecPosX', 'vecPosY', 'vecDiameter', ...
                    'scalBoxWidthX', 'scalBoxHeightY', 'K', 'P_target', 'scalPressure', ...
                    'N', 'N_original', 'packingFraction', ...
                    'vecHertzNN', 'vecHertzMM', 'vecHertzKeff');
            else
                save(strFilename, 'vecPosX', 'vecPosY', 'vecDiameter', ...
                    'scalBoxWidthX', 'scalBoxHeightY', 'K', 'P_target', 'scalPressure', ...
                    'N', 'N_original', 'packingFraction');
            end
        end
    end

    fprintf('File saved to: %s\n', strFilename);

    if x_mult ~= 1 || y_mult ~= 1
        if ~boolThreeD
            packRepeatTile(N_original, K, P_target, scalRoundedWidth, seed, x_mult, y_mult, ...
                calc_eig, save_path, save_path);
            disp("Tile saved to: " + strFilename);
        else
            warning('packRepeatTile not supported for 3D yet.');
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

function [vecPosX, vecPosY, vecPosZ, cellParticleList, scalNumCellsX, scalNumCellsY, scalNumCellsZ] = ...
        rebuildCellList3D(vecPosX, vecPosY, vecPosZ, scalBoxWidthX, scalBoxHeightY, scalBoxDepthZ, scalRawCellWidth, scalTimestep, N)

    scalNumCellsX  = round(scalBoxWidthX  / scalRawCellWidth);
    scalCellWidthX = scalBoxWidthX  / scalNumCellsX;
    scalNumCellsY  = round(scalBoxHeightY / scalRawCellWidth);
    scalCellWidthY = scalBoxHeightY / scalNumCellsY;
    scalNumCellsZ  = round(scalBoxDepthZ  / scalRawCellWidth);
    scalCellWidthZ = scalBoxDepthZ  / scalNumCellsZ;

    % Sanity check
    vecFloorX = floor(vecPosX / scalBoxWidthX);
    vecFloorY = floor(vecPosY / scalBoxHeightY);
    vecFloorZ = floor(vecPosZ / scalBoxDepthZ);
    if any(abs(vecFloorX) > 1) || any(abs(vecFloorY) > 1) || any(abs(vecFloorZ) > 1)
        error(['Particle moved more than one box length in a single timestep.\n' ...
               'Reduce scalTimestep (currently %.4e).'], scalTimestep);
    end

    % Wrap into [0, L)
    vecPosX = mod(vecPosX, scalBoxWidthX);
    vecPosY = mod(vecPosY, scalBoxHeightY);
    vecPosZ = mod(vecPosZ, scalBoxDepthZ);

    % Map to cell indices
    vecCellIdxX = min(max(ceil(vecPosX / scalCellWidthX), 1), scalNumCellsX);
    vecCellIdxY = min(max(ceil(vecPosY / scalCellWidthY), 1), scalNumCellsY);
    vecCellIdxZ = min(max(ceil(vecPosZ / scalCellWidthZ), 1), scalNumCellsZ);

    % Linear index: x + Nx*(y-1) + Nx*Ny*(z-1)
    vecCellLinearIdx = vecCellIdxX + scalNumCellsX*(vecCellIdxY-1) + scalNumCellsX*scalNumCellsY*(vecCellIdxZ-1);
    cellParticleList = accumarray(vecCellLinearIdx, (1:N)', [scalNumCellsX*scalNumCellsY*scalNumCellsZ 1], @(x){x});
    cellParticleList = reshape(cellParticleList, scalNumCellsX, scalNumCellsY, scalNumCellsZ);
end

function [vecPairNN, vecPairMM, scalNumPairs, scalMaxPairs] = findNeighbors3D( ...
        cellParticleList, scalNumCellsX, scalNumCellsY, scalNumCellsZ, ...
        vecPairNN, vecPairMM, scalMaxPairs)

    scalNumPairs = 0;

    for idxCellX = 1:scalNumCellsX
        for idxCellY = 1:scalNumCellsY
            for idxCellZ = 1:scalNumCellsZ

                scalCellLeft  = mod(idxCellX-2, scalNumCellsX)+1;
                scalCellRight = mod(idxCellX,   scalNumCellsX)+1;
                scalCellDown  = mod(idxCellY-2, scalNumCellsY)+1;
                scalCellUp    = mod(idxCellY,   scalNumCellsY)+1;
                scalCellBack  = mod(idxCellZ-2, scalNumCellsZ)+1;
                scalCellFront = mod(idxCellZ,   scalNumCellsZ)+1;

                vecCurrentCell = cellParticleList{idxCellX, idxCellY, idxCellZ};

                vecNeighborList = [ ...
                    cellParticleList{scalCellLeft,  scalCellDown, scalCellBack};  ...
                    cellParticleList{scalCellLeft,  scalCellDown, idxCellZ};      ...
                    cellParticleList{scalCellLeft,  scalCellDown, scalCellFront}; ...
                    cellParticleList{scalCellLeft,  idxCellY,     scalCellBack};  ...
                    cellParticleList{scalCellLeft,  idxCellY,     idxCellZ};      ...
                    cellParticleList{scalCellLeft,  idxCellY,     scalCellFront}; ...
                    cellParticleList{scalCellLeft,  scalCellUp,   scalCellBack};  ...
                    cellParticleList{scalCellLeft,  scalCellUp,   idxCellZ};      ...
                    cellParticleList{scalCellLeft,  scalCellUp,   scalCellFront}; ...
                    cellParticleList{idxCellX,      scalCellDown, scalCellBack};  ...
                    cellParticleList{idxCellX,      scalCellDown, idxCellZ};      ...
                    cellParticleList{idxCellX,      scalCellDown, scalCellFront}; ...
                    cellParticleList{idxCellX,      idxCellY,     scalCellBack};  ...
                    vecCurrentCell;                                                ...
                    cellParticleList{idxCellX,      idxCellY,     scalCellFront}; ...
                    cellParticleList{idxCellX,      scalCellUp,   scalCellBack};  ...
                    cellParticleList{idxCellX,      scalCellUp,   idxCellZ};      ...
                    cellParticleList{idxCellX,      scalCellUp,   scalCellFront}; ...
                    cellParticleList{scalCellRight, scalCellDown, scalCellBack};  ...
                    cellParticleList{scalCellRight, scalCellDown, idxCellZ};      ...
                    cellParticleList{scalCellRight, scalCellDown, scalCellFront}; ...
                    cellParticleList{scalCellRight, idxCellY,     scalCellBack};  ...
                    cellParticleList{scalCellRight, idxCellY,     idxCellZ};      ...
                    cellParticleList{scalCellRight, idxCellY,     scalCellFront}; ...
                    cellParticleList{scalCellRight, scalCellUp,   scalCellBack};  ...
                    cellParticleList{scalCellRight, scalCellUp,   idxCellZ};      ...
                    cellParticleList{scalCellRight, scalCellUp,   scalCellFront}];

                for idxNN = vecCurrentCell'
                    vecCandidates = vecNeighborList(vecNeighborList > idxNN);
                    scalNumCandidates = numel(vecCandidates);
                    if scalNumCandidates == 0; continue; end
                    if scalNumPairs + scalNumCandidates > scalMaxPairs
                        scalMaxPairs = 2 * scalMaxPairs;
                        vecPairNN(scalMaxPairs) = 0;
                        vecPairMM(scalMaxPairs) = 0;
                    end
                    vecPairNN(scalNumPairs+1 : scalNumPairs+scalNumCandidates) = idxNN;
                    vecPairMM(scalNumPairs+1 : scalNumPairs+scalNumCandidates) = vecCandidates;
                    scalNumPairs = scalNumPairs + scalNumCandidates;
                end
            end
        end
    end
end

function [vecPairNN, vecPairMM, scalNumPairs, scalMaxPairs] = findNeighbors2D( ...
        cellParticleList, scalNumCellsX, scalNumCellsY, ...
        vecPairNN, vecPairMM, scalMaxPairs)

    scalNumPairs = 0;

    for idxCellX = 1:scalNumCellsX
        for idxCellY = 1:scalNumCellsY

            scalCellLeft  = mod(idxCellX-2, scalNumCellsX)+1;
            scalCellRight = mod(idxCellX,   scalNumCellsX)+1;
            scalCellDown  = mod(idxCellY-2, scalNumCellsY)+1;
            scalCellUp    = mod(idxCellY,   scalNumCellsY)+1;

            vecCurrentCell = cellParticleList{idxCellX, idxCellY};

            vecNeighborList = [ ...
                cellParticleList{scalCellLeft,  scalCellDown}; ...
                cellParticleList{scalCellLeft,  idxCellY};     ...
                cellParticleList{scalCellLeft,  scalCellUp};   ...
                cellParticleList{idxCellX,      scalCellDown}; ...
                vecCurrentCell;                                ...
                cellParticleList{idxCellX,      scalCellUp};   ...
                cellParticleList{scalCellRight, scalCellDown}; ...
                cellParticleList{scalCellRight, idxCellY};     ...
                cellParticleList{scalCellRight, scalCellUp}];

            for idxNN = vecCurrentCell'
                vecCandidates = vecNeighborList(vecNeighborList > idxNN);
                scalNumCandidates = numel(vecCandidates);
                if scalNumCandidates == 0; continue; end
                if scalNumPairs + scalNumCandidates > scalMaxPairs
                    scalMaxPairs = 2 * scalMaxPairs;
                    vecPairNN(scalMaxPairs) = 0;
                    vecPairMM(scalMaxPairs) = 0;
                end
                vecPairNN(scalNumPairs+1 : scalNumPairs+scalNumCandidates) = idxNN;
                vecPairMM(scalNumPairs+1 : scalNumPairs+scalNumCandidates) = vecCandidates;
                scalNumPairs = scalNumPairs + scalNumCandidates;
            end
        end
    end
end
