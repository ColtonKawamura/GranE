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

    tic

    rng(seed)

    %% Box and particle setup
    if boolThreeD
        scalBoxWidth  = 2*N^(1/3)*D;   % Lx
        scalBoxHeight = 2*N^(1/3)*D;   % Ly
        scalBoxDepth  = 2*N^(1/3)*D;   % Lz
    else
        scalBoxWidth  = 2*sqrt(N)*D;   % Lx
        scalBoxHeight = 2*sqrt(N)*D;   % Ly
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

    % time step shoule be 1/100 of a particle oscillation period
    scalTimestep = 2*pi * sqrt(M/K) * 0.01;
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
    [vecPosX, vecPosY, vecPosZ] = ndgrid(D/2 : G*D : scalBoxWidth-D/2, ...
                                          D/2 : G*D : scalBoxHeight-D/2, ...
                                          D/2 : G*D : scalBoxDepth-D/2);
    else
        [vecPosX, vecPosY] = ndgrid(D/2 : G*D : scalBoxWidth-D/2, ...
                                    D/2 : G*D : scalBoxHeight-D/2);
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
    scalNumCellsX  = round(scalBoxWidth  / scalRawCellWidth);
    scalCellWidthX = scalBoxWidth  / scalNumCellsX;
    scalNumCellsY  = round(scalBoxHeight / scalRawCellWidth);
    scalCellWidthY = scalBoxHeight / scalNumCellsY;
    if boolThreeD
        scalNumCellsZ  = round(scalBoxDepth  / scalRawCellWidth);
        scalCellWidthZ = scalBoxDepth  / scalNumCellsZ;
    end

    % Wrap positions into [0, L) before rebuilding cell list
    % mod(x, L) == x - L*floor(x/L), handles both positive and negative overshoot
    vecPosX = mod(vecPosX, scalBoxWidth);   % [N x 1]
    vecPosY = mod(vecPosY, scalBoxHeight);  % [N x 1]
    if boolThreeD
        vecPosZ = mod(vecPosZ, scalBoxDepth);
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
        axis equal; axis([0 scalBoxWidth 0 scalBoxHeight]);
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
    scalNextLogFrac = 0.01;  % next threshold to print
    for nt = 1:scalMaxSteps

        % Progress logging
        fracDone = nt / scalMaxSteps;
        if fracDone >= scalNextLogFrac
            fprintf('  %d%% complete | step %d | P=%.4e | Ek=%.4e\n', ...
                round(scalNextLogFrac*100), nt, scalPressure, scalEk);
            scalNextLogFrac = scalNextLogFrac + 0.1;
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
                ylim([0 scalBoxHeight]); xlim([0 scalBoxWidth]);
                title(num2str(scalBoxHeight));
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

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Re-assign particles to cells %%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if boolCellUpdateNeeded || mod(nt, scalCellUpdateInterval) == 0
            if boolThreeD
                [vecPosX, vecPosY, vecPosZ, cellParticleList, scalNumCellsX, scalNumCellsY, scalNumCellsZ] = ...
                    rebuildCellList3D(vecPosX, vecPosY, vecPosZ, scalBoxWidth, scalBoxHeight, scalBoxDepth, scalRawCellWidth, scalTimestep, N);
            else
                [vecPosX, vecPosY, cellParticleList, scalNumCellsX, scalNumCellsY] = ...
                    rebuildCellList(vecPosX, vecPosY, scalBoxWidth, scalBoxHeight, scalRawCellWidth, scalTimestep, N);
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
        vecSepX = vecSepX - scalBoxWidth  * round(vecSepX / scalBoxWidth);
        vecSepY = vecPosY(vecActivePairMM) - vecPosY(vecActivePairNN);
        vecSepY = vecSepY - scalBoxHeight * round(vecSepY / scalBoxHeight);
        if boolThreeD
            vecSepZ = vecPosZ(vecActivePairMM) - vecPosZ(vecActivePairNN);
            vecSepZ = vecSepZ - scalBoxDepth * round(vecSepZ / scalBoxDepth);
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
            vecRadiiEff = (vecRadiiNN .* vecRadiiMM) ./ (vecRadiiNN + vecRadiiMM); % [scalNumContacts x 1]
        end

        vecSepDist = sqrt(vecSepDistSq); % [scalNumContacts x 1]
        vecOverlap = vecContactDist - vecSepDist; % [scalNumContacts x 1] positive when overlapping

        % Force magnitude and potential energy per contact
        if options.hertzian
            vecForceMag = -K .* vecOverlap.^(3/2); % [scalNumContacts x 1]
            vecPotentialContact = (2/5) * K .* vecOverlap.^(5/2);% [scalNumContacts x 1]
        else
            vecForceMag = -K .* vecOverlap; % [scalNumContacts x 1]
            vecPotentialContact = 0.5 * K .* vecOverlap.^2;      % [scalNumContacts x 1]
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

        boolLeftWallContact  = vecPosX < vecDiameter/2;
        boolRightWallContact = vecPosX > scalBoxWidth - vecDiameter/2;

        vecPosX = vecPosX - scalBoxWidth  .* floor(vecPosX / scalBoxWidth);
        vecPosY = vecPosY - scalBoxHeight .* floor(vecPosY / scalBoxHeight);
        if boolThreeD
            vecPosZ = vecPosZ - scalBoxDepth .* floor(vecPosZ / scalBoxDepth);
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
            vecPosZ = vecPosZ + vecVelZ*scalTimestep + vecAccelZPrev.*(scalTimestep^2/2);
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
            scalPressure = (scalEp * (5/2) / K)^(2/5);
        else
            scalPressure = sqrt(2 * scalEp / K);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% COMPRESSION DECISIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        scalEk = vecKineticEnergyHistory(nt);

        if boolFastCompressPhase
            if scalPressure < P_target/50
                scalBoxWidth= scalBoxWidth * (1-scalCompressionRateFast);
                scalBoxHeight = scalBoxHeight * (1-scalCompressionRateFast);
                vecPosX = vecPosX * (1-scalCompressionRateFast);  % [N x 1]
                vecPosY = vecPosY * (1-scalCompressionRateFast);  % [N x 1]
                if boolThreeD
                    scalBoxDepth = scalBoxDepth * (1-scalCompressionRateFast);
                    vecPosZ = vecPosZ * (1-scalCompressionRateFast);
                end
                boolConverged = true;
                boolCellUpdateNeeded = true;
                scalLastCompressStep = nt;
            elseif scalPressure < P_target && scalEk < 1e-8
                scalBoxWidth  = scalBoxWidth  * (1-scalCompressionRateFast);
                scalBoxHeight = scalBoxHeight * (1-scalCompressionRateFast);
                vecPosX = vecPosX * (1-scalCompressionRateFast);  % [N x 1]
                vecPosY = vecPosY * (1-scalCompressionRateFast);  % [N x 1]
                if boolThreeD
                    scalBoxDepth = scalBoxDepth * (1-scalCompressionRateFast);
                    vecPosZ = vecPosZ * (1-scalCompressionRateFast);
                end
                boolConverged = true;
                boolCellUpdateNeeded = true;
                scalLastCompressStep = nt;
            elseif scalPressure > P_target && scalEk < 1e-10 && nt > (scalLastCompressStep+100)
                scalBoxWidth  = scalBoxWidth  * (1+scalCompressionRateFast);
                scalBoxHeight = scalBoxHeight * (1+scalCompressionRateFast);
                vecPosX = vecPosX * (1+scalCompressionRateFast);  % [N x 1]
                vecPosY = vecPosY * (1+scalCompressionRateFast);  % [N x 1]
                if boolThreeD
                    scalBoxDepth = scalBoxDepth * (1+scalCompressionRateFast);
                    vecPosZ = vecPosZ * (1+scalCompressionRateFast);
                end
                boolConverged = true;
                boolCellUpdateNeeded = true;
                scalLastCompressStep = nt;
                boolFastCompressPhase = false;
            end
        else
            if scalPressure < scalPressureFastGrow
                scalBoxWidth  = scalBoxWidth  * (1-scalCompressionRate);
                scalBoxHeight = scalBoxHeight * (1-scalCompressionRate);
                vecPosX = vecPosX * (1-scalCompressionRate);  % [N x 1]
                vecPosY = vecPosY * (1-scalCompressionRate);  % [N x 1]
                if boolThreeD
                    scalBoxDepth = scalBoxDepth * (1-scalCompressionRate);
                    vecPosZ = vecPosZ * (1-scalCompressionRate);
                end
                boolConverged = true;
                boolCellUpdateNeeded = true;
                scalLastCompressStep = nt;
            elseif scalPressure < P_target && scalEk < 1e-8
                scalBoxWidth  = scalBoxWidth * (1-scalCompressionRate);
                scalBoxHeight = scalBoxHeight * (1-scalCompressionRate);
                vecPosX = vecPosX * (1-scalCompressionRate);  % [N x 1]
                vecPosY = vecPosY * (1-scalCompressionRate);  % [N x 1]
                if boolThreeD
                    scalBoxDepth = scalBoxDepth * (1-scalCompressionRate);
                    vecPosZ = vecPosZ * (1-scalCompressionRate);
                end
                boolConverged = true;
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
    fprintf('Running cleanRats...\n');
    if boolThreeD
        matPositions = [vecPosX, vecPosY, vecPosZ];
        vecRadii = vecDiameter ./ 2;
        [matPositions, vecRadii] = cleanRats(matPositions, vecRadii, scalBoxHeight, scalBoxWidth, scalBoxDepth);
        vecPosX = matPositions(:,1);
        vecPosY = matPositions(:,2);
        vecPosZ = matPositions(:,3);
    else
        matPositions = [vecPosX, vecPosY];
        vecRadii = vecDiameter ./ 2;
        [matPositions, vecRadii] = cleanRats(matPositions, vecRadii, scalBoxHeight, scalBoxWidth);
        vecPosX = matPositions(:,1);
        vecPosY = matPositions(:,2);
    end
    vecDiameter = vecRadii .* 2;
    N_clean = size(matPositions, 1);
    fprintf('cleanRats complete. %d particles remaining (of %d original).\n', N_clean, N);

    %% Final plot
    figure;
    hold on;
    if boolThreeD
        [sx, sy, sz] = sphere(16);

        % Main particles
        for np = 1:N
            r = vecDiameter(np)/2;
            surf(r*sx + vecPosX(np), r*sy + vecPosY(np), r*sz + vecPosZ(np), ...
                'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.6);
        end

        % Ghost particles on +x, +y, +z faces
        vecOffsets = [scalBoxWidth, 0, 0; ...
                      0, scalBoxHeight, 0; ...
                      0, 0, scalBoxDepth];  % [3 x 3] one offset per face

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
        axis([-scalBoxWidth*0.0 2*scalBoxWidth ...
              -scalBoxHeight*0.0 2*scalBoxHeight ...
              -scalBoxDepth*0.0  2*scalBoxDepth]);
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
        vecOffsets2D = [scalBoxWidth,  0; ...
                       -scalBoxWidth,  0; ...
                        0,  scalBoxHeight; ...
                        0, -scalBoxHeight; ...
                        scalBoxWidth,  scalBoxHeight; ...
                       -scalBoxWidth,  scalBoxHeight; ...
                        scalBoxWidth, -scalBoxHeight; ...
                       -scalBoxWidth, -scalBoxHeight];

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
        axis([-scalBoxWidth 2*scalBoxWidth -scalBoxHeight 2*scalBoxHeight]);
    end
    drawnow;
    hold off;
%% Save results

    if boolThreeD
        scalRoundedWidth = round(N^(1/3));
        strFilename = sprintf('%s3D_N%d_P%s_Width%d_Seed%d.mat', ...
            save_path, N, num2str(P_target), scalRoundedWidth, seed);
    else
        scalRoundedWidth = round(sqrt(N));
        strFilename = sprintf('%s2D_N%d_P%s_Width%d_Seed%d.mat', ...
            save_path, N, num2str(P_target), scalRoundedWidth, seed);
    end

    if calc_eig
        if boolThreeD
            % 3D hessian not yet implemented — skip eigenmodes
            warning('calc_eig not supported for 3D yet — saving positions only.');
            save(strFilename, 'vecPosX', 'vecPosY', 'vecPosZ', 'vecDiameter', ...
                'scalBoxWidth', 'scalBoxHeight', 'scalBoxDepth', ...
                'K', 'P_target', 'scalPressure', 'N');
        else
            matPositions = [vecPosX, vecPosY];
            vecRadii = vecDiameter ./ 2;
            [matPositions, vecRadii] = cleanRats(matPositions, vecRadii, K, scalBoxHeight, scalBoxWidth);
            matHessian = hess2d(matPositions, vecRadii, K, scalBoxHeight, scalBoxWidth);
            [matEigenVectors, matEigenValues] = eig(matHessian);
            save(strFilename, 'vecPosX', 'vecPosY', 'vecDiameter', ...
                'scalBoxWidth', 'scalBoxHeight', 'K', 'P_target', 'scalPressure', 'N', ...
                'matEigenVectors', 'matEigenValues');
        end
    else
        if boolThreeD
            save(strFilename, 'vecPosX', 'vecPosY', 'vecPosZ', 'vecDiameter', ...
                'scalBoxWidth', 'scalBoxHeight', 'scalBoxDepth', ...
                'K', 'P_target', 'scalPressure', 'N');
        else
            save(strFilename, 'vecPosX', 'vecPosY', 'vecDiameter', ...
                'scalBoxWidth', 'scalBoxHeight', 'K', 'P_target', 'scalPressure', 'N');
        end
    end

    fprintf('File saved to: %s\n', strFilename);

    if x_mult ~= 1 || y_mult ~= 1
        if ~boolThreeD
            packRepeatTile(N, K, P_target, scalRoundedWidth, seed, x_mult, y_mult, ...
                calc_eig, save_path, save_path);
            disp("Tile saved to: " + strFilename);
        else
            warning('packRepeatTile not supported for 3D yet.');
        end
    end

    toc
end

function [vecPosX, vecPosY, cellParticleList, scalNumCellsX, scalNumCellsY] = ...
        rebuildCellList(vecPosX, vecPosY, scalBoxWidth, scalBoxHeight, scalRawCellWidth, scalTimestep, N)

    scalNumCellsX  = round(scalBoxWidth  / scalRawCellWidth);
    scalCellWidthX = scalBoxWidth  / scalNumCellsX;
    scalNumCellsY  = round(scalBoxHeight / scalRawCellWidth);
    scalCellWidthY = scalBoxHeight / scalNumCellsY;

    % Sanity check: no particle should move more than one box length per step
    vecFloorX = floor(vecPosX / scalBoxWidth);
    vecFloorY = floor(vecPosY / scalBoxHeight);
    if any(abs(vecFloorX) > 1) || any(abs(vecFloorY) > 1)
        error(['Particle moved more than one box length in a single timestep.\n' ...
               'Max x overshoot: %.2f box lengths\n' ...
               'Max y overshoot: %.2f box lengths\n' ...
               'Reduce scalTimestep (currently %.4e).'], ...
               max(abs(vecFloorX)), max(abs(vecFloorY)), scalTimestep);
    end

    % Wrap positions into [0, L) before rebuilding
    vecPosX = mod(vecPosX, scalBoxWidth);   % [N x 1]
    vecPosY = mod(vecPosY, scalBoxHeight);  % [N x 1]

    % Map to cell indices, clamped to [1, scalNumCells]
    vecCellIdxX = min(max(ceil(vecPosX / scalCellWidthX), 1), scalNumCellsX);  % [N x 1]
    vecCellIdxY = min(max(ceil(vecPosY / scalCellWidthY), 1), scalNumCellsY);  % [N x 1]

    % O(N) bucketing via accumarray
    vecCellLinearIdx = vecCellIdxX + scalNumCellsX * (vecCellIdxY - 1);        % [N x 1]
    cellParticleList = accumarray(vecCellLinearIdx, (1:N)', [scalNumCellsX*scalNumCellsY 1], @(x){x});
    cellParticleList = reshape(cellParticleList, scalNumCellsX, scalNumCellsY); % [scalNumCellsX x scalNumCellsY]

end

function [vecPosX, vecPosY, vecPosZ, cellParticleList, scalNumCellsX, scalNumCellsY, scalNumCellsZ] = ...
        rebuildCellList3D(vecPosX, vecPosY, vecPosZ, scalBoxWidth, scalBoxHeight, scalBoxDepth, scalRawCellWidth, scalTimestep, N)

    scalNumCellsX  = round(scalBoxWidth  / scalRawCellWidth);
    scalCellWidthX = scalBoxWidth  / scalNumCellsX;
    scalNumCellsY  = round(scalBoxHeight / scalRawCellWidth);
    scalCellWidthY = scalBoxHeight / scalNumCellsY;
    scalNumCellsZ  = round(scalBoxDepth  / scalRawCellWidth);
    scalCellWidthZ = scalBoxDepth  / scalNumCellsZ;

    % Sanity check
    vecFloorX = floor(vecPosX / scalBoxWidth);
    vecFloorY = floor(vecPosY / scalBoxHeight);
    vecFloorZ = floor(vecPosZ / scalBoxDepth);
    if any(abs(vecFloorX) > 1) || any(abs(vecFloorY) > 1) || any(abs(vecFloorZ) > 1)
        error(['Particle moved more than one box length in a single timestep.\n' ...
               'Reduce scalTimestep (currently %.4e).'], scalTimestep);
    end

    % Wrap into [0, L)
    vecPosX = mod(vecPosX, scalBoxWidth);
    vecPosY = mod(vecPosY, scalBoxHeight);
    vecPosZ = mod(vecPosZ, scalBoxDepth);

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
