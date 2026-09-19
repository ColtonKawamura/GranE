function packRepeatTile(N, K, P_target, scalWidthFactor, seed, scalXMult, scalYMult, calc_eig, strInPath, strSavePath)
    % packRepeatTile  Load a saved 2D packing and tile it in x and/or y.
    %
    % N,              Number of particles in the base tile
    % K,              Spring constant
    % P_target,       Target pressure
    % scalWidthFactor, Width factor used in the base tile filename
    % seed,           RNG seed used in the base tile filename
    % scalXMult,      Number of times to tile in x
    % scalYMult,      Number of times to tile in y
    % calc_eig,       Boolean: compute and save eigenmodes
    % strInPath,      Path to the base tile .mat file
    % strSavePath,    Path to save the tiled output .mat file
    %
    % Example:
    % packRepeatTile(400, 100, 0.05, 22, 1, 2, 2, false, 'in/tiles/', 'out/tiles/')

    %% Load base tile
    strBaseName   = sprintf('2D_N%d_P%s_Width%d_Seed%d', N, num2str(P_target), scalWidthFactor, seed);
    strFilenameIn = strInPath + strBaseName + ".mat";

    fprintf('Loading file: %s\n', strFilenameIn);

    try
        load(strFilenameIn);
        fprintf('Load SUCCESS: %s\n', strFilenameIn);
    catch ME
        fprintf('Load FAILED: %s\n', strFilenameIn);
        fprintf('Error message: %s\n', ME.message);
    end

    %% Tile in x
    % Shift copies of vecPosX by multiples of scalBoxWidthX, keep vecPosY unchanged
    vecPosXTiled = zeros(N * scalXMult, 1);       % [N*scalXMult x 1]
    vecPosYTiled = zeros(N * scalXMult, 1);       % [N*scalXMult x 1]
    vecDiameterTiled = zeros(N * scalXMult, 1);   % [N*scalXMult x 1]

    for ii = 0:scalXMult-1
        idxStart = ii*N + 1;
        idxEnd   = ii*N + N;
        vecPosXTiled(idxStart:idxEnd)    = vecPosX + ii * scalBoxWidthX;  % shift x only
        vecPosYTiled(idxStart:idxEnd)    = vecPosY;                      % y unchanged
        vecDiameterTiled(idxStart:idxEnd) = vecDiameter;
    end

    scalBoxWidthXTiled = scalBoxWidthX * scalXMult;
    scalNTiled        = N * scalXMult;

    %% Tile in y
    % Shift copies of vecPosYTiled by multiples of scalBoxHeightY, keep x unchanged
    vecPosXFinal    = zeros(scalNTiled * scalYMult, 1);   % [N*scalXMult*scalYMult x 1]
    vecPosYFinal    = zeros(scalNTiled * scalYMult, 1);   % [N*scalXMult*scalYMult x 1]
    vecDiameterFinal = zeros(scalNTiled * scalYMult, 1);  % [N*scalXMult*scalYMult x 1]

    for ii = 0:scalYMult-1
        idxStart = ii*scalNTiled + 1;
        idxEnd   = ii*scalNTiled + scalNTiled;
        vecPosXFinal(idxStart:idxEnd)     = vecPosXTiled;                       % x unchanged
        vecPosYFinal(idxStart:idxEnd)     = vecPosYTiled + ii * scalBoxHeightY;  % shift y only
        vecDiameterFinal(idxStart:idxEnd) = vecDiameterTiled;
    end

    scalBoxHeightYFinal = scalBoxHeightY * scalYMult;
    scalNFinal         = scalNTiled * scalYMult;
    scalWidthFactorFinal = scalWidthFactor * scalYMult;
    scalDiameterAverage  = mean(vecDiameterFinal);
    scalMass             = 1;

    %% Save
    strFilenameOut = sprintf('%s2D_N%d_P%s_Width%d_Seed%d.mat', ...
        strSavePath, scalNFinal, num2str(P_target), scalWidthFactorFinal, seed);

    if calc_eig
        matPositions = [vecPosXFinal, vecPosYFinal];              % [scalNFinal x 2]
        vecRadii     = vecDiameterFinal ./ 2;                     % [scalNFinal x 1]
        [matPositions, vecRadii] = cleanRats(matPositions, vecRadii, K, scalBoxHeightYFinal, scalBoxWidthXTiled);
        matHessian = hess2d(matPositions, vecRadii, K, scalBoxHeightYFinal, scalBoxWidthXTiled);  % [2*scalNFinal x 2*scalNFinal]
        [matEigenVectors, matEigenValues] = eig(matHessian);
        save(strFilenameOut, 'vecPosXFinal', 'vecPosYFinal', 'vecDiameterFinal', ...
            'scalBoxWidthXTiled', 'scalBoxHeightYFinal', 'K', 'P_target', 'scalPressure', ...
            'scalNFinal', 'matEigenVectors', 'matEigenValues');
    else
        % Standardize variable names to match pack.m output
        vecPosX        = vecPosXFinal;
        vecPosY        = vecPosYFinal;
        vecDiameter    = vecDiameterFinal;
        scalBoxWidthX  = scalBoxWidthXTiled;
        scalBoxHeightY = scalBoxHeightYFinal;
        N              = scalNFinal;
        N_original     = N;  % no rattler removal in tiling

        % Packing fraction of the tiled packing
        scalVolumeSpheres_clean = sum(pi * (vecDiameter/2).^2);
        scalVolumeBox_clean     = scalBoxWidthX * scalBoxHeightY;
        scalPackingFraction     = scalVolumeSpheres_clean / scalVolumeBox_clean;
        scalPackingFractionFull = scalPackingFraction;  % same for tiling

        % Coordination number not recomputed here — mark as NaN
        scalMeanCoordNum = NaN;

        % Friction flags: default if not loaded from the base tile
        if ~exist('boolFrictionOn','var')
            boolFrictionOn = false;
        end
        if ~exist('scalMu','var')
            scalMu = 0.0;
        end
        if ~exist('scalKt','var')
            scalKt = 0.0;
        end

        save(strFilenameOut, 'vecPosX', 'vecPosY', 'vecDiameter', ...
            'scalBoxWidthX', 'scalBoxHeightY', 'K', 'P_target', 'scalPressure', ...
            'N', 'N_original', 'scalPackingFraction', 'scalPackingFractionFull', ...
            'scalMeanCoordNum', 'boolFrictionOn', 'scalMu', 'scalKt');
    end
        %% Plot tiled packing
        figure;
        hold on;
        for np = 1:scalNFinal
            rectangle('Position', [vecPosXFinal(np) - vecDiameterFinal(np)/2, ...
                                    vecPosYFinal(np) - vecDiameterFinal(np)/2, ...
                                    vecDiameterFinal(np), vecDiameterFinal(np)], ...
                      'Curvature', [1 1], 'FaceColor', 'b', 'EdgeColor', 'none');
        end

        % Ghost tiles around the main box (same style as pack.m)
        vecOffsets2D = [scalBoxWidthXTiled,  0; ...
                       -scalBoxWidthXTiled,  0; ...
                        0,  scalBoxHeightYFinal; ...
                        0, -scalBoxHeightYFinal; ...
                        scalBoxWidthXTiled,  scalBoxHeightYFinal; ...
                       -scalBoxWidthXTiled,  scalBoxHeightYFinal; ...
                        scalBoxWidthXTiled, -scalBoxHeightYFinal; ...
                       -scalBoxWidthXTiled, -scalBoxHeightYFinal];

        for iface = 1:8
            ox = vecOffsets2D(iface, 1);
            oy = vecOffsets2D(iface, 2);
            for np = 1:scalNFinal
                rectangle('Position', [vecPosXFinal(np) + ox - vecDiameterFinal(np)/2, ...
                                        vecPosYFinal(np) + oy - vecDiameterFinal(np)/2, ...
                                        vecDiameterFinal(np), vecDiameterFinal(np)], ...
                          'Curvature', [1 1], 'FaceColor', 'r', 'EdgeColor', 'none', ...
                          'FaceAlpha', 0.15);
            end
        end

        axis equal;
        axis([-scalBoxWidthXTiled  2*scalBoxWidthXTiled  ...
              -scalBoxHeightYFinal 2*scalBoxHeightYFinal]);
        title(sprintf('Tiled packing: N=%d, Xmult=%d, Ymult=%d', ...
            scalNFinal, scalXMult, scalYMult));
        drawnow;
        hold off;

    disp("Saved to: " + strFilenameOut);

end

