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
    strBaseName = sprintf('2D_N%d_P%s_Width%d_Seed%d', N, num2str(P_target), scalWidthFactor, seed);
    strFilenameIn = strInPath + strBaseName + ".mat";
    load(strFilenameIn, 'vecPosX', 'vecPosY', 'vecDiameter', ...
        'scalBoxWidth', 'scalBoxHeight', 'K', 'P_target', 'scalPressure', 'N');

    %% Tile in x
    % Shift copies of vecPosX by multiples of scalBoxWidth, keep vecPosY unchanged
    vecPosXTiled = zeros(N * scalXMult, 1);       % [N*scalXMult x 1]
    vecPosYTiled = zeros(N * scalXMult, 1);       % [N*scalXMult x 1]
    vecDiameterTiled = zeros(N * scalXMult, 1);   % [N*scalXMult x 1]

    for ii = 0:scalXMult-1
        idxStart = ii*N + 1;
        idxEnd   = ii*N + N;
        vecPosXTiled(idxStart:idxEnd)    = vecPosX + ii * scalBoxWidth;  % shift x only
        vecPosYTiled(idxStart:idxEnd)    = vecPosY;                      % y unchanged
        vecDiameterTiled(idxStart:idxEnd) = vecDiameter;
    end

    scalBoxWidthTiled = scalBoxWidth * scalXMult;
    scalNTiled        = N * scalXMult;

    %% Tile in y
    % Shift copies of vecPosYTiled by multiples of scalBoxHeight, keep x unchanged
    vecPosXFinal    = zeros(scalNTiled * scalYMult, 1);   % [N*scalXMult*scalYMult x 1]
    vecPosYFinal    = zeros(scalNTiled * scalYMult, 1);   % [N*scalXMult*scalYMult x 1]
    vecDiameterFinal = zeros(scalNTiled * scalYMult, 1);  % [N*scalXMult*scalYMult x 1]

    for ii = 0:scalYMult-1
        idxStart = ii*scalNTiled + 1;
        idxEnd   = ii*scalNTiled + scalNTiled;
        vecPosXFinal(idxStart:idxEnd)     = vecPosXTiled;                       % x unchanged
        vecPosYFinal(idxStart:idxEnd)     = vecPosYTiled + ii * scalBoxHeight;  % shift y only
        vecDiameterFinal(idxStart:idxEnd) = vecDiameterTiled;
    end

    scalBoxHeightFinal = scalBoxHeight * scalYMult;
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
        [matPositions, vecRadii] = cleanRats(matPositions, vecRadii, K, scalBoxHeightFinal, scalBoxWidthTiled);
        matHessian = hess2d(matPositions, vecRadii, K, scalBoxHeightFinal, scalBoxWidthTiled);  % [2*scalNFinal x 2*scalNFinal]
        [matEigenVectors, matEigenValues] = eig(matHessian);
        save(strFilenameOut, 'vecPosXFinal', 'vecPosYFinal', 'vecDiameterFinal', ...
            'scalBoxWidthTiled', 'scalBoxHeightFinal', 'K', 'P_target', 'scalPressure', ...
            'scalNFinal', 'matEigenVectors', 'matEigenValues');
    else
        save(strFilenameOut, 'vecPosXFinal', 'vecPosYFinal', 'vecDiameterFinal', ...
            'scalBoxWidthTiled', 'scalBoxHeightFinal', 'K', 'P_target', 'scalPressure', ...
            'scalNFinal', 'scalWidthFactorFinal', 'scalMass', 'scalDiameterAverage');
    end

    disp("Saved to: " + strFilenameOut);

end

