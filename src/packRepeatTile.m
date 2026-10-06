function packRepeatTile(N, K, P_target, scalWidthFactor, seed, scalXMult, scalYMult, calc_eig, strInPath, strSavePath, scalZMult, boolHertzian)
    % packRepeatTile  Load a saved 2D or 3D packing and tile it in x, y, and/or z.
    %
    % N,              Number of particles in the base tile
    % K,              Spring constant
    % P_target,       Target pressure
    % scalWidthFactor, Width factor used in the base tile filename
    % seed,           RNG seed used in the base tile filename
    % scalXMult,      Number of times to tile in x
    % scalYMult,      Number of times to tile in y
    % calc_eig,       Boolean: compute and save eigenmodes (2D only)
    % strInPath,      Path to the base tile .mat file
    % strSavePath,    Path to save the tiled output .mat file
    % scalZMult,      (optional, default 0) Number of times to tile in z;
    %                 ~= 0 selects 3D (same convention as pack.m)
    % boolHertzian,   (optional, default false) 3D only: base tile and
    %                 output use the '_Hertz' filename suffix
    %
    % 2D: loads the backbone 2D_N<N>_P<P>_Width<W>_Seed<s>.mat and saves
    %     2D_N<Ntiled>_P<P>_Width<W*scalYMult>_Seed<s>.mat.
    % 3D: loads the FULL (pre-cleanRats, rattlers included) tile source
    %     3D_N<N>_P<P>_Width<W>_Seed<s>_Full[_Hertz].mat that pack.m writes
    %     when options.saveFullState = true, and saves the superlattice under
    %     the normal packing name 3D_N<Ntiled>_P<P>_Width<W>_Seed<s>[_Hertz].mat
    %     (multipliers stored inside the file as x_mult, y_mult, z_mult).
    %
    % Example:
    % packRepeatTile(400, 100, 0.05, 22, 1, 2, 2, false, 'in/tiles/', 'out/tiles/')
    % packRepeatTile(1000, 100, 0.05, 10, 1, 1, 1, false, 'in/', 'out/', 9)   % 3D, 9x in z

    if nargin < 11 || isempty(scalZMult)
        scalZMult = 0;
    end
    if nargin < 12 || isempty(boolHertzian)
        boolHertzian = false;
    end

    if scalZMult ~= 0
        packRepeatTile3D(N, K, P_target, scalWidthFactor, seed, scalXMult, scalYMult, scalZMult, ...
            calc_eig, strInPath, strSavePath, boolHertzian);
        return
    end

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

    %% Tile in x and y (shared tiler, z_mult = 1 for 2D)
    [vecPosXFinal, vecPosYFinal, ~, vecDiameterFinal, ...
        scalBoxWidthXTiled, scalBoxHeightYFinal, ~, scalNFinal] = tileGrid( ...
        vecPosX, vecPosY, zeros(size(vecPosX)), vecDiameter, ...
        scalBoxWidthX, scalBoxHeightY, 0, scalXMult, scalYMult, 1);

    scalWidthFactorFinal = scalWidthFactor * scalYMult;
    scalDiameterAverage  = mean(vecDiameterFinal);
    scalMass             = 1;

    %% Save
    strFilenameOut = sprintf('%s2D_N%d_P%s_Width%d_Seed%d.mat', ...
        strSavePath, scalNFinal, num2str(P_target), scalWidthFactorFinal, seed);

    if calc_eig
        % Eigen path: remove rattlers first, then build the Hessian on the
        % backbone (pack.m does the same before saving eigenmodes).
        matPositions = [vecPosXFinal, vecPosYFinal];              % [scalNFinal x 2]
        vecRadii     = vecDiameterFinal ./ 2;                     % [scalNFinal x 1]
        N_original   = scalNFinal;                                % before cleanRats
        [matPositions, vecRadii] = cleanRats(matPositions, vecRadii, K, scalBoxHeightYFinal, scalBoxWidthXTiled);
        vecPosX        = matPositions(:,1);
        vecPosY        = matPositions(:,2);
        vecDiameter    = vecRadii .* 2;
        N              = size(matPositions, 1);
        matHessian = hess2d(matPositions, vecRadii, K, scalBoxHeightYFinal, scalBoxWidthXTiled);  % [2*N x 2*N]
        [matEigenVectors, matEigenValues] = eig(matHessian);

        % Recompute the saved metrics on the backbone, exactly like pack.m
        % (shared functions: single source of truth for the definitions)
        scalBoxWidthX  = scalBoxWidthXTiled;
        scalBoxHeightY = scalBoxHeightYFinal;
        scalPackingFraction = computePackingFraction(vecDiameter, scalBoxWidthX, scalBoxHeightY);
        scalPackingFractionFull = scalPackingFraction;  % same for tiling
        scalMeanCoordNum = computeMeanCoordNum(vecPosX, vecPosY, vecDiameter, scalBoxWidthX, scalBoxHeightY);
        fprintf('Tiled backbone: N=%d, PF=%.4f, mean coordination number=%.4f\n', ...
            N, scalPackingFraction, scalMeanCoordNum);

        save(strFilenameOut, 'vecPosX', 'vecPosY', 'vecDiameter', ...
            'scalBoxWidthX', 'scalBoxHeightY', 'K', 'P_target', 'scalPressure', ...
            'N', 'N_original', 'scalPackingFraction', 'scalPackingFractionFull', ...
            'scalMeanCoordNum', 'matEigenVectors', 'matEigenValues');
    else
        % Standardize variable names to match pack.m output
        vecPosX        = vecPosXFinal;
        vecPosY        = vecPosYFinal;
        vecDiameter    = vecDiameterFinal;
        scalBoxWidthX  = scalBoxWidthXTiled;
        scalBoxHeightY = scalBoxHeightYFinal;
        N              = scalNFinal;
        N_original     = N;  % no rattler removal in tiling

        % Packing fraction of the tiled packing, via the shared function
        % (identical to the base tile by construction: tiling scales box
        % area and disk area by the same factor)
        scalPackingFraction = computePackingFraction(vecDiameter, scalBoxWidthX, scalBoxHeightY);
        scalPackingFractionFull = scalPackingFraction;  % same for tiling

        % Coordination number recomputed on the tiled packing under full PBC
        % via the shared function (same definition as pack.m). Tiling
        % replicates the base contact network, so for a periodic base
        % packing this equals the base tile's scalMeanCoordNum.
        scalMeanCoordNum = computeMeanCoordNum(vecPosX, vecPosY, vecDiameter, scalBoxWidthX, scalBoxHeightY);
        fprintf('Tiled packing: N=%d, PF=%.4f, mean coordination number=%.4f\n', ...
            N, scalPackingFraction, scalMeanCoordNum);

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

function packRepeatTile3D(N, K, P_target, scalWidthFactor, seed, scalXMult, scalYMult, scalZMult, ...
        calc_eig, strInPath, strSavePath, boolHertzian)
    % packRepeatTile3D  3D branch of packRepeatTile: tile the full
    % (pre-cleanRats) 3D tile source in x, y, and z and save the superlattice.

    if scalXMult < 1 || scalYMult < 1 || scalZMult < 1
        error('packRepeatTile:BadMult', ...
            'packRepeatTile: 3D multipliers must be >= 1 (got x=%g, y=%g, z=%g).', ...
            scalXMult, scalYMult, scalZMult);
    end
    if calc_eig
        warning('calc_eig not supported for 3D yet — saving positions only.');
    end
    if boolHertzian
        strSuffix = '_Hertz';
    else
        strSuffix = '';
    end

    %% Load the full-state base tile (rattlers included)
    strFilenameIn = sprintf('%s3D_N%d_P%s_Width%d_Seed%d_Full%s.mat', ...
        strInPath, N, num2str(P_target), scalWidthFactor, seed, strSuffix);
    fprintf('Loading file: %s\n', strFilenameIn);
    if ~isfile(strFilenameIn)
        error('packRepeatTile:MissingFullState', ...
            ['3D repeat-tile requires the full (pre-cleanRats) tile source %s. ', ...
             'Run pack(..., options) with options.saveFullState = true.'], strFilenameIn);
    end
    sBase = load(strFilenameIn);
    scalNTileSrc = numel(sBase.vecPosX);   % full per-tile particle count

    %% Tile in x, y, and z
    [vecPosXFinal, vecPosYFinal, vecPosZFinal, vecDiameterFinal, ...
        scalBoxWidthXTiled, scalBoxHeightYFinal, scalBoxDepthZFinal, scalNTiled] = tileGrid( ...
        sBase.vecPosX, sBase.vecPosY, sBase.vecPosZ, sBase.vecDiameter, ...
        sBase.scalBoxWidthX, sBase.scalBoxHeightY, sBase.scalBoxDepthZ, ...
        scalXMult, scalYMult, scalZMult);

    % Metrics on the tiled packing (rattlers included, matching the stored
    % particles). Packing fraction is identical to the base tile by
    % construction; coordination number under full PBC equals the base
    % tile's full-state Zn (tiling replicates the contact network).
    scalPackingFractionTiled     = computePackingFraction(vecDiameterFinal, scalBoxWidthXTiled, scalBoxHeightYFinal, scalBoxDepthZFinal);
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
        scalMeanCoordNumTiled = computeMeanCoordNum(sBase.vecPosX, sBase.vecPosY, sBase.vecDiameter, ...
            sBase.scalBoxWidthX, sBase.scalBoxHeightY, sBase.vecPosZ, sBase.scalBoxDepthZ);
    end
    fprintf('3D tiled packing: N=%d (of %d per tile), PF=%.4f, mean coordination number=%.4f\n', ...
        scalNTiled, scalNTileSrc, scalPackingFractionTiled, scalMeanCoordNumTiled);

    % Tiled output filename: the normal packing name with N = the TILED
    % particle count, so simMD finds it like any other packing. Width and
    % Seed are the base tile's; the multipliers are stored inside the file.
    strFilenameOut = sprintf('%s3D_N%d_P%s_Width%d_Seed%d%s.mat', ...
        strSavePath, scalNTiled, num2str(P_target), scalWidthFactor, seed, strSuffix);

    % Friction flags: default if not stored in the base tile
    boolFrictionOn = false;
    scalMu = 0.0;
    scalKt = 0.0;
    if isfield(sBase, 'boolFrictionOn'), boolFrictionOn = sBase.boolFrictionOn; end
    if isfield(sBase, 'scalMu'),         scalMu = sBase.scalMu;                 end
    if isfield(sBase, 'scalKt'),         scalKt = sBase.scalKt;                 end

    % Stored variable names match the backbone .mat convention (vecPosX, N,
    % scalPackingFraction, ...) so downstream loaders are unchanged; save
    % -struct writes one variable per field, in the explicit field order.
    sTiledFile = struct();
    sTiledFile.vecPosX                 = vecPosXFinal;
    sTiledFile.vecPosY                 = vecPosYFinal;
    sTiledFile.vecPosZ                 = vecPosZFinal;
    sTiledFile.vecDiameter             = vecDiameterFinal;
    sTiledFile.scalBoxWidthX           = scalBoxWidthXTiled;
    sTiledFile.scalBoxHeightY          = scalBoxHeightYFinal;
    sTiledFile.scalBoxDepthZ           = scalBoxDepthZFinal;
    sTiledFile.K                       = K;
    sTiledFile.P_target                = P_target;
    sTiledFile.scalPressure            = sBase.scalPressure;
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
    sTiledFile.NTileSrc                = scalNTileSrc;
    sTiledFile.scalPackingFractionTiled     = scalPackingFractionTiled;
    sTiledFile.scalPackingFractionFullTiled = scalPackingFractionFullTiled;
    sTiledFile.scalMeanCoordNumTiled        = scalMeanCoordNumTiled;
    cellTiledFields = fieldnames(sTiledFile);
    save(strFilenameOut, '-struct', 'sTiledFile', cellTiledFields{:});
    fprintf('3D tiled packing saved to: %s\n', strFilenameOut);
end

function [vecPosXFinal, vecPosYFinal, vecPosZFinal, vecDiameterFinal, ...
          scalBoxWidthXTiled, scalBoxHeightYFinal, scalBoxDepthZFinal, ...
          scalNTiled] = tileGrid(vecPosX, vecPosY, vecPosZ, vecDiameter, ...
          scalBoxWidthX, scalBoxHeightY, scalBoxDepthZ, x_mult, y_mult, z_mult)
    % tileGrid  Repeat a packing x_mult/y_mult/z_mult times along x/y/z.
    %
    % Shared by the 2D (z_mult = 1, vecPosZ = 0) and 3D paths. The base-tile
    % positions must already be wrapped into [0, L) on every axis (pack()
    % wraps positions each step), so the shifted copies sit in
    % [i*L, (i+1)*L) and the periodic structure is reproduced exactly across
    % tile boundaries. Copies are appended x fastest, then y, then z.
    N = numel(vecPosX);
    if numel(vecPosY) ~= N || numel(vecPosZ) ~= N || numel(vecDiameter) ~= N
        error('packRepeatTile:SizeMismatch', ...
            ['packRepeatTile: position/diameter vectors must all have N=%d elements; ', ...
             'got y=%d, z=%d, D=%d.'], N, numel(vecPosY), numel(vecPosZ), numel(vecDiameter));
    end
    scalNTiled = N * x_mult * y_mult * z_mult;

    vecPosXFinal     = zeros(scalNTiled, 1);
    vecPosYFinal     = zeros(scalNTiled, 1);
    vecPosZFinal     = zeros(scalNTiled, 1);
    vecDiameterFinal = zeros(scalNTiled, 1);

    idxEnd = 0;
    for iz = 0:z_mult-1
        for iy = 0:y_mult-1
            for ix = 0:x_mult-1
                idxStart = idxEnd + 1;
                idxEnd   = idxStart + N - 1;
                vecPosXFinal(idxStart:idxEnd)     = vecPosX(:) + ix * scalBoxWidthX;
                vecPosYFinal(idxStart:idxEnd)     = vecPosY(:) + iy * scalBoxHeightY;
                vecPosZFinal(idxStart:idxEnd)     = vecPosZ(:) + iz * scalBoxDepthZ;
                vecDiameterFinal(idxStart:idxEnd) = vecDiameter(:);
            end
        end
    end

    scalBoxWidthXTiled  = scalBoxWidthX * x_mult;
    scalBoxHeightYFinal = scalBoxHeightY * y_mult;
    scalBoxDepthZFinal  = scalBoxDepthZ * z_mult;
end
