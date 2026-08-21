function testEigenpair(strLoadPath, options)
%TESTEIGENPAIRCONSISTENCY  Check if eigenpair matches integrator model.
%
% Usage (old eigenData file):
%   testEigenpairConsistency("...results_2D_iso_N100_P0.1_Seed5_gamma_3.00e-03.mat", ...
%                            struct('loadTestData', false, 'modeIndex', 1));
%
% Usage (NEW.mat test data):
%   testEigenpairConsistency("...NEW.mat", ...
%                            struct('loadTestData', true, 'modeIndex', 1));


        if nargin < 2
            options = struct();
        end
        if ~isfield(options, 'loadTestData')
            options.loadTestData = false;
        end
        if ~isfield(options, 'modeIndex')
            options.modeIndex = 1;
        end
        if ~isfield(options, 'scalEpsilon')
            options.scalEpsilon = 1e-6;
        end

    % ---- Load data, mirroring simEigenmode branches ---------------------
    if options.loadTestData
        % NEW.mat format, from outData
        s = load(strLoadPath);

        outData = s.outData;

        scalSpringConst = outData.springConstant;
        scalPressure    = outData.pressure;
        scalBv          = outData.damping;
        scalLx          = outData.Lx;
        scalLy          = outData.Ly;

        matEigenvectors = outData.eigenVectors{1};   % (2N) x Nmodes
        vecEigenvalues  = outData.eigenValues{1};    % Nmodes x 1

        vecRadii    = outData.radii{1};              % N x 1
        matPos      = outData.positions{1};          % N x 2
        vecX        = matPos(:,1);
        vecY        = matPos(:,2);
        vecDn       = 2*vecRadii;
        scalNumPart = numel(vecRadii);

        fprintf('[testEigenpairConsistency] NEW.mat format detected.\n');

    else
        % Old eigenData format, using same fields as loadSimParams
        s = load(strLoadPath);

        vecEigenvalues  = s.eigenvalues;
        matEigenvectors = s.eigenvectors;
        scalGamma       = s.gamma;
        sPacking        = s.packing;

        scalSpringConst = sPacking.K;
        scalNumPart     = sPacking.N;
        scalPressure    = sPacking.P_target;
        scalLx          = sPacking.Lx;
        scalLy          = sPacking.Ly;
        scalBv          = scalGamma;

        vecX  = sPacking.x(:);
        vecY  = sPacking.y(:);
        vecDn = sPacking.Dn(:);

        fprintf('[testEigenpairConsistency] Old results_2D_*.mat format detected.\n');

        % If eigenvectors cover fewer DOFs than packing, truncate packing
        scalNdofMode = size(matEigenvectors, 1);
        if mod(scalNdofMode, 2) ~= 0
            error('testEigenpairConsistency:ModeDOFParity', ...
                  'Mode DOF count (%d) is not even; expected 2*NumParticles.', scalNdofMode);
        end
        scalNumPartDyn = scalNdofMode / 2;

        if scalNumPartDyn ~= scalNumPart
            fprintf(['[testEigenpairConsistency] WARNING: packing has scalNumPart=%d, ', ...
                     'but eigenvectors describe scalNumPartDyn=%d particles. ', ...
                     'Using scalNumPartDyn (truncating packing).\n'], ...
                    scalNumPart, scalNumPartDyn);
            scalNumPart = scalNumPartDyn;
            vecX        = vecX(1:scalNumPart);
            vecY        = vecY(1:scalNumPart);
            vecDn       = vecDn(1:scalNumPart);
        end
    end

    % ---- Basic derived quantities --------------------------------------
    vecMass    = (pi/4) .* vecDn.^2;
    vecInvMass = 1 ./ vecMass;

    vecLambda = vecEigenvalues(:);
    if options.modeIndex < 1 || options.modeIndex > numel(vecLambda)
        error('testEigenpairConsistency:ModeIndex', ...
            'modeIndex=%d out of range (1..%d).', options.modeIndex, numel(vecLambda));
    end

    scalLambda = vecLambda(options.modeIndex);
    scalOmega  = imag(scalLambda);
    scalBeta   = real(scalLambda);

    fprintf('[testEigenpairConsistency] Testing mode %d: lambda = %g + %gi\n', ...
        options.modeIndex, scalBeta, scalOmega);

    % ---- Extract eigenvector and mode displacement ----------------------
    vecModeRawReal = real(matEigenvectors(:, options.modeIndex));
    scalNdof       = numel(vecModeRawReal);

    if scalNdof < 2*scalNumPart
        error('testEigenpairConsistency:ModeSize', ...
              'Mode DOF count (%d) < 2*scalNumPart (%d).', scalNdof, 2*scalNumPart);
    end

    vecMode     = vecModeRawReal(1:2*scalNumPart);
    matModeDisp = reshape(vecMode, 2, scalNumPart).';    % N x 2: [dx_i, dy_i]

    % ---- Build neighbor list with periodic boundaries -------------------
    fprintf('[testEigenpairConsistency] Building neighbor list ...\n');
    scalSkin            = 0;
    vecZnList           = zeros(scalNumPart, 1);
    cellNeighborListAll = cell(1, scalNumPart);
    cellSpringListAll   = cell(1, scalNumPart);

    for i = 1:scalNumPart
        vecNeighborList = [];
        vecSpringList   = [];
        for j = [1:i-1, i+1:scalNumPart]
            dy = vecY(j) - vecY(i);
            dy = dy - round(dy/scalLy)*scalLy;     % PBC y

            dx = vecX(j) - vecX(i);
            dx = dx - round(dx/scalLx)*scalLx;     % PBC x

            dnm = (1+scalSkin)*(vecDn(i)+vecDn(j))/2;
            dsq = dx^2 + dy^2;

            if dsq < dnm^2
                vecNeighborList = [vecNeighborList, j];    %#ok<AGROW>
                vecSpringList   = [vecSpringList, sqrt(dsq)];  %#ok<AGROW>
            end
        end
        cellNeighborListAll{i} = vecNeighborList;
        cellSpringListAll{i}   = vecSpringList;
        vecZnList(i)           = length(vecSpringList);
    end
    fprintf('[testEigenpairConsistency] Neighbor list built.\n');

    % Edge list (same as simEigenmode)
    scalMaxNb  = max(vecZnList);
    matNbrIdx  = zeros(scalNumPart, scalMaxNb, 'int32');
    matSprLen  = zeros(scalNumPart, scalMaxNb);

    for i = 1:scalNumPart
        vecNbr = cellNeighborListAll{i};
        vecSpr = cellSpringListAll{i};
        nb     = length(vecNbr);
        matNbrIdx(i, 1:nb) = int32(vecNbr);
        matSprLen(i, 1:nb) = vecSpr;
    end
    clear cellNeighborListAll cellSpringListAll;

    matValid = (matNbrIdx > 0);
    matSrc   = repmat((1:scalNumPart)', 1, scalMaxNb);

    vecSrcFlat = matSrc(matValid);            % E×1
    vecDstFlat = double(matNbrIdx(matValid)); % E×1
    vecD0Flat  = matSprLen(matValid);         % E×1

    fprintf('[testEigenpairConsistency] Edge list has %d directed edges.\n', ...
        length(vecSrcFlat));

    % ---- Apply small displacement along eigenvector ---------------------
    eps  = options.scalEpsilon;
    vecX0 = vecX;
    vecY0 = vecY;

    vecX = vecX0 + eps * matModeDisp(:, 1);
    vecY = vecY0 + eps * matModeDisp(:, 2);

    % Enforce PBC
    vecX = mod(vecX, scalLx);
    vecY = mod(vecY, scalLy);

    % Zero velocities (so damping not active in this test)
    vecVx = zeros(scalNumPart, 1);
    vecVy = zeros(scalNumPart, 1);

    % ---- Compute forces using same kernel as simEigenmode --------------
    K  = scalSpringConst;
    Bv = scalBv;

    % Source/dest positions and velocities
    xs  = vecX(vecSrcFlat);    ys  = vecY(vecSrcFlat);
    xd  = vecX(vecDstFlat);    yd  = vecY(vecDstFlat);
    vxs = vecVx(vecSrcFlat);   vys = vecVy(vecSrcFlat);
    vxd = vecVx(vecDstFlat);   vyd = vecVy(vecDstFlat);

    dx = xd - xs;
    dy = yd - ys;

    dx = dx - round(dx./scalLx).*scalLx;
    dy = dy - round(dy./scalLy).*scalLy;

    dnm  = sqrt(dx.^2 + dy.^2);
    Fmag = -K .* (vecD0Flat./dnm - 1);     % same formula as simEigenmode

    dvx = vxs - vxd;
    dvy = vys - vyd;

    FxEdge = Fmag.*dx - Bv.*dvx;
    FyEdge = Fmag.*dy - Bv.*dvy;

    Fx = accumarray(vecSrcFlat, FxEdge, [scalNumPart, 1]);
    Fy = accumarray(vecSrcFlat, FyEdge, [scalNumPart, 1]);

    % ---- Accelerations and eigenvalue ratio -----------------------------
    Ax = Fx .* vecInvMass;
    Ay = Fy .* vecInvMass;

    dispX = eps * matModeDisp(:, 1);
    dispY = eps * matModeDisp(:, 2);

    maskX = abs(dispX) > 1e-12;
    maskY = abs(dispY) > 1e-12;

    ratioX = Ax(maskX) ./ dispX(maskX);
    ratioY = Ay(maskY) ./ dispY(maskY);
    ratioAll = [ratioX; ratioY];

    scalRatioMean = mean(ratioAll);
    scalRatioStd  = std(ratioAll);

    fprintf('[testEigenpairConsistency] Mean(a/u) = %g, Std(a/u) = %g\n', ...
        scalRatioMean, scalRatioStd);

    % Expected value if eigenvalues are frequencies: a ≈ -omega^2 u
    scalOmegaSq = -scalOmega^2;

    fprintf('[testEigenpairConsistency] Expected a/u ~ %g (= -omega_eig^2)\n', ...
        scalOmegaSq);
    fprintf('[testEigenpairConsistency] Relative std = %g\n', ...
        scalRatioStd / abs(scalRatioMean));

    % Simple visual sanity check
    figure;
    subplot(2,1,1);
    plot(ratioX, 'o');
    ylabel('a_x / u_x');
    grid on;
    subplot(2,1,2);
    plot(ratioY, 'o');
    ylabel('a_y / u_y');
    xlabel('particle index');
    grid on;
    sgtitle(sprintf('Mode %d: lambda = %.4g + %.4gi', options.modeIndex, scalBeta, scalOmega));

end
