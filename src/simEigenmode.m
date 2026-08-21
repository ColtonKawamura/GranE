function simEigenmode(strLoadPath, options)

    arguments
        strLoadPath (1,1) string
        options.scalModeFreq   (1,1) double  = 0;   % target omega
        options.scalModeAtten  (1,1) double  = 0;   % target beta
        options.fullSpectrum   (1,1) logical = false;
        options.visSim         (1,1) logical = false;
        options.loadTestData    (1,1) logical = false;
    end

%% Load eigenmode data and packing
    strInPath = strLoadPath;

    if options.loadTestData
        % Expect a struct "outData" with the fields:
        % pressure, damping, eigenVectors, eigenValues,
        % Ly, Lx, radii, positions, springConstant
        load(strLoadPath)

        % Basic scalars
        scalSpringConst = outData.springConstant;
        scalPressure    = outData.pressure;
        scalBv          = outData.damping;
        scalLx          = outData.Lx;
        scalLy          = outData.Ly;
        scalWidth       = outData.Lx;


        % Eigenvalues / eigenvectors (stored in cells)
        matEigenvectors = outData.eigenVectors{1};   % (2N) x Nmodes
        vecEigenvalues  = outData.eigenValues{1};    % Nmodes x 1

        % Particle data
        vecRadii    = outData.radii{1};              % N x 1
        matPos      = outData.positions{1};          % N x 2
        vecX        = matPos(:,1);
        vecY        = matPos(:,2);
        vecDn       = 2*vecRadii;
        scalNumPart = numel(vecRadii);

        % Dummy / logging-only scalars
        scalMass        = 1.0;
        scalFreqDriving = 0.0;
        scalSeed        = 0;

        % Target pressure for logging
        scalPTarget = scalPressure;
    else
        [scalSpringConst, scalMass, scalFreqDriving, scalNumPart, ...
         scalPressure, scalWidth, scalSeed, scalBv, ...
         vecEigenvalues, matEigenvectors, scalGamma, sPacking, ...
         strPackingFile, strResultsFile] = loadSimParams(strInPath);

        % Mode DOFs determine how many particles actually have eigenvectors
        scalNdofMode = size(matEigenvectors, 1);
        if mod(scalNdofMode, 2) ~= 0
            error('simEigenmode:ModeDOFParity', ...
                  'Mode DOF count (%d) is not even; expected 2*NumParticles.', scalNdofMode);
        end
        scalNumPartDyn = scalNdofMode / 2  % number of particles with eigenvectors

        if scalNumPartDyn ~= scalNumPart
            fprintf(['[simEigenmode] WARNING: packing has scalNumPart=%d, ', ...
                     'but eigenvectors describe scalNumPartDyn=%d particles. ', ...
                     'Using scalNumPartDyn for dynamics (truncating packing).\n'], ...
                    scalNumPart, scalNumPartDyn);
            scalNumPart = scalNumPartDyn;  % override for simulation
        end


        % Packing fields (Hungarian) — truncate to dynamic particles
        vecX        = sPacking.x(:);
        vecY        = sPacking.y(:);
        vecDn       = sPacking.Dn(:);
        vecX        = vecX(1:scalNumPart);
        vecY        = vecY(1:scalNumPart);
        vecDn       = vecDn(1:scalNumPart);
        scalLx      = sPacking.Lx;
        scalLy      = sPacking.Ly;
        scalPTarget = sPacking.P_target;
    end


%% Logging
    fprintf(['[simEigenmode] START  scalSpringConst=%g scalMass=%g ', ...
             'scalBv=%g scalFreqDriving=%g scalNumPart=%d scalPressure=%g ', ...
             'scalWidth=%g scalSeed=%d\n'], ...
        scalSpringConst, scalMass, scalBv, scalFreqDriving, ...
        scalNumPart, scalPressure, scalWidth, scalSeed);
    fprintf('[simEigenmode]   strInPath  = %s\n',  char(strInPath));
    scalTStart = tic;

%% Mode selection: by omega (imag part) or attenuation beta (real part)
    vecLambda = vecEigenvalues(:);  % eigenvalues as column

    if options.scalModeFreq ~= 0
        % Select by imaginary part (positive omega)
        vecOmega        = imag(vecLambda);
        vecMaskPositive = (vecOmega > 0);
        if ~any(vecMaskPositive)
            error('simEigenmode:NoPositiveFreq', ...
                  'No eigenvalues with positive imaginary part found.');
        end
        vecOmegaSel = vecOmega(vecMaskPositive);
        vecIdxAll   = find(vecMaskPositive);
        [~, scalIdxLocal] = min(abs(vecOmegaSel - options.scalModeFreq));
        scalIdxMode = vecIdxAll(scalIdxLocal);
        fprintf(['[simEigenmode] Mode selected by omega: ', ...
                 'target omega=%g, lambda(idx)=%g + %gi (idx=%d)\n'], ...
            options.scalModeFreq, real(vecLambda(scalIdxMode)), ...
            imag(vecLambda(scalIdxMode)), scalIdxMode);

    elseif options.scalModeAtten ~= 0
        % Select by real part (attenuation beta)
        vecBetaReal = real(vecLambda);
        [~, scalIdxMode] = min(abs(vecBetaReal - options.scalModeAtten));
        fprintf(['[simEigenmode] Mode selected by attenuation: ', ...
                 'target beta=%g, lambda(idx)=%g + %gi (idx=%d)\n'], ...
            options.scalModeAtten, real(vecLambda(scalIdxMode)), ...
            imag(vecLambda(scalIdxMode)), scalIdxMode);

    else
        error('simEigenmode:ModeSelection', ...
              'Either options.scalModeFreq or options.scalModeAtten must be non-zero.');
    end

    scalLambda    = vecLambda(scalIdxMode);
    scalOmegaMode = imag(scalLambda);  % omega
    scalBetaMode  = real(scalLambda);  % attenuation beta

    % Mode vector → per-particle dx,dy with Hungarian naming
    vecModeRawReal = real(matEigenvectors(:, scalIdxMode));  % column, real part only
    scalNdof       = numel(vecModeRawReal);

    if scalNdof < 2*scalNumPart
        error('simEigenmode:ModeSize', ...
              'Mode DOF count (%d) < 2*scalNumPart (%d).', scalNdof, 2*scalNumPart);
    elseif scalNdof > 2*scalNumPart
        fprintf(['[simEigenmode] WARNING: mode has %d DOFs, using first %d ', ...
                 'entries as particle displacements.\n'], ...
            scalNdof, 2*scalNumPart);
    end

    vecMode     = vecModeRawReal(1:2*scalNumPart);
    matModeDisp = reshape(vecMode, 2, scalNumPart).';    % scalNumPart x 2: [dx_i, dy_i]

try

%% GPU setup
    logUseGPU = canUseGPU();
    if logUseGPU
        sGpuDev = gpuDevice();
        fprintf('[simEigenmode] GPU detected: %s (%.1f GB)\n', ...
            sGpuDev.Name, sGpuDev.TotalMemory/1024^3);
    else
        fprintf('[simEigenmode] No GPU found — running on CPU.\n');
    end

%% Enforce column vectors and masses
    vecX  = vecX(:);
    vecY  = vecY(:);
    vecDn = vecDn(:);

    vecMass                 = (pi/4) .* vecDn.^2;
    vecInvMass              = 1 ./ vecMass;
    scalMassParticleAverage = mean(vecMass);
    scalMassMin             = min(vecMass);

%% Simulation parameters
    scalB          = 0;                      % global background damping
    scalAmpDriving = scalPTarget / 100;      % unused in eigenmode test, kept for logging
    scalDt         = pi*sqrt(scalMassMin/scalSpringConst)*0.005;
    scalC0         = min(vecDn)*sqrt(scalSpringConst/scalMassMin);
    scalNt         = round(1000*(scalLx/scalC0)/scalDt);

    fprintf('[simEigenmode] scalNt = %d, scalDt = %g\n', scalNt, scalDt);

    scalDt2Half = scalDt^2 / 2;
    scalDtHalf  = scalDt / 2;

%% Memory estimate
    scalTrajBytes = double(scalNumPart) * double(scalNt) * 8;
    fprintf(['[simEigenmode] Avoided trajectory arrays: %.2f GB each, ', ...
             '%.2f GB total.\n'], ...
        scalTrajBytes/1024^3, 2*scalTrajBytes/1024^3);

%% Wall / bulk masks — fully periodic system: no walls
    vecLeftWallMask  = false(scalNumPart, 1);
    vecRightWallMask = false(scalNumPart, 1);
    vecBulkMask      = true(scalNumPart, 1);  %#ok<NASGU>

    vecIdxLeftWall  = find(vecLeftWallMask);
    vecIdxRightWall = find(vecRightWallMask);

%% Build neighbor lists (CPU, done once) with periodic x,y
    fprintf('[simEigenmode] Building neighbor lists ...\n');
    scalSkin            = 0;
    vecZnList           = zeros(scalNumPart, 1);
    cellNeighborListAll = cell(1, scalNumPart);
    cellSpringListAll   = cell(1, scalNumPart);

    for scalIdxPart = 1:scalNumPart
        vecNeighborList = [];
        vecSpringList   = [];
        for scalIdxNeighbor = [1:scalIdxPart-1, scalIdxPart+1:scalNumPart]
            scalDy = vecY(scalIdxNeighbor) - vecY(scalIdxPart);
            scalDy = scalDy - round(scalDy/scalLy)*scalLy;     % PBC in y

            scalDx = vecX(scalIdxNeighbor) - vecX(scalIdxPart);
            scalDx = scalDx - round(scalDx/scalLx)*scalLx;     % PBC in x

            scalDnm = (1+scalSkin)*(vecDn(scalIdxPart)+vecDn(scalIdxNeighbor))/2;
            scalDsq = scalDx^2 + scalDy^2;

            if scalDsq < scalDnm^2
                vecNeighborList = [vecNeighborList, scalIdxNeighbor]; %#ok<AGROW>
                vecSpringList   = [vecSpringList,   sqrt(scalDsq)];   %#ok<AGROW>
            end
        end
        cellNeighborListAll{scalIdxPart} = vecNeighborList;
        cellSpringListAll{scalIdxPart}   = vecSpringList;
        vecZnList(scalIdxPart)           = length(vecSpringList);
    end
    fprintf('[simEigenmode] Neighbor lists built.\n');

%% Source and destination particle indices (edge list)
    scalMaxNb  = max(vecZnList);
    matNbrIdx  = zeros(scalNumPart, scalMaxNb, 'int32');
    matSprLen  = zeros(scalNumPart, scalMaxNb);

    for scalIdxPart = 1:scalNumPart
        vecNbr = cellNeighborListAll{scalIdxPart};
        vecSpr = cellSpringListAll{scalIdxPart};
        scalNb = length(vecNbr);
        matNbrIdx(scalIdxPart, 1:scalNb) = int32(vecNbr);
        matSprLen(scalIdxPart, 1:scalNb) = vecSpr;
    end
    clear cellNeighborListAll cellSpringListAll;

    matValid = (matNbrIdx > 0);
    matSrc   = repmat((1:scalNumPart)', 1, scalMaxNb);

    vecSrcFlat = matSrc(matValid);            % E×1
    vecDstFlat = double(matNbrIdx(matValid)); % E×1
    vecD0Flat  = matSprLen(matValid);         % E×1

    clear matNbrIdx matSprLen matValid matSrc;

    fprintf('[simEigenmode] Edge list: %d directed edges (avg %.1f neighbors/particle).\n', ...
        length(vecSrcFlat), length(vecSrcFlat)/scalNumPart);

%% Initial displacement along selected eigenmode (small amplitude)
    scalAmpInit = scalPTarget * 0.01;

    % Base equilibrium positions from packing
    vecX0 = vecX;
    vecY0 = vecY;

    % Apply displacement: vecX, vecY shifted by mode shape
    vecX = vecX0 + scalAmpInit * matModeDisp(:, 1);
    vecY = vecY0 + scalAmpInit * matModeDisp(:, 2);

    % Enforce periodic boundaries in BOTH directions
    vecX = mod(vecX, scalLx);
    vecY = mod(vecY, scalLy);

%% Verlet state arrays
    vecAxOld = zeros(scalNumPart, 1);
    vecAyOld = zeros(scalNumPart, 1);
    vecVx    = zeros(scalNumPart, 1);  % start at rest
    vecVy    = zeros(scalNumPart, 1);
    vecEk    = zeros(scalNt, 1);
    vecEp    = zeros(scalNt, 1);
    scalG    = 0;

%% Probe particle time series (for frequency & attenuation check)
    scalIdxProbe   = 11;                      % probe particle index
    vecProbeDispX  = zeros(scalNt, 1);
    vecProbeDispY  = zeros(scalNt, 1);
    vecProbeTime   = (0:scalNt-1).' * scalDt;

%% Full-spectrum trajectory storage (optional)
    if options.fullSpectrum
        scalTrajBytes = double(scalNumPart) * double(scalNt) * 8;
        fprintf(['[simEigenmode] fullSpectrum enabled — allocating %.2f GB ', ...
                 'per trajectory array.\n'], ...
            scalTrajBytes/1024^3);
        matTrajDispX = zeros(scalNumPart, scalNt);
        matTrajDispY = zeros(scalNumPart, scalNt);
    end

    % Index vector for visualizer (all particles)
    vecIdxAll = (1:scalNumPart)';

%% GPU transfer (everything hot goes to GPU once)
    if logUseGPU
        % Particle state
        vecX        = gpuArray(vecX);
        vecY        = gpuArray(vecY);
        vecX0_gpu   = gpuArray(vecX0);
        vecY0_gpu   = gpuArray(vecY0);
        vecVx       = gpuArray(vecVx);
        vecVy       = gpuArray(vecVy);
        vecAxOld    = gpuArray(vecAxOld);
        vecAyOld    = gpuArray(vecAyOld);
        vecEk       = gpuArray(vecEk);
        vecEp       = gpuArray(vecEp);
        vecInvMass_gpu = gpuArray(vecInvMass);
        vecMass_gpu    = gpuArray(vecMass);

        % Edge list
        vecSrcFlat = gpuArray(vecSrcFlat);
        vecDstFlat = gpuArray(vecDstFlat);
        vecD0Flat  = gpuArray(vecD0Flat);

        % Wall index lists (empty)
        vecIdxLeftWall_gpu  = gpuArray(vecIdxLeftWall);
        vecIdxRightWall_gpu = gpuArray(vecIdxRightWall);

        fprintf('[simEigenmode] Arrays transferred to GPU.\n');
    else
        % CPU aliases
        vecX0_gpu          = vecX0;
        vecY0_gpu          = vecY0;
        vecMass_gpu        = vecMass;
        vecInvMass_gpu     = vecInvMass;
        vecIdxLeftWall_gpu = vecIdxLeftWall;
        vecIdxRightWall_gpu= vecIdxRightWall;
    end

%% Main Verlet loop
    fprintf('[simEigenmode] Starting main loop (%d timesteps) ...\n', scalNt);
    scalLoopTic     = tic;
    scalLogInterval = max(1, floor(scalNt/10));

    for scalNtStep = 1:scalNt

        if mod(scalNtStep, scalLogInterval) == 0
            fprintf(['[simEigenmode]   timestep %d / %d  (%.1f%%, elapsed %.1f s)\n'], ...
                scalNtStep, scalNt, 100*scalNtStep/scalNt, toc(scalLoopTic));
        end

        % ── Verlet step 1: update positions ──────────────────────────────
        vecX = vecX + vecVx.*scalDt + vecAxOld.*scalDt2Half;
        vecY = vecY + vecVy.*scalDt + vecAyOld.*scalDt2Half;

        % Enforce periodic boundaries in BOTH directions
        vecX = mod(vecX, scalLx);
        vecY = mod(vecY, scalLy);

        % Optional visualization
        if options.visSim
            if logUseGPU
                vecX_vis = gather(vecX);
                vecY_vis = gather(vecY);
            else
                vecX_vis = vecX;
                vecY_vis = vecY;
            end
            visSim(vecX_vis, vecX0, vecY_vis, vecY0, ...
                   vecIdxAll, scalAmpDriving, scalBv, ...
                   scalPressure, scalFreqDriving, scalSpringConst, true);
        end

        % ── Vectorized force kernel with periodic boundaries ─────────────
        vecXs  = vecX(vecSrcFlat);    vecYs  = vecY(vecSrcFlat);
        vecXd  = vecX(vecDstFlat);    vecYd  = vecY(vecDstFlat);
        vecVxs = vecVx(vecSrcFlat);   vecVys = vecVy(vecSrcFlat);
        vecVxd = vecVx(vecDstFlat);   vecVyd = vecVy(vecDstFlat);

        vecDx = vecXd - vecXs;
        vecDy = vecYd - vecYs;

        vecDx = vecDx - round(vecDx./scalLx).*scalLx;
        vecDy = vecDy - round(vecDy./scalLy).*scalLy;

        vecDnm  = sqrt(vecDx.^2 + vecDy.^2);
        vecFmag = -scalSpringConst .* (vecD0Flat./vecDnm - 1);

        vecDvx = vecVxs - vecVxd;
        vecDvy = vecVys - vecVyd;

        vecFxEdge = vecFmag.*vecDx - scalBv.*vecDvx;
        vecFyEdge = vecFmag.*vecDy - scalBv.*vecDvy;
        vecEpEdge = 0.5*scalSpringConst.*(vecD0Flat - vecDnm).^2;

        vecFx = accumarray(vecSrcFlat, vecFxEdge, [scalNumPart, 1]);
        vecFy = accumarray(vecSrcFlat, vecFyEdge, [scalNumPart, 1]);

        % Background damping
        vecFx = vecFx - scalB.*vecVx;
        vecFy = vecFy - scalB.*vecVy;

        % Wall particles (none) — zero if any
        vecFx(vecIdxLeftWall_gpu)  = 0;
        vecFx(vecIdxRightWall_gpu) = 0;

        % Energies (Ep /2 because undirected pairs appear twice)
        vecEp(scalNtStep) = sum(vecEpEdge) / (2*scalNumPart);
        vecEk(scalNtStep) = (0.5 * sum(vecMass_gpu .* (vecVx.^2 + vecVy.^2))) / scalNumPart;

        % ── Accelerations ────────────────────────────────────────────────
        vecAx = vecFx .* vecInvMass_gpu;
        vecAy = vecFy .* vecInvMass_gpu - scalG;

        % ── Verlet step 2: update velocities ─────────────────────────────
        vecVx = vecVx + (vecAxOld + vecAx) .* scalDtHalf;
        vecVy = vecVy + (vecAyOld + vecAy) .* scalDtHalf;

        vecAxOld = vecAx;
        vecAyOld = vecAy;

        % Record probe displacement (gather from GPU if needed)
        if logUseGPU
            vecProbeDispX(scalNtStep) = gather(vecX(scalIdxProbe));
            vecProbeDispY(scalNtStep) = gather(vecY(scalIdxProbe));
        else
            vecProbeDispX(scalNtStep) = vecX(scalIdxProbe);
            vecProbeDispY(scalNtStep) = vecY(scalIdxProbe);
        end

        % Optional full-spectrum storage
        if options.fullSpectrum
            if logUseGPU
                matTrajDispX(:, scalNtStep) = gather(vecX);
                matTrajDispY(:, scalNtStep) = gather(vecY);
            else
                matTrajDispX(:, scalNtStep) = vecX;
                matTrajDispY(:, scalNtStep) = vecY;
            end
        end

    end

    fprintf('[simEigenmode] Main loop finished (%.1f s).\n', toc(scalLoopTic));

%% Gather GPU arrays back to CPU
    if logUseGPU
        vecX  = gather(vecX);
        vecY  = gather(vecY);
        vecX0 = gather(vecX0_gpu);
        vecY0 = gather(vecY0_gpu);
        vecVx = gather(vecVx);
        vecVy = gather(vecVy);
        vecEk = gather(vecEk);
        vecEp = gather(vecEp);
        fprintf('[simEigenmode] GPU arrays gathered to CPU.\n');
    end

%% Frequency and attenuation analysis (omega domain, displacement magnitude)
    vecProbeDispX_cpu = vecProbeDispX;
    vecProbeDispY_cpu = vecProbeDispY;
    vecProbeTime_cpu  = vecProbeTime;

    % Equilibrium position of probe particle
    scalX0Probe = vecX0(scalIdxProbe);
    scalY0Probe = vecY0(scalIdxProbe);

    % Displacement vector u(t) = r(t) - r0
    vecProbeUx = vecProbeDispX_cpu - scalX0Probe;
    vecProbeUy = vecProbeDispY_cpu - scalY0Probe;

    % Magnitude of displacement vector |u(t)|
    % vecProbeMag = sqrt(vecProbeUx.^2 + vecProbeUy.^2);
     vecProbeMag = vecProbeUx;

    scalNtFFT   = length(vecProbeMag);
    vecOmegaFFT = (0:scalNtFFT-1).' * (2*pi/(scalNtFFT*scalDt));
    vecFFT      = fft(vecProbeMag);
    vecAmp      = abs(vecFFT) / scalNtFFT;

    scalOmegaTol = 1e-3;
    scalAmpTol   = 1e-12;
    vecMaskGood  = (abs(vecOmegaFFT) > scalOmegaTol) & (vecAmp > scalAmpTol);

    if any(vecMaskGood)
        [scalAmpMax, idxLocal] = max(vecAmp(vecMaskGood));
        idxAll       = find(vecMaskGood);
        idxPeak      = idxAll(idxLocal);
        scalOmegaMax = vecOmegaFFT(idxPeak);
    else
        scalAmpMax   = 0;
        scalOmegaMax = NaN;
    end

    fprintf(['[simEigenmode] Mode data:          beta_dat = %g --- omega=%g\n'], ...
        real(scalLambda), imag(scalLambda), scalBetaMode, scalOmegaMode);

%% Attenuation measurement from |u(t)| = sqrt(ux^2 + uy^2)
    if abs(scalBetaMode) > 0
        vecMagPos    = vecProbeMag;
        vecMaskMag   = isfinite(vecMagPos) & (vecMagPos > 0);
        vecTFit      = vecProbeTime_cpu(vecMaskMag);
        vecLogMagFit = log(vecMagPos(vecMaskMag));

        if numel(vecTFit) > 10
            vecPolyCoeff = polyfit(vecTFit, vecLogMagFit, 1);
            scalBetaFit  = vecPolyCoeff(1);   % slope of ln|u| vs t


            fprintf(['[simEigenmode] Measured mode data: beta_fit = %g --- eigen omega = %g)\n'], ...
                scalBetaFit, scalOmegaMax);
        else
            scalBetaFit = NaN;
            fprintf('[simEigenmode] WARNING: insufficient displacement points for attenuation fit.\n');
        end
    else
        scalBetaFit = NaN;
        fprintf('[simEigenmode] Mode beta (real part of lambda) is zero — skipping attenuation fit.\n');
    end

%% Zero-crossing times of x and y (for oblique oscillations)
    % x(t) zero crossings
    vecSignX  = sign(vecProbeDispX_cpu);
    vecSignX(vecProbeDispX_cpu == 0) = 0;
    vecSignXProd = vecSignX(1:end-1) .* vecSignX(2:end);
    vecIdxZeroX  = find(vecSignXProd <= 0 & diff(vecProbeDispX_cpu) ~= 0);
    vecTZeroX    = vecProbeTime_cpu(vecIdxZeroX) ...
                 - vecProbeDispX_cpu(vecIdxZeroX) .* ...
                   (scalDt ./ (vecProbeDispX_cpu(vecIdxZeroX+1) - vecProbeDispX_cpu(vecIdxZeroX)));

    % y(t) zero crossings
    vecSignY  = sign(vecProbeDispY_cpu);
    vecSignY(vecProbeDispY_cpu == 0) = 0;
    vecSignYProd = vecSignY(1:end-1) .* vecSignY(2:end);
    vecIdxZeroY  = find(vecSignYProd <= 0 & diff(vecProbeDispY_cpu) ~= 0);
    vecTZeroY    = vecProbeTime_cpu(vecIdxZeroY) ...
                 - vecProbeDispY_cpu(vecIdxZeroY) .* ...
                   (scalDt ./ (vecProbeDispY_cpu(vecIdxZeroY+1) - vecProbeDispY_cpu(vecIdxZeroY)));

%% Plots: displacement magnitude, FFT (omega), attenuation fit
    figure;

    subplot(4,1,1);
    plot(vecProbeTime_cpu, vecProbeDispX_cpu-mean(vecProbeDispX_cpu), 'w-');
    hold on;
    hold off;
    xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('$\Delta x$', 'Interpreter', 'latex', 'FontSize', 20);
    grid on;

    subplot(4,1,2);
    plot(vecProbeTime_cpu, vecProbeDispY_cpu-mean(vecProbeDispY_cpu), 'w-');
    hold on;
    hold off;
    xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('$\Delta y$', 'Interpreter', 'latex', 'FontSize', 20);
    grid on;


    % Attenuation fit on ln|u(t)|
    subplot(4,1,3);
    if abs(scalBetaMode) > 0 && ~isnan(scalBetaFit)
        vecMagPos    = vecProbeMag;
        vecMaskMag   = isfinite(vecMagPos) & (vecMagPos > 0);
        vecTPlot     = vecProbeTime_cpu(vecMaskMag);
        vecLogMag    = log(vecMagPos(vecMaskMag));
        vecPolyCoeff = polyfit(vecTPlot, vecLogMag, 1);
        vecLogFit    = polyval(vecPolyCoeff, vecTPlot);
        vecMagFit    = exp(vecLogFit);

        plot(vecProbeTime_cpu, vecProbeMag, 'Color', [0.2 0.6 0.2]); hold on;
        plot(vecTPlot, vecMagFit, 'r--', 'LineWidth', 1.5);
        hold off;
        xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 20);
        ylabel('$|\mathbf{u}|$', 'Interpreter', 'latex', 'FontSize', 20);
        title(sprintf('$\\beta_{\\mathrm{eig}} = %.4g,\\; \\beta_{\\mathrm{fit}} = %.4g$', ...
            scalBetaMode, scalBetaFit), ...
            'Interpreter', 'latex', 'FontSize', 20);
    else
        plot(vecProbeTime_cpu, vecProbeMag, 'Color', [0.2 0.6 0.2]);
        xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 20);
        ylabel('$|\mathbf{u}_{\mathrm{probe}}(t)|$', 'Interpreter', 'latex', 'FontSize', 20);
        title('Magnitude (attenuation fit skipped)', ...
              'Interpreter', 'latex', 'FontSize', 20);
    end
    grid on;

    subplot(4,1,4);
    plot(vecOmegaFFT, vecAmp, 'w-');
    hold on;
    plot([scalOmegaMode, scalOmegaMode], [0, scalAmpMax], 'r--', 'LineWidth', 1.5);
    hold off;
    xlim([0, 2*scalOmegaMode]);
    xlabel('$\omega$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('$|\tilde{u}(\omega)|$', 'Interpreter', 'latex', 'FontSize', 20);
    title(sprintf('$\\omega{\\mathrm{eig}} = %.4g,\\; \\omega{\\mathrm{max}} = %.4g$', ...
        scalOmegaMode, scalOmegaMax), ...
        'Interpreter', 'latex', 'FontSize', 20);
    grid on;

%% Save / log
    try
        memTest('memlog.txt');
    catch
        % ignore if memTest unavailable
    end

    fprintf('[simEigenmode] DONE (total %.1f s)\n', toc(scalTStart));

%% Error handler
catch sSimEigenmodeME
    fprintf(2, '[simEigenmode] ERROR: %s\n', sSimEigenmodeME.message);
    for scalK = 1:length(sSimEigenmodeME.stack)
        fprintf(2, '[simEigenmode]   in %s (line %d)\n', ...
            sSimEigenmodeME.stack(scalK).name, sSimEigenmodeME.stack(scalK).line);
    end
    rethrow(sSimEigenmodeME);
end

end  % function simEigenmode


function [scalSpringConst, scalMass, scalFreqDriving, scalNumPart, ...
          scalPressure, scalWidth, scalSeed, scalBv, ...
          vecEigenvalues, matEigenvectors, scalGamma, sPacking, ...
          strPackingFile, strResultsFile] = loadSimParams(strInPath)

    sData           = load(strInPath);
    vecEigenvalues  = sData.eigenvalues;
    matEigenvectors = sData.eigenvectors;
    scalGamma       = sData.gamma;
    sPacking        = sData.packing;
    strPackingFile  = sData.packingFile;
    strResultsFile  = sData.resultsFile;

    scalSpringConst = sPacking.K;
    scalNumPart     = sPacking.N;
    scalPressure    = sPacking.P_target;
    scalWidth       = sPacking.Lx;

    cellSeedTok = regexp(strResultsFile, 'Seed(\d+)', 'tokens', 'once');
    scalSeed    = str2double(cellSeedTok{1});

    scalMass        = 1.0;
    scalFreqDriving = 0.0;
    scalBv          = scalGamma;
end


function visSim(x, x0, y, y0, idx, A, dampingConstant, pressureValue, omega, springConstant, sameAxis)

    omegaHat        = omega/sqrt(springConstant);
    dampingConstant = dampingConstant / sqrt(springConstant);

    xColor = [0.1490    0.5490    0.8660];
    yColor = [0.9600    0.4660    0.1600];

    xLimit = max(x0);
    if sameAxis

        subplot(2,1,1)
        plot(x0(idx), x(idx) - x0(idx), 'o', 'Color', xColor)
        ylim(1.2*[-A,A])
        xlim([0,xLimit])
        xlabel('$x_0$', 'Interpreter', 'latex', 'FontSize', 20)
        ylabel('$\Delta x$', 'Interpreter', 'latex', 'FontSize', 20)
        grid on

        subplot(2,1,2)
        plot(x0(idx), y(idx) - y0(idx), 'o', 'Color', yColor)
        ylim(1.2*[-A,A])
        xlim([0,xLimit])
        xlabel('$x_0$', 'Interpreter', 'latex', 'FontSize', 20)
        ylabel('$\Delta y$', 'Interpreter', 'latex', 'FontSize', 20)
        grid on

        sgtitle(['$\hat{\omega} = ', num2str(omegaHat), ...
                 '$, $\hat{\beta} = ', num2str(dampingConstant), ...
                 '$, $\hat{P} = ', num2str(pressureValue), '$'], ...
                'Interpreter', 'latex', 'FontSize', 20);
        drawnow

    else
        subplot(1,2,1)
        plot(x0(idx), x(idx) - x0(idx), 'o')
        xlabel('$x_0$', 'Interpreter', 'latex', 'FontSize', 20)
        ylabel('$\Delta x$', 'Interpreter', 'latex', 'FontSize', 20)
        ylim(1.2*[-A,A])
        xlim([0,xLimit])
        grid on

        subplot(1,2,2)
        plot(x0(idx), y(idx) - y0(idx), 'o')
        xlabel('$x_0$', 'Interpreter', 'latex', 'FontSize', 20)
        ylabel('$\Delta y$', 'Interpreter', 'latex', 'FontSize', 20)
        ylim(1.2*[-A,A])
        xlim([0,xLimit])
        grid on

        sgtitle(['$\hat{\omega} = ', num2str(omegaHat), ...
                 '$, $\hat{\beta} = ', num2str(dampingConstant), ...
                 '$, $\hat{P} = ', num2str(pressureValue), '$'], ...
                'Interpreter', 'latex', 'FontSize', 20);
        drawnow
    end
end

