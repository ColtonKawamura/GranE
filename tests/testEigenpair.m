function testEigenpair(strLoadPath, scalModeFreq, logOldCode)
% testEigenpair(strLoadPath, scalModeFreq, logOldCode)
% Simple check of how well a selected eigenpair (λ, v) satisfies
%   (K + λ C + λ^2 M) v = 0

%% Load file (eigenData-style)
sData           = load(strLoadPath);
vecEigenvalues  = sData.eigenvalues;    % Nmodes×1 or 1×Nmodes
matEigenvectors = sData.eigenvectors;   % (2N)×Nmodes
scalGamma       = sData.gamma;
sPacking        = sData.packing;

scalSpringConst = sPacking.K;
scalNumPart     = sPacking.N;
scalPressure    = sPacking.P_target;
scalLx          = sPacking.Lx;
scalLy          = sPacking.Ly;

fprintf('[testEigenpair] Loaded %s\n', char(strLoadPath));
fprintf('[testEigenpair]   K=%g, gamma=%g, N=%d, P=%g\n', ...
    scalSpringConst, scalGamma, scalNumPart, scalPressure);

%% Optional old-code flip (fix Kxy sign convention)
if logOldCode
    scalNdof = size(matEigenvectors, 1);
    if mod(scalNdof, 2) ~= 0
        error('testEigenpair:OldCodeFlip', ...
              'oldCode=true but DOF count %d is not 2*N.', scalNdof);
    end
    % Assuming ordering [x1; y1; x2; y2; ...]
    matEigenvectors(2:2:end, :) = -matEigenvectors(2:2:end, :);
    fprintf('[testEigenpair] Applied oldCode y-flip to eigenvectors.\n');
end

%% Select mode by target frequency (imaginary part of λ)
vecLambda = vecEigenvalues(:);        % Nmodes×1
vecOmega  = imag(vecLambda);
vecMaskPos = (vecOmega > 0);

if ~any(vecMaskPos)
    error('testEigenpair:NoPositiveFreq', ...
          'No eigenvalues with positive imaginary part found.');
end

vecOmegaPos = vecOmega(vecMaskPos);
vecIdxPos   = find(vecMaskPos);

[~, scalIdxLocal] = min(abs(vecOmegaPos - scalModeFreq));
idxMode           = vecIdxPos(scalIdxLocal);

scalLambdaSel = vecLambda(idxMode);
scalOmegaSel  = imag(scalLambdaSel);
scalBetaSel   = real(scalLambdaSel);

fprintf(['[testEigenpair] Selected mode: target omega=%g, ', ...
         'lambda=%g + %gi (idx=%d)\n'], ...
    scalModeFreq, scalBetaSel, scalOmegaSel, idxMode);

%% Build K, C, M for current packing / damping
matPositions = [sPacking.x(:), sPacking.y(:)];   % N×2
vecRadii     = (sPacking.Dn(:).') / 2;          % 1×N

[matK, matC, matM] = matSpringDampMass( ...
    matPositions, vecRadii, scalLy, scalLx, scalGamma, scalSpringConst, ...
    "periodic", true);

%% Compute residual ||(K + λ C + λ^2 M) v|| / ||v||
vecVMode = matEigenvectors(:, idxMode);

vecResidual = (matK + scalLambdaSel*matC + (scalLambdaSel^2)*matM) * vecVMode;
scalResNorm = norm(vecResidual) / norm(vecVMode);

fprintf('[testEigenpair] Residual ||(K + λ C + λ^2 M)v||/||v|| = %g\n', scalResNorm);
fprintf('[testEigenpair] Done.\n');
end

