% tests matSpringDampMass.m, frictionless and frictional (2D)
% matlab should be cd into ~/repos/GranE in order to run with:
% run("tests/testMatSpringDampMass.m")

clearvars

% Resolve the repo root from this file's location so the test works no matter
% what the current folder is (run() executes scripts from their own folder).
strRepoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(strRepoRoot, 'src'));

% ---- Real 2D frictional packing from pack.m ----
% Replace the synthetic perturbed grid with a genuine jammed packing built by
% the DEM compression protocol in pack.m (Cundall-Strack friction, bidisperse
% disks). N is nominally 1000; the saved backbone count (after cleanRats
% rattler removal) can be a little lower.
% pack.m names it 2D_N1000_P0.001_Width32_Seed1_Friction.mat (Width = round(sqrt(1000)) = 32).
strPackDir = fullfile(strRepoRoot, 'data');
strPackingFile = fullfile(strPackDir, '2D_N1000_P0.001_Width32_Seed1_Friction.mat');
packOpts = struct( ...
    'flagFrictionOn',      true, ...  % Cundall-Strack friction (2D only)
    'scalFricCoef',        0.5, ...   % Coulomb coefficient mu
    'scalTangentialK',     1/3, ...   % Kt/K ratio (Cundall-Strack default)
    'scalGammaNormal',     0, ...     % normal dashpot (0 = original scheme)
    'scalGammaTangential', 0, ...     % tangential dashpot (off)
    'saveFrictionalState', false, ...
    'saveFullState',       false);
pack(1000, 100, 1, 1.4, 1, 0.001, 1, false, 1, 1, 0, false, [strPackDir, filesep], packOpts);

% Load the jammed packing and hand it to matSpringDampMass.
sPack = load(strPackingFile);
assert(sPack.N >= 500, ...
    sprintf('frictional packing did not jam: backbone N = %d', sPack.N));

packing = struct( ...
    'vecPosX',        sPack.vecPosX, ...
    'vecPosY',        sPack.vecPosY, ...
    'scalBoxWidthX',  sPack.scalBoxWidthX, ...
    'scalBoxHeightY', sPack.scalBoxHeightY, ...
    'vecDiameter',    sPack.vecDiameter);

% Geometry for the independent energy-form check below.
vecPosX = sPack.vecPosX;
vecPosY = sPack.vecPosY;
vecDiameter = sPack.vecDiameter;

% Clean up the generated file; the test only needs it in memory.
delete(strPackingFile);

N = sPack.N;
kn = 100;
kt = 30;
gn = 0.2;
gt = 0.05;
m  = 1;
tol = 1e-10;

% ----------------- frictionless -----------------
[matK0, matG0, matM0] = matSpringDampMass(packing, gn, kn, m, "periodic", true, "algorithm", "particle");
assert(isequal(size(matK0), [2*N 2*N]), 'frictionless K should be 2N x 2N');
assert(isequal(size(matM0), [2*N 2*N]), 'frictionless M should be 2N x 2N');

% ----------------- frictional -----------------
[matK, matG, matM] = matSpringDampMass(packing, gn, kn, m, "periodic", true, "algorithm", "particle", ...
    "friction", true, "tangentialSpringConstant", kt, "tangentialDampingConstant", gt);

% 1) sizes are 3N x 3N
assert(isequal(size(matK), [3*N 3*N]), 'frictional K should be 3N x 3N');
assert(isequal(size(matG), [3*N 3*N]), 'frictional Gamma should be 3N x 3N');
assert(isequal(size(matM), [3*N 3*N]), 'frictional M should be 3N x 3N');

% 2) symmetric and positive semi-definite
assert(norm(matK - matK.', 'fro') < tol*norm(matK, 'fro'), 'K should be symmetric');
assert(norm(matG - matG.', 'fro') < tol*norm(matG, 'fro'), 'Gamma should be symmetric');
assert(min(eig((matK + matK.')/2)) > -tol*norm(matK, 'fro'), 'K should be positive semi-definite');
assert(min(eig((matG + matG.')/2)) > -tol*norm(matG, 'fro'), 'Gamma should be positive semi-definite');

% 3) rigid translations are zero modes of K and Gamma
vecTransX = repmat([1; 0; 0], N, 1);
vecTransY = repmat([0; 1; 0], N, 1);
assert(norm(matK*vecTransX) < tol*norm(matK, 'fro') && norm(matK*vecTransY) < tol*norm(matK, 'fro'), ...
    'rigid translations should be zero modes of K');
assert(norm(matG*vecTransX) < tol*norm(matG, 'fro') && norm(matG*vecTransY) < tol*norm(matG, 'fro'), ...
    'rigid translations should be zero modes of Gamma');

% 4) quadratic form matches sum over contacts of k_n*delta_n^2 + k_t*zeta_t^2
vecR = randn(3*N, 1);
vecV = randn(3*N, 1);
vecRadius = vecDiameter/2;
vecBox = [packing.scalBoxWidthX ; packing.scalBoxHeightY];
scalEnergyK = 0;
scalEnergyG = 0;
for i = 1:N
    for j = (i + 1):N
        dX = [vecPosX(j) - vecPosX(i) ; vecPosY(j) - vecPosY(i)];
        dX = dX - round(dX./vecBox).*vecBox;
        r  = norm(dX);
        if r < vecRadius(i) + vecRadius(j)
            nhat = [dX/r ; 0];
            ti = [-dX(2)/r ; dX(1)/r ;  vecRadius(i)];
            tj = [-dX(2)/r ; dX(1)/r ; -vecRadius(j)];
            ri = vecR(3*i-2:3*i);  rj = vecR(3*j-2:3*j);
            vi = vecV(3*i-2:3*i);  vj = vecV(3*j-2:3*j);
            scalEnergyK = scalEnergyK + kn*(nhat.'*(rj - ri))^2 + kt*(tj.'*rj - ti.'*ri)^2;
            scalEnergyG = scalEnergyG + gn*(nhat.'*(vj - vi))^2 + gt*(tj.'*vj - ti.'*vi)^2;
        end
    end
end
assert(abs(vecR.'*matK*vecR - scalEnergyK) < 1e-9*abs(scalEnergyK), ...
    'r''*K*r should equal sum of k_n*delta_n^2 + k_t*zeta_t^2 over contacts');
assert(abs(vecV.'*matG*vecV - scalEnergyG) < 1e-9*abs(scalEnergyG), ...
    'v''*Gamma*v should equal sum of g_n*dot(delta_n)^2 + g_t*dot(zeta_t)^2 over contacts');

% 5) mass matrix is diag(m_i, m_i, I_i) with I_i = m_i*R_i^2/2
vecMass = m*pi*vecRadius.^2;
vecMassDiagExpected = reshape([vecMass.' ; vecMass.' ; (0.5*vecMass.*vecRadius.^2).'], 3*N, 1);
assert(norm(diag(matM) - vecMassDiagExpected) < tol*norm(vecMassDiagExpected), ...
    'M should be diag(m_i, m_i, m_i*R_i^2/2)');
assert(nnz(matM - diag(diag(matM))) == 0, 'M should be diagonal');

% 6) frictionless limit: k_t = 0 and g_t = g_n recover the frictionless
%    x,y blocks, and K has no theta entries
[matKLim, matGLim, matMLim] = matSpringDampMass(packing, gn, kn, m, "periodic", true, "algorithm", "particle", ...
    "friction", true, "tangentialSpringConstant", 0, "tangentialDampingConstant", gn);
vecIdxXY = sort([1:3:3*N, 2:3:3*N]);
vecIdxTheta = 3:3:3*N;
assert(norm(matKLim(vecIdxXY, vecIdxXY) - matK0, 'fro') < tol*norm(matK0, 'fro'), ...
    'k_t = 0 should recover the frictionless K');
assert(norm(matKLim(vecIdxTheta, :), 'fro') == 0, 'k_t = 0 should leave no theta entries in K');
assert(norm(matGLim(vecIdxXY, vecIdxXY) - matG0, 'fro') < tol*norm(matG0, 'fro'), ...
    'g_t = g_n should recover the frictionless Gamma on x,y');
assert(norm(matMLim(vecIdxXY, vecIdxXY) - matM0, 'fro') == 0, ...
    'friction should not change the translational masses');

% 7) default tangential constants: k_t = k_n/3 and g_t = g_n
[matKDef, matGDef] = matSpringDampMass(packing, gn, kn, m, "periodic", true, "algorithm", "particle", "friction", true);
[matKExp, matGExp] = matSpringDampMass(packing, gn, kn, m, "periodic", true, "algorithm", "particle", ...
    "friction", true, "tangentialSpringConstant", kn/3, "tangentialDampingConstant", gn);
assert(isequal(matKDef, matKExp) && isequal(matGDef, matGExp), ...
    'default tangential constants should be k_t = k_n/3 and g_t = g_n');

% 8) cell and particle algorithms agree, and sparse matches dense
[matKCell, matGCell, matMCell] = matSpringDampMass(packing, gn, kn, m, "periodic", true, "algorithm", "cell", ...
    "friction", true, "tangentialSpringConstant", kt, "tangentialDampingConstant", gt, "sparseOutput", true);
assert(issparse(matKCell) && issparse(matGCell) && issparse(matMCell), 'sparseOutput should give sparse matrices');
assert(norm(full(matKCell) - matK, 'fro') < tol*norm(matK, 'fro'), 'cell and particle K should agree');
assert(norm(full(matGCell) - matG, 'fro') < tol*norm(matG, 'fro'), 'cell and particle Gamma should agree');
assert(norm(full(matMCell) - matM) == 0, 'sparse and dense M should agree');

% 9) friction is 2D only
packing3D = struct( ...
    'vecPosX', vecPosX, 'vecPosY', vecPosY, 'vecPosZ', 0.5*ones(N, 1), ...
    'scalBoxWidthX', packing.scalBoxWidthX, 'scalBoxHeightY', packing.scalBoxHeightY, 'scalBoxDepthZ', 2, ...
    'vecDiameter', vecDiameter);
boolThrew = false;
try
    matSpringDampMass(packing3D, gn, kn, m, "periodic", true, "friction", true);
catch err
    boolThrew = strcmp(err.identifier, 'matSpringDampMass:Friction3DNotSupported');
end
assert(boolThrew, 'friction = true should error for 3D packings');

% 10) the frictional QEP solves with polyeig
[matPsi, vecLambda] = polyeig(matK, matG, matM);
assert(numel(vecLambda) == 6*N, 'frictional QEP should have 6N eigenvalues');
% rigid-body zero modes are defective (lambda = 0 twice), so polyeig only
% resolves them to ~sqrt(eps); skip them and check the rest decay
vecLambdaNonZero = vecLambda(abs(vecLambda) > 1e-4*max(abs(vecLambda)));
assert(all(real(vecLambdaNonZero) < 0), 'damped non-zero modes should decay');

fprintf('testMatSpringDampMass passed.\n');
