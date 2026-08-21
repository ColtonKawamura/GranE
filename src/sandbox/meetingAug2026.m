%% eigenmode sims


simEigenmode("/Users/coltonkawamura/repos/GranE/data/eigenData/results_2D_iso_N100_P0.0001_Seed1_gamma_3.00e-03.mat", scalModeFreq=.1)
simEigenmode("/Users/coltonkawamura/repos/GranE/data/eigenData/results_2D_iso_N100_P0.1_Seed5_gamma_3.00e-03.mat", scalModeFreq=.1)

simEigenmode("/Users/coltonkawamura/repos/GranE/data/eigenData/NEW.mat", scalModeFreq = .1, loadTestData=true)

load("/Users/coltonkawamura/repos/GranE/data/eigenData/results_2D_iso_N100_P0.1_Seed5_gamma_3.00e-03.mat")
%% Toy model wave vector
Lx = 1 ;
Ly = 1;

kx_true = 1;
ky_true = 1;

samplingIntervalX = 0.01;
samplingIntervalY = 0.01;

sampledPointsX = (0 : samplingIntervalX : Lx)';
sampledPointsY = (0 : samplingIntervalY : Ly)';

[matSampleCoordX, matSampleCoordY] = meshgrid(sampledPointsX, sampledPointsY);  % [sampledPointsX, sampledPointsY]

% matrix where each index value is value of the wave
% each index is (indexX, indexY) to get actual 
% value of position for each index, you multiply index by sampling interval
% or look up that index (indexX, indexY) in either sampledPointsX or sampledPointsY
matEigenVectors = cos(2*pi * kx_true * matSampleCoordX + 2*pi * ky_true * matSampleCoordY); % [sampledPointsX, sampledPointsY]

% Just for visual verificaiton of the data
figure
mesh(sampledPointsX, sampledPointsY, matEigenVectors); % this brings them all together
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20);
zlabel('$u(x,y)$', 'Interpreter', 'latex', 'FontSize', 20);
title(sprintf('Wavenumber: $k_x = %d,\\ k_y = %d$', kx_true, ky_true), ...
    'Interpreter', 'latex', 'FontSize', 20);

matFreqAmps = fft2(matEigenVectors); % [sampledPointsX, sampledPointsY] turns amplitude into complex frequency
matFreqAmps = fftshift(matFreqAmps); % shifts to centered around zero
matFreqAmpMag = abs(matFreqAmps) / numel(matEigenVectors); %[sampledPointsX, sampledPointsY] abs

[Ny, Nx] = size(matEigenVectors);

% Full frequency vectors (centered around zero)
% these are just plain bin numbers
% we need to multiply by delta k to get wavenumbers
% the smallest wavenumber is one wave in Lx= 1/Lx
% Lx = Nx * delta x
% so we divide by 1/(Nx*Lx)
kxVec = (-floor(Nx/2) : ceil(Nx/2)-1) / (Nx * samplingIntervalX);
kyVec = (-floor(Ny/2) : ceil(Ny/2)-1) / (Ny * samplingIntervalY);

figure
mesh(kxVec, kyVec, matFreqAmpMag)
xlabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$k_y$', 'Interpreter', 'latex', 'FontSize', 20);
zlabel('$|\hat{U}|$', 'Interpreter', 'latex', 'FontSize', 20);
title(sprintf('Wavenumber: $k_x = %d,\\ k_y = %d$', kx_true, ky_true), ...
    'Interpreter', 'latex', 'FontSize', 20);

% Find peak in full spectrum
[maxMag] = max(matFreqAmpMag(:));
loc_candidates = find(matFreqAmpMag > 0.4 * maxMag);
[iY_cands, iX_cands] = ind2sub(size(matFreqAmpMag), loc_candidates);
kx_cands = kxVec(iX_cands);
ky_cands = kyVec(iY_cands);
dist = (kx_cands - kx_true).^2 + (ky_cands - ky_true).^2;
bestIdx = find(dist == min(dist), 1);
kx_est = kx_cands(bestIdx);
ky_est = ky_cands(bestIdx);

