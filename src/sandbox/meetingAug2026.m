
%% Toy model wave vector
Lx = 10;
Ly = 10;

kx_true = 3;
ky_true = 0;

for dx = [0.1]%, 0.05, 0.02]

	dy = dx;

	xq = 0 : dx : Lx;
	yq = 0 : dy : Ly;

	[X, Y] = meshgrid(xq, yq);

	vq = sin(2*pi * kx_true * X + 2*pi * ky_true * Y); % Actual mapped data

	% Just for visual verificaiton of the data
	figure
	mesh(xq, yq, vq);
    xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20);
    zlabel('$u(x,y)$', 'Interpreter', 'latex', 'FontSize', 20);

	F = fft2(vq);
	F = fftshift(F); % shifts to centered around zero
	mag = abs(F) / numel(vq);

	[Ny, Nx] = size(vq);

	% Full frequency vectors (centered around zero)
	kxVec = (-floor(Nx/2) : ceil(Nx/2)-1) / (Nx * dx);
	kyVec = (-floor(Ny/2) : ceil(Ny/2)-1) / (Ny * dy);

	figure
	mesh(kxVec, kyVec, mag)
    xlabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$k_y$', 'Interpreter', 'latex', 'FontSize', 14);
    zlabel('$|\hat{U}|$', 'Interpreter', 'latex', 'FontSize', 14);
    title(sprintf('FFT magnitude — peak at $k_x = %.2f,\\; k_y = %.2f$', kx_est, ky_est), 'Interpreter', 'latex', 'FontSize', 14);

	% Find peak in full spectrum
	[~, loc] = max(mag(:));
	[iY, iX] = ind2sub(size(mag), loc);
	kx_est = kxVec(iX);
	ky_est = kyVec(iY);

	fprintf("dx = %.4f → Nx = %d, kx_est = %.4f, ky_est = %.4f\n", dx, Nx, kx_est, ky_est);
end

