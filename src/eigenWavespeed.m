function wavespeedVector = plotEigenWavespeed(plotData, pressureList, dampingList)

	gridSpacingX = 0.1;
	gridSpacingY = gridSpacingX;

	figure
	for i = 1:length(pressureList)
		pressureValue = pressureList(i);
		for j = 1:length(dampingList)
			dampingValue = dampingList(j);
			dataPressureDamping = filterData(plotData, 'pressure', pressureValue,  'damping', dampingValue);
			%realEigenValues = real(dataPressureDamping.eigenValues{1});
			%imagEigenValues = imag(dataPressureDamping.eigenValues{1});
			imagEigenValues = imag(dataPressureDamping.eigenValues{1}(imag(dataPressureDamping.eigenValues{1})>0));
			realEigenValues = real(dataPressureDamping.eigenValues{1}(imag(dataPressureDamping.eigenValues{1})>0));
			
			kx = zeros(size(imagEigenValues));
			for k = 1:length(imagEigenValues)
			%for k = 1:200
				freqDamped = imagEigenValues(k);

                % get the eigenvectors and particle positions
				[dx, dy, x0, y0] = plotEigenmode(plotData,...
								pressureValue, ...
								dampingValue, ...
								freqDamped, ...
								"plot", false);

				Lx = max(x0);
				Ly = max(y0);

				xq = 0 : gridSpacingX : Lx;
				yq = 0 : gridSpacingY : Ly;

				[X, Y] = meshgrid(xq, yq);

				%% For x-component of eigen vector, create a surface plot
				vq = griddata(x0,y0,dx, X,Y); % this takes
				vq(isnan(vq))= 0; % clears NaN's (query's outside of actual data)k

				% Just for visual verificaiton of the data
				F = fft2(vq);
				F = fftshift(F); % fft() repeats spectrum after f_NYQ, shift so repeats on negative side
				mag = abs(F) / numel(vq);

				[Ny, Nx] = size(vq);

				% Full frequency vectors (centered around zero)
				kx_xVec = (-floor(Nx/2) : ceil(Nx/2)-1) / (Nx * gridSpacingX);
				kx_yVec = (-floor(Ny/2) : ceil(Ny/2)-1) / (Ny * gridSpacingY);


				% Find peak in full spectrum
				[~, loc] = max(mag(:));
				[iY, iX] = ind2sub(size(mag), loc);
				kx_est = kx_xVec(iX);
				ky_est = kx_yVec(iY);

				%kx(k) = sqrt(kx_est^2 + ky_est^2);
				kx(k) = kx_est;
			end

			if isscalar(pressureList)
				markerColor = [0,0,1];
			else
				[~, markerColor] = normVarColor(pressureList, pressureValue, 1);
			end

			if isscalar(dampingList)
				markerSize = 6;
			else
				markerSize = exp(dataPressureDamping.damping/max(dampingList))*3;
			end

			springConstant = dataPressureDamping.springConstant(1);
			mass = 1;
			dampingDimensionless = dataPressureDamping.damping./sqrt(springConstant*mass);
			pressureLabel = sprintf('$ %.4f, %.4f $', dataPressureDamping.pressure, dampingDimensionless); 
			wavespeed = freqDamped./kx;
			
			
			plot(imagEigenValues,...
				wavespeed,...
				'o', ...
				'MarkerSize', markerSize , ...
				'MarkerEdgeColor', markerColor, ...
				'Color', markerColor, ...
				'DisplayName', pressureLabel);
			xlabel('$\omega_n$', 'Interpreter', 'latex', 'FontSize', 16);
			ylabel('$c_n$', 'Interpreter', 'latex', 'FontSize', 16');
			title(sprintf('$L_x$ by $L_y$: %.2f by %.2f', Lx, Ly), ...
				'Interpreter', 'latex', ...
				'FontSize', 16);

			grid on; 
			hold on;
			set(gca, 'XScale', 'log', 'YScale', 'log');
		end

		legend('show', 'Interpreter', 'latex');
		leg = legend('show', ...
			'Location', 'northeastoutside', ...
			'Interpreter', 'latex', ...
			'FontSize', 15);
		title(leg, "$  \hat{P}, \hat{\Gamma} $")
		ax = gca;
		ax.FontSize = 20;



end
