function processEigenModesDampedPara(in_path, out_path, dampingConstants, options)
% Processes eigenmodes for damped granular packings with parallel computation
%
% This function loads granular packing data from .mat files, calculates the 
% eigenmodes with various damping constants, and saves the results. Uses 
% parallel processing to handle multiple files simultaneously.
%
% Parameters:
% in_path - Directory containing .mat files with granular packing data
% out_path - Directory to save the processed results
% dampingConstants - Array of damping constants to apply to each packing
% options.periodic - Boolean flag for periodic boundary conditions (default: false)

	arguments
	    in_path (1,1) string
	    out_path (1,1) string
	    dampingConstants (1,:) double
	    options.periodic (1,1) logical = false
	    options.serial (1,1) logical = false
	    options.singleFiles (1,1) logical = false
	end

	filenameList = dir(fullfile(in_path, '*.mat'));
	numFiles = length(filenameList);
	numDamping = length(dampingConstants);
	numTotalCombinations = numFiles * numDamping;

	% --- Pre-allocation for the final structure ---
	if numFiles == 0
	    warning('No .mat files found in %s', in_path);
	    return;
	end


	idx = 0;
	if options.serial
		outData = struct();
		numTotalCombinations = 1; % added this beause was inputting empty data
		outData.pressure       = zeros(numTotalCombinations, 1);
		outData.damping        = zeros(numTotalCombinations, 1);
		outData.eigenVectors  = cell(numTotalCombinations, 1);
		outData.eigenValues   = cell(numTotalCombinations, 1);
		outData.Ly             = zeros(numTotalCombinations, 1);
		outData.Lx             = zeros(numTotalCombinations, 1);
		outData.radii              = cell(numTotalCombinations, 1);
		outData.positions              = cell(numTotalCombinations, 1);
		outData.springConstant = zeros(numTotalCombinations, 1);

	    for i = 1:numFiles

		filename = fullfile(in_path, filenameList(i).name);

		try
		    loadedVars = load(filename, 'x', 'y', 'Dn', 'K', 'Ly', 'Lx', 'P', 'N');
		    mass= 1;
		catch
		    warning("File %s does not contain the expected variables. Skipping...", filename);
		    continue;
		end

		springConstant = loadedVars.K;
		positions = [loadedVars.x', loadedVars.y'];
		radii = loadedVars.Dn' / 2;
		[positions, radii] = cleanRats(positions, radii, loadedVars.Ly, loadedVars.Lx);

		for j = 1:numDamping
			idx = idx +1;	
		    dampingConstant = dampingConstants(j);
		    fprintf('Serial processing pressure %f with damping %d \n', loadedVars.P, dampingConstant);

		    if options.periodic
			[Hessian, matDamp, matMass] = matSpringDampMass(positions, radii, loadedVars.Ly, loadedVars.Lx, dampingConstant, springConstant, "periodic", true);
		    else
			[Hessian, matDamp, matMass] = matSpringDampMass(positions, radii, loadedVars.Ly, loadedVars.Lx, dampingConstant, springConstant);
		    end

		    [eigenVectors_j, eigenValues_j] = polyeig(Hessian, matDamp, matMass);


		    outData.pressure(idx,1) = loadedVars.P;
		    outData.damping(idx,1) = dampingConstant;
		    outData.eigenVectors{idx,1} = eigenVectors_j;;
		    outData.eigenValues{idx,1} = eigenValues_j;
		    outData.radii{idx,1} = radii;
		    outData.Ly(idx,1) = loadedVars.Ly;
		    outData.Lx(idx,1) = loadedVars.Lx;
		    outData.positions{idx,1} = positions;
		    outData.springConstant(idx,1) = loadedVars.K;
		
		    if options.singleFiles
			filename_output = string(sprintf("2D_eigenData_%dby%d_P%.3f_damp%.3f_K%d_M%d.mat", round(loadedVars.Lx), round(loadedVars.Ly), loadedVars.P, dampingConstant, loadedVars.K, mass));
			save_path = fullfile(out_path, filename_output);
			save(save_path, 'outData', '-v7.3'); % need this for file sizes larger than 2G's
			display("Saved to: " + save_path);
		    end
			
		end
	    end

	    if ~options.singleFiles
		filename_output = string(sprintf("2D_eigenData_%dby%d_P%.3f_damp%.3f_K%d_M%d.mat", round(loadedVars.Lx), round(loadedVars.Ly), loadedVars.P, dampingConstant, loadedVars.K, mass));
		save_path = fullfile(out_path, filename_output);
		save(save_path, 'outData', '-v7.3'); % need this for file sizes larger than 2G's
		display("Saved to: " + save_path);
	    end

	else
		fprintf("Performing Parallel Processing\n");

		% --- Pre-allocate the final output structure ---
		outData = struct();
		outData.pressure       = zeros(numTotalCombinations, 1);
		outData.damping        = zeros(numTotalCombinations, 1);
		outData.eigenVectors   = cell(numTotalCombinations, 1);
		outData.eigenValues    = cell(numTotalCombinations, 1);
		outData.Ly             = zeros(numTotalCombinations, 1);
		outData.Lx             = zeros(numTotalCombinations, 1);
		outData.radii          = cell(numTotalCombinations, 1);
		outData.positions      = cell(numTotalCombinations, 1);
		outData.springConstant = zeros(numTotalCombinations, 1);

		% ----------------------------------------------------------------
		% MEMORY OPTIMIZATION: Pre-load all file metadata (positions, radii,
		% pressure, box sizes, spring constant) here in a serial loop.
		%
		% Why? Positions and radii are small (kilobytes), while eigenvectors
		% are huge (hundreds of MB). By separating metadata loading from the
		% eigenvalue computation, we avoid embedding redundant metadata copies
		% inside each worker's result cell. Workers only need to broadcast
		% fileData once and then each stores only the {eigenVectors, eigenValues}
		% pair — not a full struct with all fields repeated numDamping times.
		% ----------------------------------------------------------------
		fileData = cell(numFiles, 1);
		for i = 1:numFiles
			filename = fullfile(in_path, filenameList(i).name);
			try
				loadedVars = load(filename, 'x', 'y', 'Dn', 'K', 'Ly', 'Lx', 'P', 'N');
			catch
				warning("File %s does not contain expected variables. Skipping...", filename);
				fileData{i} = [];
				continue;
			end
			fd = struct();
			fd.positions = [loadedVars.x', loadedVars.y'];
			fd.radii     = loadedVars.Dn' / 2;
			[fd.positions, fd.radii] = cleanRats(fd.positions, fd.radii, loadedVars.Ly, loadedVars.Lx);
			fd.P  = loadedVars.P;
			fd.Ly = loadedVars.Ly;
			fd.Lx = loadedVars.Lx;
			fd.K  = loadedVars.K;
			fileData{i} = fd;
		end

		% ----------------------------------------------------------------
		% PARALLEL EFFICIENCY: Choose the outer (parfor) dimension to be
		% whichever is larger — numFiles or numDamping.  The outer loop
		% determines how many iterations are spread across workers, so
		% assigning the larger dimension to parfor maximises the number of
		% workers that can be busy at once.
		%
		% allEigenData stores only the computed {eigenVectors, eigenValues}
		% pairs.  Metadata (positions, Lx, etc.) is read from fileData
		% during assembly, so it never duplicated inside allEigenData.
		% ----------------------------------------------------------------

		if numFiles >= numDamping
			% ---- CASE A: more files than damping values ----
			% Outer parfor iterates over files; each worker runs all
			% damping values serially for its assigned file.
			fprintf("  Outer loop: files (%d) >= damping (%d) — parfor over files\n", numFiles, numDamping);

			% allEigenData{i}{j} = {eigenVectors, eigenValues} for file i, damping j
			allEigenData = cell(numFiles, 1);

			parfor i = 1:numFiles
				fd = fileData{i};   % sliced input: MATLAB sends only fileData{i} to this worker
				if isempty(fd)
					allEigenData{i} = cell(numDamping, 1); %#ok<PFOUS>
					continue;
				end

				% Allocate a cell to collect results for this file across all dampings
				workerEigen_i = cell(numDamping, 1);

				for j = 1:numDamping
					dampingConstant = dampingConstants(j);
					fprintf('  [parfor file %d] pressure %.4g, damping %.4g\n', i, fd.P, dampingConstant);

					% Build system matrices for this (file, damping) pair
					if options.periodic
						%                   fprintf("***OLD VERSION \n")
						% [Hessian, matDamp, matMass] = OLDmatSpringDampMass(fd.positions, fd.radii, fd.Ly, fd.Lx, dampingConstant, fd.K, "periodic", true);
						[Hessian, matDamp, matMass] = matSpringDampMass(fd.positions, fd.radii, fd.Ly, fd.Lx, 0, dampingConstant, fd.K, "periodic", true);
                        Hessian(1:5, 1:5)
					else
						[Hessian, matDamp, matMass] = matSpringDampMass(fd.positions, fd.radii, fd.Ly, fd.Lx, dampingConstant, fd.K);
					end

					% Solve the quadratic eigenvalue problem: (Hessian + lambda*matDamp + lambda^2*matMass)*v = 0
					[eigenVectors_j, eigenValues_j] = polyeig(Hessian, matDamp, matMass);

					% Store only the eigenpair — NOT the full metadata struct.
					% This is the key memory saving: positions/radii (~KB) are NOT
					% duplicated numDamping times inside the worker results.
					workerEigen_i{j} = {eigenVectors_j, eigenValues_j};
				end

				allEigenData{i} = workerEigen_i;
			end

			% --- Assemble results from Case A ---
			for i = 1:numFiles
				fd = fileData{i};
				if isempty(fd), continue; end
				for j = 1:numDamping
					% Linear index: rows are ordered by file, then by damping
					idx = (i - 1) * numDamping + j;
					eigenPair = allEigenData{i}{j};
					outData.pressure(idx)      = fd.P;
					outData.damping(idx)        = dampingConstants(j);
					outData.eigenVectors{idx}   = eigenPair{1};
					outData.eigenValues{idx}    = eigenPair{2};
					outData.radii{idx}          = fd.radii;
					outData.Ly(idx)             = fd.Ly;
					outData.Lx(idx)             = fd.Lx;
					outData.positions{idx}      = fd.positions;
					outData.springConstant(idx) = fd.K;
				end
			end

		else
			% ---- CASE B: more damping values than files ----
			% Outer parfor iterates over damping; each worker loads all
			% files serially for its assigned damping value.
			% This keeps more workers busy when numDamping > numFiles.
			fprintf("  Outer loop: damping (%d) > files (%d) — parfor over damping\n", numDamping, numFiles);

			% allEigenData{j}{i} = {eigenVectors, eigenValues} for file i, damping j
			allEigenData = cell(numDamping, 1);

			% fileData is broadcast to all workers (it only contains small
			% metadata — positions are kilobytes, not gigabytes).
			parfor j = 1:numDamping
				dampingConstant = dampingConstants(j);

				% Allocate a cell to collect results across all files for this damping
				workerEigen_j = cell(numFiles, 1);

				for i = 1:numFiles
					fd = fileData{i};   % broadcast variable — full fileData sent to each worker
					if isempty(fd), continue; end
					fprintf('  [parfor damp %d] pressure %.4g, damping %.4g\n', j, fd.P, dampingConstant);

					% Build system matrices
					if options.periodic
						[Hessian, matDamp, matMass] = matSpringDampMass(fd.positions, fd.radii, fd.Ly, fd.Lx, dampingConstant, fd.K, "periodic", true);
					else
						[Hessian, matDamp, matMass] = matSpringDampMass(fd.positions, fd.radii, fd.Ly, fd.Lx, dampingConstant, fd.K);
					end

					[eigenVectors_j, eigenValues_j] = polyeig(Hessian, matDamp, matMass);

					% Store only the eigenpair (no metadata duplication)
					workerEigen_j{i} = {eigenVectors_j, eigenValues_j};
				end

				allEigenData{j} = workerEigen_j;
			end

			% --- Assemble results from Case B ---
			% Index mapping: same as Case A — idx = (i-1)*numDamping + j
			for j = 1:numDamping
				for i = 1:numFiles
					fd = fileData{i};
					if isempty(fd), continue; end
					idx = (i - 1) * numDamping + j;
					eigenPair = allEigenData{j}{i};
					outData.pressure(idx)      = fd.P;
					outData.damping(idx)        = dampingConstants(j);
					outData.eigenVectors{idx}   = eigenPair{1};
					outData.eigenValues{idx}    = eigenPair{2};
					outData.radii{idx}          = fd.radii;
					outData.Ly(idx)             = fd.Ly;
					outData.Lx(idx)             = fd.Lx;
					outData.positions{idx}      = fd.positions;
					outData.springConstant(idx) = fd.K;
				end
			end
		end

		% --- Saving ---
		% Use the pre-loaded fileData instead of re-reading a file just to
		% get the spring constant (avoids the extra disk I/O that was
		% needed before because parfor couldn't see the last K value).
		first_valid_idx = find(outData.pressure ~= 0 | ~cellfun('isempty', outData.eigenValues), 1, 'first');
		if ~isempty(first_valid_idx)
			N_save  = length(outData.positions{first_valid_idx});
			Lx_save = outData.Lx(first_valid_idx);
			Ly_save = outData.Ly(first_valid_idx);
			K_save  = outData.springConstant(first_valid_idx);
			mass_save = 1;
			filename_output = sprintf("2D_damped_eigenstuff_N%d_%dby%d_K%d_M%d.mat", N_save, round(Lx_save), round(Ly_save), K_save, mass_save);
		else
			error('No valid data found to save.');
		end

		save_path = fullfile(out_path, filename_output);
		save(save_path, 'outData', '-v7.3');
		fprintf("Saved %d entries to: %s\n", numTotalCombinations, save_path);
	end
end
