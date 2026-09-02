function [matSpring, matDamp, matMass] = matSpringDampMass(packing, dampingConstant, springConstant, mass, options)
% Computes the system matrices K, Γ, and M for a granular packing (2D or 3D)
% with Hooke's force law described by the system of ODEs
%
%     M*x''(t) + Γ*x'(t) + K*x(t) = 0
%
% The matrices correspond to a listing of the degrees of freedom in x(t) are
% in the following order:
%
%     x(t) = [x-pos. of particle 1
%             y-pos. of particle 1
%             z-pos. of particle 1
%             ...
%             x-pos. of particle N
%             y-pos. of particle N
%             z-pos. of particle N]
%
% Arguments:
%     packing            Struct with the packing information
%     dampingConstant    Damping parameter
%     springConstant     Hooke's spring constant
%     mass               Mass of a particle of unit radius
%
% Options/Keyword Arguments:
%     algorithm          Select algorithm for checking particle contacts:
%                          "cell"     - Cell-based method (linear scaling, default, faster for large packings)
%                          "particle" - Check all pairs of particles (quadratic scaling, faster for small packings)
%     periodic           Use periodic boundary conditions
%     sparseOutput       Generate sparse matrices
%     uniformMass        Force all particles (even of different sizes) to have the same mass
%     oldPackingFormat   Set to true if packing is stored in the old format
%
% All options default to false.
%
% Returns:
%     matSpring          Stiffness matrix, K (Hessian)
%     matDamp            Damping matrix, Γ
%     matMass            Mass matrix, M

    arguments
        packing                  (1, 1) struct
        dampingConstant          (1, 1) double {mustBeReal} = 1.0
        springConstant           (1, 1) double {mustBeReal} = 1.0
        mass                     (1, 1) double {mustBeReal} = 1.0
        options.algorithm        (1, 1) string = "cell"
        options.periodic         (1, 1) logical = false % true ==> periodic boundary condition
        options.sparseOutput     (1, 1) logical = false % true ==> output sparse matrices
        options.uniformMass      (1, 1) logical = false % true ==> assume all particles are of the same mass
        options.oldPackingFormat (1, 1) logical = false % true ==> assume packing is stored in the old format
    end

    % Grab the data out of the packing.
    [dim, N, dof, positions, radii, boxDims] = getPackingData(packing, options.oldPackingFormat);

    % Which particles touch the left and right walls (2D)?
    if ((dim == 2) && ~options.periodic)
        maskLeftWall  = positions(:, 1) < radii;
        maskRightWall = positions(:, 1) > boxDims(1) - radii;
    end

    % Build the stiffness and damping matrices.
    switch (options.algorithm)
    case "cell"
        [matSpring, matDamp] = buildMatricesCell(N, dim, dof, positions, radii, boxDims, springConstant, options);
    case "particle"
        [matSpring, matDamp] = buildMatricesParticle(N, dim, dof, positions, radii, boxDims, springConstant, options);
    otherwise
        error("Not implemented.");
    end

    % Add entries for wall contacts (2D).
    if ((dim == 2) && ~options.periodic)
        for (i = 1:1:N)
            if (maskLeftWall(i) || maskRightWall(i))
                dofs_i = (dim*(i - 1) + 1):(dim*i);
                matSpring(dofs_i, dofs_i) = matSpring(dofs_i, dofs_i) + springConstant*eye(dim);
            end
        end
    end

    % Apply the damping coefficient.
    matDamp = dampingConstant*matDamp;

    % Build the mass matrix.
    if (options.uniformMass)
        if (options.sparseOutput)
            matMass = mass*speye(dof);
        else
            matMass = mass*eye(dof);
        end
    else
        V = pi^(dim/2)/gamma(dim/2 + 1); % Volume of unit ball in R^dim
        r = reshape(repmat(radii.', dim, 1), dof, 1);
        if (options.sparseOutput)
            matMass = (mass*V)*spdiags(r.^dim, 0, dof, dof);
        else
            matMass = (mass*V)*diag(r.^dim);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [matSpring, matDamp] = allocateMatrices(dof, dim, sparse)
    if (sparse)
        % TODO:  Do these over-allocate?
        matSpring = spalloc(dof, dof, ceil(0.005*(dof^2)));
        matDamp   = spalloc(dof, dof, ceil(0.005*(dof^2)));
    else
        matSpring = zeros(dof, dof);
        matDamp   = zeros(dof, dof);
    end
end

function [matSpring, matDamp] = buildMatricesCell(N, dim, dof, positions, radii, boxDims, springConstant, options)
    % Pre-allocate the outputs.
    [matSpring, matDamp] = allocateMatrices(dof, dim, options.sparseOutput);

    % Data for the cell-based method for checking particle contacts.
    baseCellWidth    = 4.0*max(radii);                 % Nominal width of each cell
    cellCounts       = round(boxDims/baseCellWidth);   % Number of cells in each dimension
    cellWidths       = boxDims./cellCounts;            % Actual cell widths
    cellCounts       = ceil(boxDims./cellWidths);      % Update cell counts to use actual widths
    particlesToCells = ceil(positions./cellWidths.');  % ith row contains cell numbers for ith particle
    cellsToParticles = cell(cellCounts.');             % (i, j, k) entry lists particles in cell (i, j, k)

    % Build the map from cell numbers the list of particles each cell contains.
    ii = cell(1, dim);
    for (I = 1:prod(cellCounts))
        [ii{:}] = ind2sub(cellCounts, I);
        cellsToParticles{ii{:}} = find(all(particlesToCells == [ii{:}], 2)).';
    end

    % Loop over all cells looking for contacts between particles in this cell
    % and in adjacent ones.
    for (I = 1:prod(cellCounts))
        [ii{:}] = ind2sub(cellCounts, I);

        % Get a list of all the particles that are in either this cell or a
        % neighboring one.  The call to unique() eliminates duplicates in the
        % event that the cell count is very low.
        %
        % TODO:  This is not very efficient.
        neighboringParticles = cell(1, 3^dim);
        S = cell(1, dim);
        for (n = 1:1:(3^dim))
            [S{:}] = ind2sub(3*ones(1, dim), n);
            s = cell2mat(S) - 2;

            if (options.periodic)
                jj = mat2cell(mod(cell2mat(ii) + s - 1, cellCounts.') + 1, 1, ones(1, dim));
            else
                error('Not implemented.');
            end

            neighboringParticles(n) = cellsToParticles(jj{:});
        end
        neighboringParticles = unique(horzcat(neighboringParticles{:}));

        % Check contacts for all particles in this cell.
        for (i = cellsToParticles{I})
            for (j = neighboringParticles)
                if (j == i)
                    continue;
                end

                % Compute distance between particles i and j.
                dX = (positions(i, :) - positions(j, :)).';
                if (options.periodic)
                    dX = dX - round(dX./boxDims).*boxDims;
                end
                r = norm(dX);

                % If the particles are in contact, add terms for the contact forces.
                if (r < radii(i) + radii(j))
                    nhat = dX/r;
                    Kblock = springConstant*(nhat*nhat');
                    Dblock = eye(dim);

                    dofs_i = (dim*(i - 1) + 1):(dim*i);
                    dofs_j = (dim*(j - 1) + 1):(dim*j);

                    matSpring(dofs_i, dofs_i) = matSpring(dofs_i, dofs_i) + Kblock;
                    matSpring(dofs_j, dofs_j) = matSpring(dofs_j, dofs_j) + Kblock;
                    matSpring(dofs_i, dofs_j) = matSpring(dofs_i, dofs_j) - Kblock;
                    matSpring(dofs_j, dofs_i) = matSpring(dofs_j, dofs_i) - Kblock;

                    matDamp(dofs_i, dofs_i) = matDamp(dofs_i, dofs_i) + Dblock;
                    matDamp(dofs_j, dofs_j) = matDamp(dofs_j, dofs_j) + Dblock;
                    matDamp(dofs_i, dofs_j) = matDamp(dofs_i, dofs_j) - Dblock;
                    matDamp(dofs_j, dofs_i) = matDamp(dofs_j, dofs_i) - Dblock;
                end
            end
        end
    end

    % This method double-counts contacts, so divide by 2.
    matSpring = matSpring/2.0;
    matDamp = matDamp/2.0;
end

function [matSpring, matDamp] = buildMatricesParticle(N, dim, dof, positions, radii, boxDims, springConstant, options)
    % Pre-allocate the outputs.
    [matSpring, matDamp] = allocateMatrices(dof, dim, options.sparseOutput);

    % Check for contact between each pair of particles to generate the
    % stiffness and damping matrices.
    for (i = 1:1:N)
        for (j = (i + 1):N)
            % Compute distance between particles i and j.
            dX = (positions(i, :) - positions(j, :)).';
            if (options.periodic)
                dX = dX - round(dX./boxDims).*boxDims;
            end
            r = norm(dX);

            % If the particles are in contact, add terms for the contact forces.
            if (r < radii(i) + radii(j))
                nhat = dX/r;
                Kblock = springConstant*(nhat*nhat');
                Dblock = eye(dim);

                dofs_i = (dim*(i - 1) + 1):(dim*i);
                dofs_j = (dim*(j - 1) + 1):(dim*j);

                matSpring(dofs_i, dofs_i) = matSpring(dofs_i, dofs_i) + Kblock;
                matSpring(dofs_j, dofs_j) = matSpring(dofs_j, dofs_j) + Kblock;
                matSpring(dofs_i, dofs_j) = matSpring(dofs_i, dofs_j) - Kblock;
                matSpring(dofs_j, dofs_i) = matSpring(dofs_j, dofs_i) - Kblock;

                matDamp(dofs_i, dofs_i) = matDamp(dofs_i, dofs_i) + Dblock;
                matDamp(dofs_j, dofs_j) = matDamp(dofs_j, dofs_j) + Dblock;
                matDamp(dofs_i, dofs_j) = matDamp(dofs_i, dofs_j) - Dblock;
                matDamp(dofs_j, dofs_i) = matDamp(dofs_j, dofs_i) - Dblock;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dim, N, dof, positions, radii, boxDims] = getPackingData(packing, oldPackingFormat)
    if (oldPackingFormat)              % Old packing format
        if (isfield(packing, 'z'))
            dim = 3;
            positions = [packing.x.' packing.y.' packing.z.'];
            boxDims   = [packing.Lx ; packing.Ly ; packing.Lz];
        else
            dim = 2;
            positions = [packing.x.' packing.y.'];
            boxDims   = [packing.Lx ; packing.Ly];
        end
        radii = packing.Dn.'/2.0;
    else                               % New packing format
        if (isfield(packing, 'vecPosZ'))
            dim = 3;
            positions = [packing.vecPosX packing.vecPosY packing.vecPosZ];
            boxDims   = [packing.scalBoxWidthX ; packing.scalBoxHeightY ; packing.scalBoxDepthZ];
        else
            dim = 2;
            positions = [packing.vecPosX packing.vecPosY];
            boxDims   = [packing.scalBoxWidthX ; packing.scalBoxHeightY];
        end
        radii = packing.vecDiameter/2.0;
    end

    N = size(positions, 1);
    dof = N*dim;
end


