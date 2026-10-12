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
% With friction (2D only) each particle also rotates, so there are f = 3
% degrees of freedom per particle and every matrix is 3N x 3N:
%
%     x(t) = [x_1, y_1, theta_1, ..., x_N, y_N, theta_N]
%
% For a contact (i, j) with unit normal n (pointing from i to j) and unit
% tangent t = z x n = (-n_y, n_x), the projection vectors are
%
%     n_hat = [n_x; n_y; 0],  t_i = [t_x; t_y; +R_i],  t_j = [t_x; t_y; -R_j]
%
% and the contact adds
%
%     K += e_i e_i' ⊗ H_ii + e_j e_j' ⊗ H_jj - e_i e_j' ⊗ H_ij - e_j e_i' ⊗ H_ij'
%     H_ii = k_n n_hat n_hat' + k_t t_i t_i'
%     H_jj = k_n n_hat n_hat' + k_t t_j t_j'
%     H_ij = k_n n_hat n_hat' + k_t t_i t_j'
%
% with Γ built the same way from γ_n and γ_t. The mass matrix is
% diag(m_i, m_i, I_i) per particle with I_i = m_i R_i^2 / 2 (uniform disk).
% Without friction the frictionless blocks k_n n n' (stiffness) and γ_n I
% (damping) are used. Neither K form includes a pre-stress term.
%
% Arguments:
%     packing            Struct with the packing information
%     dampingConstant    Damping parameter (normal damping γ_n with friction)
%     springConstant     Hooke's spring constant (normal stiffness k_n with friction)
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
%     friction           Include tangential (frictional) contacts and particle
%                        rotation, giving 3N x 3N matrices (2D only)
%     tangentialSpringConstant   Tangential stiffness k_t (friction only).
%                                Defaults to springConstant/3, the
%                                Cundall-Strack ratio used by pack.m
%     tangentialDampingConstant  Tangential damping γ_t (friction only).
%                                Defaults to dampingConstant
%
% All logical options default to false.
%
% Returns:
%     matSpring          Stiffness matrix, K (Hessian)
%     matDamp            Damping matrix, Γ
%     matMass            Mass matrix, M
%
% Example (frictional QEP, dense matrices):
%     [matK, matGamma, matM] = matSpringDampMass(packing, gamma, K, 1, ...
%         "periodic", true, "friction", true);
%     [matPsi, vecLambda] = polyeig(matK, matGamma, matM);
%     % solves (matK + lambda*matGamma + lambda^2*matM) * psi = 0

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
        options.friction         (1, 1) logical = false % true ==> tangential contacts + rotation (2D only)
        options.tangentialSpringConstant  (1, 1) double {mustBeReal} = NaN % k_t; NaN ==> springConstant/3
        options.tangentialDampingConstant (1, 1) double {mustBeReal} = NaN % γ_t; NaN ==> dampingConstant
    end

    % Grab the data out of the packing.
    [dim, N, ~, positions, radii, boxDims] = getPackingData(packing, options.oldPackingFormat);
    radii = radii(:);

    % Contact constants and degrees of freedom per particle (f).  With
    % friction each particle also has a rotation theta, so f = dim + 1.
    consts = struct('kn', springConstant, 'kt', 0, 'gn', dampingConstant, 'gt', 0, 'friction', options.friction);
    if (options.friction)
        if (dim ~= 2)
            error("matSpringDampMass:Friction3DNotSupported", "friction = true is 2D only.");
        end
        f = dim + 1;
        consts.kt = options.tangentialSpringConstant;
        if (isnan(consts.kt))
            consts.kt = springConstant/3;
        end
        consts.gt = options.tangentialDampingConstant;
        if (isnan(consts.gt))
            consts.gt = dampingConstant;
        end
    else
        f = dim;
    end
    dof = f*N;

    % Which particles touch the left and right walls (2D)?
    if ((dim == 2) && ~options.periodic)
        maskLeftWall  = positions(:, 1) < radii;
        maskRightWall = positions(:, 1) > boxDims(1) - radii;
    end

    % Build the stiffness and damping matrices.
    switch (options.algorithm)
    case "cell"
        [matSpring, matDamp] = buildMatricesCell(N, dim, f, dof, positions, radii, boxDims, consts, options);
    case "particle"
        [matSpring, matDamp] = buildMatricesParticle(N, dim, f, dof, positions, radii, boxDims, consts, options);
    otherwise
        error("Not implemented.");
    end

    % Add entries for wall contacts (2D, translational DOFs only).
    if ((dim == 2) && ~options.periodic)
        for (i = 1:1:N)
            if (maskLeftWall(i) || maskRightWall(i))
                dofs_i = (f*(i - 1) + 1):(f*(i - 1) + dim);
                matSpring(dofs_i, dofs_i) = matSpring(dofs_i, dofs_i) + springConstant*eye(dim);
            end
        end
    end

    % Build the mass matrix.
    if (options.uniformMass)
        vecMass = mass*ones(N, 1);
    else
        V = pi^(dim/2)/gamma(dim/2 + 1); % Volume of unit ball in R^dim
        vecMass = (mass*V)*radii.^dim;
    end
    if (options.friction)
        vecInertia = 0.5*vecMass.*radii.^2; % Moment of inertia of a uniform disk
        vecMassDiag = reshape([vecMass.' ; vecMass.' ; vecInertia.'], dof, 1);
    else
        vecMassDiag = reshape(repmat(vecMass.', dim, 1), dof, 1);
    end
    if (options.sparseOutput)
        matMass = spdiags(vecMassDiag, 0, dof, dof);
    else
        matMass = diag(vecMassDiag);
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

function [matSpring, matDamp] = buildMatricesCell(N, dim, f, dof, positions, radii, boxDims, consts, options)
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
                    nhat = -dX/r; % Unit normal pointing from i to j
                    [Hii, Hjj, Hij, Gii, Gjj, Gij] = contactBlocks(nhat, radii(i), radii(j), dim, consts);

                    dofs_i = (f*(i - 1) + 1):(f*i);
                    dofs_j = (f*(j - 1) + 1):(f*j);

                    matSpring(dofs_i, dofs_i) = matSpring(dofs_i, dofs_i) + Hii;
                    matSpring(dofs_j, dofs_j) = matSpring(dofs_j, dofs_j) + Hjj;
                    matSpring(dofs_i, dofs_j) = matSpring(dofs_i, dofs_j) - Hij;
                    matSpring(dofs_j, dofs_i) = matSpring(dofs_j, dofs_i) - Hij.';

                    matDamp(dofs_i, dofs_i) = matDamp(dofs_i, dofs_i) + Gii;
                    matDamp(dofs_j, dofs_j) = matDamp(dofs_j, dofs_j) + Gjj;
                    matDamp(dofs_i, dofs_j) = matDamp(dofs_i, dofs_j) - Gij;
                    matDamp(dofs_j, dofs_i) = matDamp(dofs_j, dofs_i) - Gij.';
                end
            end
        end
    end

    % This method double-counts contacts, so divide by 2.
    matSpring = matSpring/2.0;
    matDamp = matDamp/2.0;
end

function [matSpring, matDamp] = buildMatricesParticle(N, dim, f, dof, positions, radii, boxDims, consts, options)
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
                nhat = -dX/r; % Unit normal pointing from i to j
                [Hii, Hjj, Hij, Gii, Gjj, Gij] = contactBlocks(nhat, radii(i), radii(j), dim, consts);

                dofs_i = (f*(i - 1) + 1):(f*i);
                dofs_j = (f*(j - 1) + 1):(f*j);

                matSpring(dofs_i, dofs_i) = matSpring(dofs_i, dofs_i) + Hii;
                matSpring(dofs_j, dofs_j) = matSpring(dofs_j, dofs_j) + Hjj;
                matSpring(dofs_i, dofs_j) = matSpring(dofs_i, dofs_j) - Hij;
                matSpring(dofs_j, dofs_i) = matSpring(dofs_j, dofs_i) - Hij.';

                matDamp(dofs_i, dofs_i) = matDamp(dofs_i, dofs_i) + Gii;
                matDamp(dofs_j, dofs_j) = matDamp(dofs_j, dofs_j) + Gjj;
                matDamp(dofs_i, dofs_j) = matDamp(dofs_i, dofs_j) - Gij;
                matDamp(dofs_j, dofs_i) = matDamp(dofs_j, dofs_i) - Gij.';
            end
        end
    end
end

function [Hii, Hjj, Hij, Gii, Gjj, Gij] = contactBlocks(nhat, Ri, Rj, dim, consts)
% Stiffness (H) and damping (G) blocks of one contact (i, j), with nhat the
% unit normal pointing from i to j.  The contact contributes
%     [ H_ii  -H_ij ; -H_ij'  H_jj ]
% to the (i, j) rows/columns of K, and likewise for Γ.
    if (consts.friction)
        % Projection vectors on (x, y, theta).  t = z x n, and the third
        % entries are the lever arms that turn sliding into spin.
        n  = [nhat ; 0];
        ti = [-nhat(2) ; nhat(1) ;  Ri];
        tj = [-nhat(2) ; nhat(1) ; -Rj];
        nn = n*n.';

        Hii = consts.kn*nn + consts.kt*(ti*ti.');
        Hjj = consts.kn*nn + consts.kt*(tj*tj.');
        Hij = consts.kn*nn + consts.kt*(ti*tj.');

        Gii = consts.gn*nn + consts.gt*(ti*ti.');
        Gjj = consts.gn*nn + consts.gt*(tj*tj.');
        Gij = consts.gn*nn + consts.gt*(ti*tj.');
    else
        Hii = consts.kn*(nhat*nhat.');
        Hjj = Hii;
        Hij = Hii;

        Gii = consts.gn*eye(dim);
        Gjj = Gii;
        Gij = Gii;
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


