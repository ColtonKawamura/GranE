function [matSpring, matDamp, matMass, sContactInfo] = matSpringDampMass(packing, dampingConstant, springConstant, mass, options)
% Computes the system matrices K, Γ, and M for a granular packing (2D or 3D)
% with Hooke's force law (default) or the Hertz(-Mindlin) force law (see
% "contactLaw" below) described by the system of ODEs
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
% (damping) are used. By default neither K form includes a pre-stress term
% (see "stressed" below).
%
% ---------------------------------------------------------------------------
% Optional physics for frictional eigenmodes near jamming
% ---------------------------------------------------------------------------
% These options follow the dynamical matrix of Henkes, van Hecke & van
% Saarloos, EPL 90, 14003 (2010), eq. (2), whose packings are prepared as in
% Shundyak, van Hecke & van Saarloos, PRE 75, 010301(R) (2007). Henkes et al.
% use Hertz-Mindlin contacts; by default k_n and k_t here are Hookean
% constants, and contactLaw = "hertz" switches to Hertz-Mindlin (option 3).
% All options are off by default, so K, Γ and M are unchanged unless they are
% requested.
%
% 1) stressed = true: the pre-stress (geometric) term.
%    The exact Hessian of a pair potential V(r) for one contact is
%        V''(r) n n' + (V'(r)/r) (I - n n').
%    For Hooke's law V'' = k_n is a constant, but V'(r) = -f_n with
%    f_n = k_n (R_i + R_j - r) > 0 is not zero in a compressed packing.
%    Sliding the pair sideways by u lengthens the center distance to
%    r + |u|^2/(2r), so the compressed spring releases energy f_n |u|^2/(2r).
%    Each contact therefore also adds
%        S = -(f_n/r) (I - n n')     (translational DOFs only)
%    to H_ii, H_jj and H_ij. This is the -(f_n/r)(dr.t)^2 term of Henkes et
%    al. eq. (2). Its size relative to k_n is the overlap/r ~ pressure, so it
%    barely moves phonons or elastic moduli, but it sets the modes below the
%    plateau frequency omega* and, with friction, it is the only restoring
%    term for rolling motions (zeta_t = 0) on which k_t does no work.
%    Leaving it out (the default) is the "unstressed" approximation: each
%    contact is replaced by a relaxed spring of the same stiffness (Wyart,
%    Silbert, Nagel & Witten, PRE 72, 051306 (2005); Mizuno, Saitoh & Silbert,
%    arXiv:1606.06599). The matching geometric term of the tangential force,
%    of order f_t/r <= mu f_n/r, is neglected, as in Henkes et al. With the
%    term on, K can have small negative eigenvalues if the packing is not in
%    force balance, and Γ is no longer proportional to K for frictionless
%    packings.
%
% 2) slipMobilized = true (friction only, needs fricState): slipping contacts.
%    The mobilization of a contact is m = |f_t|/(mu f_n), between 0 and 1.
%    A contact with m = 1 sits on the Coulomb cone: pushed further in the
%    direction it is already loaded, it slides at constant f_t and gives no
%    tangential stiffness. Following Henkes et al., such contacts slip in both
%    directions under small vibrations and get k_t = 0; contacts with m < 1
%    keep k_t. The tangential dashpot γ_t stays on slipping contacts, because
%    pack.m applies its tangential dashpot outside the Coulomb cap. With this
%    choice the floppy (zero-frequency) modes appear at the generalized
%    isostatic contact number z = 3 + n_m in 2D, with n_m the number of fully
%    mobilized contacts per particle.
%    The contact forces come from the fricState struct in the one file pack.m
%    writes for a frictional packing (*_Friction.mat, or *_Hertz_Friction.mat).
%    Its top-level fricState is already re-indexed onto the backbone, so pass
%    that file (the backbone positions) as the packing together with its
%    fricState. (sFullState holds the pre-cleanRats positions and an
%    un-reindexed copy, sFullState.fricState; those two also match each other.)
%    A contact counts as fully mobilized when m >= 1 - mobilizationTol.
%    fricState stores the normal force of the contact law the packing was
%    made with (pack.m's Hookean or Hertzian force).
%
% 3) contactLaw = "hertz": Hertz normal force, and Hertz-Mindlin tangential
%    stiffness with friction, matching pack.m with hertzian = true (and
%    hertzMindlin = true). For a contact with overlap delta = R_i + R_j - r and
%    R_eff = R_i R_j / (R_i + R_j), pack.m uses
%        f_n = (4/3) K sqrt(R_eff) delta^(3/2),
%    so the linearized normal stiffness is
%        k_n = df_n/d(delta) = 2 K sqrt(R_eff delta),
%    and Mindlin's tangential stiffness is proportional to the contact radius
%    sqrt(R_eff delta) as well:
%        k_t = (k_t0/K) k_n,   k_t0 = tangentialSpringConstant (default K/3).
%    Both stiffnesses therefore vary from contact to contact (k ~ delta^(1/2),
%    Henkes et al.), and the pre-stress term uses the Hertz f_n. Because pack.m
%    advances the Mindlin tangential force incrementally (dF_t = -k_t dU), a
%    change of the overlap does not change F_t to first order, so eq. (2) of
%    Henkes et al. is the complete linearization. Damping stays constant.
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
%                        (pack.m uses one mass for every particle)
%     oldPackingFormat   Set to true if packing is stored in the old format
%     friction           Include tangential (frictional) contacts and particle
%                        rotation, giving 3N x 3N matrices (2D only)
%     tangentialSpringConstant   Tangential stiffness k_t (friction only).
%                                Defaults to springConstant/3, the
%                                Cundall-Strack ratio used by pack.m
%     tangentialDampingConstant  Tangential damping γ_t (friction only).
%                                Defaults to dampingConstant
%     stressed           Add the pre-stress term -(f_n/r)(I - n n') to K
%     fricState          fricState struct saved by pack.m (per-contact forces)
%     slipMobilized      Set k_t = 0 on fully mobilized contacts (needs
%                        friction = true and fricState)
%     mobilizationTol    m >= 1 - mobilizationTol counts as fully mobilized
%                        (default 1e-6)
%     contactLaw         "hooke" (default) or "hertz" (Hertz-Mindlin with
%                        friction); see option 3 above
%
% All logical options default to false.
%
% Returns:
%     matSpring          Stiffness matrix, K (Hessian)
%     matDamp            Damping matrix, Γ
%     matMass            Mass matrix, M
%     sContactInfo       Struct with the contact network that K was built from:
%                          numContacts       particle-particle contacts
%                          vecNumContacts    [N x 1] contacts per particle
%                          numMobilized      contacts with m >= 1 - mobilizationTol
%                                            (NaN without fricState)
%                          numNotInFricState contacts found here but absent from
%                                            fricState (treated as m < 1)
%                          vecMobilization   m for every contact listed in
%                                            fricState (empty without it)
%
% Example (frictional QEP, dense matrices):
%     [matK, matGamma, matM] = matSpringDampMass(packing, gamma, K, 1, ...
%         "periodic", true, "friction", true);
%     [matPsi, vecLambda] = polyeig(matK, matGamma, matM);
%     % solves (matK + lambda*matGamma + lambda^2*matM) * psi = 0
%
% Example (undamped benchmark against Henkes et al., pack.m frictional run;
% the one output file holds the backbone and its re-indexed fricState):
%     sPack = load('2D_N1000_P0.001_Width32_Seed1_Friction.mat');
%     [matK, ~, matM, info] = matSpringDampMass(sPack, 0, sPack.K, 1, ...
%         "periodic", true, "uniformMass", true, "friction", true, ...
%         "stressed", true, "fricState", sPack.fricState, "slipMobilized", true);
%     vecOmega2 = eig(full(matK), full(matM));   % Γ = 0, so lambda = ±i*omega

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
        options.stressed         (1, 1) logical = false % true ==> add pre-stress term -(f_n/r)(I - n n')
        options.fricState        (1, 1) struct = struct() % pack.m fricState (per-contact forces)
        options.slipMobilized    (1, 1) logical = false % true ==> k_t = 0 on fully mobilized contacts
        options.mobilizationTol  (1, 1) double {mustBeNonnegative} = 1e-6 % m >= 1 - tol ==> fully mobilized
        options.contactLaw       (1, 1) string {mustBeMember(options.contactLaw, ["hooke", "hertz"])} = "hooke"
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

    % Pre-stress term, slipping (fully mobilized) contacts and contact law.
    consts.boolStressed = options.stressed;
    consts.boolSlipMobilized = options.slipMobilized;
    consts.boolUseFricState = ~isempty(fieldnames(options.fricState));
    consts.boolHertz = (options.contactLaw == "hertz");
    consts.scalKtOverKn = consts.kt/consts.kn;   % Mindlin: k_t = (k_t0/K) k_n(delta)
    if (options.slipMobilized && ~(options.friction && consts.boolUseFricState))
        error("matSpringDampMass:SlipNeedsFricState", ...
            "slipMobilized = true needs friction = true and the fricState saved by pack.m.");
    end
    vecMobilization = [];
    if (consts.boolUseFricState)
        [consts.matMobilized, consts.matInFricState, vecMobilization] = ...
            mobilizedContacts(options.fricState, N, radii, options.mobilizationTol);
    end

    % Which particles touch the left and right walls (2D)?
    if ((dim == 2) && ~options.periodic)
        maskLeftWall  = positions(:, 1) < radii;
        maskRightWall = positions(:, 1) > boxDims(1) - radii;
    end

    % Build the stiffness and damping matrices.
    switch (options.algorithm)
    case "cell"
        [matSpring, matDamp, sCounts] = buildMatricesCell(N, dim, f, dof, positions, radii, boxDims, consts, options);
    case "particle"
        [matSpring, matDamp, sCounts] = buildMatricesParticle(N, dim, f, dof, positions, radii, boxDims, consts, options);
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

    % Summary of the contact network K was built from.
    sContactInfo = struct();
    sContactInfo.numContacts = sum(sCounts.vecNumContacts)/2;
    sContactInfo.vecNumContacts = sCounts.vecNumContacts;
    if (consts.boolUseFricState)
        sContactInfo.numMobilized = sCounts.numMobilized;
        sContactInfo.numNotInFricState = sCounts.numNotInFricState;
        if (sCounts.numNotInFricState > 0)
            warning("matSpringDampMass:ContactsNotInFricState", ...
                "%d contacts are not in fricState and are treated as not mobilized.", sCounts.numNotInFricState);
        end
    else
        sContactInfo.numMobilized = NaN;
        sContactInfo.numNotInFricState = NaN;
    end
    sContactInfo.vecMobilization = vecMobilization;
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

function [matSpring, matDamp, sCounts] = buildMatricesCell(N, dim, f, dof, positions, radii, boxDims, consts, options)
    % Pre-allocate the outputs.
    [matSpring, matDamp] = allocateMatrices(dof, dim, options.sparseOutput);
    sCounts = struct('vecNumContacts', zeros(N, 1), 'numMobilized', 0, 'numNotInFricState', 0);

    % Data for the cell-based method for checking particle contacts.
    % The cells must tile the box exactly (cellCounts.*cellWidths == boxDims) so
    % that the periodic cell wrap matches the minimum-image wrap. Recomputing the
    % count as ceil(boxDims./cellWidths) is a floating-point trap: boxDims./
    % (boxDims./k) can come out as k + 2e-15, which ceil turns into k + 1 cells
    % that no longer tile the box, and contacts across the boundary are dropped.
    baseCellWidth    = 4.0*max(radii);                         % Nominal width of each cell
    cellCounts       = max(round(boxDims/baseCellWidth), 1);   % Number of cells in each dimension
    cellWidths       = boxDims./cellCounts;                    % Actual cell widths (tile the box)
    particlesToCells = min(max(ceil(positions./cellWidths.'), 1), cellCounts.'); % ith row: cell numbers of particle i
    cellsToParticles = cell(cellCounts.');                     % (i, j, k) entry lists particles in cell (i, j, k)

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
                    [boolSlip, boolKnown] = lookupContact(i, j, consts);
                    [Hii, Hjj, Hij, Gii, Gjj, Gij] = contactBlocks(nhat, radii(i), radii(j), r, dim, consts, boolSlip);

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

                    % Each contact is visited once from each side, so only
                    % particle i is counted here and the totals are halved.
                    sCounts.vecNumContacts(i) = sCounts.vecNumContacts(i) + 1;
                    sCounts.numMobilized      = sCounts.numMobilized + boolSlip/2;
                    sCounts.numNotInFricState = sCounts.numNotInFricState + (consts.boolUseFricState && ~boolKnown)/2;
                end
            end
        end
    end

    % This method double-counts contacts, so divide by 2.
    matSpring = matSpring/2.0;
    matDamp = matDamp/2.0;
end

function [matSpring, matDamp, sCounts] = buildMatricesParticle(N, dim, f, dof, positions, radii, boxDims, consts, options)
    % Pre-allocate the outputs.
    [matSpring, matDamp] = allocateMatrices(dof, dim, options.sparseOutput);
    sCounts = struct('vecNumContacts', zeros(N, 1), 'numMobilized', 0, 'numNotInFricState', 0);

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
                [boolSlip, boolKnown] = lookupContact(i, j, consts);
                [Hii, Hjj, Hij, Gii, Gjj, Gij] = contactBlocks(nhat, radii(i), radii(j), r, dim, consts, boolSlip);

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

                sCounts.vecNumContacts(i) = sCounts.vecNumContacts(i) + 1;
                sCounts.vecNumContacts(j) = sCounts.vecNumContacts(j) + 1;
                sCounts.numMobilized      = sCounts.numMobilized + boolSlip;
                sCounts.numNotInFricState = sCounts.numNotInFricState + (consts.boolUseFricState && ~boolKnown);
            end
        end
    end
end

function [Hii, Hjj, Hij, Gii, Gjj, Gij] = contactBlocks(nhat, scalRadiusI, scalRadiusJ, scalDist, dim, consts, boolSlip)
% Stiffness (H) and damping (G) blocks of one contact (i, j), with nhat the
% unit normal pointing from i to j and scalDist the center distance.  The
% contact contributes
%     [ H_ii  -H_ij ; -H_ij'  H_jj ]
% to the (i, j) rows/columns of K, and likewise for Γ.  boolSlip marks a fully
% mobilized contact, which loses its tangential spring when slipMobilized is on.

    % Contact stiffnesses and normal force of the contact law.
    if (consts.boolHertz)
        % Hertz(-Mindlin), as in pack.m: f_n = (4/3) K sqrt(R_eff) delta^(3/2),
        % k_n = 2 K sqrt(R_eff delta), k_t = (k_t0/K) k_n.
        scalOverlap   = scalRadiusI + scalRadiusJ - scalDist;
        scalRadiusEff = scalRadiusI*scalRadiusJ/(scalRadiusI + scalRadiusJ);
        scalKn = 2*consts.kn*sqrt(scalRadiusEff*scalOverlap);
        scalFn = (4/3)*consts.kn*sqrt(scalRadiusEff)*scalOverlap^(3/2);
        scalKt = consts.scalKtOverKn*scalKn;
    else
        scalKn = consts.kn;
        scalFn = consts.kn*(scalRadiusI + scalRadiusJ - scalDist);
        scalKt = consts.kt;
    end

    if (consts.friction)
        % Projection vectors on (x, y, theta).  t = z x n, and the third
        % entries are the lever arms that turn sliding into spin.
        n  = [nhat ; 0];
        ti = [-nhat(2) ; nhat(1) ;  scalRadiusI];
        tj = [-nhat(2) ; nhat(1) ; -scalRadiusJ];
        nn = n*n.';

        if (consts.boolSlipMobilized && boolSlip)
            scalKt = 0; % slipping contact: no tangential spring (Henkes et al.)
        end

        Hii = scalKn*nn + scalKt*(ti*ti.');
        Hjj = scalKn*nn + scalKt*(tj*tj.');
        Hij = scalKn*nn + scalKt*(ti*tj.');

        Gii = consts.gn*nn + consts.gt*(ti*ti.');
        Gjj = consts.gn*nn + consts.gt*(tj*tj.');
        Gij = consts.gn*nn + consts.gt*(ti*tj.');
    else
        Hii = scalKn*(nhat*nhat.');
        Hjj = Hii;
        Hij = Hii;

        Gii = consts.gn*eye(dim);
        Gjj = Gii;
        Gij = Gii;
    end

    if (consts.boolStressed)
        % Pre-stress term V'(r)/r (I - n n') with V'(r) = -f_n, the normal force
        % of the contact law.  It only involves the translational DOFs.
        matPreStress = -(scalFn/scalDist)*(eye(dim) - nhat*nhat.');
        if (consts.friction)
            matPreStress = blkdiag(matPreStress, 0);
        end
        Hii = Hii + matPreStress;
        Hjj = Hjj + matPreStress;
        Hij = Hij + matPreStress;
    end
end

function [boolSlip, boolKnown] = lookupContact(i, j, consts)
% Is contact (i, j) fully mobilized, and is it listed in fricState at all?
    if (~consts.boolUseFricState)
        boolSlip  = false;
        boolKnown = false;
        return;
    end
    a = min(i, j);
    b = max(i, j);
    boolKnown = full(consts.matInFricState(a, b));
    boolSlip  = full(consts.matMobilized(a, b));
end

function [matMobilized, matInFricState, vecMobilization] = mobilizedContacts(fricState, N, radii, scalTol)
% Sparse N x N logical maps (upper triangle, i < j) of the contacts listed in
% pack.m's fricState and of the fully mobilized ones, m = |f_t|/(mu f_n) >= 1 - tol.
    if ((fricState.N ~= N) || (max(abs(fricState.vecRadius(:) - radii)) > 1e-12*max(radii)))
        error("matSpringDampMass:FricStateMismatch", ...
            "fricState does not match the packing. Pass pack.m's frictional file (backbone) with its top-level fricState, or sFullState with sFullState.fricState.");
    end
    vecI = min(fricState.vecContactNN(:), fricState.vecContactMM(:));
    vecJ = max(fricState.vecContactNN(:), fricState.vecContactMM(:));
    if (fricState.mu > 0)
        vecMobilization = abs(fricState.vecFt(:))./(fricState.mu*abs(fricState.vecFn(:)));
    else
        vecMobilization = ones(numel(vecI), 1); % mu = 0: no contact can hold a tangential force
    end
    vecBoolMobilized = vecMobilization >= 1 - scalTol;
    matInFricState = sparse(vecI, vecJ, true, N, N);
    matMobilized   = sparse(vecI(vecBoolMobilized), vecJ(vecBoolMobilized), true, N, N);
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

