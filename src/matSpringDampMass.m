function [matSpring, matDamp, matMass, sContactInfo] = matSpringDampMass(sPacking, scalDampingConstant, scalSpringConstant, scalMass, sOptions)
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
% In the code, n_hat is vecNormalProjection, t_i and t_j are
% vecTangentProjectionI/J, n_hat n_hat' is matNormalProjection, H_ii, H_jj,
% H_ij are matStiffnessBlockII/JJ/IJ and the Γ blocks matDampingBlockII/JJ/IJ.
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
%     sPacking             Struct with the packing information
%     scalDampingConstant  Damping parameter (normal damping γ_n with friction)
%     scalSpringConstant   Hooke's spring constant (normal stiffness k_n with friction)
%     scalMass             Mass of a particle of unit radius
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
%                                Defaults to scalSpringConstant/3, the
%                                Cundall-Strack ratio used by pack.m
%     tangentialDampingConstant  Tangential damping γ_t (friction only).
%                                Defaults to scalDampingConstant
%     stressed           Add the pre-stress term -(f_n/r)(I - n n') to K
%     fricState          fricState struct saved by pack.m (per-contact forces)
%     slipMobilized      Set k_t = 0 on fully mobilized contacts (needs
%                        friction = true and fricState)
%     mobilizationTol    m >= 1 - mobilizationTol counts as fully mobilized
%                        (default 1e-6)
%     contactLaw         "hooke" (default) or "hertz" (Hertz-Mindlin with
%                        friction); see option 3 above
%
% The option names are the keyword interface used by the callers and the
% tests, so they keep their names; inside the function they are read from the
% struct sOptions.
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
%     [matK, matGamma, matM] = matSpringDampMass(sPacking, scalGamma, scalK, 1, ...
%         "periodic", true, "friction", true);
%     [matPsi, vecLambda] = polyeig(matK, matGamma, matM);
%     % solves (matK + lambda*matGamma + lambda^2*matM) * psi = 0
%
% Example (undamped benchmark against Henkes et al., pack.m frictional run;
% the one output file holds the backbone and its re-indexed fricState):
%     sPack = load('2D_N1000_P0.001_Width32_Seed1_Friction.mat');
%     [matK, ~, matM, sInfo] = matSpringDampMass(sPack, 0, sPack.K, 1, ...
%         "periodic", true, "uniformMass", true, "friction", true, ...
%         "stressed", true, "fricState", sPack.fricState, "slipMobilized", true);
%     vecOmega2 = eig(full(matK), full(matM));   % Γ = 0, so lambda = ±i*omega

    arguments
        sPacking                 (1, 1) struct
        scalDampingConstant      (1, 1) double {mustBeReal} = 1.0
        scalSpringConstant       (1, 1) double {mustBeReal} = 1.0
        scalMass                 (1, 1) double {mustBeReal} = 1.0
        sOptions.algorithm        (1, 1) string = "cell"
        sOptions.periodic         (1, 1) logical = false % true ==> periodic boundary condition
        sOptions.sparseOutput     (1, 1) logical = false % true ==> output sparse matrices
        sOptions.uniformMass      (1, 1) logical = false % true ==> assume all particles are of the same mass
        sOptions.oldPackingFormat (1, 1) logical = false % true ==> assume packing is stored in the old format
        sOptions.friction         (1, 1) logical = false % true ==> tangential contacts + rotation (2D only)
        sOptions.tangentialSpringConstant  (1, 1) double {mustBeReal} = NaN % k_t; NaN ==> scalSpringConstant/3
        sOptions.tangentialDampingConstant (1, 1) double {mustBeReal} = NaN % γ_t; NaN ==> scalDampingConstant
        sOptions.stressed         (1, 1) logical = false % true ==> add pre-stress term -(f_n/r)(I - n n')
        sOptions.fricState        (1, 1) struct = struct() % pack.m fricState (per-contact forces)
        sOptions.slipMobilized    (1, 1) logical = false % true ==> k_t = 0 on fully mobilized contacts
        sOptions.mobilizationTol  (1, 1) double {mustBeNonnegative} = 1e-6 % m >= 1 - tol ==> fully mobilized
        sOptions.contactLaw       (1, 1) string {mustBeMember(sOptions.contactLaw, ["hooke", "hertz"])} = "hooke"
    end

    % Grab the data out of the packing.
    [scalDim, scalNumParticles, ~, matPositions, vecRadii, vecBoxDims] = getPackingData(sPacking, sOptions.oldPackingFormat);
    vecRadii = vecRadii(:);

    % Contact constants and degrees of freedom per particle.  With friction
    % each particle also has a rotation theta, so there are scalDim + 1.
    sConsts = struct('scalKn', scalSpringConstant, 'scalKt', 0, 'scalGn', scalDampingConstant, 'scalGt', 0, ...
        'boolFriction', sOptions.friction);
    if (sOptions.friction)
        if (scalDim ~= 2)
            error("matSpringDampMass:Friction3DNotSupported", "friction = true is 2D only.");
        end
        scalDofPerParticle = scalDim + 1;
        sConsts.scalKt = sOptions.tangentialSpringConstant;
        if (isnan(sConsts.scalKt))
            sConsts.scalKt = scalSpringConstant/3;
        end
        sConsts.scalGt = sOptions.tangentialDampingConstant;
        if (isnan(sConsts.scalGt))
            sConsts.scalGt = scalDampingConstant;
        end
    else
        scalDofPerParticle = scalDim;
    end
    scalNumDof = scalDofPerParticle*scalNumParticles;

    % Pre-stress term, slipping (fully mobilized) contacts and contact law.
    sConsts.boolStressed = sOptions.stressed;
    sConsts.boolSlipMobilized = sOptions.slipMobilized;
    sConsts.boolUseFricState = ~isempty(fieldnames(sOptions.fricState));
    sConsts.boolHertz = (sOptions.contactLaw == "hertz");
    sConsts.scalKtOverKn = sConsts.scalKt/sConsts.scalKn;   % Mindlin: k_t = (k_t0/K) k_n(delta)
    if (sOptions.slipMobilized && ~(sOptions.friction && sConsts.boolUseFricState))
        error("matSpringDampMass:SlipNeedsFricState", ...
            "slipMobilized = true needs friction = true and the fricState saved by pack.m.");
    end
    vecMobilization = [];
    if (sConsts.boolUseFricState)
        [sConsts.matMobilized, sConsts.matInFricState, vecMobilization] = ...
            mobilizedContacts(sOptions.fricState, scalNumParticles, vecRadii, sOptions.mobilizationTol);
    end

    % Which particles touch the left and right walls (2D)?
    if ((scalDim == 2) && ~sOptions.periodic)
        vecBoolLeftWall  = matPositions(:, 1) < vecRadii;
        vecBoolRightWall = matPositions(:, 1) > vecBoxDims(1) - vecRadii;
    end

    % Build the stiffness and damping matrices.
    switch (sOptions.algorithm)
    case "cell"
        [matSpring, matDamp, sCounts] = buildMatricesCell(scalNumParticles, scalDim, scalDofPerParticle, scalNumDof, ...
            matPositions, vecRadii, vecBoxDims, sConsts, sOptions);
    case "particle"
        [matSpring, matDamp, sCounts] = buildMatricesParticle(scalNumParticles, scalDim, scalDofPerParticle, scalNumDof, ...
            matPositions, vecRadii, vecBoxDims, sConsts, sOptions);
    otherwise
        error("Not implemented.");
    end

    % Add entries for wall contacts (2D, translational DOFs only).
    if ((scalDim == 2) && ~sOptions.periodic)
        for (idxParticle = 1:1:scalNumParticles)
            if (vecBoolLeftWall(idxParticle) || vecBoolRightWall(idxParticle))
                vecDofsParticle = (scalDofPerParticle*(idxParticle - 1) + 1):(scalDofPerParticle*(idxParticle - 1) + scalDim);
                matSpring(vecDofsParticle, vecDofsParticle) = matSpring(vecDofsParticle, vecDofsParticle) + scalSpringConstant*eye(scalDim);
            end
        end
    end

    % Build the mass matrix.
    if (sOptions.uniformMass)
        vecMass = scalMass*ones(scalNumParticles, 1);
    else
        scalUnitBallVolume = pi^(scalDim/2)/gamma(scalDim/2 + 1); % Volume of unit ball in R^dim
        vecMass = (scalMass*scalUnitBallVolume)*vecRadii.^scalDim;
    end
    if (sOptions.friction)
        vecInertia = 0.5*vecMass.*vecRadii.^2; % Moment of inertia of a uniform disk
        vecMassDiag = reshape([vecMass.' ; vecMass.' ; vecInertia.'], scalNumDof, 1);
    else
        vecMassDiag = reshape(repmat(vecMass.', scalDim, 1), scalNumDof, 1);
    end
    if (sOptions.sparseOutput)
        matMass = spdiags(vecMassDiag, 0, scalNumDof, scalNumDof);
    else
        matMass = diag(vecMassDiag);
    end

    % Summary of the contact network K was built from.  The field names are
    % read by runFrictionalDOS, runUnresolvedDOS and the tests.
    sContactInfo = struct();
    sContactInfo.numContacts = sum(sCounts.vecNumContacts)/2;
    sContactInfo.vecNumContacts = sCounts.vecNumContacts;
    if (sConsts.boolUseFricState)
        sContactInfo.numMobilized = sCounts.scalNumMobilized;
        sContactInfo.numNotInFricState = sCounts.scalNumNotInFricState;
        if (sCounts.scalNumNotInFricState > 0)
            warning("matSpringDampMass:ContactsNotInFricState", ...
                "%d contacts are not in fricState and are treated as not mobilized.", sCounts.scalNumNotInFricState);
        end
    else
        sContactInfo.numMobilized = NaN;
        sContactInfo.numNotInFricState = NaN;
    end
    sContactInfo.vecMobilization = vecMobilization;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [matSpring, matDamp] = allocateMatrices(scalNumDof, boolSparse)
    if (boolSparse)
        % TODO:  Do these over-allocate?
        matSpring = spalloc(scalNumDof, scalNumDof, ceil(0.005*(scalNumDof^2)));
        matDamp   = spalloc(scalNumDof, scalNumDof, ceil(0.005*(scalNumDof^2)));
    else
        matSpring = zeros(scalNumDof, scalNumDof);
        matDamp   = zeros(scalNumDof, scalNumDof);
    end
end

function [matSpring, matDamp, sCounts] = buildMatricesCell(scalNumParticles, scalDim, scalDofPerParticle, scalNumDof, ...
        matPositions, vecRadii, vecBoxDims, sConsts, sOptions)
    % Pre-allocate the outputs.
    [matSpring, matDamp] = allocateMatrices(scalNumDof, sOptions.sparseOutput);
    sCounts = struct('vecNumContacts', zeros(scalNumParticles, 1), 'scalNumMobilized', 0, 'scalNumNotInFricState', 0);

    % Data for the cell-based method for checking particle contacts.
    % The cells must tile the box exactly (vecCellCounts.*vecCellWidths == vecBoxDims)
    % so that the periodic cell wrap matches the minimum-image wrap. Recomputing
    % the count as ceil(vecBoxDims./vecCellWidths) is a floating-point trap:
    % vecBoxDims./(vecBoxDims./k) can come out as k + 2e-15, which ceil turns into
    % k + 1 cells that no longer tile the box, and contacts across the boundary
    % are dropped.
    scalBaseCellWidth = 4.0*max(vecRadii);                              % Nominal width of each cell
    vecCellCounts     = max(round(vecBoxDims/scalBaseCellWidth), 1);    % Number of cells in each dimension
    vecCellWidths     = vecBoxDims./vecCellCounts;                      % Actual cell widths (tile the box)
    matParticleToCell = min(max(ceil(matPositions./vecCellWidths.'), 1), vecCellCounts.'); % row p: cell subscripts of particle p
    cellCellToParticles = cell(vecCellCounts.');                        % (a, b, c) entry lists the particles in cell (a, b, c)

    % Build the map from cell subscripts to the list of particles each cell contains.
    cellCellSubscript = cell(1, scalDim);
    for (idxCell = 1:prod(vecCellCounts))
        [cellCellSubscript{:}] = ind2sub(vecCellCounts, idxCell);
        cellCellToParticles{cellCellSubscript{:}} = find(all(matParticleToCell == [cellCellSubscript{:}], 2)).';
    end

    % Loop over all cells looking for contacts between particles in this cell
    % and in adjacent ones.
    for (idxCell = 1:prod(vecCellCounts))
        [cellCellSubscript{:}] = ind2sub(vecCellCounts, idxCell);

        % Get a list of all the particles that are in either this cell or a
        % neighboring one.  The call to unique() eliminates duplicates in the
        % event that the cell count is very low.
        %
        % TODO:  This is not very efficient.
        cellNeighborLists = cell(1, 3^scalDim);
        cellStencilSubscript = cell(1, scalDim);
        for (idxStencil = 1:1:(3^scalDim))
            [cellStencilSubscript{:}] = ind2sub(3*ones(1, scalDim), idxStencil);
            vecStencilOffset = cell2mat(cellStencilSubscript) - 2;

            if (sOptions.periodic)
                cellNeighborCellSubscript = mat2cell(mod(cell2mat(cellCellSubscript) + vecStencilOffset - 1, vecCellCounts.') + 1, ...
                    1, ones(1, scalDim));
            else
                error('Not implemented.');
            end

            cellNeighborLists(idxStencil) = cellCellToParticles(cellNeighborCellSubscript{:});
        end
        vecNeighborParticles = unique(horzcat(cellNeighborLists{:}));

        % Check contacts for all particles in this cell.
        for (idxParticleI = cellCellToParticles{idxCell})
            for (idxParticleJ = vecNeighborParticles)
                if (idxParticleJ == idxParticleI)
                    continue;
                end

                % Separation vector from particle j to particle i.
                vecSeparation = (matPositions(idxParticleI, :) - matPositions(idxParticleJ, :)).';
                if (sOptions.periodic)
                    vecSeparation = vecSeparation - round(vecSeparation./vecBoxDims).*vecBoxDims;
                end
                scalCenterDistance = norm(vecSeparation);

                % If the particles are in contact, add terms for the contact forces.
                if (scalCenterDistance < vecRadii(idxParticleI) + vecRadii(idxParticleJ))
                    vecUnitNormal = -vecSeparation/scalCenterDistance; % Unit normal pointing from i to j
                    [boolSlip, boolKnown] = lookupContact(idxParticleI, idxParticleJ, sConsts);
                    [matStiffnessBlockII, matStiffnessBlockJJ, matStiffnessBlockIJ, ...
                        matDampingBlockII, matDampingBlockJJ, matDampingBlockIJ] = contactBlocks(vecUnitNormal, ...
                        vecRadii(idxParticleI), vecRadii(idxParticleJ), scalCenterDistance, scalDim, sConsts, boolSlip);

                    vecDofsI = (scalDofPerParticle*(idxParticleI - 1) + 1):(scalDofPerParticle*idxParticleI);
                    vecDofsJ = (scalDofPerParticle*(idxParticleJ - 1) + 1):(scalDofPerParticle*idxParticleJ);

                    matSpring(vecDofsI, vecDofsI) = matSpring(vecDofsI, vecDofsI) + matStiffnessBlockII;
                    matSpring(vecDofsJ, vecDofsJ) = matSpring(vecDofsJ, vecDofsJ) + matStiffnessBlockJJ;
                    matSpring(vecDofsI, vecDofsJ) = matSpring(vecDofsI, vecDofsJ) - matStiffnessBlockIJ;
                    matSpring(vecDofsJ, vecDofsI) = matSpring(vecDofsJ, vecDofsI) - matStiffnessBlockIJ.';

                    matDamp(vecDofsI, vecDofsI) = matDamp(vecDofsI, vecDofsI) + matDampingBlockII;
                    matDamp(vecDofsJ, vecDofsJ) = matDamp(vecDofsJ, vecDofsJ) + matDampingBlockJJ;
                    matDamp(vecDofsI, vecDofsJ) = matDamp(vecDofsI, vecDofsJ) - matDampingBlockIJ;
                    matDamp(vecDofsJ, vecDofsI) = matDamp(vecDofsJ, vecDofsI) - matDampingBlockIJ.';

                    % Each contact is visited once from each side, so only
                    % particle i is counted here and the totals are halved.
                    sCounts.vecNumContacts(idxParticleI) = sCounts.vecNumContacts(idxParticleI) + 1;
                    sCounts.scalNumMobilized      = sCounts.scalNumMobilized + boolSlip/2;
                    sCounts.scalNumNotInFricState = sCounts.scalNumNotInFricState + (sConsts.boolUseFricState && ~boolKnown)/2;
                end
            end
        end
    end

    % This method double-counts contacts, so divide by 2.
    matSpring = matSpring/2.0;
    matDamp = matDamp/2.0;
end

function [matSpring, matDamp, sCounts] = buildMatricesParticle(scalNumParticles, scalDim, scalDofPerParticle, scalNumDof, ...
        matPositions, vecRadii, vecBoxDims, sConsts, sOptions)
    % Pre-allocate the outputs.
    [matSpring, matDamp] = allocateMatrices(scalNumDof, sOptions.sparseOutput);
    sCounts = struct('vecNumContacts', zeros(scalNumParticles, 1), 'scalNumMobilized', 0, 'scalNumNotInFricState', 0);

    % Check for contact between each pair of particles to generate the
    % stiffness and damping matrices.
    for (idxParticleI = 1:1:scalNumParticles)
        for (idxParticleJ = (idxParticleI + 1):scalNumParticles)
            % Separation vector from particle j to particle i.
            vecSeparation = (matPositions(idxParticleI, :) - matPositions(idxParticleJ, :)).';
            if (sOptions.periodic)
                vecSeparation = vecSeparation - round(vecSeparation./vecBoxDims).*vecBoxDims;
            end
            scalCenterDistance = norm(vecSeparation);

            % If the particles are in contact, add terms for the contact forces.
            if (scalCenterDistance < vecRadii(idxParticleI) + vecRadii(idxParticleJ))
                vecUnitNormal = -vecSeparation/scalCenterDistance; % Unit normal pointing from i to j
                [boolSlip, boolKnown] = lookupContact(idxParticleI, idxParticleJ, sConsts);
                [matStiffnessBlockII, matStiffnessBlockJJ, matStiffnessBlockIJ, ...
                    matDampingBlockII, matDampingBlockJJ, matDampingBlockIJ] = contactBlocks(vecUnitNormal, ...
                    vecRadii(idxParticleI), vecRadii(idxParticleJ), scalCenterDistance, scalDim, sConsts, boolSlip);

                vecDofsI = (scalDofPerParticle*(idxParticleI - 1) + 1):(scalDofPerParticle*idxParticleI);
                vecDofsJ = (scalDofPerParticle*(idxParticleJ - 1) + 1):(scalDofPerParticle*idxParticleJ);

                matSpring(vecDofsI, vecDofsI) = matSpring(vecDofsI, vecDofsI) + matStiffnessBlockII;
                matSpring(vecDofsJ, vecDofsJ) = matSpring(vecDofsJ, vecDofsJ) + matStiffnessBlockJJ;
                matSpring(vecDofsI, vecDofsJ) = matSpring(vecDofsI, vecDofsJ) - matStiffnessBlockIJ;
                matSpring(vecDofsJ, vecDofsI) = matSpring(vecDofsJ, vecDofsI) - matStiffnessBlockIJ.';

                matDamp(vecDofsI, vecDofsI) = matDamp(vecDofsI, vecDofsI) + matDampingBlockII;
                matDamp(vecDofsJ, vecDofsJ) = matDamp(vecDofsJ, vecDofsJ) + matDampingBlockJJ;
                matDamp(vecDofsI, vecDofsJ) = matDamp(vecDofsI, vecDofsJ) - matDampingBlockIJ;
                matDamp(vecDofsJ, vecDofsI) = matDamp(vecDofsJ, vecDofsI) - matDampingBlockIJ.';

                sCounts.vecNumContacts(idxParticleI) = sCounts.vecNumContacts(idxParticleI) + 1;
                sCounts.vecNumContacts(idxParticleJ) = sCounts.vecNumContacts(idxParticleJ) + 1;
                sCounts.scalNumMobilized      = sCounts.scalNumMobilized + boolSlip;
                sCounts.scalNumNotInFricState = sCounts.scalNumNotInFricState + (sConsts.boolUseFricState && ~boolKnown);
            end
        end
    end
end

function [matStiffnessBlockII, matStiffnessBlockJJ, matStiffnessBlockIJ, matDampingBlockII, matDampingBlockJJ, matDampingBlockIJ] = ...
        contactBlocks(vecUnitNormal, scalRadiusI, scalRadiusJ, scalCenterDistance, scalDim, sConsts, boolSlip)
% Stiffness and damping blocks of one contact (i, j), with vecUnitNormal the
% unit normal pointing from i to j and scalCenterDistance the center distance.
% The contact contributes
%     [ matStiffnessBlockII  -matStiffnessBlockIJ ; -matStiffnessBlockIJ'  matStiffnessBlockJJ ]
% to the (i, j) rows/columns of K, and likewise for Γ.  boolSlip marks a fully
% mobilized contact, which loses its tangential spring when slipMobilized is on.

    % Contact stiffnesses and normal force of the contact law.
    if (sConsts.boolHertz)
        % Hertz(-Mindlin), as in pack.m: f_n = (4/3) K sqrt(R_eff) delta^(3/2),
        % k_n = 2 K sqrt(R_eff delta), k_t = (k_t0/K) k_n.
        scalOverlap   = scalRadiusI + scalRadiusJ - scalCenterDistance;
        scalRadiusEff = scalRadiusI*scalRadiusJ/(scalRadiusI + scalRadiusJ);
        scalKn = 2*sConsts.scalKn*sqrt(scalRadiusEff*scalOverlap);
        scalFn = (4/3)*sConsts.scalKn*sqrt(scalRadiusEff)*scalOverlap^(3/2);
        scalKt = sConsts.scalKtOverKn*scalKn;
    else
        scalKn = sConsts.scalKn;
        scalFn = sConsts.scalKn*(scalRadiusI + scalRadiusJ - scalCenterDistance);
        scalKt = sConsts.scalKt;
    end

    if (sConsts.boolFriction)
        % Projection vectors on (x, y, theta).  t = z x n, and the third
        % entries are the lever arms that turn sliding into spin.
        vecNormalProjection    = [vecUnitNormal ; 0];
        vecTangentProjectionI  = [-vecUnitNormal(2) ; vecUnitNormal(1) ;  scalRadiusI];
        vecTangentProjectionJ  = [-vecUnitNormal(2) ; vecUnitNormal(1) ; -scalRadiusJ];
        matNormalProjection    = vecNormalProjection*vecNormalProjection.';

        if (sConsts.boolSlipMobilized && boolSlip)
            scalKt = 0; % slipping contact: no tangential spring (Henkes et al.)
        end

        matStiffnessBlockII = scalKn*matNormalProjection + scalKt*(vecTangentProjectionI*vecTangentProjectionI.');
        matStiffnessBlockJJ = scalKn*matNormalProjection + scalKt*(vecTangentProjectionJ*vecTangentProjectionJ.');
        matStiffnessBlockIJ = scalKn*matNormalProjection + scalKt*(vecTangentProjectionI*vecTangentProjectionJ.');

        matDampingBlockII = sConsts.scalGn*matNormalProjection + sConsts.scalGt*(vecTangentProjectionI*vecTangentProjectionI.');
        matDampingBlockJJ = sConsts.scalGn*matNormalProjection + sConsts.scalGt*(vecTangentProjectionJ*vecTangentProjectionJ.');
        matDampingBlockIJ = sConsts.scalGn*matNormalProjection + sConsts.scalGt*(vecTangentProjectionI*vecTangentProjectionJ.');
    else
        matStiffnessBlockII = scalKn*(vecUnitNormal*vecUnitNormal.');
        matStiffnessBlockJJ = matStiffnessBlockII;
        matStiffnessBlockIJ = matStiffnessBlockII;

        matDampingBlockII = sConsts.scalGn*eye(scalDim);
        matDampingBlockJJ = matDampingBlockII;
        matDampingBlockIJ = matDampingBlockII;
    end

    if (sConsts.boolStressed)
        % Pre-stress term V'(r)/r (I - n n') with V'(r) = -f_n, the normal force
        % of the contact law.  It only involves the translational DOFs.
        matPreStress = -(scalFn/scalCenterDistance)*(eye(scalDim) - vecUnitNormal*vecUnitNormal.');
        if (sConsts.boolFriction)
            matPreStress = blkdiag(matPreStress, 0);
        end
        matStiffnessBlockII = matStiffnessBlockII + matPreStress;
        matStiffnessBlockJJ = matStiffnessBlockJJ + matPreStress;
        matStiffnessBlockIJ = matStiffnessBlockIJ + matPreStress;
    end
end

function [boolSlip, boolKnown] = lookupContact(idxParticleI, idxParticleJ, sConsts)
% Is contact (i, j) fully mobilized, and is it listed in fricState at all?
    if (~sConsts.boolUseFricState)
        boolSlip  = false;
        boolKnown = false;
        return;
    end
    idxLow  = min(idxParticleI, idxParticleJ);
    idxHigh = max(idxParticleI, idxParticleJ);
    boolKnown = full(sConsts.matInFricState(idxLow, idxHigh));
    boolSlip  = full(sConsts.matMobilized(idxLow, idxHigh));
end

function [matMobilized, matInFricState, vecMobilization] = mobilizedContacts(sFricState, scalNumParticles, vecRadii, scalTol)
% Sparse N x N logical maps (upper triangle, i < j) of the contacts listed in
% pack.m's fricState and of the fully mobilized ones, m = |f_t|/(mu f_n) >= 1 - tol.
    if ((sFricState.N ~= scalNumParticles) || (max(abs(sFricState.vecRadius(:) - vecRadii)) > 1e-12*max(vecRadii)))
        error("matSpringDampMass:FricStateMismatch", ...
            "fricState does not match the packing. Pass pack.m's frictional file (backbone) with its top-level fricState, or sFullState with sFullState.fricState.");
    end
    vecI = min(sFricState.vecContactNN(:), sFricState.vecContactMM(:));
    vecJ = max(sFricState.vecContactNN(:), sFricState.vecContactMM(:));
    if (sFricState.mu > 0)
        vecMobilization = abs(sFricState.vecFt(:))./(sFricState.mu*abs(sFricState.vecFn(:)));
    else
        vecMobilization = ones(numel(vecI), 1); % mu = 0: no contact can hold a tangential force
    end
    vecBoolMobilized = vecMobilization >= 1 - scalTol;
    matInFricState = sparse(vecI, vecJ, true, scalNumParticles, scalNumParticles);
    matMobilized   = sparse(vecI(vecBoolMobilized), vecJ(vecBoolMobilized), true, scalNumParticles, scalNumParticles);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [scalDim, scalNumParticles, scalNumDof, matPositions, vecRadii, vecBoxDims] = getPackingData(sPacking, boolOldPackingFormat)
    if (boolOldPackingFormat)          % Old packing format
        if (isfield(sPacking, 'z'))
            scalDim = 3;
            matPositions = [sPacking.x.' sPacking.y.' sPacking.z.'];
            vecBoxDims   = [sPacking.Lx ; sPacking.Ly ; sPacking.Lz];
        else
            scalDim = 2;
            matPositions = [sPacking.x.' sPacking.y.'];
            vecBoxDims   = [sPacking.Lx ; sPacking.Ly];
        end
        vecRadii = sPacking.Dn.'/2.0;
    else                               % New packing format
        if (isfield(sPacking, 'vecPosZ'))
            scalDim = 3;
            matPositions = [sPacking.vecPosX sPacking.vecPosY sPacking.vecPosZ];
            vecBoxDims   = [sPacking.scalBoxWidthX ; sPacking.scalBoxHeightY ; sPacking.scalBoxDepthZ];
        else
            scalDim = 2;
            matPositions = [sPacking.vecPosX sPacking.vecPosY];
            vecBoxDims   = [sPacking.scalBoxWidthX ; sPacking.scalBoxHeightY];
        end
        vecRadii = sPacking.vecDiameter/2.0;
    end

    scalNumParticles = size(matPositions, 1);
    scalNumDof = scalNumParticles*scalDim;
end
