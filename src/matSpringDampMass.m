function [matSpring, matDamp, matMass] = matSpringDampMass(positions, radii, Ly, Lx, Lz, damping_constant, springConstant, options)

    % example script:

    % load("/Users/coltonkawamura/repos/GranMA/in/3d/15by15by15/3D_N3375_P0.1_Width15_Seed1.mat")
    %
    % positions = [x(:), y(:), z(:)]; % N×3
    % radii = Dn(:)' / 2; % 1×N  (radii, not diameters)
    %
    %
    % [matSpring, matDamp, matMass] = matSpringDampMass( ...
    % positions, radii, Ly, Lx, Lz, 0, K, ...
    % 'periodic', true, 'threeD', true);

    arguments
        positions (:,:) double {mustBeReal}
        radii (1,:) double {mustBeReal}
        Ly (1,1) double {mustBeReal} = 1
        Lx (1,1) double {mustBeReal} = 1
        Lz (1,1) double {mustBeReal} = 1
        damping_constant (1,1) double {mustBeReal} = 1
        springConstant (1,1) double {mustBeReal} = 1
        options.periodic (1,1) logical = false
        options.threeD (1,1) logical = false
    end

    k = springConstant;
    if options.threeD && size(positions, 2) ~= 3
        error('options.threeD is true but positions have %d columns, expected 3.', size(positions, 2))
    end
    dim = size(positions, 2);
    mass = 1;
    N = size(positions,1);
    dof = dim * N;
    Zn = zeros(N,1);
    left_wall_list = (positions(:,1)<radii);
    right_wall_list = (positions(:,1)>Lx-radii);
    Zn(left_wall_list|right_wall_list) = 2;
    matSpring = zeros(dof, dof);
    matDamp = zeros(dof, dof);

    % Loop over all pairs of particles
    for i = 1:N
        for j = i+1:N

            % distnaces 
            dx = positions(i, 1) - positions(j, 1);
            dy = positions(i, 2) - positions(j, 2);
            dz = 0;
            if options.threeD
                 dz = positions(i, 3) - positions(j, 3);
                 if options.periodic
                     dz = dz  - round(dz/Lz)*Lz;
                 end
            end

            % Aapply periodic boundary conditions
            if options.periodic
                dy = dy - round(dy / Ly) * Ly;
                dx= dx - round(dx / Lx) * Lx;
            end
            r = sqrt(dx^2 + dy^2 +dz^2); 
            
            overlap = radii(i) + radii(j) - r;
            
            if overlap > 0 

                % keep tally of contacts for each partcile
                Zn(i) = Zn(i)+1;
                Zn(j) = Zn(j)+1;

                nvec = [dx; dy; dz] / r; % unit vector pointing from particle j to i for a signal contact pair
                nvec = nvec(1:dim); % for dz, just cuts it off so 3D and 2D are the same

                Kblock = k * (nvec * nvec'); % dim x dim, this is the block matrix
                Dblock = eye(dim);  % dim x dim (scaled later)

                % each row pair (and column pair) correspond to a particles x and y
                %  particle 1 → dofs_i = (2*0+1):(2*1) = 1:2   (rows 1,2)
                %   particle 2 → dofs_j = (2*1+1):(2*2) = 3:4   (rows 3,4)
                dofs_i = (dim*(i-1)+1):(dim*i);
                dofs_j = (dim*(j-1)+1):(dim*j);

                %  this assembles the laplacian
                %  replaces rows and columns dofs_i and dof_j with Kblock and Dblock
                matSpring(dofs_i, dofs_i) = matSpring(dofs_i, dofs_i) + Kblock; % Diagonal, coordination number (+1) for each contact for this particle
                matSpring(dofs_j, dofs_j) = matSpring(dofs_j, dofs_j) + Kblock; % same as above, but for other particle
                matSpring(dofs_i, dofs_j) = matSpring(dofs_i, dofs_j) - Kblock; %  gets a (-1) for each particle in contact
                matSpring(dofs_j, dofs_i) = matSpring(dofs_j, dofs_i) - Kblock; % same as above, but for particle j

                matDamp(dofs_i, dofs_i) = matDamp(dofs_i, dofs_i) + Dblock;
                matDamp(dofs_j, dofs_j) = matDamp(dofs_j, dofs_j) + Dblock;
                matDamp(dofs_i, dofs_j) = matDamp(dofs_i, dofs_j) - Dblock;
                matDamp(dofs_j, dofs_i) = matDamp(dofs_j, dofs_i) - Dblock;
            end
        end
    end
        if ~options.periodic && ~options.threeD
            for i = 1:N
                if left_wall_list(i) || right_wall_list(i)
                    dofs_i = (dim*(i-1)+1):(dim*i);
                    matSpring(dofs_i, dofs_i) = matSpring(dofs_i, dofs_i) + k * eye(dim);
                end
            end
        end
    matMass = mass * eye(dof);
    matDamp = damping_constant * matDamp;
end


