function [vecPosXFinal, vecPosYFinal, vecPosZFinal, vecDiameterFinal, ...
          scalBoxWidthXTiled, scalBoxHeightYFinal, scalBoxDepthZFinal, ...
          scalNTiled] = tile3D(vecPosX, vecPosY, vecPosZ, vecDiameter, ...
          scalBoxWidthX, scalBoxHeightY, scalBoxDepthZ, ...
          x_mult, y_mult, z_mult, N)
% tile3D  Repeat a 3D packing in x, y, and/or z to build a superlattice.
%
% Makes x_mult copies shifted along x, y_mult along y, and z_mult along z, so
% the output is an x_mult-by-y_mult-by-z_mult superlattice of the base tile.
% Mirrors the 2D packRepeatTile and GranMA's pack3dRepeatTile tiling.
%
% The base-tile positions must already be wrapped into [0, L) on every axis
% (pack() wraps positions each step), so the shifted copies sit in
% [i*L, (i+1)*L) and the periodic structure is reproduced exactly across tile
% boundaries: a particle crossing a tile face continues in the neighboring
% copy at the same phase, and every base contact (minimum-image over the base
% box) is reproduced in the superlattice. The packing fraction is invariant by
% construction (particle volume and box volume scale by the same factor), and
% the mean coordination number under full PBC over the superlattice box equals
% the base tile's (tiling replicates the contact network).
%
% Inputs:
%   vecPosX, vecPosY, vecPosZ   [N x 1] base-tile positions (in [0, L))
%   vecDiameter                 [N x 1] base-tile diameters
%   scalBoxWidthX/HeightY/DepthZ  base box dimensions (Lx, Ly, Lz)
%   x_mult, y_mult, z_mult      tile counts along each axis (>= 1)
%   N                           number of particles in the base tile
%
% Outputs:
%   vecPosXFinal / vecPosYFinal / vecPosZFinal  [scalNTiled x 1] tiled positions
%   vecDiameterFinal                            [scalNTiled x 1] tiled diameters
%   scalBoxWidthXTiled / scalBoxHeightYFinal / scalBoxDepthZFinal  superlattice box
%   scalNTiled                                       total particle count
%
% Example (stack a 1000-particle tile 9 times in z -> 9000 particles):
%   [x,y,z,D,Lx,Ly,Lz,Nt] = tile3D(x,y,z,D,Lx,Ly,Lz, 1, 1, 9, 1000);

    if numel(vecPosX) ~= N || numel(vecPosY) ~= N || numel(vecPosZ) ~= N || numel(vecDiameter) ~= N
        error('tile3D:SizeMismatch', ...
            ['tile3D: each position/diameter vector must have N=%d elements; ', ...
             'got x=%d, y=%d, z=%d, D=%d.'], ...
            N, numel(vecPosX), numel(vecPosY), numel(vecPosZ), numel(vecDiameter));
    end
    scalNTiled = N * x_mult * y_mult * z_mult;

    vecPosXFinal     = zeros(scalNTiled, 1);
    vecPosYFinal     = zeros(scalNTiled, 1);
    vecPosZFinal     = zeros(scalNTiled, 1);
    vecDiameterFinal = zeros(scalNTiled, 1);

    % Fill the superlattice: for each (z, y, x) tile offset, append one shifted
    % copy of the base tile. The particle ordering in the array is immaterial
    % to the physics and to the metrics (PF, Zn).
    idxEnd = 0;
    for iz = 0:z_mult-1
        for iy = 0:y_mult-1
            for ix = 0:x_mult-1
                idxStart = idxEnd + 1;
                idxEnd   = idxStart + N - 1;
                vecPosXFinal(idxStart:idxEnd)    = vecPosX + ix * scalBoxWidthX;
                vecPosYFinal(idxStart:idxEnd)    = vecPosY + iy * scalBoxHeightY;
                vecPosZFinal(idxStart:idxEnd)    = vecPosZ + iz * scalBoxDepthZ;
                vecDiameterFinal(idxStart:idxEnd) = vecDiameter;
            end
        end
    end

    scalBoxWidthXTiled  = scalBoxWidthX * x_mult;
    scalBoxHeightYFinal = scalBoxHeightY * y_mult;
    scalBoxDepthZFinal  = scalBoxDepthZ * z_mult;
end
