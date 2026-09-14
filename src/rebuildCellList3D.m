% rebuildCellList3D -- Re-bucket 3D particle positions into a periodic cell list.
% Returns the wrapped positions, the cell-of-particles cell array, and the
% number of cells per axis.
%
% Extracted from pack.m (was a local subfunction) so it can be unit-tested
% directly with synthetic positions.

function [vecPosX, vecPosY, vecPosZ, cellParticleList, scalNumCellsX, scalNumCellsY, scalNumCellsZ] = ...
        rebuildCellList3D(vecPosX, vecPosY, vecPosZ, scalBoxWidthX, scalBoxHeightY, scalBoxDepthZ, scalRawCellWidth, scalTimestep, N)

    scalNumCellsX  = round(scalBoxWidthX  / scalRawCellWidth);
    scalCellWidthX = scalBoxWidthX  / scalNumCellsX;
    scalNumCellsY  = round(scalBoxHeightY / scalRawCellWidth);
    scalCellWidthY = scalBoxHeightY / scalNumCellsY;
    scalNumCellsZ  = round(scalBoxDepthZ  / scalRawCellWidth);
    scalCellWidthZ = scalBoxDepthZ  / scalNumCellsZ;

    % Sanity check
    vecFloorX = floor(vecPosX / scalBoxWidthX);
    vecFloorY = floor(vecPosY / scalBoxHeightY);
    vecFloorZ = floor(vecPosZ / scalBoxDepthZ);
    if any(abs(vecFloorX) > 1) || any(abs(vecFloorY) > 1) || any(abs(vecFloorZ) > 1)
        error(['Particle moved more than one box length in a single timestep.\n' ...
               'Reduce scalTimestep (currently %.4e).'], scalTimestep);
    end

    % Wrap into [0, L)
    vecPosX = mod(vecPosX, scalBoxWidthX);
    vecPosY = mod(vecPosY, scalBoxHeightY);
    vecPosZ = mod(vecPosZ, scalBoxDepthZ);

    % Map to cell indices
    vecCellIdxX = min(max(ceil(vecPosX / scalCellWidthX), 1), scalNumCellsX);
    vecCellIdxY = min(max(ceil(vecPosY / scalCellWidthY), 1), scalNumCellsY);
    vecCellIdxZ = min(max(ceil(vecPosZ / scalCellWidthZ), 1), scalNumCellsZ);

    % Linear index: x + Nx*(y-1) + Nx*Ny*(z-1)
    vecCellLinearIdx = vecCellIdxX + scalNumCellsX*(vecCellIdxY-1) + scalNumCellsX*scalNumCellsY*(vecCellIdxZ-1);
    cellParticleList = accumarray(vecCellLinearIdx, (1:N)', [scalNumCellsX*scalNumCellsY*scalNumCellsZ 1], @(x){x});
    cellParticleList = reshape(cellParticleList, scalNumCellsX, scalNumCellsY, scalNumCellsZ);
end
