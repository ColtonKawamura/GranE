function scalPackingFraction = computePackingFraction(vecDiameter, scalBoxWidthX, scalBoxHeightY, scalBoxDepthZ)
% computePackingFraction  Packing fraction = total particle volume / box
% volume.
%
% SINGLE SOURCE OF TRUTH for the packing-fraction definition shared by
% pack.m and packRepeatTile.m. Update the definition here (e.g. a new
% particle shape) and both algorithms follow.
%
% Inputs:
%   vecDiameter        [N x 1] particle diameters
%   scalBoxWidthX, scalBoxHeightY   box dimensions (Lx, Ly)
%   scalBoxDepthZ      Lz — omit for 2D packings
%
% 2D (no Lz): disk area pi*(d/2)^2 over box area Lx*Ly.
% 3D (Lz given): sphere volume (4/3)*pi*(d/2)^3 over box volume Lx*Ly*Lz.

    if nargin < 4 || isempty(scalBoxDepthZ)
        scalVolumeSpheres = sum(pi * (vecDiameter/2).^2);
        scalVolumeBox     = scalBoxWidthX * scalBoxHeightY;
    else
        scalVolumeSpheres = sum((4/3) * pi * (vecDiameter/2).^3);
        scalVolumeBox     = scalBoxWidthX * scalBoxHeightY * scalBoxDepthZ;
    end
    scalPackingFraction = scalVolumeSpheres / scalVolumeBox;
end
