function scalMeanCoordNum = computeMeanCoordNum(vecPosX, vecPosY, vecDiameter, scalBoxWidthX, scalBoxHeightY, vecPosZ, scalBoxDepthZ)
% computeMeanCoordNum  Mean particle-particle coordination number.
%
% Thin wrapper over computeCoordNum so pack.m and packRepeatTile.m (and
% any future packing code) share ONE definition of the mean coordination
% number:  mean of the per-particle contact counts, all particles
% included. The contact convention itself lives in computeCoordNum —
% update it there and both algorithms follow.
%
% Arguments: same as computeCoordNum. Omit vecPosZ/scalBoxDepthZ for 2D.

    if nargin < 6 || isempty(vecPosZ)
        vecCoordNum = computeCoordNum(vecPosX, vecPosY, vecDiameter, scalBoxWidthX, scalBoxHeightY);
    else
        vecCoordNum = computeCoordNum(vecPosX, vecPosY, vecDiameter, scalBoxWidthX, scalBoxHeightY, vecPosZ, scalBoxDepthZ);
    end
    scalMeanCoordNum = mean(vecCoordNum);
end
