function scalMeanCoordNum = computeMeanCoordNum(vecPosX, vecPosY, vecDiameter, scalBoxWidthX, scalBoxHeightY, vecPosZ, scalBoxDepthZ)
% computeMeanCoordNum  Mean particle-particle coordination number.
%
% SINGLE SOURCE OF TRUTH for the mean-coordination-number definition
% shared by pack.m and packRepeatTile.m. A pair (i,j) is in contact when
% the minimum-image center distance is STRICTLY LESS THAN the sum of the
% radii, (d_i + d_j)/2; the mean is over ALL particles (rattlers included
% at zero). If the contact rule ever changes, change it here and both the
% packing and the tiling algorithms pick it up.
%
% Inputs:
%   vecPosX, vecPosY, vecDiameter   [N x 1] positions and diameters
%   scalBoxWidthX, scalBoxHeightY   box dimensions (Lx, Ly)
%   vecPosZ, scalBoxDepthZ          [N x 1] z-positions and Lz — omit both
%                                   for 2D packings
%
% Output:
%   scalMeanCoordNum   [1 x 1] mean per-particle contact count
%
% Vectorized O(N^2) over unique pairs (i < j); fine for the packing sizes
% used here (a cell list can replace this for much larger N).

    if nargin < 6 || isempty(vecPosZ)
        boolThreeD = false;
    else
        boolThreeD = true;
    end

    N = numel(vecPosX);
    [i, j] = ndgrid(1:N, 1:N);
    boolUpper = i < j;
    i = i(boolUpper);
    j = j(boolUpper);

    % Minimum-image separations
    vecSepX = vecPosX(j) - vecPosX(i);
    vecSepX = vecSepX - scalBoxWidthX * round(vecSepX / scalBoxWidthX);
    vecSepY = vecPosY(j) - vecPosY(i);
    vecSepY = vecSepY - scalBoxHeightY * round(vecSepY / scalBoxHeightY);
    if boolThreeD
        vecSepZ = vecPosZ(j) - vecPosZ(i);
        vecSepZ = vecSepZ - scalBoxDepthZ * round(vecSepZ / scalBoxDepthZ);
    end

    vecContactDist = (vecDiameter(i) + vecDiameter(j)) / 2;   % r_i + r_j
    vecSepDistSq   = vecSepX.^2 + vecSepY.^2;
    if boolThreeD
        vecSepDistSq = vecSepDistSq + vecSepZ.^2;
    end
    boolContact = vecSepDistSq < vecContactDist.^2;          % strictly touching

    % Per-particle contact counts; each contact contributes 1 to both
    % partners. The ~any() guard avoids accumarray's empty-input error on
    % a packing with zero contacts.
    if ~any(boolContact)
        vecCoordNum = zeros(N, 1);
    else
        vecCoordNum = accumarray(i(boolContact), 1, [N 1]) ...
                    + accumarray(j(boolContact), 1, [N 1]);
    end
    scalMeanCoordNum = mean(vecCoordNum);
end
