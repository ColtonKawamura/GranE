function scalOmegaStar = calcOmegaStar(vecBinCenters, vecDOS, scalOmegaRef, scalThreshold)
% calcOmegaStar  Crossover frequency omega* below which the DOS plateau is lost.
%
%   scalOmegaStar = calcOmegaStar(vecBinCenters, vecDOS, scalOmegaRef)
%   scalOmegaStar = calcOmegaStar(vecBinCenters, vecDOS, scalOmegaRef, scalThreshold)
%
%   Near jamming the DOS has a plateau down to a crossover frequency omega*,
%   below which it falls towards the Debye form D ~ omega (2D). Henkes, van
%   Hecke & van Saarloos, EPL 90, 14003 (2010) define omega* as the point
%   where D(omega), normalized by the plateau height D(omega_ref) at a
%   reference frequency inside the plateau, first reaches 0.2:
%
%       D(omega*) / D(omega_ref) = scalThreshold  (default 0.2).
%
%   The DOS is scanned upward from the lowest bin; the crossing is linearly
%   interpolated between the two bins that bracket it. NaN is returned if
%   D/D(omega_ref) is already above the threshold in the first bin (omega*
%   below the resolution) or never reaches it below omega_ref.
%
%   Arguments:
%       vecBinCenters  bin centers of the DOS (ascending)
%       vecDOS         D at the bin centers (calcDOS)
%       scalOmegaRef   reference frequency inside the plateau; D(omega_ref)
%                      is interpolated from the bins
%       scalThreshold  fraction of the plateau height (default 0.2)
%
%   Example:
%       [vecDOS, vecW] = calcDOS(vecOmegaHat, scalNumModes, linspace(0, 3, 61));
%       scalOmegaStar = calcOmegaStar(vecW, vecDOS, 1);

    arguments
        vecBinCenters (:, 1) double
        vecDOS        (:, 1) double
        scalOmegaRef  (1, 1) double {mustBePositive}
        scalThreshold (1, 1) double {mustBePositive} = 0.2
    end

    scalOmegaStar = NaN;
    scalPlateau = interp1(vecBinCenters, vecDOS, scalOmegaRef, 'linear');
    if ~(scalPlateau > 0)
        return;
    end
    vecRatio = vecDOS/scalPlateau;
    vecBoolBelowRef = vecBinCenters < scalOmegaRef;
    idxFirst = find(vecBoolBelowRef & vecRatio >= scalThreshold, 1, 'first');
    if isempty(idxFirst) || idxFirst == 1
        return;
    end
    % Linear interpolation between the bins bracketing the crossing.
    scalW1 = vecBinCenters(idxFirst - 1);
    scalW2 = vecBinCenters(idxFirst);
    scalR1 = vecRatio(idxFirst - 1);
    scalR2 = vecRatio(idxFirst);
    scalOmegaStar = scalW1 + (scalThreshold - scalR1)*(scalW2 - scalW1)/(scalR2 - scalR1);
end
