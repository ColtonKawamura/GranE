function structStats = meanPhaseDev(vecXValues, vecYValues, scalBinWidth, structOptions)
%MEANPHASEDEV Compute the mean of the per-bin standard deviations of y vs x.
%Ported from GranMA src/matlab_functions/meanPhaseDev.m (hungarian naming).
%
%   Args:
%       vecXValues  (1,:) double
%       vecYValues  (1,:) double
%       scalBinWidth (1,1) double
%       structOptions.cV (1,1) logical = false % Optional: flag for coefficient of variation

    arguments
        vecXValues (1,:) double
        vecYValues (1,:) double
        scalBinWidth (1,1) double
        structOptions.cV (1,1) logical = false % Optional: flag for coefficient of variation
    end

    % Ensure xValues and yValues are sorted by xValues.
    [vecXSorted, vecSortedIdx] = sort(vecXValues);
    vecYSorted = vecYValues(vecSortedIdx);

    % Make the bins.
    scalMinX = min(vecXSorted);
    scalMaxX = max(vecXSorted);
    vecBinEdges = scalMinX:scalBinWidth:scalMaxX;

    % Ensure the last edge includes scalMaxX if it falls exactly on a bin edge.
    if vecBinEdges(end) < scalMaxX
        vecBinEdges = [vecBinEdges, scalMaxX];
    end
    vecBinStddevs = [];
    vecBinMeans = [];



    % Iterate over each bin to calculate the standard deviation of yValues
    % within the bin.
    for scalI = 1:length(vecBinEdges) - 1
        scalBinStart = vecBinEdges(scalI);
        scalBinEnd = vecBinEdges(scalI + 1);

        % Find indices of xValues within the current bin range. Include the
        % start edge, exclude the end edge, except for the last bin.
        if scalI < length(vecBinEdges) - 1
            vecBinIdx = find(vecXSorted >= scalBinStart & vecXSorted < scalBinEnd);
        else
            % Include the end edge for the last bin.
            vecBinIdx = find(vecXSorted >= scalBinStart & vecXSorted <= scalBinEnd);
        end

        if ~isempty(vecBinIdx)
            vecBinStddevs(end+1) = std(vecYSorted(vecBinIdx));
            vecBinMeans(end+1) = mean(vecYSorted(vecBinIdx));
        end
    end

    structStats = mean(vecBinStddevs(~isnan(vecBinStddevs)));
    if structOptions.cV
        vecMeanMeans = mean(vecBinMeans(~isnan(vecBinMeans)));
        structStats = structStats / vecMeanMeans;
    end

    % If all bins were empty or resulted in NaN, set as NaN.
    if isempty(structStats)
        structStats = NaN;
    end
end
