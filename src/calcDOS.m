function [vecDOS, vecBinCenters] = calcDOS(vecOmega, scalNumModes, vecBinEdges)
% calcDOS  Vibrational density of states D(omega) as a normalized histogram.
%
%   [vecDOS, vecBinCenters] = calcDOS(vecOmega, scalNumModes, vecBinEdges)
%
%   D(omega) d(omega) is the fraction of all modes with frequency in
%   [omega, omega + d(omega)]. Following Henkes, van Hecke & van Saarloos,
%   EPL 90, 14003 (2010), the spectrum is the histogram of the eigenvalues
%   omega_k^2 of the dynamical matrix as a function of omega_k:
%
%       D(omega_b) = (number of modes in bin b) / (scalNumModes * width_b).
%
%   The histogram is normalized by the TOTAL number of modes, including the
%   zero and unstable ones, so the integral of D over the bins is the
%   fraction of finite-frequency modes. Spectra of several packings are
%   averaged by passing their pooled frequencies and the summed mode count.
%
%   Arguments:
%       vecOmega       frequencies (e.g. sSpectrum.vecOmega from
%                      calcEigenFrequencies, possibly pooled over packings
%                      and rescaled). NaN entries (unstable modes) and
%                      frequencies outside the bins are not counted.
%       scalNumModes   total number of modes the frequencies come from
%       vecBinEdges    bin edges (linear or logarithmic)
%
%   Returns:
%       vecDOS         [numBins x 1] density of states
%       vecBinCenters  [numBins x 1] bin centers (arithmetic mean of the edges)
%
%   Example:
%       sSpectrum = calcEigenFrequencies(matK, matM);
%       [vecDOS, vecW] = calcDOS(sSpectrum.vecOmega, sSpectrum.scalNumModes, linspace(0, 30, 61));
%       plot(vecW, vecDOS)

    arguments
        vecOmega      (:, 1) double
        scalNumModes  (1, 1) double {mustBePositive}
        vecBinEdges   (:, 1) double
    end

    vecOmega = vecOmega(~isnan(vecOmega));
    vecCounts = histcounts(vecOmega, vecBinEdges);
    vecBinWidths = diff(vecBinEdges);
    vecDOS = vecCounts(:)./(scalNumModes*vecBinWidths(:));
    vecBinCenters = (vecBinEdges(1:end-1) + vecBinEdges(2:end))/2;
end
