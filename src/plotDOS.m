function plotDOS(resultsDir, options)
%PLODOS Plot the vibrational density of states (DOS) of eigenmode data.
%
%   plotDOS(resultsDir)
%   plotDOS(resultsDir, options)
%
%   Scans resultsDir for eigenmode results files of the form
%       results_3D_N*_P*_Width*_Seed*_gamma_*.mat
%       results_2D_iso_N*_P*_Seed*_gamma_*.mat
%   (same naming convention as plotModes.m), computes the normalized
%   vibrational density of states for each file, and overlays the curves on
%   one log-log panel per damping value (gamma), grouped into columns.
%
%   The density of states g(omega) is the distribution of the mode
%   frequencies extracted from the complex eigenvalues lambda_i saved with
%   each result:
%       omega_i = |Im(lambda_i)|            (reduced angular frequency)
%
%   Default mode (flagScaleByCutoff = true): every spectrum is rescaled by
%   its own maximum (Debye-cutoff) frequency, x = omega / omega_max, so that
%   all files collapse onto the common axis x in (0, 1] regardless of their
%   absolute frequency scale, and g(x) is normalized to unit area.  This
%   makes the DOS of packings with very different frequency scales directly
%   comparable.
%
%   Absolute mode (flagScaleByCutoff = false): the DOS is plotted in the raw
%   reduced frequency omega with a shared axis across all files.
%
%   Degenerate spectra (runs whose omega_max is far below the largest in the
%   data, i.e. with no real modes) are skipped automatically.
%
%   For each panel (fixed gamma) every matching file is overlaid, colored by
%   its reduced pressure, together with a mean +/- 1 standard-deviation band
%   (translucent) and a bold mean curve.
%
%   options : name-value struct with fields:
%            argGamma        : target gammas (closest-match selection)
%            argPressure     : target pressures (closest-match selection)
%            argN            : target particle numbers (closest-match)
%            argSeed         : target seeds (closest-match selection)
%            numBins        : number of histogram bins (default 256)
%            flagScaleByCutoff : rescale each spectrum by its own omega_max
%                               (default true); when false plot the absolute
%                               reduced frequency instead
%            lightMode      : use a light figure theme (default false)
%            showBand       : display the mean +/- 1 std band and mean line
%                               (default true)
%
%   Example:
%       plotDOS('data/eigenData/3d/hooke/iso3d-uniform-mass/')
%       opts.argPressure = [0.01 0.1]; opts.numBins = 512;
%       plotDOS('data/eigenData/3d/hooke/iso3d-uniform-mass/', opts)
%
%         o = struct();
%         o.argN = 1000;
%         o.argGamma = 0.001;
%         o.flagScaleByCutoff = false;
%         o.numBins = 40;
%         plotDOS('data/eigenData/3d/hooke/iso3d-uniform-mass/', o);

    if nargin < 2 || isempty(options)
        options = struct();
    end

    try

%% Parameters
        scalNumBins = getfield_or_empty(options, 'numBins');
        if isempty(scalNumBins) || scalNumBins < 2
            scalNumBins = 256;
        end

        flagScaleByCutoff = logical(getfield_or_empty(options, 'flagScaleByCutoff'));
        if isempty(flagScaleByCutoff)
            flagScaleByCutoff = true;
        end

        flagLightMode = logical(getfield_or_empty(options, 'lightMode'));
        if isempty(flagLightMode)
            flagLightMode = false;
        end

        flagShowBand = logical(getfield_or_empty(options, 'showBand'));
        if isempty(flagShowBand)
            flagShowBand = true;
        end

        axisFontSize = 16;

%% Scan directory and collect metadata (same convention as plotModes.m)
        files3D = dir(fullfile(resultsDir, 'results_3D_N*_gamma_*.mat'));
        files2D = dir(fullfile(resultsDir, 'results_2D_iso_N*_gamma_*.mat'));

        meta = struct('file', {}, 'dim', {}, 'N', {}, 'P', {}, 'Width', {}, ...
                      'seed', {}, 'gamma', {});

        for k = 1:numel(files3D)
            name = files3D(k).name;
            tok = regexp(name, ...
                '^results_3D_N(\d+)_P([0-9\.eE+-]+)_Width([0-9\.eE+-]+)_Seed(\d+)_gamma_([0-9\.eE+-]+)\.mat$', ...
                'tokens', 'once');
            if isempty(tok), continue; end
            meta(end + 1) = struct('file',  fullfile(files3D(k).folder, name), ...
                                   'dim',   '3D', ...
                                   'N',     str2double(tok{1}), ...
                                   'P',     str2double(tok{2}), ...
                                   'Width', str2double(tok{3}), ...
                                   'seed',  str2double(tok{4}), ...
                                   'gamma', str2double(tok{5}));
        end

        for k = 1:numel(files2D)
            name = files2D(k).name;
            tok = regexp(name, ...
                '^results_2D_iso_N(\d+)_P([0-9\.eE+-]+)_Seed(\d+)_gamma_([0-9\.eE+-]+)\.mat$', ...
                'tokens', 'once');
            if isempty(tok), continue; end
            meta(end + 1) = struct('file',  fullfile(files2D(k).folder, name), ...
                                   'dim',   '2D', ...
                                   'N',     str2double(tok{1}), ...
                                   'P',     str2double(tok{2}), ...
                                   'Width', NaN, ...
                                   'seed',  str2double(tok{4}), ...
                                   'gamma', str2double(tok{5}));
        end

        if isempty(meta)
            error('plotDOS:NoFiles', ...
                  'No matching results files found in directory: %s', resultsDir);
        end

%% Filtering (closest-match selection, reused from plotModes.m)
        argGamma    = getfield_or_empty(options, 'argGamma');
        argPressure = getfield_or_empty(options, 'argPressure');
        argN        = getfield_or_empty(options, 'argN');
        argSeed     = getfield_or_empty(options, 'argSeed');

        allGamma = unique([meta.gamma]);
        allP     = unique([meta.P]);
        allN     = unique([meta.N]);
        allSeed  = unique([meta.seed]);

        gammaValsSel = selectClosest(allGamma, argGamma);
        PList        = selectClosest(allP,   argPressure);
        NList        = selectClosest(allN,   argN);
        seedList     = selectClosest(allSeed, argSeed);

        % Keep only the files that match ALL four (closest-match) selections,
        % exactly as plotModes.m filters per file.
        keep = false(numel(meta), 1);
        for k = 1:numel(meta)
            if ismember(meta(k).gamma, gammaValsSel) && ...
                   ismember(meta(k).P,   PList)   && ...
                   ismember(meta(k).N,   NList)   && ...
                   ismember(meta(k).seed, seedList)
                keep(k) = true;
            end
        end
        meta = meta(keep);
        if isempty(meta)
            error('plotDOS:NoFiles', ...
                'No result files match the requested selection.');
        end

%% Per-file reduced frequencies (drop near-zero numerical modes)
        metaValid     = meta([]);
        cellOmegaMax  = {};
        cellFreq      = {};
        numFiles      = 0;

        for k = 1:numel(meta)
            results = load(meta(k).file);
            if ~isfield(results, 'eigenvalues') || isempty(results.eigenvalues)
                warning('plotDOS:NoEigenvalues', ...
                        'No eigenvalues in %s; skipping.', meta(k).file);
                continue;
            end

            vecOmega = abs(imag(results.eigenvalues(:)));   % reduced frequency >= 0
            vecOmega = vecOmega(vecOmega > 0);
            if isempty(vecOmega)
                continue;
            end

            scalOmegaMax = max(vecOmega);
            % Drop the near-zero numerical modes (static / solver-noise modes);
            % the physical spectrum starts well above 1e-3 * omega_max.
            vecOmega = vecOmega(vecOmega > 1e-3 * scalOmegaMax);
            if isempty(vecOmega)
                continue;
            end

            numFiles  = numFiles + 1;
            metaValid(end + 1) = meta(k);
            cellOmegaMax{numFiles} = scalOmegaMax;
            cellFreq{numFiles} = vecOmega;
        end

        if numFiles == 0
            error('plotDOS:NoData', 'No result files contained usable eigenvalues.');
        end

%% Drop degenerate spectra (no real modes) — omega_max far below the largest
        scalOmegaGlobalMax = max([cellOmegaMax{:}]);
        scalOmegaFloor     = 1e-3 * scalOmegaGlobalMax;
        keepIdx            = find(cellfun(@(v) v > scalOmegaFloor, cellOmegaMax));
        if isempty(keepIdx)
            error('plotDOS:NoData', 'No result files contained a real spectrum.');
        end
        metaValid    = metaValid(keepIdx);
        cellFreq     = cellFreq(keepIdx);
        numFiles     = numel(keepIdx);

%% Build the frequency axis and the normalized DOS for each file
        if flagScaleByCutoff
            % Rescale each spectrum by its own Debye cutoff: x = omega/omega_max.
            xMin = 1e-3;
            vecEdges = logspace(log10(xMin), 0, scalNumBins + 1);
            vecGridMin = xMin;
            vecGridMax = 1;
            strXLabel = '$\hat{\omega}/\hat{\omega}_{\max}$';
            for k = 1:numFiles
                cellFreq{k} = cellFreq{k} / cellOmegaMax{k};
            end
        else
            % Absolute reduced frequency over the shared range present in data.
            vecAllFreq = cat(1, cellFreq{:});
            vecGridMin = min(vecAllFreq);
            vecGridMax = max(vecAllFreq);
            vecEdges = logspace(log10(vecGridMin), log10(vecGridMax), scalNumBins + 1);
            strXLabel = '$\hat{\omega}$';
        end

        cellCenters = {};
        cellG       = {};
        for k = 1:numFiles
            [counts, ~] = histcounts(cellFreq{k}, vecEdges, 'normalization', 'pdf');
            % geometric bin centers (true centers in log space)
            cellCenters{k} = sqrt(vecEdges(1:end - 1) .* vecEdges(2:end));
            cellG{k} = counts(:);
        end

%% Common frequency grid for the mean +/- std band
        vecOmegaGrid = logspace(log10(vecGridMin), log10(vecGridMax), 400);
        cellGI = cell(numFiles, 1);
        for k = 1:numFiles
            gI = interp1(cellCenters{k}, cellG{k}, vecOmegaGrid, 'linear', NaN);
            valid = isfinite(gI);
            if ~all(valid)
                firstF = find(valid, 1);
                lastF  = find(valid, 1, 'last');
                % No modes below the lowest frequency -> g = 0 there.
                if firstF > 1
                    gI(1:firstF - 1) = 0;
                end
                % Smooth constant tail above the highest populated bin.
                if lastF < numel(gI)
                    gI(lastF + 1:end) = gI(lastF);
                end
            end
            cellGI{k} = gI;
        end

%% Build the figure
        gammaVals = unique([metaValid.gamma]);
        numGamma  = numel(gammaVals);
        fig = figure('Name', 'plotDOS');
        currVersion = version();
        if currVersion(2) ~= '3'
            try
                if flagLightMode
                    theme(fig, 'light');
                else
                    theme(fig, 'dark');
                end
            catch
                % theme() unavailable on this MATLAB version; fall back.
            end
        end
        t = tiledlayout(1, numGamma, 'TileSpacing', 'tight', 'Padding', 'tight'); %#ok<NASGU>

        minP = min(PList);
        maxP = max(PList);
        pDenom = log10(maxP) - log10(minP);
        if ~isfinite(pDenom) || pDenom == 0
            pDenom = 1;
        end

        decLo = floor(log10(vecGridMin));
        decHi = ceil(log10(vecGridMax));
        tickDecs = 10.^(decLo:decHi);
        xticksSet = tickDecs(tickDecs >= vecGridMin & tickDecs <= vecGridMax);
        if isempty(xticksSet)
            xticksSet = [vecGridMin, vecGridMax];
        end

        for gidx = 1:numGamma
            gamma = gammaVals(gidx);

            colIdx = find([metaValid.gamma] == gamma);

            ax = nexttile(gidx);
            hold(ax, 'on');
            set(ax, 'XScale', 'log', 'YScale', 'log', 'Box', 'on');
            grid(ax, 'on');

            % Overlay every matching file, colored by reduced pressure
            for ci = 1:numel(colIdx)
                kk = colIdx(ci);
                pVal = (log10(metaValid(kk).P) - log10(minP)) / pDenom;
                pCol = [pVal, 0, 1 - pVal];
                loglog(ax, cellCenters{kk}, cellG{kk}, ...
                       'Color', pCol, 'LineWidth', 1, ...
                       'Marker', 'o', 'MarkerSize', 3, 'MarkerFaceColor', pCol);
            end

            % Mean +/- 1 std band and mean curve
            if flagShowBand && numel(colIdx) > 1
                gMat = zeros(numel(vecOmegaGrid), numel(colIdx));
                for ci = 1:numel(colIdx)
                    gMat(:, ci) = cellGI{colIdx(ci)};
                end
                gMean = mean(gMat, 2);
                gStd  = std(gMat, 0, 2);
                fill([vecOmegaGrid(:); fliplr(vecOmegaGrid(:))], ...
                     [(gMean + gStd); fliplr(gMean - gStd)], ...
                     [0.55 0.55 0.55], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
                loglog(ax, vecOmegaGrid(:), gMean, 'k-', 'LineWidth', 2.5);
            end

            % Tighten the y-axis so the curves fill the panel (log scale
            % cannot display zero, so floor the axis at 1e-3).
            gVals = [];
            for ci = 1:numel(colIdx)
                gVals = [gVals; cellG{colIdx(ci)}];
            end
            if flagShowBand && numel(colIdx) > 1
                gVals = [gVals; gMean];
            end
            gVals = gVals(isfinite(gVals));
            if ~isempty(gVals)
                set(ax, 'Ylim', [1e-3, 10^ceil(log10(max(gVals) * 2.5))]);
            end

            hold(ax, 'off');

            title(ax, ['$\gamma=' num2str(gamma) '$'], ...
                  'Interpreter', 'latex', 'FontSize', axisFontSize);

            set(ax, 'XTick', xticksSet);
            set(ax, 'Xlim', [xticksSet(1), xticksSet(end)]);

            xlabel(ax, strXLabel, 'Interpreter', 'latex', 'FontSize', axisFontSize);
            if gidx == 1
                ylabel(ax, '$g(\hat{\omega})$', 'Interpreter', 'latex', 'FontSize', axisFontSize);
            else
                ax.YLabel.String = '';
                ax.YTickLabel = [];
            end
        end

%% Pressure colorbar
        if numel(unique([metaValid.P])) > 1
            cb = colorbar;
            cbTicks = 10.^(floor(log10(minP)):ceil(log10(maxP)));
            cbTicks = cbTicks(cbTicks >= minP & cbTicks <= maxP);
            if isempty(cbTicks)
                cbTicks = [minP, maxP];
            end
            set(cb, 'YTick', cbTicks);
            ylabel(cb, '$\hat{P}$', 'Interpreter', 'latex', 'FontSize', 14);
        end

        fprintf('[plotDOS] %d file(s) plotted across %d gamma column(s) [scaleByCutoff=%d].\n', ...
                numFiles, numGamma, flagScaleByCutoff);

    catch sPlotDOSME
        fprintf(2, '[plotDOS] ERROR: %s\n', sPlotDOSME.message);
        for scalK = 1:length(sPlotDOSME.stack)
            fprintf(2, '[plotDOS]   in %s (line %d)\n', ...
                    sPlotDOSME.stack(scalK).name, sPlotDOSME.stack(scalK).line);
        end
        rethrow(sPlotDOSME);
    end
end

%==============================================================
% Helper functions (nested)
%==============================================================
function valsOut = selectClosest(allVals, targetVals)
    % selectClosest(allVals, targetVals) -> closest member of allVals for
    % each target (reused from plotModes.m). If targetVals is empty, returns
    % allVals unchanged.
    if isempty(targetVals)
        valsOut = allVals;
        return;
    end
    valsOut = zeros(size(targetVals));
    for i = 1:numel(targetVals)
        [~, idx] = min(abs(allVals - targetVals(i)));
        valsOut(i) = allVals(idx);
    end
    valsOut = unique(valsOut);
end

function val = getfield_or_empty(s, fieldname)
    % getfield_or_empty(s, fieldname) returns s.fieldname, or [] if s is not
    % a struct or does not have that field (reused from plotModes.m).
    if isstruct(s) && isfield(s, fieldname)
        val = s.(fieldname);
    else
        val = [];
    end
end
