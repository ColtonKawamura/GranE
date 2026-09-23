function plotModes(resultsDir, avg_mass, options)
    % plotModes(resultsDir, avg_mass, options)
    %
    % Default behavior:
    %   - Scan resultsDir for:
    %       results_3D_N*_P*_Width*_Seed*_gamma_*.mat
    %       results_2D_iso_N*_P*_Seed*_gamma_*.mat
    %   - Plot ALL files found (2D + 3D).
    %
    % Optional filtering (closest matches, not exact):
    %   options.argGamma    = [0.01 0.001];      % target gammas
    %   options.argPressure = [0.1 0.02 0.03];   % target pressures
    %   options.argN        = [216 1000];        % target N values
    %   options.argSeed     = [1 3];             % target seeds
    %
    % Optional pressure scaling:
    %   options.scalPressure = 0.25;             % multiply y by P^0.25
    %                                            % (0 or empty -> no scaling)
    %
    % Optional coordination-number scaling row:
    %   options.flagCoordNumRow = true;         % add a 3rd row whose y is the
    %                                            % (beta/gamma) ratio normalized by
    %                                            % the packing's mean coordination
    %                                            % number (scalMeanCoordNum). Default
    %                                            % OFF; if the saved packing lacks
    %                                            % scalMeanCoordNum it is recomputed
    %                                            % from the saved positions, exactly
    %                                            % as pack.m does.
    %
    % Example:
    %   plotModes('3dUniformMass/', 1.0);
    %   plotModes("/Users/coltonkawamura/repos/GranE/data/eigenData/iso3d-uniform-mass/", 1)
    %
    %   opts = struct()
    %   opts.flagCoordNumRow = true;
    %   plotModes("~/repos/GranE/data/eigenData/iso3d-uniform-mass/", 1, opts)
    %
    %   opts.argGamma    = [0.01 0.001];
    %   opts.argPressure = [0.1 0.02];
    %   opts.scalPressure = 0.25;   % use P^(1/4) scaling
    %   plotModes('3dUniformMass/', 1.0, opts);

    if nargin < 3 || isempty(options)
        options = struct();
    end

    %--------------------------------------------------------------
    % Scan directory and collect metadata for all result files
    %--------------------------------------------------------------
    files3D = dir(fullfile(resultsDir,'results_3D_N*_gamma_*.mat'));
    files2D = dir(fullfile(resultsDir,'results_2D_iso_N*_gamma_*.mat'));

    meta = struct('file',{},'dim',{},'N',{},'P',{},'Width',{},'seed',{},'gamma',{});

    % --- 3D files ---
    for k = 1:numel(files3D)
        name = files3D(k).name;
        tok = regexp(name, ...
            '^results_3D_N(\d+)_P([0-9\.eE+-]+)_Width([0-9\.eE+-]+)_Seed(\d+)_gamma_([0-9\.eE+-]+)\.mat$', ...
            'tokens','once');
        if isempty(tok), continue; end
        N     = str2double(tok{1});
        P     = str2double(tok{2});
        Width = str2double(tok{3});
        seed  = str2double(tok{4});
        gamma = str2double(tok{5});
        meta(end+1) = struct( ...
            'file',  fullfile(files3D(k).folder,name), ...
            'dim',   '3D', ...
            'N',     N, ...
            'P',     P, ...
            'Width', Width, ...
            'seed',  seed, ...
            'gamma', gamma);
    end

    % --- 2D files ---
    for k = 1:numel(files2D)
        name = files2D(k).name;
        tok = regexp(name, ...
            '^results_2D_iso_N(\d+)_P([0-9\.eE+-]+)_Seed(\d+)_gamma_([0-9\.eE+-]+)\.mat$', ...
            'tokens','once');
        if isempty(tok), continue; end
        N     = str2double(tok{1});
        P     = str2double(tok{2});
        seed  = str2double(tok{3});
        gamma = str2double(tok{4});
        meta(end+1) = struct( ...
            'file',  fullfile(files2D(k).folder,name), ...
            'dim',   '2D', ...
            'N',     N, ...
            'P',     P, ...
            'Width', NaN, ...
            'seed',  seed, ...
            'gamma', gamma);
    end

    if isempty(meta)
        error('plotModes:NoFiles',...
              'No matching results files found in directory: %s',resultsDir);
    end

    %--------------------------------------------------------------
    % Build lists of all values actually present
    %--------------------------------------------------------------
    allGamma = unique([meta.gamma]);
    allP     = unique([meta.P]);
    allN     = unique([meta.N]);
    allSeed  = unique([meta.seed]);

    % Apply "closest match" selection if options are given
    argGamma    = getfield_or_empty(options,'argGamma');
    argPressure = getfield_or_empty(options,'argPressure');
    argN        = getfield_or_empty(options,'argN');
    argSeed     = getfield_or_empty(options,'argSeed');

    % Pressure scaling exponent (0 -> no scaling)
    scalPressure = getfield_or_empty(options,'scalPressure');
    if isempty(scalPressure)
        scalPressure = 0;
    end

    % Third row: additionally normalize y by the packing's mean coordination
    % number (scalMeanCoordNum). Default off.
    flagCoordNumRow = logical(getfield_or_empty(options,'flagCoordNumRow'));
    if isempty(flagCoordNumRow)
        flagCoordNumRow = false;
    end

    gamma_vals = selectClosest(allGamma, argGamma);    % columns in tiledlayout
    P_list     = selectClosest(allP,     argPressure);
    N_list     = selectClosest(allN,     argN);
    seed_list  = selectClosest(allSeed,  argSeed);

    %--------------------------------------------------------------
    % Figure + layout
    %--------------------------------------------------------------
    AxisFontSize = 16;
    SYMLIST = 's*>po.^vx+hd';   % enough distinct symbols for many N values
    curr_version = version();

    numGamma = numel(gamma_vals);

    numRows = 2 + flagCoordNumRow;
    fig = figure;
    t = tiledlayout(numRows,numGamma,'TileSpacing','tight','Padding','tight'); %#ok<NASGU>

    if curr_version(2) ~= '3'
        theme(fig,'light');
    end

    % Color scaling uses the P_list actually used
    minP = min(P_list);
    maxP = max(P_list);

    %--------------------------------------------------------------
    % Prepare dynamic y-label string based on scalPressure
    %--------------------------------------------------------------
    if scalPressure == 0
        yLabelStr = '$\hat{\beta}_i/\hat{\gamma}$';
    elseif abs(scalPressure - 0.25) < 1e-8
        % Special case: show 1/4 nicely
        yLabelStr = '$\left(\hat{\beta}_i/\hat{\gamma}\right)\hat{P}^{1/4}$';
    else
        % Generic exponent, e.g. P^0.3
        yLabelStr = ['$\left(\hat{\beta}_i/\hat{\gamma}\right)\hat{P}^{' num2str(scalPressure) '}$'];
    end

    % Y-label for the coordination-number-normalized (third) row: the same
    % quantity as yLabelStr, divided by the mean coordination number z-bar
    % (kept inside one math span).
    zLabelStr = [yLabelStr(1:end-1) ' / \bar{z}$'];

    %--------------------------------------------------------------
    % Loop over gamma columns
    %--------------------------------------------------------------
    for gidx = 1:numGamma

        gamma = gamma_vals(gidx);

        %----------------------------------------------------------
        % TOP ROW: raw scaling
        %----------------------------------------------------------
        ax1 = nexttile(gidx);

        hold(ax1,'on')
        set(ax1,...
            'XScale','log',...
            'YScale','log',...
            'Box','on')

        grid(ax1,'on')

        %----------------------------------------------------------
        % BOTTOM ROW: P-scaled frequency
        %----------------------------------------------------------
        ax2 = nexttile(numGamma + gidx);

        hold(ax2,'on')
        set(ax2,...
            'XScale','log',...
            'YScale','log',...
            'Box','on')

        grid(ax2,'on')

        %----------------------------------------------------------
        % THIRD ROW: y additionally normalized by the packing's mean
        % coordination number (only when the option is enabled)
        %----------------------------------------------------------
        ax3 = [];
        if flagCoordNumRow
            ax3 = nexttile(2*numGamma + gidx);
            hold(ax3,'on')
            set(ax3,...
                'XScale','log',...
                'YScale','log',...
                'Box','on')
            grid(ax3,'on')
        end

        %----------------------------------------------------------
        % Loop over all files and plot those matching this gamma
        % and the selected P/N/seed lists
        %----------------------------------------------------------
        for k = 1:numel(meta)
            m = meta(k);

            % Match gamma exactly (gamma_vals is subset of allGamma)
            if m.gamma ~= gamma
                continue;
            end
            % Filter by closest-selected P, N, seed
            if ~ismember(m.P,    P_list),    continue; end
            if ~ismember(m.N,    N_list),    continue; end
            if ~ismember(m.seed, seed_list), continue; end

            % Load and reactively check 2D vs 3D via packing struct
            results = load(m.file);

            is3D = false;
            if isfield(results,'packing') && ...
                    isfield(results.packing,'scalBoxDepthZ') && ...
                    (results.packing.scalBoxDepthZ ~= 0)
                is3D = true;
            end

            if is3D
                disp([m.file ' [3D]']);
            else
                disp([m.file ' [2D]']);
            end

            eigvals = results.eigenvalues;

            % Color from pressure
            pDenom = log10(maxP) - log10(minP);
            if pDenom == 0 || ~isfinite(pDenom)
                % Single (or degenerate) pressure: no color ramp to span.
                P_color_val = 0.5;
            else
                P_color_val = (log10(m.P) - log10(minP)) / pDenom;
            end
            P_color = [P_color_val, 0, 1-P_color_val];

            % Symbol from N
            idxN = find(N_list == m.N, 1);
            if isempty(idxN)
                idxN = 1;
            end
            idxN = min(idxN, length(SYMLIST));
            SYM = SYMLIST(idxN);

            w = abs(imag(eigvals));
            a = -real(eigvals);

            %------------------------------------------------------
            % Apply pressure scaling to attenuation (if requested)
            %------------------------------------------------------
            yval = a*sqrt(avg_mass)/(gamma/sqrt(avg_mass)); % original \hat{\beta}_i/\hat{\gamma}
            if scalPressure ~= 0
                yval = yval .* (m.P.^scalPressure);
            end

            %------------------------------------------------------
            % Top row
            %------------------------------------------------------
            loglog(ax1,...
                w*sqrt(avg_mass),...
                yval,...
                SYM,...
                'Color',P_color)

            %------------------------------------------------------
            % Bottom row
            %------------------------------------------------------
            loglog(ax2,...
                w*sqrt(avg_mass)/sqrt(m.P),...
                yval,...
                SYM,...
                'Color',P_color)

            %------------------------------------------------------
            % Third row: same reduced-frequency x, y additionally
            % normalized by the packing's mean coordination number.
            % If scalMeanCoordNum is not saved with the packing, it is
            % recomputed from the saved positions, exactly as pack.m does.
            %------------------------------------------------------
            if flagCoordNumRow
                % Zn = getfield_or_empty(results.packing,'scalMeanCoordNum');
                Zn = computeMeanCoordNum(results.packing);
                if isempty(Zn) || ~isfinite(Zn) || Zn <= 0
                    Zn = computeMeanCoordNum(results.packing);
                end
                if isfinite(Zn) && Zn > 0
                    loglog(ax3,...
                        w*sqrt(avg_mass)/sqrt(m.P),...
                        yval / Zn,...
                        SYM,...
                        'Color',P_color)
                end
            end

        end

        %==========================================================
        % TOP ROW FORMATTING
        %==========================================================
        xlim(ax1,[3e-3 10])
        ylim(ax1,[3e-3 10])

        plot(ax1,[3e-2 1],[3e-2 1].^2,'k-')

        text(ax1,...
            0.2,...
            0.02,...
            '2',...
            'FontSize',16,...
            'Interpreter','latex')

        title(ax1,...
            ['$\hat{\gamma}=' num2str(gamma) '$'],...
            'Interpreter','latex',...
            'FontSize',AxisFontSize)

        set(ax1,'XTick',10.^(-3:1))

        xlabel(ax1,...
            '$\hat{\omega}_i$',...
            'FontSize',AxisFontSize,...
            'Interpreter','latex')

        % Only first column gets y-axis labels/ticks
        if gidx == 1

            ylabel(ax1,...
                yLabelStr,...
                'FontSize',AxisFontSize,...
                'Interpreter','latex')

            % Reference 1/3 line stays in original coordinates
            plot(ax1,[3e-2 1],3*[3e-2 1].^(1/3),'k-')
            text(ax1,...
                0.15,...
                4,...
                '1/3',...
                'FontSize',16,...
                'Interpreter','latex')

        else

            ax1.YTickLabel = [];

        end

        %==========================================================
        % BOTTOM ROW FORMATTING
        %==========================================================
        xlim(ax2,[3e-2 100])
        ylim(ax2,[3e-3 10])

        plot(ax2,[3e-2 1],[3e-2 1].^2,'k-')

        text(ax2,...
            0.1,...
            0.05,...
            '2',...
            'FontSize',16,...
            'Interpreter','latex')

        set(ax2,'XTick',10.^(-3:2))

        xlabel(ax2,...
            '$\hat{\omega}_i/\hat{\omega}_c$',...
            'FontSize',AxisFontSize,...
            'Interpreter','latex')

        % Only first column gets y-axis labels/ticks
        if gidx == 1

            ylabel(ax2,...
                yLabelStr,...
                'FontSize',AxisFontSize,...
                'Interpreter','latex')

        else

            ax2.YTickLabel = [];

        end

        %==========================================================
        % THIRD ROW FORMATTING (coordination-number normalized)
        %==========================================================
        if flagCoordNumRow
            xlim(ax3,[3e-2 100])
            ylim(ax3,[3e-3 10])

            plot(ax3,[3e-2 1],[3e-2 1].^2,'k-')

            set(ax3,'XTick',10.^(-3:2))

            xlabel(ax3,...
                '$\hat{\omega}_i/\hat{\omega}_c$',...
                'FontSize',AxisFontSize,...
                'Interpreter','latex')

            if gidx == 1
                ylabel(ax3,...
                    zLabelStr,...
                    'FontSize',AxisFontSize,...
                    'Interpreter','latex')
            else
                ax3.YTickLabel = [];
            end
        end

    end

end

%==============================================================
% Helper functions (nested)
%==============================================================
function valsOut = selectClosest(allVals, targetVals)
    % If no targets given, use allVals.
    if isempty(targetVals)
        valsOut = allVals;
        return;
    end
    valsOut = zeros(size(targetVals));
    for i = 1:numel(targetVals)
        [~,idx] = min(abs(allVals - targetVals(i)));
        valsOut(i) = allVals(idx);
    end
    valsOut = unique(valsOut);
end

function val = getfield_or_empty(s, fieldname)
    if isstruct(s) && isfield(s, fieldname)
        val = s.(fieldname);
    else
        val = [];
    end
end

function Zn = computeMeanCoordNum(packing)
    % computeMeanCoordNum(packing) -- mean particle-particle coordination
    % number (Zn) of a saved packing, recomputed from its positions.
    %
    % Mirrors pack.m: vecCoordNum is the per-particle count of OTHER
    % particles overlapping under the periodic minimum-image distance
    % (r_ij < r_i + r_j), and scalMeanCoordNum = mean(vecCoordNum). Wall
    % contacts are NOT counted (pack.m's vecCoordNum excludes them).
    %
    % Results files store the packing in one of two schemas, both handled:
    %   3D: vecPosX/vecPosY/vecPosZ (N x 1), vecDiameter (N x 1),
    %       scalBoxWidthX / scalBoxHeightY / scalBoxDepthZ
    %   2D: x / y / Dn (1 x N), Lx / Ly
    %
    % Returns NaN when the required fields are absent so the caller can
    % fall back (and leave the row blank for that point).
    if isfield(packing,'vecPosX') && isfield(packing,'vecPosY') && ...
            isfield(packing,'vecPosZ') && isfield(packing,'vecDiameter')
        % 3D schema
        vecX  = packing.vecPosX;
        vecY  = packing.vecPosY;
        vecZ  = packing.vecPosZ;
        vecD  = packing.vecDiameter;
        Lx    = packing.scalBoxWidthX;
        Ly    = packing.scalBoxHeightY;
        Lz    = packing.scalBoxDepthZ;
        bool3D = true;
    elseif isfield(packing,'x') && isfield(packing,'y') && isfield(packing,'Dn')
        % 2D schema
        vecX  = packing.x(:);
        vecY  = packing.y(:);
        vecD  = packing.Dn(:);
        Lx    = packing.Lx;
        Ly    = packing.Ly;
        bool3D = false;
    else
        Zn = NaN;
        return;
    end

    N = numel(vecX);
    if N < 2 || any(~isfinite([Lx; Ly])) || (bool3D && ~isfinite(Lz))
        Zn = NaN;
        return;
    end

    vecCoordNum = zeros(N, 1);
    for i = 1:N
        dx = vecX - vecX(i);
        dx = dx - Lx * round(dx / Lx);
        dy = vecY - vecY(i);
        dy = dy - Ly * round(dy / Ly);
        if bool3D
            dz = vecZ - vecZ(i);
            dz = dz - Lz * round(dz / Lz);
            r  = sqrt(dx.^2 + dy.^2 + dz.^2);
        else
            r  = sqrt(dx.^2 + dy.^2);
        end
        rSum = vecD + vecD(i);              % contact distance r_i + r_j
        vecCoordNum(i) = sum(r < rSum & (1:N)' ~= i);
    end
    Zn = mean(vecCoordNum);
end

