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
    % Example:
    %   plotModes('3dUniformMass/', 1.0);
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

    fig = figure;
    t = tiledlayout(2,numGamma,'TileSpacing','tight','Padding','tight'); %#ok<NASGU>

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
            P_color_val = ...
                (log10(m.P) - log10(minP)) / ...
                (log10(maxP) - log10(minP));
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

