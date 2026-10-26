function hFig = plotUncertContacts(sUncert, options)
% plotUncertContacts  How much the sub-noise contacts matter.
%
%   plotUncertContacts(sUncert)
%   plotUncertContacts(sUncert, 'boolLightMode', false, 'strSavePath', "uncertContacts.png")
%   hFig = plotUncertContacts(...)
%
%   A contact whose normal force is below scalForceFloor*<f> (pack.m's
%   force-balance tolerance, 1e-4) is inside the equilibration noise: it is
%   not resolved whether it exists, yet cleanRats and matSpringDampMass count
%   it like any other contact. This figure shows how impactful those
%   contacts are, one row per contact law:
%     (a) how deep the weak tail of f/<f> goes on the backbone,
%     (b) how many weak contacts there are (backbone, and all contacts),
%     (c) their share of the backbone's normal stiffness sum k_n,
%     (d) how much cutting them would move z, as a share of z - z_iso
%         (the distance from isostaticity the vibrational spectrum follows).
%
%   sUncert is the 4th output of cleanRats (one struct per packing; pass a
%   struct array to compare packings). cleanRats calls this function itself
%   when its sUncertOptions.boolPlot is true.
%
%   Colour: pressure, red (high) -> blue (low) on a log scale, as in
%   plotAttenuation. Squares; larger squares = larger friction coefficient.
%
%   Options:
%       boolLightMode  light (true, default) or dark (false) figure
%       strSavePath    export a PNG there ("" = no export)

    arguments
        sUncert (1, :) struct
        options.boolLightMode (1, 1) logical = true
        options.strSavePath (1, 1) string = ""
    end

    % ---- theme -------------------------------------------------------------
    if options.boolLightMode
        vecBg = [1 1 1];  vecFg = [0 0 0];  vecGrid = [0.15 0.15 0.15];  vecShade = [0.85 0.85 0.85];
    else
        vecBg = [0.118 0.118 0.118];  vecFg = [0.88 0.88 0.88];  vecGrid = [0.7 0.7 0.7];  vecShade = [0.25 0.25 0.25];
    end

    % ---- colour by pressure (normVarColor: [n, 0, 1-n], n = normalized log P) ----
    vecP = [sUncert.scalPressure];
    vecPList = unique(vecP(~isnan(vecP)));
    vecLogP = log(vecPList);
    if numel(vecPList) > 1
        vecNormP = (vecLogP - min(vecLogP)) / (max(vecLogP) - min(vecLogP));
    else
        vecNormP = zeros(size(vecPList));
    end
    matColourP = [vecNormP(:), zeros(numel(vecNormP), 1), 1 - vecNormP(:)];
    if options.boolLightMode
        matColourP = min(matColourP + 0.2*(1 - matColourP), 1);    % as plotAttenuation's lightMode
    else
        matColourP = min(matColourP + 0.35*(1 - matColourP), 1);   % lifted for a dark background
    end

    % ---- square size by friction (normalized log mu) ----
    vecMu = [sUncert.scalMu];
    vecMuList = unique(vecMu);
    vecLogMu = log(max(vecMuList, eps));
    if numel(vecMuList) > 1
        vecNormMu = (vecLogMu - min(vecLogMu)) / (max(vecLogMu) - min(vecLogMu));
    else
        vecNormMu = 0.5*ones(size(vecMuList));
    end
    vecSizeMu = 5 + 9*vecNormMu;

    cellLaw = unique({sUncert.strContactLaw}, 'stable');
    scalNumRows = numel(cellLaw);
    scalFloor = sUncert(1).scalForceFloor;

    hFig = figure('Color', vecBg, 'Position', [60 60 1750 430*scalNumRows]);
    try
        if options.boolLightMode
            theme(hFig, 'light');
        else
            theme(hFig, 'dark');
        end
    catch
        % theme() exists only in R2025a+
    end

    for idxRow = 1:scalNumRows
        vecIdxLaw = find(strcmp({sUncert.strContactLaw}, cellLaw{idxRow}));
        strLaw = 'Hooke';
        if strcmp(cellLaw{idxRow}, 'hertz')
            strLaw = 'Hertz';
        end
        cellAx = cell(1, 4);
        for idxCol = 1:4
            cellAx{idxCol} = subplot(scalNumRows, 4, 4*(idxRow - 1) + idxCol);
            hold(cellAx{idxCol}, 'on');
        end
        [axTail, axCount, axStiff, axDz] = cellAx{:};

        % (a) cumulative distribution of f/<f> on the backbone
        patch(axTail, [1e-9 scalFloor scalFloor 1e-9], [1e-4 1e-4 1 1], vecShade, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        for idxK = vecIdxLaw
            vecF = sort(sUncert(idxK).vecForceRel(sUncert(idxK).vecBoolBackbone));
            vecCdf = (1:numel(vecF))'/numel(vecF);
            plot(axTail, vecF, vecCdf, '-', 'Color', colourOf(sUncert(idxK).scalPressure), 'LineWidth', 0.9, 'HandleVisibility', 'off');
        end
        xline(axTail, scalFloor, '--', 'Color', vecFg, 'LineWidth', 1.3, 'HandleVisibility', 'off');
        text(axTail, 2.5e-9, 0.45, {'$\mathrm{below\ the}$', '$\mathrm{noise\ floor}$'}, 'Interpreter', 'latex', 'FontSize', 12, 'Color', vecFg);
        set(axTail, 'XScale', 'log', 'YScale', 'log', 'XLim', [1e-9 10], 'YLim', [1e-4 1], 'XTick', 10.^(-8:2:0));
        xlabel(axTail, '$f_c/\langle f\rangle$', 'Interpreter', 'latex');
        ylabel(axTail, '$\mathrm{fraction\ of\ backbone\ contacts}\ \le f_c$', 'Interpreter', 'latex');
        title(axTail, sprintf('$\\mathrm{(a)\\ %s:\\ the\\ weak\\ tail}$', strLaw), 'Interpreter', 'latex');

        % (b)-(d): one square per packing, x = pressure (spread by mu)
        for idxK = vecIdxLaw
            s = sUncert(idxK);
            vecCol = colourOf(s.scalPressure);
            scalSize = vecSizeMu(vecMuList == s.scalMu);
            scalX = s.scalPressure * 10^(0.07*(find(vecMuList == s.scalMu) - (numel(vecMuList) + 1)/2));
            plot(axCount, scalX, s.scalNumWeak, 's', 'MarkerSize', scalSize, 'Color', vecCol, 'LineWidth', 1.1, 'HandleVisibility', 'off');
            plot(axCount, scalX, s.scalNumWeakBackbone, 's', 'MarkerSize', scalSize, 'Color', vecCol, 'MarkerFaceColor', vecCol, 'HandleVisibility', 'off');
            plot(axStiff, scalX, 100*s.scalStiffShareWeak, 's', 'MarkerSize', scalSize, 'Color', vecCol, 'MarkerFaceColor', vecCol, 'HandleVisibility', 'off');
            plot(axDz, scalX, 100*s.scalDeltaZRel, 's', 'MarkerSize', scalSize, 'Color', vecCol, 'MarkerFaceColor', vecCol, 'HandleVisibility', 'off');
        end
        ylabel(axCount, '$N_{\rm weak}\ \mathrm{per\ packing}$', 'Interpreter', 'latex');
        title(axCount, sprintf('$\\mathrm{(b)\\ how\\ many:}\\ f_c<%s\\langle f\\rangle$', latexPow(scalFloor)), 'Interpreter', 'latex');
        ylabel(axStiff, '$\sum_{\rm weak} k_n \,/\, \sum_{\rm backbone} k_n\ (\%)$', 'Interpreter', 'latex');
        title(axStiff, '$\mathrm{(c)\ share\ of\ the\ backbone\ stiffness}$', 'Interpreter', 'latex');
        ylabel(axDz, '$(z - z_{\rm cut})\,/\,(z - z_{\rm iso})\ (\%)$', 'Interpreter', 'latex');
        title(axDz, '$\mathrm{(d)\ cutting\ them:\ shift\ of}\ z - z_{\rm iso}$', 'Interpreter', 'latex');
        % y ranges from the data: (d) is signed (cutting can also RAISE z, when
        % it removes low-z disks), and leave headroom above (b) for its legend
        sLaw = sUncert(vecIdxLaw);
        cellY = {[sLaw.scalNumWeak], 100*[sLaw.scalStiffShareWeak], 100*[sLaw.scalDeltaZRel]};
        vecHeadroom = [1.75, 1.15, 1.15];
        cellAxY = {axCount, axStiff, axDz};
        for idxAx = 1:3
            ax = cellAxY{idxAx};
            vecY = cellY{idxAx};
            scalLo = min([0, vecY]);
            scalHi = max([vecY, 0]);
            scalSpan = max(scalHi - scalLo, eps);
            ylim(ax, [scalLo - 0.08*scalSpan*(scalLo < 0), scalLo + vecHeadroom(idxAx)*scalSpan]);
            yline(ax, 0, '-', 'Color', vecGrid, 'LineWidth', 0.6, 'HandleVisibility', 'off');
            vecTick = unique(10.^floor(log10(vecPList) + 1e-12));   % decades only: 3e-5 and 1e-4 labels would collide
            set(ax, 'XScale', 'log', 'XLim', [min(vecPList)/3, max(vecPList)*3], 'XTick', vecTick, ...
                'XTickLabel', arrayfun(@(x) ['$' latexPow(x) '$'], vecTick, 'UniformOutput', false));
            xlabel(ax, '$\mathrm{target\ pressure}\ P$', 'Interpreter', 'latex');
        end

        % legends: pressure (colour) in (a), friction (size) and fill in (b)
        for idxP = numel(vecPList):-1:1
            plot(axTail, NaN, NaN, '-', 'Color', matColourP(idxP, :), 'LineWidth', 2, ...
                'DisplayName', sprintf('$P = %s$', latexPow(vecPList(idxP))));
        end
        legend(axTail, 'Location', 'southeast', 'Interpreter', 'latex', 'FontSize', 10, 'TextColor', vecFg, 'Color', vecBg, 'EdgeColor', vecGrid);
        plot(axCount, NaN, NaN, 's', 'MarkerSize', 8, 'Color', vecFg, 'MarkerFaceColor', vecFg, 'DisplayName', '$\mathrm{backbone}$');
        plot(axCount, NaN, NaN, 's', 'MarkerSize', 8, 'Color', vecFg, 'DisplayName', '$\mathrm{all\ (incl.\ rattlers)}$');
        for idxMu = 1:numel(vecMuList)
            plot(axCount, NaN, NaN, 's', 'MarkerSize', vecSizeMu(idxMu), 'Color', vecGrid, 'DisplayName', sprintf('$\\mu = %g$', vecMuList(idxMu)));
        end
        legend(axCount, 'Location', 'north', 'Interpreter', 'latex', 'FontSize', 9, 'TextColor', vecFg, 'Color', vecBg, 'EdgeColor', vecGrid, 'NumColumns', 3);

        for ax = [axTail, axCount, axStiff, axDz]
            set(ax, 'Color', vecBg, 'XColor', vecFg, 'YColor', vecFg, 'GridColor', vecGrid, ...
                'TickLabelInterpreter', 'latex', 'FontSize', 12, 'Box', 'on');
            ax.Title.Color = vecFg;
            grid(ax, 'on');
            axis(ax, 'square');
        end
    end

    if strlength(options.strSavePath) > 0
        exportgraphics(hFig, options.strSavePath, 'Resolution', 150, 'BackgroundColor', vecBg);
    end

    function vecCol = colourOf(scalP)
        if isnan(scalP)
            vecCol = vecGrid;
        else
            vecCol = matColourP(vecPList == scalP, :);
        end
    end
end

function str = latexPow(scalX)
% 1e-4 -> '10^{-4}', 3e-5 -> '3\times10^{-5}' (for labels inside $...$)
    scalExp = floor(log10(scalX) + 1e-12);
    scalMant = scalX / 10^scalExp;
    if abs(scalMant - 1) < 1e-9
        str = sprintf('10^{%d}', scalExp);
    else
        str = sprintf('%g\\times10^{%d}', scalMant, scalExp);
    end
end
