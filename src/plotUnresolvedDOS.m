function sSummary = plotUnresolvedDOS(options)
% plotUnresolvedDOS  Figures: how unresolved particles bias the frictional DOS.
%
%   plotUnresolvedDOS()
%   plotUnresolvedDOS("boolLightMode", false, "strFigureDir", "figs")
%   sSummary = plotUnresolvedDOS("vecStrCompareDir", [dirSeeds dirTol dirN2000], ...
%       "vecStrCompareLabel", ["6 seeds", "tolerance 1e-5", "N = 2000"])
%
%   Reads the runUnresolvedDOS output (<strDataDir>/<law>/unresolved_<rule>.mat)
%   and writes, into strFigureDir:
%
%     fig0_tail_before.png, fig0_tail_after.png
%                          plotUncertContacts of every packing before and
%                          after cleanRats' unresolved-particle rule
%                          (sUncertOptions.boolRemoveUnresolved): the
%                          uncertain tail of the contact forces disappears
%                          (boolTail; reads the packings, about 1 min)
%     fig1_dos.png         DOS with and without the unresolved particles,
%                          Hooke (top) and Hertz (bottom), one column per
%                          (mu, P) of matDosGroups, Frictional solve
%     fig2_collapse.png    omega* against z - 3 with and without them, with
%                          linear fits: the Hertz collapse tightens, Hooke
%                          does not move
%     fig3_estimators.png  the omega* shift under four estimators (20% and
%                          50% of D(omega_ref), equal-count bins, knee fit)
%     fig4_modes.png       every mode of the original Frictional matrix:
%                          share on the removed particles against omega_hat,
%                          coloured by how many particles it moves
%     fig5_mechanism.png   (a,b) k_n/<k_n> against f_n/<f> (Hooke, Hertz),
%                          contacts of removed particles highlighted; (c) the
%                          Coulomb cone test tan(beta) <= mu of two-contact
%                          particles, kept and removed
%     fig6_robustness.png  the Hertz omega* shift per (mu, P) for the base data
%                          set and every vecStrCompareDir (only if given)
%     fig7_rules.png       the narrow rule ("unresolved") against the broad
%                          rule ("twoContact", every z <= 2 particle): particles
%                          removed, collapsed packings, omega* shift (only if
%                          unresolved_twoContact.mat exists)
%
%   omega_hat = omega/sqrt(<k_n>/m); seeds of the same (mu, P) are pooled.
%   omega* estimators (sGroups fields):
%     T20  D(omega*)/D(omega_ref) = 0.2, bins scalBinWidth (Henkes et al.,
%          plotFrictionalDOS' default)
%     T50  the same at 0.5, bins scalBinWidth/2
%     Eq   T20 on equal-count bins of scalModesPerBin modes
%     Knee least-squares fit of D = D_p omega/(omega + omega_k) below
%          omega_ref, bins scalBinWidth/2; omega* = omega_k (NaN when the
%          best omega_k is the largest trial value: no knee below omega_ref)
%
%   Options:
%       strDataDir        base data set (default <GranE>/data/frictionalDOS)
%       strFigureDir      output (default <strDataDir>/figures_unresolved)
%       boolLightMode     light (default) or dark figures
%       boolTail          make fig0 (default true)
%       matDosGroups      [mu P] rows for fig1 (default [0.01 1e-4; 0.1 1e-3; 0.3 1e-4])
%       vecStrCompareDir  extra data sets for fig6 (hertz only); a data set
%                         whose (mu, P, seed) rows repeat the base set's is
%                         compared on its own, otherwise its seeds are pooled
%                         with the base set's ("seeds" comparisons)
%       vecStrCompareLabel labels for them
%       vecBoolComparePool true: pool that data set with the base set
%       scalBinWidth      (default 0.08), scalOmegaRef (default 1),
%       scalModesPerBin   (default 40)
%
%   Returns sSummary.(law): sGroups (per-(mu, P) omega* with and without,
%   every estimator) and the fits of fig2.
%
%   The sideways jitter that separates overlapping points (fig5 panel c,
%   fig7 panels a-b) is drawn from a private random stream with a fixed
%   seed, so every run writes identical PNGs and the caller's global random
%   state is left untouched.

    arguments
        options.strDataDir   (1, 1) string = ""
        options.strFigureDir (1, 1) string = ""
        options.boolLightMode (1, 1) logical = true
        options.boolTail      (1, 1) logical = true
        options.matDosGroups (:, 2) double = [0.01 1e-4; 0.1 1e-3; 0.3 1e-4]
        options.vecStrCompareDir   (1, :) string = string.empty(1, 0)
        options.vecStrCompareLabel (1, :) string = string.empty(1, 0)
        options.vecBoolComparePool (1, :) logical = false(1, 0)
        options.scalBinWidth    (1, 1) double {mustBePositive} = 0.08
        options.scalOmegaRef    (1, 1) double {mustBePositive} = 1
        options.scalModesPerBin (1, 1) double {mustBePositive, mustBeInteger} = 40
    end

    % Fixed-seed stream for the plotting jitter (reproducible figures).
    sJitterStream = RandStream('mt19937ar', 'Seed', 0);

    strRepoRoot = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(strRepoRoot, 'src'));
    strDataDir = options.strDataDir;
    if strDataDir == ""
        strDataDir = fullfile(strRepoRoot, 'data', 'frictionalDOS');
    end
    strFigureDir = options.strFigureDir;
    if strFigureDir == ""
        strFigureDir = fullfile(strDataDir, 'figures_unresolved');
    end
    if ~isfolder(strFigureDir)
        mkdir(strFigureDir);
    end
    sTheme = makeTheme(options.boolLightMode);
    vecLaws = ["hooke", "hertz"];
    cellLawName = {'Hooke', 'Hertz'};

    sSummary = struct();
    for strLaw = vecLaws
        sRes = loadResults(strDataDir, strLaw, "unresolved");
        sSummary.(strLaw).sResults = sRes;
        sSummary.(strLaw).sGroups = groupResults(sRes, options);
    end

    % ---- fig0: the uncertain tail before and after the rule -------------------
    if options.boolTail
        [sBefore, sAfter] = uncertainTail(strDataDir, vecLaws);
        hFig = plotUncertContacts(sBefore, 'boolLightMode', options.boolLightMode, 'strTitle', "before removing unresolved particles", ...
            'strSavePath', string(fullfile(strFigureDir, 'fig0_tail_before.png')));
        close(hFig);
        hFig = plotUncertContacts(sAfter, 'boolLightMode', options.boolLightMode, 'strTitle', "after removing unresolved particles", ...
            'strSavePath', string(fullfile(strFigureDir, 'fig0_tail_after.png')));
        close(hFig);
        sSummary.sUncertBefore = sBefore;
        sSummary.sUncertAfter = sAfter;
    end

    % ---- fig1: DOS with and without ------------------------------------------
    scalNumCols = size(options.matDosGroups, 1);
    hFig = newFigure(sTheme, [60 60 430*scalNumCols 820]);
    for idxLaw = 1:2
        sG = sSummary.(vecLaws(idxLaw)).sGroups;
        for idxCol = 1:scalNumCols
            ax = subplot(2, scalNumCols, (idxLaw - 1)*scalNumCols + idxCol); hold(ax, 'on');
            idxG = find([sG.scalMu] == options.matDosGroups(idxCol, 1) & [sG.scalPressureTarget] == options.matDosGroups(idxCol, 2), 1);
            if isempty(idxG)
                axis(ax, 'off'); continue;
            end
            g = sG(idxG);
            scalRef = interp1(g.vecBinCenters, g.vecDOSStick, options.scalOmegaRef);
            scalRefCut = interp1(g.vecBinCenters, g.vecDOSStickCut, options.scalOmegaRef);
            plot(ax, g.vecBinCenters, g.vecDOSStick/scalRef, '-o', 'Color', sTheme.vecKept, 'MarkerSize', 4, 'LineWidth', 1.3, ...
                'DisplayName', '$\mathrm{all\ particles}$');
            plot(ax, g.vecBinCenters, g.vecDOSStickCut/scalRefCut, '-s', 'Color', sTheme.vecCut, 'MarkerSize', 4, 'LineWidth', 1.3, ...
                'DisplayName', '$\mathrm{unresolved\ removed}$');
            xline(ax, g.scalT20, '--', 'Color', sTheme.vecKept, 'LineWidth', 1.1, 'HandleVisibility', 'off');
            xline(ax, g.scalT20Cut, '--', 'Color', sTheme.vecCut, 'LineWidth', 1.1, 'HandleVisibility', 'off');
            yline(ax, 0.2, ':', 'Color', sTheme.vecFg, 'HandleVisibility', 'off');
            set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', [0.03 3], 'YLim', [0.02 3]);
            xlabel(ax, '$\hat\omega$', 'Interpreter', 'latex');
            ylabel(ax, '$D(\hat\omega)/D(\hat\omega_{\mathrm{ref}})$', 'Interpreter', 'latex');
            title(ax, sprintf('$\\mathrm{(%s)\\ %s},\\ \\mu=%g,\\ P=%s$', char('a' + (idxLaw - 1)*scalNumCols + idxCol - 1), ...
                cellLawName{idxLaw}, g.scalMu, latexPow(g.scalPressureTarget)), 'Interpreter', 'latex');
            subtitle(ax, sprintf('$\\omega^*: %.3f \\rightarrow %.3f\\ (%d\\ \\mathrm{removed})$', g.scalT20, g.scalT20Cut, g.scalNumRemoved), ...
                'Interpreter', 'latex', 'Color', sTheme.vecFg);
            if idxCol == 1
                legend(ax, 'Location', 'southeast', 'Interpreter', 'latex', 'FontSize', 10, 'TextColor', sTheme.vecFg, 'Color', sTheme.vecBg, 'EdgeColor', sTheme.vecGrid);
            end
            styleAxes(ax, sTheme);
        end
    end
    saveFigure(hFig, fullfile(strFigureDir, 'fig1_dos.png'), sTheme);

    % ---- fig2: omega* against z - 3 ------------------------------------------
    hFig = newFigure(sTheme, [60 60 1100 500]);
    for idxLaw = 1:2
        sG = sSummary.(vecLaws(idxLaw)).sGroups;
        ax = subplot(1, 2, idxLaw); hold(ax, 'on');
        [vecFitKept, scalR2Kept] = linearFit([sG.scalDeltaZ], [sG.scalT20]);
        [vecFitCut, scalR2Cut] = linearFit([sG.scalDeltaZCut], [sG.scalT20Cut]);
        for idxG = 1:numel(sG)   % arrows: how far each group moves
            plot(ax, [sG(idxG).scalDeltaZ sG(idxG).scalDeltaZCut], [sG(idxG).scalT20 sG(idxG).scalT20Cut], '-', 'Color', sTheme.vecGrid, 'HandleVisibility', 'off');
        end
        plot(ax, [sG.scalDeltaZ], [sG.scalT20], 'o', 'Color', sTheme.vecKept, 'MarkerSize', 7, 'LineWidth', 1.2, ...
            'DisplayName', sprintf('$\\mathrm{all:}\\ R^2 = %.3f,\\ \\mathrm{slope}\\ %.3f$', scalR2Kept, vecFitKept(1)));
        plot(ax, [sG.scalDeltaZCut], [sG.scalT20Cut], 's', 'Color', sTheme.vecCut, 'MarkerFaceColor', sTheme.vecCut, 'MarkerSize', 7, ...
            'DisplayName', sprintf('$\\mathrm{removed:}\\ R^2 = %.3f,\\ \\mathrm{slope}\\ %.3f$', scalR2Cut, vecFitCut(1)));
        vecX = linspace(0, max([sG.scalDeltaZ sG.scalDeltaZCut], [], 'omitnan')*1.05, 20);
        plot(ax, vecX, polyval(vecFitKept, vecX), '--', 'Color', sTheme.vecKept, 'HandleVisibility', 'off');
        plot(ax, vecX, polyval(vecFitCut, vecX), '-', 'Color', sTheme.vecCut, 'HandleVisibility', 'off');
        xlabel(ax, '$\Delta z = z - 3$', 'Interpreter', 'latex');
        ylabel(ax, '$\omega^*\ (\hat\omega)$', 'Interpreter', 'latex');
        title(ax, sprintf('$\\mathrm{(%s)\\ %s,\\ Frictional\\ solve}$', char('a' + idxLaw - 1), cellLawName{idxLaw}), 'Interpreter', 'latex');
        legend(ax, 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 10, 'TextColor', sTheme.vecFg, 'Color', sTheme.vecBg, 'EdgeColor', sTheme.vecGrid);
        styleAxes(ax, sTheme);
        sSummary.(vecLaws(idxLaw)).vecFitKept = vecFitKept; sSummary.(vecLaws(idxLaw)).scalR2Kept = scalR2Kept;
        sSummary.(vecLaws(idxLaw)).vecFitCut = vecFitCut;   sSummary.(vecLaws(idxLaw)).scalR2Cut = scalR2Cut;
    end
    saveFigure(hFig, fullfile(strFigureDir, 'fig2_collapse.png'), sTheme);

    % ---- fig3: the shift under four estimators -------------------------------
    cellEst = {'T20', 'T50', 'Eq', 'Knee'};
    cellEstName = {'$20\%\ \mathrm{of}\ D(\hat\omega_{\mathrm{ref}})$', '$50\%$', '$\mathrm{equal\ count}$', '$\mathrm{knee\ fit}$'};
    cellEstMk = {'o', 's', '^', 'd'};
    hFig = newFigure(sTheme, [60 60 1100 500]);
    for idxLaw = 1:2
        sG = sSummary.(vecLaws(idxLaw)).sGroups;
        ax = subplot(1, 2, idxLaw); hold(ax, 'on');
        for idxE = 1:numel(cellEst)
            vecShift = 100*([sG.(['scal' cellEst{idxE} 'Cut'])]./[sG.(['scal' cellEst{idxE}])] - 1);
            plot(ax, [sG.scalDeltaZ] + 0.006*(idxE - 2.5), vecShift, cellEstMk{idxE}, 'Color', sTheme.matSeries(idxE, :), ...
                'MarkerFaceColor', sTheme.matSeries(idxE, :), 'MarkerSize', 7, 'DisplayName', cellEstName{idxE});
        end
        yline(ax, 0, '-', 'Color', sTheme.vecGrid, 'HandleVisibility', 'off');
        xlabel(ax, '$\Delta z = z - 3\ \mathrm{(all\ particles)}$', 'Interpreter', 'latex');
        ylabel(ax, '$\omega^*_{\mathrm{removed}}/\omega^*_{\mathrm{all}} - 1\ (\%)$', 'Interpreter', 'latex');
        title(ax, sprintf('$\\mathrm{(%s)\\ %s:\\ shift\\ of}\\ \\omega^*\\ \\mathrm{per}\\ (\\mu, P)$', char('a' + idxLaw - 1), cellLawName{idxLaw}), 'Interpreter', 'latex');
        legend(ax, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 10, 'TextColor', sTheme.vecFg, 'Color', sTheme.vecBg, 'EdgeColor', sTheme.vecGrid);
        ax.YAxis.Exponent = 0;
        styleAxes(ax, sTheme);
    end
    saveFigure(hFig, fullfile(strFigureDir, 'fig3_estimators.png'), sTheme);

    % ---- fig4: mode localization ----------------------------------------------
    hFig = newFigure(sTheme, [60 60 1100 500]);
    for idxLaw = 1:2
        sRes = sSummary.(vecLaws(idxLaw)).sResults;
        sRes = sRes([sRes.scalNumRemoved] > 0 & isfield(sRes, 'vecModeOmegaHat'));
        ax = subplot(1, 2, idxLaw); hold(ax, 'on');
        vecW = vertcat(sRes.vecModeOmegaHat);
        vecShare = vertcat(sRes.vecModeShareRemoved);
        vecNumMoved = vertcat(sRes.vecModePR).*repelem([sRes.scalNi].', arrayfun(@(s) numel(s.vecModeOmegaHat), sRes));
        vecBool = vecW > 1e-3;
        scatter(ax, vecW(vecBool), vecShare(vecBool), 10, log10(vecNumMoved(vecBool)), 'filled', 'MarkerFaceAlpha', 0.6);
        colormap(ax, turbo);
        clim(ax, [0 log10(max([sRes.scalNi]))]);
        hBar = colorbar(ax);
        hBar.Label.String = '$\log_{10}(\mathrm{particles\ moved} = N_i \times \mathrm{PR})$';
        hBar.Label.Interpreter = 'latex'; hBar.TickLabelInterpreter = 'latex'; hBar.Color = sTheme.vecFg;
        yline(ax, 0.5, '--', 'Color', sTheme.vecFg, 'HandleVisibility', 'off');
        scalNumLoc = nnz(vecBool & vecShare > 0.5);
        scalNumRem = sum([sRes.scalNumRemoved]);
        set(ax, 'XScale', 'log', 'XLim', [1e-2 5], 'YLim', [0 1]);
        xlabel(ax, '$\hat\omega\ \mathrm{(original\ Frictional\ matrix)}$', 'Interpreter', 'latex');
        ylabel(ax, '$\mathrm{share\ of\ the\ mode\ on\ removed\ particles}$', 'Interpreter', 'latex');
        title(ax, sprintf('$\\mathrm{(%s)\\ %s}$', char('a' + idxLaw - 1), cellLawName{idxLaw}), 'Interpreter', 'latex');
        subtitle(ax, {sprintf('$%d\\ \\mathrm{modes\\ with\\ share} > 0.5,\\ %d\\ \\mathrm{removed\\ particles}$', scalNumLoc, scalNumRem), ...
            sprintf('$\\mathrm{their\\ median}\\ \\hat\\omega = %.2f$', median(vecW(vecBool & vecShare > 0.5)))}, ...
            'Interpreter', 'latex', 'Color', sTheme.vecFg, 'FontSize', 10);
        styleAxes(ax, sTheme);
        sSummary.(vecLaws(idxLaw)).scalNumLocalizedModes = scalNumLoc;
        sSummary.(vecLaws(idxLaw)).scalNumRemoved = scalNumRem;
    end
    saveFigure(hFig, fullfile(strFigureDir, 'fig4_modes.png'), sTheme);

    % ---- fig5: mechanism (stiffness of weak contacts; cone test) -------------
    hFig = newFigure(sTheme, [60 60 1500 480]);
    for idxLaw = 1:2
        sRes = sSummary.(vecLaws(idxLaw)).sResults;
        ax = subplot(1, 3, idxLaw); hold(ax, 'on');
        vecF = vertcat(sRes.vecForceRel); vecK = vertcat(sRes.vecKnRel); vecBoolRem = vertcat(sRes.vecBoolContactRemoved);
        plot(ax, vecF(~vecBoolRem), vecK(~vecBoolRem), '.', 'Color', sTheme.vecKept, 'MarkerSize', 4, 'DisplayName', '$\mathrm{kept\ contacts}$');
        plot(ax, vecF(vecBoolRem), vecK(vecBoolRem), 'o', 'Color', sTheme.vecCut, 'MarkerSize', 5, 'LineWidth', 1.1, 'DisplayName', '$\mathrm{contacts\ of\ removed\ particles}$');
        if idxLaw == 2
            vecX = logspace(-9, 1, 50);
            plot(ax, vecX, vecX.^(1/3)*median(vecK(vecF > 0.5 & vecF < 2)), '--', 'Color', sTheme.vecFg, 'DisplayName', '$k_n \propto f_n^{1/3}$');
        end
        xline(ax, 1e-4, ':', 'Color', sTheme.vecFg, 'HandleVisibility', 'off');
        set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', [1e-9 10], 'YLim', [1e-4 10], 'XTick', 10.^(-8:2:0));
        xlabel(ax, '$f_n/\langle f\rangle$', 'Interpreter', 'latex');
        ylabel(ax, '$k_n/\langle k_n\rangle$', 'Interpreter', 'latex');
        title(ax, sprintf('$\\mathrm{(%s)\\ %s:\\ contact\\ stiffness}$', char('a' + idxLaw - 1), cellLawName{idxLaw}), 'Interpreter', 'latex');
        legend(ax, 'Location', 'southeast', 'Interpreter', 'latex', 'FontSize', 9, 'TextColor', sTheme.vecFg, 'Color', sTheme.vecBg, 'EdgeColor', sTheme.vecGrid);
        styleAxes(ax, sTheme);
    end
    ax = subplot(1, 3, 3); hold(ax, 'on');
    sAll = [sSummary.hooke.sResults, sSummary.hertz.sResults];
    vecMuList = unique([sAll.scalMu]);
    for idxMu = 1:numel(vecMuList)
        sMu = sAll([sAll.scalMu] == vecMuList(idxMu));
        vecKept = vertcat(sMu.vecTanBetaKept); vecRem = vertcat(sMu.vecTanBetaRemoved);
        scalX = idxMu;
        plot(ax, scalX - 0.15 + 0.1*(rand(sJitterStream, size(vecKept)) - 0.5), vecKept, '.', 'Color', sTheme.vecKept, 'MarkerSize', 5, 'HandleVisibility', 'off');
        plot(ax, scalX + 0.15 + 0.1*(rand(sJitterStream, size(vecRem)) - 0.5), vecRem, 'o', 'Color', sTheme.vecCut, 'MarkerSize', 5, 'LineWidth', 1.1, 'HandleVisibility', 'off');
        plot(ax, scalX + [-0.4 0.4], [vecMuList(idxMu) vecMuList(idxMu)], '-', 'Color', sTheme.vecFg, 'LineWidth', 2, 'HandleVisibility', 'off');
        text(ax, scalX, 1.4e-3, sprintf('$%d/%d$', nnz(vecRem > vecMuList(idxMu)), numel(vecRem)), 'Interpreter', 'latex', ...
            'HorizontalAlignment', 'center', 'Color', sTheme.vecCut, 'FontSize', 10);
    end
    plot(ax, NaN, NaN, '.', 'Color', sTheme.vecKept, 'MarkerSize', 12, 'DisplayName', '$\mathrm{kept\ two\mbox{-}contact}$');
    plot(ax, NaN, NaN, 'o', 'Color', sTheme.vecCut, 'DisplayName', '$\mathrm{removed\ (unresolved)}$');
    plot(ax, NaN, NaN, '-', 'Color', sTheme.vecFg, 'LineWidth', 2, 'DisplayName', '$\tan\beta = \mu$');
    set(ax, 'YScale', 'log', 'XLim', [0.4 numel(vecMuList) + 0.6], 'YLim', [1e-3 1e3], 'XTick', 1:numel(vecMuList), ...
        'XTickLabel', arrayfun(@(x) sprintf('$%g$', x), vecMuList, 'UniformOutput', false));
    xlabel(ax, '$\mu$', 'Interpreter', 'latex');
    ylabel(ax, '$\tan\beta\ \mathrm{of\ two\mbox{-}contact\ particles}$', 'Interpreter', 'latex');
    title(ax, '$\mathrm{(c)\ Coulomb\ cone\ test:\ load\ needs}\ \tan\beta \le \mu$', 'Interpreter', 'latex');
    subtitle(ax, '$\mathrm{numbers:\ removed\ particles\ failing\ the\ test}$', 'Interpreter', 'latex', 'Color', sTheme.vecFg);
    legend(ax, 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 9, 'TextColor', sTheme.vecFg, 'Color', sTheme.vecBg, 'EdgeColor', sTheme.vecGrid);
    styleAxes(ax, sTheme);
    saveFigure(hFig, fullfile(strFigureDir, 'fig5_mechanism.png'), sTheme);

    % ---- fig6: robustness (Hertz) --------------------------------------------
    if ~isempty(options.vecStrCompareDir)
        hFig = newFigure(sTheme, [60 60 1500 520]);
        sBase = sSummary.hertz.sResults;
        cellSets = {sSummary.hertz.sGroups};
        cellLabel = {'$\mathrm{base}$'};
        for idxC = 1:numel(options.vecStrCompareDir)
            sExtra = loadResults(options.vecStrCompareDir(idxC), "hertz", "unresolved");
            strLabel = char(replace(options.vecStrCompareLabel(idxC), " ", "\ "));
            if numel(options.vecBoolComparePool) >= idxC && options.vecBoolComparePool(idxC)
                sExtra = [sBase, sExtra]; %#ok<AGROW>
            else
                % a fair baseline: the base packings with the same (mu, P, seed)
                % when the comparison covers only some of the base seeds
                vecBoolSame = arrayfun(@(b) any([sExtra.scalMu] == b.scalMu & [sExtra.scalPressureTarget] == b.scalPressureTarget & ...
                    [sExtra.scalSeed] == b.scalSeed), sBase);
                if any(vecBoolSame) && ~all(ismember(unique([sBase.scalSeed]), unique([sExtra.scalSeed])))
                    cellSets{end + 1} = groupResults(sBase(vecBoolSame), options); %#ok<AGROW>
                    cellLabel{end + 1} = ['$\mathrm{base,\ same\ packings\ as\ ' strLabel '}$']; %#ok<AGROW>
                end
            end
            cellSets{end + 1} = groupResults(sExtra, options); %#ok<AGROW>
            cellLabel{end + 1} = ['$\mathrm{' strLabel '}$']; %#ok<AGROW>
        end
        for idxE = 1:2
            ax = subplot(1, 2, idxE); hold(ax, 'on');
            strE = cellEst{1 + 2*(idxE - 1)};   % T20, Eq
            matKeys = unique(cell2mat(cellfun(@(g) [[g.scalMu].', [g.scalPressureTarget].'], cellSets, 'UniformOutput', false).'), 'rows');
            vecBoolHasData = false(size(matKeys, 1), 1);   % drop groups where omega* is unmeasurable everywhere
            for idxS = 1:numel(cellSets)
                g = cellSets{idxS};
                vecBoolOk = ~isnan([g.(['scal' strE 'Cut'])]./[g.(['scal' strE])]);
                matOk = [reshape([g(vecBoolOk).scalMu], [], 1), reshape([g(vecBoolOk).scalPressureTarget], [], 1)];
                vecBoolHasData = vecBoolHasData | ismember(matKeys, matOk, 'rows');
            end
            matKeys = matKeys(vecBoolHasData, :);
            for idxS = 1:numel(cellSets)
                g = cellSets{idxS};
                [vecBoolIn, vecX] = ismember([[g.scalMu].', [g.scalPressureTarget].'], matKeys, 'rows');
                vecShift = 100*([g.(['scal' strE 'Cut'])]./[g.(['scal' strE])] - 1);
                vecX = vecX(vecBoolIn).'; vecShift = vecShift(vecBoolIn.');
                plot(ax, vecX + 0.1*(idxS - (numel(cellSets) + 1)/2), vecShift, cellEstMk{1 + mod(idxS - 1, 4)}, 'Color', sTheme.matSeries(idxS, :), ...
                    'MarkerFaceColor', sTheme.matSeries(idxS, :), 'MarkerSize', 8, 'DisplayName', cellLabel{idxS});
            end
            yline(ax, 0, '-', 'Color', sTheme.vecGrid, 'HandleVisibility', 'off');
            set(ax, 'XTick', 1:size(matKeys, 1), 'XTickLabel', arrayfun(@(i) sprintf('$%g,\\ %s$', matKeys(i, 1), latexPow(matKeys(i, 2))), ...
                1:size(matKeys, 1), 'UniformOutput', false), 'XLim', [0.4 size(matKeys, 1) + 0.6]);
            xtickangle(ax, 45);
            xlabel(ax, '$(\mu,\ P)$', 'Interpreter', 'latex');
            ylabel(ax, '$\omega^*\ \mathrm{shift}\ (\%)$', 'Interpreter', 'latex');
            title(ax, sprintf('$\\mathrm{(%s)\\ Hertz,\\ %s}$', char('a' + idxE - 1), ...
                char(replace(string(cellEstName{1 + 2*(idxE - 1)}), '$', ''))), 'Interpreter', 'latex');
            if idxE == 2
                legend(ax, 'Location', 'eastoutside', 'Interpreter', 'latex', 'FontSize', 10, 'TextColor', sTheme.vecFg, 'Color', sTheme.vecBg, 'EdgeColor', sTheme.vecGrid);
            end
            ax.YAxis.Exponent = 0;
            styleAxes(ax, sTheme);
        end
        saveFigure(hFig, fullfile(strFigureDir, 'fig6_robustness.png'), sTheme);
        sSummary.cellCompareGroups = cellSets;
    end

    % ---- fig7: narrow rule against the broad one -----------------------------
    if isfile(fullfile(strDataDir, 'hertz', 'unresolved_twoContact.mat'))
        hFig = newFigure(sTheme, [60 60 1500 480]);
        for idxLaw = 1:2
            sNarrow = sSummary.(vecLaws(idxLaw)).sResults;
            sBroad = loadResults(strDataDir, vecLaws(idxLaw), "twoContact");
            sSummary.(vecLaws(idxLaw)).sResultsTwoContact = sBroad;
            ax = subplot(1, 3, idxLaw); hold(ax, 'on');
            plot(ax, [sNarrow.scalPressureTarget].*10.^(0.05*randn(sJitterStream, 1, numel(sNarrow))), 100*[sNarrow.scalNumRemoved]./[sNarrow.scalNi], 's', ...
                'Color', sTheme.vecCut, 'MarkerFaceColor', sTheme.vecCut, 'MarkerSize', 6, 'DisplayName', '$\mathrm{unresolved\ (narrow)}$');
            plot(ax, [sBroad.scalPressureTarget].*10.^(0.05*randn(sJitterStream, 1, numel(sBroad))), 100*[sBroad.scalNumRemoved]./[sBroad.scalNi], 'o', ...
                'Color', sTheme.vecKept, 'MarkerSize', 6, 'LineWidth', 1.1, 'DisplayName', '$z_i \le 2\ \mathrm{(broad)}$');
            vecPTick = 10.^(floor(log10(min([sBroad.scalPressureTarget]))):ceil(log10(max([sBroad.scalPressureTarget]))));
            set(ax, 'XScale', 'log', 'YScale', 'log', 'YLim', [0.05 150], 'XLim', [vecPTick(1)/1.5 vecPTick(end)*1.5], 'XTick', vecPTick, ...
                'XTickLabel', arrayfun(@(x) ['$' latexPow(x) '$'], vecPTick, 'UniformOutput', false));
            xlabel(ax, '$\mathrm{target\ pressure}\ P$', 'Interpreter', 'latex');
            ylabel(ax, '$\mathrm{particles\ removed\ (\%\ of\ backbone)}$', 'Interpreter', 'latex');
            title(ax, sprintf('$\\mathrm{(%s)\\ %s}$', char('a' + idxLaw - 1), cellLawName{idxLaw}), 'Interpreter', 'latex');
            subtitle(ax, sprintf('$\\mathrm{broad\\ rule\\ dismantles}\\ %d/%d\\ \\mathrm{packings}$', nnz([sBroad.boolCollapsed]), numel(sBroad)), ...
                'Interpreter', 'latex', 'Color', sTheme.vecFg);
            legend(ax, 'Location', 'southwest', 'Interpreter', 'latex', 'FontSize', 9, 'TextColor', sTheme.vecFg, 'Color', sTheme.vecBg, 'EdgeColor', sTheme.vecGrid);
            styleAxes(ax, sTheme);
        end
        ax = subplot(1, 3, 3); hold(ax, 'on');
        gN = sSummary.hertz.sGroups;
        gB = groupResults(sSummary.hertz.sResultsTwoContact, options);
        plot(ax, [gN.scalDeltaZ], 100*([gN.scalT20Cut]./[gN.scalT20] - 1), 's', 'Color', sTheme.vecCut, 'MarkerFaceColor', sTheme.vecCut, ...
            'MarkerSize', 7, 'DisplayName', '$\mathrm{unresolved\ (narrow)}$');
        plot(ax, [gB.scalDeltaZ], 100*([gB.scalT20Cut]./[gB.scalT20] - 1), 'o', 'Color', sTheme.vecKept, 'MarkerSize', 7, 'LineWidth', 1.1, ...
            'DisplayName', '$z_i \le 2\ \mathrm{(broad)}$');
        yline(ax, 0, '-', 'Color', sTheme.vecGrid, 'HandleVisibility', 'off');
        xlabel(ax, '$\Delta z = z - 3\ \mathrm{(all\ particles)}$', 'Interpreter', 'latex');
        ylabel(ax, '$\omega^*\ \mathrm{shift}\ (\%)$', 'Interpreter', 'latex');
        title(ax, '$\mathrm{(c)\ Hertz:}\ \omega^*\ \mathrm{shift,\ surviving\ packings}$', 'Interpreter', 'latex');
        legend(ax, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 9, 'TextColor', sTheme.vecFg, 'Color', sTheme.vecBg, 'EdgeColor', sTheme.vecGrid);
        ax.YAxis.Exponent = 0;
        styleAxes(ax, sTheme);
        saveFigure(hFig, fullfile(strFigureDir, 'fig7_rules.png'), sTheme);
        sSummary.hertz.sGroupsTwoContact = gB;
    end
end

function [sBefore, sAfter] = uncertainTail(strDataDir, vecLaws)
% cleanRats on every packing's full pre-cleanRats state, with the
% unresolved-particle rule on: its 4th/5th outputs are the uncertain-contact
% statistics after and before the rule (plotUncertContacts' input).
    sBefore = struct([]); sAfter = struct([]);
    for strLaw = vecLaws
        sSpectra = load(fullfile(strDataDir, strLaw, 'spectra.mat'), 'sRuns');
        for sRun = sSpectra.sRuns
            scalSeed = sRun.scalSeed;
            if isfield(sRun, 'scalSeedUsed') && ~isempty(sRun.scalSeedUsed)
                scalSeed = sRun.scalSeedUsed;
            end
            strSuffix = '_Friction';
            if strLaw == "hertz"
                strSuffix = '_Hertz_Friction';
            end
            sFile = load(fullfile(strDataDir, strLaw, 'packings', sprintf('mu%g', sRun.scalMu), ...
                sprintf('2D_N%d_P%s_Width%d_Seed%d%s.mat', sRun.scalN, num2str(sRun.scalPressureTarget), ...
                round(sqrt(sRun.scalN)), scalSeed, strSuffix)), 'sFullState');
            sF = sFile.sFullState; %#ok<NASGU> used inside evalc (cleanRats prints per packing)
            sOpt = struct('strContactLaw', char(strLaw), 'scalPressure', sRun.scalPressureTarget, 'boolRemoveUnresolved', true); %#ok<NASGU>
            [~, ~, ~, ~, sA, sB] = evalc(['cleanRats([sF.vecPosX(:), sF.vecPosY(:)], sF.vecDiameter(:)/2, sF.scalBoxHeightY, ' ...
                'sF.scalBoxWidthX, [], true, sRun.scalMu, sOpt)']);
            if isempty(sBefore)
                sBefore = sB; sAfter = sA;
            else
                sBefore(end + 1) = sB; sAfter(end + 1) = sA; %#ok<AGROW>
            end
        end
    end
end

function sRes = loadResults(strDataDir, strLaw, strRule)
    strFile = fullfile(strDataDir, strLaw, sprintf('unresolved_%s.mat', strRule));
    if ~isfile(strFile)
        error('plotUnresolvedDOS:NoResults', 'No %s. Run runUnresolvedDOS("%s", "strDataDir", "%s", "strRule", "%s") first.', ...
            strFile, strLaw, strDataDir, strRule);
    end
    sLoaded = load(strFile, 'sResults');
    sRes = sLoaded.sResults;
end

function sGroups = groupResults(sRes, options)
% Pool the seeds of every (mu, P); DOS and omega* (four estimators) of the
% Frictional solve with all particles and with the removed ones gone.
% Collapsed packings are left out of both versions.
    sRes = sRes(~[sRes.boolCollapsed]);
    matKeys = unique([[sRes.scalMu].', [sRes.scalPressureTarget].'], 'rows');
    scalOmegaMax = 0;
    for s = sRes
        scalOmegaMax = max([scalOmegaMax; s.sStick.vecOmega/sqrt(s.scalMeanKn)], [], 'omitnan');
    end
    vecEdges = (0:options.scalBinWidth:(scalOmegaMax + options.scalBinWidth)).';
    vecEdgesFine = (0:options.scalBinWidth/2:(scalOmegaMax + options.scalBinWidth)).';
    sGroups = struct([]);
    for idxKey = 1:size(matKeys, 1)
        sIn = sRes([sRes.scalMu] == matKeys(idxKey, 1) & [sRes.scalPressureTarget] == matKeys(idxKey, 2));
        g = struct();
        g.scalMu = matKeys(idxKey, 1);
        g.scalPressureTarget = matKeys(idxKey, 2);
        g.scalNumSeeds = numel(sIn);
        g.scalNumRemoved = sum([sIn.scalNumRemoved]);
        g.scalDeltaZ = mean([sIn.scalZ]) - 3;
        g.scalDeltaZCut = mean([sIn.scalZCut]) - 3;
        for strVer = ["", "Cut"]
            vecW = []; scalModes = 0;
            for s = sIn
                sSpec = s.("sStick" + strVer);
                vecW = [vecW; sSpec.vecOmega/sqrt(s.("scalMeanKn" + strVer))]; %#ok<AGROW>
                scalModes = scalModes + sSpec.scalNumModes;
            end
            vecW = vecW(vecW > 0);
            [vecDOS, vecCenters] = calcDOS(vecW, scalModes, vecEdges);
            [vecDOSFine, vecCentersFine] = calcDOS(vecW, scalModes, vecEdgesFine);
            vecSorted = sort(vecW);
            vecEdgesEq = unique([0; vecSorted(options.scalModesPerBin:options.scalModesPerBin:end)]);
            [vecDOSEq, vecCentersEq] = calcDOS(vecW, scalModes, vecEdgesEq);
            g.("vecDOSStick" + strVer) = vecDOS;
            g.("scalT20" + strVer) = calcOmegaStar(vecCenters, vecDOS, options.scalOmegaRef);
            g.("scalT50" + strVer) = calcOmegaStar(vecCentersFine, vecDOSFine, options.scalOmegaRef, 0.5);
            g.("scalEq" + strVer) = calcOmegaStar(vecCentersEq, vecDOSEq, options.scalOmegaRef);
            g.("scalKnee" + strVer) = kneeFit(vecCentersFine, vecDOSFine, options.scalOmegaRef);
        end
        g.vecBinCenters = vecCenters;
        if isempty(sGroups)
            sGroups = g;
        else
            sGroups(end + 1) = g; %#ok<AGROW>
        end
    end
end

function scalOmegaK = kneeFit(vecCenters, vecDOS, scalOmegaRef)
% Least squares D = D_p w/(w + w_k) on 0 < w <= omega_ref: for each trial
% w_k (400 log-spaced values, 1e-3 to 2) D_p is linear and solved exactly;
% the best w_k is returned.
    vecBool = vecCenters > 0 & vecCenters <= scalOmegaRef;
    vecX = vecCenters(vecBool); vecY = vecDOS(vecBool);
    scalBest = Inf; scalOmegaK = NaN;
    vecWk = logspace(-3, 0.3, 400);
    for scalWk = vecWk
        vecG = vecX./(vecX + scalWk);
        scalDp = (vecG.'*vecY)/(vecG.'*vecG);
        scalErr = sum((vecY - scalDp*vecG).^2);
        if scalErr < scalBest
            scalBest = scalErr; scalOmegaK = scalWk;
        end
    end
    if scalOmegaK >= vecWk(end - 1)
        scalOmegaK = NaN;   % pinned at the upper end of the grid: no knee below omega_ref
    end
end

function [vecFit, scalR2] = linearFit(vecX, vecY)
    vecBool = ~isnan(vecX) & ~isnan(vecY);
    vecFit = [NaN NaN]; scalR2 = NaN;
    if nnz(vecBool) < 3
        return;
    end
    vecFit = polyfit(vecX(vecBool), vecY(vecBool), 1);
    vecRes = vecY(vecBool) - polyval(vecFit, vecX(vecBool));
    scalR2 = 1 - sum(vecRes.^2)/sum((vecY(vecBool) - mean(vecY(vecBool))).^2);
end

function sTheme = makeTheme(boolLightMode)
    if boolLightMode
        sTheme.vecBg = [1 1 1]; sTheme.vecFg = [0 0 0]; sTheme.vecGrid = [0.6 0.6 0.6];
        sTheme.vecKept = [0.25 0.25 0.25]; sTheme.vecCut = [0.85 0.33 0.10];
        sTheme.matSeries = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19; 0.49 0.18 0.56; 0.30 0.75 0.93; 0.64 0.08 0.18; 0.93 0.69 0.13];
    else
        sTheme.vecBg = [0.118 0.118 0.118]; sTheme.vecFg = [0.88 0.88 0.88]; sTheme.vecGrid = [0.5 0.5 0.5];
        sTheme.vecKept = [0.36 0.55 0.94]; sTheme.vecCut = [1.00 0.44 0.26];
        sTheme.matSeries = [0.36 0.55 0.94; 1.00 0.44 0.26; 0.40 0.73 0.42; 0.80 0.60 0.95; 0.95 0.85 0.35; 0.95 0.45 0.65; 0.55 0.85 0.85];
    end
    sTheme.boolLight = boolLightMode;
end

function hFig = newFigure(sTheme, vecPosition)
    hFig = figure('Color', sTheme.vecBg, 'Position', vecPosition);
    try
        if sTheme.boolLight
            theme(hFig, 'light');
        else
            theme(hFig, 'dark');
        end
    catch
        % theme() exists only in R2025a+
    end
end

function styleAxes(ax, sTheme)
    set(ax, 'Color', sTheme.vecBg, 'XColor', sTheme.vecFg, 'YColor', sTheme.vecFg, 'GridColor', sTheme.vecGrid, ...
        'TickLabelInterpreter', 'latex', 'FontSize', 12, 'Box', 'on');
    ax.Title.Color = sTheme.vecFg;
    ax.Title.Interpreter = 'latex';
    grid(ax, 'on');
end

function saveFigure(hFig, strPath, sTheme)
    exportgraphics(hFig, strPath, 'Resolution', 150, 'BackgroundColor', sTheme.vecBg);
    fprintf('[plotUnresolvedDOS] %s\n', strPath);
end

function str = latexPow(scalX)
% 1e-4 -> '10^{-4}', 3e-5 -> '3\times10^{-5}' (for labels inside $...$)
    scalExp = floor(log10(scalX) + 1e-12);
    scalMant = scalX/10^scalExp;
    if abs(scalMant - 1) < 1e-9
        str = sprintf('10^{%d}', scalExp);
    else
        str = sprintf('%g\\times10^{%d}', scalMant, scalExp);
    end
end
