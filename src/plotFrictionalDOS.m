function sSummary = plotFrictionalDOS(options)
% plotFrictionalDOS  Henkes et al. (2010) style figures from runFrictionalDOS output.
%
%   plotFrictionalDOS()
%   sSummary = plotFrictionalDOS("vecStrContactLaw", ["hooke", "hertz"], "strFigureDir", "figs")
%
%   Reads <strDataDir>/<law>/spectra.mat written by runFrictionalDOS and makes
%   the analogs of the figures of Henkes, van Hecke & van Saarloos, EPL 90,
%   14003 (2010), one panel per contact law:
%
%     fig1_z_nm.png        Fig. 1: packings in the (z, n_m) plane with the
%                          generalized isostaticity line n_m = z - 3
%     fig2_dos_mu.png      Fig. 3a/3b: DOS at the lowest pressure for every mu,
%                          fully mobilized contacts slipping (top) and sticking
%                          (bottom)
%     fig3_dos_pressure.png Fig. 3c: DOS for mu = 0.3 (slipping) at every pressure
%     fig4_all_slip.png    Fig. 3d: k_t -> 0 on every contact, compared with the
%                          frictionless DOS of the same backbone
%     fig5_omega_star.png  Fig. 4: crossover frequency omega* against the
%                          distance to isostaticity, z - (3 + n_m) (slipping,
%                          top) and z - 3 (sticking, bottom), with omega* from
%                          fixed-width bins (Henkes) and from equal-count bins
%     fig6_z_vs_p.png      z against p^(1/2) (Hooke) or p^(1/3) (Hertz): the
%                          excess contact number vanishes as a power of the
%                          pressure; the intercepts are z_c(mu) at p -> 0
%     fig7_vs_mu.png       z_c(mu) and the packing fraction against mu, with the
%                          mu -> 0 and mu -> infinity limits of Shundyak et al.
%     fig8_mobilization.png fraction of contacts with mobilization above m: the
%                          jump at m = 1 is the delta peak of fully mobilized
%                          contacts (Shundyak et al. Fig. 2b)
%     fig9_exact_checks.png two exact statements checked on EVERY packing:
%                          slipping only lowers modes (Courant-Fischer,
%                          omega_k^2(slip) <= omega_k^2(stick)), and k_t = 0 on
%                          every contact gives the frictionless spectrum plus
%                          N_i zero-frequency rotations
%
%   Frequencies are measured in units of the mean contact frequency,
%       omega_hat = omega / sqrt(<k_n>/m),
%   with <k_n> the mean normal contact stiffness of the packing. For Hookean
%   contacts <k_n> = K. For Hertz contacts <k_n> ~ p^(1/3), so omega_hat is
%   Henkes' omega/p^(1/6) up to a constant. Spectra of the seeds at the same
%   (mu, P) are pooled. omega* is where D(omega_hat)/D(omega_ref) first
%   reaches 0.2 (calcOmegaStar), with omega_ref = scalOmegaRef.
%
%   Options:
%       vecStrContactLaw  laws to plot (default ["hooke", "hertz"]; a law
%                         without spectra.mat is skipped with a warning)
%       strDataDir        runFrictionalDOS output (default <GranE>/data/frictionalDOS)
%       strFigureDir      where the PNGs go (default <strDataDir>/figures)
%       scalOmegaRef      reference frequency for omega* (default 1)
%       scalBinWidth      DOS bin width in omega_hat (default 0.08)
%       scalMuPressure    mu of the pressure series (default 0.3)
%       scalModesPerBin   modes per bin for the equal-count omega* estimate
%                         (default 40)
%
%   Returns sSummary(law) with the per-(mu, P) table: vecMu, vecPressure,
%   vecZ, vecNm, vecOmegaStarSlip, vecOmegaStarStick (and the fits).

    arguments
        options.vecStrContactLaw (1, :) string = ["hooke", "hertz"]
        options.strDataDir       (1, 1) string = ""
        options.strFigureDir     (1, 1) string = ""
        options.scalOmegaRef     (1, 1) double {mustBePositive} = 1
        options.scalBinWidth     (1, 1) double {mustBePositive} = 0.08
        options.scalMuPressure   (1, 1) double = 0.3
        options.scalModesPerBin  (1, 1) double {mustBePositive, mustBeInteger} = 40
    end

    strRepoRoot = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(strRepoRoot, 'src'));
    strDataDir = options.strDataDir;
    if strDataDir == ""
        strDataDir = fullfile(strRepoRoot, 'data', 'frictionalDOS');
    end
    strFigureDir = options.strFigureDir;
    if strFigureDir == ""
        strFigureDir = fullfile(strDataDir, 'figures');
    end
    if ~isfolder(strFigureDir)
        mkdir(strFigureDir);
    end

    % Load the spectra of every law that has been run.
    cellRuns = {};
    vecStrLaw = strings(1, 0);
    for strLaw = options.vecStrContactLaw
        strFile = fullfile(strDataDir, strLaw, 'spectra.mat');
        if ~isfile(strFile)
            warning('plotFrictionalDOS:NoData', 'No %s; run runFrictionalDOS("%s") first.', strFile, strLaw);
            continue;
        end
        sLoaded = load(strFile, 'sRuns');
        cellRuns{end + 1} = sLoaded.sRuns; %#ok<AGROW>
        vecStrLaw(end + 1) = strLaw; %#ok<AGROW>
    end
    if isempty(cellRuns)
        error('plotFrictionalDOS:NoData', 'No spectra found in %s.', strDataDir);
    end
    scalNumLaws = numel(cellRuns);
    cellLawTitle = containers.Map({'hooke', 'hertz'}, {'Hookean + Cundall-Strack', 'Hertz-Mindlin'});

    % Pool the seeds of every (mu, P) group.
    cellGroups = cell(1, scalNumLaws);
    for idxLaw = 1:scalNumLaws
        cellGroups{idxLaw} = groupRuns(cellRuns{idxLaw}, options.scalBinWidth, options.scalOmegaRef, options.scalModesPerBin);
    end

    sSummary = struct([]);
    hFig1 = newFigure([0 0 560*scalNumLaws 520]);
    hTile1 = tiledlayout(hFig1, 1, scalNumLaws, 'TileSpacing', 'compact');
    hFig2 = newFigure([0 0 560*scalNumLaws 820]);
    hTile2 = tiledlayout(hFig2, 2, scalNumLaws, 'TileSpacing', 'compact');
    hFig3 = newFigure([0 0 560*scalNumLaws 450]);
    hTile3 = tiledlayout(hFig3, 1, scalNumLaws, 'TileSpacing', 'compact');
    hFig4 = newFigure([0 0 560*scalNumLaws 450]);
    hTile4 = tiledlayout(hFig4, 1, scalNumLaws, 'TileSpacing', 'compact');
    hFig5 = newFigure([0 0 600*scalNumLaws 880]);
    hTile5 = tiledlayout(hFig5, 2, scalNumLaws, 'TileSpacing', 'compact');
    hFig6 = newFigure([0 0 560*scalNumLaws 480]);
    hTile6 = tiledlayout(hFig6, 1, scalNumLaws, 'TileSpacing', 'compact');
    hFig7 = newFigure([0 0 1100 460]);
    hTile7 = tiledlayout(hFig7, 1, 2, 'TileSpacing', 'compact');
    hAxZc = nexttile(hTile7, 1); hold(hAxZc, 'on'); box(hAxZc, 'on'); grid(hAxZc, 'on');
    hAxPhi = nexttile(hTile7, 2); hold(hAxPhi, 'on'); box(hAxPhi, 'on'); grid(hAxPhi, 'on');
    hFig8 = newFigure([0 0 560*scalNumLaws 450]);
    hTile8 = tiledlayout(hFig8, 1, scalNumLaws, 'TileSpacing', 'compact');
    hFig9 = newFigure([0 0 1100 460]);
    hTile9 = tiledlayout(hFig9, 1, 2, 'TileSpacing', 'compact');
    hAxCF = nexttile(hTile9, 1); hold(hAxCF, 'on'); box(hAxCF, 'on'); grid(hAxCF, 'on');
    hAxAS = nexttile(hTile9, 2); hold(hAxAS, 'on'); box(hAxAS, 'on'); grid(hAxAS, 'on');
    matColorLaw = [0.85 0.33 0.10; 0.00 0.45 0.74];

    for idxLaw = 1:scalNumLaws
        sGroups = cellGroups{idxLaw};
        strTitle = cellLawTitle(char(vecStrLaw(idxLaw)));
        vecMuAll = unique([sGroups.scalMu]);
        vecPAll = unique([sGroups.scalPressureTarget]);
        matColorMu = parula(numel(vecMuAll) + 1);
        cellMarkerP = {'o', 's', '^', 'd', 'v', 'p'};

        % ---------- Fig. 1: (z, n_m) plane ----------
        hAx = nexttile(hTile1); hold(hAx, 'on'); box(hAx, 'on'); grid(hAx, 'on');
        plot(hAx, [3 4], [0 1], 'k-', 'LineWidth', 2, 'DisplayName', 'n_m = z - 3 (generalized isostaticity)');
        for sRun = cellRuns{idxLaw}
            idxMu = find(vecMuAll == sRun.scalMu);
            idxP = find(vecPAll == sRun.scalPressureTarget);
            plot(hAx, sRun.scalZ, sRun.scalNm, cellMarkerP{1 + mod(idxP - 1, numel(cellMarkerP))}, ...
                'MarkerFaceColor', matColorMu(idxMu, :), 'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'HandleVisibility', 'off');
        end
        addLegendEntries(hAx, vecMuAll, vecPAll, matColorMu, cellMarkerP);
        xlabel(hAx, 'z'); ylabel(hAx, 'n_m');
        xlim(hAx, [2.9 4.3]); ylim(hAx, [-0.02 1]);
        title(hAx, {strTitle, 'stable packings lie right of the line'});
        legend(hAx, 'Location', 'northwest', 'FontSize', 8);

        % ---------- Fig. 2: DOS at the lowest pressure, every mu ----------
        scalPLow = lowestCommonPressure(sGroups, vecMuAll);
        for idxRow = 1:2
            hAx = nexttile(hTile2, (idxRow - 1)*scalNumLaws + idxLaw); hold(hAx, 'on'); box(hAx, 'on'); grid(hAx, 'on');
            for idxMu = 1:numel(vecMuAll)
                sG = sGroups([sGroups.scalMu] == vecMuAll(idxMu) & [sGroups.scalPressureTarget] == scalPLow);
                if isempty(sG), continue; end
                if idxRow == 1
                    vecD = sG.vecDOSSlip;
                else
                    vecD = sG.vecDOSStick;
                end
                plot(hAx, sG.vecBinCenters, vecD, '-', 'Color', matColorMu(idxMu, :), 'LineWidth', 1.6, ...
                    'DisplayName', sprintf('\\mu = %g', vecMuAll(idxMu)));
            end
            xline(hAx, options.scalOmegaRef, 'k:', 'HandleVisibility', 'off');
            xlabel(hAx, '\omega / (<k_n>/m)^{1/2}'); ylabel(hAx, 'D(\omega)');
            if idxRow == 1
                title(hAx, {strTitle, sprintf('P = %g, m = 1 contacts SLIP (Henkes Fig. 3a)', scalPLow)});
            else
                title(hAx, sprintf('P = %g, all contacts STICK (Henkes Fig. 3b)', scalPLow));
            end
            legend(hAx, 'Location', 'northeast', 'FontSize', 8);
        end

        % ---------- Fig. 3: pressure series at mu = scalMuPressure ----------
        hAx = nexttile(hTile3); hold(hAx, 'on'); box(hAx, 'on'); grid(hAx, 'on');
        sGmu = sGroups([sGroups.scalMu] == options.scalMuPressure);
        [~, vecOrder] = sort([sGmu.scalPressureTarget], 'descend');
        matColorP = cool(max(numel(sGmu), 2));
        for idxG = 1:numel(vecOrder)
            sG = sGmu(vecOrder(idxG));
            plot(hAx, sG.vecBinCenters, sG.vecDOSSlip, '-', 'Color', matColorP(idxG, :), 'LineWidth', 1.6, ...
                'DisplayName', sprintf('P = %g (z - 3 - n_m = %.2f)', sG.scalPressureTarget, sG.scalDeltaZSlip));
        end
        xlabel(hAx, '\omega / (<k_n>/m)^{1/2}'); ylabel(hAx, 'D(\omega)');
        xlim(hAx, [0 1.2]);
        title(hAx, {strTitle, sprintf('\\mu = %g, slipping, low-frequency end (Henkes Fig. 3c)', options.scalMuPressure)});
        legend(hAx, 'Location', 'northeast', 'FontSize', 8);

        % ---------- Fig. 4: all contacts slipping vs frictionless ----------
        hAx = nexttile(hTile4); hold(hAx, 'on'); box(hAx, 'on'); grid(hAx, 'on');
        sG = sGroups([sGroups.scalMu] == min(vecMuAll) & [sGroups.scalPressureTarget] == scalPLow);
        plot(hAx, sG.vecBinCenters, (2/3)*sG.vecDOSFrictionless, 'k-', 'LineWidth', 3.5, ...
            'DisplayName', 'frictionless backbone \times 2/3 (2 of 3 DOF)');
        plot(hAx, sG.vecBinCenters, sG.vecDOSAllSlip, 'r--', 'LineWidth', 1.8, ...
            'DisplayName', sprintf('frictional, k_t = 0 on every contact (zero modes: %.3f of all)', sG.scalZeroFractionAllSlip));
        plot(hAx, sG.vecBinCenters, sG.vecDOSStick, '-', 'Color', [0.3 0.6 1], 'LineWidth', 1.2, ...
            'DisplayName', 'frictional, sticking');
        xlabel(hAx, '\omega / (<k_n>/m)^{1/2}'); ylabel(hAx, 'D(\omega)');
        title(hAx, {strTitle, sprintf('\\mu = %g, P = %g: k_t \\rightarrow 0 (Henkes Fig. 3d)', min(vecMuAll), scalPLow)});
        legend(hAx, 'Location', 'northeast', 'FontSize', 8);

        % ---------- Fig. 5: omega* vs distance to isostaticity ----------
        % One point per PACKING: seeds at the same (mu, P) can have very
        % different n_m (e.g. 0.05 and 0.86), so their Delta z must not be averaged.
        sPerRun = perRunOmegaStar(cellRuns{idxLaw}, options.scalBinWidth, options.scalOmegaRef, options.scalModesPerBin);
        vecMuRun = [cellRuns{idxLaw}.scalMu];
        vecX1 = sPerRun.vecDeltaZSlip; vecY1 = sPerRun.vecOmegaStarSlip; vecY1b = sPerRun.vecOmegaStarSlipEq;
        vecX2 = sPerRun.vecDeltaZStick; vecY2 = sPerRun.vecOmegaStarStick; vecY2b = sPerRun.vecOmegaStarStickEq;
        cellRows = {{vecX1, vecY1, vecY1b, 'slipping (Henkes): \Deltaz = z - (3 + n_m)'}, ...
                    {vecX2, vecY2, vecY2b, 'sticking: \Deltaz = z - 3'}};
        for idxRow = 1:2
            vecX = cellRows{idxRow}{1}; vecYa = cellRows{idxRow}{2}; vecYb = cellRows{idxRow}{3};
            hAx = nexttile(hTile5, (idxRow - 1)*scalNumLaws + idxLaw); hold(hAx, 'on'); box(hAx, 'on'); grid(hAx, 'on');
            for idxG = 1:numel(vecX)
                idxMu = find(vecMuAll == vecMuRun(idxG));
                plot(hAx, vecX(idxG), vecYa(idxG), 'o', 'MarkerFaceColor', matColorMu(idxMu, :), 'MarkerEdgeColor', 'k', ...
                    'MarkerSize', 8, 'HandleVisibility', 'off');
                if isnan(vecYa(idxG))
                    % D/D(omega_ref) >= 0.2 already in the first bin: omega* is
                    % below the first bin, shown as an upper limit at its edge.
                    plot(hAx, vecX(idxG), options.scalBinWidth, 'v', 'MarkerFaceColor', matColorMu(idxMu, :), ...
                        'MarkerEdgeColor', 'k', 'MarkerSize', 7, 'HandleVisibility', 'off');
                end
                plot(hAx, vecX(idxG), vecYb(idxG), 'd', 'MarkerEdgeColor', matColorMu(idxMu, :), 'LineWidth', 1.5, ...
                    'MarkerSize', 8, 'HandleVisibility', 'off');
            end
            [vecFitA, scalR2A] = linearFit(vecX, vecYa);
            [vecFitB, scalR2B] = linearFit(vecX, vecYb);
            vecXLine = linspace(0, max(vecX, [], 'omitnan')*1.05, 50);
            plot(hAx, vecXLine, polyval(vecFitA, vecXLine), 'k-', 'DisplayName', ...
                sprintf('fixed bins: \\omega* = %.2f \\Deltaz %+.3f (R^2 = %.2f)', vecFitA(1), vecFitA(2), scalR2A));
            plot(hAx, vecXLine, polyval(vecFitB, vecXLine), 'k--', 'DisplayName', ...
                sprintf('equal-count bins: \\omega* = %.2f \\Deltaz %+.3f (R^2 = %.2f)', vecFitB(1), vecFitB(2), scalR2B));
            for idxMu = 1:numel(vecMuAll)
                plot(hAx, NaN, NaN, 'o', 'MarkerFaceColor', matColorMu(idxMu, :), 'MarkerEdgeColor', 'k', ...
                    'DisplayName', sprintf('\\mu = %g', vecMuAll(idxMu)));
            end
            plot(hAx, NaN, NaN, 'ko', 'MarkerFaceColor', [0.6 0.6 0.6], 'DisplayName', 'o: fixed-width bins (Henkes)');
            plot(hAx, NaN, NaN, 'kd', 'DisplayName', sprintf('\\diamondsuit: %d modes per bin', modesPerBin(options)));
            plot(hAx, NaN, NaN, 'kv', 'MarkerFaceColor', [0.6 0.6 0.6], 'DisplayName', '\nabla: \omega* below the first bin (upper limit)');
            xlabel(hAx, cellRows{idxRow}{4}); ylabel(hAx, '\omega* / (<k_n>/m)^{1/2}');
            xlim(hAx, [0, max(vecXLine)]); ylim(hAx, [0, max([vecYa vecYb options.scalBinWidth], [], 'omitnan')*1.2]);
            if idxRow == 1
                title(hAx, {strTitle, 'omega* vs distance to isostaticity, one point per packing (Henkes Fig. 4)'});
            end
            legend(hAx, 'Location', 'northwest', 'FontSize', 8);
            if idxRow == 1
                vecFitSlip = vecFitA; scalR2Slip = scalR2A;
            else
                vecFitStick = vecFitA; scalR2Stick = scalR2A;
            end
        end

        % ---------- Fig. 6: z vs p^beta, intercept z_c(mu) ----------
        if vecStrLaw(idxLaw) == "hertz"
            scalBeta = 1/3; strBeta = '1/3';
        else
            scalBeta = 1/2; strBeta = '1/2';
        end
        hAx = nexttile(hTile6); hold(hAx, 'on'); box(hAx, 'on'); grid(hAx, 'on');
        vecZc = NaN(size(vecMuAll));
        for idxMu = 1:numel(vecMuAll)
            sRunsMu = cellRuns{idxLaw}([cellRuns{idxLaw}.scalMu] == vecMuAll(idxMu));
            vecXp = [sRunsMu.scalPressureVirial].^scalBeta;
            vecYz = [sRunsMu.scalZ];
            plot(hAx, vecXp, vecYz, 'o', 'MarkerFaceColor', matColorMu(idxMu, :), 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
            vecFitZ = polyfit(vecXp, vecYz, 1);
            vecZc(idxMu) = vecFitZ(2);
            vecXLine = linspace(0, max(vecXp), 20);
            plot(hAx, vecXLine, polyval(vecFitZ, vecXLine), '-', 'Color', matColorMu(idxMu, :), ...
                'DisplayName', sprintf('\\mu = %g: z_c = %.2f', vecMuAll(idxMu), vecZc(idxMu)));
        end
        yline(hAx, 4, 'k--', 'z = 4 (frictionless isostatic)', 'HandleVisibility', 'off');
        yline(hAx, 3, 'k--', 'z = 3 (d + 1)', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
        xlabel(hAx, sprintf('p^{%s}  (virial pressure)', strBeta)); ylabel(hAx, 'z (backbone)');
        title(hAx, {strTitle, sprintf('z - z_c \\propto p^{%s}: lines extrapolate to z_c(\\mu)', strBeta)});
        legend(hAx, 'Location', 'southeast', 'FontSize', 8);

        % ---------- Fig. 7: z_c and phi vs mu (both laws on one axis) ----------
        vecPhiLow = arrayfun(@(scalMu) mean([sGroups([sGroups.scalMu] == scalMu & [sGroups.scalPressureTarget] == scalPLow).scalPhi]), vecMuAll);
        plot(hAxZc, vecMuAll, vecZc, 'o-', 'Color', matColorLaw(idxLaw, :), 'MarkerFaceColor', matColorLaw(idxLaw, :), ...
            'DisplayName', [strTitle ': z_c (p \rightarrow 0 fit)']);
        plot(hAxPhi, vecMuAll, vecPhiLow, 'o-', 'Color', matColorLaw(idxLaw, :), 'MarkerFaceColor', matColorLaw(idxLaw, :), ...
            'DisplayName', sprintf('%s: P = %g', strTitle, scalPLow));

        % ---------- Fig. 8: mobilization distribution ----------
        hAx = nexttile(hTile8); hold(hAx, 'on'); box(hAx, 'on'); grid(hAx, 'on');
        for idxMu = 1:numel(vecMuAll)
            sRunsMu = cellRuns{idxLaw}([cellRuns{idxLaw}.scalMu] == vecMuAll(idxMu) & [cellRuns{idxLaw}.scalPressureTarget] == scalPLow);
            vecM = sort(vertcat(sRunsMu.vecMobilization));
            vecFracAbove = 1 - (0:numel(vecM) - 1).'/numel(vecM);
            stairs(hAx, [0; vecM], [1; vecFracAbove], '-', 'Color', matColorMu(idxMu, :), 'LineWidth', 1.6, ...
                'DisplayName', sprintf('\\mu = %g (fraction at m = 1: %.2f)', vecMuAll(idxMu), mean(vecM >= 1 - 1e-6)));
        end
        xlim(hAx, [0 1.02]); ylim(hAx, [0 1]);
        xlabel(hAx, 'mobilization m = |f_t| / (\mu f_n)'); ylabel(hAx, 'fraction of contacts with mobilization \geq m');
        title(hAx, {strTitle, sprintf('P = %g: jump at m = 1 = fully mobilized contacts', scalPLow)});
        legend(hAx, 'Location', 'southwest', 'FontSize', 8);

        % ---------- Fig. 9: exact checks on every packing ----------
        sChecks = exactChecks(cellRuns{idxLaw});
        plot(hAxCF, sChecks.vecDeltaZ, max(sChecks.vecCourantFischer, 1e-18), 'o', 'Color', matColorLaw(idxLaw, :), ...
            'MarkerFaceColor', matColorLaw(idxLaw, :), 'DisplayName', sprintf('%s (%d packings)', strTitle, numel(sChecks.vecDeltaZ)));
        plot(hAxAS, sChecks.vecDeltaZ, max(sChecks.vecAllSlipError, 1e-18), 'o', 'Color', matColorLaw(idxLaw, :), ...
            'MarkerFaceColor', matColorLaw(idxLaw, :), 'DisplayName', sprintf('%s: zero-mode count = N_i + frictionless in %d of %d', ...
            strTitle, sum(sChecks.vecBoolZeroCount), numel(sChecks.vecBoolZeroCount)));
        fprintf('   exact checks: max_k [omega_k^2(slip) - omega_k^2(stick)]/max omega^2 <= %.1e over all packings; all-slip vs frictionless max rel. diff %.1e; zero-mode count N_i + frictionless in %d of %d\n', ...
            max(sChecks.vecCourantFischer), max(sChecks.vecAllSlipError), sum(sChecks.vecBoolZeroCount), numel(sChecks.vecBoolZeroCount));

        % ---------- summary table ----------
        sLaw = struct();
        sLaw.strContactLaw = char(vecStrLaw(idxLaw));
        sLaw.vecMu = [sGroups.scalMu];
        sLaw.vecPressure = [sGroups.scalPressureTarget];
        sLaw.vecNumSeeds = [sGroups.scalNumSeeds];
        sLaw.vecZ = [sGroups.scalZ];
        sLaw.vecNm = [sGroups.scalNm];
        sLaw.vecOmegaStarSlip = [sGroups.scalOmegaStarSlip];
        sLaw.vecOmegaStarStick = [sGroups.scalOmegaStarStick];
        sLaw.sPerRun = sPerRun;
        sLaw.vecFitSlip = vecFitSlip; sLaw.scalR2Slip = scalR2Slip;
        sLaw.vecFitStick = vecFitStick; sLaw.scalR2Stick = scalR2Stick;
        sLaw.vecOmegaStarSlipEq = [sGroups.scalOmegaStarSlipEq];
        sLaw.vecOmegaStarStickEq = [sGroups.scalOmegaStarStickEq];
        sLaw.vecMuZc = vecMuAll;
        sLaw.vecZc = vecZc;
        sLaw.vecPhiLowP = vecPhiLow;
        if isempty(sSummary)
            sSummary = sLaw;
        else
            sSummary(end + 1) = sLaw; %#ok<AGROW>
        end
        fprintf('\n[plotFrictionalDOS] %s\n   mu        P   seeds     z     n_m   z-3-n_m   phi    w*slip  w*stick (equal-count: slip stick)\n', strTitle);
        for idxG = 1:numel(sGroups)
            fprintf('%6g %8g %5d %7.3f %7.3f %8.3f %7.4f %7.3f %7.3f   %7.3f %7.3f\n', sGroups(idxG).scalMu, sGroups(idxG).scalPressureTarget, ...
                sGroups(idxG).scalNumSeeds, sGroups(idxG).scalZ, sGroups(idxG).scalNm, sGroups(idxG).scalDeltaZSlip, ...
                sGroups(idxG).scalPhi, sGroups(idxG).scalOmegaStarSlip, sGroups(idxG).scalOmegaStarStick, ...
                sGroups(idxG).scalOmegaStarSlipEq, sGroups(idxG).scalOmegaStarStickEq);
        end
        fprintf('   z_c(mu) from z vs p^beta: %s at mu = %s\n', mat2str(vecZc, 4), mat2str(vecMuAll));
        fprintf('   per-packing fits (fixed bins), %d of %d slip and %d of %d stick omega* resolved:\n', ...
            sum(isfinite(vecY1)), numel(vecY1), sum(isfinite(vecY2)), numel(vecY2));
        fprintf('   fit slip:  omega* = %.3f dz + %.3f (R^2 = %.3f)\n   fit stick: omega* = %.3f dz + %.3f (R^2 = %.3f)\n', ...
            vecFitSlip(1), vecFitSlip(2), scalR2Slip, vecFitStick(1), vecFitStick(2), scalR2Stick);
    end

    title(hTile1, 'Packings in the (z, n_m) plane (Henkes et al. 2010, Fig. 1)');
    title(hTile2, 'Density of states at the lowest pressure (Henkes et al. 2010, Fig. 3a/3b)');
    exportgraphics(hFig1, fullfile(strFigureDir, 'fig1_z_nm.png'), 'Resolution', 130);
    exportgraphics(hFig2, fullfile(strFigureDir, 'fig2_dos_mu.png'), 'Resolution', 130);
    exportgraphics(hFig3, fullfile(strFigureDir, 'fig3_dos_pressure.png'), 'Resolution', 130);
    exportgraphics(hFig4, fullfile(strFigureDir, 'fig4_all_slip.png'), 'Resolution', 130);
    exportgraphics(hFig5, fullfile(strFigureDir, 'fig5_omega_star.png'), 'Resolution', 130);
    yline(hAxZc, 3.98, 'b--', '\mu \rightarrow 0 (Shundyak 2007)', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
    yline(hAxZc, 3.00, 'r--', '\mu \rightarrow \infty (Shundyak 2007)', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
    yline(hAxPhi, 0.8395, 'b--', '\mu \rightarrow 0 (Shundyak 2007, p \rightarrow 0)', 'HandleVisibility', 'off');
    yline(hAxPhi, 0.758, 'r--', '\mu \rightarrow \infty (Shundyak 2007, p \rightarrow 0)', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
    set([hAxZc hAxPhi], 'XScale', 'log');
    xlabel(hAxZc, '\mu'); ylabel(hAxZc, 'z_c');
    xlabel(hAxPhi, '\mu'); ylabel(hAxPhi, '\phi (rattlers included)');
    ylim(hAxZc, [2.9 4.1]); ylim(hAxPhi, [0.74 0.86]);
    title(hAxZc, 'Contact number at jamming vs friction');
    title(hAxPhi, 'Packing fraction vs friction (lowest pressure)');
    legend(hAxZc, 'Location', 'east', 'FontSize', 8);
    legend(hAxPhi, 'Location', 'northeast', 'FontSize', 8);
    exportgraphics(hFig6, fullfile(strFigureDir, 'fig6_z_vs_p.png'), 'Resolution', 130);
    exportgraphics(hFig7, fullfile(strFigureDir, 'fig7_vs_mu.png'), 'Resolution', 130);
    exportgraphics(hFig8, fullfile(strFigureDir, 'fig8_mobilization.png'), 'Resolution', 130);
    set([hAxCF hAxAS], 'YScale', 'log');
    yline(hAxCF, eps, 'k--', 'machine \epsilon', 'HandleVisibility', 'off');
    xlabel(hAxCF, 'z - 3 (packing)'); ylabel(hAxCF, 'max_k [\omega_k^2(slip) - \omega_k^2(stick)] / max \omega^2');
    title(hAxCF, {'Slipping never raises a mode (Courant-Fischer)', 'positive values would violate it; floor 10^{-18} = exactly 0'});
    xlabel(hAxAS, 'z - 3 (packing)'); ylabel(hAxAS, 'max |\omega^2_{all slip} - \omega^2_{frictionless}| / max \omega^2');
    title(hAxAS, {'k_t = 0 everywhere: frictionless spectrum + N_i rotations', '(N_i modes closest to 0 removed before comparing)', ''});
    legend(hAxCF, 'Location', 'northeast', 'FontSize', 8);
    legend(hAxAS, 'Location', 'northeast', 'FontSize', 8);
    exportgraphics(hFig9, fullfile(strFigureDir, 'fig9_exact_checks.png'), 'Resolution', 130);
    fprintf('[plotFrictionalDOS] figures saved to %s\n', strFigureDir);
end

% =========================================================================
function sGroups = groupRuns(sRuns, scalBinWidth, scalOmegaRef, scalModesPerBin)
% Pool the seeds of every (mu, P) and compute the DOS of each stiffness version.
    matKeys = unique([[sRuns.scalMu].', [sRuns.scalPressureTarget].'], 'rows');
    % Common bins in omega_hat for all groups of this law.
    scalOmegaMax = 0;
    for sRun = sRuns
        scalOmegaMax = max([scalOmegaMax; sRun.sStick.vecOmega/sqrt(sRun.scalMeanKn)], [], 'omitnan');
    end
    vecEdges = 0:scalBinWidth:(scalOmegaMax + scalBinWidth);
    sGroups = struct([]);
    for idxKey = 1:size(matKeys, 1)
        vecBoolIn = [sRuns.scalMu] == matKeys(idxKey, 1) & [sRuns.scalPressureTarget] == matKeys(idxKey, 2);
        sIn = sRuns(vecBoolIn);
        sG = struct();
        sG.scalMu = matKeys(idxKey, 1);
        sG.scalPressureTarget = matKeys(idxKey, 2);
        sG.scalNumSeeds = numel(sIn);
        sG.scalZ = mean([sIn.scalZ]);
        sG.scalNm = mean([sIn.scalNm]);
        sG.scalDeltaZSlip = mean([sIn.scalZ] - 3 - [sIn.scalNm]);
        sG.scalDeltaZStick = mean([sIn.scalZ] - 3);
        sG.scalPressureVirial = mean([sIn.scalPressureVirial]);
        sG.scalPhi = mean([sIn.scalPhi]);
        cellVersions = {'sStick', 'sSlip', 'sAllSlip', 'sFrictionless'};
        cellNames = {'Stick', 'Slip', 'AllSlip', 'Frictionless'};
        for idxV = 1:numel(cellVersions)
            vecOmegaPooled = [];
            scalModes = 0;
            scalZero = 0;
            for sRun = sIn
                sSpec = sRun.(cellVersions{idxV});
                vecOmegaPooled = [vecOmegaPooled; sSpec.vecOmega/sqrt(sRun.scalMeanKn)]; %#ok<AGROW>
                scalModes = scalModes + sSpec.scalNumModes;
                scalZero = scalZero + sSpec.scalNumZero;
            end
            vecOmegaPooled = vecOmegaPooled(vecOmegaPooled > 0);   % zero modes are not in a bin
            [vecDOS, vecCenters] = calcDOS(vecOmegaPooled, scalModes, vecEdges);
            sG.(['vecDOS' cellNames{idxV}]) = vecDOS;
            % Equal-count bins (scalModesPerBin modes each): a bin-width-free
            % estimate of omega*, finest where the modes are densest.
            vecSorted = sort(vecOmegaPooled);
            vecEdgesEq = unique([0; vecSorted(scalModesPerBin:scalModesPerBin:end)]);
            [vecDOSEq, vecCentersEq] = calcDOS(vecOmegaPooled, scalModes, vecEdgesEq);
            sG.(['vecDOSEq' cellNames{idxV}]) = vecDOSEq;
            sG.(['vecCentersEq' cellNames{idxV}]) = vecCentersEq;
            sG.(['scalZeroFraction' cellNames{idxV}]) = scalZero/scalModes;
        end
        sG.vecBinCenters = vecCenters;
        sG.scalOmegaStarSlip = calcOmegaStar(vecCenters, sG.vecDOSSlip, scalOmegaRef);
        sG.scalOmegaStarStick = calcOmegaStar(vecCenters, sG.vecDOSStick, scalOmegaRef);
        sG.scalOmegaStarSlipEq = calcOmegaStar(sG.vecCentersEqSlip, sG.vecDOSEqSlip, scalOmegaRef);
        sG.scalOmegaStarStickEq = calcOmegaStar(sG.vecCentersEqStick, sG.vecDOSEqStick, scalOmegaRef);
        if isempty(sGroups)
            sGroups = sG;
        else
            sGroups(end + 1) = sG; %#ok<AGROW>
        end
    end
end

function scalPLow = lowestCommonPressure(sGroups, vecMuAll)
% Lowest target pressure that every mu was run at.
    vecP = unique([sGroups.scalPressureTarget]);
    scalPLow = max(vecP);
    for scalP = sort(vecP, 'descend')
        if all(arrayfun(@(scalMu) any([sGroups.scalMu] == scalMu & [sGroups.scalPressureTarget] == scalP), vecMuAll))
            scalPLow = scalP;
        end
    end
end

function [vecFit, scalR2] = linearFit(vecX, vecY)
% Least-squares line through the finite points, and its R^2.
    vecBool = isfinite(vecX) & isfinite(vecY);
    vecFit = [NaN NaN];
    scalR2 = NaN;
    if sum(vecBool) < 2
        return;
    end
    vecFit = polyfit(vecX(vecBool), vecY(vecBool), 1);
    vecRes = vecY(vecBool) - polyval(vecFit, vecX(vecBool));
    scalR2 = 1 - sum(vecRes.^2)/sum((vecY(vecBool) - mean(vecY(vecBool))).^2);
end

function addLegendEntries(hAx, vecMuAll, vecPAll, matColorMu, cellMarkerP)
    for idxMu = 1:numel(vecMuAll)
        plot(hAx, NaN, NaN, 'o', 'MarkerFaceColor', matColorMu(idxMu, :), 'MarkerEdgeColor', 'k', ...
            'DisplayName', sprintf('\\mu = %g', vecMuAll(idxMu)));
    end
    for idxP = 1:numel(vecPAll)
        plot(hAx, NaN, NaN, cellMarkerP{1 + mod(idxP - 1, numel(cellMarkerP))}, 'MarkerEdgeColor', 'k', ...
            'DisplayName', sprintf('P = %g', vecPAll(idxP)));
    end
end

function hFig = newFigure(vecPosition)
    hFig = figure('Position', vecPosition, 'Visible', 'off');
    try
        theme(hFig, 'light');
    catch
        % theme() is not available before R2025a
    end
end

function scalModes = modesPerBin(sOptions)
    scalModes = sOptions.scalModesPerBin;
end

function sChecks = exactChecks(sRuns)
% Per packing: (1) Courant-Fischer, omega_k^2(slip) <= omega_k^2(stick) for
% every sorted k, as the largest positive excess relative to max omega^2;
% (2) the all-slip spectrum minus its N_i modes closest to zero equals the
% frictionless spectrum of the same backbone (relative max difference), and
% the zero-mode counts add up.
    scalNumRuns = numel(sRuns);
    sChecks = struct('vecDeltaZ', zeros(1, scalNumRuns), 'vecCourantFischer', zeros(1, scalNumRuns), ...
        'vecAllSlipError', zeros(1, scalNumRuns), 'vecBoolZeroCount', false(1, scalNumRuns));
    for idxRun = 1:scalNumRuns
        sRun = sRuns(idxRun);
        vecStick = sRun.sStick.vecOmegaSq;
        vecSlip = sRun.sSlip.vecOmegaSq;
        scalMax = max(vecStick);
        sChecks.vecDeltaZ(idxRun) = sRun.scalZ - 3;
        sChecks.vecCourantFischer(idxRun) = max(0, max(vecSlip - vecStick))/scalMax;
        vecAll = sRun.sAllSlip.vecOmegaSq;
        [~, vecOrder] = sort(abs(vecAll));
        vecAllNonRot = sort(vecAll(vecOrder(sRun.scalNi + 1:end)));
        vecFrictionless = sRun.sFrictionless.vecOmegaSq;
        sChecks.vecAllSlipError(idxRun) = max(abs(vecAllNonRot - vecFrictionless))/max(abs(vecFrictionless));
        sChecks.vecBoolZeroCount(idxRun) = (sRun.sAllSlip.scalNumZero == sRun.scalNi + sRun.sFrictionless.scalNumZero);
    end
end

function sPerRun = perRunOmegaStar(sRuns, scalBinWidth, scalOmegaRef, scalModesPerBin)
% omega* of every packing on its own (fixed-width and equal-count bins).
    scalOmegaMax = 0;
    for sRun = sRuns
        scalOmegaMax = max([scalOmegaMax; sRun.sStick.vecOmega/sqrt(sRun.scalMeanKn)], [], 'omitnan');
    end
    vecEdges = 0:scalBinWidth:(scalOmegaMax + scalBinWidth);
    scalNumRuns = numel(sRuns);
    sPerRun = struct('vecDeltaZSlip', [sRuns.scalZ] - 3 - [sRuns.scalNm], 'vecDeltaZStick', [sRuns.scalZ] - 3, ...
        'vecOmegaStarSlip', NaN(1, scalNumRuns), 'vecOmegaStarStick', NaN(1, scalNumRuns), ...
        'vecOmegaStarSlipEq', NaN(1, scalNumRuns), 'vecOmegaStarStickEq', NaN(1, scalNumRuns));
    cellVersions = {'sSlip', 'sStick'};
    cellNames = {'Slip', 'Stick'};
    for idxRun = 1:scalNumRuns
        sRun = sRuns(idxRun);
        for idxV = 1:2
            sSpec = sRun.(cellVersions{idxV});
            vecOmega = sSpec.vecOmega/sqrt(sRun.scalMeanKn);
            vecOmega = vecOmega(vecOmega > 0);
            [vecDOS, vecCenters] = calcDOS(vecOmega, sSpec.scalNumModes, vecEdges);
            sPerRun.(['vecOmegaStar' cellNames{idxV}])(idxRun) = calcOmegaStar(vecCenters, vecDOS, scalOmegaRef);
            vecSorted = sort(vecOmega);
            vecEdgesEq = unique([0; vecSorted(scalModesPerBin:scalModesPerBin:end)]);
            [vecDOSEq, vecCentersEq] = calcDOS(vecOmega, sSpec.scalNumModes, vecEdgesEq);
            sPerRun.(['vecOmegaStar' cellNames{idxV} 'Eq'])(idxRun) = calcOmegaStar(vecCentersEq, vecDOSEq, scalOmegaRef);
        end
    end
end
