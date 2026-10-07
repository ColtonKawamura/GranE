function plotAttenuation(combinedFile, options)
%PLOTATTENUATION Plot attenuation coefficient alpha vs angular frequency omega
% from the combined simMD data produced by processData.m.

    arguments
        combinedFile (1,1) string
        options (1,1) struct
    end
    if ~isfield(options, 'pressure'), options.pressure = []; end

    D = load(combinedFile);
    data = D;

    atten = data.attenuation_x(:);
    atten(atten <= 0) = 1e-6;
    omega = data.omega(:);
    P     = sort(unique(data.pressure_actual(:)));
    if ~isempty(options.pressure) && isvector(options.pressure)
        P = options.pressure(:)';
    end

    fig = figure('Name', 'Attenuation vs Angular Frequency', 'Color', 'w');
    ax  = axes('Parent', fig);
    hold(ax, 'on');
    grid(ax, 'on'); box(ax, 'on');
    set(ax, 'XScale', 'log', 'YScale', 'log');

    cmap = lines(numel(P));
    for p = 1:numel(P)
        pv = P(p);
        mask = data.pressure_actual == pv;
        if any(mask)
            plot(ax, omega(mask), atten(mask), '-o', 'MarkerSize', 6, ...
                'LineWidth', 1.2, 'Color', cmap(p,:), ...
                'DisplayName', sprintf('P = %.4f', pv));
        end
    end
    hold(ax, 'off');

    xlabel(ax, '$\hat{\omega}$', 'FontSize', 16, 'Interpreter', 'latex');
    ylabel(ax, '$\hat{\alpha}$', 'FontSize', 16, 'Interpreter', 'latex');
    title(ax, 'GranE - Attenuation vs Angular Frequency');
    legend(ax, 'show', 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 12);

    outName = fullfile(pwd, 'attenuation_vs_omega.png');
    exportgraphics(ax, outName);
end
