function vecScatterY = plotPhase(structData, structOptions)
%PLOPHASE Plot the (wrapped or unwrapped) phase vs distance for a given
%simulation data set. Ports GranMA src/matlab_functions/plotPhase.m into
%GranE (hungarian variable naming).
%
%   Args:
%       structData : scalar struct with simulation data fields, e.g.
%                    initial_distance_from_oscillation_output_x/y/z_fft,
%                    unwrapped_phase_vector_x/y/z, wavespeed_x/y, omega_gamma,
%                    pressure_actual, gamma, omega, seed.
%       structOptions : name-value pairs
%            plot        (logical, default true): generate plots.
%            shear       (logical, default false): shear vs compression.
%            z           (logical, default false): plot the z direction instead.
%            wrapped     (logical, default false): use mod 2pi instead of unwrapped.
%            lightMode   (logical, default false): light theme.
%            showFitLine (logical, default true): show the linear fit line.

    arguments
        structData (1,1) struct % Required: scalar struct with simulation data
        structOptions.plot (1,1) logical = true % Optional: flag to generate plots
        structOptions.shear (1,1) logical = false % Optional: shear vs compression for parallel component
        structOptions.z (1,1) logical = false % Optional: plots z direction instead
        structOptions.wrapped (1,1) logical = false % uses wrapped phase instead of mod 2pi
        structOptions.lightMode (1,1) logical = false % Optional: flag for light mode
        structOptions.showFitLine (1,1) logical = true % Optional: flag to show fit line
    end

    if structOptions.z
        vecDistY = structData.initial_distance_from_oscillation_output_z_fft{1};
        vecDistX = structData.initial_distance_from_oscillation_output_x_fft{1};
        vecPhaseY = structData.unwrapped_phase_vector_z{1};
        vecPhaseX = structData.unwrapped_phase_vector_x{1};
        if structOptions.wrapped
            [vecDistSortedX, vecIdx] = sort(vecDistY);
            vecPhaseY = mod(vecPhaseY(vecIdx), 2*pi);
            [vecDistSortedX, vecIdx] = sort(vecDistX);
            vecPhaseX = mod(vecPhaseX(vecIdx), 2*pi);
            vecScatterX = meanPhaseDev(vecDistX, vecPhaseX, 1);
            vecScatterY = meanPhaseDev(vecDistY, vecPhaseY, 1);
        end
    else
        vecDistY = structData.initial_distance_from_oscillation_output_y_fft{1};
        vecDistX = structData.initial_distance_from_oscillation_output_x_fft{1};
        vecPhaseY = structData.unwrapped_phase_vector_y{1};
        vecPhaseX = structData.unwrapped_phase_vector_x{1};
        if structOptions.wrapped
            [vecDistY, vecIdx] = sort(vecDistY);
            vecPhaseY = mod(vecPhaseY(vecIdx), 2*pi);
            [vecDistX, vecIdx] = sort(vecDistX);
            vecPhaseX = mod(vecPhaseX(vecIdx), 2*pi);
            vecScatterX = meanPhaseDev(vecDistX, vecPhaseX, 1);
            vecScatterY = meanPhaseDev(vecDistY, vecPhaseY, 1);
        end
    end

    % shear would swap scatter_x/scatter_y and distance_x/distance_y; left as a
    % no-op placeholder, matching the original GranMA behaviour.

    if structOptions.plot
        figure
        if structOptions.shear % this does not actually do anything, might need to change
            scatter(vecDistX, vecPhaseX, "DisplayName", "$ \hat{x}")
            hold on
            scatter(vecDistY, vecPhaseY, "DisplayName", "$ \hat{y}")
            structCoeffs = polyfit(vecDistY, vecPhaseY, 1);
            vecFitLine = linspace(min(vecDistY), max(vecDistY), 100);
            if structOptions.showFitLine
                plot(vecFitLine, polyval(structCoeffs, vecFitLine), "DisplayName", sprintf("fit $\\hat{y}$: slope = %g", structCoeffs(1)))
            end
            hold on
        else
            scatter(vecDistX, vecPhaseX, "DisplayName", "$ \hat{x}")
            hold on
            scatter(vecDistY, vecPhaseY, "DisplayName", "$ \hat{y}")
            structCoeffs = polyfit(vecDistX, vecPhaseX, 1);
            vecFitLine = linspace(min(vecDistX), max(vecDistX), 100);
            if structOptions.showFitLine
                plot(vecFitLine, polyval(structCoeffs, vecFitLine), "DisplayName", sprintf("fit $\\hat{x}$: slope = %g", structCoeffs(1)))
            end
            hold on
        end
        hold on
        box on
        grid on
        if structOptions.wrapped
            set(gca,'YTick', [0, pi, 2*pi], 'YTickLabel', {'0', ' $\pi $', '$ 2\pi $'}, 'TickLabelInterpreter', 'latex');
        end
        ylabel("$ \Delta \phi $", "Interpreter", 'latex', "FontSize", 20)
        xlabel("$ x ", "Interpreter", 'latex', "FontSize", 20)
        set(get(gca, 'ylabel'), 'rotation', 0);
        legend('Interpreter', 'latex', 'FontSize', 20);
        if structOptions.wrapped
            ylim([0,2*pi])
        end

        % Wave speed ĉ: prefer a processed wavespeed_x/y; GranE sim outputs
        % store wavenumber instead, so derive c = omega / wavenumber when the
        % wavespeed field is absent (never throws on missing data fields).
        if structOptions.shear
            vecDir = 'y';
        else
            vecDir = 'x';
        end
        scalC = NaN;
        if isfield(structData, ['wavespeed_' vecDir]) && numel(structData.(['wavespeed_' vecDir])) >= 1
            scalC = structData.(['wavespeed_' vecDir])(1);
        elseif isfield(structData, 'omega') && isfield(structData, ['wavenumber_' vecDir]) && numel(structData.(['wavenumber_' vecDir])) >= 1 && structData.(['wavenumber_' vecDir])(1) ~= 0
            scalC = structData.omega(1) / structData.(['wavenumber_' vecDir])(1);
        end
        if structOptions.shear
            title({sprintf('$\\hat{P}$: %g, $\\hat{\\gamma}$: %g, $\\hat{\\omega}$: %g, $\\textrm{seed}$: %g', structData.pressure_actual(1), structData.gamma(1), structData.omega(1), structData.seed(1)), sprintf('$ \\hat{\\gamma}$: %g, $\\hat{\\omega}\\hat{\\gamma}$: %g, $\\hat{c}$: %g', structData.gamma(1), structData.omega_gamma(1), scalC)}, 'Interpreter', 'latex', 'FontSize', 20);
        else
            title({sprintf('$\\hat{P}$: %g, $\\hat{\\gamma}$: %g, $\\hat{\\omega}$: %g, $\\textrm{seed}$: %g', structData.pressure_actual(1), structData.gamma(1), structData.omega(1), structData.seed(1)), sprintf('$ \\hat{\\gamma}$: %g, $\\hat{\\omega}\\hat{\\gamma}$: %g, $\\hat{c}$: %g', structData.gamma(1), structData.omega_gamma(1), scalC)}, 'Interpreter', 'latex', 'FontSize', 20);
        end

        if structOptions.lightMode
            try
                theme(gcf, 'light')
            catch
            end
        end
    end

end
