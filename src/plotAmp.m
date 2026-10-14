function ampRatio = plotAmp(structData, structOptions)
%PLOAMP Plot amplitude decay and calculate the ratio of perpendicular to
%parallel amplitude intercepts from exponential fits. Ports GranMA
%src/matlab_functions/plotAmp.m into GranE (hungarian variable naming).
%
%   Args:
%       structData : scalar struct with simulation data fields, e.g.
%                    amplitude_vector_x/y/z, initial_distance_from_oscillation_output_x/y/z_fft,
%                    unwrapped_phase_vector_x/y/z, pressure_actual, gamma, omega, omega_gamma,
%                    attenuation_x/y, alphaoveromega_x, seed.
%       structOptions : name-value pairs
%            plot      (logical, default true): generate plots.
%            shear     (logical, default false): use shear (y) data instead of
%                                                compression (x) for the parallel component.
%            z         (logical, default false): plot the z direction instead.
%            slice     (logical, default false): filter particles by initial y/z.
%            showFit   (logical, default true): overlay the fitted decay line.
%            lightMode (logical, default false): light theme.
%
%   Returns:
%       ampRatio : exp(yIntercept_amp_perp) / exp(yIntercept_amp_parra)

    arguments
        structData (1,1) struct % Required: scalar struct with simulation data
        structOptions.plot (1,1) logical = true % Optional: flag to generate plots
        structOptions.shear (1,1) logical = false % Optional: shear vs compression for parallel component
        structOptions.z (1,1) logical = false % Optional: plots z direction instead
        structOptions.trimDist (1,1) logical = false % trim the distance
        structOptions.slice (1,1) logical = false % filter by initial y/z
        structOptions.showFit (1,1) logical = true
        structOptions.lightMode (1,1) logical = false
    end

    colYColor = [0.8500, 0.3250, 0.0980];
    colXColor = [0.0000, 0.4470, 0.7410];

    if structOptions.slice
        % FILTER SECTION ----------------------------------------------

        scalLowerY = 0;
        scalUpperY = 40;
        scalLowerZ = 0;
        scalUpperZ = 2;

        % Filter out particles whose initial y is outside [lowerY, upperY].
        vecFFTinitialY = structData(1).y_fft_initial_y{1}; % Extract contents of cell array
        boolMask = (vecFFTinitialY >= scalLowerY) & (vecFFTinitialY <= scalUpperY);

        structData(1).y_fft_initial_y{1} = structData(1).y_fft_initial_y{1}(boolMask);
        structData(1).y_fft_initial_z{1} = structData(1).y_fft_initial_z{1}(boolMask);
        structData(1).initial_distance_from_oscillation_output_y_fft{1} = structData(1).initial_distance_from_oscillation_output_y_fft{1}(boolMask);
        structData(1).amplitude_vector_y{1} = structData(1).amplitude_vector_y{1}(boolMask);
        structData(1).unwrapped_phase_vector_y{1} = structData(1).unwrapped_phase_vector_y{1}(boolMask);

        % Filter based on z.
        vecFFTinitialZ = structData(1).y_fft_initial_z{1}; % Extract contents of cell array
        boolMask = (vecFFTinitialZ >= scalLowerZ) & (vecFFTinitialZ <= scalUpperZ);

        structData(1).y_fft_initial_y{1} = structData(1).y_fft_initial_y{1}(boolMask);
        structData(1).y_fft_initial_z{1} = structData(1).y_fft_initial_z{1}(boolMask);
        structData(1).initial_distance_from_oscillation_output_y_fft{1} = structData(1).initial_distance_from_oscillation_output_y_fft{1}(boolMask);
        structData(1).amplitude_vector_y{1} = structData(1).amplitude_vector_y{1}(boolMask);
        structData(1).unwrapped_phase_vector_y{1} = structData(1).unwrapped_phase_vector_y{1}(boolMask);

        % Filter particles whose initial z is outside of [lowerZ, upperZ] for z_fft.
        vecFFTinitialYz = structData(1).z_fft_initial_y{1}; % Extract contents of cell array
        boolMask = (vecFFTinitialYz >= scalLowerY) & (vecFFTinitialYz <= scalUpperY);

        structData(1).z_fft_initial_y{1} = structData(1).z_fft_initial_y{1}(boolMask);
        structData(1).z_fft_initial_z{1} = structData(1).z_fft_initial_z{1}(boolMask);
        structData(1).initial_distance_from_oscillation_output_z_fft{1} = structData(1).initial_distance_from_oscillation_output_z_fft{1}(boolMask);
        structData(1).amplitude_vector_z{1} = structData(1).amplitude_vector_z{1}(boolMask);
        structData(1).unwrapped_phase_vector_z{1} = structData(1).unwrapped_phase_vector_z{1}(boolMask);

        vecFFTinitialZz = structData(1).z_fft_initial_z{1}; % Extract contents of cell array
        boolMask = (vecFFTinitialZz >= scalLowerZ) & (vecFFTinitialZz <= scalUpperZ);

        structData(1).z_fft_initial_y{1} = structData(1).z_fft_initial_y{1}(boolMask);
        structData(1).z_fft_initial_z{1} = structData(1).z_fft_initial_z{1}(boolMask);
        structData(1).initial_distance_from_oscillation_output_z_fft{1} = structData(1).initial_distance_from_oscillation_output_z_fft{1}(boolMask);
        structData(1).amplitude_vector_z{1} = structData(1).amplitude_vector_z{1}(boolMask);
        structData(1).unwrapped_phase_vector_z{1} = structData(1).unwrapped_phase_vector_z{1}(boolMask);

        % Filter particles whose initial z is outside of [lowerZ, upperZ] for x_fft.
        vecFFTinitialXy = structData(1).x_fft_initial_y{1}; % Extract contents of cell array
        boolMask = (vecFFTinitialXy >= scalLowerY) & (vecFFTinitialXy <= scalUpperY);

        structData(1).x_fft_initial_y{1} = structData(1).x_fft_initial_y{1}(boolMask);
        structData(1).x_fft_initial_z{1} = structData(1).x_fft_initial_z{1}(boolMask);
        structData(1).initial_distance_from_oscillation_output_x_fft{1} = structData(1).initial_distance_from_oscillation_output_x_fft{1}(boolMask);
        structData(1).amplitude_vector_x{1} = structData(1).amplitude_vector_x{1}(boolMask);
        structData(1).unwrapped_phase_vector_x{1} = structData(1).unwrapped_phase_vector_x{1}(boolMask);

        vecFFTinitialXz = structData(1).x_fft_initial_z{1}; % Extract contents of cell array
        boolMask = (vecFFTinitialXz >= scalLowerZ) & (vecFFTinitialXz <= scalUpperZ);

        structData(1).x_fft_initial_y{1} = structData(1).x_fft_initial_y{1}(boolMask);
        structData(1).x_fft_initial_z{1} = structData(1).x_fft_initial_z{1}(boolMask);
        structData(1).initial_distance_from_oscillation_output_x_fft{1} = structData(1).initial_distance_from_oscillation_output_x_fft{1}(boolMask);
        structData(1).amplitude_vector_x{1} = structData(1).amplitude_vector_x{1}(boolMask);
        structData(1).unwrapped_phase_vector_x{1} = structData(1).unwrapped_phase_vector_x{1}(boolMask);

        % FILTER SECTION END ----------------------------------------------
    end

    if structOptions.shear
        vecXPara = structData.initial_distance_from_oscillation_output_y_fft{1};
        vecYPara = structData.amplitude_vector_y{1};
    else
        vecXPara = structData.initial_distance_from_oscillation_output_x_fft{1};
        vecYPara = structData.amplitude_vector_x{1};
    end

    structCoeffsParra = fitLogLine(vecXPara, vecYPara);
    scalIntParra = structCoeffsParra(1);
    scalSlopeParra = structCoeffsParra(2);

    if structOptions.plot
        if structOptions.shear
            figure
            scatter(vecXPara, vecYPara,[], colYColor, "DisplayName", "$\\hat{y}$")
            hold on
            set(gca, 'YScale', 'log')
            grid on
            xlabel("$ x $", "Interpreter", 'latex', "FontSize", 30)
            ylabel("$A$", "Interpreter", 'latex', "FontSize", 30)
            set(get(gca, 'ylabel'), 'rotation', 0);
            box on
            hold on
        else
            figure
            scatter(vecXPara, vecYPara,[], colXColor,"DisplayName", "$\\hat{x}$")
            hold on
            set(gca, 'YScale', 'log')
            grid on
            xlabel("$ x $", "Interpreter", 'latex', "FontSize", 30)
            ylabel("$A$", "Interpreter", 'latex', "FontSize", 30)
            set(get(gca, 'ylabel'), 'rotation', 0);
            box on
            hold on
        end
    end

    if structOptions.z
        vecXPerp = structData.initial_distance_from_oscillation_output_z_fft{1};
        vecYPerp = structData.amplitude_vector_z{1};
        structCoeffsPerp = fitLogLine(vecXPerp, vecYPerp);
        scalIntPerp = structCoeffsPerp(1);
        scalSlopePerp = structCoeffsPerp(2);
    elseif structOptions.shear
        vecXPerp = structData.initial_distance_from_oscillation_output_x_fft{1};
        vecYPerp = structData.amplitude_vector_x{1};
        structCoeffsPerp = fitLogLine(vecXPerp, vecYPerp);
        scalIntPerp = structCoeffsPerp(1);
        scalSlopePerp = structCoeffsPerp(2);
    else
        vecXPerp = structData.initial_distance_from_oscillation_output_y_fft{1};
        vecYPerp = structData.amplitude_vector_y{1};
        structCoeffsPerp = fitLogLine(vecXPerp, vecYPerp);
        scalIntPerp = structCoeffsPerp(1);
        scalSlopePerp = structCoeffsPerp(2);
    end

    ampRatio = exp(scalIntPerp) / exp(scalIntParra);

    if structOptions.plot
        if structOptions.shear
            scatter(vecXPerp, vecYPerp,[], colXColor, "DisplayName", "$ \\hat{x}")
        else
            scatter(vecXPerp, vecYPerp,[], colYColor, "DisplayName", "$ \\hat{y}")
        end

        if structOptions.showFit
            plot(vecXPara, exp(scalIntParra + scalSlopeParra .* vecXPara), 'HandleVisibility', 'off', 'DisplayName', sprintf("$Slope: %g$", scalSlopeParra), 'Color', 'k', 'LineWidth', 5)
        end

        set(gca, 'YScale', 'log')
        grid on
        legend('show', 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 20);
        if structOptions.shear
            title(sprintf("$\\hat{P}$: %g, $\\hat{\\gamma}$: %g, $\\hat{\\omega}$: %g, Ratio: %g, \\n $\\hat{\\omega}\\hat{\\gamma}$: %g, $\\hat{\\alpha}$: %g, $\\hat{\\alpha} /\\hat{\\omega}$: %g, seed: %g", structData.pressure_actual(1), structData.gamma(1), structData.omega(1), ampRatio, structData.omega_gamma(1), structData.attenuation_y, structData.alphaoveromega_x, structData.seed(1)), 'Interpreter', 'latex','FontSize', 20);
        else
            title(sprintf("$\\hat{P}$: %g, $\\hat{\\gamma}$: %g, $\\hat{\\omega}$: %g, Ratio: %g, \\n $\\hat{\\omega}\\hat{\\gamma}$: %g, $\\hat{\\alpha}$: %g, $\\hat{\\alpha} /\\hat{\\omega}$: %g, seed: %g", structData.pressure_actual(1), structData.gamma(1), structData.omega(1), ampRatio, structData.omega_gamma(1), structData.attenuation_x, structData.alphaoveromega_x, structData.seed(1)), 'Interpreter', 'latex','FontSize', 20);
        end
        if structOptions.lightMode
            try
                theme(gcf, 'light')
            catch
            end
        end
    end

end

function structCoeffs = fitLogLine(vecX, vecY)
    % Fit log(y) = p(1) + p(2)*x and return coeffs as [intercept, slope],
    % consistent with the GranMA/Julia convention.
    polyCoeffs = polyfit(vecX, log(vecY), 1);
    structCoeffs = [polyCoeffs(2), polyCoeffs(1)];
end
