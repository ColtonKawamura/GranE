function [fitted_attenuation, wavenumber, attenuation_fit_line, initial_distance_from_oscillation_output, amplitude_vector, unwrapped_phase_vector, cleaned_particle_index, initial_position_y_out, initial_position_z_out] = processDFT(dft_coefficients, particle_variance, driving_amplitude, index_particles, index_oscillating_wall, initial_distance_from_oscillation, initial_position_y, initial_position_z)

% Define minimum peak amplitude
min_peak_amplitude = driving_amplitude * 1E-10;

% Minimum variance to consider particle as having moved
min_variance = 0;

iskip = 1;

% Debug
fprintf('processDFT: N=%d, non-wall=%d, nonzero variance=%d, nonzero dft=%d\n', ...
    length(dft_coefficients), sum(~index_oscillating_wall), ...
    sum(particle_variance > min_variance), sum(abs(dft_coefficients) > 0));
fprintf('driving_amplitude = %g, min_peak_amplitude = %g\n', driving_amplitude, min_peak_amplitude);
fprintf('sample dft mag: %g %g %g\n', abs(dft_coefficients(1)), abs(dft_coefficients(2)), abs(dft_coefficients(3)));

for nn = index_particles(:)'
    if ~index_oscillating_wall(nn)
        amp = 2 * abs(dft_coefficients(nn));
        fprintf('first non-wall: nn=%d, var=%.6e, amp=%.6e, min_peak=%.6e\n', nn, particle_variance(nn), amp, min_peak_amplitude);
        break
    end
end

% Ensure figures directory exists
fig_dir = 'figures';
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end
fig_label = datestr(now, 'yyyymmdd_HHMMSS_FFF');

% Initialize output vectors
initial_distance_from_oscillation_output = [];
initial_position_y_out = [];
initial_position_z_out = [];
amplitude_vector = [];
phase_vector = [];
cleaned_particle_index = [];

has_yz = nargin >= 8 && ~isempty(initial_position_y);

for nn = index_particles(1:iskip:end)'
    if ~index_oscillating_wall(nn)
        if particle_variance(nn) > min_variance
            amp = 2 * abs(dft_coefficients(nn));
            if amp > min_peak_amplitude
                amplitude_vector = [amplitude_vector, amp];
                initial_distance_from_oscillation_output = [initial_distance_from_oscillation_output, initial_distance_from_oscillation(nn)];
                phase_vector = [phase_vector, angle(dft_coefficients(nn))];
                cleaned_particle_index = [cleaned_particle_index, nn];
                if has_yz
                    initial_position_y_out = [initial_position_y_out, initial_position_y(nn)];
                    initial_position_z_out = [initial_position_z_out, initial_position_z(nn)];
                end
            end
        end
    end
end

fprintf('processDFT: particles passing all filters = %d\n', length(amplitude_vector));

n = length(amplitude_vector);
cut_off_index = floor(1 * n);

amplitude_vector = amplitude_vector(1:cut_off_index);
initial_distance_from_oscillation_output = initial_distance_from_oscillation_output(1:cut_off_index);
phase_vector = phase_vector(1:cut_off_index);
cleaned_particle_index = cleaned_particle_index(1:cut_off_index);
if has_yz
    initial_position_y_out = initial_position_y_out(1:cut_off_index);
    initial_position_z_out = initial_position_z_out(1:cut_off_index);
end

% Attenuation Fitting and Plotting
coefficients = polyfit(initial_distance_from_oscillation_output, log(abs(amplitude_vector)), 1);
fitted_attenuation = coefficients(1);
intercept_attenuation = coefficients(2);
attenuation_fit_line = exp(intercept_attenuation) * exp(initial_distance_from_oscillation_output .* fitted_attenuation);

fig1 = figure('Visible', 'off');
semilogy(initial_distance_from_oscillation_output, abs(amplitude_vector), 'bo', 'DisplayName', 'Data');
hold on;
semilogy(initial_distance_from_oscillation_output, attenuation_fit_line, 'r-', 'DisplayName', 'Linear Fit');
xlabel('Distance from Oscillation (Particle Diameters)', 'Interpreter', 'latex');
ylabel('Particle Oscillation Amplitude', 'Interpreter', 'latex');
legend('show', 'Interpreter', 'latex');
grid on;
box on;
hold off;
saveas(fig1, fullfile(fig_dir, sprintf('amplitude_%s.png', fig_label)));
close(fig1);

% Wavenumber and Phase Fitting and Plotting
unwrapped_phase_vector = unwrap(phase_vector);

p = polyfit(initial_distance_from_oscillation_output, unwrapped_phase_vector, 1);
fitted_line = polyval(p, initial_distance_from_oscillation_output);
wavenumber = p(1);

fig2 = figure('Visible', 'off');
scatter(initial_distance_from_oscillation_output, unwrapped_phase_vector, 'o');
grid on;
box on;
hold on;
plot(initial_distance_from_oscillation_output, fitted_line, '-r');
xlabel('Distance from Oscillation (Particle Diameters)', 'Interpreter', 'latex');
ylabel('$\Delta\phi$', 'Interpreter', 'latex');

y_max = max(unwrapped_phase_vector);
y_min = min(unwrapped_phase_vector);
yticks = ceil(y_min/pi) * pi:pi:floor(y_max/pi) * pi;
yticklabels = arrayfun(@(x) sprintf('%.2f\\pi', x/pi), yticks, 'UniformOutput', false);
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);
hold off;
saveas(fig2, fullfile(fig_dir, sprintf('phase_%s.png', fig_label)));
close(fig2);

end


