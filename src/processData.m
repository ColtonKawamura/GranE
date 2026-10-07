function outData = processData(inputDir, outputFile, options)
%PROCESSData Load every .mat produced by simMD in inputDir and assemble them
% into a single combined data structure (one row per simulation), mirroring
% the GranMA script src/matlab_functions/dataCrunch.m.
%
%   outData = processData(inputDir, outputFile, options)
%
%   inputDir    - folder containing the per-run simMD .mat files
%   outputFile  - path of the single combined .mat to write (with -struct)
%   options.threeD - logical, include the z-direction fields (default true)
%
% The combined struct uses the same field names expected by the plot*
% helpers (attenuation_x, omega, gamma, pressure_actual, seed, ...).

    arguments
        inputDir (1,1) string
        outputFile (1,1) string
        options (1,1) struct
    end
    if ~isfield(options, 'threeD'), options.threeD = true; end

    % Read every .mat in the directory.
    files = dir(fullfile(inputDir, '*.mat'));
    nFiles = length(files);
    outData = struct();
    outData.attenuation_x        = zeros(nFiles, 1);
    outData.attenuation_y        = zeros(nFiles, 1);
    outData.amplitude_vector_x   = cell(nFiles, 1);
    outData.amplitude_vector_y   = cell(nFiles, 1);
    outData.gamma                = zeros(nFiles, 1);
    outData.omega                = zeros(nFiles, 1);
    outData.omega_gamma          = zeros(nFiles, 1);
    outData.pressure             = zeros(nFiles, 1);
    outData.pressure_actual      = zeros(nFiles, 1);
    outData.seed                 = zeros(nFiles, 1);
    outData.alphaoveromega_x     = zeros(nFiles, 1);
    outData.alphaoveromega_y     = zeros(nFiles, 1);
    outData.wavenumber_x         = zeros(nFiles, 1);
    outData.wavenumber_y         = zeros(nFiles, 1);
    outData.unwrapped_phase_vector_x   = cell(nFiles, 1);
    outData.unwrapped_phase_vector_y   = cell(nFiles, 1);
    outData.initial_distance_from_oscillation_output_x_fft = cell(nFiles, 1);
    outData.initial_distance_from_oscillation_output_y_fft = cell(nFiles, 1);

    if options.threeD
        outData.attenuation_z        = zeros(nFiles, 1);
        outData.alphaoveromega_z     = zeros(nFiles, 1);
        outData.wavenumber_z         = zeros(nFiles, 1);
    end

    for i = 1:nFiles
        % Load the run. load returns a struct whose fields are the saved vars.
        d = load(fullfile(inputDir, files(i).name));

        % Rename / derive the fields to the plotting convention. The simMD
        % outputs are dimensionless; attenuation is stored with the opposite
        % sign of the decay we want to show, so negate it.
        d.attenuation_x = -d.attenuation_x_dimensionless';
        d.attenuation_y = -d.attenuation_y_dimensionless';
        d.gamma         = d.gamma_dimensionless';
        d.omega         = d.driving_angular_frequency_dimensionless';
        d.alphaoveromega_x = d.attenuation_x ./ d.omega;
        d.alphaoveromega_y = d.attenuation_y ./ d.omega;
        d.omega_gamma  = d.omega .* d.gamma;
        d.pressure     = d.input_pressure';
        d.pressure_actual = d.pressure_dimensionless';
        d.seed         = d.seed;

        if isfield(d, 'wavenumber_x_dimensionless')
            d.wavenumber_x = -d.wavenumber_x_dimensionless';
            d.wavenumber_y = -d.wavenumber_y_dimensionless';
        end

        % 3D fields are only present for full 3D runs.
        if options.threeD && isfield(d, 'attenuation_z_dimensionless')
            d.attenuation_z   = -d.attenuation_z_dimensionless';
            d.alphaoveromega_z = d.attenuation_z ./ d.omega;
            if isfield(d, 'wavenumber_z_dimensionless')
                d.wavenumber_z = -d.wavenumber_z_dimensionless';
            end
        end

        % Collect into the combined struct.
        outData.attenuation_x(i)          = d.attenuation_x;
        outData.attenuation_y(i)          = d.attenuation_y;
        outData.amplitude_vector_x{i}     = d.amplitude_vector_x;
        outData.amplitude_vector_y{i}     = d.amplitude_vector_y;
        outData.gamma(i)                  = d.gamma;
        outData.omega(i)                  = d.omega;
        outData.omega_gamma(i)            = d.omega_gamma;
        outData.pressure(i)               = d.pressure;
        outData.pressure_actual(i)        = d.pressure_actual;
        outData.seed(i)                   = d.seed;
        outData.alphaoveromega_x(i)       = d.alphaoveromega_x;
        outData.alphaoveromega_y(i)       = d.alphaoveromega_y;
        outData.wavenumber_x(i)           = d.wavenumber_x;
        outData.wavenumber_y(i)           = d.wavenumber_y;
        outData.unwrapped_phase_vector_x{i}     = d.unwrapped_phase_vector_x;
        outData.unwrapped_phase_vector_y{i}     = d.unwrapped_phase_vector_y;
        outData.initial_distance_from_oscillation_output_x_fft{i} = d.initial_distance_from_oscillation_output_x_fft;
        outData.initial_distance_from_oscillation_output_y_fft{i} = d.initial_distance_from_oscillation_output_y_fft;

        if options.threeD && isfield(d, 'attenuation_z_dimensionless')
            outData.attenuation_z(i)     = d.attenuation_z;
            outData.alphaoveromega_z(i)  = d.alphaoveromega_z;
            outData.wavenumber_z(i)      = d.wavenumber_z;
        end
    end

    % Save the combined data as a single struct.
    fprintf('processData: %d files -> %s\n', nFiles, outputFile);
    save(outputFile, '-struct', 'outData', '-v7.3');
end
