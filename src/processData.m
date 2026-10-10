function outData = processData(inputDir, outputFile, options)
%PROCESSDATA Load every .mat produced by simMD in inputDir and assemble them
% into a single combined data structure (one row per simulation), mirroring
% the GranMA script src/matlab_functions/dataCrunch.m.
%
%   outData = processData(inputDir, outputFile, options)
%
%   inputDir    - folder containing the per-run simMD .mat files
%   outputFile  - path of the single combined .mat to write (with -struct)
%   options.threeD - logical, include the z-direction fields (default true)
%
% Files in inputDir that are not complete simMD outputs (e.g. a previously
% written combined file, or spectrum-only files from an early exit) are
% skipped with a warning.
%
% The combined struct uses the same field names expected by the plot*
% helpers (attenuation_x, omega, gamma, pressure, pressure_actual, seed, ...).

    arguments
        inputDir (1,1) string
        outputFile (1,1) string
        options (1,1) struct
    end
    if ~isfield(options, 'threeD'), options.threeD = true; end

    requiredVars = {'attenuation_x_dimensionless', 'attenuation_y_dimensionless', ...
        'wavenumber_x_dimensionless', 'wavenumber_y_dimensionless', ...
        'gamma_dimensionless', 'driving_angular_frequency_dimensionless', ...
        'input_pressure', 'pressure_dimensionless', 'seed', 'wavespeed_x', ...
        'amplitude_vector_x', 'amplitude_vector_y', ...
        'unwrapped_phase_vector_x', 'unwrapped_phase_vector_y', ...
        'initial_distance_from_oscillation_output_x_fft', ...
        'initial_distance_from_oscillation_output_y_fft'};

    % Collect the simMD outputs in inputDir, skipping the output file itself
    % and anything that lacks the simMD result variables.
    allFiles = dir(fullfile(inputDir, '*.mat'));
    [~, outName, outExt] = fileparts(outputFile);
    keep = false(numel(allFiles), 1);
    for i = 1:numel(allFiles)
        if strcmp(allFiles(i).name, [char(outName) char(outExt)])
            continue;
        end
        vars = who('-file', fullfile(allFiles(i).folder, allFiles(i).name));
        missing = setdiff(requiredVars, vars);
        if isempty(missing)
            keep(i) = true;
        else
            warning('processData:skipFile', 'Skipping %s (missing %s)', ...
                allFiles(i).name, strjoin(missing, ', '));
        end
    end
    files = allFiles(keep);
    nFiles = numel(files);

    % Detect max-amp tracking data from the first file (as GranMA does).
    hasMaxAmp = false;
    if nFiles > 0
        vars = who('-file', fullfile(files(1).folder, files(1).name));
        hasMaxAmp = all(ismember({'vecPosX0', 'maxAmpXAfter2', 'maxAmpYAfter2', ...
            'attenuationMaxAmpX_dimensionless', 'attenuationMaxAmpY_dimensionless'}, vars));
    end

    outData = struct();
    outData.alphaoveromega_x     = zeros(nFiles, 1);
    outData.alphaoveromega_y     = zeros(nFiles, 1);
    outData.amplitude_vector_x   = cell(nFiles, 1);
    outData.amplitude_vector_y   = cell(nFiles, 1);
    outData.attenuation_x        = zeros(nFiles, 1);
    outData.attenuation_y        = zeros(nFiles, 1);
    outData.gamma                = zeros(nFiles, 1);
    outData.initial_distance_from_oscillation_output_x_fft = cell(nFiles, 1);
    outData.initial_distance_from_oscillation_output_y_fft = cell(nFiles, 1);
    outData.omega                = zeros(nFiles, 1);
    outData.omega_gamma          = zeros(nFiles, 1);
    outData.pressure             = zeros(nFiles, 1);
    outData.pressure_actual      = zeros(nFiles, 1);
    outData.seed                 = zeros(nFiles, 1);
    outData.unwrapped_phase_vector_x = cell(nFiles, 1);
    outData.unwrapped_phase_vector_y = cell(nFiles, 1);
    outData.wavenumber_x         = zeros(nFiles, 1);
    outData.wavenumber_y         = zeros(nFiles, 1);
    outData.wavespeed_x          = zeros(nFiles, 1);

    if hasMaxAmp
        outData.x0                 = cell(nFiles, 1);
        outData.maxAmpXAfter2      = cell(nFiles, 1);
        outData.maxAmpYAfter2      = cell(nFiles, 1);
        outData.attenuationMaxAmpX = zeros(nFiles, 1);
        outData.attenuationMaxAmpY = zeros(nFiles, 1);
    end

    if options.threeD
        outData.alphaoveromega_z   = NaN(nFiles, 1);
        outData.amplitude_vector_z = cell(nFiles, 1);
        outData.attenuation_z      = NaN(nFiles, 1);
        outData.initial_distance_from_oscillation_output_z_fft = cell(nFiles, 1);
        outData.unwrapped_phase_vector_z = cell(nFiles, 1);
        outData.wavenumber_z       = NaN(nFiles, 1);
        outData.x_fft_initial_y    = cell(nFiles, 1);
        outData.x_fft_initial_z    = cell(nFiles, 1);
        outData.y_fft_initial_y    = cell(nFiles, 1);
        outData.y_fft_initial_z    = cell(nFiles, 1);
        outData.z_fft_initial_y    = cell(nFiles, 1);
        outData.z_fft_initial_z    = cell(nFiles, 1);
    end

    for i = 1:nFiles
        d = load(fullfile(files(i).folder, files(i).name));

        % Rename / derive the fields to the plotting convention (same
        % conversions as GranMA dataCrunch.m). The fitted slopes from
        % processDFT are negative for a decaying wave, so negate them.
        attenuation_x = -d.attenuation_x_dimensionless;
        attenuation_y = -d.attenuation_y_dimensionless;
        omega         = d.driving_angular_frequency_dimensionless;
        gamma         = d.gamma_dimensionless;

        outData.alphaoveromega_x(i)   = attenuation_x / omega;
        outData.alphaoveromega_y(i)   = attenuation_y / omega;
        outData.amplitude_vector_x{i} = d.amplitude_vector_x;
        outData.amplitude_vector_y{i} = d.amplitude_vector_y;
        outData.attenuation_x(i)      = attenuation_x;
        outData.attenuation_y(i)      = attenuation_y;
        outData.gamma(i)              = gamma;
        outData.initial_distance_from_oscillation_output_x_fft{i} = d.initial_distance_from_oscillation_output_x_fft;
        outData.initial_distance_from_oscillation_output_y_fft{i} = d.initial_distance_from_oscillation_output_y_fft;
        outData.omega(i)              = omega;
        outData.omega_gamma(i)        = omega * gamma;
        outData.pressure(i)           = d.input_pressure;
        outData.pressure_actual(i)    = d.pressure_dimensionless;
        outData.seed(i)               = d.seed;
        outData.unwrapped_phase_vector_x{i} = d.unwrapped_phase_vector_x;
        outData.unwrapped_phase_vector_y{i} = d.unwrapped_phase_vector_y;
        outData.wavenumber_x(i)       = -d.wavenumber_x_dimensionless;
        outData.wavenumber_y(i)       = -d.wavenumber_y_dimensionless;
        outData.wavespeed_x(i)        = -d.wavespeed_x;

        if hasMaxAmp
            outData.x0{i}                 = d.vecPosX0;
            outData.maxAmpXAfter2{i}      = d.maxAmpXAfter2;
            outData.maxAmpYAfter2{i}      = d.maxAmpYAfter2;
            outData.attenuationMaxAmpX(i) = -d.attenuationMaxAmpX_dimensionless;
            outData.attenuationMaxAmpY(i) = -d.attenuationMaxAmpY_dimensionless;
        end

        % z-direction fields are only present for 3D runs; 2D rows stay NaN/empty.
        if options.threeD && isfield(d, 'attenuation_z_dimensionless')
            attenuation_z = -d.attenuation_z_dimensionless;
            outData.alphaoveromega_z(i)   = attenuation_z / omega;
            outData.amplitude_vector_z{i} = d.amplitude_vector_z;
            outData.attenuation_z(i)      = attenuation_z;
            outData.initial_distance_from_oscillation_output_z_fft{i} = d.initial_distance_from_oscillation_output_z_fft;
            outData.unwrapped_phase_vector_z{i} = d.unwrapped_phase_vector_z;
            outData.wavenumber_z(i)       = -d.wavenumber_z_dimensionless;
            outData.x_fft_initial_y{i}    = d.x_fft_initial_y;
            outData.x_fft_initial_z{i}    = d.x_fft_initial_z;
            outData.y_fft_initial_y{i}    = d.y_fft_initial_y;
            outData.y_fft_initial_z{i}    = d.y_fft_initial_z;
            outData.z_fft_initial_y{i}    = d.z_fft_initial_y;
            outData.z_fft_initial_z{i}    = d.z_fft_initial_z;
        end
    end

    fprintf('processData: %d files -> %s\n', nFiles, outputFile);
    save(outputFile, '-struct', 'outData', '-v7.3');
end
