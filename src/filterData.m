function filtered_data_struct = filterData(data_struct, varargin)
    % Filters a scalar struct containing vector/cell array fields based on closest values.
    %
    % Args:
    %   data_struct: A scalar struct where each field contains a vector or cell array.
    %                All fields intended for filtering or output must have the same
    %                number of elements (e.g., same number of rows or total elements
    %                if vectors).
    %   varargin: Pairs of (field_name_string, value) for filtering. Filtering fields
    %             must contain numeric vector data.
    %
    % Returns:
    %   filtered_data_struct: A scalar struct containing filtered data in its fields.
    % Example:
    %      plotData = filterData(outData, 'pressure', .1, 'damping', 2)
    %      plotData = filterData(outData, 'pressure', .1, 'damping', 2, 'omega', [0.5, 1.0]);


    if ~isscalar(data_struct) || ~isstruct(data_struct)
        error('Input data must be a scalar struct.');
    end

    if mod(length(varargin), 2) ~= 0
        error('Filtering arguments must be provided in field_name/value pairs.');
    end

    all_field_names = fieldnames(data_struct);
    if isempty(all_field_names)
        filtered_data_struct = data_struct; % Return empty/original if no fields
        return;
    end

    % Determine the number of data points from the first field
    % This assumes all relevant fields have the same first dimension size
    first_field_name = all_field_names{1};
    try
        first_field_data = data_struct.(first_field_name);
        % Use size(data, 1) as the primary way to get length
        if isempty(first_field_data)
             num_elements = 0;
        else
            num_elements = size(first_field_data, 1);
            % Handle case where it might be a row vector initially
            if num_elements == 1 && ~iscell(first_field_data) && ndims(first_field_data) == 2 && size(first_field_data, 2) > 1
                num_elements = size(first_field_data, 2);
            end
        end
    catch ME
        error('Could not determine size of data from field "%s": %s', first_field_name, ME.message);
    end

    valid_indices = (1:num_elements)'; % Start with all indices valid (as column vector)

    for i = 1:2:length(varargin)
        field_name = varargin{i};
        value = varargin{i+1};
        if length(value) == 2
            if field_name == "y"
                lower_y = value(1);
                upper_y = value(2);
                for j = 1:length(data_struct.alphaoveromega_x)
                    initialY = data_struct.x_fft_initial_y{j}; % all the initial y for sim j
                    mask = (initialY >= lower_y) & (initialY <= upper_y);
                    data_struct.x_fft_initial_y{j} = initialY(mask);
                    initialZ = data_struct.x_fft_initial_z{j};
                    data_struct.x_fft_initial_z{j} = initialZ(mask);
                    amp = data_struct.amplitude_vector_x{j};
                    data_struct.amplitude_vector_x{j} = amp(mask);
                    phase = data_struct.unwrapped_phase_vector_x{j};
                    data_struct.unwrapped_phase_vector_x{j} = phase(mask);
                    initialX = data_struct.initial_distance_from_oscillation_output_x_fft{j};
                    data_struct.initial_distance_from_oscillation_output_x_fft{j} = initialX(mask); 

                    initialY = data_struct.y_fft_initial_y{j}; % all the initial y for sim j
                    mask = (initialY >= lower_y) & (initialY <= upper_y);
                    data_struct.y_fft_initial_y{j} = initialY(mask);
                    initialZ = data_struct.y_fft_initial_z{j};
                    data_struct.y_fft_initial_z{j} = initialZ(mask);
                    amp = data_struct.amplitude_vector_y{j};
                    data_struct.amplitude_vector_y{j} = amp(mask);
                    phase = data_struct.unwrapped_phase_vector_y{j};
                    data_struct.unwrapped_phase_vector_y{j} = phase(mask);
                    initialX = data_struct.initial_distance_from_oscillation_output_y_fft{j};
                    data_struct.initial_distance_from_oscillation_output_y_fft{j} = initialX(mask); 
                
                    initialY = data_struct.z_fft_initial_y{j}; % all the initial y for sim j
                    mask = (initialY >= lower_y) & (initialY <= upper_y);
                    data_struct.z_fft_initial_y{j} = initialY(mask);
                    initialZ = data_struct.z_fft_initial_z{j};
                    data_struct.z_fft_initial_z{j} = initialZ(mask);
                    amp = data_struct.amplitude_vector_z{j};
                    data_struct.amplitude_vector_z{j} = amp(mask);
                    phase = data_struct.unwrapped_phase_vector_z{j};
                    data_struct.unwrapped_phase_vector_z{j} = phase(mask);
                    initialX = data_struct.initial_distance_from_oscillation_output_z_fft{j};
                    data_struct.initial_distance_from_oscillation_output_z_fft{j} = initialX(mask); 
                end
            elseif field_name == "z"
                lower_z = value(1);
                upper_z = value(2);
                for j = 1:length(data_struct.alphaoveromega_x)
                    initialY = data_struct.x_fft_initial_y{j}; % all the initial y for sim j
                    mask = (initialY >= lower_z) & (initialY <= upper_z);
                    data_struct.x_fft_initial_y{j} = initialY(mask);
                    initialZ = data_struct.x_fft_initial_z{j};
                    data_struct.x_fft_initial_z{j} = initialZ(mask);
                    amp = data_struct.amplitude_vector_x{j};
                    data_struct.amplitude_vector_x{j} = amp(mask);
                    phase = data_struct.unwrapped_phase_vector_x{j};
                    data_struct.unwrapped_phase_vector_x{j} = phase(mask);
                    initialX = data_struct.initial_distance_from_oscillation_output_x_fft{j};
                    data_struct.initial_distance_from_oscillation_output_x_fft{j} = initialX(mask); 

                    initialY = data_struct.y_fft_initial_y{j}; % all the initial y for sim j
                    mask = (initialY >= lower_z) & (initialY <= upper_z);
                    data_struct.y_fft_initial_y{j} = initialY(mask);
                    initialZ = data_struct.y_fft_initial_z{j};
                    data_struct.y_fft_initial_z{j} = initialZ(mask);
                    amp = data_struct.amplitude_vector_y{j};
                    data_struct.amplitude_vector_y{j} = amp(mask);
                    phase = data_struct.unwrapped_phase_vector_y{j};
                    data_struct.unwrapped_phase_vector_y{j} = phase(mask);
                    initialX = data_struct.initial_distance_from_oscillation_output_y_fft{j};
                    data_struct.initial_distance_from_oscillation_output_y_fft{j} = initialX(mask); 
                
                    initialY = data_struct.z_fft_initial_y{j}; % all the initial y for sim j
                    mask = (initialY >= lower_z) & (initialY <= upper_z);
                    data_struct.z_fft_initial_y{j} = initialY(mask);
                    initialZ = data_struct.z_fft_initial_z{j};
                    data_struct.z_fft_initial_z{j} = initialZ(mask);
                    amp = data_struct.amplitude_vector_z{j};
                    data_struct.amplitude_vector_z{j} = amp(mask);
                    phase = data_struct.unwrapped_phase_vector_z{j};
                    data_struct.unwrapped_phase_vector_z{j} = phase(mask);
                    initialX = data_struct.initial_distance_from_oscillation_output_z_fft{j};
                    data_struct.initial_distance_from_oscillation_output_z_fft{j} = initialX(mask); 
                end
            elseif field_name == "omega"
                % Handle range filtering for 'omega' by updating valid_indices
                if isempty(valid_indices)
                    continue; % Skip if no valid indices left from previous filters
                end
                lower_omega = value(1);
                upper_omega = value(2);

                % --- Get and validate the full 'omega' data ---
                try
                    if ~isfield(data_struct, 'omega')
                         error('Field "omega" does not exist in the data structure.');
                    end
                    full_omega_data = data_struct.omega;

                    % Check consistency of the number of elements (using first dimension)
                    current_field_size = size(full_omega_data, 1);
                    if current_field_size == 1 && ~iscell(full_omega_data) && ndims(full_omega_data) == 2 && size(full_omega_data, 2) > 1
                        current_field_size = size(full_omega_data, 2);
                        % If it looks like a row vector, transpose it for consistency
                        if current_field_size == num_elements
                            full_omega_data = full_omega_data';
                        end
                    end

                    if ~isempty(full_omega_data) && current_field_size ~= num_elements
                        error('Field "omega" has inconsistent number of elements (%d) in the first dimension compared to expected (%d).', current_field_size, num_elements);
                    end

                    % Ensure it's numeric for comparison
                    if ~isnumeric(full_omega_data)
                        error('Filtering field "omega" must be numeric.');
                    end
                    % Ensure it's a column vector if it's a vector
                    if isvector(full_omega_data)
                        full_omega_data = full_omega_data(:);
                    elseif ~isempty(full_omega_data)
                        % If not a vector, filtering logic might need adjustment
                        error('Filtering field "omega" must be a vector (or empty).');
                    end

                catch ME
                    error('Could not extract or validate data from field "omega". Original error: %s', ME.message);
                end
                % --- End validation ---

                if isempty(full_omega_data)
                    valid_indices = []; % No data left to filter
                    continue; % Skip to next filter criterion
                end

                % Extract the subset of omega values corresponding to currently valid indices
                omega_subset = full_omega_data(valid_indices);

                % Find which indices *within the current subset* satisfy the range condition
                matching_indices_in_subset = find((omega_subset >= lower_omega) & (omega_subset <= upper_omega));

                % Update valid_indices: keep only those original indices that correspond
                % to the matching indices found in the subset.
                valid_indices = valid_indices(matching_indices_in_subset);

                continue; % Skip the default closest-value logic for this filter pair
            end
            continue
        else

            if ~ischar(field_name) && ~isstring(field_name)
                error('Field name argument must be a string.');
            end
            field_name = char(field_name);

            if isempty(valid_indices)
                break; % Stop if no valid indices remain
            end

            if ~isfield(data_struct, field_name)
                error('Field "%s" does not exist in the data structure.', field_name);
            end

            % Extract the full data for the field being filtered
            try
                full_field_data = data_struct.(field_name);
                % Check consistency of the number of elements (using first dimension)
                current_field_size = size(full_field_data, 1);
                if current_field_size == 1 && ~iscell(full_field_data) && ndims(full_field_data) == 2 && size(full_field_data, 2) > 1
                    current_field_size = size(full_field_data, 2);
                    % If it looks like a row vector, transpose it for consistency
                    if current_field_size == num_elements
                        full_field_data = full_field_data';
                    end
                end

                if ~isempty(full_field_data) && current_field_size ~= num_elements
                    error('Field "%s" has inconsistent number of elements (%d) in the first dimension compared to expected (%d).', field_name, current_field_size, num_elements);
                end

                % Ensure it's numeric for comparison
                if ~isnumeric(full_field_data)
                    error('Filtering field "%s" must be numeric.', field_name);
                end
                % Ensure it's a column vector if it's a vector
                if isvector(full_field_data)
                    full_field_data = full_field_data(:);
                elseif ~isempty(full_field_data)
                    % If not a vector, filtering logic might need adjustment
                    % For now, assume we compare based on the first column if it's a matrix?
                    % Or error? Let's error for now if it's not a vector.
                    error('Filtering field "%s" must be a vector (or empty).', field_name);
                end

            catch ME
                error('Could not extract or validate data from field "%s". Original error: %s', field_name, ME.message);
            end
        end

        % Extract the subset of values corresponding to currently valid indices
        if isempty(full_field_data)
            current_field_values = [];
        else
            current_field_values = full_field_data(valid_indices);
        end

        if isempty(current_field_values)
             valid_indices = []; % No data left to filter
             break;
        end

        if ~isnumeric(value) || ~isscalar(value)
             error('Value provided for field "%s" must be a numeric scalar for closest match filtering.', field_name);
        end

        % Calculate absolute differences for the subset
        abs_diffs = abs(current_field_values - value);

        % Find the minimum difference within the subset
        min_diff = min(abs_diffs);

        % Find which indices *within the current subset* match this minimum difference
        % Use a tolerance for floating point comparison
        tolerance = min_diff * 1e-9; % Small tolerance relative to the min difference
        matching_indices_in_subset = find(abs_diffs <= min_diff + tolerance);

        % Update valid_indices: keep only those original indices that correspond
        % to the matching indices found in the subset.
        valid_indices = valid_indices(matching_indices_in_subset);

    end % End loop through filter criteria

    % Create the final filtered struct
    filtered_data_struct = struct();
    if isempty(valid_indices)
        warning('Filtering resulted in no matching data points.');
        % Create struct with same fields but empty data
        for fn = 1:length(all_field_names)
            field = all_field_names{fn};
            original_field_data = data_struct.(field);
            if isnumeric(original_field_data)
                filtered_data_struct.(field) = [];
            elseif iscell(original_field_data)
                filtered_data_struct.(field) = {};
            else
                filtered_data_struct.(field) = []; % Default empty
            end
        end
    else
        % Select data for the final valid indices for all fields
        for fn = 1:length(all_field_names)
             field = all_field_names{fn};
             original_field_data = data_struct.(field);

             % Determine how to index based on original data dimensions
             if isempty(original_field_data)
                 filtered_data_struct.(field) = original_field_data;
             else
                 original_size = size(original_field_data);
                 elements_dim1 = original_size(1);

                 % Check if the first dimension matches the expected number of elements
                 if elements_dim1 == num_elements
                     if iscell(original_field_data)
                         filtered_data_struct.(field) = original_field_data(valid_indices, :, :, :, :, :); % Index rows
                     elseif isnumeric(original_field_data) || islogical(original_field_data)
                         filtered_data_struct.(field) = original_field_data(valid_indices, :, :, :, :, :); % Index rows
                     else
                         warning('Cannot filter field "%s" of type %s by row; copying original value.', field, class(original_field_data));
                         filtered_data_struct.(field) = original_field_data;
                     end
                 % Check if it was treated as a row vector initially
                 elseif num_elements > 1 && elements_dim1 == 1 && length(original_field_data) == num_elements
                      if iscell(original_field_data)
                         filtered_data_struct.(field) = original_field_data(:, valid_indices); % Index columns
                     elseif isnumeric(original_field_data) || islogical(original_field_data)
                         filtered_data_struct.(field) = original_field_data(:, valid_indices); % Index columns
                     else
                         warning('Cannot filter field "%s" of type %s by column; copying original value.', field, class(original_field_data));
                         filtered_data_struct.(field) = original_field_data;
                     end
                 else
                     % If dimensions don't match, cannot safely filter, copy original
                     warning('Field "%s" dimensions [%s] do not match expected element count %d for filtering; copying original value.', field, num2str(original_size), num_elements);
                     filtered_data_struct.(field) = original_field_data;
                 end
             end
        end
    end

end % End function

