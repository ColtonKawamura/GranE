function plotAttenuation(varargin)
%PLOTATTENUATION Plot attenuation coefficient alpha vs angular frequency omega,
%coloured by applied pressure. Mirrors GranMA/src/matlab_functions/
%plotAttenuationOmega.m in both logic and aesthetics.
%
%   Simplest call (file path):
%       o.pressureArray = [0.001 0.01 0.1];   % optional
%       o.lightMode   = true;                  % optional (default false)
%       plotAttenuation('data/simMD/3d/hooke/combined_simData.mat', o);
%
%   Direct-struct call (GranMA style):
%       plotAttenuation(data, gammaValues, o);
%
%   data        - scalar struct (fields: attenuation_x, omega, gamma,
%                 pressure) produced by processData.m
%   gammaValues - gamma values to plot; also sets marker sizes (default: from data)
%   options.plotFlag    - render a figure (default true)
%   options.lightMode   - lighter palette + light theme (default false)
%   options.pressureArray - pressures to plot (default: from data)

	% --- Resolve inputs: allow file-path, struct, or mixed forms ---
	nArgs = numel(varargin);
	if nArgs == 1
		if ischar(varargin{1}) || isstring(varargin{1})
			data = load(varargin{1});
		else
			data = varargin{1};
		end
		gammaValues = [];
		options = struct();
	elseif nArgs == 2
		if ischar(varargin{1}) || isstring(varargin{1})
			data = load(varargin{1});
			options = varargin{2};
			gammaValues = [];
		elseif isstruct(varargin{2})
			data = varargin{1};
			options = varargin{2};
			gammaValues = [];
		else
			data = varargin{1};
			gammaValues = varargin{2};
			options = struct();
		end
	elseif nArgs == 3
		data = varargin{1};
		gammaValues = varargin{2};
		options = varargin{3};
	else
		error('plotAttenuation:wrongNumArgs', 'Expected 1-3 arguments.');
	end

	if ~isfield(options,'plotFlag'), options.plotFlag = true; end
	if ~isfield(options,'lightMode'), options.lightMode = false; end
	if ~isfield(options,'pressureArray'), options.pressureArray = []; end

	% If no gamma list passed, derive from the data
	if isempty(gammaValues)
		gammaValues = unique(data.gamma(:));
	end
	gammaValues = sort(gammaValues(:))';

	attenField = 'attenuation_x';

	% --- Pressure list setup (GranMA filters on the input pressure) ---
	if ~isempty(options.pressureArray)
		pressureList = options.pressureArray;
	else
		pressureList = unique(data.pressure(:));
	end
	pressureList = sort(pressureList(:))';

	if options.plotFlag
		figure_attenuation = figure('Name','Attenuation vs Angular Frequency','Color','w');
		ax = axes('Parent', figure_attenuation);
		hold(ax, 'on');

		ylabel(ax, '$ \hat{\alpha} $', 'FontSize', 20, 'Interpreter', 'latex');
		xlabel(ax, '$\hat{\omega}$', 'FontSize', 20, 'Interpreter', 'latex');
		set(ax, 'XScale', 'log');
		set(ax, 'YScale', 'log');
		set(get(ax, 'ylabel'), 'rotation', 0);
		grid(ax, 'on');
		box(ax, 'on');
	end

	% --- Normalise log-pressure to a blue->red RGB ramp (normVarColor convention) ---
	% colour = [n, 0, 1-n]; lightMode blends the colours toward white.
	logP = log(pressureList);
	logRange = max(logP) - min(logP);
	if logRange > 0
		normPressure = (logP - min(logP)) ./ logRange;
	else
		normPressure = zeros(size(logP));   % single pressure: avoid 0/0
	end
	normPressure = normPressure(:);
	colourList = [normPressure, zeros(numel(normPressure), 1), (1 - normPressure)];
	if options.lightMode
		colourList = colourList + 0.2*(1 - colourList);
		colourList = min(colourList, 1);
	end

	% --- One curve per (pressure, gamma) pair, looping over the requested
	% lists exactly as GranMA's plotAttenuationOmega does.
	for iP = 1:numel(pressureList)
		pressureValue = pressureList(iP);
		maskP = isClose(data.pressure, pressureValue);
		if ~any(maskP)
			continue;
		end
		colourP = colourList(iP,:);

		for gammaValue = gammaValues
			markerSize = exp(gammaValue/max(gammaValues))*3;

			idx = find(maskP & isClose(data.gamma, gammaValue));
			if isempty(idx)
				continue;
			end
			gammaValueActual = data.gamma(idx(1));

			attVals = data.(attenField)(idx);
			wVals   = data.omega(idx);

			% GranMA default view: discard unphysical fits with alpha > 1.
			% Non-positive values are left in; the log axis omits them.
			keep = ~(attVals > 1);
			attVals = attVals(keep);
			wVals   = wVals(keep);
			if isempty(wVals)
				continue;
			end

			% Sort by omega so the line connects in frequency order.
			[wVals, sortIdx] = sort(wVals);
			attVals = attVals(sortIdx);

			if options.plotFlag
				plot(ax, wVals, attVals, '-o', ...
					'MarkerSize', markerSize, ...
					'LineWidth', 1.2, ...
					'Color', colourP, ...
					'DisplayName', sprintf(' %.4f, %.4f', pressureValue, gammaValueActual));
			end
		end
	end

	if options.plotFlag
		leg = legend(ax, 'show', 'Location', 'eastoutside', 'Interpreter', 'latex', 'FontSize', 15);
		leg.Location = 'eastoutside';
		title(leg, '$  \hat{P}, \hat{\gamma} $');
		axis square;
		grid(ax, 'on');
		box(ax, 'on');
		if options.lightMode
			try
				theme(gcf, 'light');   % theme() only exists in R2025a+
			catch
			end
		end
	end
end

function tf = isClose(values, target)
	% Tolerant equality for floating-point parameter matching.
	tf = abs(values - target) <= 1e-9 * max(1, abs(target));
end
