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
%   options.scaleYGamma   - exponent p: plot alpha / gamma^p on the y-axis
%                           (default 0 = unscaled). Any real p, e.g. 0.5, 1.4;
%                           negative p multiplies. GranMA convention (divide).
%   options.scaleXCoordNum - exponent q: plot omega / Z^q on the x-axis, where
%                           Z is the packing's mean coordination number
%                           (data.coord_num from processData; default 0).
%   options.scaleXOmegaC  - exponent r: plot omega / omega_c^r on the x-axis,
%                           with omega_c = sqrt(pressure_actual) per curve
%                           (default 0; r = 1 is GranMA's xOmegaOverSqrtP).
%                           Combines with scaleXCoordNum: x = omega/(omega_c^r Z^q).
%   options.seed          - plot only this packing seed (default [] = all seeds).
%                           Must match a seed in data.seed exactly; the figure
%                           is titled "Seed = <seed>" (GranMA singleSeed style).

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
	if ~isfield(options,'scaleYGamma'), options.scaleYGamma = 0; end
	if ~isfield(options,'scaleXCoordNum'), options.scaleXCoordNum = 0; end
	if ~isfield(options,'scaleXOmegaC'), options.scaleXOmegaC = 0; end
	if ~isfield(options,'seed'), options.seed = []; end
	validateattributes(options.scaleYGamma, {'numeric'}, {'scalar', 'real', 'finite'});
	validateattributes(options.scaleXCoordNum, {'numeric'}, {'scalar', 'real', 'finite'});
	validateattributes(options.scaleXOmegaC, {'numeric'}, {'scalar', 'real', 'finite'});

	% --- Seed selection (before any other filtering) ---
	if ~isempty(options.seed)
		validateattributes(options.seed, {'numeric'}, {'scalar', 'real', 'finite'});
		% filterData picks the closest value, so check for an exact match
		% first to avoid silently plotting a different seed.
		if ~any(data.seed(:) == options.seed)
			error('plotAttenuation:noSeed', 'Seed %g not found. Available seeds: %s', ...
				options.seed, mat2str(unique(data.seed(:))'));
		end
		data = filterData(data, 'seed', options.seed);
	end

	if options.scaleXCoordNum ~= 0
		if ~isfield(data, 'coord_num') || all(isnan(data.coord_num(:)))
			error('plotAttenuation:noCoordNum', ...
				['scaleXCoordNum needs data.coord_num. Re-run processData on simMD ' ...
				 'outputs that store mean_coord_num, or pass options.packingDir to processData.']);
		elseif any(isnan(data.coord_num(:)))
			warning('plotAttenuation:partialCoordNum', ...
				'%d rows have no coord_num and are omitted from the Z-scaled plot.', ...
				sum(isnan(data.coord_num(:))));
		end
	end

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

		ylabel(ax, scaledLabel('\hat{\alpha}', {'\hat{\gamma}'}, options.scaleYGamma), ...
			'FontSize', 20, 'Interpreter', 'latex');
		xlabel(ax, scaledLabel('\hat{\omega}', {'\hat{\omega}_c', 'Z'}, ...
			[options.scaleXOmegaC, options.scaleXCoordNum]), ...
			'FontSize', 20, 'Interpreter', 'latex');
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
	% lists and selecting simulations with filterData exactly as GranMA's
	% plotAttenuationOmega does (closest-value match on each field).
	for iP = 1:numel(pressureList)
		pressureValue = pressureList(iP);
		pressureData = filterData(data, 'pressure', pressureValue);
		if isempty(pressureData.omega)
			continue;
		end
		colourP = colourList(iP,:);
		omegaC = sqrt(pressureData.pressure_actual(1));   % omega_c = sqrt(P_actual)

		for gammaValue = gammaValues
			markerSize = exp(gammaValue/max(gammaValues))*3;

			gammaData = filterData(pressureData, 'gamma', gammaValue);
			if isempty(gammaData.omega)
				continue;
			end
			gammaValueActual = gammaData.gamma(1);

			attVals = gammaData.(attenField)(:);
			wVals   = gammaData.omega(:);
			if options.scaleXCoordNum ~= 0
				zVals = gammaData.coord_num(:);
			else
				zVals = ones(size(wVals));
			end

			% GranMA default view: discard unphysical fits with alpha > 1.
			% Non-positive values are left in; the log axis omits them.
			% The filter acts on the raw alpha, before any axis scaling.
			keep = ~(attVals > 1);
			attVals = attVals(keep);
			wVals   = wVals(keep);
			zVals   = zVals(keep);
			if isempty(wVals)
				continue;
			end

			% Sort by omega so the line connects in frequency order.
			[wVals, sortIdx] = sort(wVals);
			attVals = attVals(sortIdx);
			zVals   = zVals(sortIdx);

			% Axis scaling: y = alpha / gamma^p, x = omega / (omega_c^r Z^q).
			plotY = attVals ./ gammaValueActual.^options.scaleYGamma;
			plotX = wVals ./ (omegaC.^options.scaleXOmegaC .* zVals.^options.scaleXCoordNum);

			if options.plotFlag
				plot(ax, plotX, plotY, '-o', ...
					'MarkerSize', markerSize, ...
					'LineWidth', 1.2, ...
					'Color', colourP, ...
					'DisplayName', sprintf(' %.4f, %.4f', pressureValue, gammaValueActual));
			end
		end
	end

	if options.plotFlag
		if ~isempty(options.seed)
			title(ax, sprintf('Seed = %g', options.seed), 'Interpreter', 'latex', 'FontSize', 20);
		end
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

function str = scaledLabel(base, syms, powers)
	% LaTeX label for base / prod(syms{k}^powers(k)), e.g.
	% '$\frac{\hat{\alpha}}{\hat{\gamma}^{\frac{1}{2}}}$'. Positive powers go in
	% the denominator, negative powers multiply the numerator.
	num = base;
	den = '';
	for k = 1:numel(syms)
		if powers(k) > 0
			if ~isempty(den)
				den = [den '\,'];
			end
			den = [den powerTerm(syms{k}, powers(k))];
		elseif powers(k) < 0
			num = [num '\,' powerTerm(syms{k}, -powers(k))];
		end
	end
	if isempty(den)
		str = ['$' num '$'];
	else
		str = ['$\frac{' num '}{' den '}$'];
	end
end

function str = powerTerm(sym, p)
	% sym^p with p > 0; proper fractions (|p| < 1) are written as \frac{n}{d}.
	if p == 1
		str = sym;
		return;
	end
	[n, d] = rat(p, 1e-9);
	if d == 1
		expStr = sprintf('%d', n);
	elseif n < d && d <= 16 && abs(n/d - p) < 1e-9
		expStr = sprintf('\\frac{%d}{%d}', n, d);
	else
		expStr = sprintf('%g', p);
	end
	str = [sym '^{' expStr '}'];
end
