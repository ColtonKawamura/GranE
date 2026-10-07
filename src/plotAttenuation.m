function out = plotAttenuation(varargin)
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
%                 pressure_actual) produced by processData.m
%   gammaValues - set of gamma values for marker sizes (default: from data)
%   options.plotFlag    - render a figure (default true)
%   options.lightMode   - white palette (default false)
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
		gammaValues = sort(unique(data.gamma(:)));
	end

	attenField   = 'attenuation_x';
	attenR2Field = 'attenuation_r_squared_x';

	% --- Pressure list setup ---
	if ~isempty(options.pressureArray)
		pressureList = options.pressureArray;
	else
		pressureList = unique(data.pressure_actual(:))';
	end
	pressureList = sort(pressureList);

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

	% --- Normalise pressure to a blue->red RGB ramp (normVarColor convention) ---
	% colour = [n, 0, 1-n]; dark-mode blends toward white.
	normPressure = (log(pressureList) - min(log(pressureList))) ./ ...
		(max(log(pressureList)) - min(log(pressureList)));
	normPressure(isempty(normPressure)) = 0;
	normPressure = normPressure(:);
	colourList = [normPressure, zeros(numel(normPressure), 1), (1 - normPressure)];
	if options.lightMode
		colourList = colourList + 0.2*(1 - colourList);
		colourList = min(colourList, 1);
	end

	markerSizes = exp(gammaValues/max(gammaValues))*3;

	% --- One point per unique (pressure, gamma) pair ---
	usedP = unique(data.pressure_actual(:));
	usedG = unique(data.gamma(:));
	hLines = [];

	for pi = 1:numel(usedP)
		pv = usedP(pi);
		maskP = (data.pressure_actual == pv);
		attP = data.(attenField)(maskP);
		wP   = data.omega(maskP);
		gP   = data.gamma(maskP);
		colourP = colourList(pi,:);

		for gi = 1:numel(usedG)
			gv = usedG(gi);
			% maskP indexes the FULL arrays; build a local index vector first so
			% we can subset attP/wP with the gamma condition locally.
			idxAll = find(maskP);
			subG   = data.gamma(idxAll) == gv;
			idxLocal = find(subG);
			if isempty(idxLocal)
				continue;
			end
			markerSize = markerSizes(gi);
			attVal = mean(attP(idxLocal));
			wVal   = mean(wP(idxLocal));

			% Guard the log plot: clamp non-positive attenuation.
			if attVal <= 0
				attVal = 1e-6;
			end

			if options.plotFlag
				hLines(end+1) = plot(ax, wVal, attVal, '-o', ...
					'MarkerSize', markerSize, ...
					'LineWidth', 1.2, ...
					'Color', colourP, ...
					'DisplayName', sprintf(' %.4f, %.4f', pv, gv));
			end
		end
	end

	if options.plotFlag
		leg = legend(ax, 'show', 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 15);
		leg.Location = 'eastoutside';
		title(leg, '$  \hat{P}, \hat{\gamma} $');
		axis square;
		legend(ax,'Location','northeast','Interpreter','latex','FontSize',15);
		grid(ax, 'on');
		box(ax, 'on');
		if options.lightMode
			theme(gcf, 'light');
		end
	end
end
