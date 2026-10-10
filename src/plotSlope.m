function plotSlope()
	dcm = datacursormode(gcf);
	datacursormode on;
	info = getCursorInfo(dcm);

	[p1, p2] = info.Position;
	x = [p1(1), p2(1)];
	y = [p1(2), p2(2)];

	% Take log of x and y for slope calculation in log-log space
	logx = log10(x);
	logy = log10(y);

	% Fit line in log-log space
	coeffs = polyfit(logx, logy, 1);

	slope = coeffs(1);           % this is the log-log slope
	intercept = coeffs(2);       % log10(intercept) of power-law

	% Label
	label = sprintf(' slope = %.4f', slope);

	% Plot original (non-log) data using loglog
	hold on
	loglog(x, y, '-o', 'DisplayName', label, 'LineWidth', 2)
	legend show
end


