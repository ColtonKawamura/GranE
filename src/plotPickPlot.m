function hFigNew = plotPickPlot(scalRowStart, scalRowEnd, ...
                                      scalColStart, scalColEnd, ...
                                      hFigSource)
% plotPickPlot  Copy a block of subplots to a new figure and re-layout.
%
% Usage:
%   % copy just row 3, col 1 from current figure
%   plotPickPlot(3, 3, 1, 1);
%
%   % copy rows 2–3, cols 1–3 from a specific figure
%   plotPickPlot(2, 3, 1, 3, hFig);

    %-------------------------------
    % Defaults / inputs
    %-------------------------------
    if nargin < 5 || isempty(hFigSource)
        hFigSource = gcf;
    end

    %-------------------------------
    % Get all axes in source figure
    %-------------------------------
    hAxAll = findall(hFigSource, 'Type', 'axes');

    % If you have legends/colorbars you want to ignore, you can filter them here
    % hAxAll = hAxAll(~strcmp(get(hAxAll, 'Tag'), 'legend'));

    scalNAx = numel(hAxAll);
    if scalNAx == 0
        error('plotPickPlot:NoAxes', 'No axes found in source figure.');
    end

    cellPosition = get(hAxAll, 'Position');
    matPosition  = cell2mat(cellPosition);   % [left bottom width height] rows

    vecCenterX = matPosition(:,1) + matPosition(:,3)/2;
    vecCenterY = matPosition(:,2) + matPosition(:,4)/2;

    %----------------------------------------
    % Sort axes into a regular grid (row/col)
    %----------------------------------------
    % Sort top-to-bottom (vecCenterY descending), then left-to-right (vecCenterX ascending)
    [~, vecOrder] = sortrows([-vecCenterY, vecCenterX]);
    hAxOrdered    = hAxAll(vecOrder);
    vecCenterX    = vecCenterX(vecOrder);
    vecCenterY    = vecCenterY(vecOrder);

    % Detect unique rows and columns (using a small tolerance)
    scalTol   = 1e-4;
    vecRowY   = uniquetol(vecCenterY, scalTol);
    vecColX   = uniquetol(vecCenterX, scalTol);

    scalNRows = numel(vecRowY);
    scalNCols = numel(vecColX);

    if scalNRows * scalNCols ~= scalNAx
        warning('plotPickPlot:GridGuess', ...
            'Guessed grid is %d x %d but there are %d axes.', ...
             scalNRows, scalNCols, scalNAx);
    end

    % Put into [row, col] grid, row 1 = top
    matHAxGrid = reshape(hAxOrdered, scalNCols, scalNRows).';

    %-------------------------------
    % Bounds checking for selection
    %-------------------------------
    if scalRowStart < 1 || scalRowEnd > scalNRows || scalRowStart > scalRowEnd
        error('Invalid row range [%d, %d] for %d rows.', ...
              scalRowStart, scalRowEnd, scalNRows);
    end
    if scalColStart < 1 || scalColEnd > scalNCols || scalColStart > scalColEnd
        error('Invalid column range [%d, %d] for %d columns.', ...
              scalColStart, scalColEnd, scalNCols);
    end

    %-------------------------------
    % Select the block of axes
    %-------------------------------
    matHAxSel = matHAxGrid(scalRowStart:scalRowEnd, ...
                           scalColStart:scalColEnd);

    % Row-major flatten: row1 all cols, then row2, ...
    vecHAxSel = matHAxSel.';
    vecHAxSel = vecHAxSel(:);

    %-------------------------------
    % Copy to new figure
    %-------------------------------
    hFigNew  = figure;
    vecHAxNew = copyobj(vecHAxSel, hFigNew);

    %-------------------------------
    % Re-layout to fill new figure
    %-------------------------------
    scalNRowsNew = scalRowEnd - scalRowStart + 1;
    scalNColsNew = scalColEnd - scalColStart + 1;

    scalW = 1 / scalNColsNew;
    scalH = 1 / scalNRowsNew;

    for scalIdx = 1:numel(vecHAxNew)
        scalRowRel = ceil(scalIdx / scalNColsNew);               % 1..scalNRowsNew
        scalColRel = scalIdx - (scalRowRel-1)*scalNColsNew;      % 1..scalNColsNew

        % "subplot"-like positions: no margins
        vecPosition = [ (scalColRel-1)*scalW, ...
                        1 - scalRowRel*scalH, ...
                        scalW, ...
                        scalH ];

        set(vecHAxNew(scalIdx), 'Position', vecPosition);
    end
end
