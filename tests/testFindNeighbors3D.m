% testFindNeighbors3D.m
% Unit test for findNeighbors3D / rebuildCellList3D (extracted from pack.m).
%
% Run from anywhere with:
%   matlab -batch "run('GranE/tests/testFindNeighbors3D.m')"   (from ~/repos)
% or from the project root:
%   matlab -batch "run('tests/testFindNeighbors3D.m')"
%
% Checks:
%   1. THE BUG CASE: scalNumCellsX=Y=Z = 2, one particle per cell (8 total).
%      The old scalCellLeft/Right/Up/Down/Back/Front stencil double-counted
%      neighbor cells here (a 27-reference stencil collapses to 8 distinct
%      cells, each referenced multiple times). The new 3x3x3 periodic
%      stencil + unique(...,'rows') must not.
%      Asserts: no duplicate (source,dest) pairs, source < dest for every
%      pair, and the pair set matches a hand-constructed expectation
%      (in 2x2x2 every cell pair is a stencil neighbor, so all C(8,2)=28
%      pairs must appear exactly once).
%   2. 2x2x2 with TWO particles per cell (16 total): all C(16,2)=120 pairs
%      appear exactly once, including same-cell pairs.
%   3. Mixed grid 2x2x4 (bug case on two axes, longer in z): the pair set
%      matches an independent O(scalNumPartTot^2) hand-constructed expectation built
%      from periodic minimum-image cell distances (88 pairs).
%   4. Regression at a typical size 5x5x5: every cell sees its full
%      27-cell stencil (26 other cells) -> 125*26/2 = 1625 pairs, all unique.
%   5. rebuildCellList3D buckets synthetic positions into the correct cells
%      and reports the correct cell counts.

% --- resolve src/ relative to this file so the test works from any cwd ---
thisFile = mfilename('fullpath');
[thisDir, ~, ~] = fileparts(thisFile);
srcDir = fullfile(fileparts(thisDir), 'src');   % <repo>/src
addpath(srcDir);
assert(exist(fullfile(srcDir,'findNeighbors3D.m'),'file') == 2, ...
    'findNeighbors3D.m not found in src/');

%% helper: run findNeighbors3D on a cell list and return the active pairs
function matPairs = runFinder(cellParticleList, scalNumCellsX, scalNumCellsY, scalNumCellsZ, scalNumPartTot)
    vecPairIdxSource = zeros(scalNumPartTot*12, 1);
    vecPairIdxDest   = zeros(scalNumPartTot*12, 1);
    [vecPairIdxSource, vecPairIdxDest, scalNumPairs, ~] = findNeighbors3D( ...
        cellParticleList, scalNumCellsX, scalNumCellsY, scalNumCellsZ, ...
        vecPairIdxSource, vecPairIdxDest, scalNumPartTot*12);
    % (findNeighbors3D internally doubles its buffers when the scalNumPartTot*12
    %  heuristic underestimates, so scalNumPairs may legitimately exceed scalNumPartTot*12)
    assert(scalNumPairs >= 0, 'negative pair count');
    % active entries must be valid particle indices 1..scalNumPartTot
    assert(all(vecPairIdxSource(1:scalNumPairs) >= 1) && ...
           all(vecPairIdxSource(1:scalNumPairs) <= scalNumPartTot) && ...
           all(vecPairIdxDest(1:scalNumPairs) >= 1) && ...
           all(vecPairIdxDest(1:scalNumPairs) <= scalNumPartTot), ...
           'pair index out of range');
    matPairs = [vecPairIdxSource(1:scalNumPairs), vecPairIdxDest(1:scalNumPairs)];
end

%% helper: build a 3D cell array with one particle per cell, particle k in
% cell k (column-major, x fastest) — the same layout rebuildCellList3D uses.
function cellParticleList = onePerCell(scalNumCellsX, scalNumCellsY, scalNumCellsZ)
    cellParticleList = cell(scalNumCellsX, scalNumCellsY, scalNumCellsZ);
    for k = 1:(scalNumCellsX*scalNumCellsY*scalNumCellsZ)
        [cx, cy, cz] = ind2sub([scalNumCellsX, scalNumCellsY, scalNumCellsZ], k);
        cellParticleList{cx, cy, cz} = k;
    end
end

%% helper: hand-constructed expectation from the SPEC, not the stencil code.
% One particle per cell, cell k (column-major, x fastest) holds particle k.
% A pair (a,b) is a candidate iff the two cells are within 1 (periodic
% minimum image) along EVERY axis.
function matExpPairs = expectedPairs(scalNumCellsX, scalNumCellsY, scalNumCellsZ)
    matExpPairs = zeros(0, 2);
    for a = 1:(scalNumCellsX*scalNumCellsY*scalNumCellsZ)-1
        [ax, ay, az] = ind2sub([scalNumCellsX, scalNumCellsY, scalNumCellsZ], a);
        for b = a+1:(scalNumCellsX*scalNumCellsY*scalNumCellsZ)
            [bx, by, bz] = ind2sub([scalNumCellsX, scalNumCellsY, scalNumCellsZ], b);
            dx = min(abs(ax-bx), scalNumCellsX-abs(ax-bx));
            dy = min(abs(ay-by), scalNumCellsY-abs(ay-by));
            dz = min(abs(az-bz), scalNumCellsZ-abs(az-bz));
            if dx <= 1 && dy <= 1 && dz <= 1
                matExpPairs = [matExpPairs; a, b]; %#ok<AGROW>
            end
        end
    end
end

%% Scenario 1: THE BUG CASE — 2x2x2 cells, one particle per cell (8 particles)
scalNumCellsX = 2;
scalNumCellsY = 2;
scalNumCellsZ = 2;
scalNumPartTot = scalNumCellsX*scalNumCellsY*scalNumCellsZ;

% create cell and drop only one particle in that cell
% that particle has an index = cell's linear index
cellParticleList = onePerCell(scalNumCellsX, scalNumCellsY, scalNumCellsZ);


% call findNieghbors3D and return mat of all active contact candidate pairs
% each row is [source, dest] particle index
matPairs = runFinder(cellParticleList, scalNumCellsX, scalNumCellsY, scalNumCellsZ, scalNumPartTot);

% no duplicate (source,dest) pairs
[~, ia] = unique(matPairs, 'rows');
assert(numel(ia) == size(matPairs, 1), ...
    sprintf('2x2x2: %d duplicate (source,dest) pairs found', ...
        size(matPairs,1) - numel(ia)));

% every pair satisfies source < dest
assert(all(matPairs(:,1) < matPairs(:,2)), ...
    '2x2x2: found a pair with source >= dest');

% in 2x2x2 the periodic 3x3x3 stencil covers ALL cells, so the expected
% pair set is every one of the C(8,2) = 28 pairs, exactly once
matExp = expectedPairs(scalNumCellsX, scalNumCellsY, scalNumCellsZ);
assert(size(matExp, 1) == 28, 'hand-check: 2x2x2 expectation should be 28 pairs');
assert(size(matPairs, 1) == 28, ...
    sprintf('2x2x2: got %d pairs, expected 28 (every particle pair)', size(matPairs,1)));
assert(isempty(setdiff(matPairs, matExp, 'rows')) && isempty(setdiff(matExp, matPairs, 'rows')), ...
    '2x2x2: pair set does not match hand-constructed expectation');
disp('Scenario 1 (2x2x2 bug case, 8 particles): PASSED — 28 unique pairs, none duplicated');

%% Scenario 2: 2x2x2, TWO particles per cell (16 particles)
scalNumCellsX = 2;
scalNumCellsY = 2;
scalNumCellsZ = 2;
scalNumPartTot = 16;

% put 2 particles in each cell. 
cellParticleList = cell(2, 2, 2);
k = 0;
for iz = 1:2
    for iy = 1:2
        for ix = 1:2
            k = k + 1;
            cellParticleList{ix, iy, iz} = [2*(k-1)+1; 2*(k-1)+2];
        end
    end
end
matPairs = runFinder(cellParticleList, scalNumCellsX, scalNumCellsY, scalNumCellsZ, scalNumPartTot);

[~, ia] = unique(matPairs, 'rows');
assert(numel(ia) == size(matPairs, 1), ...
    sprintf('2x2x2 dense: %d duplicate pairs found', size(matPairs,1) - numel(ia)));
assert(all(matPairs(:,1) < matPairs(:,2)), '2x2x2 dense: source >= dest pair found');

% every particle pair is a stencil neighbor (all cells touch), so all
% C(16,2) = 120 pairs must appear exactly once, including same-cell pairs
matAllPairs = zeros(0, 2);
for a = 1:15
    for b = a+1:16
        matAllPairs = [matAllPairs; a, b]; %#ok<AGROW>
    end
end
assert(size(matPairs, 1) == 120, ...
    sprintf('2x2x2 dense: got %d pairs, expected 120', size(matPairs,1)));
assert(isempty(setdiff(matPairs, matAllPairs, 'rows')) && ...
       isempty(setdiff(matAllPairs, matPairs, 'rows')), ...
    '2x2x2 dense: pair set does not match all C(16,2) pairs');
disp('Scenario 2 (2x2x2 dense, 16 particles): PASSED — 120 unique pairs incl. same-cell');

%% Scenario 3: 2x2x4 grid — bug case on x,y; long axis in z
scalNumCellsX = 2; scalNumCellsY = 2; scalNumCellsZ = 4; scalNumPartTot = scalNumCellsX*scalNumCellsY*scalNumCellsZ;
cellParticleList = onePerCell(scalNumCellsX, scalNumCellsY, scalNumCellsZ);
matPairs = runFinder(cellParticleList, scalNumCellsX, scalNumCellsY, scalNumCellsZ, scalNumPartTot);

[~, ia] = unique(matPairs, 'rows');
assert(numel(ia) == size(matPairs, 1), '2x2x4: duplicate pairs found');
assert(all(matPairs(:,1) < matPairs(:,2)), '2x2x4: source >= dest pair found');

% hand-constructed expectation from periodic minimum-image cell distances
matExp = expectedPairs(scalNumCellsX, scalNumCellsY, scalNumCellsZ);
assert(size(matExp, 1) == 88, ...
    sprintf('hand-check: 2x2x4 expectation should be 88 pairs, got %d', size(matExp,1)));
assert(size(matPairs, 1) == size(matExp, 1), ...
    sprintf('2x2x4: got %d pairs, expected %d', size(matPairs,1), size(matExp,1)));
assert(isempty(setdiff(matPairs, matExp, 'rows')) && isempty(setdiff(matExp, matPairs, 'rows')), ...
    '2x2x4: pair set does not match hand-constructed expectation');
disp('Scenario 3 (2x2x4 mixed grid, 16 particles): PASSED — 88 pairs, coverage matches');

%% Scenario 4: regression at typical size 5x5x5 (no duplicates expected)
%% verifies that old code works on larger case (regression)
scalNumCellsX = 5; scalNumCellsY = 5; scalNumCellsZ = 5; scalNumPartTot = scalNumCellsX*scalNumCellsY*scalNumCellsZ;
cellParticleList = onePerCell(scalNumCellsX, scalNumCellsY, scalNumCellsZ);
matPairs = runFinder(cellParticleList, scalNumCellsX, scalNumCellsY, scalNumCellsZ, scalNumPartTot);

[~, ia] = unique(matPairs, 'rows');
assert(numel(ia) == size(matPairs, 1), '5x5x5: duplicate pairs found');
assert(all(matPairs(:,1) < matPairs(:,2)), '5x5x5: source >= dest pair found');

% with >= 3 cells per axis every cell sees exactly 26 other stencil cells
assert(size(matPairs, 1) == scalNumPartTot*26/2, ...
    sprintf('5x5x5: got %d pairs, expected %d', size(matPairs,1), scalNumPartTot*26/2));
disp('Scenario 4 (5x5x5 regression, 125 particles): PASSED — 1625 unique pairs');

%% Scenario 5: rebuildCellList3D buckets synthetic positions correctly
%% builds a simple 4x4 cell grid, drops 6 parttyicles in known positions
%% and checks that  rebuildCellList works in pack.m, but this is hardcodeded in
scalBox = 12; scalRawCellWidth = 3;   % -> 4 cells per axis, cell width 3
% 6 particles in known corner/edge/interior cells
vecX = [1.5, 4.5, 7.5, 10.5, 1.5, 7.5]';
vecY = [1.5, 1.5, 4.5, 7.5, 10.5, 1.5]';
vecZ = [1.5, 4.5, 10.5, 1.5, 7.5, 1.5]';
scalNumPartTot = 6;
[vecXw, vecYw, vecZw, cellParticleList, nX, nY, nZ] = rebuildCellList3D( ...
    vecX, vecY, vecZ, scalBox, scalBox, scalBox, scalRawCellWidth, 0.01, scalNumPartTot);

assert(nX == 4 && nY == 4 && nZ == 4, ...
    sprintf('rebuildCellList3D: cell counts %dx%dx%d, expected 4x4x4', nX, nY, nZ));
assert(isequal(size(cellParticleList), [4, 4, 4]), 'cellParticleList shape wrong');

% cell index = ceil(pos / cellWidth), cellWidth = 12/4 = 3
expCellX = ceil(vecX / 3); expCellY = ceil(vecY / 3); expCellZ = ceil(vecZ / 3);
for k = 1:scalNumPartTot
    matFound = find(cellParticleList{expCellX(k), expCellY(k), expCellZ(k)} == k);
    assert(~isempty(matFound), ...
        sprintf('particle %d not in its expected cell (%d,%d,%d)', k, expCellX(k), expCellY(k), expCellZ(k)));
end
% wrapped positions stay in [0, L)
assert(all(vecXw >= 0) && all(vecXw < scalBox) && ...
       all(vecYw >= 0) && all(vecYw < scalBox) && ...
       all(vecZw >= 0) && all(vecZw < scalBox), 'wrapped positions out of [0, L)');
disp('Scenario 5 (rebuildCellList3D bucketing): PASSED — 4x4x4 cells, correct buckets');

disp('testFindNeighbors3D: ALL PASSED');
