% findNeighbors3D -- Build the candidate contact-pair list for a 3D cell list.
% For every particle in every cell, candidate neighbors are the particles in
% the 27 cells of the periodic 3x3x3 stencil around that cell. Each (source,
% dest) pair is emitted once, with source < dest.
%
% Extracted from pack.m (was a local subfunction) so it can be unit-tested
% directly with a synthetic cellParticleList.

function [vecPairIdxSource, vecPairIdxDest, scalNumPairs, scalMaxPairs] = findNeighbors3D( ...
        cellParticleList, scalNumCellsX, scalNumCellsY, scalNumCellsZ, ...
        vecPairIdxSource, vecPairIdxDest, scalMaxPairs)

    scalNumPairs = 0;

    % go throgh each combination of cells
    for idxCellX = 1:scalNumCellsX
        for idxCellY = 1:scalNumCellsY
            for idxCellZ = 1:scalNumCellsZ

                % ----------------------------- START WRONG ----------------------------
                % % wrap the cells
                % % MUDSUCK : if there is only 2 Cells per axis, then the same cell will appear twince
                % % after wraped.
                % % example: cellIndex = [1,2], so 2's neighbors  will be [1,2(this cell), 1]
                % % which will double -count neighbors for each particles.
                % scalCellLeft  = mod(idxCellX-2, scalNumCellsX)+1;
                % scalCellRight = mod(idxCellX,   scalNumCellsX)+1;
                % scalCellDown  = mod(idxCellY-2, scalNumCellsY)+1;
                % scalCellUp    = mod(idxCellY,   scalNumCellsY)+1;
                % scalCellBack  = mod(idxCellZ-2, scalNumCellsZ)+1;
                % scalCellFront = mod(idxCellZ,   scalNumCellsZ)+1;
                %
                % % define "this cell" in this loop
                % % vector of  all particles index in this cell
                % vecCurrentCellPartIdx = cellParticleList{idxCellX, idxCellY, idxCellZ};
                %
                % % (at this point) a vector that contains the particle index of every 
                % % particle that COULD be a neighbor beacuse it's in a cell that is 
                % % adjacent "this" cell 
                % vecNeighborList = [ ...
                %     cellParticleList{scalCellLeft,  scalCellDown, scalCellBack};  ...
                %     cellParticleList{scalCellLeft,  scalCellDown, idxCellZ};      ...
                %     cellParticleList{scalCellLeft,  scalCellDown, scalCellFront}; ...
                %     cellParticleList{scalCellLeft,  idxCellY,     scalCellBack};  ...
                %     cellParticleList{scalCellLeft,  idxCellY,     idxCellZ};      ...
                %     cellParticleList{scalCellLeft,  idxCellY,     scalCellFront}; ...
                %     cellParticleList{scalCellLeft,  scalCellUp,   scalCellBack};  ...
                %     cellParticleList{scalCellLeft,  scalCellUp,   idxCellZ};      ...
                %     cellParticleList{scalCellLeft,  scalCellUp,   scalCellFront}; ...
                %     cellParticleList{idxCellX,      scalCellDown, scalCellBack};  ...
                %     cellParticleList{idxCellX,      scalCellDown, idxCellZ};      ...
                %     cellParticleList{idxCellX,      scalCellDown, scalCellFront}; ...
                %     cellParticleList{idxCellX,      idxCellY,     scalCellBack};  ...
                %     vecCurrentCellPartIdx;                                                ...
                %     cellParticleList{idxCellX,      idxCellY,     scalCellFront}; ...
                %     cellParticleList{idxCellX,      scalCellUp,   scalCellBack};  ...
                %     cellParticleList{idxCellX,      scalCellUp,   idxCellZ};      ...
                %     cellParticleList{idxCellX,      scalCellUp,   scalCellFront}; ...
                %     cellParticleList{scalCellRight, scalCellDown, scalCellBack};  ...
                %     cellParticleList{scalCellRight, scalCellDown, idxCellZ};      ...
                %     cellParticleList{scalCellRight, scalCellDown, scalCellFront}; ...
                %     cellParticleList{scalCellRight, idxCellY,     scalCellBack};  ...
                %     cellParticleList{scalCellRight, idxCellY,     idxCellZ};      ...
                %     cellParticleList{scalCellRight, idxCellY,     scalCellFront}; ...
                %     cellParticleList{scalCellRight, scalCellUp,   scalCellBack};  ...
                %     cellParticleList{scalCellRight, scalCellUp,   idxCellZ};      ...
                %     cellParticleList{scalCellRight, scalCellUp,   scalCellFront}];
                % ----------------------------- END WRONG ----------------------------

                % --- Build unique neighbor-cell list (3D periodic) ---
                % 3x3x3 stencil: all (scalIdxOffsetX,scalIdxOffsetY,scalIdxOffsetZ) in {-1,0,1} around this cell.
                % Use mod(...) to wrap, then unique(...,'rows') to remove
                % duplicate (ix,iy,iz) triplets that occur when scalNumCells
                % is small (e.g. 2 per axis).
                matIdxNeighborCells = zeros(27, 3); % each row is a neighbor cell [ix, iy,  iz] 
                idxNeighbor = 0; % start at zero, and increment each loop. After [-1, 0, +1] times will be 27 (all cells inclusive)
                for scalIdxOffsetX = -1:1
                    for scalIdxOffsetY = -1:1
                        for scalIdxOffsetZ = -1:1
                            idxNeighbor = idxNeighbor + 1; % increments the row. each row = neighbor cell
                            ix = mod(idxCellX - 1 + scalIdxOffsetX, scalNumCellsX) + 1;
                            iy = mod(idxCellY - 1 + scalIdxOffsetY, scalNumCellsY) + 1;
                            iz = mod(idxCellZ - 1 + scalIdxOffsetZ, scalNumCellsZ) + 1;
                            matIdxNeighborCells(idxNeighbor, :) = [ix, iy, iz];
                        end
                    end
                end

                % Remove duplicate neighbor cells (same [ix,iy,iz])
                matIdxNeighborCells = unique(matIdxNeighborCells, 'rows');

                % Particles in the current cell
                vecCurrentCellPartIdx = cellParticleList{idxCellX, idxCellY, idxCellZ};

                % Build neighbor particle (not nessearily contacts) list by concatenating each unique cell
                vecNeighborList = [];
                for k = 1:size(matIdxNeighborCells, 1)
                    ix = matIdxNeighborCells(k, 1);
                    iy = matIdxNeighborCells(k, 2);
                    iz = matIdxNeighborCells(k, 3);
                    vecNeighborList = [vecNeighborList; cellParticleList{ix, iy, iz}];
                end

                % go through every particle index in this cell
                for idxNN = vecCurrentCellPartIdx'

                    % pick out only the index of particles that are greater than this one
                    % prevents from counting contact twince
                    vecContactCandidatePartIdx = vecNeighborList(vecNeighborList > idxNN);
                    scalNumCandidates = numel(vecContactCandidatePartIdx);
                    if scalNumCandidates == 0; continue; end

                    % expand in case we undercounted contacts (see intialization of
                    % these vectors at the begging of the file
                    % I don't think this should ever happen
                    if scalNumPairs + scalNumCandidates > scalMaxPairs
                        scalMaxPairs = 2 * scalMaxPairs;
                        vecPairIdxSource(scalMaxPairs) = 0;
                        vecPairIdxDest(scalMaxPairs) = 0;
                    end

                    % for each pair, set the source particle index
                    % example: 
                    % vecPairIdxSource(1:3) = [5; 5; 5];
                    vecPairIdxSource(scalNumPairs+1 : scalNumPairs+scalNumCandidates) = idxNN;

                    % set the destination:
                    % example:
                    % vecPairIdxDest(1:3)   = [9; 10; 12];
                    vecPairIdxDest(scalNumPairs+1 : scalNumPairs+scalNumCandidates) = vecContactCandidatePartIdx;

                    % example contacts: (5,9), (5,10,) (5, 12)


                    scalNumPairs = scalNumPairs + scalNumCandidates;
                end
            end
        end
    end
end
