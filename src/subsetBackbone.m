function [sPackCut, sFricCut] = subsetBackbone(sPack, sFric, vecBoolKeep)
% subsetBackbone  The backbone without some particles, fricState re-indexed onto it.
%
%   [sPackCut, sFricCut] = subsetBackbone(sPack, sFric, vecBoolKeep)
%
%   sPack is a backbone packing struct (vecPosX, vecPosY, vecDiameter, box
%   dims) and sFric its re-indexed fricState; vecBoolKeep marks the
%   particles that stay. Contacts that touch a removed particle are dropped,
%   and the remaining contacts are re-numbered onto the kept particles (as
%   pack.m's reindexFricState), so the pair can go straight to
%   matSpringDampMass. Used by runFrictionalDOS and runUnresolvedDOS.
    vecIdxNew = zeros(numel(vecBoolKeep), 1);
    vecIdxNew(vecBoolKeep) = 1:nnz(vecBoolKeep);
    sPackCut = sPack;
    sPackCut.vecPosX = sPack.vecPosX(vecBoolKeep);
    sPackCut.vecPosY = sPack.vecPosY(vecBoolKeep);
    sPackCut.vecDiameter = sPack.vecDiameter(vecBoolKeep);
    vecBoolLive = vecBoolKeep(sFric.vecContactNN) & vecBoolKeep(sFric.vecContactMM);
    sFricCut = sFric;
    for strField = ["vecOverlap", "vecFn", "vecFt", "vecUnitTanX", "vecUnitTanY"]
        sFricCut.(strField) = sFric.(strField)(vecBoolLive);
    end
    sFricCut.vecContactNN = vecIdxNew(sFric.vecContactNN(vecBoolLive));
    sFricCut.vecContactMM = vecIdxNew(sFric.vecContactMM(vecBoolLive));
    sFricCut.vecRadius = sFric.vecRadius(vecBoolKeep);
    sFricCut.vecInertia = sFric.vecInertia(vecBoolKeep);
    sFricCut.N = nnz(vecBoolKeep);
end
