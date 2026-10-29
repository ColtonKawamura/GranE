function [vecBoolKeep, vecBoolHit] = unresolvedParticles(sFric, scalNi, scalForceMean, strRule, scalForceFloor)
% unresolvedParticles  Which backbone particles survive a removal rule.
%
%   [vecBoolKeep, vecBoolHit] = unresolvedParticles(sFric, scalNi, scalForceMean)
%   [vecBoolKeep, vecBoolHit] = unresolvedParticles(sFric, scalNi, scalForceMean, "twoContact")
%
%   sFric is pack.m's fricState re-indexed onto a backbone of scalNi
%   particles; scalForceMean is <f>, the mean |f_n| over all contacts of the
%   full pre-cleanRats state (the normalization of cleanRats and
%   plotUncertContacts).
%
%   Rules (strRule):
%     "unresolved"  (default) Agnolin & Roux, PRE 76, 061302 (2007): a
%                   particle whose backbone contacts ALL carry
%                   f_n < scalForceFloor*<f> is held only by forces inside
%                   pack.m's equilibration noise (force balance tolerance
%                   1e-4 <|f_n|>). It is removed, then cleanRats' z_i < 2
%                   rule is re-applied until nothing changes (the cascade).
%     "twoContact"  Shimada et al., arXiv:2608.30151 (2026) and Ishima et
%                   al., arXiv:2207.06632 (2023): every particle with
%                   z_i <= 2 is removed, iterated.
%
%   Returns
%     vecBoolKeep   [scalNi x 1] true for particles that stay
%     vecBoolHit    [scalNi x 1] true for particles the rule itself removes
%                   (the rest of ~vecBoolKeep is the cascade)
%
%   Used by runFrictionalDOS (spectra with and without the unresolved
%   particles) and runUnresolvedDOS (the paper figures).

    arguments
        sFric          (1, 1) struct
        scalNi         (1, 1) double {mustBeInteger, mustBeNonnegative}
        scalForceMean  (1, 1) double {mustBePositive}
        strRule        (1, 1) string {mustBeMember(strRule, ["unresolved", "twoContact"])} = "unresolved"
        scalForceFloor (1, 1) double {mustBePositive} = 1e-4
    end

    vecI = sFric.vecContactNN(:);
    vecJ = sFric.vecContactMM(:);
    vecZ = accumarray([vecI; vecJ], 1, [scalNi 1]);
    if strRule == "unresolved"
        vecForceRel = abs(sFric.vecFn(:))/scalForceMean;
        vecBoolWeak = vecForceRel < scalForceFloor;
        vecZWeak = accumarray([vecI(vecBoolWeak); vecJ(vecBoolWeak)], 1, [scalNi 1]);
        vecBoolHit = vecZ > 0 & vecZWeak == vecZ;
        scalZMin = 2;                                             % cleanRats' rattler rule
    else
        vecBoolHit = vecZ <= 2;
        scalZMin = 3;                                             % z_i <= 2 removed, iterated
    end
    vecBoolKeep = cascade(~vecBoolHit, vecI, vecJ, scalNi, scalZMin);
end

function vecBoolKeep = cascade(vecBoolKeep, vecI, vecJ, scalNi, scalZMin)
% Remove kept particles with fewer than scalZMin contacts among kept
% particles, repeatedly, until nothing changes.
    while true
        vecBoolLive = vecBoolKeep(vecI) & vecBoolKeep(vecJ);
        vecZ = zeros(scalNi, 1);
        if any(vecBoolLive)
            vecZ = accumarray([vecI(vecBoolLive); vecJ(vecBoolLive)], 1, [scalNi 1]);
        end
        vecBoolDrop = vecBoolKeep & vecZ < scalZMin;
        if ~any(vecBoolDrop)
            return;
        end
        vecBoolKeep(vecBoolDrop) = false;
    end
end
