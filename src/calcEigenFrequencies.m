function sSpectrum = calcEigenFrequencies(matSpring, matMass, options)
% calcEigenFrequencies  Undamped vibrational frequencies of a packing.
%
%   sSpectrum = calcEigenFrequencies(matSpring, matMass)
%   sSpectrum = calcEigenFrequencies(matSpring, matMass, "zeroTol", 1e-9)
%
%   With no damping the equation of motion M x'' + K x = 0 has solutions
%   x ~ exp(+-i omega t) with K v = omega^2 M v. M is diagonal (masses and
%   moments of inertia), so v = M^(-1/2) w turns this into the symmetric
%   problem
%       D w = omega^2 w,   D = M^(-1/2) K M^(-1/2)   (the dynamical matrix),
%   whose eigenvalues omega^2 are real and are computed here with eig.
%   This is the matrix whose spectrum Henkes, van Hecke & van Saarloos,
%   EPL 90, 14003 (2010) histogram to get the density of states.
%
%   Arguments:
%       matSpring   stiffness matrix K (dense or sparse), e.g. from
%                   matSpringDampMass
%       matMass     diagonal mass matrix M from matSpringDampMass
%
%   Options:
%       zeroTol     |omega^2| <= zeroTol * max(omega^2) counts as a zero mode
%                   (default 1e-9)
%
%   Returns struct sSpectrum with fields:
%       vecOmegaSq     [dof x 1] sorted eigenvalues omega^2 of D
%       vecOmega       [dof x 1] sqrt(omega^2) of the positive modes, 0 for
%                      zero modes, NaN for unstable (omega^2 < 0) modes
%       scalNumModes   dof, the total number of modes
%       scalNumZero    number of zero modes (|omega^2| <= tol)
%       scalNumNegative number of unstable modes (omega^2 < -tol)
%
%   Example:
%       sPack = load('2D_N500_P0.001_Width22_Seed1_Friction.mat');   % pack.m output
%       [matK, ~, matM] = matSpringDampMass(sPack, 0, sPack.K, 1, "periodic", true, ...
%           "uniformMass", true, "friction", true);
%       sSpectrum = calcEigenFrequencies(matK, matM);

    arguments
        matSpring (:, :) double
        matMass   (:, :) double
        options.zeroTol (1, 1) double {mustBePositive} = 1e-9
    end

    vecMassDiag = full(diag(matMass));
    vecInvSqrtMass = 1./sqrt(vecMassDiag);
    matDynamical = full(matSpring).*(vecInvSqrtMass*vecInvSqrtMass.');
    matDynamical = (matDynamical + matDynamical.')/2;   % remove round-off asymmetry

    vecOmegaSq = sort(eig(matDynamical));
    scalTol = options.zeroTol*max(abs(vecOmegaSq));
    vecBoolZero = abs(vecOmegaSq) <= scalTol;
    vecBoolNegative = vecOmegaSq < -scalTol;

    vecOmega = sqrt(max(vecOmegaSq, 0));
    vecOmega(vecBoolZero) = 0;
    vecOmega(vecBoolNegative) = NaN;

    sSpectrum = struct();
    sSpectrum.vecOmegaSq = vecOmegaSq;
    sSpectrum.vecOmega = vecOmega;
    sSpectrum.scalNumModes = numel(vecOmegaSq);
    sSpectrum.scalNumZero = sum(vecBoolZero);
    sSpectrum.scalNumNegative = sum(vecBoolNegative);
end
