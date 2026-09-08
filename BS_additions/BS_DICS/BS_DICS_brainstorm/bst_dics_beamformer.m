function [Filters, SourcePower, Orientations, Meta] = bst_dics_beamformer(Cx, Leadfield, varargin)
% BST_DICS_BEAMFORMER Compute scalar DICS filters without FieldTrip.
%
% Inputs:
%   Cx        - Hermitian sensor cross-spectral matrix [nChannels x nChannels]
%   Leadfield - Fixed [nChannels x nSources], free
%               [nChannels x 3*nSources], or [nChannels x 3 x nSources]
%
% Name-value options:
%   Regularization        - Fraction of mean sensor power (default: 0.05).
%                           Values >= 1 are interpreted as an absolute ridge.
%   Orientation           - 'max-power', 'fixed', or 'nai'.
%   NoiseCov              - Sensor noise covariance/CSD for 'nai'.
%   FixedOrient           - [3 x nSources] real unit orientations for 'fixed'.
%   LeadfieldFormat       - 'auto', 'fixed', or 'free'.
%   ReturnAllOrientations - Store vector filters in Meta (default: false).
%
% Outputs:
%   Filters      - Scalar unit-gain filters [nChannels x nSources]
%   SourcePower  - DICS power, or NAI, [nSources x 1]
%   Orientations - Selected orientations [3 x nSources]
%   Meta         - Solver diagnostics
%
% Author:
%   Vahab YoussofZadeh, 2026
%
% Notes:
%   Native Brainstorm implementation of scalar DICS beamforming.
%   Uses a common sensor CSD to estimate unit-gain spatial filters and
%   applies the filters to condition-specific CSD matrices.

% Parse and validate options
p = inputParser();
p.FunctionName = 'bst_dics_beamformer';
p.addRequired('Cx', @(x) isnumeric(x) && ismatrix(x));
p.addRequired('Leadfield', @(x) isnumeric(x) && (ismatrix(x) || ndims(x) == 3));
p.addParameter('Regularization', 0.05, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
p.addParameter('Orientation', 'max-power', @(x) ischar(x) || isstring(x));
p.addParameter('NoiseCov', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('FixedOrient', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('LeadfieldFormat', 'auto', @(x) ischar(x) || isstring(x));
p.addParameter('ReturnAllOrientations', false, @(x) islogical(x) && isscalar(x));
p.parse(Cx, Leadfield, varargin{:});
opts = p.Results;

if size(Cx, 1) ~= size(Cx, 2)
    error('Cx must be square.');
end
if any(~isfinite(Cx(:)))
    error('Cx contains NaN or Inf values.');
end

orientationMode = lower(char(opts.Orientation));
if ~ismember(orientationMode, {'max-power', 'fixed', 'nai'})
    error('Unknown orientation mode: %s.', orientationMode);
end
leadfieldFormat = lower(char(opts.LeadfieldFormat));
if ~ismember(leadfieldFormat, {'auto', 'fixed', 'free'})
    error('LeadfieldFormat must be auto, fixed, or free.');
end

% Hermitian sensor matrices
Cx = (Cx + Cx') ./ 2;
nChannels = size(Cx, 1);
[Leadfield3D, nSources, nOrient] = local_parse_leadfield(Leadfield, leadfieldFormat);
if size(Leadfield3D, 1) ~= nChannels
    error('Leadfield has %d channels, but Cx has %d.', size(Leadfield3D, 1), nChannels);
end
if any(~isfinite(Leadfield3D(:)))
    error('Leadfield contains NaN or Inf values.');
end

if strcmp(orientationMode, 'nai')
    if isempty(opts.NoiseCov)
        error('NoiseCov is required for NAI orientation.');
    end
    if ~isequal(size(opts.NoiseCov), size(Cx)) || any(~isfinite(opts.NoiseCov(:)))
        error('NoiseCov must be finite and the same size as Cx.');
    end
    Cn = (opts.NoiseCov + opts.NoiseCov') ./ 2;
else
    Cn = [];
end

if strcmp(orientationMode, 'fixed') && nOrient == 3
    if ~isequal(size(opts.FixedOrient), [3, nSources])
        error('FixedOrient must be [3 x nSources] for a free-orientation leadfield.');
    end
    opts.FixedOrient = real(opts.FixedOrient);
    orientNorm = sqrt(sum(opts.FixedOrient .^ 2, 1));
    if any(orientNorm <= 0)
        error('FixedOrient contains one or more zero vectors.');
    end
    opts.FixedOrient = bsxfun(@rdivide, opts.FixedOrient, orientNorm);
end

% Regularize only in sensor space
sensorScale = real(trace(Cx)) ./ nChannels;
if ~isfinite(sensorScale) || sensorScale <= 0
    sensorScale = mean(abs(diag(Cx)));
end
if opts.Regularization < 1
    regAbs = opts.Regularization .* sensorScale;
else
    regAbs = opts.Regularization;
end
Creg = Cx + regAbs .* eye(nChannels, 'like', Cx);
Ci = pinv(Creg);
Ci = (Ci + Ci') ./ 2;

Filters = complex(zeros(nChannels, nSources, 'like', Cx));
SourcePower = zeros(nSources, 1, 'like', real(Cx));
Orientations = zeros(3, nSources, 'like', real(Cx));
if opts.ReturnAllOrientations
    allOrientFilters = complex(zeros(nChannels, nOrient, nSources, 'like', Cx));
else
    allOrientFilters = [];
end

for iSource = 1:nSources
    Li = Leadfield3D(:, :, iSource);
    Gi = Li' * Ci * Li;
    Gi = (Gi + Gi') ./ 2;

    if opts.ReturnAllOrientations
        allOrientFilters(:, :, iSource) = Ci * Li * pinv(Gi);
    end

    if nOrient == 1
        orient = 1;
        Orientations(:, iSource) = [1; 0; 0];
    else
        switch orientationMode
            case 'fixed'
                orient = opts.FixedOrient(:, iSource);

            case 'max-power'
                % For a real dipole orientation, maximum unit-gain power is
                % obtained from the smallest eigenvalue of Re(L' C^-1 L).
                orientMetric = real(Gi);
                orientMetric = (orientMetric + orientMetric') ./ 2;
                [vectors, values] = eig(orientMetric);
                [~, iBest] = min(real(diag(values)));
                orient = real(vectors(:, iBest));

            case 'nai'
                % The scalar-filter denominator cancels in this generalized
                % source-to-noise ratio, leaving a stable orientation problem.
                signalMetric = Li' * Ci' * Cx * Ci * Li;
                noiseMetric = Li' * Ci' * Cn * Ci * Li;
                signalMetric = real((signalMetric + signalMetric') ./ 2);
                noiseMetric = real((noiseMetric + noiseMetric') ./ 2);
                noiseScale = max(abs(diag(noiseMetric)));
                if isempty(noiseScale) || ~isfinite(noiseScale) || noiseScale <= 0
                    noiseScale = 1;
                end
                noiseMetric = noiseMetric + eps(class(noiseMetric)) .* noiseScale .* eye(nOrient, 'like', noiseMetric);
                [vectors, values] = eig(signalMetric, noiseMetric);
                scores = real(diag(values));
                scores(~isfinite(scores)) = -Inf;
                [~, iBest] = max(scores);
                orient = real(vectors(:, iBest));
        end

        orientNorm = norm(orient);
        if ~isfinite(orientNorm) || orientNorm <= eps(class(orientNorm))
            Filters(:, iSource) = NaN;
            SourcePower(iSource) = NaN;
            Orientations(:, iSource) = NaN;
            continue;
        end
        orient = orient ./ orientNorm;
        Orientations(:, iSource) = orient;
    end

    % Build one scalar, unit-gain beamformer for the selected orientation.
    projectedLeadfield = Li * orient;
    unitGainDenom = real(projectedLeadfield' * Ci * projectedLeadfield);
    denomTol = eps(class(unitGainDenom)) .* max(1, norm(projectedLeadfield) .* norm(Ci * projectedLeadfield));
    if ~isfinite(unitGainDenom) || unitGainDenom <= denomTol
        Filters(:, iSource) = NaN;
        SourcePower(iSource) = NaN;
        continue;
    end

    filter = (Ci * projectedLeadfield) ./ unitGainDenom;
    response = filter' * projectedLeadfield;
    if abs(response) > denomTol
        filter = filter ./ conj(response);
    end
    Filters(:, iSource) = filter;

    signalPower = real(filter' * Cx * filter);
    if strcmp(orientationMode, 'nai')
        noisePower = real(filter' * Cn * filter);
        if noisePower > 0
            SourcePower(iSource) = signalPower ./ noisePower;
        else
            SourcePower(iSource) = NaN;
        end
    else
        SourcePower(iSource) = signalPower;
    end
end

Meta = struct();
Meta.nChannels = nChannels;
Meta.nSources = nSources;
Meta.nOrientations = nOrient;
Meta.Regularization = opts.Regularization;
Meta.RegularizationAbs = regAbs;
Meta.OrientationMode = orientationMode;
Meta.LeadfieldFormat = leadfieldFormat;
Meta.NoiseAware = strcmp(orientationMode, 'nai');
Meta.UnitGain = true;
if opts.ReturnAllOrientations
    Meta.AllOrientFilters = allOrientFilters;
end
end


function [leadfield3D, nSources, nOrient] = local_parse_leadfield(leadfieldIn, format)
if ndims(leadfieldIn) == 3
    [nChannels, nOrient, nSources] = size(leadfieldIn);
    if ~ismember(nOrient, [1, 3])
        error('A 3-D leadfield must be [nChannels x 1|3 x nSources].');
    end
    leadfield3D = reshape(leadfieldIn, [nChannels, nOrient, nSources]);
    return;
end

[nChannels, nColumns] = size(leadfieldIn);
if strcmp(format, 'fixed')
    nOrient = 1;
    nSources = nColumns;
    leadfield3D = reshape(leadfieldIn, [nChannels, 1, nSources]);
elseif strcmp(format, 'free') || (strcmp(format, 'auto') && mod(nColumns, 3) == 0)
    if mod(nColumns, 3) ~= 0
        error('A free-orientation leadfield must have 3*nSources columns.');
    end
    nOrient = 3;
    nSources = nColumns ./ 3;
    leadfield3D = zeros(nChannels, 3, nSources, 'like', leadfieldIn);
    for iSource = 1:nSources
        columns = (iSource - 1) * 3 + (1:3);
        leadfield3D(:, :, iSource) = leadfieldIn(:, columns);
    end
else
    nOrient = 1;
    nSources = nColumns;
    leadfield3D = reshape(leadfieldIn, [nChannels, 1, nSources]);
end
end
