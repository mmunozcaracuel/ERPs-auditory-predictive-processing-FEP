function time = erp_time_axis(data, matrixField, kind)
% Return and validate the S4-relative time axis of a study ERP matrix.
% Long epochs span a nominal 2000 ms, from -1000 to +1000 ms at 1000 Hz.
% The final sample may be +999 or +1000 ms, depending on endpoint inclusion.
% S1, S2, S3 and S4 occur at -900, -600, -300 and 0 ms, respectively.
if ~isfield(data, matrixField)
    error('ERPTime:MissingMatrix', 'Missing matrix: %s.', matrixField);
end
switch kind
    case 'long'
        timeField = 'timeLongS4'; startField = 'cnvPinvStartMs';
        endField = 'cnvPinvEndMs'; first = -1000; last = 1000;
    case 'pinv'
        timeField = 'timePinvS4'; startField = ''; endField = '';
        first = -1000; last = 1000;
    case 'mmn'
        timeField = 'timeMmnS4'; startField = 'mmnStartMs';
        endField = 'mmnEndMs'; first = -50; last = 300;
    otherwise
        error('ERPTime:Kind', 'Unknown epoch type: %s.', kind);
end
n = size(data.(matrixField), 2);
tol = 1e-6;
if isfield(data, 'samplingRateHz') && ...
        (~isscalar(data.samplingRateHz) || ~isfinite(data.samplingRateHz) || ...
        abs(data.samplingRateHz - 1000) > tol)
    error('ERPTime:Rate', 'The study matrices must be sampled at 1000 Hz.');
end
if isfield(data, timeField)
    time = double(data.(timeField)(:))';
else
    % Reconstruct only from explicit epoch bounds, never from matrix length alone.
    if isempty(startField) || ~all(isfield(data, {startField, endField})) || ...
            ~isequal(data.(startField), first) || ~isequal(data.(endField), last)
        error('ERPTime:Metadata', ...
            'Missing time metadata for %s. Export the axis from the processed EEG.', matrixField);
    end
    time = double(data.(startField)) + (0:n-1);
end
if numel(time) ~= n || isempty(time) || any(~isfinite(time)) || ...
        abs(time(1)-first) > tol || any(abs(diff(time)-1) > tol) || ...
        ~(abs(time(end)-(last-1)) < tol || abs(time(end)-last) < tol)
    error('ERPTime:Axis', ...
        'Incompatible time axis for %s: expected %g to %g ms at 1000 Hz (final endpoint optional).', ...
        matrixField, first, last);
end
end
