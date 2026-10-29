function y = comparators(action, x, fs)
%COMPARATORS  MSPTDfast run exactly as the ppg-beats benchmark runs it.
%   Requires the ppg-beats toolbox on the MATLAB path (see README).
%
%   beats = comparators('msptd_pb', ppg, fs)  pre-filter + MSPTDfast + tidying
%                                             + mid-amplitude fiducial (MSPTD_PB)
%   beats = comparators('msptd',    ppg, fs)  the same without the pre-filter
%   ppg_f = comparators('filter',   ppg, fs)  the ppg-beats pre-filter only
%
%   The pre-filter reproduces ppg-beats preprocess_ppg_signal.m: a Kaiser-window
%   FIR low-pass (stop 6.93 Hz, pass 10 Hz), then baseline removal by
%   subtracting a low-pass (stop 0.6 Hz, pass 1.02 Hz) of the Tukey-windowed
%   signal, both zero-phase. Filter design is cached per sampling rate and so
%   excluded from timing after the first (warm-up) call, as in the paper.
switch action
    case 'filter'
        y = charlton_filter(x, fs);
    case 'msptd'
        y = msptd_core(x, fs);
    case 'msptd_pb'
        y = msptd_core(charlton_filter(x, fs), fs);
    otherwise
        error('comparators: unknown action "%s"', action);
end
end


% =========================================================================
function mid = msptd_core(x, fs)
% MSPTDfast v2 via ppg-beats detect_ppg_beats: detector, onset-peak tidying,
% mid-amplitude points (the fiducial the benchmark scores).
S.v = x(:); S.fs = fs;
[~, ~, mid] = detect_ppg_beats(S, 'MSPTDfastv2');
mid = unique(mid(~isnan(mid)));
end

% =========================================================================
function y = charlton_filter(x, fs)
[b_high, b_low] = charlton_coeffs(fs);
x = x(:); n = numel(x); y = x;
if n > 3*(length(b_high)-1)
    y = filtfilt(b_high, 1, y);                    % remove > ~8 Hz
end
if n > 3*(length(b_low)-1)
    prop_win = 0.5*2 / ((n-1)/fs);                 % 0.5 s Tukey edges
    y = y - filtfilt(b_low, 1, y .* tukeywin(n, prop_win));   % remove < ~0.67 Hz
end
end

% =========================================================================
function [b_high, b_low] = charlton_coeffs(fs)
persistent cache
if isempty(cache), cache = containers.Map('KeyType', 'double', 'ValueType', 'any'); end
if ~isKey(cache, fs)
    cache(fs) = {kaiser_lpf(6.93, 10, fs), kaiser_lpf(0.6, 1.02, fs)};
end
c = cache(fs); b_high = c{1}; b_low = c{2};
end

function b = kaiser_lpf(Fstop, Fpass, fs)
% ppg-beats create_lpf(): Kaiser-window FIR, deviations 0.01 (stop), 0.05 (pass).
[N, Wn, BETA, TYPE] = kaiserord([Fstop Fpass]/(fs/2), [1 0], [0.01 0.05]);
b = fir1(N, Wn, TYPE, kaiser(N+1, BETA), 'scale');
end
