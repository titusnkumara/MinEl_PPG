%EXAMPLE  Minimal use of minel_ppg on a synthetic PPG.
fs = 125;                                    % Hz
t  = (0:1/fs:60)';                           % one minute
hr = 72;                                     % bpm
phase = mod(t * hr/60, 1);
ppg = exp(-((phase - 0.20)/0.08).^2) + 0.35*exp(-((phase - 0.45)/0.06).^2);   % pulse + dicrotic wave
ppg = ppg + 0.02*randn(size(t));

beats = minel_ppg(ppg, fs);                  % native sample indices

fprintf('%d beats, mean heart rate %.1f bpm\n', numel(beats), 60/mean(diff(beats)/fs));
plot(t, ppg); hold on; plot(t(beats), ppg(beats), 'ro'); xlabel('Time (s)'); ylabel('PPG');
