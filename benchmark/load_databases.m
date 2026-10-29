function [subj, Fs, is_ecg_ref] = load_databases(DB, cfg)
%LOAD_DATABASES  Load one benchmark database into a common format.
%   [subj, Fs, is_ecg_ref] = load_databases(DB, cfg)
%
%   DB   : 'mimic_large_test' | 'mimic_test' | 'mimic_train' | 'mimic_af' |
%          'bidmc' | 'capnobase' | 'dalia'
%   subj : struct array, one element per recording, with fields
%            ppg    raw PPG
%            ecg    raw ECG (ECG-referenced databases) or reference beat
%                   indices (CapnoBase PPG annotations, PPG-DaLiA R peaks)
%            group  0 none | 1 adult / non-AF | 2 neonate / AF
%            pid    patient ID (held-out set only; empty elsewhere = one
%                   recording per patient)
%            act    activity labels at 4 Hz (PPG-DaLiA only)
%   File locations are set in run_paper.m (cfg.files).
switch DB
    case 'mimic_af'
        S1 = load(cfg.files.mimic_af_nonaf);
        S2 = load(cfg.files.mimic_af_af);
        subj = [from_ppg_beats(S1.data, 1); from_ppg_beats(S2.data, 2)];
        Fs = 125; is_ecg_ref = true;

    case 'mimic_test'
        S = load(cfg.files.mimic_test);
        subj = from_ppg_beats(S.data, []);
        Fs = 125; is_ecg_ref = true;

    case 'mimic_train'
        S = load(cfg.files.mimic_train);
        subj = from_ppg_beats(S.data, []);
        Fs = 125; is_ecg_ref = true;

    case 'bidmc'
        S = load(cfg.files.bidmc);
        subj = from_ppg_beats(S.data, []);
        Fs = 125; is_ecg_ref = true;

    case 'mimic_large_test'
        % Held-out set. Records from patients who appear in any other MIMIC
        % PERform dataset are removed (patient_overlap.m; 796 of 1,257
        % recordings kept, 111 patients).
        S = load(cfg.files.mimic_large_test);
        keep = patient_overlap(cfg.files, S);
        idx = find(keep);
        subj = from_ppg_beats(S.data(idx), 0);
        for j = 1:numel(idx)
            subj(j).pid = string(S.data(idx(j)).fix.mimic_details.subj_id);  % MIMIC-III subject ID
        end
        Fs = S.data(1).ppg.fs; is_ecg_ref = true;

    case 'capnobase'
        S = load(cfg.files.capnobase);
        D = S.data;
        valid = find(~cellfun(@isempty, {D.ref_peaks}));
        if isempty(valid), error('No CapnoBase records with reference peaks.'); end
        Fs = D(valid(1)).fs;
        subj = struct('group', {D(valid).group}', 'ecg', {D(valid).ref_peaks}', ...
                      'ppg', {D(valid).ppg}', 'pid', [], 'act', []);
        is_ecg_ref = false;

    case 'dalia'
        % 64 Hz wrist PPG with ECG R peaks on the PPG timeline and activity labels.
        subj = struct('group', {}, 'ecg', {}, 'ppg', {}, 'pid', {}, 'act', {});
        for s = cfg.dalia.subjects
            f = fullfile(cfg.dalia.data_dir, sprintf('S%d.mat', s));
            if ~exist(f, 'file'), fprintf('missing %s\n', f); continue; end
            D = load(f);
            p = double(D.ppg(:)); p(~isfinite(p)) = 0;
            r = round(sort(double(D.rpeaks_ppg(:))));
            r = unique(r(r >= 1 & r <= numel(p)));
            subj(end+1, 1) = struct('group', 0, 'ecg', r, 'ppg', p, ...
                                    'pid', [], 'act', double(D.activity(:))); %#ok<AGROW>
        end
        Fs = 64; is_ecg_ref = false;

    otherwise
        error('Unknown database: %s', DB);
end
end


% =========================================================================
function subj = from_ppg_beats(data, group)
% ppg-beats MIMIC PERform / BIDMC format. If group is empty, the adult /
% neonate label is read from data(i).fix.group ('a' -> 1, 'n' -> 2).
n = numel(data);
subj = struct('group', cell(n,1), 'ecg', cell(n,1), 'ppg', cell(n,1), ...
              'pid', cell(n,1), 'act', cell(n,1));
for i = 1:n
    if isempty(group)
        g = 0;
        if isfield(data(i), 'fix') && isstruct(data(i).fix) && isfield(data(i).fix, 'group')
            c = lower(char(data(i).fix.group));
            if ~isempty(c) && c(1) == 'a', g = 1; elseif ~isempty(c) && c(1) == 'n', g = 2; end
        end
        subj(i).group = g;
    else
        subj(i).group = group;
    end
    subj(i).ecg = data(i).ekg.v(:);
    subj(i).ppg = data(i).ppg.v(:);
end
end
