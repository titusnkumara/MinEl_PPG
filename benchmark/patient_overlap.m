function keep = patient_overlap(files, large)
%PATIENT_OVERLAP  Held-out mask: Large Testing records from patients in no other dataset.
%   keep = patient_overlap(files, large) matches MIMIC-III subject IDs across
%   MIMIC PERform Large Testing (large = its loaded struct) and the Training,
%   Testing and AF / non-AF datasets (files.mimic_train, .mimic_test,
%   .mimic_af_af, .mimic_af_nonaf), and returns a logical mask over Large
%   Testing records: true = patient appears in no other dataset.
%   Called automatically by load_databases; prints the overlap counts.
path = find_subject_field(large.data(1).fix);
if isempty(path), error('No subject-ID field found in Large Testing data.fix.'); end

lt = arrayfun(@(d) string(get_path(d.fix, path)), large.data)';
id = @(f) arrayfun(@(d) string(field_or_path(d.fix, 'subj_id', path)), getfield(load(f), 'data'))';
tr = id(files.mimic_train); te = id(files.mimic_test);
af = [id(files.mimic_af_af); id(files.mimic_af_nonaf)];

ov = @(o) numel(intersect(unique(lt), unique(o)));
others = unique([tr; te; af]);
keep = ~ismember(lt, others);
fprintf('  [overlap] patients shared with Large Testing: Training %d | Testing %d | AF %d | any %d\n', ...
    ov(tr), ov(te), ov(af), ov(others));
fprintf('  [overlap] keeping %d of %d records (%d patients)\n', sum(keep), numel(keep), numel(unique(lt(keep))));
end

function field = find_subject_field(s)
% Dotted path to the first field named like a subject ID ('' if none).
field = '';
if ~isstruct(s), return; end
fn = fieldnames(s);
for k = 1:numel(fn)
    nm = lower(fn{k});
    if contains(nm, 'subject') || strcmp(nm, 'subj_id') || strcmp(nm, 'subjid')
        field = fn{k}; return;
    end
end
for k = 1:numel(fn)
    v = s.(fn{k});
    if isstruct(v) && isscalar(v)
        sub = find_subject_field(v);
        if ~isempty(sub), field = [fn{k} '.' sub]; return; end
    end
end
end

function v = get_path(s, path)
parts = strsplit(path, '.');
v = s;
for p = 1:numel(parts), v = v.(parts{p}); end
end

function v = field_or_path(s, field, path)
if isfield(s, field) && ~isempty(s.(field)), v = s.(field); else, v = get_path(s, path); end
end
