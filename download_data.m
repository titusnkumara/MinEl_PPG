function download_data(root)
%DOWNLOAD_DATA  Fetch the benchmark databases from Zenodo into the layout run_paper.m expects.
%
%   download_data          downloads into the parent of this file's folder,
%                          i.e. DATA = '.' as set in run_paper.m
%   download_data(root)    downloads into root
%
%   Source: Jayarathna T (2026) MinEl_PPG: Benchmarking Dataset, Zenodo,
%   https://doi.org/10.5281/zenodo.23721153
%
%   Files already present with the correct MD5 are skipped. Each download
%   shows a progress bar (via curl, which ships with macOS, Linux and
%   Windows 10+; without curl, websave is used and no bar is shown). Every
%   file is checked against its published MD5 after download. PPG-DaLiA is
%   unzipped into dalia/dalia_mat.
%
%   Resulting layout:
%     mimic-large-test/mimic_perform_large_test_a_data.mat
%     mimic_test/mimic_perform_test_all_data.mat
%     mimic_train/mimic_perform_train_all_data.mat
%     mimic_perform_af/mimic_perform_af_data.mat
%     mimic_perform_af/mimic_perform_non_af_data.mat
%     bidmc/bidmc_data.mat
%     capnobase/capnobase_db.mat
%     dalia/dalia_mat/...                (from dalia_mat.zip)

if nargin < 1 || isempty(root), root = fileparts(fileparts(mfilename('fullpath'))); end   % DATA = '../'
REC = '23721153';
F = { % file name                             folder               md5
    'mimic_perform_large_test_a_data.mat', 'mimic-large-test',  'f75a87885ae1dffc0a33d3ce6446abeb'
    'mimic_perform_test_all_data.mat',     'mimic_test',        '4cf66ff255ab369433ae2d733444edf6'
    'mimic_perform_train_all_data.mat',    'mimic_train',       '1003ac99d4a4c01fbbe8c54e3f981e43'
    'mimic_perform_af_data.mat',           'mimic_perform_af',  '38adc07e4420226a0edfe9769a23284f'
    'mimic_perform_non_af_data.mat',       'mimic_perform_af',  'd1c4f93442a3b4c9cb95718ca0a57013'
    'bidmc_data.mat',                      'bidmc',             'c0294ec2bb5494c756e48044e4ca93b5'
    'capnobase_db.mat',                    'capnobase',         '9d2c6b5c2764b5d1a1e8f35fdc610fc4'
    'dalia_mat.zip',                       'dalia',             '9d5b91adb421a4ff9101e2f47f903ff0' };

sizes = zenodo_sizes(REC);
[st, ~] = system('curl --version'); use_curl = st == 0;
if ~use_curl, fprintf('curl not found: downloading with websave (no progress bar).\n'); end
fprintf('\nDownloading to %s\n\n', root);

for k = 1:size(F, 1)
    [name, folder, md5] = F{k, :};
    dest = fullfile(root, folder, name);
    is_zip = strcmp(name, 'dalia_mat.zip');
    dalia_dir = fullfile(root, 'dalia', 'dalia_mat');
    if is_zip && isfile(fullfile(dalia_dir, '.complete'))
        fprintf('%-38s already unzipped, skipped\n', name); continue;
    end
    if isfile(dest) && strcmp(md5_file(dest), md5)
        fprintf('%-38s present, MD5 ok, skipped\n', name);
    else
        if ~isfolder(fileparts(dest)), mkdir(fileparts(dest)); end
        url = sprintf('https://zenodo.org/records/%s/files/%s?download=1', REC, name);
        total = 0; if isKey(sizes, name), total = sizes(name); end
        fetch(url, dest, name, total, use_curl);
        got = md5_file(dest);
        if ~strcmp(got, md5)
            delete(dest);
            error('download_data:md5', '%s: MD5 mismatch (got %s, expected %s). Run download_data again.', name, got, md5);
        end
        fprintf('%-38s MD5 ok\n', '');
    end
    if is_zip, unzip_dalia(dest, dalia_dir); end
end
fprintf('\nAll databases ready in %s\n', root);
end

% =========================================================================
function sizes = zenodo_sizes(rec)
% File sizes from the Zenodo API, for the progress bar (empty if unavailable).
sizes = containers.Map('KeyType', 'char', 'ValueType', 'double');
try
    r = webread(sprintf('https://zenodo.org/api/records/%s', rec), weboptions('Timeout', 30));
    f = r.files; if isstruct(f), f = num2cell(f); end
    for i = 1:numel(f), sizes(f{i}.key) = double(f{i}.size); end
catch
    fprintf('(Zenodo API not reachable: progress shown in MB only)\n');
end
end

% =========================================================================
function fetch(url, dest, name, total, use_curl)
% Download url to dest with a progress bar. Interrupted transfers are
% resumed (curl -C -) until the file reaches the size Zenodo reports, so a
% dropped connection costs only a retry; a partial .part file left by an
% earlier run is resumed too.
tmp = [dest '.part']; done = [dest '.done'];
if ~use_curl
    fprintf('%-38s downloading ...\n', name);
    websave(tmp, url, weboptions('Timeout', 600));
    movefile(tmp, dest); return;
end
MAX_TRY = 30; t0 = tic; last = ''; b0 = fsize(tmp);
for attempt = 1:MAX_TRY
    if total > 0 && fsize(tmp) >= total, break; end
    if isfile(done), delete(done); end
    if ispc
        cmd = sprintf('start "" /B cmd /c "curl -L -f -s -C - -o "%s" "%s" && (echo ok>"%s") || (echo fail>"%s")"', tmp, url, done, done);
    else
        cmd = sprintf('(curl -L -f -s -C - -o "%s" "%s" && echo ok > "%s" || echo fail > "%s") > /dev/null 2>&1 &', tmp, url, done, done);
    end
    system(cmd);
    b_prev = -1; t_prev = toc(t0);
    while true
        pause(0.25);
        b = fsize(tmp); el = toc(t0); rate = (b - b0) / max(el, 1e-3) / 2^20;
        if total > 0
            frac = min(1, b / total); n = round(30 * frac);
            line = sprintf('%-38s [%s%s] %5.1f%%  %6.1f/%6.1f MB  %5.2f MB/s', name, ...
                repmat('#', 1, n), repmat('.', 1, 30 - n), 100 * frac, b / 2^20, total / 2^20, rate);
        else
            line = sprintf('%-38s %6.1f MB  %5.2f MB/s', name, b / 2^20, rate);
        end
        if attempt > 1, line = sprintf('%s  (resume %d)', line, attempt - 1); end
        fprintf([repmat('\b', 1, numel(last)) '%s'], line); last = line;
        if isfile(done), break; end
        if b ~= b_prev, b_prev = b; t_prev = el;
        elseif el - t_prev > 120                       % stalled: abandon this attempt
            if ~ispc, system('pkill -f "curl -L -f -s -C - -o"'); end
            break;
        end
    end
    ok = isfile(done) && strcmp(strtrim(fileread(done)), 'ok');
    if isfile(done), delete(done); end
    if total == 0 && ok, break; end                     % size unknown: trust curl
    if total > 0 && fsize(tmp) >= total, break; end
    pause(2);                                          % dropped or stalled: resume
end
fprintf('\n');
if total > 0 && fsize(tmp) < total
    error('download_data:failed', ['%s: still incomplete after %d attempts (%.1f of %.1f MB). ' ...
        'Run download_data again; it resumes from where it stopped.'], name, MAX_TRY, fsize(tmp) / 2^20, total / 2^20);
end
movefile(tmp, dest);
end

function b = fsize(f)
d = dir(f); b = 0; if ~isempty(d), b = d.bytes; end
end

% =========================================================================
function unzip_dalia(zipfile, dalia_dir)
% Unzip PPG-DaLiA into dalia/dalia_mat, whether or not the zip has a top folder.
fprintf('%-38s unzipping ...\n', '');
tmpd = tempname; mkdir(tmpd);
unzip(zipfile, tmpd);
if isfolder(fullfile(tmpd, '__MACOSX')), rmdir(fullfile(tmpd, '__MACOSX'), 's'); end
top = dir(tmpd); top = top(~ismember({top.name}, {'.', '..'}));
if isfolder(dalia_dir), rmdir(dalia_dir, 's'); end
if numel(top) == 1 && top(1).isdir           % zip contains a single top folder
    movefile(fullfile(tmpd, top(1).name), dalia_dir);
else                                          % files at the zip root
    mkdir(dalia_dir); movefile(fullfile(tmpd, '*'), dalia_dir);
end
rmdir(tmpd, 's');
fid = fopen(fullfile(dalia_dir, '.complete'), 'w'); fclose(fid);
fprintf('%-38s unzipped to %s\n', '', dalia_dir);
end

% =========================================================================
function h = md5_file(f)
% MD5 of a file, read in 16 MB chunks.
md = java.security.MessageDigest.getInstance('MD5');
fid = fopen(f, 'r'); c = onCleanup(@() fclose(fid));
while true
    b = fread(fid, 2^24, '*uint8');
    if isempty(b), break; end
    md.update(typecast(b, 'int8'));
end
h = lower(sprintf('%02x', typecast(int8(md.digest()), 'uint8')));
end