function [qrs_amp_raw,qrs_i_raw,delay]=pan_tompkin(ecg,fs,gr)
%PAN_TOMPKIN  Pan-Tompkins QRS detector (Pan and Tompkins 1985).
%% ============== Licensce ========================================== %%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
% FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
% OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
% TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% Author :
% Hooman Sedghamiz, Feb, 2018
% MSc. Biomedical Engineering, Linkoping University
% Email : Hooman.sedghamiz@gmail.com
%% ============ Update History ================== %%
% Feb 2018 : 
%           1- Cleaned up the code and added more comments
%           2- Added to BioSigKit Toolbox
%% ================= Now Part of BioSigKit ==================== %%
%
%   Local changes: for fs <= 30 Hz the band-pass is skipped and
%   findpeaks_minpeakdistance_fast replaces findpeaks. Neither branch is
%   used in the benchmark (all ECGs are >= 64 Hz).
if ~isvector(ecg)
  error('ecg must be a row or column vector');
end
if nargin < 3
    gr = 1;
end
ecg = ecg(:);

delay = 0;
skip = 0;
m_selected_RR = 0;
mean_RR = 0;
ser_back = 0; 
ax = zeros(1,6);

if fs == 200
  ecg = ecg - mean(ecg);
  Wn = 12*2/fs;
  N = 3;
  [a,b] = butter(N,Wn,'low');
  ecg_l = filtfilt(a,b,ecg); 
  ecg_l = ecg_l/ max(abs(ecg_l));
  if gr
    figure;
    ax(1) = subplot(321);plot(ecg);axis tight;title('Raw signal');
    ax(2)=subplot(322);plot(ecg_l);axis tight;title('Low pass filtered');
  end
  Wn = 5*2/fs;
  N = 3;
  [a,b] = butter(N,Wn,'high');
  ecg_h = filtfilt(a,b,ecg_l); 
  ecg_h = ecg_h/ max(abs(ecg_h));
  if gr
    ax(3)=subplot(323);plot(ecg_h);axis tight;title('High Pass Filtered');
  end
else
    if(fs<=30)
            ecg_h = ecg;
    else
            f1=5; f2=15;
          Wn=[f1 f2]*2/fs;
          N = 3;
          [a,b] = butter(N,Wn);
          ecg_h = filtfilt(a,b,ecg);
          ecg_h = ecg_h/ max( abs(ecg_h));
          if gr
            ax(1) = subplot(3,2,[1 2]);plot(ecg);axis tight;title('Raw Signal');
            ax(3)=subplot(323);plot(ecg_h);axis tight;title('Band Pass Filtered');
          end
    end  
end

if fs ~= 200
  int_c = (5-1)/(fs*1/40);
  b = interp1(1:5,[1 2 0 -2 -1].*(1/8)*fs,1:int_c:5);
else
  b = [1 2 0 -2 -1].*(1/8)*fs;   
end
ecg_d = filtfilt(b,1,ecg_h);
ecg_d = ecg_d/max(ecg_d);
if gr
  ax(4)=subplot(324);plot(ecg_d); axis tight; title('Derivative');
end
ecg_s = ecg_d.^2;
if gr
  ax(5)=subplot(325); plot(ecg_s); axis tight; title('Squared');
end
ecg_m = conv(ecg_s ,ones(1 ,round(0.150*fs))/round(0.150*fs));
delay = delay + round(0.150*fs)/2;
if gr
  ax(6)=subplot(326);plot(ecg_m); axis tight; title('Moving average');
end

if(fs<=30)
    [pks,locs] = findpeaks_minpeakdistance_fast(ecg_m,round(0.2*fs));
else
    [pks,locs] = findpeaks(ecg_m,'MINPEAKDISTANCE',round(0.2*fs));
end

LLp = length(pks);
qrs_c = zeros(1,LLp);
qrs_i = zeros(1,LLp);
qrs_i_raw = zeros(1,LLp);
qrs_amp_raw= zeros(1,LLp);
nois_c = zeros(1,LLp);
nois_i = zeros(1,LLp);
SIGL_buf = zeros(1,LLp);
NOISL_buf = zeros(1,LLp);
SIGL_buf1 = zeros(1,LLp);
NOISL_buf1 = zeros(1,LLp);
THRS_buf1 = zeros(1,LLp);
THRS_buf = zeros(1,LLp);

THR_SIG = max(ecg_m(1:2*fs))*1/3;
THR_NOISE = mean(ecg_m(1:2*fs))*1/2;
SIG_LEV= THR_SIG;
NOISE_LEV = THR_NOISE;
THR_SIG1 = max(ecg_h(1:2*fs))*1/3;
THR_NOISE1 = mean(ecg_h(1:2*fs))*1/2;
SIG_LEV1 = THR_SIG1;
NOISE_LEV1 = THR_NOISE1;
Beat_C = 0; Beat_C1 = 0; Noise_Count = 0;

for i = 1 : LLp  
    if locs(i)-round(0.150*fs)>= 1 && locs(i)<= length(ecg_h)
        [y_i,x_i] = max(ecg_h(locs(i)-round(0.150*fs):locs(i)));
    else
        if i == 1
            [y_i,x_i] = max(ecg_h(1:locs(i)));
            ser_back = 1;
        elseif locs(i)>= length(ecg_h)
            [y_i,x_i] = max(ecg_h(locs(i)-round(0.150*fs):end));
        end       
    end       
    if Beat_C >= 9        
        diffRR = diff(qrs_i(Beat_C-8:Beat_C));
        mean_RR = mean(diffRR);
        comp = qrs_i(Beat_C)-qrs_i(Beat_C-1);
        if comp <= 0.92*mean_RR || comp >= 1.16*mean_RR
            THR_SIG = 0.5*(THR_SIG);
            THR_SIG1 = 0.5*(THR_SIG1);               
        else
            m_selected_RR = mean_RR;
        end          
    end
    if m_selected_RR
        test_m = m_selected_RR;
    elseif mean_RR && m_selected_RR == 0
        test_m = mean_RR;
    else
        test_m = 0;
    end
    if test_m
        if (locs(i) - qrs_i(Beat_C)) >= round(1.66*test_m)
            [pks_temp,locs_temp] = max(ecg_m(qrs_i(Beat_C)+ round(0.200*fs):locs(i)-round(0.200*fs)));
            locs_temp = qrs_i(Beat_C)+ round(0.200*fs) + locs_temp -1;
            if pks_temp > THR_NOISE
                Beat_C = Beat_C + 1;
                qrs_c(Beat_C) = pks_temp;
                qrs_i(Beat_C) = locs_temp;
                if locs_temp <= length(ecg_h)
                    [y_i_t,x_i_t] = max(ecg_h(locs_temp-round(0.150*fs):locs_temp));
                else
                    [y_i_t,x_i_t] = max(ecg_h(locs_temp-round(0.150*fs):end));
                end
                if y_i_t > THR_NOISE1 
                    Beat_C1 = Beat_C1 + 1;
                    qrs_i_raw(Beat_C1) = locs_temp-round(0.150*fs)+ (x_i_t - 1);
                    qrs_amp_raw(Beat_C1) = y_i_t;
                    SIG_LEV1 = 0.25*y_i_t + 0.75*SIG_LEV1;
                end
                not_nois = 1;
                SIG_LEV = 0.25*pks_temp + 0.75*SIG_LEV;
            end             
        else
            not_nois = 0;         
        end
    end
    if pks(i) >= THR_SIG      
        if Beat_C >= 3
            if (locs(i)-qrs_i(Beat_C)) <= round(0.3600*fs)
                Slope1 = mean(diff(ecg_m(locs(i)-round(0.075*fs):locs(i))));
                Slope2 = mean(diff(ecg_m(qrs_i(Beat_C)-round(0.075*fs):qrs_i(Beat_C))));
                if abs(Slope1) <= abs(0.5*(Slope2))
                    Noise_Count = Noise_Count + 1;
                    nois_c(Noise_Count) = pks(i);
                    nois_i(Noise_Count) = locs(i);
                    skip = 1;
                    NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1;
                    NOISE_LEV = 0.125*pks(i) + 0.875*NOISE_LEV; 
                else
                    skip = 0;
                end
            end
        end
        if skip == 0    
            Beat_C = Beat_C + 1;
            qrs_c(Beat_C) = pks(i);
            qrs_i(Beat_C) = locs(i);
            if y_i >= THR_SIG1  
                Beat_C1 = Beat_C1 + 1;
                if ser_back 
                    qrs_i_raw(Beat_C1) = x_i;
                else
                    qrs_i_raw(Beat_C1)= locs(i)-round(0.150*fs)+ (x_i - 1);
                end
                qrs_amp_raw(Beat_C1) =  y_i;
                SIG_LEV1 = 0.125*y_i + 0.875*SIG_LEV1;
            end
            SIG_LEV = 0.125*pks(i) + 0.875*SIG_LEV;
        end
    elseif (THR_NOISE <= pks(i)) && (pks(i) < THR_SIG)
        NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1;
        NOISE_LEV = 0.125*pks(i) + 0.875*NOISE_LEV;
    elseif pks(i) < THR_NOISE
        Noise_Count = Noise_Count + 1;
        nois_c(Noise_Count) = pks(i);
        nois_i(Noise_Count) = locs(i);    
        NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1;
        NOISE_LEV = 0.125*pks(i) + 0.875*NOISE_LEV;
    end
    if NOISE_LEV ~= 0 || SIG_LEV ~= 0
        THR_SIG = NOISE_LEV + 0.25*(abs(SIG_LEV - NOISE_LEV));
        THR_NOISE = 0.5*(THR_SIG);
    end
    if NOISE_LEV1 ~= 0 || SIG_LEV1 ~= 0
        THR_SIG1 = NOISE_LEV1 + 0.25*(abs(SIG_LEV1 - NOISE_LEV1));
        THR_NOISE1 = 0.5*(THR_SIG1);
    end
    SIGL_buf(i) = SIG_LEV;
    NOISL_buf(i) = NOISE_LEV;
    THRS_buf(i) = THR_SIG;
    SIGL_buf1(i) = SIG_LEV1;
    NOISL_buf1(i) = NOISE_LEV1;
    THRS_buf1(i) = THR_SIG1;
    skip = 0; not_nois = 0; ser_back = 0;    
end
qrs_i_raw = qrs_i_raw(1:Beat_C1);
qrs_amp_raw = qrs_amp_raw(1:Beat_C1);
qrs_c = qrs_c(1:Beat_C);
qrs_i = qrs_i(1:Beat_C);
if gr
    hold on,scatter(qrs_i,qrs_c,'m');
    hold on,plot(locs,NOISL_buf,'--k','LineWidth',2);
    hold on,plot(locs,SIGL_buf,'--r','LineWidth',2);
    hold on,plot(locs,THRS_buf,'--g','LineWidth',2);
    if any(ax), ax(~ax) = []; linkaxes(ax,'x'); zoom on; end
    figure;
    az(1)=subplot(311); plot(ecg_h); axis tight; hold on;
    scatter(qrs_i_raw,qrs_amp_raw,'m');
    az(2)=subplot(312); plot(ecg_m); axis tight; hold on;
    scatter(qrs_i,qrs_c,'m');
    az(3)=subplot(313); plot(ecg-mean(ecg)); axis tight; hold on;
    line(repmat(qrs_i_raw,[2 1]), repmat([min(ecg-mean(ecg))/2; max(ecg-mean(ecg))/2],size(qrs_i_raw)),...
         'LineWidth',2.5,'LineStyle','-.','Color','r');
    linkaxes(az,'x'); zoom on;
end
end


function [pks, locs] = findpeaks_minpeakdistance_fast(ecg, minDistSamp)
% Local maxima with minimum separation (greedy, highest first);
% equivalent to findpeaks(ecg, 'MinPeakDistance', minDistSamp).
if nargin < 2 || isempty(minDistSamp) || minDistSamp < 1
    minDistSamp = 1;
end
ecg = ecg(:);
isLocalMax = (ecg(2:end-1) > ecg(1:end-2)) & (ecg(2:end-1) > ecg(3:end));
allLocs = find(isLocalMax) + 1;
allPks  = ecg(allLocs);
if isempty(allLocs)
    pks = []; locs = [];
    return;
end
[~, sortIdx] = sort(allPks, 'descend');
sortedLocs = allLocs(sortIdx);
sortedPks  = allPks(sortIdx);
keep = true(size(sortedLocs));
for i = 1:length(sortedLocs)
    if keep(i)
        tooClose = (abs(sortedLocs - sortedLocs(i)) < minDistSamp);
        tooClose(1:i) = false;
        keep = keep & ~tooClose;
    end
end
locs = sortedLocs(keep);
pks  = sortedPks(keep);
[locs, order] = sort(locs);
pks = pks(order);
end
