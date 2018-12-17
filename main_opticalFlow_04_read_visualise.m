%% Script for quantifying zebrafish heartbeat from LSM images
% Reads back processed data

%% User Inputs
% Process scan: #'s or (-1 for all)
% -- 1:20 (before drug); 21:80 (drug); 81:140 (wash)
% process_scan = [10,51,111];
% process_scan = process_scan-5;
% process_scan = 1:20;
% process_scan = 21:80;
% process_scan = 81:140;
% process_scan = 1;
process_scan = [20,50,110];
% process_scan = -1;

% Sample Number
sample_no = 7;

% Path to output folder
% output_folder = 'D:\LSM\sample8\sample8-output';
output_folder = sprintf('F:\\LSM\\sample%d\\sample%d-output',sample_no,sample_no);

% Strain Locations
if 1  % --Atrium
    switch sample_no
        case 1
            % sample1
            x_ = 563+90;
            y_ = 665-50;
        case 4
            % sample4
            x_ = 438-30;
            y_ = 768-10;
        case 5
            % sample5
            x_ = 644-20;
            y_ = 674-20;
        case 6
            % sample6
            x_ = 510-40;
            y_ = 596+20;
        case 7
            % sample7
            x_ = 584;
            y_ = 551;
        case 8
            % sample8
            x_ = 530-70;
            y_ = 700+20;
    end
else % --Ventricle
    switch sample_no
        case 1
            % sample1
            x_ = 405;
            y_ = 350;
        case 4
            % sample4
            x_ = 253+90;
            y_ = 528+00;
        case 5
            % sample5
            x_ = 512;
            y_ = 448;
        case 6
            % sample6
            x_ = 250;
            y_ = 349;
        case 7
            % sample7
            x_ = 441;
            y_ = 297;
        case 8
            % sample8
            x_ = 280-10;
            y_ = 481-10;
    end
end


% Fourier analysis
Fs = 20; % Sampling frequency (Hz)
T = 1/Fs; % Sampling period
L = 300; % Length of signal
% t = (0:L-1)*T; % Time vector
f = Fs*(0:(L/2))/L; % Freq space

%% MAIN
%% Read back dataID file structure
load(fullfile(output_folder,'dataID.mat'));

%% Reorganise order: Matlabs thinks 15.10 is before 15.8
t_hr = zeros([1 numel(data.img_timepoints)],'double');
t_min = zeros([1 numel(data.img_timepoints)],'double');
for id_ = 1:numel(data.img_timepoints)
    tmp_split = strsplit(data.img_timepoints{id_},'.');
    t_hr(id_) = str2double(tmp_split{1});
    t_min(id_) = str2double(tmp_split{2});
end
if sum(t_hr==23)~=0
    t_hr(t_hr==1)=25;
    t_hr(t_hr==0)=24;
end
[t_sorted, sort_idx] = sortrows([t_hr', t_min'],[1 2]);

data.img_timepoints_sorted = data.img_timepoints(sort_idx);
data.img_paths_sorted = data.img_paths(sort_idx);
data.img_scan_names_sorted = data.img_scan_names(sort_idx);

%% Display LSM
if ~exist('lsm_img','var')
    raw_images = matfile(fullfile(output_folder,data.img_timepoints_sorted{100},'raw_images.mat'));
    lsm_img = mean(raw_images.raw_images(:,:,1:55),3);
end
figure(1);
subplot(1,2,1);
imagesc(lsm_img);
colormap(gray);
subplot(1,2,2);
imagesc(lsm_img(mean(x_)+(-200:200),mean(y_)+(-200:200)));
colormap(gray);

%% Read back each timeseries
% Which timepoints to process
if process_scan ~= -1
    id = process_scan;
else
    id = 1:numel(data.img_timepoints_sorted);
end

tmp_strain = matfile(fullfile(output_folder,data.img_timepoints_sorted{id_},'strain.mat'));
tmp_strain_t = squeeze(tmp_strain.strain(x_,y_,:));
strain_t = zeros([length(tmp_strain_t), numel(id)]);
strain_tz = strain_t;
P1_max_amp = zeros([1 numel(id)]);
P1_max_freq = zeros([1 numel(id)]);

c_ = 0;
for id_ = id
    c_ = c_+1;
    % >>console
    fprintf('Processing %d of %d sequences: \n', c_, numel(id));
    fprintf('-- %s \n',data.img_timepoints_sorted{id_});
    % >>end
    tmp_strain = matfile(fullfile(output_folder,data.img_timepoints_sorted{id_},'strain.mat'));
    tmp_strain_t = squeeze(mean(mean(tmp_strain.strain(x_,y_,:),1),2));
    
    strain_t(:,c_) = tmp_strain_t;%/T; % Strain velocity val/period
    
    for t_ = 1:size(tmp_strain_t,1)
        strain_tz(t_,c_) = sum(tmp_strain_t(1:t_));%/T);
    end
    strain_tz_gs = smooth(strain_tz(:,c_),71,'loess')';
    strain_tz(:,c_) = strain_tz(:,c_)-strain_tz_gs';
    
    % Fourier analysis
    F_strain_tz = fft(strain_tz(:,c_),L);
    P2 = abs(F_strain_tz/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
%     [P1_max_psd(c_), P1_max_idx] = max(P1);
    
%     P1_max_freq(c_) = f(P1_max_idx);
    
    [pks, loc] = findpeaks(P1);
    [pks_s, s_id_] = sort(pks,'descend');
    % min loc of 2 peaks
    P1_max_idx = min(loc(s_id_(1:2)));
    P1_max_amp(c_) = sqrt( sum(P1((1:2)*P1_max_idx) .* f((1:2)*P1_max_idx)' ) );
    P1_max_freq(c_) = f(P1_max_idx);
    
    
    
%     figure;
%     plot(f,P1)
%     title('Single-Sided Amplitude Spectrum of X(t)')
%     xlabel('f (Hz)')
%     ylabel('|P1(f)|')
 
end

%% Visualisation
% Timeseries
t_s = linspace(0,length(id),length(strain_t(:)));

figure(2);
amplitude_p2p = 2*smooth(P1_max_amp,7,'lowess');
plot(amplitude_p2p);
title('Amplitude (Peak to peak)');

figure;
beat_rate = 60*smooth(medfilt1(P1_max_freq),7,'lowess');
plot(beat_rate);
title('Beat rate (bpm)');

figure;
plot(t_s,strain_t(:));
figure;
plot(t_s,strain_tz(:)/T/2);
figure;
strain_tz_all = strain_tz(:);
plot(strain_tz_all([1:150,301:450,601:750])/T/2);


while 0
    
    save(fullfile(output_folder,'atrium_amplitude_beat.mat'),'amplitude_p2p','beat_rate');
    
    save(fullfile(output_folder,'ventricle_amplitude_beat.mat'),'amplitude_p2p','beat_rate');
    
end
