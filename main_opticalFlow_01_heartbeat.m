%% Script for quantifying zebrafish heartbeat from LSM images

%% User Inputs
% Path to data folder
% data.path = 'E:\pwijesinghe\sample1\sample1finalVerapomil';
data.path = 'E:\pwijesinghe\sample4\sample 4';

% Process scan: #'s or (-1 for all)
process_scan = -1;

% Path to output folder
data.output_folder = 'E:\pwijesinghe\sample4\sample4-output';

% Save outputs
data.save = 1;
% Save opticalFlow?
data.save_flow = 0;

%% Load in folder structure
if ~exist(data.output_folder,'dir')
    mkdir(data.output_folder);
end
% Read the folder structure in the data.path
data.dir = dir(data.path); %read dir structure
data.dir = data.dir(3:end); %drop '.' and '..'
data.dir = data.dir([data.dir.isdir]); %keep folders only
% Extract all sample names
data.scan_names = {data.dir.name};
% Create a cell array of all sample folders
data.scan_folders = cellfun(@(x) fullfile(data.path,x),data.scan_names,'UniformOutput',0);
% Parse each folder for more folders and append to a structure
data.img_paths = {};
data.img_scan_names = {};
data.img_timepoints = {};
for iScan_ = 1:numel(data.scan_names)
    scanDir = dir( data.scan_folders{iScan_} ); %read dir structure
    scanDir = scanDir(3:end); %drop '.' and '..'
    scanDir = scanDir([scanDir.isdir]); %keep folders only
    timepointNames = {scanDir.name}; %extract sample names
    timepointFolders = cellfun(@(x) fullfile(data.scan_folders{iScan_},x),timepointNames,'UniformOutput',0); %cell array of scan folders
    for iTime_ = 1:numel(timepointNames)
        %octName = dir( [timepointFolders{iTime_} '\*-oct-snr-db.mat'] );
        data.img_paths{end+1} = timepointFolders{iTime_};
        data.img_scan_names{end+1} = data.scan_names{iScan_};
        data.img_timepoints{end+1} = timepointNames{iTime_};
    end
end
clear timepointNames timepointFolders scanDir iScan_ iTime_

%% Save data ids
save(fullfile(data.output_folder,'dataID.mat'),'data','-v7.3');

%% Process data
% Which timepoints to process
if process_scan ~= -1
    iImg = process_scan;
else
    iImg = 1:numel(data.img_paths);
end

% Process all requested timepoints
c_ = 0;
for iImg_ = iImg
    c_ = c_+1;
    tmp_path = data.img_paths{iImg_};
    tmp_dir = dir(tmp_path); %read dir structure
    tmp_dir = tmp_dir(3:end); %drop '.' and '..'
    tmp_dir = tmp_dir(~[tmp_dir.isdir]); %keep files only
    
    % >>console
    fprintf('Processing %d of %d sequences: \n', c_, numel(iImg));
    fprintf('-- %s \n',data.img_timepoints{iImg_});
    % >>end
    
    % Initialise optical flow
%     opticFlow = opticalFlowLK;
    opticFlow = opticalFlowHS('Smoothness',3);
%     opticFlow = opticalFlowFarneback;
    
    % Preallocate resources
    tmp_img_u16 = imread(fullfile(tmp_path, tmp_dir(1).name));
    strain = zeros([size(tmp_img_u16,1),size(tmp_img_u16,2),numel(tmp_dir)],'single');
    raw_images = zeros([size(tmp_img_u16,1),size(tmp_img_u16,2),numel(tmp_dir)],'double');
    
    % Load image stack
    for iImgNo_ = 1:numel(tmp_dir)
        % >>console
        if mod(iImgNo_-1,50)==0
            fprintf('---- %d of %d \n',iImgNo_, numel(tmp_dir));
        end
        % >>end
        
        % Read image out of sequence
        tmp_img_path = fullfile(tmp_path, tmp_dir(iImgNo_).name);
        
        % Scale to a sensible range [0 255]
        tmp_img_u16 = imread(tmp_img_path);
        tmp_img = double(tmp_img_u16)/65535*12000;
        
        % Store images
        raw_images(:,:,iImgNo_) = tmp_img;
        
        % Estimate optical flow
        if data.save_flow
            flow{iImgNo_} = estimateFlow(opticFlow,tmp_img); %#ok<SAGROW>
            % Estimate strain
            tmp_strain = divergence( flow{iImgNo_}.Vx, flow{iImgNo_}.Vy );
        else
            flow = estimateFlow(opticFlow,tmp_img);
            % Estimate strain
            tmp_strain = divergence( flow.Vx, flow.Vy );
        end
        
        % Record strain
        strain(:,:,iImgNo_) = gaussian_smooth(tmp_strain, [50 50]);
    end
    % >>console
    fprintf('---- %d of %d ... done \n', numel(tmp_dir), numel(tmp_dir));
    fprintf('-- saving \n');
    % >>end
    
    if data.save
        % Save flow data
        if ~exist(fullfile(data.output_folder,data.img_timepoints{iImg_}),'dir')
            mkdir(fullfile(data.output_folder,data.img_timepoints{iImg_}));
        end
        
        save(fullfile(data.output_folder,data.img_timepoints{iImg_},'strain.mat'),'strain','-v7.3');
        save(fullfile(data.output_folder,data.img_timepoints{iImg_},'raw_images.mat'),'raw_images','-v7.3');
        if data.save_flow
            save(fullfile(data.output_folder,data.img_timepoints{iImg_},'opticalFlowStruct.mat'),'flow','-v7.3');
        end
    end
    
    % >>console
    fprintf('------- \n');
    % >>end
end




%% END
if false
    % ---TEMPS
    % Plot output
    figure;
    imagesc(raw_images(:,:,10));
    colormap(gray);
    figure;
    imagesc(squeeze(strain(:,:,10)));
    
    x_ = 516;
    y_ = 507;
    t = 1:size(strain,3);
    figure;
    plot(t,squeeze(strain(x_,y_,:)));
    
    for t_ = 1:size(strain,3)
        strain_tz(t_) = sum(squeeze(strain(x_,y_,1:t_)));
    end
    strain_tz_gs = smooth(strain_tz,71,'loess')';
    figure;
    plot(t,strain_tz_gs);
    
    figure;
    plot(t,strain_tz-strain_tz_gs);
    
    figure;
    plot(t,squeeze(strain(x_,y_,:)), t,strain_tz-strain_tz_gs);
    
end

