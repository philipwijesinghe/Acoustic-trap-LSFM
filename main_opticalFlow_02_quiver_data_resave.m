%% Displays quiver plots of a sample 

%% Sample 6
% 15 processed to flow
for process_scan = 1:15

% Sample Number
sample_no = 6;

% Path to output folder
% output_folder = 'D:\LSM\sample8\sample8-output';
output_folder = sprintf('F:\\LSM\\sample%d\\sample%d-output-quiver',sample_no,sample_no);


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

% Processed id = 10:10:150, so make sure it's a multiple of 10
% raw_images = matfile(fullfile(output_folder,data.img_timepoints_sorted{process_scan*10},'raw_images.mat'));
% flow_data = matfile(fullfile(output_folder,data.img_timepoints_sorted{process_scan*10},'opticalFlowStruct.mat'));
load(fullfile(output_folder,data.img_timepoints_sorted{process_scan*10},'opticalFlowStruct.mat'));

%% Resave data
for id_ = 1:numel(flow)
    tmp = flow{id_};
    tmp2.Vx = tmp.Vx;
    tmp2.Vy = tmp.Vy;
    flow{id_} = tmp2;
end
save(fullfile(output_folder,data.img_timepoints_sorted{process_scan*10},'opticalFlowStruct.mat'),'flow','-v7.3');
clear flow

process_scan
end


% lsm_img = raw_images.raw_images(:,:,1);
% figure(1);
% imagesc(lsm_img);
% colormap(gray);






