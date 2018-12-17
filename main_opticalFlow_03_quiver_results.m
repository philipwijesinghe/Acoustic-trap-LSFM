%% Displays quiver plots of a sample 

%% Sample 6
% 15 processed to flow
process_scan = 1;

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
raw_images = matfile(fullfile(output_folder,data.img_timepoints_sorted{process_scan*10},'raw_images.mat'));
flow_data = matfile(fullfile(output_folder,data.img_timepoints_sorted{process_scan*10},'opticalFlowStruct.mat'));

%% Visualise
% View a series of frames as intensity image with quiver plot over it

% Read in the current frame from disk (lsm image and flow vectors
lsm_img_all = raw_images.raw_images;
flow_img_all = flow_data.flow;
% flow_img_vx = flow_img{1,1}.Vx;
% flow_img_vy = flow_img{1,1}.Vy;

for frame_id = 9:18
    
%     % Read in the current frame from disk (lsm image and flow vectors
%     lsm_img = raw_images.raw_images(:,:,frame_id);
%     flow_img = flow_data.flow(1,frame_id);
%     flow_img_vx = flow_img{1,1}.Vx;
%     flow_img_vy = flow_img{1,1}.Vy;

    lsm_img = lsm_img_all(:,:,frame_id);
    flow_img = flow_img_all(1,frame_id);
    flow_img_vx = flow_img{1,1}.Vx;
    flow_img_vy = flow_img{1,1}.Vy;

    
    % Calculate the X,Y coordinates for the quiver
    [X,Y] = meshgrid(linspace(1,size(lsm_img,2),size(flow_img_vx,2)),linspace(1,size(lsm_img,1),size(flow_img_vx,1)));
    
    % figure;imagesc(flow_img_vx);colorbar;caxis([-4 4])
    
    % Perform complex form smoothing weighted by LSM 4th root of intensity
    flow_img_cplx = (lsm_img).^(1/4).*exp(1i*angle(complex(flow_img_vx,flow_img_vy)));
    % flow_img_cplx = complex(flow_img_vx,flow_img_vy);
    flow_img_cplx_mean = gaussian_smooth(flow_img_cplx,[40 40]);
    flow_img_vx_mean = real(flow_img_cplx_mean);
    flow_img_vy_mean = imag(flow_img_cplx_mean);
    
%     figure;imagesc(flow_img_vx_mean);colorbar;caxis([-4 4])
%     figure;imagesc(flow_img_vy_mean);colorbar;caxis([-4 4])
    
    % Smooth the lsm image to allow for thresholding
    lsm_img_mean = gaussian_smooth(lsm_img,[5 5]);
    lsm_img_mean_block = gaussian_smooth(lsm_img,[50 50]);
    lsm_img_mean_diff = lsm_img_mean-lsm_img_mean_block;
    lsm_img_msk_adapt = (lsm_img_mean_diff>-1)&(lsm_img_mean>30);
    lsm_img_msk_adapt = imclose(lsm_img_msk_adapt,strel('disk', 95, 0));
%     figure;imagesc(lsm_img_msk_adapt);
%     figure;imagesc(lsm_img_mean_diff);
    
    % figure;imagesc(lsm_img_mean);colorbar;caxis([0 100])
%     figure;imagesc(thresholdLocally(lsm_img_mean/100))
    
    % Threshold quiver to not plot in low intensity areas
    flow_img_vx_msk = flow_img_vx_mean.*lsm_img_msk_adapt;%(lsm_img_mean>30);
    flow_img_vy_msk = flow_img_vy_mean.*lsm_img_msk_adapt;%(lsm_img_mean>30);
    
    % Subsample quiver by skip #
    skip = 15;
    
    % FIGURE
    %     figure(200+frame_id);
    figure(200);
    set(gcf,'Position',[0 0 1024 1024])
    imagesc(lsm_img);caxis([0 700]);
    colormap(gray);
    colorbar
    hold on
    quiver(X(1:skip:end,1:skip:end), Y(1:skip:end,1:skip:end),...
        flow_img_vx_msk(1:skip:end,1:skip:end), flow_img_vy_msk(1:skip:end,1:skip:end),...
        1,'LineWidth',1,'Color','r')
    axis image
    hold off
    set(gca,'Position',[0 0 1 1])
    set(gcf,'PaperPositionMode','auto')
    
    saveas(gcf,['./out_images/heart_quiver_s6f' num2str(frame_id) '.png'])
    %     close(200+frame_id);
    
    
    %% Strain
    strain_img = divergence( flow_img_vx_mean, flow_img_vx_mean );
    
    im = repmat(imgMat2Ratio(lsm_img,[0 700]),[1 1 3]);
    
    st = imgMat2Ratio(strain_img.*lsm_img_msk_adapt,[-0.025 0.025]);
    st = ind2rgb(uint16(st*255+1),gen_coolwarm(256));
        
    fuse = 1*(im.*st)+0*(st) ;
    
    figure(300)
    set(gcf,'Position',[0 0 1024 1024])
	%     imagesc(strain_img.*lsm_img_msk_adapt);caxis([-0.025 0.025]);colorbar;
	image(fuse)
    colormap(gen_coolwarm(256));
    axis image
    set(gca,'Position',[0 0 1 1])
    set(gcf,'PaperPositionMode','auto')
    saveas(gcf,['./out_images/heart_strain_s6f' num2str(frame_id) '.png'])
    
    
    
end



