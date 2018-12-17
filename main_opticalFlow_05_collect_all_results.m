% Read all datasets

for id_ = [1,4,5,6,7,8]
    
    output_folder = ['F:\LSM\sample' num2str(id_) '\sample' num2str(id_) '-output'];
    
    load(fullfile(output_folder,'atrium_amplitude_beat.mat'));
    all_atrium_amp{id_} = amplitude_p2p;
    all_atrium_beat{id_} = beat_rate;
    
    load(fullfile(output_folder,'ventricle_amplitude_beat.mat'));
    all_ventricle_amp{id_} = amplitude_p2p;
    all_ventricle_beat{id_} = beat_rate;
    
        drug_offset_min(id_) = 150 - length(amplitude_p2p);
    x{id_} = (1:length(amplitude_p2p)) + drug_offset_min(id_);
    
    % Normalise beats
    beat_rate = all_ventricle_beat{id_};
    all_ventricle_beat{id_} = beat_rate ./ mean(beat_rate(21-drug_offset_min:11-drug_offset_min+10));
    beat_rate = all_atrium_beat{id_};
    all_atrium_beat{id_} = beat_rate ./ mean(beat_rate(21-drug_offset_min:11-drug_offset_min+10));
end


%% DISPLAY PROPERTIES
XLIMIT = [10 150];
BEATAMP = [0 1.5];
BEATLABELY = 'Normalised Beat Rate';
% BEATAMP = [40 150];
% BEATLABELY = 'Beat Rate (bpm)';

%% VISUALISE

figure(1);
hold on;
for id_ = [1,4,5,6,7,8]
    plot(x{id_}, all_atrium_beat{id_},'Color',[0.7 0.7 0.7],...
        'LineStyle',':');
end
% xlabel('Time (min)');
% ylabel('Beat Rate (bpm)');


figure(2);
hold on;
for id_ = [1,4,5,6,7,8]
    plot(x{id_}, all_atrium_amp{id_},'Color',[0.7 0.7 0.7],...
        'LineStyle',':');
end
% xlabel('Time (min)');
% ylabel('Beat Amplitude (atrium)');


figure(3);
hold on;
for id_ = [1,4,5,6,7,8]
    plot(x{id_}, all_ventricle_beat{id_},'Color',[0.7 0.7 0.7],...
        'LineStyle',':');
end
% xlabel('Time (min)');
% ylabel('Beat Rate (bpm)');

figure(4);
hold on;
for id_ = [1,4,5,6,7,8]
    plot(x{id_}, all_ventricle_amp{id_},'Color',[0.7 0.7 0.7],...
        'LineStyle',':');
end
% xlabel('Time (min)');
% ylabel('Beat Amplitude (ventricle)');






%% MEANS
figure(1);
hold on;
c_ = 0;
clear mean_beat_a
for id_ = [1,4,5,6,7,8]
    c_ = c_+1;
    tmp = all_atrium_beat{id_};
    tmp = [NaN([drug_offset_min(id_),1]); tmp];offst = 0;
%     tmp = tmp(11-drug_offset_min(id_):150-drug_offset_min(id_));offst = 10;
    mean_beat_a(:,c_) = tmp;
end
env_max = nanmean(mean_beat_a,2)+nanstd(mean_beat_a,[],2);
env_min = nanmean(mean_beat_a,2)-nanstd(mean_beat_a,[],2);
f_x = [(1:size(mean_beat_a,1)), flip((1:size(mean_beat_a,1)))]+offst;
f_y = [env_max', flip(env_min)'];

fill(f_x,f_y,[0.9 0.9 0.9],'LineStyle','none','FaceAlpha',0.5);
plot((1:size(mean_beat_a,1))+offst,nanmean(mean_beat_a,2),'k','LineWidth',1.5);
plot((1:size(mean_beat_a,1))+offst,env_max,'--k');
plot((1:size(mean_beat_a,1))+offst,env_min,'--k');
line([30,30],[0 150],'LineWidth',2,'LineStyle',':','Color',[0 0 0]);
line([90,90],[0 150],'LineWidth',2,'LineStyle',':','Color',[0 0 0]);

ylim(BEATAMP)
xlim(XLIMIT);

xlabel('Time (min)');
ylabel(BEATLABELY);
title('Beat Rate - Atrium')






figure(2);
hold on;
c_ = 0;
clear mean_amp_a
for id_ = [1,4,5,6,7,8]
    c_ = c_+1;
    tmp = all_atrium_amp{id_};
    tmp = [NaN([drug_offset_min(id_),1]); tmp];offst = 0;
%     tmp = tmp(11-drug_offset_min(id_):150-drug_offset_min(id_));offst = 10;
    mean_amp_a(:,c_) = tmp;
end
env_max = nanmean(mean_amp_a,2)+nanstd(mean_amp_a,[],2);
env_min = nanmean(mean_amp_a,2)-nanstd(mean_amp_a,[],2);
f_x = [(1:size(mean_amp_a,1)), flip((1:size(mean_amp_a,1)))]+offst;
f_y = [env_max', flip(env_min)'];

fill(f_x,f_y,[0.9 0.9 0.9],'LineStyle','none','FaceAlpha',0.5);
plot((1:size(mean_amp_a,1))+offst,nanmean(mean_amp_a,2),'k','LineWidth',1.5);
plot((1:size(mean_amp_a,1))+offst,env_max,'--k');
plot((1:size(mean_amp_a,1))+offst,env_min,'--k');
line([30,30],[0 1],'LineWidth',2,'LineStyle',':','Color',[0 0 0]);
line([90,90],[0 1],'LineWidth',2,'LineStyle',':','Color',[0 0 0]);

ylim([0 1]);
xlim(XLIMIT);

xlabel('Time (min)');
ylabel('Beat Amplitude (\DeltaV / V)');
title('Beat Amplitude - Atrium')




figure(3);
hold on;
c_ = 0;
clear mean_beat_v
for id_ = [1,4,5,6,7,8]
    c_ = c_+1;
    tmp = all_ventricle_beat{id_};
    tmp = [NaN([drug_offset_min(id_),1]); tmp];offst = 0;
%     tmp = tmp(11-drug_offset_min(id_):150-drug_offset_min(id_));offst = 10;
    mean_beat_v(:,c_) = tmp;
end
env_max = nanmean(mean_beat_v,2)+nanstd(mean_beat_v,[],2);
env_min = nanmean(mean_beat_v,2)-nanstd(mean_beat_v,[],2);
f_x = [(1:size(mean_beat_v,1)), flip((1:size(mean_beat_v,1)))]+offst;
f_y = [env_max', flip(env_min)'];

fill(f_x,f_y,[0.9 0.9 0.9],'LineStyle','none','FaceAlpha',0.5);
plot((1:size(mean_beat_v,1))+offst,nanmean(mean_beat_v,2),'k','LineWidth',1.5);
plot((1:size(mean_beat_v,1))+offst,env_max,'--k');
plot((1:size(mean_beat_v,1))+offst,env_min,'--k');
line([30,30],[0 150],'LineWidth',2,'LineStyle',':','Color',[0 0 0]);
line([90,90],[0 150],'LineWidth',2,'LineStyle',':','Color',[0 0 0]);

ylim(BEATAMP)
xlim(XLIMIT);

xlabel('Time (min)');
ylabel(BEATLABELY);
title('Beat Rate - Ventricle');


figure(4);
hold on;
c_ = 0;
clear mean_amp_v
for id_ = [1,4,5,6,7,8]
    c_ = c_+1;
    tmp = all_ventricle_amp{id_};
    tmp = [NaN([drug_offset_min(id_),1]); tmp];offst = 0;
%     tmp = tmp(11-drug_offset_min(id_):150-drug_offset_min(id_));offst = 10;
    mean_amp_v(:,c_) = tmp;
end
env_max = nanmean(mean_amp_v,2)+nanstd(mean_amp_v,[],2);
env_min = nanmean(mean_amp_v,2)-nanstd(mean_amp_v,[],2);
f_x = [(1:size(mean_amp_v,1)), flip((1:size(mean_amp_v,1)))]+offst;
f_y = [env_max', flip(env_min)'];

fill(f_x,f_y,[0.9 0.9 0.9],'LineStyle','none','FaceAlpha',0.5);
plot((1:size(mean_amp_v,1))+offst,nanmean(mean_amp_v,2),'k','LineWidth',1.5);
plot((1:size(mean_amp_v,1))+offst,env_max,'--k');
plot((1:size(mean_amp_v,1))+offst,env_min,'--k');
line([30,30],[0 1],'LineWidth',2,'LineStyle',':','Color',[0 0 0]);
line([90,90],[0 1],'LineWidth',2,'LineStyle',':','Color',[0 0 0]);

ylim([0.0 0.7]);
xlim(XLIMIT);

xlabel('Time (min)');
ylabel('Beat Amplitude (\DeltaV / V)');
title('Beat Amplitude - Ventricle')



