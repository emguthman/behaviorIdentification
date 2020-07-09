%clustering of resident intruder data
%created 06-25-20, e mae guthman phd

close all
clearvars

%% init & load
disp('loading data')

%inits
sleapDir = '/Volumes/falkner/Mae/loopProject/aggroExp/sleapData/70777_377/tracked_labels';
cd(sleapDir)

sleapFiles = {'70777b_slpTracked_200224_extracted.mat';
    '70777b_slpTracked_200225_extracted.mat';
    '70777b_slpTracked_200226_extracted.mat';
    '70777b_slpTracked_200227_extracted.mat';
    '70777b_slpTracked_200228_extracted.mat';
    '70777b_slpTracked_200229_extracted.mat'};

im_dist = [];
home_mouse_length = [];
intruder_mouse_length = [];
home_dist_from_center = [];
intruder_dist_from_center = [];
nose_nose_dist = [];
nose_neck_dist = [];
nose_trunk_dist = [];
nose_tush_dist = [];
home_body_angle = [];
intruder_body_angle = [];
home_right_forelimb_angle = [];
home_left_forelimb_angle = [];
intruder_right_forelimb_angle = [];
intruder_left_forelimb_angle = [];
home_to_intruder_angle = [];
intruder_to_home_angle = [];
home_tail_angle_one = [];
home_tail_angle_two = [];
home_tail_angle_three = [];
home_tailbase_angle = [];
intruder_tail_angle_one = [];
intruder_tail_angle_two = [];
intruder_tail_angle_three = [];
intruder_tailbase_angle = [];
home_right_forelimb_velo = [];
intruder_right_forelimb_velo = [];
home_left_forelimb_velo = [];
intruder_left_forelimb_velo = [];
home_mouse_velo = [];
intruder_mouse_velo = [];
home_nose_velo = [];
intruder_nose_velo = [];
im_velo = [];
home_tailtip_velo = [];
intruder_tailtip_velo = [];
home_other_nose_angle = [];
intruder_other_nose_angle = [];
home_other_tush_angle = [];
intruder_other_tush_angle = [];

for ii = 1:length(sleapFiles)
    disp(['loading file ' num2str(ii)])

    load(sleapFiles{ii},'dist_im','dist_mouse','dist_cc','dist_nn','dist_nneck','dist_nt','dist_ntrunk','orientation_body','orientation_forelimb_left','orientation_forelimb_right','orientation_mouse',...
        'orientation_tail1','orientation_tail2','orientation_tail3','orientation_tailbase','velo_forelimb_left','velo_forelimb_right','velo_im','velo_nose','velo_trunk','velo_tailtip','orientation_other_nose',...
        'orientation_other_tush')
    
    %load features
    im_dist = [im_dist; dist_im.norm]; clear dist_im
    home_mouse_length = [home_mouse_length; dist_mouse.home.norm]; 
    intruder_mouse_length = [intruder_mouse_length; dist_mouse.intruder.norm]; clear dist_mouse
    home_dist_from_center = [home_dist_from_center; dist_cc.home.norm];
    intruder_dist_from_center = [intruder_dist_from_center; dist_cc.intruder.norm]; clear dist_cc
    nose_nose_dist = [nose_nose_dist; dist_nn.norm]; clear dist_nn
    nose_neck_dist = [nose_neck_dist; dist_nneck.norm]; clear dist_nneck
    nose_trunk_dist = [nose_trunk_dist; dist_ntrunk.norm]; clear dist_ntrunk
    nose_tush_dist = [nose_tush_dist; dist_nt.norm]; clear dist_nt
    home_body_angle = [home_body_angle; orientation_body.home.norm];
    intruder_body_angle = [intruder_body_angle; orientation_body.intruder.norm]; clear dist_cc
    home_right_forelimb_angle = [home_right_forelimb_angle; orientation_forelimb_right.home.norm];
    home_left_forelimb_angle = [home_left_forelimb_angle; orientation_forelimb_left.home.norm];
    intruder_right_forelimb_angle = [intruder_right_forelimb_angle; orientation_forelimb_right.intruder.norm]; clear orientation_forelimb_right
    intruder_left_forelimb_angle = [intruder_left_forelimb_angle; orientation_forelimb_left.intruder.norm]; clear orientation_forelimb_left
    home_to_intruder_angle = [home_to_intruder_angle; orientation_mouse.home.norm];
    intruder_to_home_angle = [intruder_to_home_angle; orientation_mouse.intruder.norm]; clear orientation_mouse
    home_tail_angle_one = [home_tail_angle_one; orientation_tail1.home.norm];
    home_tail_angle_two = [home_tail_angle_two; orientation_tail2.home.norm];
    home_tail_angle_three = [home_tail_angle_three; orientation_tail3.home.norm];
    home_tailbase_angle = [home_tailbase_angle; orientation_tailbase.home.norm];
    intruder_tail_angle_one = [intruder_tail_angle_one; orientation_tail1.intruder.norm]; clear orientation_tail1
    intruder_tail_angle_two = [intruder_tail_angle_two; orientation_tail2.intruder.norm]; clear orientation_tail2
    intruder_tail_angle_three = [intruder_tail_angle_three; orientation_tail3.intruder.norm]; clear orientation_tail3
    intruder_tailbase_angle = [intruder_tailbase_angle; orientation_tailbase.intruder.norm]; clear orientation_tailbase
    home_right_forelimb_velo = [home_right_forelimb_velo; velo_forelimb_right.home.norm];
    home_left_forelimb_velo = [home_left_forelimb_velo; velo_forelimb_left.home.norm];
    intruder_right_forelimb_velo = [intruder_right_forelimb_velo; velo_forelimb_right.intruder.norm]; clear velo_forelimb_right
    intruder_left_forelimb_velo = [intruder_left_forelimb_velo; velo_forelimb_left.intruder.norm]; clear velo_forelimb_left
    home_mouse_velo = [home_mouse_velo; velo_trunk.home.norm];
    intruder_mouse_velo = [intruder_mouse_velo; velo_trunk.intruder.norm]; clear velo_trunk
    home_nose_velo = [home_nose_velo; velo_nose.home.norm];
    intruder_nose_velo = [intruder_nose_velo; velo_nose.intruder.norm]; clear velo_nose
    im_velo = [im_velo; velo_im.norm]; clear velo_im
    home_tailtip_velo = [home_tailtip_velo; velo_tailtip.home.norm];
    intruder_tailtip_velo = [intruder_tailtip_velo; velo_tailtip.intruder.norm]; clear velo_tailtip
    home_other_nose_angle = [home_other_nose_angle; orientation_other_nose.home.norm];
    intruder_other_nose_angle = [intruder_other_nose_angle; orientation_other_nose.intruder.norm]; clear orientation_other_nose
    home_other_tush_angle = [home_other_tush_angle; orientation_other_tush.home.norm];
    intruder_other_tush_angle = [intruder_other_tush_angle; orientation_other_tush.intruder.norm]; clear orientation_other_tush
    
end

 %% tSNE
sFig = figure(1);
sFig.Position = [1010 350 2000 900];

%full tSNE
posture_tSNE.full = tsne([im_dist, nose_nose_dist, nose_tush_dist, home_mouse_velo, intruder_mouse_velo, im_velo, home_nose_velo, intruder_nose_velo, ...
    home_dist_from_center, intruder_dist_from_center, home_left_forelimb_velo, home_right_forelimb_velo, intruder_left_forelimb_velo, intruder_right_forelimb_velo, ...
    home_mouse_length, intruder_mouse_length, home_to_intruder_angle, intruder_to_home_angle, home_left_forelimb_angle, intruder_left_forelimb_angle ...
    home_body_angle, intruder_body_angle, home_right_forelimb_angle, intruder_right_forelimb_angle, nose_neck_dist, nose_trunk_dist, ...
    home_tailbase_angle, intruder_tailbase_angle, home_tail_angle_one, intruder_tail_angle_one, home_tail_angle_two, intruder_tail_angle_two, ...
    home_tail_angle_three, intruder_tail_angle_three, home_tailtip_velo, intruder_tailtip_velo, home_other_nose_angle, intruder_other_nose_angle, ...
    home_other_tush_angle, intruder_other_tush_angle]);

subplot(2,4,1)
scatter(posture_tSNE.full(:,1),posture_tSNE.full(:,2),7.5,'filled')

%plot tSNE
low_x.full = floor(min(posture_tSNE.full(:,1)));
high_x.full = ceil(max(posture_tSNE.full(:,1)));
low_y.full = floor(min(posture_tSNE.full(:,2)));
high_y.full = ceil(max(posture_tSNE.full(:,2)));

xlim([low_x.full high_x.full])
ylim([low_y.full high_y.full])
setAx(gca)

subplot(2,4,5)
tsne_x.full = low_x.full:0.25:high_x.full;
tsne_y.full = low_y.full:0.25:high_y.full;
[mesh_x.full,mesh_y.full] = meshgrid(tsne_x.full,tsne_y.full);
mesh_x.full = mesh_x.full(:);
mesh_y.full = mesh_y.full(:);
mesh_mat.full = [mesh_x.full mesh_y.full];
density.full = ksdensity(posture_tSNE.full,mesh_mat.full);
density.full = reshape(density.full, [length(tsne_y.full) length(tsne_x.full)]);
contourf(tsne_x.full,tsne_y.full,density.full,100,'edgecolor','none')
colormap(viridis)

%based on lindsay's
subplot(2,4,2)
posture_tSNE.lw = tsne([im_dist, nose_nose_dist, nose_tush_dist, home_mouse_velo, intruder_mouse_velo, im_velo, ...
    home_dist_from_center, intruder_dist_from_center, home_to_intruder_angle, intruder_to_home_angle]);
scatter(posture_tSNE.lw(:,1),posture_tSNE.lw(:,2),7.5,'filled');

%plot tSNE
low_x.lw = floor(min(posture_tSNE.lw(:,1)));
high_x.lw = ceil(max(posture_tSNE.lw(:,1)));
low_y.lw = floor(min(posture_tSNE.lw(:,2)));
high_y.lw = ceil(max(posture_tSNE.lw(:,2)));

xlim([low_x.lw high_x.lw])
ylim([low_y.lw high_y.lw])
setAx(gca)

subplot(2,4,6)
tsne_x.lw = low_x.lw:1:high_x.lw;
tsne_y.lw = low_y.lw:1:high_y.lw;
[mesh_x.lw,mesh_y.lw] = meshgrid(tsne_x.lw,tsne_y.lw);
mesh_x.lw = mesh_x.lw(:);
mesh_y.lw = mesh_y.lw(:);
mesh_mat.lw = [mesh_x.lw mesh_y.lw];
density.lw = ksdensity(posture_tSNE.lw,mesh_mat.lw);
density.lw = reshape(density.lw, [length(tsne_y.lw) length(tsne_x.lw)]);
contourf(tsne_x.lw,tsne_y.lw,density.lw,100,'edgecolor','none')
colormap(viridis)

%variable set 1
subplot(2,4,3)
posture_tSNE.musculus = tsne([im_dist, home_mouse_velo, intruder_mouse_velo, home_dist_from_center, intruder_dist_from_center, home_tailtip_velo, intruder_tailtip_velo ...
    home_other_nose_angle, intruder_other_nose_angle, home_other_tush_angle, intruder_other_tush_angle, home_mouse_length, intruder_mouse_length, home_tailbase_angle, intruder_tailbase_angle]);
scatter(posture_tSNE.musculus(:,1),posture_tSNE.musculus(:,2),7.5,'filled');

%plot tSNE
low_x.musculus = floor(min(posture_tSNE.musculus(:,1)));
high_x.musculus = ceil(max(posture_tSNE.musculus(:,1)));
low_y.musculus = floor(min(posture_tSNE.musculus(:,2)));
high_y.musculus = ceil(max(posture_tSNE.musculus(:,2)));

xlim([low_x.musculus high_x.musculus])
ylim([low_y.musculus high_y.musculus])
setAx(gca)

subplot(2,4,7)
tsne_x.musculus = low_x.musculus:1:high_x.musculus;
tsne_y.musculus = low_y.musculus:1:high_y.musculus;
[mesh_x.musculus,mesh_y.musculus] = meshgrid(tsne_x.musculus,tsne_y.musculus);
mesh_x.musculus = mesh_x.musculus(:);
mesh_y.musculus = mesh_y.musculus(:);
mesh_mat.musculus = [mesh_x.musculus mesh_y.musculus];
density.musculus = ksdensity(posture_tSNE.musculus,mesh_mat.musculus);
density.musculus = reshape(density.musculus, [length(tsne_y.musculus) length(tsne_x.musculus)]);
contourf(tsne_x.musculus,tsne_y.musculus,density.musculus,100,'edgecolor','none')
colormap(viridis)

%variable set 2
subplot(2,4,4)
posture_tSNE.vulgaris = tsne([im_dist, home_mouse_velo, intruder_mouse_velo, home_dist_from_center, intruder_dist_from_center, home_tailtip_velo, ...
    home_tailbase_angle, home_other_nose_angle, intruder_other_nose_angle, home_mouse_length, intruder_mouse_length]);
scatter(posture_tSNE.vulgaris(:,1),posture_tSNE.vulgaris(:,2),7.5,'filled');

%plot tSNE
low_x.vulgaris = floor(min(posture_tSNE.vulgaris(:,1)));
high_x.vulgaris = ceil(max(posture_tSNE.vulgaris(:,1)));
low_y.vulgaris = floor(min(posture_tSNE.vulgaris(:,2)));
high_y.vulgaris = ceil(max(posture_tSNE.vulgaris(:,2)));

xlim([low_x.vulgaris high_x.vulgaris])
ylim([low_y.vulgaris high_y.vulgaris])
setAx(gca)

subplot(2,4,8)
tsne_x.vulgaris = low_x.vulgaris:1:high_x.vulgaris;
tsne_y.vulgaris = low_y.vulgaris:1:high_y.vulgaris;
[mesh_x.vulgaris,mesh_y.vulgaris] = meshgrid(tsne_x.vulgaris,tsne_y.vulgaris);
mesh_x.vulgaris = mesh_x.vulgaris(:);
mesh_y.vulgaris = mesh_y.vulgaris(:);
mesh_mat.vulgaris = [mesh_x.vulgaris mesh_y.vulgaris];
density.vulgaris = ksdensity(posture_tSNE.vulgaris,mesh_mat.vulgaris);
density.vulgaris = reshape(density.vulgaris, [length(tsne_y.vulgaris) length(tsne_x.vulgaris)]);
contourf(tsne_x.vulgaris,tsne_y.vulgaris,density.vulgaris,100,'edgecolor','none')
colormap(viridis)

disp('tSNE complete')
