%clustering of resident intruder data
%created 06-25-20, e mae guthman phd

close all
clearvars

%% init & load
disp('loading data')
%inits
sleapDir = '/Volumes/falkner/Mae/loopProject/aggroExp/sleapData/70777_377/tracked_labels';
cd(sleapDir)
fs = 30; %in Hz

%nodes
nose = 1;
trunk = 5;
tti = 6;
tailtip = 7;
tail0 = 8;
tail1 = 9;
tail2 = 10;
neck = 11;
forelimb_left = 14;
forelimb_right = 15;

noMice = 2; %number of mice tracked
home_mouse = 1; %track number for homecage mouse
intruder = 2; %track number for intruder mouse
noFeatures = 3; %number of features tracked
w = gausswin(fs/4); %.25 s window for gaussian filter
w = w/sum(w); %normalizes
gaussCut = 1:fs/8; %shifts filtered signal to match raw data
cc = colormap(viridis);

%load
occupancy_matrix = h5read('70777b_slpTracked_200229.analysis.h5','/track_occupancy');
tracks_matrix = h5read('70777b_slpTracked_200229.analysis.h5','/tracks');
nodeLabel = h5read('70777b_slpTracked_200229.analysis.h5','/node_names');

%crop to only when resident intruder is in
intruder_frames = find(occupancy_matrix(2,:) == 1);
intruderIn = intruder_frames(1); intruderOut = intruder_frames(end);
ri_tracks_matrix = tracks_matrix(intruderIn:intruderOut,:,:,1:noMice);
noFrames_intruder = size(ri_tracks_matrix,1);
t = 1:noFrames_intruder; t = t';
t_sec = t./fs;

disp('sleap tracks loaded')
disp('///////////////////////')

%% feature extraction
disp('begin feature extraction')

%inter-mouse distance
% defined as trunk to trunk distance

%get trunk x-y
trunk_pos.home.raw(:,:) = reshape(ri_tracks_matrix(:,trunk,:,home_mouse),[size(ri_tracks_matrix,1) 2]);
trunk_pos.intruder.raw(:,:) = reshape(ri_tracks_matrix(:,trunk,:,intruder),[size(ri_tracks_matrix,1) 2]);
trunk_pos.home.interp = -1.*ones(size(trunk_pos.home.raw));
trunk_pos.intruder.interp = -1.*ones(size(trunk_pos.intruder.raw));

%interpolate nans
for ii = 1:2
    nans.trunk.home(:,ii) = isnan(trunk_pos.home.raw(:,ii));
    nans.trunk.intruder(:,ii) = isnan(trunk_pos.intruder.raw(:,ii));
    trunk_pos.home.interp(nans.trunk.home(:,ii),ii) = interp1(t(~nans.trunk.home(:,ii)), trunk_pos.home.raw(~nans.trunk.home(:,ii),ii), t(nans.trunk.home(:,ii)));
    trunk_pos.home.interp(~nans.trunk.home(:,ii),ii) = trunk_pos.home.raw(~nans.trunk.home(:,ii),ii);
    trunk_pos.intruder.interp(nans.trunk.intruder(:,ii),ii) = interp1(t(~nans.trunk.intruder(:,ii)), trunk_pos.intruder.raw(~nans.trunk.intruder(:,ii),ii), t(nans.trunk.intruder(:,ii)));
    trunk_pos.intruder.interp(~nans.trunk.intruder(:,ii),ii) = trunk_pos.intruder.raw(~nans.trunk.intruder(:,ii),ii);
end

%inter-mouse distance
dist_im.raw = -1.*ones(noFrames_intruder,1);
dist_im.norm = -1.*ones(noFrames_intruder,1);
for ii = 1:noFrames_intruder
    dist_im.raw(ii) = norm(trunk_pos.home.interp(ii,:)-trunk_pos.intruder.interp(ii,:)); %intermouse distance in pixels
end

%interpolate extreme values
upper.dist_im = quantile(dist_im.raw,.97);
lower.dist_im = quantile(dist_im.raw,.03);
highs.dist_im = (dist_im.raw>upper.dist_im);
lows.dist_im = (dist_im.raw<lower.dist_im);
dist_im.raw(highs.dist_im) = interp1(t(~highs.dist_im), dist_im.raw(~highs.dist_im), t(highs.dist_im));
dist_im.raw(lows.dist_im) = interp1(t(~lows.dist_im), dist_im.raw(~lows.dist_im), t(lows.dist_im));

%filter
dist_im.filt = medfilt1(dist_im.raw,10);
dist_im.gauss = filter(w,1,dist_im.filt);
dist_im.gauss(1:gaussCut(end)+1) = [];
dist_im.norm = rescale(dist_im.gauss,-1,1);

%plot
timeseriesFig = figure(1);
timeseriesFig.Position = [40 40 1250 900];
% subplot(3,4,[1:3; 5:7; 9:11])
hold on
plot(t_sec((gaussCut(end)+2):end),dist_im.norm,'linewidth',2,'color',cc(1,:))
xlabel('time (s)')
setAx(gca)

disp('inter-mouse distance complete')

%distance from cage-center

%get cage center
% defined as middle point of mouse tracks
cageCenter = mean([max(trunk_pos.intruder.interp); min(trunk_pos.intruder.interp)]);

%distance from cc
dist_cc.home.raw = -1.*ones(noFrames_intruder,1);
dist_cc.home.norm = -1.*ones(noFrames_intruder,1);
dist_cc.intruder.raw = -1.*ones(noFrames_intruder,1);
dist_cc.intruder.norm = -1.*ones(noFrames_intruder,1);
for ii = 1:noFrames_intruder
    dist_cc.home.raw(ii) = norm(trunk_pos.home.interp(ii,:)-cageCenter); %distance in pixels
    dist_cc.intruder.raw(ii) = norm(trunk_pos.intruder.interp(ii,:)-cageCenter); %distance in pixels
end

%interpolate extreme values
upper.dist_cc.home = quantile(dist_cc.home.raw,.97);
lower.dist_cc.home = quantile(dist_cc.home.raw,.03);
highs.dist_cc.home = (dist_cc.home.raw>upper.dist_cc.home);
lows.dist_cc.home = (dist_cc.home.raw<lower.dist_cc.home);
dist_cc.home.raw(highs.dist_cc.home) = interp1(t(~highs.dist_cc.home), dist_cc.home.raw(~highs.dist_cc.home), t(highs.dist_cc.home));
dist_cc.home.raw(lows.dist_cc.home) = interp1(t(~lows.dist_cc.home), dist_cc.home.raw(~lows.dist_cc.home), t(lows.dist_cc.home));
upper.dist_cc.intruder = quantile(dist_cc.intruder.raw,.97);
lower.dist_cc.intruder = quantile(dist_cc.intruder.raw,.03);
highs.dist_cc.intruder = (dist_cc.intruder.raw>upper.dist_cc.intruder);
lows.dist_cc.intruder = (dist_cc.intruder.raw<lower.dist_cc.intruder);
dist_cc.intruder.raw(highs.dist_cc.intruder) = interp1(t(~highs.dist_cc.intruder), dist_cc.intruder.raw(~highs.dist_cc.intruder), t(highs.dist_cc.intruder));
dist_cc.intruder.raw(lows.dist_cc.intruder) = interp1(t(~lows.dist_cc.intruder), dist_cc.intruder.raw(~lows.dist_cc.intruder), t(lows.dist_cc.intruder));

%filter
dist_cc.home.filt = medfilt1(dist_cc.home.raw,10);
dist_cc.home.gauss = filter(w,1,dist_cc.home.filt);
dist_cc.home.gauss(1:gaussCut(end)+1) = [];
dist_cc.home.norm = rescale(dist_cc.home.gauss,-1,1);
dist_cc.intruder.filt = medfilt1(dist_cc.intruder.raw,10);
dist_cc.intruder.gauss = filter(w,1,dist_cc.intruder.filt);
dist_cc.intruder.gauss(1:gaussCut(end)+1) = [];
dist_cc.intruder.norm = rescale(dist_cc.intruder.gauss,-1,1);

%plot
plot(t_sec((gaussCut(end)+2):end),dist_cc.home.norm-2,'linewidth',2,'color',cc(1,:))
plot(t_sec((gaussCut(end)+2):end),dist_cc.intruder.norm-4,'linewidth',2,'color',cc(1,:))

disp('distance from cage center complete')

%nose to nose distance
%get x-y
nose_pos.home.raw(:,:) = reshape(ri_tracks_matrix(:,nose,:,home_mouse),[size(ri_tracks_matrix,1) 2]);
nose_pos.intruder.raw(:,:) = reshape(ri_tracks_matrix(:,nose,:,intruder),[size(ri_tracks_matrix,1) 2]);
nose_pos.home.interp = -1.*ones(size(nose_pos.home.raw));
nose_pos.intruder.interp = -1.*ones(size(nose_pos.intruder.raw));

%interpolate nans
for ii = 1:2
    nans.nose.home(:,ii) = isnan(nose_pos.home.raw(:,ii));
    nans.nose.intruder(:,ii) = isnan(nose_pos.intruder.raw(:,ii));
    nose_pos.home.interp(nans.nose.home(:,ii),ii) = interp1(t(~nans.nose.home(:,ii)), nose_pos.home.raw(~nans.nose.home(:,ii),ii), t(nans.nose.home(:,ii)));
    nose_pos.home.interp(~nans.nose.home(:,ii),ii) = nose_pos.home.raw(~nans.nose.home(:,ii),ii);
    nose_pos.intruder.interp(nans.nose.intruder(:,ii),ii) = interp1(t(~nans.nose.intruder(:,ii)), nose_pos.intruder.raw(~nans.nose.intruder(:,ii),ii), t(nans.nose.intruder(:,ii)));
    nose_pos.intruder.interp(~nans.nose.intruder(:,ii),ii) = nose_pos.intruder.raw(~nans.nose.intruder(:,ii),ii);
end

%nose to nose distance
dist_nn.raw = -1.*ones(noFrames_intruder,1);
dist_nn.norm = -1.*ones(noFrames_intruder,1);
for ii = 1:noFrames_intruder
    dist_nn.raw(ii) = norm(nose_pos.home.interp(ii,:)-nose_pos.intruder.interp(ii,:)); %nose-nose distance in pixels
end

%interpolate extreme values
upper.dist_nn = quantile(dist_nn.raw,.97);
lower.dist_nn = quantile(dist_nn.raw,.03);
highs.dist_nn = (dist_nn.raw>upper.dist_nn);
lows.dist_nn = (dist_nn.raw<lower.dist_nn);
dist_nn.raw(highs.dist_nn) = interp1(t(~highs.dist_nn), dist_nn.raw(~highs.dist_nn), t(highs.dist_nn));
dist_nn.raw(lows.dist_nn) = interp1(t(~lows.dist_nn), dist_nn.raw(~lows.dist_nn), t(lows.dist_nn));

%filter
dist_nn.filt = medfilt1(dist_nn.raw,10);
dist_nn.gauss = filter(w,1,dist_nn.filt);
dist_nn.gauss(1:gaussCut(end)+1) = [];
dist_nn.norm = rescale(dist_nn.gauss,-1,1);

plot(t_sec((gaussCut(end)+2):end),dist_nn.norm-6,'linewidth',2,'color',cc(1,:))
disp('nose to nose distance complete')

%nose to tti distance
% home_nose, intruder_tti

%get x-y
tti_pos.intruder.raw(:,:) = reshape(ri_tracks_matrix(:,tti,:,intruder),[size(ri_tracks_matrix,1) 2]);
tti_pos.intruder.interp = -1.*ones(size(tti_pos.intruder.raw));

%interpolate nans
for ii = 1:2
    nans.tti.intruder(:,ii) = isnan(tti_pos.intruder.raw(:,ii));
    tti_pos.intruder.interp(nans.tti.intruder(:,ii),ii) = interp1(t(~nans.tti.intruder(:,ii)), tti_pos.intruder.raw(~nans.tti.intruder(:,ii),ii), t(nans.tti.intruder(:,ii)));
    tti_pos.intruder.interp(~nans.tti.intruder(:,ii),ii) = tti_pos.intruder.raw(~nans.tti.intruder(:,ii),ii);
end

%nose to tush distance
dist_nt.raw = -1.*ones(noFrames_intruder,1);
dist_nt.norm = -1.*ones(noFrames_intruder,1);
for ii = 1:noFrames_intruder
    dist_nt.raw(ii) = norm(nose_pos.home.interp(ii,:)-tti_pos.intruder.interp(ii,:)); %nose-tti distance in pixels
end

%interpolate extreme values
upper.dist_nt = quantile(dist_nt.raw,.97);
lower.dist_nt = quantile(dist_nt.raw,.03);
highs.dist_nt = (dist_nt.raw>upper.dist_nt);
lows.dist_nt = (dist_nt.raw<lower.dist_nt);
dist_nt.raw(highs.dist_nt) = interp1(t(~highs.dist_nt), dist_nt.raw(~highs.dist_nt), t(highs.dist_nt));
dist_nt.raw(lows.dist_nt) = interp1(t(~lows.dist_nt), dist_nt.raw(~lows.dist_nt), t(lows.dist_nt));

%filter
dist_nt.filt = medfilt1(dist_nt.raw,10);
dist_nt.gauss = filter(w,1,dist_nt.filt);
dist_nt.gauss(1:gaussCut(end)+1) = [];
dist_nt.norm = rescale(dist_nt.gauss,-1,1);

plot(t_sec((gaussCut(end)+2):end),dist_nt.norm-8,'linewidth',2,'color',cc(1,:))
disp('nose to tush distance complete')

%mouse length
% defined as nose to tti distance

%get home tti x-y
tti_pos.home.raw(:,:) = reshape(ri_tracks_matrix(:,tti,:,home_mouse),[size(ri_tracks_matrix,1) 2]);
tti_pos.home.interp = -1.*ones(size(tti_pos.home.raw));

%interpolate nans
for ii = 1:2
    nans.tti.home(:,ii) = isnan(tti_pos.home.raw(:,ii));
    tti_pos.home.interp(nans.tti.home(:,ii),ii) = interp1(t(~nans.tti.home(:,ii)), tti_pos.home.raw(~nans.tti.home(:,ii),ii), t(nans.tti.home(:,ii)));
    tti_pos.home.interp(~nans.tti.home(:,ii),ii) = tti_pos.home.raw(~nans.tti.home(:,ii),ii);
end

%distance
dist_mouse.home.raw = -1.*ones(noFrames_intruder,1);
dist_mouse.home.norm = -1.*ones(noFrames_intruder,1);
dist_mouse.intruder.raw = -1.*ones(noFrames_intruder,1);
dist_mouse.intruder.norm = -1.*ones(noFrames_intruder,1);
for ii = 1:noFrames_intruder
    dist_mouse.home.raw(ii) = norm(nose_pos.home.interp(ii,:)-tti_pos.home.interp(ii,:)); %mouse length in pixels
    dist_mouse.intruder.raw(ii) = norm(nose_pos.intruder.interp(ii,:)-tti_pos.intruder.interp(ii,:)); %mouse length in pixels
end

%interpolate extreme values
upper.dist_mouse.home = quantile(dist_mouse.home.raw,.97);
lower.dist_mouse.home = quantile(dist_mouse.home.raw,.03);
highs.dist_mouse.home = (dist_mouse.home.raw>upper.dist_mouse.home);
lows.dist_mouse.home = (dist_mouse.home.raw<lower.dist_mouse.home);
dist_mouse.home.raw(highs.dist_mouse.home) = interp1(t(~highs.dist_mouse.home), dist_mouse.home.raw(~highs.dist_mouse.home), t(highs.dist_mouse.home));
dist_mouse.home.raw(lows.dist_mouse.home) = interp1(t(~lows.dist_mouse.home), dist_mouse.home.raw(~lows.dist_mouse.home), t(lows.dist_mouse.home));
upper.dist_mouse.intruder = quantile(dist_mouse.intruder.raw,.97);
lower.dist_mouse.intruder = quantile(dist_mouse.intruder.raw,.03);
highs.dist_mouse.intruder = (dist_mouse.intruder.raw>upper.dist_mouse.intruder);
lows.dist_mouse.intruder = (dist_mouse.intruder.raw<lower.dist_mouse.intruder);
dist_mouse.intruder.raw(highs.dist_mouse.intruder) = interp1(t(~highs.dist_mouse.intruder), dist_mouse.intruder.raw(~highs.dist_mouse.intruder), t(highs.dist_mouse.intruder));
dist_mouse.intruder.raw(lows.dist_mouse.intruder) = interp1(t(~lows.dist_mouse.intruder), dist_mouse.intruder.raw(~lows.dist_mouse.intruder), t(lows.dist_mouse.intruder));

%filter
dist_mouse.home.filt = medfilt1(dist_mouse.home.raw,10);
dist_mouse.home.gauss = filter(w,1,dist_mouse.home.filt);
dist_mouse.home.gauss(1:gaussCut(end)+1) = [];
dist_mouse.home.norm = rescale(dist_mouse.home.gauss,-1,1);
dist_mouse.intruder.filt = medfilt1(dist_mouse.intruder.raw,10);
dist_mouse.intruder.gauss = filter(w,1,dist_mouse.intruder.filt);
dist_mouse.intruder.gauss(1:gaussCut(end)+1) = [];
dist_mouse.intruder.norm = rescale(dist_mouse.intruder.gauss,-1,1);

%plot
plot(t_sec((gaussCut(end)+2):end),dist_mouse.home.norm-10,'linewidth',2,'color',cc(128,:))
plot(t_sec((gaussCut(end)+2):end),dist_mouse.intruder.norm-12,'linewidth',2,'color',cc(128,:))

disp('mouse length complete')

%mouse velocity
% defined as trunk velocity
velo_trunk.home.raw = -1.*ones(noFrames_intruder-1,1);
velo_trunk.intruder.raw = -1.*ones(noFrames_intruder-1,1);
for ii = 1:noFrames_intruder-1
    velo_trunk.home.raw(ii) = norm(trunk_pos.home.interp(ii+1,:)-trunk_pos.home.interp(ii,:)); %velocity in pixels/frame
    velo_trunk.intruder.raw(ii) = norm(trunk_pos.intruder.interp(ii+1,:)-trunk_pos.intruder.interp(ii,:)); %velocity in pixels/frame
end

%interpolate extreme values
upper.velo_trunk.home = quantile(velo_trunk.home.raw,.97);
lower.velo_trunk.home = quantile(velo_trunk.home.raw,.03);
highs.velo_trunk.home = (velo_trunk.home.raw>upper.velo_trunk.home);
lows.velo_trunk.home = (velo_trunk.home.raw<lower.velo_trunk.home);
velo_trunk.home.raw(highs.velo_trunk.home) = interp1(t(~highs.velo_trunk.home), velo_trunk.home.raw(~highs.velo_trunk.home), t(highs.velo_trunk.home));
velo_trunk.home.raw(lows.velo_trunk.home) = interp1(t(~lows.velo_trunk.home), velo_trunk.home.raw(~lows.velo_trunk.home), t(lows.velo_trunk.home));
upper.velo_trunk.intruder = quantile(velo_trunk.intruder.raw,.97);
lower.velo_trunk.intruder = quantile(velo_trunk.intruder.raw,.03);
highs.velo_trunk.intruder = (velo_trunk.intruder.raw>upper.velo_trunk.intruder);
lows.velo_trunk.intruder = (velo_trunk.intruder.raw<lower.velo_trunk.intruder);
velo_trunk.intruder.raw(highs.velo_trunk.intruder) = interp1(t(~highs.velo_trunk.intruder), velo_trunk.intruder.raw(~highs.velo_trunk.intruder), t(highs.velo_trunk.intruder));
velo_trunk.intruder.raw(lows.velo_trunk.intruder) = interp1(t(~lows.velo_trunk.intruder), velo_trunk.intruder.raw(~lows.velo_trunk.intruder), t(lows.velo_trunk.intruder));

%filter
velo_trunk.home.filt = medfilt1(velo_trunk.home.raw,10);
velo_trunk.home.gauss = filter(w,1,velo_trunk.home.filt);
velo_trunk.home.gauss(1:gaussCut(end)) = [];
velo_trunk.home.norm = rescale(velo_trunk.home.gauss,-1,1);
velo_trunk.intruder.filt = medfilt1(velo_trunk.intruder.raw,10);
velo_trunk.intruder.gauss = filter(w,1,velo_trunk.intruder.filt);
velo_trunk.intruder.gauss(1:gaussCut(end)) = [];
velo_trunk.intruder.norm = rescale(velo_trunk.intruder.gauss,-1,1);

%plot
plot(t_sec((gaussCut(end)+2):end),velo_trunk.home.norm-14,'linewidth',2,'color',cc(1,:))
plot(t_sec((gaussCut(end)+2):end),velo_trunk.intruder.norm-16,'linewidth',2,'color',cc(1,:))

disp('mouse velocity complete')

%intermouse velocity
velo_im.raw = -1.*ones(noFrames_intruder-1,1);
for ii = 1:noFrames_intruder-1
    velo_im.raw(ii) = norm(dist_im.raw(ii+1,:)-dist_im.raw(ii,:)); %velocity in pixels/frame
end

%interpolate extreme values
upper.velo_im = quantile(velo_im.raw,.97);
lower.velo_im = quantile(velo_im.raw,.03);
highs.velo_im = (velo_im.raw>upper.velo_im);
lows.velo_im = (velo_im.raw<lower.velo_im);
velo_im.raw(highs.velo_im) = interp1(t(~highs.velo_im), velo_im.raw(~highs.velo_im), t(highs.velo_im));
velo_im.raw(lows.velo_im) = interp1(t(~lows.velo_im), velo_im.raw(~lows.velo_im), t(lows.velo_im));

%filter
velo_im.filt = medfilt1(velo_im.raw,10);
velo_im.gauss = filter(w,1,velo_im.filt);
velo_im.gauss(1:gaussCut(end)) = [];
velo_im.norm = rescale(velo_im.gauss,-1,1);

%plot
plot(t_sec((gaussCut(end)+2):end),velo_im.norm-18,'linewidth',2,'color',cc(1,:))

disp('inter-mouse velocity complete')

%nose velocity
velo_nose.home.raw = -1.*ones(noFrames_intruder-1,1);
velo_nose.intruder.raw = -1.*ones(noFrames_intruder-1,1);
for ii = 1:noFrames_intruder-1
    velo_nose.home.raw(ii) = norm(nose_pos.home.interp(ii+1,:)-nose_pos.home.interp(ii,:)); %velocity in pixels/frame
    velo_nose.intruder.raw(ii) = norm(nose_pos.intruder.interp(ii+1,:)-nose_pos.intruder.interp(ii,:)); %velocity in pixels/frame
end

%interpolate extreme values
upper.velo_nose.home = quantile(velo_nose.home.raw,.97);
lower.velo_nose.home = quantile(velo_nose.home.raw,.03);
highs.velo_nose.home = (velo_nose.home.raw>upper.velo_nose.home);
lows.velo_nose.home = (velo_nose.home.raw<lower.velo_nose.home);
velo_nose.home.raw(highs.velo_nose.home) = interp1(t(~highs.velo_nose.home), velo_nose.home.raw(~highs.velo_nose.home), t(highs.velo_nose.home));
velo_nose.home.raw(lows.velo_nose.home) = interp1(t(~lows.velo_nose.home), velo_nose.home.raw(~lows.velo_nose.home), t(lows.velo_nose.home));
upper.velo_nose.intruder = quantile(velo_nose.intruder.raw,.97);
lower.velo_nose.intruder = quantile(velo_nose.intruder.raw,.03);
highs.velo_nose.intruder = (velo_nose.intruder.raw>upper.velo_nose.intruder);
lows.velo_nose.intruder = (velo_nose.intruder.raw<lower.velo_nose.intruder);
velo_nose.intruder.raw(highs.velo_nose.intruder) = interp1(t(~highs.velo_nose.intruder), velo_nose.intruder.raw(~highs.velo_nose.intruder), t(highs.velo_nose.intruder));
velo_nose.intruder.raw(lows.velo_nose.intruder) = interp1(t(~lows.velo_nose.intruder), velo_nose.intruder.raw(~lows.velo_nose.intruder), t(lows.velo_nose.intruder));

%filter
velo_nose.home.filt = medfilt1(velo_nose.home.raw,10);
velo_nose.home.gauss = filter(w,1,velo_nose.home.filt);
velo_nose.home.gauss(1:gaussCut(end)) = [];
velo_nose.home.norm = rescale(velo_nose.home.gauss,-1,1);
velo_nose.intruder.filt = medfilt1(velo_nose.intruder.raw,10);
velo_nose.intruder.gauss = filter(w,1,velo_nose.intruder.filt);
velo_nose.intruder.gauss(1:gaussCut(end)) = [];
velo_nose.intruder.norm = rescale(velo_nose.intruder.gauss,-1,1);

%plot
plot(t_sec((gaussCut(end)+2):end),velo_nose.home.norm-20,'linewidth',2,'color',cc(16,:))
plot(t_sec((gaussCut(end)+2):end),velo_nose.intruder.norm-22,'linewidth',2,'color',cc(32,:))
disp('nose velocity complete')


%left paw velocity
%get left paw x-y
forelimb_left_pos.home.raw(:,:) = reshape(ri_tracks_matrix(:,forelimb_left,:,home_mouse),[size(ri_tracks_matrix,1) 2]);
forelimb_left_pos.intruder.raw(:,:) = reshape(ri_tracks_matrix(:,forelimb_left,:,intruder),[size(ri_tracks_matrix,1) 2]);
forelimb_left_pos.home.interp = -1.*ones(size(forelimb_left_pos.home.raw));
forelimb_left_pos.intruder.interp = -1.*ones(size(forelimb_left_pos.intruder.raw));

%interpolate nans
for ii = 1:2
    nans.forelimb_left.home(:,ii) = isnan(forelimb_left_pos.home.raw(:,ii));
    nans.forelimb_left.intruder(:,ii) = isnan(forelimb_left_pos.intruder.raw(:,ii));
    forelimb_left_pos.home.interp(nans.forelimb_left.home(:,ii),ii) = interp1(t(~nans.forelimb_left.home(:,ii)), forelimb_left_pos.home.raw(~nans.forelimb_left.home(:,ii),ii), t(nans.forelimb_left.home(:,ii)));
    forelimb_left_pos.home.interp(~nans.forelimb_left.home(:,ii),ii) = forelimb_left_pos.home.raw(~nans.forelimb_left.home(:,ii),ii);
    forelimb_left_pos.intruder.interp(nans.forelimb_left.intruder(:,ii),ii) = interp1(t(~nans.forelimb_left.intruder(:,ii)), forelimb_left_pos.intruder.raw(~nans.forelimb_left.intruder(:,ii),ii), t(nans.forelimb_left.intruder(:,ii)));
    forelimb_left_pos.intruder.interp(~nans.forelimb_left.intruder(:,ii),ii) = forelimb_left_pos.intruder.raw(~nans.forelimb_left.intruder(:,ii),ii);
end

%velocity
velo_forelimb_left.home.raw = -1.*ones(noFrames_intruder-1,1);
velo_forelimb_left.intruder.raw = -1.*ones(noFrames_intruder-1,1);
for ii = 1:noFrames_intruder-1
    velo_forelimb_left.home.raw(ii) = norm(forelimb_left_pos.home.interp(ii+1,:)-forelimb_left_pos.home.interp(ii,:)); %velocity in pixels/frame
    velo_forelimb_left.intruder.raw(ii) = norm(forelimb_left_pos.intruder.interp(ii+1,:)-forelimb_left_pos.intruder.interp(ii,:)); %velocity in pixels/frame
end

%interpolate extreme values
upper.velo_forelimb_left.home = quantile(velo_forelimb_left.home.raw,.97);
lower.velo_forelimb_left.home = quantile(velo_forelimb_left.home.raw,.03);
highs.velo_forelimb_left.home = (velo_forelimb_left.home.raw>upper.velo_forelimb_left.home);
lows.velo_forelimb_left.home = (velo_forelimb_left.home.raw<lower.velo_forelimb_left.home);
velo_forelimb_left.home.raw(highs.velo_forelimb_left.home) = interp1(t(~highs.velo_forelimb_left.home), velo_forelimb_left.home.raw(~highs.velo_forelimb_left.home), t(highs.velo_forelimb_left.home));
velo_forelimb_left.home.raw(lows.velo_forelimb_left.home) = interp1(t(~lows.velo_forelimb_left.home), velo_forelimb_left.home.raw(~lows.velo_forelimb_left.home), t(lows.velo_forelimb_left.home));
upper.velo_forelimb_left.intruder = quantile(velo_forelimb_left.intruder.raw,.97);
lower.velo_forelimb_left.intruder = quantile(velo_forelimb_left.intruder.raw,.03);
highs.velo_forelimb_left.intruder = (velo_forelimb_left.intruder.raw>upper.velo_forelimb_left.intruder);
lows.velo_forelimb_left.intruder = (velo_forelimb_left.intruder.raw<lower.velo_forelimb_left.intruder);
velo_forelimb_left.intruder.raw(highs.velo_forelimb_left.intruder) = interp1(t(~highs.velo_forelimb_left.intruder), velo_forelimb_left.intruder.raw(~highs.velo_forelimb_left.intruder), t(highs.velo_forelimb_left.intruder));
velo_forelimb_left.intruder.raw(lows.velo_forelimb_left.intruder) = interp1(t(~lows.velo_forelimb_left.intruder), velo_forelimb_left.intruder.raw(~lows.velo_forelimb_left.intruder), t(lows.velo_forelimb_left.intruder));

%filter
velo_forelimb_left.home.filt = medfilt1(velo_forelimb_left.home.raw,10);
velo_forelimb_left.home.gauss = filter(w,1,velo_forelimb_left.home.filt);
velo_forelimb_left.home.gauss(1:gaussCut(end)) = [];
velo_forelimb_left.home.norm = rescale(velo_forelimb_left.home.gauss,-1,1);
velo_forelimb_left.intruder.filt = medfilt1(velo_forelimb_left.intruder.raw,10);
velo_forelimb_left.intruder.gauss = filter(w,1,velo_forelimb_left.intruder.filt);
velo_forelimb_left.intruder.gauss(1:gaussCut(end)) = [];
velo_forelimb_left.intruder.norm = rescale(velo_forelimb_left.intruder.gauss,-1,1);

%plot
plot(t_sec((gaussCut(end)+2):end),velo_forelimb_left.home.norm-24,'linewidth',2,'color',cc(48,:))
plot(t_sec((gaussCut(end)+2):end),velo_forelimb_left.intruder.norm-26,'linewidth',2,'color',cc(64,:))

disp('left paw velocity complete')

%right paw velocity
%get right paw x-y
forelimb_right_pos.home.raw(:,:) = reshape(ri_tracks_matrix(:,forelimb_right,:,home_mouse),[size(ri_tracks_matrix,1) 2]);
forelimb_right_pos.intruder.raw(:,:) = reshape(ri_tracks_matrix(:,forelimb_right,:,intruder),[size(ri_tracks_matrix,1) 2]);
forelimb_right_pos.home.interp = -1.*ones(size(forelimb_right_pos.home.raw));
forelimb_right_pos.intruder.interp = -1.*ones(size(forelimb_right_pos.intruder.raw));

%interpolate nans
for ii = 1:2
    nans.forelimb_right.home(:,ii) = isnan(forelimb_right_pos.home.raw(:,ii));
    nans.forelimb_right.intruder(:,ii) = isnan(forelimb_right_pos.intruder.raw(:,ii));
    forelimb_right_pos.home.interp(nans.forelimb_right.home(:,ii),ii) = interp1(t(~nans.forelimb_right.home(:,ii)), forelimb_right_pos.home.raw(~nans.forelimb_right.home(:,ii),ii), t(nans.forelimb_right.home(:,ii)));
    forelimb_right_pos.home.interp(~nans.forelimb_right.home(:,ii),ii) = forelimb_right_pos.home.raw(~nans.forelimb_right.home(:,ii),ii);
    forelimb_right_pos.intruder.interp(nans.forelimb_right.intruder(:,ii),ii) = interp1(t(~nans.forelimb_right.intruder(:,ii)), forelimb_right_pos.intruder.raw(~nans.forelimb_right.intruder(:,ii),ii), t(nans.forelimb_right.intruder(:,ii)));
    forelimb_right_pos.intruder.interp(~nans.forelimb_right.intruder(:,ii),ii) = forelimb_right_pos.intruder.raw(~nans.forelimb_right.intruder(:,ii),ii);
end

%velocity
velo_forelimb_right.home.raw = -1.*ones(noFrames_intruder-1,1);
velo_forelimb_right.intruder.raw = -1.*ones(noFrames_intruder-1,1);
for ii = 1:noFrames_intruder-1
    velo_forelimb_right.home.raw(ii) = norm(forelimb_right_pos.home.interp(ii+1,:)-forelimb_right_pos.home.interp(ii,:)); %velocity in pixels/frame
    velo_forelimb_right.intruder.raw(ii) = norm(forelimb_right_pos.intruder.interp(ii+1,:)-forelimb_right_pos.intruder.interp(ii,:)); %velocity in pixels/frame
end

%interpolate extreme values
upper.velo_forelimb_right.home = quantile(velo_forelimb_right.home.raw,.97);
lower.velo_forelimb_right.home = quantile(velo_forelimb_right.home.raw,.03);
highs.velo_forelimb_right.home = (velo_forelimb_right.home.raw>upper.velo_forelimb_right.home);
lows.velo_forelimb_right.home = (velo_forelimb_right.home.raw<lower.velo_forelimb_right.home);
velo_forelimb_right.home.raw(highs.velo_forelimb_right.home) = interp1(t(~highs.velo_forelimb_right.home), velo_forelimb_right.home.raw(~highs.velo_forelimb_right.home), t(highs.velo_forelimb_right.home));
velo_forelimb_right.home.raw(lows.velo_forelimb_right.home) = interp1(t(~lows.velo_forelimb_right.home), velo_forelimb_right.home.raw(~lows.velo_forelimb_right.home), t(lows.velo_forelimb_right.home));
upper.velo_forelimb_right.intruder = quantile(velo_forelimb_right.intruder.raw,.97);
lower.velo_forelimb_right.intruder = quantile(velo_forelimb_right.intruder.raw,.03);
highs.velo_forelimb_right.intruder = (velo_forelimb_right.intruder.raw>upper.velo_forelimb_right.intruder);
lows.velo_forelimb_right.intruder = (velo_forelimb_right.intruder.raw<lower.velo_forelimb_right.intruder);
velo_forelimb_right.intruder.raw(highs.velo_forelimb_right.intruder) = interp1(t(~highs.velo_forelimb_right.intruder), velo_forelimb_right.intruder.raw(~highs.velo_forelimb_right.intruder), t(highs.velo_forelimb_right.intruder));
velo_forelimb_right.intruder.raw(lows.velo_forelimb_right.intruder) = interp1(t(~lows.velo_forelimb_right.intruder), velo_forelimb_right.intruder.raw(~lows.velo_forelimb_right.intruder), t(lows.velo_forelimb_right.intruder));

%filter
velo_forelimb_right.home.filt = medfilt1(velo_forelimb_right.home.raw,10);
velo_forelimb_right.home.gauss = filter(w,1,velo_forelimb_right.home.filt);
velo_forelimb_right.home.gauss(1:gaussCut(end)) = [];
velo_forelimb_right.home.norm = rescale(velo_forelimb_right.home.gauss,-1,1);
velo_forelimb_right.intruder.filt = medfilt1(velo_forelimb_right.intruder.raw,10);
velo_forelimb_right.intruder.gauss = filter(w,1,velo_forelimb_right.intruder.filt);
velo_forelimb_right.intruder.gauss(1:gaussCut(end)) = [];
velo_forelimb_right.intruder.norm = rescale(velo_forelimb_right.intruder.gauss,-1,1);

%plot
plot(t_sec((gaussCut(end)+2):end),velo_forelimb_right.home.norm-28,'linewidth',2,'color',cc(80,:))
plot(t_sec((gaussCut(end)+2):end),velo_forelimb_right.intruder.norm-30,'linewidth',2,'color',cc(96,:))

disp('right paw velocity complete')

% orientation_mouse (nose-to-trunk angle)
% defined as dot product of [trunk_to_nose dot nose_to_othertrunk] divided by product of vector length [trunk_to_nose x nose_to_othertrunk]
% 1 is facing towards, -1 is facing away, 0 is right angle
orientation_mouse.home.raw = -1.*ones(noFrames_intruder,1); %orientation_mouse of home mouse to intruder mouse
orientation_mouse.intruder.raw = -1.*ones(noFrames_intruder,1); %orientation_mouse of intruder mouse to mouse mouse

%neck x-y
neck_pos.home.raw(:,:) = reshape(ri_tracks_matrix(:,neck,:,home_mouse),[size(ri_tracks_matrix,1) 2]);
neck_pos.intruder.raw(:,:) = reshape(ri_tracks_matrix(:,neck,:,intruder),[size(ri_tracks_matrix,1) 2]);
neck_pos.home.interp = -1.*ones(size(neck_pos.home.raw));
neck_pos.intruder.interp = -1.*ones(size(neck_pos.intruder.raw));

%interpolate nans
for ii = 1:2
    nans.neck.home(:,ii) = isnan(neck_pos.home.raw(:,ii));
    nans.neck.intruder(:,ii) = isnan(neck_pos.intruder.raw(:,ii));
    neck_pos.home.interp(nans.neck.home(:,ii),ii) = interp1(t(~nans.neck.home(:,ii)), neck_pos.home.raw(~nans.neck.home(:,ii),ii), t(nans.neck.home(:,ii)));
    neck_pos.home.interp(~nans.neck.home(:,ii),ii) = neck_pos.home.raw(~nans.neck.home(:,ii),ii);
    neck_pos.intruder.interp(nans.neck.intruder(:,ii),ii) = interp1(t(~nans.neck.intruder(:,ii)), neck_pos.intruder.raw(~nans.neck.intruder(:,ii),ii), t(nans.neck.intruder(:,ii)));
    neck_pos.intruder.interp(~nans.neck.intruder(:,ii),ii) = neck_pos.intruder.raw(~nans.neck.intruder(:,ii),ii);
end

%get orientation_mouse
for ii = 1:noFrames_intruder
    orientation_mouse.home.raw(ii) =dot(trunk_pos.home.interp(ii,:)-nose_pos.home.interp(ii,:),nose_pos.home.interp(ii,:)-trunk_pos.intruder.interp(ii,:))...
        /(norm(trunk_pos.home.interp(ii,:)-nose_pos.home.interp(ii,:))*norm(nose_pos.home.interp(ii,:)-trunk_pos.intruder.interp(ii,:)));
    orientation_mouse.intruder.raw(ii) = dot(trunk_pos.intruder.interp(ii,:)-nose_pos.intruder.interp(ii,:),nose_pos.intruder.interp(ii,:)-trunk_pos.home.interp(ii,:))...
        /(norm(trunk_pos.intruder.interp(ii,:)-nose_pos.intruder.interp(ii,:))*norm(nose_pos.intruder.interp(ii,:)-trunk_pos.home.interp(ii,:)));
end

%interpolate extreme values
upper.orientation_mouse.home = quantile(orientation_mouse.home.raw,.97);
lower.orientation_mouse.home = quantile(orientation_mouse.home.raw,.03);
highs.orientation_mouse.home = (orientation_mouse.home.raw>upper.orientation_mouse.home);
lows.orientation_mouse.home = (orientation_mouse.home.raw<lower.orientation_mouse.home);
orientation_mouse.home.raw(highs.orientation_mouse.home) = interp1(t(~highs.orientation_mouse.home), orientation_mouse.home.raw(~highs.orientation_mouse.home), t(highs.orientation_mouse.home));
orientation_mouse.home.raw(lows.orientation_mouse.home) = interp1(t(~lows.orientation_mouse.home), orientation_mouse.home.raw(~lows.orientation_mouse.home), t(lows.orientation_mouse.home));
upper.orientation_mouse.intruder = quantile(orientation_mouse.intruder.raw,.97);
lower.orientation_mouse.intruder = quantile(orientation_mouse.intruder.raw,.03);
highs.orientation_mouse.intruder = (orientation_mouse.intruder.raw>upper.orientation_mouse.intruder);
lows.orientation_mouse.intruder = (orientation_mouse.intruder.raw<lower.orientation_mouse.intruder);
orientation_mouse.intruder.raw(highs.orientation_mouse.intruder) = interp1(t(~highs.orientation_mouse.intruder), orientation_mouse.intruder.raw(~highs.orientation_mouse.intruder), t(highs.orientation_mouse.intruder));
orientation_mouse.intruder.raw(lows.orientation_mouse.intruder) = interp1(t(~lows.orientation_mouse.intruder), orientation_mouse.intruder.raw(~lows.orientation_mouse.intruder), t(lows.orientation_mouse.intruder));

%filter
orientation_mouse.home.filt = medfilt1(orientation_mouse.home.raw,10);
orientation_mouse.home.norm = filter(w,1,orientation_mouse.home.filt);
orientation_mouse.home.norm(1:gaussCut(end)+1) = [];
orientation_mouse.intruder.filt = medfilt1(orientation_mouse.intruder.raw,10);
orientation_mouse.intruder.norm = filter(w,1,orientation_mouse.intruder.filt);
orientation_mouse.intruder.norm(1:gaussCut(end)+1) = [];

%plot
plot(t_sec((gaussCut(end)+2):end),orientation_mouse.home.norm-32,'linewidth',2,'color',cc(1,:))
plot(t_sec((gaussCut(end)+2):end),orientation_mouse.intruder.norm-34,'linewidth',2,'color',cc(1,:))

disp('mouse orientation complete')

% forelimb orientation_mouse (tti-neck-to-forelimb angle)
% 1 is facing towards, -1 backwards, 0 is right angle
orientation_forelimb_left.home.raw = -1.*ones(noFrames_intruder,1); %orientation_mouse of home mouse to intruder mouse
orientation_forelimb_left.intruder.raw = -1.*ones(noFrames_intruder,1); %orientation_mouse of intruder mouse to mouse mouse

%get orientation
for ii = 1:noFrames_intruder
    orientation_forelimb_left.home.raw(ii) =dot(tti_pos.home.interp(ii,:)-neck_pos.home.interp(ii,:),neck_pos.home.interp(ii,:)-forelimb_left_pos.home.interp(ii,:))...
        /(norm(tti_pos.home.interp(ii,:)-neck_pos.home.interp(ii,:))*norm(neck_pos.home.interp(ii,:)-forelimb_left_pos.home.interp(ii,:)));
    orientation_forelimb_left.intruder.raw(ii) =dot(tti_pos.intruder.interp(ii,:)-neck_pos.intruder.interp(ii,:),neck_pos.intruder.interp(ii,:)-forelimb_left_pos.intruder.interp(ii,:))...
        /(norm(tti_pos.intruder.interp(ii,:)-neck_pos.intruder.interp(ii,:))*norm(neck_pos.intruder.interp(ii,:)-forelimb_left_pos.intruder.interp(ii,:)));
end

%interpolate extreme values
upper.orientation_forelimb_left.home = quantile(orientation_forelimb_left.home.raw,.97);
lower.orientation_forelimb_left.home = quantile(orientation_forelimb_left.home.raw,.03);
highs.orientation_forelimb_left.home = (orientation_forelimb_left.home.raw>upper.orientation_forelimb_left.home);
lows.orientation_forelimb_left.home = (orientation_forelimb_left.home.raw<lower.orientation_forelimb_left.home);
orientation_forelimb_left.home.raw(highs.orientation_forelimb_left.home) = interp1(t(~highs.orientation_forelimb_left.home), orientation_forelimb_left.home.raw(~highs.orientation_forelimb_left.home), t(highs.orientation_forelimb_left.home));
orientation_forelimb_left.home.raw(lows.orientation_forelimb_left.home) = interp1(t(~lows.orientation_forelimb_left.home), orientation_forelimb_left.home.raw(~lows.orientation_forelimb_left.home), t(lows.orientation_forelimb_left.home));
upper.orientation_forelimb_left.intruder = quantile(orientation_forelimb_left.intruder.raw,.97);
lower.orientation_forelimb_left.intruder = quantile(orientation_forelimb_left.intruder.raw,.03);
highs.orientation_forelimb_left.intruder = (orientation_forelimb_left.intruder.raw>upper.orientation_forelimb_left.intruder);
lows.orientation_forelimb_left.intruder = (orientation_forelimb_left.intruder.raw<lower.orientation_forelimb_left.intruder);
orientation_forelimb_left.intruder.raw(highs.orientation_forelimb_left.intruder) = interp1(t(~highs.orientation_forelimb_left.intruder), orientation_forelimb_left.intruder.raw(~highs.orientation_forelimb_left.intruder), t(highs.orientation_forelimb_left.intruder));
orientation_forelimb_left.intruder.raw(lows.orientation_forelimb_left.intruder) = interp1(t(~lows.orientation_forelimb_left.intruder), orientation_forelimb_left.intruder.raw(~lows.orientation_forelimb_left.intruder), t(lows.orientation_forelimb_left.intruder));

%filter
orientation_forelimb_left.home.filt = medfilt1(orientation_forelimb_left.home.raw,10);
orientation_forelimb_left.home.norm = filter(w,1,orientation_forelimb_left.home.filt);
orientation_forelimb_left.home.norm(1:gaussCut(end)+1) = [];
orientation_forelimb_left.intruder.filt = medfilt1(orientation_forelimb_left.intruder.raw,10);
orientation_forelimb_left.intruder.norm = filter(w,1,orientation_forelimb_left.intruder.filt);
orientation_forelimb_left.intruder.norm(1:gaussCut(end)+1) = [];

%plot
plot(t_sec((gaussCut(end)+2):end),orientation_forelimb_left.home.norm-36,'linewidth',2,'color',cc(112,:))
plot(t_sec((gaussCut(end)+2):end),orientation_forelimb_left.intruder.norm-38,'linewidth',2,'color',cc(144,:))

disp('left forelimb orientation complete')

% body orientation (tti-trunk-to-nose angle)
% 1 is facing towards, -1 backwards, 0 is right angle
orientation_body.home.raw = -1.*ones(noFrames_intruder,1); %orientation_mouse of home mouse to intruder mouse
orientation_body.intruder.raw = -1.*ones(noFrames_intruder,1); %orientation_mouse of intruder mouse to mouse mouse

%get orientation
for ii = 1:noFrames_intruder
    orientation_body.home.raw(ii) =dot(tti_pos.home.interp(ii,:)-trunk_pos.home.interp(ii,:),trunk_pos.home.interp(ii,:)-nose_pos.home.interp(ii,:))...
        /(norm(tti_pos.home.interp(ii,:)-trunk_pos.home.interp(ii,:))*norm(trunk_pos.home.interp(ii,:)-nose_pos.home.interp(ii,:)));
    orientation_body.intruder.raw(ii) =dot(tti_pos.intruder.interp(ii,:)-trunk_pos.intruder.interp(ii,:),trunk_pos.intruder.interp(ii,:)-nose_pos.intruder.interp(ii,:))...
        /(norm(tti_pos.intruder.interp(ii,:)-trunk_pos.intruder.interp(ii,:))*norm(trunk_pos.intruder.interp(ii,:)-nose_pos.intruder.interp(ii,:)));
end

%interpolate extreme values
upper.orientation_body.home = quantile(orientation_body.home.raw,.97);
lower.orientation_body.home = quantile(orientation_body.home.raw,.03);
highs.orientation_body.home = (orientation_body.home.raw>upper.orientation_body.home);
lows.orientation_body.home = (orientation_body.home.raw<lower.orientation_body.home);
orientation_body.home.raw(highs.orientation_body.home) = interp1(t(~highs.orientation_body.home), orientation_body.home.raw(~highs.orientation_body.home), t(highs.orientation_body.home));
orientation_body.home.raw(lows.orientation_body.home) = interp1(t(~lows.orientation_body.home), orientation_body.home.raw(~lows.orientation_body.home), t(lows.orientation_body.home));
upper.orientation_body.intruder = quantile(orientation_body.intruder.raw,.97);
lower.orientation_body.intruder = quantile(orientation_body.intruder.raw,.03);
highs.orientation_body.intruder = (orientation_body.intruder.raw>upper.orientation_body.intruder);
lows.orientation_body.intruder = (orientation_body.intruder.raw<lower.orientation_body.intruder);
orientation_body.intruder.raw(highs.orientation_body.intruder) = interp1(t(~highs.orientation_body.intruder), orientation_body.intruder.raw(~highs.orientation_body.intruder), t(highs.orientation_body.intruder));
orientation_body.intruder.raw(lows.orientation_body.intruder) = interp1(t(~lows.orientation_body.intruder), orientation_body.intruder.raw(~lows.orientation_body.intruder), t(lows.orientation_body.intruder));

%filter
orientation_body.home.filt = medfilt1(orientation_body.home.raw,10);
orientation_body.home.norm = filter(w,1,orientation_body.home.filt);
orientation_body.home.norm(1:gaussCut(end)+1) = [];
orientation_body.intruder.filt = medfilt1(orientation_body.intruder.raw,10);
orientation_body.intruder.norm = filter(w,1,orientation_body.intruder.filt);
orientation_body.intruder.norm(1:gaussCut(end)+1) = [];

%plot
plot(t_sec((gaussCut(end)+2):end),orientation_body.home.norm-40,'linewidth',2,'color',cc(128,:))
plot(t_sec((gaussCut(end)+2):end),orientation_body.intruder.norm-42,'linewidth',2,'color',cc(128,:))

disp('body orientation complete')

% forelimb orientation_mouse (tti-neck-to-forelimb angle)
% 1 is facing towards, -1 backwards, 0 is right angle
orientation_forelimb_right.home.raw = -1.*ones(noFrames_intruder,1); %orientation_mouse of home mouse to intruder mouse
orientation_forelimb_right.intruder.raw = -1.*ones(noFrames_intruder,1); %orientation_mouse of intruder mouse to mouse mouse

%get orientation
for ii = 1:noFrames_intruder
    orientation_forelimb_right.home.raw(ii) =dot(tti_pos.home.interp(ii,:)-neck_pos.home.interp(ii,:),neck_pos.home.interp(ii,:)-forelimb_right_pos.home.interp(ii,:))...
        /(norm(tti_pos.home.interp(ii,:)-neck_pos.home.interp(ii,:))*norm(neck_pos.home.interp(ii,:)-forelimb_right_pos.home.interp(ii,:)));
    orientation_forelimb_right.intruder.raw(ii) =dot(tti_pos.intruder.interp(ii,:)-neck_pos.intruder.interp(ii,:),neck_pos.intruder.interp(ii,:)-forelimb_right_pos.intruder.interp(ii,:))...
        /(norm(tti_pos.intruder.interp(ii,:)-neck_pos.intruder.interp(ii,:))*norm(neck_pos.intruder.interp(ii,:)-forelimb_right_pos.intruder.interp(ii,:)));
end

%interpolate extreme values
upper.orientation_forelimb_right.home = quantile(orientation_forelimb_right.home.raw,.97);
lower.orientation_forelimb_right.home = quantile(orientation_forelimb_right.home.raw,.03);
highs.orientation_forelimb_right.home = (orientation_forelimb_right.home.raw>upper.orientation_forelimb_right.home);
lows.orientation_forelimb_right.home = (orientation_forelimb_right.home.raw<lower.orientation_forelimb_right.home);
orientation_forelimb_right.home.raw(highs.orientation_forelimb_right.home) = interp1(t(~highs.orientation_forelimb_right.home), orientation_forelimb_right.home.raw(~highs.orientation_forelimb_right.home), t(highs.orientation_forelimb_right.home));
orientation_forelimb_right.home.raw(lows.orientation_forelimb_right.home) = interp1(t(~lows.orientation_forelimb_right.home), orientation_forelimb_right.home.raw(~lows.orientation_forelimb_right.home), t(lows.orientation_forelimb_right.home));
upper.orientation_forelimb_right.intruder = quantile(orientation_forelimb_right.intruder.raw,.97);
lower.orientation_forelimb_right.intruder = quantile(orientation_forelimb_right.intruder.raw,.03);
highs.orientation_forelimb_right.intruder = (orientation_forelimb_right.intruder.raw>upper.orientation_forelimb_right.intruder);
lows.orientation_forelimb_right.intruder = (orientation_forelimb_right.intruder.raw<lower.orientation_forelimb_right.intruder);
orientation_forelimb_right.intruder.raw(highs.orientation_forelimb_right.intruder) = interp1(t(~highs.orientation_forelimb_right.intruder), orientation_forelimb_right.intruder.raw(~highs.orientation_forelimb_right.intruder), t(highs.orientation_forelimb_right.intruder));
orientation_forelimb_right.intruder.raw(lows.orientation_forelimb_right.intruder) = interp1(t(~lows.orientation_forelimb_right.intruder), orientation_forelimb_right.intruder.raw(~lows.orientation_forelimb_right.intruder), t(lows.orientation_forelimb_right.intruder));

%filter
orientation_forelimb_right.home.filt = medfilt1(orientation_forelimb_right.home.raw,10);
orientation_forelimb_right.home.norm = filter(w,1,orientation_forelimb_right.home.filt);
orientation_forelimb_right.home.norm(1:gaussCut(end)+1) = [];
orientation_forelimb_right.intruder.filt = medfilt1(orientation_forelimb_right.intruder.raw,10);
orientation_forelimb_right.intruder.norm = filter(w,1,orientation_forelimb_right.intruder.filt);
orientation_forelimb_right.intruder.norm(1:gaussCut(end)+1) = [];

%plot
plot(t_sec((gaussCut(end)+2):end),orientation_forelimb_right.home.norm-44,'linewidth',2,'color',cc(160,:))
plot(t_sec((gaussCut(end)+2):end),orientation_forelimb_right.intruder.norm-46,'linewidth',2,'color',cc(176,:))

disp('right forelimb orientation complete')

%nose to neck distance
% home_nose, intruder_neck

%nose to neck distance
dist_nneck.raw = -1.*ones(noFrames_intruder,1);
dist_nneck.norm = -1.*ones(noFrames_intruder,1);
for ii = 1:noFrames_intruder
    dist_nneck.raw(ii) = norm(nose_pos.home.interp(ii,:)-neck_pos.intruder.interp(ii,:)); %nose-neck distance in pixels
end

%interpolate extreme values
upper.dist_nneck = quantile(dist_nneck.raw,.97);
lower.dist_nneck = quantile(dist_nneck.raw,.03);
highs.dist_nneck = (dist_nneck.raw>upper.dist_nneck);
lows.dist_nneck = (dist_nneck.raw<lower.dist_nneck);
dist_nneck.raw(highs.dist_nneck) = interp1(t(~highs.dist_nneck), dist_nneck.raw(~highs.dist_nneck), t(highs.dist_nneck));
dist_nneck.raw(lows.dist_nneck) = interp1(t(~lows.dist_nneck), dist_nneck.raw(~lows.dist_nneck), t(lows.dist_nneck));

%filter
dist_nneck.filt = medfilt1(dist_nneck.raw,10);
dist_nneck.gauss = filter(w,1,dist_nneck.filt);
dist_nneck.gauss(1:gaussCut(end)+1) = [];
dist_nneck.norm = rescale(dist_nneck.gauss,-1,1);

plot(t_sec((gaussCut(end)+2):end),dist_nneck.norm-48,'linewidth',2,'color',cc(128,:))
disp('nose to neck distance complete')

%nose to trunk distance
% home_nose, intruder_trunk

%nose to trunk distance
dist_ntrunk.raw = -1.*ones(noFrames_intruder,1);
dist_ntrunk.norm = -1.*ones(noFrames_intruder,1);
for ii = 1:noFrames_intruder
    dist_ntrunk.raw(ii) = norm(nose_pos.home.interp(ii,:)-trunk_pos.intruder.interp(ii,:)); %nose-trunk distance in pixels
end

%interpolate extreme values
upper.dist_ntrunk = quantile(dist_ntrunk.raw,.97);
lower.dist_ntrunk = quantile(dist_ntrunk.raw,.03);
highs.dist_ntrunk = (dist_ntrunk.raw>upper.dist_ntrunk);
lows.dist_ntrunk = (dist_ntrunk.raw<lower.dist_ntrunk);
dist_ntrunk.raw(highs.dist_ntrunk) = interp1(t(~highs.dist_ntrunk), dist_ntrunk.raw(~highs.dist_ntrunk), t(highs.dist_ntrunk));
dist_ntrunk.raw(lows.dist_ntrunk) = interp1(t(~lows.dist_ntrunk), dist_ntrunk.raw(~lows.dist_ntrunk), t(lows.dist_ntrunk));

%filter
dist_ntrunk.filt = medfilt1(dist_ntrunk.raw,10);
dist_ntrunk.gauss = filter(w,1,dist_ntrunk.filt);
dist_ntrunk.gauss(1:gaussCut(end)+1) = [];
dist_ntrunk.norm = rescale(dist_ntrunk.gauss,-1,1);

plot(t_sec((gaussCut(end)+2):end),dist_ntrunk.norm-50,'linewidth',2,'color',cc(128,:))
disp('nose to trunk distance complete')

% tail_base orientation (trunk-tti-to-tail0 angle)
% 1 is facing towards, -1 backwards, 0 is right angle
orientation_tailbase.home.raw = -1.*ones(noFrames_intruder,1); %orientation_mouse of home mouse to intruder mouse
orientation_tailbase.intruder.raw = -1.*ones(noFrames_intruder,1); %orientation_mouse of intruder mouse to mouse mouse

%tail0 x-y
tail0_pos.home.raw(:,:) = reshape(ri_tracks_matrix(:,tail0,:,home_mouse),[size(ri_tracks_matrix,1) 2]);
tail0_pos.intruder.raw(:,:) = reshape(ri_tracks_matrix(:,tail0,:,intruder),[size(ri_tracks_matrix,1) 2]);
tail0_pos.home.interp = -1.*ones(size(tail0_pos.home.raw));
tail0_pos.intruder.interp = -1.*ones(size(tail0_pos.intruder.raw));

%interpolate nans
for ii = 1:2
    nans.tail0.home(:,ii) = isnan(tail0_pos.home.raw(:,ii));
    nans.tail0.intruder(:,ii) = isnan(tail0_pos.intruder.raw(:,ii));
    tail0_pos.home.interp(nans.tail0.home(:,ii),ii) = interp1(t(~nans.tail0.home(:,ii)), tail0_pos.home.raw(~nans.tail0.home(:,ii),ii), t(nans.tail0.home(:,ii)));
    tail0_pos.home.interp(~nans.tail0.home(:,ii),ii) = tail0_pos.home.raw(~nans.tail0.home(:,ii),ii);
    tail0_pos.intruder.interp(nans.tail0.intruder(:,ii),ii) = interp1(t(~nans.tail0.intruder(:,ii)), tail0_pos.intruder.raw(~nans.tail0.intruder(:,ii),ii), t(nans.tail0.intruder(:,ii)));
    tail0_pos.intruder.interp(~nans.tail0.intruder(:,ii),ii) = tail0_pos.intruder.raw(~nans.tail0.intruder(:,ii),ii);
end

%get orientation
for ii = 1:noFrames_intruder
    orientation_tailbase.home.raw(ii) =dot(trunk_pos.home.interp(ii,:)-tti_pos.home.interp(ii,:),tti_pos.home.interp(ii,:)-tail0_pos.home.interp(ii,:))...
        /(norm(trunk_pos.home.interp(ii,:)-tti_pos.home.interp(ii,:))*norm(tti_pos.home.interp(ii,:)-tail0_pos.home.interp(ii,:)));
    orientation_tailbase.intruder.raw(ii) =dot(trunk_pos.intruder.interp(ii,:)-tti_pos.intruder.interp(ii,:),tti_pos.intruder.interp(ii,:)-tail0_pos.intruder.interp(ii,:))...
        /(norm(trunk_pos.intruder.interp(ii,:)-tti_pos.intruder.interp(ii,:))*norm(tti_pos.intruder.interp(ii,:)-tail0_pos.intruder.interp(ii,:)));
end

%interpolate extreme values
upper.orientation_tailbase.home = quantile(orientation_tailbase.home.raw,.97);
lower.orientation_tailbase.home = quantile(orientation_tailbase.home.raw,.03);
highs.orientation_tailbase.home = (orientation_tailbase.home.raw>upper.orientation_tailbase.home);
lows.orientation_tailbase.home = (orientation_tailbase.home.raw<lower.orientation_tailbase.home);
orientation_tailbase.home.raw(highs.orientation_tailbase.home) = interp1(t(~highs.orientation_tailbase.home), orientation_tailbase.home.raw(~highs.orientation_tailbase.home), t(highs.orientation_tailbase.home));
orientation_tailbase.home.raw(lows.orientation_tailbase.home) = interp1(t(~lows.orientation_tailbase.home), orientation_tailbase.home.raw(~lows.orientation_tailbase.home), t(lows.orientation_tailbase.home));
upper.orientation_tailbase.intruder = quantile(orientation_tailbase.intruder.raw,.97);
lower.orientation_tailbase.intruder = quantile(orientation_tailbase.intruder.raw,.03);
highs.orientation_tailbase.intruder = (orientation_tailbase.intruder.raw>upper.orientation_tailbase.intruder);
lows.orientation_tailbase.intruder = (orientation_tailbase.intruder.raw<lower.orientation_tailbase.intruder);
orientation_tailbase.intruder.raw(highs.orientation_tailbase.intruder) = interp1(t(~highs.orientation_tailbase.intruder), orientation_tailbase.intruder.raw(~highs.orientation_tailbase.intruder), t(highs.orientation_tailbase.intruder));
orientation_tailbase.intruder.raw(lows.orientation_tailbase.intruder) = interp1(t(~lows.orientation_tailbase.intruder), orientation_tailbase.intruder.raw(~lows.orientation_tailbase.intruder), t(lows.orientation_tailbase.intruder));

%filter
orientation_tailbase.home.filt = medfilt1(orientation_tailbase.home.raw,10);
orientation_tailbase.home.norm = filter(w,1,orientation_tailbase.home.filt);
orientation_tailbase.home.norm(1:gaussCut(end)+1) = [];
orientation_tailbase.intruder.filt = medfilt1(orientation_tailbase.intruder.raw,10);
orientation_tailbase.intruder.norm = filter(w,1,orientation_tailbase.intruder.filt);
orientation_tailbase.intruder.norm(1:gaussCut(end)+1) = [];

%plot
plot(t_sec((gaussCut(end)+2):end),orientation_tailbase.home.norm-52,'linewidth',2,'color',cc(192,:))
plot(t_sec((gaussCut(end)+2):end),orientation_tailbase.intruder.norm-60,'linewidth',2,'color',cc(208,:))

disp('tail_base orientation complete')

% tail1 orientation (tti-to-tail0-to-tail1 angle)
% 1 is facing towards, -1 backwards, 0 is right angle
orientation_tail1.home.raw = -1.*ones(noFrames_intruder,1); %orientation_mouse of home mouse to intruder mouse
orientation_tail1.intruder.raw = -1.*ones(noFrames_intruder,1); %orientation_mouse of intruder mouse to mouse mouse

%tail1 x-y
tail1_pos.home.raw(:,:) = reshape(ri_tracks_matrix(:,tail1,:,home_mouse),[size(ri_tracks_matrix,1) 2]);
tail1_pos.intruder.raw(:,:) = reshape(ri_tracks_matrix(:,tail1,:,intruder),[size(ri_tracks_matrix,1) 2]);
tail1_pos.home.interp = -1.*ones(size(tail1_pos.home.raw));
tail1_pos.intruder.interp = -1.*ones(size(tail1_pos.intruder.raw));

%interpolate nans
for ii = 1:2
    nans.tail1.home(:,ii) = isnan(tail1_pos.home.raw(:,ii));
    nans.tail1.intruder(:,ii) = isnan(tail1_pos.intruder.raw(:,ii));
    tail1_pos.home.interp(nans.tail1.home(:,ii),ii) = interp1(t(~nans.tail1.home(:,ii)), tail1_pos.home.raw(~nans.tail1.home(:,ii),ii), t(nans.tail1.home(:,ii)));
    tail1_pos.home.interp(~nans.tail1.home(:,ii),ii) = tail1_pos.home.raw(~nans.tail1.home(:,ii),ii);
    tail1_pos.intruder.interp(nans.tail1.intruder(:,ii),ii) = interp1(t(~nans.tail1.intruder(:,ii)), tail1_pos.intruder.raw(~nans.tail1.intruder(:,ii),ii), t(nans.tail1.intruder(:,ii)));
    tail1_pos.intruder.interp(~nans.tail1.intruder(:,ii),ii) = tail1_pos.intruder.raw(~nans.tail1.intruder(:,ii),ii);
end

%get orientation
for ii = 1:noFrames_intruder
    orientation_tail1.home.raw(ii) =dot(tti_pos.home.interp(ii,:)-tail0_pos.home.interp(ii,:),tail0_pos.home.interp(ii,:)-tail1_pos.home.interp(ii,:))...
        /(norm(tti_pos.home.interp(ii,:)-tail0_pos.home.interp(ii,:))*norm(tail0_pos.home.interp(ii,:)-tail1_pos.home.interp(ii,:)));
    orientation_tail1.intruder.raw(ii) =dot(tti_pos.intruder.interp(ii,:)-tail0_pos.intruder.interp(ii,:),tail0_pos.intruder.interp(ii,:)-tail1_pos.intruder.interp(ii,:))...
        /(norm(tti_pos.intruder.interp(ii,:)-tail0_pos.intruder.interp(ii,:))*norm(tail0_pos.intruder.interp(ii,:)-tail1_pos.intruder.interp(ii,:)));
end

%interpolate extreme values
upper.orientation_tail1.home = quantile(orientation_tail1.home.raw,.97);
lower.orientation_tail1.home = quantile(orientation_tail1.home.raw,.03);
highs.orientation_tail1.home = (orientation_tail1.home.raw>upper.orientation_tail1.home);
lows.orientation_tail1.home = (orientation_tail1.home.raw<lower.orientation_tail1.home);
orientation_tail1.home.raw(highs.orientation_tail1.home) = interp1(t(~highs.orientation_tail1.home), orientation_tail1.home.raw(~highs.orientation_tail1.home), t(highs.orientation_tail1.home));
orientation_tail1.home.raw(lows.orientation_tail1.home) = interp1(t(~lows.orientation_tail1.home), orientation_tail1.home.raw(~lows.orientation_tail1.home), t(lows.orientation_tail1.home));
upper.orientation_tail1.intruder = quantile(orientation_tail1.intruder.raw,.97);
lower.orientation_tail1.intruder = quantile(orientation_tail1.intruder.raw,.03);
highs.orientation_tail1.intruder = (orientation_tail1.intruder.raw>upper.orientation_tail1.intruder);
lows.orientation_tail1.intruder = (orientation_tail1.intruder.raw<lower.orientation_tail1.intruder);
orientation_tail1.intruder.raw(highs.orientation_tail1.intruder) = interp1(t(~highs.orientation_tail1.intruder), orientation_tail1.intruder.raw(~highs.orientation_tail1.intruder), t(highs.orientation_tail1.intruder));
orientation_tail1.intruder.raw(lows.orientation_tail1.intruder) = interp1(t(~lows.orientation_tail1.intruder), orientation_tail1.intruder.raw(~lows.orientation_tail1.intruder), t(lows.orientation_tail1.intruder));

%filter
orientation_tail1.home.filt = medfilt1(orientation_tail1.home.raw,10);
orientation_tail1.home.norm = filter(w,1,orientation_tail1.home.filt);
orientation_tail1.home.norm(1:gaussCut(end)+1) = [];
orientation_tail1.intruder.filt = medfilt1(orientation_tail1.intruder.raw,10);
orientation_tail1.intruder.norm = filter(w,1,orientation_tail1.intruder.filt);
orientation_tail1.intruder.norm(1:gaussCut(end)+1) = [];

%plot
plot(t_sec((gaussCut(end)+2):end),orientation_tail1.home.norm-54,'linewidth',2,'color',cc(192,:))
plot(t_sec((gaussCut(end)+2):end),orientation_tail1.intruder.norm-62,'linewidth',2,'color',cc(208,:))

disp('tail_node_one orientation complete')

% tail2 orientation (tail0-to-tail1-to-tail2 angle)
% 1 is facing towards, -1 backwards, 0 is right angle
orientation_tail2.home.raw = -1.*ones(noFrames_intruder,1); %orientation_mouse of home mouse to intruder mouse
orientation_tail2.intruder.raw = -1.*ones(noFrames_intruder,1); %orientation_mouse of intruder mouse to mouse mouse

%tail2 x-y
tail2_pos.home.raw(:,:) = reshape(ri_tracks_matrix(:,tail2,:,home_mouse),[size(ri_tracks_matrix,1) 2]);
tail2_pos.intruder.raw(:,:) = reshape(ri_tracks_matrix(:,tail2,:,intruder),[size(ri_tracks_matrix,1) 2]);
tail2_pos.home.interp = -1.*ones(size(tail2_pos.home.raw));
tail2_pos.intruder.interp = -1.*ones(size(tail2_pos.intruder.raw));

%interpolate nans
for ii = 1:2
    nans.tail2.home(:,ii) = isnan(tail2_pos.home.raw(:,ii));
    nans.tail2.intruder(:,ii) = isnan(tail2_pos.intruder.raw(:,ii));
    tail2_pos.home.interp(nans.tail2.home(:,ii),ii) = interp1(t(~nans.tail2.home(:,ii)), tail2_pos.home.raw(~nans.tail2.home(:,ii),ii), t(nans.tail2.home(:,ii)));
    tail2_pos.home.interp(~nans.tail2.home(:,ii),ii) = tail2_pos.home.raw(~nans.tail2.home(:,ii),ii);
    tail2_pos.intruder.interp(nans.tail2.intruder(:,ii),ii) = interp1(t(~nans.tail2.intruder(:,ii)), tail2_pos.intruder.raw(~nans.tail2.intruder(:,ii),ii), t(nans.tail2.intruder(:,ii)));
    tail2_pos.intruder.interp(~nans.tail2.intruder(:,ii),ii) = tail2_pos.intruder.raw(~nans.tail2.intruder(:,ii),ii);
end

%get orientation
for ii = 1:noFrames_intruder
    orientation_tail2.home.raw(ii) =dot(tail0_pos.home.interp(ii,:)-tail1_pos.home.interp(ii,:),tail1_pos.home.interp(ii,:)-tail2_pos.home.interp(ii,:))...
        /(norm(tail0_pos.home.interp(ii,:)-tail1_pos.home.interp(ii,:))*norm(tail1_pos.home.interp(ii,:)-tail2_pos.home.interp(ii,:)));
    orientation_tail2.intruder.raw(ii) =dot(tail0_pos.intruder.interp(ii,:)-tail1_pos.intruder.interp(ii,:),tail1_pos.intruder.interp(ii,:)-tail2_pos.intruder.interp(ii,:))...
        /(norm(tail0_pos.intruder.interp(ii,:)-tail1_pos.intruder.interp(ii,:))*norm(tail1_pos.intruder.interp(ii,:)-tail2_pos.intruder.interp(ii,:)));
end

%interpolate extreme values
upper.orientation_tail2.home = quantile(orientation_tail2.home.raw,.97);
lower.orientation_tail2.home = quantile(orientation_tail2.home.raw,.03);
highs.orientation_tail2.home = (orientation_tail2.home.raw>upper.orientation_tail2.home);
lows.orientation_tail2.home = (orientation_tail2.home.raw<lower.orientation_tail2.home);
orientation_tail2.home.raw(highs.orientation_tail2.home) = interp1(t(~highs.orientation_tail2.home), orientation_tail2.home.raw(~highs.orientation_tail2.home), t(highs.orientation_tail2.home));
orientation_tail2.home.raw(lows.orientation_tail2.home) = interp1(t(~lows.orientation_tail2.home), orientation_tail2.home.raw(~lows.orientation_tail2.home), t(lows.orientation_tail2.home));
upper.orientation_tail2.intruder = quantile(orientation_tail2.intruder.raw,.97);
lower.orientation_tail2.intruder = quantile(orientation_tail2.intruder.raw,.03);
highs.orientation_tail2.intruder = (orientation_tail2.intruder.raw>upper.orientation_tail2.intruder);
lows.orientation_tail2.intruder = (orientation_tail2.intruder.raw<lower.orientation_tail2.intruder);
orientation_tail2.intruder.raw(highs.orientation_tail2.intruder) = interp1(t(~highs.orientation_tail2.intruder), orientation_tail2.intruder.raw(~highs.orientation_tail2.intruder), t(highs.orientation_tail2.intruder));
orientation_tail2.intruder.raw(lows.orientation_tail2.intruder) = interp1(t(~lows.orientation_tail2.intruder), orientation_tail2.intruder.raw(~lows.orientation_tail2.intruder), t(lows.orientation_tail2.intruder));

%filter
orientation_tail2.home.filt = medfilt1(orientation_tail2.home.raw,10);
orientation_tail2.home.norm = filter(w,1,orientation_tail2.home.filt);
orientation_tail2.home.norm(1:gaussCut(end)+1) = [];
orientation_tail2.intruder.filt = medfilt1(orientation_tail2.intruder.raw,10);
orientation_tail2.intruder.norm = filter(w,1,orientation_tail2.intruder.filt);
orientation_tail2.intruder.norm(1:gaussCut(end)+1) = [];

%plot
plot(t_sec((gaussCut(end)+2):end),orientation_tail2.home.norm-56,'linewidth',2,'color',cc(192,:))
plot(t_sec((gaussCut(end)+2):end),orientation_tail2.intruder.norm-64,'linewidth',2,'color',cc(208,:))

disp('tail_node_two orientation complete')

% tail3 orientation (tail1-to-tail2-to-tailtip angle)
% 1 is facing towards, -1 backwards, 0 is right angle
orientation_tail3.home.raw = -1.*ones(noFrames_intruder,1); %orientation_mouse of home mouse to intruder mouse
orientation_tail3.intruder.raw = -1.*ones(noFrames_intruder,1); %orientation_mouse of intruder mouse to mouse mouse

%tailtip x-y
tailtip_pos.home.raw(:,:) = reshape(ri_tracks_matrix(:,tailtip,:,home_mouse),[size(ri_tracks_matrix,1) 2]);
tailtip_pos.intruder.raw(:,:) = reshape(ri_tracks_matrix(:,tailtip,:,intruder),[size(ri_tracks_matrix,1) 2]);
tailtip_pos.home.interp = -1.*ones(size(tailtip_pos.home.raw));
tailtip_pos.intruder.interp = -1.*ones(size(tailtip_pos.intruder.raw));

%interpolate nans
for ii = 1:2
    nans.tailtip.home(:,ii) = isnan(tailtip_pos.home.raw(:,ii));
    nans.tailtip.intruder(:,ii) = isnan(tailtip_pos.intruder.raw(:,ii));
    tailtip_pos.home.interp(nans.tailtip.home(:,ii),ii) = interp1(t(~nans.tailtip.home(:,ii)), tailtip_pos.home.raw(~nans.tailtip.home(:,ii),ii), t(nans.tailtip.home(:,ii)));
    tailtip_pos.home.interp(~nans.tailtip.home(:,ii),ii) = tailtip_pos.home.raw(~nans.tailtip.home(:,ii),ii);
    tailtip_pos.intruder.interp(nans.tailtip.intruder(:,ii),ii) = interp1(t(~nans.tailtip.intruder(:,ii)), tailtip_pos.intruder.raw(~nans.tailtip.intruder(:,ii),ii), t(nans.tailtip.intruder(:,ii)));
    tailtip_pos.intruder.interp(~nans.tailtip.intruder(:,ii),ii) = tailtip_pos.intruder.raw(~nans.tailtip.intruder(:,ii),ii);
end

%get orientation
for ii = 1:noFrames_intruder
    orientation_tail3.home.raw(ii) =dot(tail1_pos.home.interp(ii,:)-tail2_pos.home.interp(ii,:),tail2_pos.home.interp(ii,:)-tailtip_pos.home.interp(ii,:))...
        /(norm(tail1_pos.home.interp(ii,:)-tail2_pos.home.interp(ii,:))*norm(tail2_pos.home.interp(ii,:)-tailtip_pos.home.interp(ii,:)));
    orientation_tail3.intruder.raw(ii) =dot(tail1_pos.intruder.interp(ii,:)-tail2_pos.intruder.interp(ii,:),tail2_pos.intruder.interp(ii,:)-tailtip_pos.intruder.interp(ii,:))...
        /(norm(tail1_pos.intruder.interp(ii,:)-tail2_pos.intruder.interp(ii,:))*norm(tail2_pos.intruder.interp(ii,:)-tailtip_pos.intruder.interp(ii,:)));
end

%interpolate extreme values
upper.orientation_tail3.home = quantile(orientation_tail3.home.raw,.97);
lower.orientation_tail3.home = quantile(orientation_tail3.home.raw,.03);
highs.orientation_tail3.home = (orientation_tail3.home.raw>upper.orientation_tail3.home);
lows.orientation_tail3.home = (orientation_tail3.home.raw<lower.orientation_tail3.home);
orientation_tail3.home.raw(highs.orientation_tail3.home) = interp1(t(~highs.orientation_tail3.home), orientation_tail3.home.raw(~highs.orientation_tail3.home), t(highs.orientation_tail3.home));
orientation_tail3.home.raw(lows.orientation_tail3.home) = interp1(t(~lows.orientation_tail3.home), orientation_tail3.home.raw(~lows.orientation_tail3.home), t(lows.orientation_tail3.home));
upper.orientation_tail3.intruder = quantile(orientation_tail3.intruder.raw,.97);
lower.orientation_tail3.intruder = quantile(orientation_tail3.intruder.raw,.03);
highs.orientation_tail3.intruder = (orientation_tail3.intruder.raw>upper.orientation_tail3.intruder);
lows.orientation_tail3.intruder = (orientation_tail3.intruder.raw<lower.orientation_tail3.intruder);
orientation_tail3.intruder.raw(highs.orientation_tail3.intruder) = interp1(t(~highs.orientation_tail3.intruder), orientation_tail3.intruder.raw(~highs.orientation_tail3.intruder), t(highs.orientation_tail3.intruder));
orientation_tail3.intruder.raw(lows.orientation_tail3.intruder) = interp1(t(~lows.orientation_tail3.intruder), orientation_tail3.intruder.raw(~lows.orientation_tail3.intruder), t(lows.orientation_tail3.intruder));

%filter
orientation_tail3.home.filt = medfilt1(orientation_tail3.home.raw,10);
orientation_tail3.home.norm = filter(w,1,orientation_tail3.home.filt);
orientation_tail3.home.norm(1:gaussCut(end)+1) = [];
orientation_tail3.intruder.filt = medfilt1(orientation_tail3.intruder.raw,10);
orientation_tail3.intruder.norm = filter(w,1,orientation_tail3.intruder.filt);
orientation_tail3.intruder.norm(1:gaussCut(end)+1) = [];

%plot
plot(t_sec((gaussCut(end)+2):end),orientation_tail3.home.norm-58,'linewidth',2,'color',cc(192,:))
plot(t_sec((gaussCut(end)+2):end),orientation_tail3.intruder.norm-66,'linewidth',2,'color',cc(208,:))

disp('tail_node_three orientation complete')

% orientation to other nose (trunk-nose-other_nose angle)
% 1 is facing towards, -1 backwards, 0 is right angle
orientation_other_nose.home.raw = -1.*ones(noFrames_intruder,1); %orientation_mouse of home mouse to intruder mouse
orientation_other_nose.intruder.raw = -1.*ones(noFrames_intruder,1); %orientation_mouse of intruder mouse to mouse mouse

%get orientation
for ii = 1:noFrames_intruder
    orientation_other_nose.home.raw(ii) =dot(trunk_pos.home.interp(ii,:)-nose_pos.home.interp(ii,:),nose_pos.home.interp(ii,:)-nose_pos.intruder.interp(ii,:))...
        /(norm(trunk_pos.home.interp(ii,:)-nose_pos.home.interp(ii,:))*norm(nose_pos.home.interp(ii,:)-nose_pos.intruder.interp(ii,:)));
    orientation_other_nose.intruder.raw(ii) =dot(trunk_pos.intruder.interp(ii,:)-nose_pos.intruder.interp(ii,:),nose_pos.intruder.interp(ii,:)-nose_pos.home.interp(ii,:))...
        /(norm(trunk_pos.intruder.interp(ii,:)-nose_pos.intruder.interp(ii,:))*norm(nose_pos.intruder.interp(ii,:)-nose_pos.home.interp(ii,:)));
end

%interpolate extreme values
upper.orientation_other_nose.home = quantile(orientation_other_nose.home.raw,.97);
lower.orientation_other_nose.home = quantile(orientation_other_nose.home.raw,.03);
highs.orientation_other_nose.home = (orientation_other_nose.home.raw>upper.orientation_other_nose.home);
lows.orientation_other_nose.home = (orientation_other_nose.home.raw<lower.orientation_other_nose.home);
orientation_other_nose.home.raw(highs.orientation_other_nose.home) = interp1(t(~highs.orientation_other_nose.home), orientation_other_nose.home.raw(~highs.orientation_other_nose.home), t(highs.orientation_other_nose.home));
orientation_other_nose.home.raw(lows.orientation_other_nose.home) = interp1(t(~lows.orientation_other_nose.home), orientation_other_nose.home.raw(~lows.orientation_other_nose.home), t(lows.orientation_other_nose.home));
upper.orientation_other_nose.intruder = quantile(orientation_other_nose.intruder.raw,.97);
lower.orientation_other_nose.intruder = quantile(orientation_other_nose.intruder.raw,.03);
highs.orientation_other_nose.intruder = (orientation_other_nose.intruder.raw>upper.orientation_other_nose.intruder);
lows.orientation_other_nose.intruder = (orientation_other_nose.intruder.raw<lower.orientation_other_nose.intruder);
orientation_other_nose.intruder.raw(highs.orientation_other_nose.intruder) = interp1(t(~highs.orientation_other_nose.intruder), orientation_other_nose.intruder.raw(~highs.orientation_other_nose.intruder), t(highs.orientation_other_nose.intruder));
orientation_other_nose.intruder.raw(lows.orientation_other_nose.intruder) = interp1(t(~lows.orientation_other_nose.intruder), orientation_other_nose.intruder.raw(~lows.orientation_other_nose.intruder), t(lows.orientation_other_nose.intruder));

%filter
orientation_other_nose.home.filt = medfilt1(orientation_other_nose.home.raw,10);
orientation_other_nose.home.norm = filter(w,1,orientation_other_nose.home.filt);
orientation_other_nose.home.norm(1:gaussCut(end)+1) = [];
orientation_other_nose.intruder.filt = medfilt1(orientation_other_nose.intruder.raw,10);
orientation_other_nose.intruder.norm = filter(w,1,orientation_other_nose.intruder.filt);
orientation_other_nose.intruder.norm(1:gaussCut(end)+1) = [];

%plot
plot(t_sec((gaussCut(end)+2):end),orientation_other_nose.home.norm-68,'linewidth',2,'color',cc(224,:))
plot(t_sec((gaussCut(end)+2):end),orientation_other_nose.intruder.norm-70,'linewidth',2,'color',cc(240,:))

disp('angle to other nose complete')

% orientation to other tush (trunk-nose-other_tti angle)
% 1 is facing towards, -1 backwards, 0 is right angle
orientation_other_tush.home.raw = -1.*ones(noFrames_intruder,1); %orientation_mouse of home mouse to intruder mouse
orientation_other_tush.intruder.raw = -1.*ones(noFrames_intruder,1); %orientation_mouse of intruder mouse to mouse mouse

%get orientation
for ii = 1:noFrames_intruder
    orientation_other_tush.home.raw(ii) =dot(trunk_pos.home.interp(ii,:)-nose_pos.home.interp(ii,:),nose_pos.home.interp(ii,:)-tti_pos.intruder.interp(ii,:))...
        /(norm(trunk_pos.home.interp(ii,:)-nose_pos.home.interp(ii,:))*norm(nose_pos.home.interp(ii,:)-tti_pos.intruder.interp(ii,:)));
    orientation_other_tush.intruder.raw(ii) =dot(trunk_pos.intruder.interp(ii,:)-nose_pos.intruder.interp(ii,:),nose_pos.intruder.interp(ii,:)-tti_pos.home.interp(ii,:))...
        /(norm(trunk_pos.intruder.interp(ii,:)-nose_pos.intruder.interp(ii,:))*norm(nose_pos.intruder.interp(ii,:)-tti_pos.home.interp(ii,:)));
end

%interpolate extreme values
upper.orientation_other_tush.home = quantile(orientation_other_tush.home.raw,.97);
lower.orientation_other_tush.home = quantile(orientation_other_tush.home.raw,.03);
highs.orientation_other_tush.home = (orientation_other_tush.home.raw>upper.orientation_other_tush.home);
lows.orientation_other_tush.home = (orientation_other_tush.home.raw<lower.orientation_other_tush.home);
orientation_other_tush.home.raw(highs.orientation_other_tush.home) = interp1(t(~highs.orientation_other_tush.home), orientation_other_tush.home.raw(~highs.orientation_other_tush.home), t(highs.orientation_other_tush.home));
orientation_other_tush.home.raw(lows.orientation_other_tush.home) = interp1(t(~lows.orientation_other_tush.home), orientation_other_tush.home.raw(~lows.orientation_other_tush.home), t(lows.orientation_other_tush.home));
upper.orientation_other_tush.intruder = quantile(orientation_other_tush.intruder.raw,.97);
lower.orientation_other_tush.intruder = quantile(orientation_other_tush.intruder.raw,.03);
highs.orientation_other_tush.intruder = (orientation_other_tush.intruder.raw>upper.orientation_other_tush.intruder);
lows.orientation_other_tush.intruder = (orientation_other_tush.intruder.raw<lower.orientation_other_tush.intruder);
orientation_other_tush.intruder.raw(highs.orientation_other_tush.intruder) = interp1(t(~highs.orientation_other_tush.intruder), orientation_other_tush.intruder.raw(~highs.orientation_other_tush.intruder), t(highs.orientation_other_tush.intruder));
orientation_other_tush.intruder.raw(lows.orientation_other_tush.intruder) = interp1(t(~lows.orientation_other_tush.intruder), orientation_other_tush.intruder.raw(~lows.orientation_other_tush.intruder), t(lows.orientation_other_tush.intruder));

%filter
orientation_other_tush.home.filt = medfilt1(orientation_other_tush.home.raw,10);
orientation_other_tush.home.norm = filter(w,1,orientation_other_tush.home.filt);
orientation_other_tush.home.norm(1:gaussCut(end)+1) = [];
orientation_other_tush.intruder.filt = medfilt1(orientation_other_tush.intruder.raw,10);
orientation_other_tush.intruder.norm = filter(w,1,orientation_other_tush.intruder.filt);
orientation_other_tush.intruder.norm(1:gaussCut(end)+1) = [];

%plot
plot(t_sec((gaussCut(end)+2):end),orientation_other_tush.home.norm-72,'linewidth',2,'color',cc(224,:))
plot(t_sec((gaussCut(end)+2):end),orientation_other_tush.intruder.norm-74,'linewidth',2,'color',cc(240,:))

disp('angle to other tush complete')

%tailtip velocity
velo_tailtip.home.raw = -1.*ones(noFrames_intruder-1,1);
velo_tailtip.intruder.raw = -1.*ones(noFrames_intruder-1,1);
for ii = 1:noFrames_intruder-1
    velo_tailtip.home.raw(ii) = norm(tailtip_pos.home.interp(ii+1,:)-tailtip_pos.home.interp(ii,:)); %velocity in pixels/frame
    velo_tailtip.intruder.raw(ii) = norm(tailtip_pos.intruder.interp(ii+1,:)-tailtip_pos.intruder.interp(ii,:)); %velocity in pixels/frame
end

%interpolate extreme values
upper.velo_tailtip.home = quantile(velo_tailtip.home.raw,.97);
lower.velo_tailtip.home = quantile(velo_tailtip.home.raw,.03);
highs.velo_tailtip.home = (velo_tailtip.home.raw>upper.velo_tailtip.home);
lows.velo_tailtip.home = (velo_tailtip.home.raw<lower.velo_tailtip.home);
velo_tailtip.home.raw(highs.velo_tailtip.home) = interp1(t(~highs.velo_tailtip.home), velo_tailtip.home.raw(~highs.velo_tailtip.home), t(highs.velo_tailtip.home));
velo_tailtip.home.raw(lows.velo_tailtip.home) = interp1(t(~lows.velo_tailtip.home), velo_tailtip.home.raw(~lows.velo_tailtip.home), t(lows.velo_tailtip.home));
upper.velo_tailtip.intruder = quantile(velo_tailtip.intruder.raw,.97);
lower.velo_tailtip.intruder = quantile(velo_tailtip.intruder.raw,.03);
highs.velo_tailtip.intruder = (velo_tailtip.intruder.raw>upper.velo_tailtip.intruder);
lows.velo_tailtip.intruder = (velo_tailtip.intruder.raw<lower.velo_tailtip.intruder);
velo_tailtip.intruder.raw(highs.velo_tailtip.intruder) = interp1(t(~highs.velo_tailtip.intruder), velo_tailtip.intruder.raw(~highs.velo_tailtip.intruder), t(highs.velo_tailtip.intruder));
velo_tailtip.intruder.raw(lows.velo_tailtip.intruder) = interp1(t(~lows.velo_tailtip.intruder), velo_tailtip.intruder.raw(~lows.velo_tailtip.intruder), t(lows.velo_tailtip.intruder));

%filter
velo_tailtip.home.filt = medfilt1(velo_tailtip.home.raw,10);
velo_tailtip.home.gauss = filter(w,1,velo_tailtip.home.filt);
velo_tailtip.home.gauss(1:gaussCut(end)) = [];
velo_tailtip.home.norm = rescale(velo_tailtip.home.gauss,-1,1);
velo_tailtip.intruder.filt = medfilt1(velo_tailtip.intruder.raw,10);
velo_tailtip.intruder.gauss = filter(w,1,velo_tailtip.intruder.filt);
velo_tailtip.intruder.gauss(1:gaussCut(end)) = [];
velo_tailtip.intruder.norm = rescale(velo_tailtip.intruder.gauss,-1,1);

%plot
plot(t_sec((gaussCut(end)+2):end),velo_tailtip.home.norm-76,'linewidth',2,'color',cc(224,:))
plot(t_sec((gaussCut(end)+2):end),velo_tailtip.intruder.norm-78,'linewidth',2,'color',cc(240,:))

disp('tailtip velocity complete')

disp('feature extraction complete')
disp('///////////////////////')

ylim([-79 1])

%% SAVE OUTPUT
appendThis = input('file name for save (in quotes):');
save([appendThis,'_extracted'])

 %% tSNE
% subplot(3,4,4)
% posture_tSNE.full = tsne([dist_im.norm, dist_nn.norm, dist_nt.norm, velo_trunk.home.norm, velo_trunk.intruder.norm, velo_im.norm, velo_nose.home.norm, velo_nose.intruder.norm, ...
%     dist_cc.home.norm, dist_cc.intruder.norm, velo_forelimb_left.home.norm, velo_forelimb_right.home.norm, velo_forelimb_left.intruder.norm, velo_forelimb_right.intruder.norm, ...
%     dist_mouse.home.norm, dist_mouse.intruder.norm, orientation_mouse.home.norm, orientation_mouse.intruder.norm, orientation_forelimb_left.home.norm, orientation_forelimb_left.intruder.norm ...
%     orientation_body.home.norm, orientation_body.intruder.norm, orientation_forelimb_right.home.norm, orientation_forelimb_right.intruder.norm, dist_nneck.norm, dist_ntrunk.norm, ...
%     orientation_tailbase.home.norm, orientation_tailbase.intruder.norm, orientation_tail1.home.norm, orientation_tail1.intruder.norm, orientation_tail2.home.norm, orientation_tail2.intruder.norm, ...
%     orientation_tail3.home.norm, orientation_tail3.intruder.norm]);
% gscatter(posture_tSNE.full(:,1),posture_tSNE.full(:,2));
% lostFrames.full = noFrames_intruder-length(posture_tSNE.full) 
% lostFrames.full = 100*(lostFrames.full/noFrames_intruder) %#ok<*NOPTS>
% 
% subplot(3,4,8)
% posture_tSNE.lw = tsne([dist_im.norm, dist_nn.norm, dist_nt.norm, velo_trunk.home.norm, velo_trunk.intruder.norm, velo_im.norm, ...
%     dist_cc.home.norm, dist_cc.intruder.norm, orientation_mouse.home.norm, orientation_mouse.intruder.norm]);
% gscatter(posture_tSNE.lw(:,1),posture_tSNE.lw(:,2));
% lostFrames.lw = noFrames_intruder-length(posture_tSNE.lw)
% lostFrames.lw = 100*(lostFrames.lw/noFrames_intruder)
% 
% subplot(3,4,12)
% posture_tSNE.lw_plus = tsne([dist_im.norm, dist_nn.norm, dist_nt.norm, velo_trunk.home.norm, velo_trunk.intruder.norm, velo_im.norm, ...
%     dist_mouse.home.norm, dist_mouse.intruder.norm, dist_cc.home.norm, dist_cc.intruder.norm, orientation_mouse.home.norm, orientation_mouse.intruder.norm, ...
%     orientation_body.home.norm, orientation_body.intruder.norm, dist_nneck.norm, dist_ntrunk.norm, orientation_tailbase.home.norm, orientation_tailbase.intruder.norm]);
% gscatter(posture_tSNE.lw_plus(:,1),posture_tSNE.lw_plus(:,2));
% 
% lostFrames.lw_plus = noFrames_intruder-length(posture_tSNE.lw_plus)
% lostFrames.lw_plus = 100*(lostFrames.lw_plus/noFrames_intruder)
% disp('tSNE complete')

%% sub functions
function frameNumber = findFrame(ri_frame,intruderIn) %#ok<DEFNU>
frameNumber = ri_frame + intruderIn - 1;

end