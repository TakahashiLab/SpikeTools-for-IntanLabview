function [tuning_curve, theta, wu2,ang_hd, mr,mvl,occupancy,num_spikes,sig,spk,headdir] = DirectionalTuningFnc(Spks, Traj, msT, varargin)
% original from https://github.com/hasselmonians/CMBHOME
%modified 2024/05/16

shuffleLag = 20000; %20sec
SpksShuffle = [];
sig=0;
spk=[];
headdir=[];

p = inputParser;
p.addParamValue('fstart', 0, @isnumeric);
p.addParamValue('animal', 'rat', @ischar);
p.addParamValue('maxdist', 400, @isnumeric);
p.addParamValue('khz', 31.25, @isnumeric);
p.addParamValue('spon', 0, @isnumeric);
p.addParamValue('speed', 2.5, @isnumeric);
p.addParamValue('shuffle', 0, @isnumeric);
p.addParamValue('shufflen', 100, @isnumeric);
p.addParamValue('shuffletype', 'shift', @ischar);
p.addParamValue('verbose', 0, @isnumeric);
p.addParamValue('posture', [3 4 5 6], @isvector);
p.addParamValue('dirtype', 'head', @ischar);

p.parse(varargin{:});
fstart = p.Results.fstart;
animal = p.Results.animal;
kHz = p.Results.khz;
maxdist = p.Results.maxdist;
spON = p.Results.spon;
ThS = p.Results.speed;
shuffle = p.Results.shuffle;
shuffleN = p.Results.shufflen;
shuffleType = p.Results.shuffletype;
verbose = p.Results.verbose;
posture = p.Results.posture;
dirType = p.Results.dirtype;

binsize = .5; %.5
theta = -180 + binsize / 2:binsize:180 - binsize / 2;
theta = theta(:);

BinWidthCm = 2.5;

%ThS=2;
yLen = max(Traj(:, 2)) - min(Traj(:, 2));
xLen = max(Traj(:, 1)) - min(Traj(:, 1));
Len = yLen;

if yLen < xLen
    Len = xLen;
end

if nargin == 4
    cmPerPixel = 120 / Len; %120cm circle;
    Bin = BinWidthCm / cmPerPixel;
    spON = 0;
    shuffle = 0;
end

switch lower(animal)
    case 'bird'
        cmPerPixel = maxdist / Len; %120cm circle;
        Bin = BinWidthCm / cmPerPixel;
        spON = 0;
    case 'rat'
        cmPerPixel = maxdist / Len; %
        Bin = BinWidthCm / cmPerPixel;
        kHz = 25;
        spON = 1;
        ThS = 2.5;
        msT = msT / kHz;
    case 'fish'
        cmPerPixel = maxdist / Len; %120cm circle;
        Bin = BinWidthCm / cmPerPixel;
        %spON=0;
    case 'seaturtle',
        cmPerPixel=20/Len;%120cm circle;
        Bin=BinWidthCm/cmPerPixel;
        kHz=30;
        spON=0;
        ThS=2.5;
        binsize = 6;%0.5
        theta=-180+binsize/2:binsize:180-binsize/2;
        theta = theta(:);
        posture=[1:4];
end

FPS = floor(1 / (median(diff(msT)) / 1000));
fs_video = FPS;
msFPS = floor(1 / FPS * 1000);

movement = sqrt(sum(diff(Traj(:, 1:2)) .^ 2, 2)) * cmPerPixel; %
speed = movement ./ diff(msT);
speed = smooth(speed, FPS, 'moving') * 1000;

if spON
    Good = find(speed >= ThS);
    Traj = Traj(Good, :);
    msT = msT(Good);
end

Spks = ceil(Spks / kHz) + fstart; %msec

Traj = Traj';
x = Traj(posture(1), :);
y = Traj(posture(2), :);
x2 = Traj(posture(3), :);
y2 = Traj(posture(4), :);
Traj = Traj';

sx = smooth(x, 30)';
sy = smooth(y, 30)';
sx2 = smooth(x2, 30)';
sy2 = smooth(y2, 30)';

headdir = atan2d(y - y2, x - x2)';

%moving direction
pX = [sx(30:end) sx(1:29)];
pY = [sy(30:end) sy(1:29)];

diffX = pX - sx;
diffY = pY - sy;
hovering = find(speed < 5)'; %5cm/s

if ~isempty(hovering)
    if hovering(1) == 1
        hovering(1) = [];
    end
end

movedir = atan2d(diffY, diffX)';

%hovering
movedir(hovering) = headdir(hovering);

switch (dirType)
    case 'move',
        headdir = movedir;

    case 'diff',
        headdir = movedir - headdir;
        headdir(headdir > 180) = headdir(headdir > 180) - 360;
        headdir(headdir < -180) = headdir(headdir <- 180) + 360;
end

StartTraj = msT(1);
EndTraj = msT(end);

MRs = [];
WU2s = [];
MVLs = [];

if shuffle

    beginSpks = msT(1);
    endSpks = msT(end);
    entireLength = ceil(endSpks - beginSpks);
    rp = randperm(entireLength - shuffleLag * 2 + 1);
    rp = rp(1:shuffleN);
    orgRp = shuffleLag:(entireLength - shuffleLag);
    rp = orgRp(rp);
    lenOrgRp = length(orgRp);
    lenSpks = length(Spks);
    SpksShuffle = zeros(shuffleN, lenSpks);

    switch (shuffleType)
        case 'shift',

            for i = 1:shuffleN
                spkss = Spks + rp(i);
                topSpks = find(spkss >= endSpks);
                spkss(topSpks) = spkss(topSpks) - endSpks + beginSpks;
                SpksShuffle(i, :) = sort(spkss);
            end

        case 'bootstrap',

            for i = 1:shuffleN
                bootsamp = sort(randsample(length(Spks), length(Spks), true));
                SpksShuffle(i, :) = Spks(bootsamp);
            end

    end

    parfor i = 1:size(SpksShuffle, 1)

        Spks = SpksShuffle(i, :);
        Spks = Spks(find(Spks > StartTraj & Spks < EndTraj));
        spk_headdir = [];

        for j = 1:(size(headdir, 1) - 1)
            if msT(j+1)-msT(j) < msFPS*2
                spk_headdir = [spk_headdir ones(1, sum(Spks >= msT(j) & Spks < msT(j + 1))) .* headdir(j)];
            end
        end

        [~, ~, mr, wu2, mvl] = DirectionTuningCore(headdir, spk_headdir, binsize, fs_video, theta);

        MRs = [MRs mr];
        WU2s = [WU2s wu2];
        MVLs = [MVLs mvl];

    end

    tuning_curve = [];
    ang_hd = [];
    mr = MRs;
    wu2 = WU2s;
    mvl = MVLs;
    occupancy = [];
    num_spikes = [];

else
    Spks = Spks(find(Spks > StartTraj & Spks < EndTraj));

    spk_headdir = [];
    spk_headdir_t = [];

    for i = 1:(size(headdir, 1) - 1)
        if msT(i+1)-msT(i) < msFPS*2
            spk_headdir=[spk_headdir ones(1,sum(Spks >= msT(i) & Spks < msT(i+1))).*headdir(i)];
            spk_headdir_t=[spk_headdir_t ones(1,sum(Spks >= msT(i) & Spks < msT(i+1))).*msT(i)];


        end
    end

    spk=[spk_headdir; spk_headdir_t];

    switch lower(animal)
        case 'seaturtle',
            [tuning_curve,ang_hd,mr,wu2,mvl,occupancy,num_spikes,sig]=DirectionTuningCore2(headdir,spk_headdir,binsize,fs_video,theta,spk_headdir_t);
        otherwise,
            [tuning_curve, ang_hd, mr, wu2, mvl, occupancy, num_spikes] = DirectionTuningCore(headdir, spk_headdir, binsize, fs_video, theta);
    end
end

return;
%%%%%%%%%%%%%%%%%%%%%%%%
function [tuning_curve, ang_hd, mr, wu2, mvl, angle_occupancy, num_spikes] = DirectionTuningCore(headdir, spk_headdir, binsize, fs_video, theta)

spk_headdir = spk_headdir';
%%%%%%%%%%%%%%%%%%%%%%%%head direction
angle_occupancy = DirectionalOccupancy(binsize, headdir, fs_video);
num_spikes = hist(spk_headdir, theta);

window = 29;

num_spikes = [fliplr(num_spikes(end - window + 1:end)) num_spikes fliplr(num_spikes(1:window))];

angle_occupancy = [fliplr(angle_occupancy(end - window + 1:end)) angle_occupancy fliplr(angle_occupancy(1:window))];
num_spikes = smooth(num_spikes, window, 'moving')';
angle_occupancy = smooth(angle_occupancy, window, 'moving')';

tuning_curve = num_spikes ./ angle_occupancy;

tuning_curve = tuning_curve(window + 1:end - window);
angle_occupancy = angle_occupancy(window + 1:end - window);
num_spikes = num_spikes(window + 1:end - window);

[ang_hd, mr, mvl] = GetOtherStats(tuning_curve, theta);
wu2 = WatsonsU2(spk_headdir, headdir);

return;
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
function [tuning_curve,ang_hd,mr,wu2,mvl,angle_occupancy,num_spikes,sig]=DirectionTuningCore2(headdir,spk_headdir,binsize,fs_video,theta,spk_headdir_t)
sig=0;ang_hd=NaN;
gaussianSD = 3; % Standard deviation for Gaussian smoothing
significantPValue = 0.0001; % Threshold for Rayleigh test
minPeakRate = 1; % Minimum peak firing rate in Hz
minZScore = 50; % Minimum Z-score

tuning_curve=NaN*ones(1,60);
spk_headdir=spk_headdir';
%%%%%%%%%%%%%%%%%%%%%%%%head direction
angle_occupancy = DirectionalOccupancy(binsize, headdir,fs_video);

if isempty(spk_headdir)

    ang_hd=NaN;
    mr=0;
    wu2=0;
    mvl=0;
    num_spikes=0;

else
    num_spikes = hist(spk_headdir,theta);
    
    window=floor(size(num_spikes,2)/2)-1;
    num_spikes=[fliplr(num_spikes(end-window+1:end)) num_spikes fliplr(num_spikes(1:window))];

    angle_occupancy=[fliplr(angle_occupancy(end-window+1:end)) angle_occupancy fliplr(angle_occupancy(1:window))];
    %num_spikes = smooth(num_spikes,window,'moving')';
    %angle_occupancy = smooth(angle_occupancy,window,'moving')';

    tuning_curve = num_spikes ./ angle_occupancy;

    % plot(num_spikes,'r');
    % hold on;
    % plot(tuning_curve,'g');
    % plot(angle_occupancy,'k');

    tuning_curve(isnan(tuning_curve))=0;

    % Apply Gaussian smoothing
    kernelSize = ceil(gaussianSD * 6); % 3 standard deviations on either side
    gaussianKernel = fspecial('gaussian', [1, kernelSize], gaussianSD);
    tuning_curve = conv2(tuning_curve, gaussianKernel, 'same');


    % Rayleigh test for uniformity
    [pValue, z] = circ_rtest(deg2rad(spk_headdir), spk_headdir_t);



    %tuning_curve=smooth(tuning_curve,window,'moving')';

    tuning_curve = tuning_curve(window+1:end-window);
    angle_occupancy = angle_occupancy(window+1:end-window);
    num_spikes = num_spikes(window+1:end-window);

    % Find peak firing rate
    peakRate = max(tuning_curve);

    % Check conditions
    if pValue < significantPValue && z > minZScore && peakRate > minPeakRate
        sig=1;
    end

    [ang_hd, mr,mvl] = GetOtherStats(tuning_curve, theta);
    wu2 = WatsonsU2(spk_headdir, headdir);
end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ang_hd, mr, mvl] = GetOtherStats(tuning_curve, theta)

tf = ~isnan(tuning_curve); % where we arent nan

theta = theta * unitsratio('rad', 'deg');

theta = theta(tf); % remove nans

tuning_curve = tuning_curve(tf)';

xs = tuning_curve .* cos(theta); % average
ys = tuning_curve .* sin(theta);

ang_hd = atan2(mean(ys), mean(xs)); % mean direction

mr = (cos(ang_hd) * sum(xs) + sin(ang_hd) * sum(ys)) / sum(tuning_curve); % mean resultant length

mvl = sqrt(mean(ys) ^ 2 + mean(xs) ^ 2); %mean vector length

return;
