function [Out,msT,st,Pos,States,vt,StatesAll] = readSpksSGLX(baseDir,inputPrefix,varargin)
debug=0;
StatesAll=[];
p = inputParser;
p.addParamValue('animal', 'rat', @ischar);
p.addParamValue('khz', 30, @isnumeric);
p.parse(varargin{:});
animal = p.Results.animal;
kHz = p.Results.khz;
Hz=1000*kHz;

if startsWith(lower(animal),'sea')
    animal=animal(4:end);
end

% ユーザーが指定する部分
%inputPrefix = '2024092501C'; % ユーザーが入力する変更可能な部分

% ベースディレクトリの設定
%baseDir = 'P:\seabird\2024\ORIGINAL\SUBJECTR\TEST';


% 'catgt_' + inputPrefix で始まるディレクトリを取得
dirPattern = fullfile(baseDir, ['catgt_' inputPrefix '_g0'])
dirs = dir(dirPattern);

% フォルダが見つかった場合に処理を実行
if isempty(dirs)
    error('No directory matching pattern "catgt_%s_g0" found.', inputPrefix);
end

catgtDir = fullfile(baseDir, dirPattern)


imec0Dir = fullfile(catgtDir, [inputPrefix '_g0_imec0']);

% 'imec????_ks4' フォルダを対象にする（ディレクトリのみ）
imecFolders = dir(fullfile(imec0Dir, 'imec*'));
imecFolders = imecFolders([imecFolders.isdir]); % ディレクトリのみを抽出
Out=[];

for i = 1:length(imecFolders)
    folderName = imecFolders(i).name;

    % imec????_ks4_org フォルダは除外
    if contains(folderName, '_ks4_orig')
        continue;
    end

    ks4Dir = fullfile(imec0Dir, folderName);
    spikeTimesFile = fullfile(ks4Dir, 'spike_times_sec_adj.npy')
    spikeClustersFile = fullfile(ks4Dir, 'spike_clusters.npy');
    clusterGroupsFile = fullfile(ks4Dir, 'cluster_group.tsv');

    % データの読み込み
    t = readNPY(spikeTimesFile);
    c = readNPY(spikeClustersFile);
    [cid, grp] = readClusterGroupsCSV(clusterGroupsFile);

    % クラスターのグループが '2' (good) であるクラスターIDを抽出
    goodClusterIDs = cid(grp == 2);

    % 出力用セル配列を初期化
    bufOut = cell(length(goodClusterIDs), 1);

    % 各 'good' クラスターについて処理
    for i = 1:length(goodClusterIDs)
        % 現在のクラスターIDに対応する発火時間を取得
        currentClusterTimes = t(c == goodClusterIDs(i))';
        % セル配列に発火時間を格納
        bufOut{i, 3} = currentClusterTimes*Hz;
    end

    Out=[Out;bufOut];
    
end

% 動画信号として xa_1_0.txt を読み込む
videoSigFile = fullfile(catgtDir, '*xa_1_0.txt');
videoFiles = dir(videoSigFile);
if ~isempty(videoFiles)
    videoSigPath = fullfile(catgtDir, videoFiles(1).name);
    videosig = load(videoSigPath);
    %fprintf('Read video signal from: %s\n', videoSigPath);
    vt=videosig*Hz;

end

% イベントアップとして xa_2_0.txt を読み込む
eventUpFile = fullfile(catgtDir, '*xa_2_0.txt');
eventUpFiles = dir(eventUpFile);
if ~isempty(eventUpFiles)
    eventUpPath = fullfile(catgtDir, eventUpFiles(1).name);
    event_up = load(eventUpPath);
    %fprintf('Read event up from: %s\n', eventUpPath);
    eu=event_up*Hz;
end

% イベントダウンとして xia_2_0.txt を読み込む
eventDownFile = fullfile(catgtDir, '*xia_2_0.txt');
eventDownFiles = dir(eventDownFile);
if ~isempty(eventDownFiles)
    eventDownPath = fullfile(catgtDir, eventDownFiles(1).name);
    event_down = load(eventDownPath);
    %fprintf('Read event down from: %s\n', eventDownPath);
    ed=event_down*Hz;
end


return;
%%%%%%%%%%%%%%%%%%%


% motion tracking
% 指定されたディレクトリパターンに従ってサブディレクトリを検索
vDir = dir(fullfile(directoryPath, 'video', '**', '*.mat')); % Include '**' to search all subdirectories for .mat files

Pos=[];
PosLen=[];
axisV=[];
for i = 1:length(vDir)
    matFile = fullfile(vDir(i).folder, vDir(i).name);

    load(matFile,'xydims');

    if isempty(strfind(matFile,'axis'))
        disp(matFile)
        PosLen=[PosLen size(xydims,2)];
        Pos=[Pos; xydims']; % Display the path of each .mat file
    else%get axis
        axisV=[axisV; xydims'];
    end
end

size(vt)
PosLen
sum(PosLen)
if   debug
    plot(Pos(:,1),Pos(:,2),'k.');
    hold on;
    plot(Pos(1,7),Pos(1,8),'ro');
    plot(Pos(1,9),Pos(1,10),'go');
    plot(Pos(1,11),Pos(1,12),'bo');
end

switch(animal)
    case 'turtle',
        if size(Pos,2)==12
            Pos(:,[11 12])=Pos(:,[9 10]);
        end
        Pos=correctNorth(Pos);

    case 'bird',
        if size(Pos,2)==12
            Pos(:,[11 12])=Pos(:,[9 10]);
            Pos=correctNorth(Pos);
        else
            Pos=correctNorth(Pos,axisV);
        end
end

States=[];
if ismember(animal, {'bird', 'turtle'})
    s=1;
    slen=8;
  
    if size(st,2) > 16
        slen=16;
    elseif size(st,2)==4
        slen=3;
    end

    switch(animal)
        case 'turtle',
            %correction
            msT=vt;
            msTC=[];
            msTC=[msTC msT(1:1+PosLen(1)-1)];
            id=find(vt>st(5));
            msTC=[msTC; msT(id(1):id(1)+PosLen(2)-1)];
            if length(PosLen)>2
                id=find(vt>st(9));
                msTC=[msTC; msT(id(1):id(1)+PosLen(3)-1)];
                id=find(vt>st(13));
                id(1)
                msTC=[msTC; msT(id(1):id(1)+PosLen(4)-1)];
            end
            vt=msTC;

            for i=1:4:slen
                %divide states
                id=find(vt > st(i) & vt <st(i+4));
                States{s}=id;
                s=s+1;
            end
            if size(msT,1)>size(Pos,1)
                msT=msT(1:size(Pos,1));
            end
            Pos(:,1:2:end)=-Pos(:,1:2:end);%myPolar coordinate west-east transpose

        case 'bird',
            %correction
           
            msT=vt;
            msTC=[];
            msTC=[msTC msT(1:1+PosLen(1)-1)];

            if length(PosLen)>1
                id=find(vt>st(5));
                msTC=[msTC; msT(id(1):id(1)+PosLen(2)-1)];
                id=find(vt>st(9));
                msTC=[msTC; msT(id(1):id(1)+PosLen(3)-1)];
                id=find(vt>st(13));
                msTC=[msTC; msT(id(1):id(1)+PosLen(4)-1)];
            end
            vt=msTC;
            size(st)
            States=cell(1,2);
            StatesAll=cell(1,slen);
            for i=1:slen
                %divide states
                id=find(vt > st(i) & vt <st(i+1));
                
                StatesAll{i}=id;
                if mod(i,2)
                    States{1}=[States{1}; id];%normal
                else
                    States{2}=[States{2}; id];%counterbalance / indonesia
                end
            end
        

            %debug
            if debug
                plot(Pos(:,1),Pos(:,2),'k.');
                hold on;
                plot(Pos(1,7),Pos(1,8),'ro');
                plot(Pos(1,11),Pos(1,12),'bo');
            end

            %change trackpoint
            Pos=Pos(:,[3:6 1:2 7 8]);
            %Pos(:,2:2:end)=-Pos(:,2:2:end);%myPolar coordinate north-south transpose
            Pos(:,1:2:end)=-Pos(:,1:2:end);%myPolar coordinate west-east transpose

    end
end

msT=double(vt)./kHz;
return;
%%
function Pos = correctNorth(Pos,axisV)
debug =0;

if nargin==1
    % 各行に対して処理を行う
    % 北向きベクトル（その行のデータに基づいて計算）
    north_vector = Pos(1, 11:12) - Pos(1, 7:8);
    if debug
        plot(Pos(1,11),Pos(1,12),'ro');
        hold on;
        plot(Pos(1,7),Pos(1,8),'go');
        plot(Pos(1,9),Pos(1,10),'bo');
    end

elseif nargin==2
    if debug
        plot(axisV(1,1),axisV(1,2),'ro');
        hold on;
        plot(axisV(1,3),axisV(1,4),'go');
        plot(axisV(1,5),axisV(1,6),'bo');
        plot(Pos(:,1),Pos(:,2))
    end
    north_vector = axisV(1, 3:4) - axisV(1, 1:2);
end

% 北向きベクトルとy軸の正の方向との角度を計算
angle_to_north =  -atan2(north_vector(2), north_vector(1))+pi/2
for row = 1:size(Pos, 1)

    % 2D回転行列を作成
    R = [cos(angle_to_north) -sin(angle_to_north); sin(angle_to_north) cos(angle_to_north)];

    % Posの各座標点を回転
    for i = 1:2:size(Pos, 2)
        rotated_point = R * Pos(row, i:i+1)';
        Pos(row, i:i+1) = rotated_point';
    end
end
return;


