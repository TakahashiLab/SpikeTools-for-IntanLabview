function [Out,msT,st,Pos,States,vt,StatesAll] = readSpks(directoryPath,varargin)
debug=0;
StatesAll=[];
p = inputParser;
p.addParamValue('animal', 'rat', @ischar);
p.addParamValue('khz', 30, @isnumeric);
p.parse(varargin{:});
animal = p.Results.animal;
kHz = p.Results.khz;

%kilosort timing signal
% 指定されたディレクトリパターンに従ってサブディレクトリを検索
ksDir = dir(fullfile(directoryPath, '*.kilosort'));
ksDir = dir(fullfile(ksDir(1).folder, ksDir(1).name));
if isempty(ksDir)
    error('kilosort.timestampファイルが見つかりません。');
end

% 各.kilosortファイルについて*timestamps.datにマッチするか調べる
for i = 1:length(ksDir)
    currentFileName = fullfile(ksDir(i).folder, ksDir(i).name);
    % ファイル名が*Din1.datにマッチするか確認
    if contains(currentFileName, 'timestamps.dat')
        ksTime = currentFileName;
    end
end
%state timing
tmp=readTrodesExtractedDataFile(ksTime);
[t]=tmp.fields.data;
ksT=t';
ksT=double(ksT);


% 指定されたディレクトリパターンに従ってサブディレクトリを検索
kilosortDir = fullfile(directoryPath, '*.kilosort');
stFolderDir = fullfile(kilosortDir, 'stFolder');
kilosort4Dir = fullfile(stFolderDir, 'kilosort4');
mPhyDir = dir(fullfile(kilosort4Dir, 'mPhy*'));

% mPhyで始まるディレクトリが存在するか確認
if isempty(mPhyDir)
    error('指定されたパターンに一致するディレクトリが見つかりません。');
end

% 最初に見つかったmPhyディレクトリのフルパスを取得
targetDir = fullfile(mPhyDir(1).folder, mPhyDir(1).name);
fprintf('loading spikes from %s\n',targetDir);

% ファイルパスを生成
spikeTimesFile = fullfile(targetDir, 'spike_times.npy');
spikeClustersFile = fullfile(targetDir, 'spike_clusters.npy');
clusterGroupsFile = fullfile(targetDir, 'cluster_group.tsv');

% データの読み込み
t = readNPY(spikeTimesFile);
c = readNPY(spikeClustersFile);
[cid, grp] = readClusterGroupsCSV(clusterGroupsFile);

% クラスターのグループが '2' (good) であるクラスターIDを抽出
goodClusterIDs = cid(grp == 2);

% 出力用セル配列を初期化
Out = cell(length(goodClusterIDs), 1);

% 各 'good' クラスターについて処理
for i = 1:length(goodClusterIDs)
    % 現在のクラスターIDに対応する発火時間を取得
    currentClusterTimes = t(c == goodClusterIDs(i))';

    % セル配列に発火時間を格納
    Out{i, 3} = ksT(currentClusterTimes);
end


% 指定されたディレクトリパターンに従ってサブディレクトリを検索
dioDir = dir(fullfile(directoryPath, '*.DIO'));
dioDir = dir(fullfile(dioDir(1).folder, dioDir(1).name));
if isempty(dioDir)
    error('DIOファイルが見つかりません。');
end

% 各.DIOファイルについて*Din1.datにマッチするか調べる
for i = 1:length(dioDir)
    currentFileName = fullfile(dioDir(i).folder, dioDir(i).name);
    % ファイル名が*Din1.datにマッチするか確認
    if contains(currentFileName, 'Din1.dat')
        VideoTime = currentFileName;
        % マッチしたファイルが見つかったら終了
    elseif contains(currentFileName, 'Din2.dat')
        StateTime = currentFileName;
    end
end

%video signal
tmp=readTrodesExtractedDataFile(VideoTime);
[t,state]=tmp.fields.data;
vt=t(state==1);

%state timing
tmp=readTrodesExtractedDataFile(StateTime);
[t,state]=tmp.fields.data;
st=t';

if strmatch(animal,'bird')
    id=find((diff(st)/30000/60)<1);
    st(id+1)=[];
    id=find((diff(st)/30000/60)>8);
end

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
    if size(st,2)>16
        slen=16;
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
            id=find(vt>st(5));
            msTC=[msTC; msT(id(1):id(1)+PosLen(2)-1)];
            id=find(vt>st(9));
            msTC=[msTC; msT(id(1):id(1)+PosLen(3)-1)];
            id=find(vt>st(13));
            msTC=[msTC; msT(id(1):id(1)+PosLen(4)-1)];
            vt=msTC;

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


