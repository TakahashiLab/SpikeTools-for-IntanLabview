function phase=WhitePhase(sig,fs)
N=length(sig);
X=fft(sig,N);
f=(0:N-1)*(fs/N);

% 周波数帯の分離
df = 2; % 周波数刻み
fmin = 2; % 最小周波数
fmax = 120; % 最大周波数
nbands = (fmax-fmin)/df + 1; % 周波数帯の数

phase=zeros(nbands,N);

% 各周波数帯ごとに処理
for i=1:nbands
    
    % 周波数帯の範囲を決める
    fstart = fmin + (i-1)*df;
    fend = fstart + df;
    
    % 周波数帯に対応するインデックスを求める
    istart = find(f>=fstart,1);
    iend = find(f>=fend,1)-1;
    
    % 周波数帯以外の成分を0にする
    Y = X;
   
    Y(1:istart-1) = 0;
    Y(iend+1:end) = 0;
    X2=ifft(Y);
        
    % phaseを計算する（ラジアン）
    phase(i,:) = angle(X2);
end