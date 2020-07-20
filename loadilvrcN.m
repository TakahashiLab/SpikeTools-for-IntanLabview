% 
%load data on the assumption 
% that 8 tetrodes, and 8 events are fully assigned
% tetrodes: 1-8
% Event:        9
%
function [out,ref,Hz]=loadilvrcN(filename,elecNum,numOfElec,elecType,realConv,atOnce)

byte=2;
step=1000000;

numOfEvent=8;%default
if nargin==2
    numOfElec=8;
    numOfMicrowires=4;
    realConv=0;
    atOnce=1;
else
    switch lower(elecType)
      case 'tetrode',
        numOfMicrowires=4;
      case 'stereotrode',
        numOfMicrowires=2;
    end
end

fp=fopen(filename,'r','b');

%header
numOfCh=fread(fp,1,'uchar');
Hz=fread(fp,1,'float32');
ref=fread(fp,1,'float32');

if atOnce
  numOfMicrowires=numOfCh-numOfEvent;
  numOfElec=1;
  numOfEvent=8;
end


headerPos=ftell(fp);
fseek(fp,0,1);
pos=ftell(fp);

size=pos-headerPos;

loop=floor(size/(numOfCh*byte*step));
lastStep=(size-loop*numOfCh*byte*step)/(numOfCh*byte);


if elecNum<=numOfElec
  initialPos=(elecNum-1)*numOfMicrowires*byte;%initial # of bytes for skip
  numOfMLE=numOfMicrowires;% # of microwires, LFPs or events
elseif elecNum==(numOfElec+1)
  initialPos=numOfElec*numOfMicrowires*byte;
  numOfMLE=numOfEvent;
else
  fprintf('Error elecNum=%d\n',elecNum);
end


out=[];

Skip=(numOfCh-numOfMLE)*byte;
blocksize=[numOfMLE step];
lastblocksize=[numOfMLE lastStep];
fseek(fp,headerPos+initialPos,-1);

Precision=sprintf('%d*int16=>int16',numOfMLE);

for i=1:loop
    tmp=fread(fp,blocksize,Precision,Skip);
    out=[out tmp];
end

tmp=fread(fp,lastblocksize,Precision,Skip);
out=[out tmp];

if realConv
  %conversion intan max amplitude:5mV or 3.V, int16: 2^15 or 2^16
  if elecNum<=numOfElec
    out=double(out)./2^15.*5*1000;%uV
  else
    out=double(out)/2^15.*3.3;%V
  end
end
fclose(fp);

return;
