% 
%load data on the assumption 
% that 3 aux inputs exist
%
function [out,ref]=loadauxilvrcN(filename)

byte=2;
step=1000000;
numOfElec=1;
numOfEvent=3;%default

fp=fopen(filename,'r','b');

%header
numOfCh=fread(fp,1,'uchar');
Hz=fread(fp,1,'float32');
ref=fread(fp,1,'float32');

headerPos=ftell(fp);
fseek(fp,0,1);
pos=ftell(fp);

size=pos-headerPos;

loop=floor(size/(numOfCh*byte*step));
lastStep=(size-loop*numOfCh*byte*step)/(numOfCh*byte);

numOfMLE=numOfEvent;

out=[];

Skip=(numOfCh-numOfMLE)*byte;
blocksize=[numOfMLE step];
lastblocksize=[numOfMLE lastStep];
fseek(fp,headerPos,-1);

Precision=sprintf('%d*int16=>int16',numOfMLE);

for i=1:loop
    tmp=fread(fp,blocksize,Precision,Skip);
    out=[out tmp];
end

tmp=fread(fp,lastblocksize,Precision,Skip);
out=[out tmp];

fclose(fp);

return;
