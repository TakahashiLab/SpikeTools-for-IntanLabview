%%%%%%%%%%%%%%
function tmp=makeMeanWF(wf,step,stepmode)

[num,len]=size(wf);

tmp=reshape(wf,num,step,len/step);
tmp=mean(tmp,3);
%

if nargin==3
  range=8:27;
  step=20;
  tmp=tmp(:,range);
end

tmp=reshape(tmp',1,num*step);
return;
%%%%%%%%%%%%%old


for i=1:loop
    tmp=tmp+wf(:,1+(i-1)*step:i*step);
end

TMP=zeros(1,step);
for j=1:num
%  tmp(j,:)=tmp(j,:)/loop;
  TMP=TMP+tmp(j,:)/loop;
end

TMP=TMP/num;

tmp2=zeros(1,num*step);
for j=1:num
  tmp2(1+(j-1)*step:j*step)=tmp(j,:);
end
%tmp= remmean(TMP);

%tmp=TMP;

%maxv=max(tmp);
%minv=min(tmp);


%if maxv < -minv
%  maxv=-minv;
%end

%tmp=tmp/maxv;

return;
