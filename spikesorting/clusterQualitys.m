function [is,lr]=clusterQualitys(o,step,kHz,loop)
if nargin==1
    step=32;
    kHz=25;
elseif nargin==3
    loop=size(o,1);
end


is=zeros(loop,2);
lr=zeros(loop,2);

if size(o,2)>3
    ics=1;
else
    ics=0;
end

for i=1:loop
    %for i=1
    %    fprintf('%d/%d\n',i,loop);
  [is(i,1),lr(i,1)]=MahaDist(o,i,1);
  oLen=size(o,1);
  %spike overlapping
  if 0
  for j=1:(oLen-1)
    tmp=floor(o{j,3}/kHz);
    for k=j+1:oLen
      [~,id]=intersect(floor(o{k,3}/kHz),tmp);
      aid=setdiff(1:size(o{k,3},2),id);
      if isinteger(o{k,1})
          if ics
              o{k,4}=getSpikesS(o{k,4},step,aid);
          end
          o{k,1}=getSpikesS(o{k,1},step,aid);
      else
          if ics
              o{k,4}=getSpikesD(o{k,4},step,aid);
          end
	o{k,1}=getSpikesD(o{k,1},step,aid);
      end
      o{k,3}(id)=[];
    end
  end
  end
  %  [is(i,2),lr(i,2)]=MahaDist(o,i,4);%ic

  if size(o,1)>=2
      if ics 
          [is(i,2),lr(i,2)]=MahaDist(o,i,4);% 4->1 ?         
      else
           [is(i,2),lr(i,2)]=MahaDist(o,i,1);% 4->1 ?
      end
  else
    is(1,2)=NaN;
    lr(1,2)=NaN;
  end
end


return;