function ind=plotMTS(Cs,phis,fs,pyr,int,confCs,elecNum,N)
%sm=floor(10/1);
sm=1;
%xlabel=[5 10 15 20 25 30 35 40];
xlabel=[ 5 10 15 20 25 30 35 40 45:5:130];
cellNum=[];
if nargin>6
  for i=1:length(N)
    cellNum=[cellNum; find(elecNum==N(i))];
  end
  cellNum=sort(cellNum);
  pyr=intersect(pyr,cellNum);
  int=intersect(int,cellNum);
end

subplot(2,2,1);
indPyr=coreMTS(Cs,fs,pyr,sm,xlabel,confCs);
title('pyramidal cell');

subplot(2,2,2);
indInt=coreMTS(Cs,fs,int,sm,xlabel,confCs);
title('interneuron');


subplot(2,2,3);
coreMTSphi(phis,fs,pyr,sm,xlabel,indPyr);
title('pyramidal cell');

subplot(2,2,4);
coreMTSphi(phis,fs,int,sm,xlabel,indInt);
title('interneuron');

return;
%%%%%%%%%%%%%%%
function indPyr=coreMTS(Cs,fs,celltype,sm,xlabel,confCs)

loop=length(celltype);
Coh=[];
loop2=0;
indPyr=[];
for i=1:loop
  Cohtmp=Cs{celltype(i)};
  indPyr{i}=find(Cohtmp<confCs{celltype(i)});
  Cohtmp(indPyr{i})=NaN;
  Coh=[Coh Cohtmp];
  loop2=loop2+1;
end


%Coh(1:200,:)=0;

%Coh=Coh(1:4000,:);


[M,ind]=max(Coh);
nCoh=Coh./repmat(M,size(Coh,1),1);
[~,ind2]=sort(ind,'descend');

imagesc(nCoh(:,ind2)');
colormap(jet);
hold on;

tmp=nanmean(Coh,2);

%out=tmp./0.3;

%out=tmp./max(tmp).*loop2;
%out=tmp./max(tmp);
%plot(tmp./max(tmp).*loop2,'k');
%out=smooth(out,sm);

sx=fs{1};
xl=[];
%xlabel=0:10:130;
%xlabel=[5 10 15 20 25 30 35 40];
for i=xlabel
  [~,ind]=min(abs(sx-i));
  xl=[xl ind];
end
out=tmp;
minOut=min(out(xl(1):end));
out=out-minOut;
out=out./max(out);
out=out.*loop;

%out=tmp.*length(celltype);
plot(out,'k');


axis([xl(1) xl(end) 1 size(nCoh,2)]);

set(gca,'xtick',xl,'xticklabel',xlabel);
set(gca,'ydir','normal');
return;
%%%%%%%%%%%%%%%
function coreMTSphi(phis,fs,celltype,sm,xlabel,indPyr)
loop=length(celltype);
Phis=[];
loop2=0;
hold on;

for i=1:loop
  tmp=phis{celltype(i)};
  tmp(indPyr{i})=NaN;
  
  Phis=[Phis tmp];
  loop2=loop2+1;
end

%Phis(1:200,:)=0;

%tmp=(nanmean(Phis,2)./pi)*180;
tmp=nanmedian(Phis,2);
%tmp=Phis./pi.*180;

if ~isempty(tmp)
  
%  [~,ind]=max(Phis);
%  [~,ind2]=sort(ind,'descend');
%  imagesc(Phis(:,ind2)');
%colorbar;
%  hold on;
%  plot(fs{1},tmp,'k');

sm=Phis;
for i=1:size(Phis,2)
%  sm(:,i)=smooth(Phis(:,i),100);
%  plot(fs{i},smooth(Phis(:,i),100));
end
plot(fs{i},nanmean(sm,2));
%  plot(fs{1},smooth(tmp,sm),'k');
%  axis([0 40 -180 180]);

end


%xlabel=[ 5 10 15 20 25 30 35 40];

axis([xlabel(1) xlabel(end) -pi pi]);
set(gca,'xtick',xlabel,'xticklabel',xlabel);
set(gca,'ydir','normal');

%save test.mat Phis

return;
sx=fs{1};
xl=[];
xlabel=0:10:130;
for i=xlabel
  [~,ind]=min(abs(sx-i));
  xl=[xl ind];
end
set(gca,'xtick',xl,'xticklabel',xlabel);
set(gca,'ydir','normal');
return;