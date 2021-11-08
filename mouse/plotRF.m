function [vars,labels]=plotRF(segPara)
num=2;%maxDist
loop=size(segPara,1);

mm=[];
for i=1:loop
  mm=[mm nanmedian(segPara{i,num})];
end

plot(mm);

xticklabels={'pre','peak','s','trough','s','rising','s','falling','post'};
set(gca,'xticklabel',xticklabels);


for i=1:loop
  segPara{i,num}(isnan(segPara{i,num}))=[];
end



figure;
%stats
labels=[];
vars=[];
for i=1:loop
  initC=num2str(i);
  labels=[labels repmat({initC},1,size(segPara{i,num},2))];
  vars=[vars segPara{i,num}];
end

boxplot(vars,labels,'notch','on');
xticklabels={'pre','peak','s','trough','s','rising','s','falling','post'};
set(gca,'xticklabel',xticklabels);

%[p,tbl,stats]=kruskalwallis(vars,labels);
[p,tbl,stats]=anova1(vars,labels);
r=multcompare(stats);

%anova1(vars,labels);
return;