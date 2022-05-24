function plotGainMap(phaseHist)

sd=10;
edges=[-3*sd:5:3*sd];
kernel=normpdf(edges,0,sd);sd=10;
edges=[-3*sd:5:3*sd];
kernel=normpdf(edges,0,sd)*sd;
%imagesc(conv2(kernel,kernel,mean(phaseHist(:,[1:40 80:120],:),3)','same'));
%mph=conv2(kernel,kernel,mean(phaseHist(:,[1:40 81:120],:),3)','same');
subplot(1,2,1);
mph=conv2(kernel,kernel,mean(phaseHist(:,[1:40],:),3)','same');

imagesc([mph(:,11:20) mph mph(:,1:10)]);
title('ascending');
set(gca,'ydir','normal');
set(gca,'ytick',[1:5:40],'yticklabel',[1:5:40]*2);
set(gca,'xtick',[1 11 21 31 41],'xticklabel',{'-2\pi','-\pi','0','\pi','2\pi'});
hold on;
t=0:40;
plot(t+1,cos(t/pi)*20+20,'w');
colorbar;
xlabel('Phase [radian]');
ylabel('Frequency [Hz]');

subplot(1,2,2);
%mph=conv2(kernel,kernel,mean(phaseHist(:,[81:120],:),3)','same');
mph=conv2(kernel,kernel,mean(phaseHist(:,[120:-1:81],:),3)','same');
imagesc([mph(:,11:20) mph mph(:,1:10)]);
title('descending');
set(gca,'ydir','normal');
%set(gca,'ytick',[1:5:40],'yticklabel',[(40:-5:1)-4]*2);
set(gca,'ytick',[1:5:40],'yticklabel',[1:5:40]*2);
set(gca,'xtick',[1 11 21 31 41],'xticklabel',{'-2\pi','-\pi','0','\pi','2\pi'});
hold on;
t=0:40;
plot(t+1,cos(t/pi)*20+20,'w');
xlabel('Phase [radian]');
ylabel('Frequency [Hz]');
colormap(jet);
colorbar;

return;
