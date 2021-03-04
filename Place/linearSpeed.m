function [speed,LTraj]=linearSpeed(LTraj,msT)
[speed,cmPerPixel]=createSpeed(LTraj, msT,0,'animal','rat');
LTraj=ceil(LTraj*cmPerPixel);
[~,ind]=sort(LTraj);

dim=0:400;
occupancy=histcounts(LTraj,dim);

for i=dim(2:end)
    speedC(i)=sum(speed(LTraj==i));
end

%speedC=histcounts(speedC,dim);
speed=speedC./occupancy;

speed=SmoothMat(speed,[2 2],2.5);
%a=histogram(speed(ind),LTraj(ind));
return;
%%%%%%%%%%%%%%%%%%%%%%%%%5555555555
function mat = SmoothMat(mat, kernel_size, std)
if nargin<3
    std=1;
end

if std == 0, return; end


Rgrid=-kernel_size/2: kernel_size/2;
kernel=normpdf(Rgrid,0,std);

kernel = kernel./sum(sum(kernel));

mat = conv(mat, kernel, 'same');
return;