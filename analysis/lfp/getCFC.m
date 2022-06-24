function out=getCFC(xphase,xamp)

%xphase=xphase([2 3],:);
%xamp=xamp([2 3],:); 

nbins=18;
nsurrogates=1;
randtype=2;
%out = get_mi(xphase,xamp,nbins,nsurrogates,randtype);
out = get_mi(xphase,xamp,nbins);

return;


sampl=1000;
step=25;
unfilteredtrace=decimate(double(unfilteredtrace),step);

%delta
ft(1,:)=filterX(unfilteredtrace,0.1,4,sampl);
%theta
ft(2,:)=filterX(unfilteredtrace,4,12,sampl);
%beta
ft(3,:)=filterX(unfilteredtrace,13,30,sampl);
%lowgamma
ft(4,:)=filterX(unfilteredtrace,30,60,sampl);
%highgamma
ft(5,:)=filterX(unfilteredtrace,60,120,sampl);

for i=1:5
  Hx(i,:)=hilbert(ft(i,:));
end

xphase=angle(Hx);
xamp=abs(Hx);


