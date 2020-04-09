function ks=convertKS(fn,sampl)
if nargin==1
    sampl=25000;
end
ks=loadKSdir(fn);
ks.st=ceil(ks.st*sampl);
return;