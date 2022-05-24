function region=batchSpectrum(dlfp,Eps)
params.Fs=1000;
params.fpass=[0 50];
params.tapers=[5 9];
tetrode=4;
step=25;


for j=1:size(Eps,2)
  figure;
  fprintf('Eps:%d\n',j);
  loop=size(Eps{j},1);
  region=[];
  Eps{j}=floor(Eps{j}/step);
  for i=1:loop
    region=[region Eps{j}(i,1):Eps{j}(i,2)];
  end
  
%  region=region(1000000:1500000);
  len=(size(dlfp,1)/tetrode);
  for i=1:len
    fprintf('plot:%d\n',i);
    subplot(len/4,4,i);
    [S,f]=mtspectrumc(dlfp(tetrode*(i-1)+1,region),params);
    plot_vector(S,f);
  end
end


return;