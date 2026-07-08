function parsave(fn,varargin)
loop=size(varargin,2);
for i=1:2:loop
  eval(  [varargin{i} '= varargin{i+1};']);
  if i==1
    save(fn,varargin{i},'-v7.3');
  else
    save(fn,varargin{i},'-append','-v7.3');
  end
end


return;