%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FRs,pyr,interneuron]=checkFR(Out)
Th=5;
loop=size(Out,1);

if size(Out,2)>1
  inputType=1;
else
  inputType=0;
end
FRs=[];
interneuron=[];
kHz=25;

E=realmin;
S=realmax;
for i=1:loop
    if inputType
    tmps=Out{i,3};
  else
    tmps=Out{i};
  end
  if ~isempty(tmps)
    if tmps(end) > E
      E=tmps(end);
    end
    if tmps(1)<S
      S=tmps(1);
    end
  end
end  

len=(E-S)/(kHz*1000);

for i=1:loop
  if inputType
    tmps=Out{i,3};
  else
    tmps=Out{i};
  end
    
  if ~isempty(tmps)

    FR=size(tmps,2)/double(len);

    if FR >Th
      fprintf('deleting #%d cell %2.1fHz\n',i,FR);
      interneuron=[interneuron i];
    end
  else
    FR=0;
  end
  
  FRs=[FRs; FR];
end

pyr=setdiff(1:loop,interneuron);

return;