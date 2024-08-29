%bank-0: 0-383
%bank-1: 384-768
%bank-2: 769-960 0-191
%
function np1IMRO(suffix,bankNum,map)
if nargin>2
    fname=sprintf('%s-map.imro',suffix);
else    
    fname=sprintf('%s-%d.imro',suffix,bankNum);
end
fid=fopen(fname,'w');
fprintf(fid,'(0,384)');

%manual mapping 
if nargin>2
    for i=1:size(map,1)
        for j=map(i,1):map(i,2)
            fprintf(fid,'(%d %d 0 500 250 0)',i,map(i,3));            
        end
    end
else%systematic 
    if bankNum~=2
        for i=0:383
            fprintf(fid,'(%d %d 0 500 250 0)',i,bankNum);
        end
    elseif bankNum==2
        for i=0:191
            fprintf(fid,'(%d %d 0 500 250 0)',i,bankNum);
        end
        for i=[192:383]
            fprintf(fid,'(%d %d 0 500 250 0)',i,bankNum-1);    
        end
    end
end
    

fprintf(fid,'\n');
fclose(fid);
return;