function np1IMRO(suffix,bankNum)
fname=sprintf('%s-%d.imro',suffix,bankNum);
fid=fopen(fname,'w');
fprintf(fid,'(0,384)');

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

fprintf(fid,'\n');
fclose(fid);
return;