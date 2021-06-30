function [out,b]=countPlaceFields(rate_map,binside,ThF)

if nargin==2
    ThF=50;    
end
out=[];
xy=[];
verbose=0;

rm = rate_map;
binside2=binside^2;

peakRate=max(rm(:));
if peakRate>1.0
    sp_thr = peakRate/5;

    field = (rm>=sp_thr);

    a = regionprops(field,{'Area' 'PixelList'});

    b = [];

    c=1;
    for i = 1:length(a)
        if ((a(i).Area * binside2) > ThF)
            b{c} = a(i).PixelList;
            out=[out a(i).Area * binside2];
            c=c+1;
        end
    end

    if verbose
        imagesc(rate_map);colormap jet;
        for k=1:(c-1)
            xy=b{k};
            hold on;
            plot(xy(:,1),xy(:,2),'.');
        end
    end
end

return;