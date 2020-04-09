function occupancy = OmitIslands(occupancy)
% Takes matrix of occupancy values, calculates center of mass of pixels>0,
% and then finds all disconnected pixels and sets them = 0

%subplot(1, 2, 1), imagesc(occupancy)

s = regionprops(occupancy>0, {'FilledArea', 'PixelIdxList'});

l = numel(s);

areas = vertcat(s.FilledArea);

[~, inds] = sort(areas);

for i = 1:length(inds)-1
    
    occupancy(s(inds(i)).PixelIdxList) = 0;
    
end 

%subplot(1, 2, 2), imagesc(occupancy)

return;