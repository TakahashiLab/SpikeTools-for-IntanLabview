function datr=gpuFiltFilt(b1,a1,buff)

dataRAW = gpuArray(buff); % move int16 data to GPU
dataRAW = dataRAW';
dataRAW = single(dataRAW); % convert to float32 so GPU operations are fast

datr = filter(b1, a1, dataRAW); % causal forward filter
datr = flipud(datr); % reverse time
datr = filter(b1, a1, datr); % causal forward filter again
datr = flipud(datr); % reverse time back

datr=gather(int16(datr));