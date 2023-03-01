function [Cs, phis, fs, confCs, phistds, zerosps] = batchMTS(Event, ensemble, seq, tetNum, trialave, preSilentSeq,SeqPoint)
    Hz = 25000;
    params.Fs = 1000; %2500
    params.fpass = [0 120];
    %params.tapers=[20 49];
    %params.tapers=[3 5];
    %params.tapers=[10 19];
    params.tapers = [3.5 6];
    %params.tapers=[15 29];
    params.pad = 1;
    params.err = [1 0.95];
    params.trialave = trialave;

    loop = size(ensemble, 1);
    Cs = cell(loop, 1);
    phis = cell(loop, 1);
    fs = cell(loop, 1);
    confCs = cell(loop, 1);
    phistds = cell(loop, 1);
    zerosps = cell(loop, 1);

    step = Hz ./ params.Fs; 

    Event = decimate(double(Event((tetNum - 1) * 4 + 1, :)), step);
   
    seq = floor(seq ./ step);

    %preSilentPeriod
    pSeq = floor(preSilentSeq./step);
    

    for i = 1:loop
        %for i=1
        fprintf('loop=%d/%d\n', i, loop);
        [Cs{i}, phis{i}, fs{i}, confCs{i}, phistds{i}, zerosps{i}] = mtanalysiscpt(Event, tetNum, ensemble{i}, params, seq, step, pSeq,SeqPoint);
    end

    return;
