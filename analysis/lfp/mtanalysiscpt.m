function [Cs, phi, f, confC, phistd, zerosp] = mtanalysiscpt(data, tetNum, unit, params, seq, step)
    Hz = 25000;
    movingwin = [1 0.5];
    unit = double(unit) ./ step;

    loop = size(seq, 2);

    if loop == 2
        loop = 1;
    end

    duration = min(floor(diff(seq) / 1000) * 1000);
    unitDuration = duration * step;

    if 0 %trial average
        Data1 = zeros(duration, loop);

        for i = 1:loop
            Data1(:, i) = data(seq(i) + 1:seq(i) + duration);
            %  spk=unit(unit>seq(i) & unit<seq(i)+unitDuration);
            spk = unit(unit > seq(i) & unit < seq(i) + duration);

            %  Spk(i).times=(double(spk')-seq(i))./Hz;
            Spk(i).times = (double(spk') - seq(i)) ./ params.Fs;
        end

        [C, phi, S12, S1, S2, f, zerosp, confC, phistd] = coherencycpt(Data1, Spk, params, 1);
    else %segmentation every 1sec
        win = 1; %segment duration in sec
        Data1 = data(seq(1):seq(end) + duration)';
        spk = unit(unit > seq(1) & unit < seq(end) + duration);
        Spk.times = (double(spk') - seq(1)) ./ params.Fs; %convert from kHz to Hz

        [C, phi, S12, S1, S2, f, zerosp, confC, phistd] = coherencysegcpt(Data1, Spk, win, params, 0);
    end

    [Cs,phi] = calcCohere(S12,S1,S2, 1); %ascending
    [Ctmp,phiTmp]=calcCohere(S12,S1,S2,2);%descending
    Cs=cat(3,Cs,Ctmp);
    phi=cat(3,phi,phiTmp);
    [Ctmp,phiTmp]=calcCohere(S12,S1,S2,3);%both
    Cs=cat(3,Cs,Ctmp);
    phi=cat(3,phi,phiTmp);

    %[C,phi,S12,S1,S2,t,f,zerosp,confC,phistd]=cohgramcpt(Data1,Spk,movingwin,params,1);

    for i = 1:size(phistd, 2)
        nanRange = find(phistd(:, i) > 2);
        %  phi(nanRange,i)=NaN;
        %  C(nanRange,i)=NaN;
    end

end

%%%%
function [C, phi] = calcCohere(S12,S1,S2, div)

    S12 = squeeze(Coh(:, :, 2));
    S1 = squeeze(Coh(:, :, 3));
    S2 = squeeze(Coh(:, :, 4));
    S12 = getDir(S12, div);
    S1 = getDir(S1, div);
    S2 = getDir(S2, div);
    S12 = squeeze(nanmean(S12, 2));
    S1 = squeeze(nanmean(S1, 2));
    S2 = squeeze(nanmean(S2, 2));

    Coh = S12 ./ sqrt(S1 .* S2);
    C = abs(Coh);
    phi = angle(Coh);
end

%%%%%%
function S = getDir(S, div)
    segment = 60;
    repeat = 20;

    if div > 0 % 1:ascending, 2:descending
        S = reshape(S, size(S, 1), segment, repeat);

        switch (div)
            case 1,
                dir = 1:2:20;
            case 2,
                dir = 2:2:20;
            case 3,
                dir = 1:20;
        end

        S = S(:, :, dir);
        S = reshape(S, size(S, 1), segment * length(dir));

    end

end
