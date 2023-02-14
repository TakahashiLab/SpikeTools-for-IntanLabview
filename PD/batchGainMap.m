function [phaseHistPyr, phaseHistInt, phaseHistPyrCtrl, phaseHistIntCtrl] = batchGainMap(Event, ensemble, seq, tetNum, normSeq, Pyr, Int)
    %function [Data1]=batchGainMap(Event,ensemble,seq,tetNum)
    Hz = 25000;
    params.Fs = 1000; %2500

   step = Hz ./ params.Fs;
  
    Event = decimate(double(Event((tetNum - 1) * 4 + 1, :)), step);
    %Event=double(Event((tetNum-1)*4+1,1:step:end));

    seq = floor(seq ./ step);
    normSeq = floor(normSeq ./ step);

    duration = min(floor(diff(seq) / 1000) * 1000);
  
    Event = Event(seq(1):seq(end) + duration)';

    if size(seq, 2) == 30 %for noise
        Data4phase = hilbert(Event);
        xphase = angle(Data4phase);
    elseif size(seq, 2) == 20 %chirp
        
        xphase = chirpPhase(Event);
        
    elseif size(seq, 2) == 10 %pulse

    end
    

    %xphase0=analogPhase(Data1(1:60000),1.5)';
    phaseCnt = 20; %18 degree?
    phaseHistPyr = calcGM(xphase, params, ensemble(Pyr), step, seq, normSeq, duration, phaseCnt);
    phaseHistInt = calcGM(xphase, params, ensemble(Int), step, seq, normSeq, duration, phaseCnt);
    phaseHistPyrCtrl = calcGM4Ctrl(xphase, params, ensemble(Pyr), step, normSeq, phaseCnt);
    phaseHistIntCtrl = calcGM4Ctrl(xphase, params, ensemble(Int), step, normSeq, phaseCnt);

    return;
    %%%%%%%%%%%%%%%%%%%%%%%%%
    function phaseHist = calcGM(xphase, params, ensemble, step, seq, normSeq, duration, phaseCnt)

        segmentSec = 0.5;
        segment = segmentSec * params.Fs; %0.5sec

        edges = -pi:pi / (phaseCnt / 2):pi;

        loop = size(ensemble, 1);
        phaseHist = zeros(phaseCnt, 120, loop);

        nFr = zeros(loop, 1);

        for i = 1:loop %unit
            fprintf('loop=%d/%d\n', i, loop);
            unit = ensemble{i};
            unit = round(unit ./ step);

            %calculate regular firing rate during pre and post silent periods
            for k = 1 %pre silent preriod only, if both pre and post periods, then 1:2
                nFr(i) = nFr(i) + sum(unit > normSeq(k, 1) & unit < normSeq(k, 2));
            end

            dnFr = sum(diff(normSeq')) ./ params.Fs;
            nFr(i) = nFr(i) ./ dnFr; %Hz

            for j = 1:length(seq)
                %    spk=unit(unit>seq(j) & unit<seq(j)+duration);
                %    spk=(double(spk')-seq(j))./params.Fs;
                c = 0;
                alpha = 100;

                while 1 %frequency
                    loc = seq(j) + c * segment;
                    phaseLoc = 1 + c * segment:(c + 1) * segment;

                    %if loc+segment+alpha > seq(j)+duration
                    if loc + segment > seq(j) + duration
                        break;
                    end

                    spk = unit(unit >= loc & unit < loc + segment);
                    spk = spk - seq(1) + 1;
                    xphasePos = xphase;

                    if ~isempty(spk)

                        phaseHist(:, c + 1, i) = phaseHist(:, c + 1, i) + (histcounts(xphasePos(spk), edges)' ./ (segmentSec / phaseCnt)); %Hz
                        %phaseHist(:, c + 1, i) = hist(xphasePos(spk), -pi:pi / 9.5:pi) ./ (segmentSec / 20); %Hz
                    end

                    c = c + 1;
                end

            end

            %average
            %phaseHist(:, :, i) = phaseHist(:, :, i) ./ length(seq) ./ nFr(i);
            phaseHist(:, :, i) = phaseHist(:, :, i) ./ length(seq); %length(seq)= optostim cycles

        end

        return;
        %%%%%%%%%%%%%%%%%%%%%%%%%
        function phaseHist = calcGM4Ctrl(xphase, params, ensemble, step, normSeq, phaseCnt)

            segmentSec = 0.5;
            segment = segmentSec * params.Fs; %0.5sec

            edges = -pi:pi / (phaseCnt / 2):pi;

            loop = size(ensemble, 1);
            phaseHist = zeros(phaseCnt, 120, loop);

            nFr = zeros(loop, 1);

            for i = 1:loop %unit
                fprintf('loop=%d/%d\n', i, loop);
                unit = ensemble{i};
                unit = floor(unit ./ step);

                %calculate regular firing rate during pre and post silent periods
                for k = 1 %pre silent preriod only, if both pre and post periods, then 1:2
                    nFr(i) = nFr(i) + sum(unit > normSeq(k, 1) & unit < normSeq(k, 2));
                end

                dnFr = sum(diff(normSeq')) ./ params.Fs;
                nFr(i) = nFr(i) ./ dnFr; %Hz

                loc0 = normSeq(1);
                xphasePos = xphase;

                for j = 1:10

                    for c = 0:(120 - 1) %frequency
                        loc = loc0 + c * segment;
                        phaseLoc = 1 + c * segment:(c + 1) * segment;

                        spk = unit(unit >= loc & unit < loc + segment);
                        spk = spk - loc0 + 1;

                        if ~isempty(spk)
                            phaseHist(:, c + 1, i) = phaseHist(:, c + 1, i) + (histcounts(xphasePos(spk), edges)' ./ (segmentSec / phaseCnt)); %Hz
                            %phaseHist(:, c + 1, i) = basehist; %Hz
                        end

                    end

                    loc0 = 120 * segment;
                end

                %average
                %phaseHist(:, :, i) = phaseHist(:, :, i) ./ 10 ./ nFr(i);
                phaseHist(:, :, i) = phaseHist(:, :, i) ./ 10; %10 phase cycles

            end

            return;
