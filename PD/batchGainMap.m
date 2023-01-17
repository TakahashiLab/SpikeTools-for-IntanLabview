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
    Data1 = Event(seq(1):seq(end) + duration)';

    %if duration>20000
    if 0 %for noise?
        xphase0 = analogPhase(Data1(1:duration), 1.5)';
        xphase = repmat(xphase0, 1, 20);
    else
        xphase = chirpPhase(Data1);
        %    xphase = mod(xphase,2*pi);%covert to 0-2pi rather than -pi:pi
    end

    %xphase0=analogPhase(Data1(1:60000),1.5)';

    phaseHistPyr = calcGM(xphase, params, ensemble(Pyr), step, seq, normSeq, duration);
    phaseHistInt = calcGM(xphase, params, ensemble(Int), step, seq, normSeq, duration);
    phaseHistPyrCtrl = calcGM4Ctrl(xphase, params, ensemble(Pyr), step, normSeq);
    phaseHistIntCtrl = calcGM4Ctrl(xphase, params, ensemble(Int), step, normSeq);

    return;
    %%%%%%%%%%%%%%%%%%%%%%%%%
    function phaseHist = calcGM(xphase, params, ensemble, step, seq, normSeq, duration)

        segmentSec = 0.5;
        segment = segmentSec * params.Fs; %0.5sec

        loop = size(ensemble, 1);
        phaseHist = zeros(20, 120, loop);

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

            for j = 1:length(seq)
                %    spk=unit(unit>seq(j) & unit<seq(j)+duration);
                %    spk=(double(spk')-seq(j))./params.Fs;
                c = 0;
                alpha = 100;

                while 1 %frequency
                    loc = seq(j) + c * segment;
                    %if loc+segment+alpha > seq(j)+duration
                    if loc + segment > seq(j) + duration
                        break;
                    end

                    spk = unit(unit > loc & unit < loc + segment);
                    spk = spk - seq(1);
                    xphasePos = xphase;
                    %    phaseHist(:,c+1,i)=hist(xphasePos(spk),-pi:pi/10:pi);
                    %phase

                    phaseHist(:, c + 1, i) = hist(xphasePos(spk), -pi:pi / 9.5:pi) ./ (segmentSec / 20); %Hz
                    c = c + 1;
                end

            end

            phaseHist(:, :, i) = phaseHist(:, :, i) ./ nFr(i);

        end

        return;
        %%%%%%%%%%%%%%%%%%%%%%%%%
        function phaseHist = calcGM4Ctrl(xphase, params, ensemble, step, normSeq)
     
            segmentSec = 0.5;
            segment = segmentSec * params.Fs; %0.5sec

            loop = size(ensemble, 1);
            phaseHist = zeros(20, 120, loop);

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
                for j = 1:10
               
                    c = 0;
                  
                    for c=0:(120-1)%frequency
                        loc =  loc0 + c * segment;
                        
                        spk = unit(unit > loc & unit < loc + segment);
                        spk = spk - normSeq(1);
                        xphasePos = xphase;

                        %phase
                        phaseHist(:, c + 1, i) = hist(xphasePos(spk), -pi:pi / 9.5:pi) ./ (segmentSec / 20); %Hz
                        c = c + 1;
                    end
                    loc0 = loc; 
                end

                phaseHist(:, :, i) = phaseHist(:, :, i) ./ nFr(i);

            end

            return;
