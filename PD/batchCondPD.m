function batchCondPD(server, varargin)
    p = inputParser;
    p.addParamValue('method', 'gainmap', @ischar);
    p.addParamValue('wavetype', 'chirp', @ischar);
    p.addParamValue('rawlfp', 0, @isnumeric); %0: not used, 1: left lfp, 2: right lfp

    p.parse(varargin{:});
    method = p.Results.method;
    waveType = p.Results.wavetype;
    rawlfp = p.Results.rawlfp;


    switch (server)
        case 'deepMachine',
            dataParing = {
                          '2018', 'TK18052901'
                          '2018', 'TK18052902'
                          '2018', 'TK18121401'
                          '2018', 'TK19011801'
                          '2018', 'TK19012501'
                          '2019', 'TK19040101'
                          '2019', 'TK19050801'
                          '2019', 'TK19060301'
                          };
        case 'windows',
            dataParing = {
                          '2019f', 'AZ19120501'
                          '2019f', 'AZ20013101'
                          '2019f', 'AZ20031201'
                          '2020h', 'AZ20030601'
                          '2020h', 'AZ20031901'
                          '2020h', 'AZ20070101'
                          '2020e', 'AZ20121501'
                          '2020e', 'AZ21010701'
                          };
    end

    %for i=7

    for i = 1:size(dataParing, 1)
        [phaseHistPyr, phaseHistInt, PyrIntList, PyrIntListStim, fr, tp, sw, pi, phaseHistPyrCtrl, phaseHistIntCtrl, cq, pyr, interneuron] = contPD(dataParing{i, 1}, dataParing{i, 2}, 'cellClass', 0, 'localcell', 1, 'method', method, 'waveType', waveType, 'rawlfp', rawlfp);

        if strcmp(method, 'cellclassify')
            save([dataParing{i, 2} '.mat'], 'tp', 'sw', 'pi', 'cq');
        else
            save([dataParing{i, 2} '.mat'], '-v7.3', 'phaseHistPyr', 'phaseHistInt', 'PyrIntList', 'PyrIntListStim', 'phaseHistPyrCtrl', 'phaseHistIntCtrl', 'pyr', 'interneuron', 'cq');
        end

    end
