function batchCondPD(server,method)
if nargin==1
    method='gainmap';
end

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

    for i = 1:size(dataParing, 1)
        [phaseHistPyr, phaseHistInt, PyrIntList, PyrIntListStim, ~, ~, ~, ~, phaseHistPyrCtrl, phaseHistIntCtrl] = contPD(dataParing{i, 1}, dataParing{i, 2}, 'cellClass', 0, 'localcell', 1,'method',method);
        save([dataParing{i,2} '.mat'],'phaseHistPyr','phaseHistInt','PyrIntList','PyrIntListStim','phaseHistPyrCtrl','phaseHistIntCtrl');
    end
