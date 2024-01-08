%
%pairs=[1 2;2 1];
%
function [ensemble, msT, fstart, Traj] = jointTrout(msT1, fstart1, Traj1, ensemble1, msT2, fstart2, Traj2, ensemble2, pairs)
   [x,y]=size(msT1);
   if x < y
        msT1=msT1';
   end
   [x,y]=size(msT2);
   if x < y
        msT2=msT2';
   end
    msT = [msT1; msT2 - fstart2 + msT1(end)];
   
    fstart = fstart1;
    Traj = [Traj1; Traj2];

    ensemble = cell(size(pairs, 1), 3);

    for i = 1:size(pairs, 1)
        ensemble{i, 1} = [ensemble1{pairs(i, 1), 1} ensemble2{pairs(i, 2), 1}];
        ensemble{i, 3} = [ensemble1{pairs(i, 1), 3} ensemble2{pairs(i, 2), 3} + msT1(end)];
    end

    ensemble = [ensemble; ensemble1(setdiff(1:size(ensemble1, 1), sort(pairs(:, 1))'), :)];
    ensemble = [ensemble; ensemble2(setdiff(1:size(ensemble2, 1), sort(pairs(:, 2))'), :)];
    return;
