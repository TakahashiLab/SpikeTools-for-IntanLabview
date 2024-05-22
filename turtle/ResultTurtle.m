function ResultTurtle(cc,sb,Sn,HD,Sn2,HD2,FR,FR2)
if nargin==0
    fprintf('loading data directly\n');
    % Load the first set of variables
    load('O:\seaturtle\2023\combine\alpha\20230820\behav.mat');

    % Store the variables and their values from the first dataset
    data_alpha = struct;
    vars_alpha = who;
    for i = 1:length(vars_alpha)
        data_alpha.(vars_alpha{i}) = eval(vars_alpha{i});
    end

    % Load the second set of variables
    load('O:\seaturtle\2023\combine\gamma\20230821-1\behav.mat');

    % Concatenate each variable from both datasets if they have the same name, using transpose if needed
    for i = 1:length(vars_alpha)
        var_name = vars_alpha{i};



        if eval(['exist(''' var_name ''', ''var'')']) % Check if variable exists in the second dataset
            alpha_var = data_alpha.(var_name);
            gamma_var = eval(var_name);


            if strmatch(var_name,'cc')

                gamma_var = gamma_var'; % Transpose gamma_var

                alpha_var = alpha_var'; % Transpose alpha_var if gamma_var cannot be transposed

            end

            if strmatch(var_name,'data_alpha')
            else
                % Concatenate variables
                eval([var_name ' = [alpha_var; gamma_var];']);
            end
        end
    end
elseif nargin==1
    switch (cc)
        case 1,
            load('O:\seaturtle\2023\combine\alpha\20230820\behav.mat');
        case 2,
            load('O:\seaturtle\2023\combine\gamma\20230821-1\behav.mat');
    end
end

hdcells=find(any(Sn'));
%hdcells=find(sum(Sn')==2);

subplot(4,4,1)
histogram(cc(hdcells),-1:0.1:1);
c=cc(hdcells);
title(sprintf('skew=%1.2f',skewness(c)));
xlabel('correlation')

subplot(4,4,2)
histogram(rad2deg(HD(hdcells,1)),-180:6:180)
set(gca,'xtick',[-180 -90 0 90 180],'xticklabel',{'W','S','E','N','W'})
title('ang HD State norm')
rad2deg(circ_mean(HD(hdcells,1)))

subplot(4,4,3)
histogram(rad2deg(HD(hdcells,2)),-180:6:180)
set(gca,'xtick',[-180 -90 0 90 180],'xticklabel',{'W','S','E','N','W'})
title('ang HD State zero')
rad2deg(circ_mean(HD(hdcells,2)))

subplot(4,4,4)
fr=FR(hdcells,:)';
histogram(diff(fr)./sum(fr),-1:0.1:1);
c=diff(fr)./sum(fr);
title(sprintf('skew=%1.2f',skewness(c)));
xlabel('diff/sum (+:norm < zero)')

if ~isempty(sb)
    subplot(4,4,5)
    histogram(sb(hdcells,1),-1:0.1:1);
    c=sb(hdcells,1);
    title(sprintf('skew=%1.2f',skewness(c)));
    xlabel('cc states 1-2');

    subplot(4,4,6)
    histogram(sb(hdcells,2),-1:0.1:1);
    c=sb(hdcells,2);
    title(sprintf('skew=%1.2f',skewness(c)));
    xlabel('cc states 3-4');
end

if size(sb,2)>2
    subplot(4,4,7)
    histogram(sb(hdcells,3),-1:0.1:1)
    c=sb(hdcells,3);
    title(sprintf('skew=%1.2f',skewness(c)));
    xlabel('cc states 1-3');

    subplot(4,4,8)
    histogram(sb(hdcells,4),-1:0.1:1)
    c=sb(hdcells,4);
    title(sprintf('skew=%1.2f',skewness(c)));
    xlabel('cc states 2-4');
end

if ~isempty(HD2)
    subplot(4,4,9)
    histogram(rad2deg(HD2(find(any(Sn2(:,[1 3])')),1)),-180:6:180)
    title('ang HD State norm 1')
    set(gca,'xtick',[-180 -90 0 90 180],'xticklabel',{'W','S','E','N','W'})

    subplot(4,4,10)
    histogram(rad2deg(HD2(find(any(Sn2(:,[1 3])')),3)),-180:6:180)
    title('ang HD State norm 2')
    set(gca,'xtick',[-180 -90 0 90 180],'xticklabel',{'W','S','E','N','W'})

    subplot(4,4,11)
    histogram(rad2deg(HD2(find(any(Sn2(:,[2 4])')),2)),-180:6:180)
    title('ang HD State zero 1')
    set(gca,'xtick',[-180 -90 0 90 180],'xticklabel',{'W','S','E','N','W'})

    subplot(4,4,12)
    histogram(rad2deg(HD2(find(any(Sn2(:,[2 4])')),4)),-180:6:180)
    title('ang HD State zero 2')
    set(gca,'xtick',[-180 -90 0 90 180],'xticklabel',{'W','S','E','N','W'})


    subplot(4,4,13)
    fr=FR2(find(any(Sn2(:,[1 3])')),[1 3])';
    histogram(diff(fr)./sum(fr),-1:0.1:1);
    c=diff(fr)./sum(fr);
    title(sprintf('skew=%1.2f',skewness(c)));
    xlabel('diff/sum (+:norm < zero)')

    subplot(4,4,14)
    fr=FR2(find(any(Sn2(:,[2 4])')),[2 4])';
    histogram(diff(fr)./sum(fr),-1:0.1:1);
    c=diff(fr)./sum(fr);
    title(sprintf('skew=%1.2f',skewness(c)));
    xlabel('diff/sum (+:norm < zero)')
end

return;