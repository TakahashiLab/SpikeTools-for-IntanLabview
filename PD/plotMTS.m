function [thetaIdx, peakCoh, peakFreq] = plotMTS(Cs, phis, fs, pyr, int, confCs, Xlabel, topFreq, div, elecNum, N)
    %sm=floor(10/1);
    sm = 1;
    %xlabel=[5 10 15 20 25 30 35 40];
    %xlabel = [5 10 15 20 25 30 35 40 45:5:130];
    cellNum = [];

    if nargin > 9
        hoge

        for i = 1:length(N)
            cellNum = [cellNum; find(elecNum == N(i))];
        end

        cellNum = sort(cellNum);
        pyr = intersect(pyr, cellNum);
        int = intersect(int, cellNum);
    end

    %pyramidal cell
    subplot(5, 2, 1);
    [out, xl, indPyr, thetaIdx{1}, spN, peakCoh{1}, peakFreq{1}, Rc{1}] = coreMTS(Cs, fs, pyr, sm, Xlabel, confCs, topFreq, div);
    title('pyramidal cell');
    fprintf('pyramidal cell: theta index=%f\n', nanmean(thetaIdx{1}));
    subplot(5, 2, 3);
    plot(out, 'k');
    hold on;
    spN(spN > length(out)) = [];
    plot(spN, out(spN), 'r');
    axis([xl(1) xl(end) 0 max(out(:))])
    set(gca, 'xtick', xl, 'xticklabel', Xlabel);
    set(gca, 'ydir', 'normal');

    %interneuron
    subplot(5, 2, 2);
    [out, ~, indInt, thetaIdx{2}, spP, peakCoh{2}, peakFreq{2}, Rc{2}] = coreMTS(Cs, fs, int, sm, Xlabel, confCs, topFreq, div);
    title('interneuron');
    fprintf('interneuron: theta index=%f\n', nanmean(thetaIdx{2}));
    spP(spP > length(out)) = [];
    subplot(5, 2, 4);
    plot(out, 'k');
    hold on;
    plot(spP, out(spP), 'r')
    axis([xl(1) xl(end) 0 max(out(:))])
    set(gca, 'xtick', xl, 'xticklabel', Xlabel);
    set(gca, 'ydir', 'normal');

    subplot(5, 2, 5);
    coreMTSphi(phis, fs, pyr, sm, Xlabel, indPyr, div, xl, spN);
    title('pyramidal cell');

    subplot(5, 2, 6);
    coreMTSphi(phis, fs, int, sm, Xlabel, indInt, div, xl, spP);
    title('interneuron');

    subplot(5, 2, 7);
    plot(Rc{1}{1}, Rc{1}{2}, '.');
    hold on;
    xlabel('Time correlation');
    ylabel('Frequency correlation');
    plot([-1 1], [0 0], [0 0], [-1 1], '-.');
    axis([-1 1 -1 1]);
    p = signrank(Rc{1}{1}, Rc{1}{2});
    fprintf('Pyramidal cell signed rank=%f\n', p);

    subplot(5, 2, 9);
    plot((cell2mat(Rc{1})'));
    set(gca, 'xtick', [1 2], 'xticklabel', {'Time', 'Freq'});
    ylabel('Rank correlation');
    axis([0.5 2.5 -1 1]);

    subplot(5, 2, 8);
    plot(Rc{2}{1}, Rc{2}{2}, '.');
    hold on;
    xlabel('Time correlation');
    ylabel('Frequency correlation');
    plot([-1 1], [0 0], [0 0], [-1 1], '-.');
    axis([-1 1 -1 1]);
    p = signrank(Rc{2}{1}, Rc{2}{2});
    fprintf('Interneuron signed rank%f\n', p);

    subplot(5, 2, 10);
    plot((cell2mat(Rc{2})'));
    set(gca, 'xtick', [1 2], 'xticklabel', {'Time', 'Freq'});
    ylabel('Rank correlation');
    axis([0.5 2.5 -1 1]);
    set(gcf, 'position', [100 100 400 800]);
end

%%%%%%%%%%%%%%%
function [out, xl, indPyr, thetaIdx, sp, peakCoh, peakFreq, Rc] = coreMTS(Cs, fs, celltype, sm, xlabel, confCs, topFreq, div)

    loop = length(celltype);
    Coh = [];
    Ctrl = [];
    loop2 = 0;
    indPyr = [];
    CohAS = [];
    CohDS = [];

    for i = 1:loop
        Cohtmp = Cs{celltype(i)};

        %CtrlTmp = squeeze(Cohtmp(:, :, div + 3));
        CtrlTmp = squeeze(Cohtmp(:, :, 6));
        CohtmpAS = squeeze(Cohtmp(:, :, 1));
        CohtmpDS = squeeze(Cohtmp(:, :, 2));
        Cohtmp = squeeze(Cohtmp(:, :, div));

        Coh = [Coh Cohtmp];
        Ctrl = [Ctrl CtrlTmp];
        %rank correlation
        CohAS = [CohAS CohtmpAS];
        CohDS = [CohDS CohtmpDS];
    end

    sx = fs{1};
    xl = [];

    for i = xlabel
        [~, ind] = min(abs(sx - i));
        xl = [xl ind - 1];
    end

    endA = find(xlabel == topFreq);

    [peakCoh, peakFreq] = max(Coh(1:xl(endA), :));

    peakFreq = sx(peakFreq);

    %rank correlation
    CohAS = CohAS(1:xl(endA), :);
    CohDS = CohDS(1:xl(endA), :);

    Frc = diag(corr(CohAS, CohDS, 'type', 'kendall')); %frequency
    Trc = diag(corr(CohAS, CohDS(end:-1:1, :), 'type', 'kendall')); %time

    Rc{1} = Trc;
    Rc{2} = Frc;

    [M, ind] = max(Coh(1:xl(endA), :));
    nCoh = Coh ./ repmat(M, size(Coh, 1), 1);
    [~, ind2] = sort(ind, 'ascend');
    tmp = nCoh(1:xl(endA), ind2);
    imagesc(tmp');
    colormap(jet);
    hold on;

    tmp = nanmean(Coh(1:xl(endA), :), 2);

    out = tmp;
    startT = find(xlabel == 4);
    endT = find(xlabel == 12);
    thetaIdx = xl(startT):xl(endT);

    startA = find(xlabel == 2);
    allIdx = xl(startA):xl(endA);
    thetaIdx = nanmean(Coh(thetaIdx, :), 1) ./ nanmean(Coh(allIdx, :), 1);

    if ~isempty(nCoh)
        axis([xl(1) xl(end) 1 size(nCoh, 2)]);
    end

    set(gca, 'xtick', xl, 'xticklabel', xlabel);
    set(gca, 'ydir', 'normal');

    %permutation test
    startA = find(xlabel == 2);
    allIdx = xl(startA):xl(endA);
    alpha = 0.05;

    Coh = Coh(allIdx, :);
    Ctrl = Ctrl(allIdx, :);
    [clusters, p] = permutest(Coh, Ctrl, false, alpha, 5000, false);

    sp = cell2mat(clusters(p < alpha));

end

%%%%%%%%%%%%%%%
function coreMTSphi(phis, fs, celltype, sm, xlabel, indPyr, div, xl, sp)
    loop = length(celltype);
    Phis = [];
    loop2 = 0;

    for i = 1:loop
        tmp = phis{celltype(i)};

        tmp = squeeze(tmp(:, :, div));

        Phis = [Phis tmp];
        loop2 = loop2 + 1;
    end

    %[~,Neu]=find(isnan(Phis));
    %Phis(:,Neu)=[];
    tmp = (circ_mean(Phis'));

    sx = fs{1};

    hold on;
    plot(tmp, 'k');

    plot(sp, tmp(sp), 'r');
    plot([xl(1) xl(end)], [-pi -pi], 'k');
    plot([xl(1) xl(end)], [0 0], 'k');
    plot([xl(1) xl(end)], [pi pi], 'k');
    axis([xl(1) xl(end) -2 * pi 2 * pi]);
    set(gca, 'xtick', xl, 'xticklabel', xlabel, 'ytick', [-pi 0 pi], 'yticklabel', {'-\pi' '0' '\pi'});
    set(gca, 'ydir', 'normal');

    %save test.mat Phis

end
