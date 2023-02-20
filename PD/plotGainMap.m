function outMph = plotGainMap(phaseHist, varargin)

    p = inputParser;
    p.addParamValue('display', 'both', @ischar);
    p.addParamValue('band', 'all', @ischar);
    p.addParamValue('hist', 'off', @ischar);
    p.addParamValue('topfreq', 120, @isnumeric);
    p.addParamValue('control', cell(1, 1), @iscell);
    p.addParamValue('xyaxis', 'off', @ischar);
    p.parse(varargin{:});

    disp = p.Results.display;
    band = p.Results.band;
    histon = p.Results.hist;
    topFreq = p.Results.topfreq;
    C = p.Results.control;
    phaseHistCtrl = C{1};
    xyaxis = p.Results.xyaxis;

    sd = 2; %5
    edges = [-3 * sd:1:3 * sd];
    kernel = normpdf(edges, 0, sd) * sd;
    kernel = kernel ./ sum(kernel(:));

    outMph = [];
    %imagesc(conv2(kernel,kernel,mean(phaseHist(:,[1:40 80:120],:),3)','same'));
    %mph=conv2(kernel,kernel,mean(phaseHist(:,[1:40 81:120],:),3)','same');

    switch band
        case 'theta',
            frange = 4:12;
        case 'beta',
            frange = 12:30;
        case 'gamma',
            frange = 31:120-2;
        case 'slowgamma',
            frange = 31:60;
        case 'fastgamma',
            frange = 60:80;
        case 'all',
            frange = 2:topFreq-2;
        otherwise,
            frange = [];
    end

    %topFreq=80;
    %topFreq = 120;
    topFreq = topFreq / 2;

    switch lower(disp)
        case 'div',

            mph = mean(phaseHist(:, [1:topFreq], :), 3)';
            mphC = mean(phaseHistCtrl(:, [1:topFreq], :), 3)';
            nC = nanmean(mphC, 2);
            mph = mph ./ nC;
            mph(nC == 0, :) = 0;

            mph = repeatMph(mph, 0, [3 3], 1);
            endMph = size(mph, 2);

            subplot(1, 2, 1);
            imagesc(mph);
            title('ascending');
            set(gca, 'ydir', 'normal');
            set(gca, 'ytick', [1:5:topFreq], 'yticklabel', [0:5:topFreq] * 2);
            set(gca, 'xtick', [1 10 20 30 40], 'xticklabel', {'-2\pi', '-\pi', '0', '\pi', '2\pi'});
            hold on;
            t = 0:endMph - 1;
            plot(t + 1, cos(t / pi) * topFreq / 2 + topFreq / 2, 'w');
            plot([1 endMph], [4/2 4/2], 'k');
            plot([1 endMph], [12/2 12/2], 'k');
            plot([1 endMph], [30/2 30/2], 'k');
            plot([1 endMph], [80/2 80/2], 'k');
            plot([endMph / 2 endMph / 2], [0 topFreq], 'k');

            colorbar;
            xlabel('Phase [radian]');
            ylabel('Frequency [Hz]');

            subplot(1, 2, 2);
            mph = mean(phaseHist(:, [120:-1:120 - topFreq + 1], :), 3)';
            mphC = mean(phaseHistCtrl(:, [120:-1:120 - topFreq + 1], :), 3)';

            nC = nanmean(mphC, 2);
            mph = mph ./ nC;
            mph(nC == 0, :) = 0;
            mph = repeatMph(mph, 0, [3 3], 1);

            imagesc(mph);
            title('descending');
            set(gca, 'ydir', 'normal');
            %set(gca,'ytick',[1:5:40],'yticklabel',[(40:-5:1)-4]*2);
            set(gca, 'ytick', [1:5:topFreq], 'yticklabel', [0:5:topFreq] * 2);
            set(gca, 'xtick', [1 10 20 30 40], 'xticklabel', {'-2\pi', '-\pi', '0', '\pi', '2\pi'});
            hold on;

            xlabel('Phase [radian]');
            ylabel('Frequency [Hz]');
            plot(t + 1, cos(t / pi) * topFreq / 2 + topFreq / 2, 'w');
            plot([1 endMph], [4/2 4/2], 'k');
            plot([1 endMph], [12/2 12/2], 'k');
            plot([1 endMph], [30/2 30/2], 'k');
            plot([1 endMph], [80/2 80/2], 'k');
            plot([endMph / 2 endMph / 2], [0 topFreq], 'k');

            colormap(jet);
            colorbar;
            caxis([0.5 max(mph(:))]);

        case 'ascend',
            xr = [1:topFreq];
            xrange = [0:5:topFreq];

            outMph = plotfnc(phaseHist, phaseHistCtrl, xr, xrange, kernel, frange, band, xyaxis, histon, topFreq);

        case 'descend',
           
            xr = [1:topFreq];
            phaseHist=phaseHist(:,120:-1:61,:);
            phaseHistCtrl=phaseHistCtrl(:,120:-1:61,:);
            xrange = [0:5:topFreq];

            outMph = plotfnc(phaseHist, phaseHistCtrl, xr, xrange, kernel, frange, band, xyaxis, histon, topFreq);
            %%title('descending');

        case 'both',
        
            xr = [1:topFreq];
            xrange = [0:5:topFreq];
      
            phaseHist=(phaseHist(:,1:60,:)+phaseHist(:,120:-1:61,:))/2;
            phaseHistCtrl=(phaseHistCtrl(:,1:60,:)+phaseHistCtrl(:,120:-1:61,:))/2;

            outMph = plotfnc(phaseHist, phaseHistCtrl, xr, xrange, kernel, frange, band, xyaxis, histon, topFreq);
    end

    return;
    %%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
function mphV = plotfnc(phaseHist, phaseHistCtrl, xr, xrange, kernel, frange, band, xyaxis, histon, topFreq)
    binsize = 1;

    outMph = [];
    
    %set frequency band
    frange=[frange frange(end)+2];%offset 
    frange = unique(round(frange / 2));
   
    xr = xr(frange);
    
    %normalization
  
    nC=nanmean(squeeze(nanmean(phaseHistCtrl,1)),1);
   
    nC(nC==0)=NaN;
  
    nCs=zeros(size(phaseHistCtrl));
    nC=repmat(nC,size(nCs,2),1);
   
    for q=1:size(nCs,1)
        nCs(q,:,:)=nC;
    end 
  
    phaseHistN=phaseHist./nCs;
    phaseHistCtrlN=phaseHistCtrl./nCs;
   

     if strcmp(xyaxis, 'neuronphase') %neuron x phase
        xr=xr(1:end-1);
        mph = squeeze(nanmean(phaseHist(:, xr, :), 2))';
        mphC = squeeze(nanmean(phaseHistCtrlN(:, xr, :), 2))';
        mphN = squeeze(nanmean(phaseHistN(:, xr, :), 2))';
        mphStat = mphN;%statistical test
     
        mphV=mphN;%for bar  
        
        [maxMph, maxMphId] = max(mphV, [], 2);
        [~, id] = sort(maxMph, 'descend');
        %mph = mphV(id, :);
        
        mph = repeatMph(mph,0);
        mph = log10(mph);
        imagesc(mph);
        set(gca, 'ytick', [1 size(mph, 1)], 'yticklabel', [1 size(mph, 1)]);
        clim([-1 2]);
        ylabel('Neuron number');

    elseif strcmp(xyaxis, 'freqphase') %freq x phase
    
        mph = nanmean(phaseHist(:, xr, :), 3)';
        mphC = nanmean(phaseHistCtrlN(:, xr, :), 3)';
        mphN = nanmean(phaseHistN(:, xr, :), 3)';
  

        mphStat = mphN;%statistical test
        
        mphV=mphN;%for bar   
        
        mph = repeatMph(mph, 1, [3 3], 1);
       %mph=repeatMph(mph,0);
        
        imagesc(mph);
    
        lenMph=size(mph,1);
        if lenMph > 20
            ytick = 1:floor(lenMph / 10):lenMph;
        else
            ytick = 1:lenMph;
        end
       
        yticklabel = xr(ytick) * 2;
        set(gca, 'ytick', ytick, 'yticklabel', yticklabel);
        ylabel('Frequency [Hz]');
    elseif strcmp(xyaxis, 'gainphase')
        xr=xr(1:end-1);
        mph = squeeze(nanmean(phaseHist(:, xr, :), 2))';
        mphC = squeeze(nanmean(phaseHistCtrlN(:, xr, :), 2))';
        mphN = squeeze(nanmean(phaseHistN(:, xr, :), 2))';
        mphStat = mphN;%statistical test
       
        mphV=mphN;%for bar  
        outMph = mph;
        
        seMph = (std(mph,0,1) ./ sqrt(size(mph,1))); %sem
        mph = nanmean(mph, 1);
        
        mph = repeatMph(mph,0);
        seMph = repeatMph(seMph,0);
        errP=mph+seMph;
        errN=mph-seMph;
        mph=log10(mph);
        errP=log10(errP);
        errN=log10(errN);
        
        plot(mph);

        hold on; plot(errP, 'r'); plot(errN, 'r');
        endMph = size(mph, 2);
        t = 1:endMph;
        set(gca, 'xtick', [1 endMph / 4 endMph / 4 * 2 endMph / 4 * 3 endMph], 'xticklabel', {'-2\pi', '-\pi', '0', '\pi', '2\pi'});

        xlabel('phase');
        ylabel('Gain');

    elseif strcmp(xyaxis, 'gainfreq')
        xr=xr(1:end-1);
        mph = squeeze(mean(phaseHist(:, xr, :), 1));
        mphC = squeeze(mean(phaseHistCtrlN(:, xr, :), 1));
        mphN = squeeze(mean(phaseHistN(:, xr, :), 1));
        mphStat = mphN;%statistical test
        %nC = nanmean(mphC, 1);%1:neuron,2:freq
        %mph = mph ./ nC;
        mphV=mphN;%for bar  
        outMph = mph;
        mph=mph(2:end-1,:);

        mph = mean(mph, 2)';
        mph=log10(mph);

        plot(mph);

        set(gca, 'xtick', xrange, 'xticklabel', xrange * 2 * binsize);
        xlabel('Frequency [Hz]');
        ylabel('Gain');

    end

    if (strcmp(xyaxis, 'gainfreq') | strcmp(xyaxis, 'gainphase'))

    elseif strcmp(histon, 'off')

        %%title('ascending');
        set(gca, 'ydir', 'normal');
        mphY = size(mph, 1);
        endMph = size(mph, 2);
        t = 1:endMph;
        set(gca, 'xtick', [1 endMph / 4 endMph / 4 * 2 endMph / 4 * 3 endMph], 'xticklabel', {'-2\pi', '-\pi', '0', '\pi', '2\pi'});
        hold on;
        plot(t, cos(t / max(t) * 4 * pi) * mphY / 2 + mphY / 2, 'w');
        xlabel('Phase [ radian]');
        colormap(jet);
        colorbar;

        if ~isempty(mph)

            if ~isnan(mph(1))
                %  caxis([0.5 max(mph(:))]);
            end

        end

 else
       
            alpha=0.05;
                [x,y]=find(isnan(mphStat));
            mphStat(x,:)=[];
            %[clusters,p]=permutest(mphStat',mphC',false,alpha,1000,false);
            [clusters,p]=permutest(mphStat',ones(size(mphStat))',false,alpha,1000,false);
            clusters
            p


            sp=cell2mat(clusters(p<alpha));

            data=mean(mphStat);
            pd=makedist('uniform','lower',0,'upper',max(data));
                       
            [h,p]=kstest(data,'cdf',[data', cdf(pd,data')]);
            fprintf('ks test=%f\n',p);
            [h,p]=chi2gof(data,'Expected',ones(size(data)));
            fprintf('chi2gof=%f\n',p);
       
        
        zShift = size(mphV, 2)/2;
    
        mph = repeatMph(mphV,0);
        mph = nanmean(mph);
        
       bar(mph);
        
        hold on;
        for z = sp
            plot([z + zShift], mph([z + zShift]), '*');
        end

        endMph = size(mph, 2);
        t = 1:endMph;
        set(gca, 'xtick', [1 endMph / 4 endMph / 4 * 2 endMph / 4 * 3 endMph], 'xticklabel', {'-2\pi', '-\pi', '0', '\pi', '2\pi'});

    end

    
end

%%%%%%%%%%%%%%
function mph = repeatMph(mph, offset, grid, sd)

    hf = size(mph, 2) / 2;
    mph = [mph mph mph];
    
    if nargin > 2
        mph = SmoothMat(mph, grid, sd);
    end
    
    mph = mph(1:end-offset, hf + 1:end -hf);%cut offset

end
