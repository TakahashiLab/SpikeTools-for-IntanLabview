function ResultTrutle2(results,varargin)
p = inputParser;
p.addParamValue('opt', 'polar', @ischar);
p.addParamValue('sense', 'uni', @ischar);
p.addParamValue('animal', 'seaturtle', @ischar);
p.addParamValue('folded', 'unfolded', @ischar);
p.addParamValue('score', 'MRL', @ischar);

p.parse(varargin{:});
opt = p.Results.opt;
bothSense = p.Results.sense;
animal=p.Results.animal;
folded=p.Results.folded;
scoreType=p.Results.score;

declination=0;

switch lower(animal)
    case 'seaturtle',
        conv=1;
        declination=0;
    case 'seabird',
        switch (folded)
            case 'folded',
                conv=1;
            case 'unfolded'
                conv=0;%convolution applied in the batchSeaturtleGPT process
        end
        declination=-8;%
    otherwise,
        fprintf('animal type error\n');
end

declination=deg2rad(declination);

if nargin==1
    opt=1;%polar plot
    bothSense=0;%uidiretional 1: bi-directional, 2: both uni- and bi-directional
end

switch lower(opt)
    case 'polar',
        fprintf('polar plot\n');
    case 'rate',
        fprintf('normalized rate map\n');
end

switch lower(bothSense)
    case 'uni',
        fprintf('unidirectional sensitivity\n');
    case 'bi',
        fprintf('bidirectional sensitivity\n');
    case 'both',
        fprintf('both uni- and bi-directional sensitivity\n');
end

FrTh=1;%1Hz
fprintf('Firing rate threshold=%dHz\n',FrTh);

if ~isstruct(results)
    switch (animal)
        case 'seaturtle',
            switch (folded)
                case 'folded',
                    matfile='behav4.mat';
                case 'unfolded'
                    matfile='behav5.mat';
            end

            switch(results)
                %seaturtle
                case 0,
                    fprintf('cortex (alpha+beta+gamma)\n');

                    % Load the first set of variables
                    load(['O:\seaturtle\2023\combine\alpha\20230821\' matfile],'results0');
                    results=results0;
                    % Load the second set of variables
                    if 1
                        load(['O:\seaturtle\2023\combine\gamma\20230821-1\' matfile],'results0');
                        results.individual=[results.individual;results0.individual];
                        results.normal_field=[results.normal_field;results0.normal_field];
                        results.zero_field=[results.zero_field;results0.zero_field];
                        results.peakfr=[results.peakfr;results0.peakfr];
                        results.ratemap=[results.ratemap;results0.ratemap];
                    end
                    % Load the second set of variables
                    load(['O:\seaturtle\2023\combine\beta\20230820\' matfile],'results0');
                    results.individual=[results.individual;results0.individual];
                    results.normal_field=[results.normal_field;results0.normal_field];
                    results.zero_field=[results.zero_field;results0.zero_field];
                    results.peakfr=[results.peakfr;results0.peakfr];
                    results.ratemap=[results.ratemap;results0.ratemap];

                case 1,%Cortex (alpha 0821)
                    fprintf('cortex (alpha 0821)\n');
                    % Load the first set of variables
                    load(['O:\seaturtle\2023\combine\alpha\20230821\' matfile],'results0');
                    results=results0;
                case 2,%Cortex (gamma 0821)
                    fprintf('cortex (gamma 0821)\n');
                    load(['O:\seaturtle\2023\combine\gamma\20230821-1\' matfile],'results0');
                    results=results0;
                case 3,%normal-normal (alpha normal-normal)
                    fprintf('cortex (alpha 0819 normal vs normal)\n');
                    load(['O:\seaturtle\2023\combine\alpha\20230819\' matfile],'results0');
                    results=results0;
                case 4,%Hippocampus DVR
                    fprintf('DVR (alpha 0820)\n');
                    load(['O:\seaturtle\2023\combine\alpha\20230820\' matfile],'results0');
                    results=results0;
                case 5,%Cortex (beta 0820)
                    fprintf('cortex (beta 0820)\n');
                    load(['O:\seaturtle\2023\combine\beta\20230820\' matfile],'results0');
                    results=results0;
            end


        case 'seabird',
            %%%%%%%%%%%%%%%%%%%%%%%
            %seabird
            %%%%%%%%%%%%%%%%%%%
            switch (folded)
                case 'folded',
                    matfile='behav3.mat';
                case 'unfolded'
                    matfile='behav4.mat';
            end
      
            switch(results)
                %
                case 0,%indonesia, HP
                    fprintf('seabird indonesia HP (D and F)\n');
                    % Load the first set of variables
                    load(['P:\seabird\2022\D_SF_indonesia.kilosort\' matfile],'results0');
                    results=results0;
                    % Load the second set of variables
                    load(['P:\seabird\2022\F_DP_indonesia.kilosort\' matfile],'results0');

                    results.individual=[results.individual;results0.individual];
                    results.normal_field=[results.normal_field;results0.normal_field];
                    results.zero_field=[results.zero_field;results0.zero_field];
                    results.peakfr=[results.peakfr;results0.peakfr];
                    results.ratemap=[results.ratemap;results0.ratemap];
                case 1,%counterbalance, HP
                    fprintf('seabird counterbalance HP(D and F)\n');
                    % Load the first set of variables
                    load(['P:\seabird\2022\D_SF_counterbalance.kilosort\'  matfile],'results0');
                    results=results0;
                    % Load the second set of variables
                    load(['P:\seabird\2022\F_DP_counterbalance.kilosort\'  matfile],'results0');
                    results.individual=[results.individual;results0.individual];
                    results.normal_field=[results.normal_field;results0.normal_field];
                    results.zero_field=[results.zero_field;results0.zero_field];
                    results.peakfr=[results.peakfr;results0.peakfr];
                    results.ratemap=[results.ratemap;results0.ratemap];
                case 2,%indonesia PF
                    fprintf('seabird indonesia PF nidopallium (F)\n');
                    load(['P:\seabird\2022\F_SF_indonesia.kilosort\'  matfile],'results0');
                    results=results0;
                case 3,%counterbalance PF
                    fprintf('seabird counterbalance PF nidopallium ( F)\n');
                    load(['P:\seabird\2022\F_SF_counterbalance.kilosort\'  matfile],'results0');
                    results=results0;
                case 4,%indonesia, HP
                    fprintf('seabird indonesia HP (D)\n');
                    % Load the first set of variables
                    load(['P:\seabird\2022\D_SF_indonesia.kilosort\' matfile],'results0');
                    results=results0;
                case 5,
                     fprintf('seabird indonesia HP (F)\n');
                    % Load the second set of variables
                    load(['P:\seabird\2022\F_DP_indonesia.kilosort\' matfile],'results0');
                    results=results0;
                 
                case 6,%counterbalance, HP
                    fprintf('seabird counterbalance HP(D)\n');
                    % Load the first set of variables
                    load(['P:\seabird\2022\D_SF_counterbalance.kilosort\'  matfile],'results0');
                    results=results0;

                case 7,
                    fprintf('seabird counterbalance HP(F)\n');
                    % Load the second set of variables
                    load(['P:\seabird\2022\F_DP_counterbalance.kilosort\'  matfile],'results0');
                     results=results0;
            end
    end
end

r=results;

rInd=r.individual;
rPeak=r.peakfr;
rRM=r.ratemap;

%magnetic declination correction 
rInd(:,5:6)=rInd(:,5:6)+declination;
fprintf('%d neurons\n',size(rInd,1));

if size(rInd,2)==10
    conditions=2;
else
    conditions=4;
end

switch lower(scoreType)
    case 'mrl',
        st=3;
        fprintf('mean resultantn length\n');
    case 'mvl',%mean vector length, Rayleigh vector score 
        st=4;
        fprintf('mean vector length\n');
end

NeuronId=[];
%mean resultant length
hdcells=[];
for i=1:conditions
    
    switch lower(bothSense)
        case 'uni',
            neuronId=find(rInd(:,i+conditions*(st-3)));
        case 'bi',
            neuronId=find(rInd(:,i+conditions*st));
        case 'both',
            neuronId=find(rInd(:,i+conditions*(st-3)));
            neuronId=union(neuronId,find(rInd(:,i+conditions*st)));
    end
    neuronId=neuronId(rPeak(neuronId,i) > FrTh)';

    NeuronId=union(NeuronId,neuronId);
    if ~isempty(neuronId)
        p=circ_rtest(rInd(neuronId,i+conditions*2));
        fprintf('p=%f\n',p);
        p=circ_rtest(rInd(neuronId,i+conditions*2)*2);
        fprintf('double p=%f\n',p);
    end
   
    subplot(4,4,i);
    switch lower(opt)
        case 'polar',
            polarhistogram(rInd(neuronId,i+conditions*2), 'BinEdges', linspace(-pi, pi, 8));
            thetaticks([0 90 180 270]);
            thetaticklabels({'E','N','W','S'});
        case 'polarscat',
            polarscatter(rInd(neuronId,i+conditions*2),ones(size(neuronId)),'.');
        case 'rate',
            dispRM(rRM,i,neuronId,conv);
    end
    title(sprintf('condition #%d,%d',i,length(neuronId)));
    hdcells=union(hdcells,neuronId);

    subplot(4,4,i+8)
    neuronId=1:size(rInd,1);
   
    switch lower(opt)
        case 'polar',
            polarhistogram(rInd(neuronId,i+conditions*2), 'BinEdges', linspace(-pi, pi, 8));
            thetaticks([0 90 180 270]);
            thetaticklabels({'E','N','W','S'});
        case 'polarscat',
            polarscatter(rInd(neuronId,i+conditions*2),ones(size(neuronId)),'.');
        case 'rate',
            dispRM(rRM,i,neuronId,conv);
            
    end
    title(sprintf('condition #%d, all neurons',i));
end
fprintf('%d/%d neurons, %f%%\n',length(hdcells),size(rInd,1),length(hdcells)/size(rInd,1)*100);
neuronId=NeuronId';

%identified in either condition
for i=1:conditions
    subplot(4,4,i+4);

    switch lower(opt)
        case 'polar',
            polarhistogram(rInd(neuronId,i+conditions*2), 'BinEdges', linspace(-pi, pi, 8));
            thetaticks([0 90 180 270]);
            thetaticklabels({'E','N','W','S'});
        case 'polarscat',
            polarscatter(rInd(neuronId,i+conditions*2),ones(size(neuronId)),'.');
        case 'rate',
            dispRM(rRM,i,neuronId,conv);
            
    end
    title(sprintf('either condition #%d,%d',i,length(neuronId)));
end

%firing rate comparison
[mRM1,mRM2,mxRM1,mxRM2]=frComp(rRM,conditions);

subplot(4,4,16)
plot([ones(size(mRM1)) ones(size(mRM2))*2]',[mRM1 mRM2]');
[h,p]=ttest2(mRM1',mRM2');
fprintf('mean firing rate difference: p=%f\n',p);

[h,p]=ttest2(mxRM1',mxRM2');
fprintf('max firing rate difference: p=%f\n',p);


%firing rate comparison
[mRM1,mRM2,mxRM1,mxRM2]=frComp(rRM(neuronId,:,:),conditions);

subplot(4,4,16)
plot([ones(size(mRM1)) ones(size(mRM2))*2]',[mRM1 mRM2]');
[h,p]=ttest2(mRM1',mRM2');
fprintf('mean firing rate difference(selected): p=%f\n',p);

[h,p]=ttest2(mxRM1',mxRM2');
fprintf('max firing rate difference(selected): p=%f\n',p);

%circular distance between normal and manipulation  
subplot(4,4,15);
polarhistogram(circ_dist(rInd(neuronId,conditions*2+1),rInd(neuronId,conditions*2+2)), 'BinEdges', linspace(-pi, pi, 8));
title('diff');

%%%%%manipulation
rNF=r.normal_field;
rZF=r.zero_field;
%magnetic declination correction
rNF(:,3)=rNF(:,3)+declination;
rZF(:,3)=rZF(:,3)+declination;

if size(rNF,2)>3
    subplot(4,4,13);

    switch lower(bothSense)
        case 'uni',
            neuronId=find(rNF(:,1+(st-3)) & rNF(:,4)>FrTh);
        case 'bi',
            neuronId=find(rNF(:,5+(st-3)) & rNF(:,4)>FrTh);
        case 'both',
            neuronId=find(rNF(:,1+(st-3)) & rNF(:,4)>FrTh);
            neuronId=union(neuronId,find(rNF(:,5+(st-3)) & rNF(:,4)>FrTh));
    end
    
    switch lower(opt)
        case 'polar',
            polarhistogram(rNF(neuronId,3), 'BinEdges', linspace(-pi, pi, 8));
            thetaticks([0 90 180 270]);
            thetaticklabels({'E','N','W','S'});
        case 'polarscat',
            polarscatter(rNF(neuronId,3),ones(size(neuronId)),'.');
        case 'rate',
            i=5;
            dispRM(rRM,i,neuronId,conv);
        
    end
    title('normal');

    fprintf('normal all\n');
    if ~isempty(neuronId)
        p=circ_rtest(rNF(neuronId,3));
        fprintf('p=%f\n',p);
        p=circ_rtest(rNF(neuronId,3)*2);
        fprintf('double p=%f\n',p);
    end

  
    subplot(4,4,14);
    switch lower(bothSense)
        case 'uni',
            neuronId=find(rZF(:,1+(st-3)) & rZF(:,4)>FrTh);
        case 'bi',
            neuronId=find(rZF(:,5+(st-3)) & rZF(:,4)>FrTh);
        case 'both',
            neuronId=find(rZF(:,1+(st-3)) & rZF(:,4)>FrTh);
            neuronId=union(neuronId,find(rZF(:,5+(st-3)) & rZF(:,4)>FrTh));
    end
   
     switch lower(opt)
        case 'polar',
             polarhistogram(rZF(neuronId,3), 'BinEdges', linspace(-pi, pi, 8));
             thetaticks([0 90 180 270]);
             thetaticklabels({'E','N','W','S'});
        case 'polarscat',
            polarscatter(rZF(neuronId,3),ones(size(neuronId)),'.');
         case 'rate',
            i=6;
            dispRM(rRM,i,neuronId,conv);
    end
   
    title('zero / indo');

    fprintf('zero /indo all\n');
    if ~isempty(neuronId)
        p=circ_rtest(rZF(neuronId,3));
        fprintf('p=%f\n',p);
        p=circ_rtest(rZF(neuronId,3)*2);
        fprintf('double p=%f\n',p);
    end
end



return;
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mRM1,mRM2,mxRM1,mxRM2]=frComp(rRM,conditions)
%firing rate comparison
switch (conditions)
    case 2,
        rRM1=squeeze(rRM(:,1,:));
        rRM2=squeeze(rRM(:,2,:));

        mRM1=mean(rRM1,2);
        mRM2=mean(rRM2,2);
        mxRM1=max(rRM1,[],2);
        mxRM2=max(rRM2,[],2);
    case 4,
        rRM1=squeeze(sum(rRM(:,[1 3],:),2));
        rRM2=squeeze(sum(rRM(:,[2 4],:),2));
             
        mRM1=mean(rRM1,2);
        mRM2=mean(rRM2,2);
        mxRM1=max(rRM1,[],2);
        mxRM2=max(rRM2,[],2);
end
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispRM(rRM,num,neuronId,conv)
rRMd=squeeze(rRM(:,num,:));

if conv
    for k=1:size(rRMd,1)
        rRMd(k,:)=kSmooth(rRMd(k,:));
    end
end
rRMd=rRMd./max(rRMd,[],2);

dispRM=rRMd(neuronId,:);
[~,id]=max(dispRM,[],2);
[~,id]=sort(id);
dispRM=dispRM(id,:);
imagesc(dispRM);
xticks([1 15 30 45 60]);
xticklabels({'west','south','east','north','west'});
colormap("jet")
return;

function hoge()
%%%%%obsolete
NeuronId=[];
%mean vector length
hdcells=[];
for i=1:conditions
    neuronId=find(rInd(:,i+conditions));
   
    switch lower(bothSense)
        case 'uni',
            %nothiing
        case 'bi',
            neuronId=find(rInd(:,i+conditions*4));
        case 'both',
            neuronId=union(neuronId,find(rInd(:,i+conditions*4)));
    end
    neuronId=neuronId(rPeak(neuronId,i) > FrTh)';
    NeuronId=union(NeuronId,neuronId);
    if ~isempty(neuronId)
        p=circ_rtest(rInd(neuronId,i+conditions*2));
        fprintf('p=%f\n',p);
        subplot(4,4,i+8);
        switch lower(opt)
            case 'polar',
                polarhistogram(rInd(neuronId,i+conditions*2), 'BinEdges', linspace(-pi, pi, 8));
            case 'rate',
                rRMd=squeeze(rRM(:,i,:));
                if conv
                    for k=1:size(rRMd,1)
                        rRMd(k,:)=kSmooth(rRMd(k,:));
                    end
                end

                rRMd=rRMd./max(rRMd,[],2);
                dispRM=rRMd(neuronId,:);
                [~,id]=max(dispRM,[],2);
                [~,id]=sort(id);
                dispRM=dispRM(id,:);
                imagesc(dispRM);
                xticks([1 15 30 45 60]);
                xticklabels({'west','south','east','north','west'});
                colormap("jet")
        end
        title(sprintf('condition #%d,%d',i,length(neuronId)));
        hdcells=union(hdcells,neuronId);
    end
end

fprintf('%d/%d neurons, %f%%\n',length(hdcells),size(rInd,1),length(hdcells)/size(rInd,1)*100);
neuronId=NeuronId';

%identified in either condition
for i=1:conditions
    subplot(4,4,i+8+2);
    
    switch lower(opt)
        case 'polar',
            polarhistogram(rInd(neuronId,i+conditions*2), 'BinEdges', linspace(-pi, pi, 8));
        case 'rate',
            rRMd=squeeze(rRM(:,i,:));
            if conv
                for k=1:size(rRMd,1)
                    rRMd(k,:)=kSmooth(rRMd(k,:));
                end
            end
            rRMd=rRMd./max(rRMd,[],2);
            dispRM=rRMd(neuronId,:);
            [~,id]=max(dispRM,[],2);
            [~,id]=sort(id);
            dispRM=dispRM(id,:);
            imagesc(dispRM);
            xticks([1 15 30 45 60]);
            xticklabels({'west','south','east','north','west'});
            colormap("jet")
    end
    title(sprintf('condition #%d,%d',i,length(neuronId)));
end

%neuronId=find(any(rInd(:,1+conditions:conditions*2)'))
subplot(4,4,14);
polarhistogram(circ_dist(rInd(neuronId,conditions*2+1),rInd(neuronId,conditions*2+2)), 'BinEdges', linspace(-pi, pi, 8));
title('diff');

rNF=r.normal_field;
if size(rNF,2)>3
    subplot(4,4,5+8);
    neuronId=find(rNF(:,2) & rNF(:,4)>FrTh)';;
    polarhistogram(rNF(neuronId,3), 'BinEdges', linspace(-pi, pi, 8));
    title('normal')

    rZF=r.zero_field;
    subplot(4,4,7+8);
    neuronId=find(rZF(:,2) & rZF(:,4)>FrTh)';
    polarhistogram(rZF(neuronId,3), 'BinEdges', linspace(-pi, pi, 8));
    title('zero / indo')
end
return;
%%%%
function tc=kSmooth(tc)
gaussianSD = 3; %
tcLen=length(tc);
tc=[tc tc tc];
kernelSize = ceil(gaussianSD * 6); % 3 standard deviations on either side
gaussianKernel = fspecial('gaussian', [1, kernelSize], gaussianSD);
tc=conv2(tc, gaussianKernel, 'same');
tc=tc(tcLen+1:tcLen*2);
return;