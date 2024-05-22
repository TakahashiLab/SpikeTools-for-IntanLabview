function ResultTrutle2(results)
if ~isstruct(results)
    switch(results)
        case 0,
            fprintf('loading data directly\n');
            % Load the first set of variables
            load('O:\seaturtle\2023\combine\alpha\20230820\behav2.mat','results0');
            results=results0;
            % Load the second set of variables
            load('O:\seaturtle\2023\combine\gamma\20230821-1\behav2.mat','results0');
            results.individual=[results.individual;results0.individual];
            results.normal_field=[results.normal_field;results0.normal_field];
            results.zero_field=[results.zero_field;results0.zero_field];

        case 1,
            fprintf('loading data directly\n');
            % Load the first set of variables
            load('O:\seaturtle\2023\combine\alpha\20230820\behav2.mat','results0');
            results=results0;
        case 2,
            load('O:\seaturtle\2023\combine\gamma\20230821-1\behav2.mat','results0');
            results=results0;
    end
end

r=results;

rInd=r.individual;
size(rInd)
if size(rInd,2)==6
    conditions=2;
else
    conditions=4;
end

%mean resultant length
for i=1:conditions
    neuronId=find(rInd(:,i));

    subplot(4,4,i);
    polarhistogram(rInd(neuronId,i+conditions*2), 'BinEdges', linspace(-pi, pi, 8));

end

rNF=r.normal_field;
subplot(4,4,5);
neuronId=find(rNF(:,1));
polarhistogram(rNF(neuronId,3), 'BinEdges', linspace(-pi, pi, 8));


rZF=r.zero_field;
subplot(4,4,7);
neuronId=find(rZF(:,1));
polarhistogram(rZF(neuronId,3), 'BinEdges', linspace(-pi, pi, 8));

%mean vector length
for i=1:conditions
    neuronId=find(rInd(:,i+4));

    subplot(4,4,i+8);
    polarhistogram(rInd(neuronId,i+conditions*2), 'BinEdges', linspace(-pi, pi, 8));

end

rNF=r.normal_field;
subplot(4,4,5+8);
neuronId=find(rNF(:,2));
polarhistogram(rNF(neuronId,3), 'BinEdges', linspace(-pi, pi, 8));


rZF=r.zero_field;
subplot(4,4,7+8);
neuronId=find(rZF(:,2));
polarhistogram(rZF(neuronId,3), 'BinEdges', linspace(-pi, pi, 8));
return;