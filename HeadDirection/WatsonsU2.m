%watsonsu2 = WatsonsU2(spk_headdir, headdir);
function watsonsu2 = WatsonsU2(x, y)
        % watsonsu2 = root.HDWatsonsU2(cel)
        %
        % Returns U2 score for vectors of directional data x and y
        %
        % watsonsu2: Watson's U^2 test for uniformity in circular data. See
        % redish,
        % hippocampus 2005 to learn how it applies to head direction cells
        %
        % see
        % http://www.ncbi.nlm.nih.gov/pmc/articles/mid/NIHMS71281/table/T2/
        % for algorithm
        %
        % also see merriam and davis 2008 pg 10
        %
        % a common way to apply this to head direction/spiking data is to
        % pass the array of head direction samples and spike head direction
        % angles. the null hypothesis is thus that spikes head directions come 
        % from the same distribution of all head direction angles. the null 
        % hypothesis is rejected if U^2 is 'large'. >10 seems to be the 
        % neuroscience HD cell literature standard. - andrew

        % ehren newman march 2010
        
% remove any nan's from x and y
              x(isnan(x)) = [];
              y(isnan(y)) = [];

              if isempty(x) || isempty(y), watsonsu2 = NaN; return; end
              

              
              z = cat(1, x, y);
              [z,I] = sort(z,'ascend');

              N = length(z);
              n1= length(x);
              n2= length(y);

              z2 = cat(1, ones(n1,1)/n1, -1*ones(n2,1)/n2);
              z2 = z2(I); % sort in same order as z

              dk = cumsum(z2);

              correctIt = 0;
              if correctIt
                % correct for only counting each repeat as coming from one dist
                if n1<n2
                  overlap = sum(ismember(x,y));
                else
                  overlap = sum(ismember(y,x));
                end      
                correction = (1/n2) * ones(overlap,1);

                dbar = (sum(dk) - sum(correction))/N;    
                watsonsu2 = ((n1*n2) / N^2) * (sum( (dk-dbar) .^ 2 ) - sum(correction.^2));

              else

                dbar = sum(dk./N);
                watsonsu2 = ((n1*n2) / N^2) * sum( (dk-dbar) .^ 2 );

              end    
        end            
        
