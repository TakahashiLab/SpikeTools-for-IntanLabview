function ensemble=makeEnsemble64(kkOuts,an,en)

ensemble=[];
tetNum=16;
for i=1:tetNum
    tmp=kkOuts{i}(an{i},:);
    
    if ~isempty(tmp)
        step=size(tmp{1,1},2)/size(tmp{1,3},2);
    end
    
    if ~isempty(en{i})
        en{i}
        rInd=[];
        for j=1:size(en{i},1)
            sameCell=en{i}{j};
            
            oind=find(an{i}==sameCell(1));
            n=tmp(oind,:);

            for k=2:length(sameCell)
                ind=find(an{i}==sameCell(k));
                rInd=[rInd ind];
                n{1,1}=[n{1,1} tmp{ind,1}];
                n{1,3}=[n{1,3} tmp{ind,3}];
                
                %reorder spike timings and waveform
                [n{1,3},ind]=sort(n{1,3});
                n{1,1}=getSpikesD(n{1,1},step,ind);
            end
            %            tmp=[tmp;n];
                tmp(oind,:)=n;
        end
        tmp(rInd,:)=[];
        end

        ensemble=[ensemble;tmp];
end

return;
