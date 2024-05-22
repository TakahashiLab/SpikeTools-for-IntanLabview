function [cc,sb,Sn,HD,Sn2,HD2,FR,FR2]=batchSeaturtle(Out,msT,Pos,Sts)

n=size(Out,1);

if size(Sts,2)==4
    s2{1}=union(Sts{1},Sts{3});
    s2{2}=union(Sts{2},Sts{4});
else
    s2=Sts;
end

cc=[];
sb=[];
Sn=[];
HD=[];
Sn2=[];
HD2=[];
FR=[];
FR2=[];
for j=1:n
    rm=[];
    %correlation between states
    Sig=[];
    Hd=[];
    for i=1:2
        %[sig95,sig99]=makeSigShuffle(Out{j,3},Pos(s2{i},:),msT(s2{i}),'animal','seaturtle','shufflen',100,'shuffletype','shift');
        [ratemap,hd,mr,sig]=plot_polar_rate_map(Out{j,3},Pos(s2{i},:),msT(s2{i}),'animal','seaturtle');
      
        rm=[rm; ratemap];
        Sig=[Sig sig];
        %Sig=[Sig sig95];
        Hd=[Hd hd];
    end
    fr=max(rm');
 
    Sn=[Sn; Sig];
    cc=[cc corr(rm(1,:)',rm(2,:)')];
    HD=[HD; Hd];
    FR=[FR;fr];

    if size(Sts,2)>=4
        %stability
        rm=[];
        Sig=[];
        Hd=[];
        for i=1:size(Sts,2)
             %[sig95,sig99]=makeSigShuffle(Out{j,3},Pos(Sts{i},:),msT(Sts{i}),'animal','seaturtle','shufflen',100,'shuffletype','shift');
            [ratemap,hd,mr,sig]=plot_polar_rate_map(Out{j,3},Pos(Sts{i},:),msT(Sts{i}),'animal','seaturtle');
            rm=[rm; ratemap];
            Sig=[Sig sig];
            %Sig=[Sig sig95];
            Hd=[Hd hd];
        end

        fr=max(rm');
        Sn2=[Sn2; Sig];
        HD2=[HD2;Hd];
        FR2=[FR2;fr];
        sb=[sb; corr(rm(1,:)',rm(3,:)') corr(rm(2,:)',rm(4,:)') corr(rm(1,:)',rm(2,:)') corr(rm(3,:)',rm(4,:)')];
    end

end
return;
