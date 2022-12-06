%%%%%%%%%%%%%%
function [resPost,resDuring,Data,Lists,S]=processData2(AA,D,S,verbose)

c={'r.','g.','b.','k.'};

Data=AA;

%beta
betaList=find(sum(str2mat(D{:,2})=='beta ',2)==5);
%theta
thetaList=find(sum(str2mat(D{:,2})=='theta',2)==5);
%gamma
gammaList=find(sum(str2mat(D{:,2})=='gamma',2)==5);

Lists{1}=betaList;
Lists{2}=thetaList;
Lists{3}=gammaList;

resPost=[];
resDuring=[];

if verbose
    fprintf('###################################\n');
    fprintf('%%%%%%%%%%%%%%%beta stimuli\n');
    fprintf('###################################\n');
end
[resP,resD]=plotPL(betaList,Data,c{1},D,S,verbose);
resPost=[resPost;resP];
resDuring=[resDuring;resD];


%save test.mat betaList Data
%return;
if verbose
    fprintf('###################################\n');
    fprintf('%%%%%%%%%%%%%%theta stimuli\n');
    fprintf('###################################\n');
end
[resP,resD]=plotPL(thetaList,Data,c{2},D,S,verbose);
resPost=[resPost;resP];
resDuring=[resDuring;resD];

if verbose
    fprintf('###################################\n');
    fprintf('%%%%%%%%%%%%%%gamma stimuli\n');
    fprintf('###################################\n');
end

[resP,resD]=plotPL(gammaList,Data,c{3},D,S,verbose);
resPost=[resPost;resP];
resDuring=[resDuring;resD];

return;