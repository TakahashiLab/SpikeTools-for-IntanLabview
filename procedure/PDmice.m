%set working directory
fn=pwd;
%
%set reference
%split hemisphere
forceRef{1,1}=1;forceRef{1,2}=1:4;forceRef{2,2}=[5:8];
%%
%load images....mat
Traj=FixPos(Pos);
save DLC.mat Traj;
%batch process for constructing neural and event signals
%1fn:working directory
%2: process option, 0:all
%3: reference assignment: forceRef
%4: verion: 2018 or new
%5: LED position: tetrode number
%6: realtime feedback?: 
%7  gpu to use?: 
%example for continuous stimulation
preProcessIlvrc4Mouse(fn,0,forceRef,0,3,0,0);
%example for realtime feedback
preProcessIlvrc4Mouse(fn,0,forceRef,0,3,1,0);
