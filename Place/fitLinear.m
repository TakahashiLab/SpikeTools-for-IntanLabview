function [LTraj,LTrajS,GI]=fitLinear(Traj)
h=figure;
plot(Traj(:,3),Traj(:,4));
fprintf(['Check the four corners of the rectangule with the mouse\n' ...
         'in a clockwise direction\n']);
fprintf('Please start at the upper left corner');
[x,y]=ginput(8);

x1=[x(1):x(2)];
y1=[y(1):((y(2)-y(1))./(length(x1)-1)):y(2)];

x2=[x(2):x(3)];
%y2=[ones(1,length(x2))*y(2)];
y2=[y(2):((y(3)-y(2))./(length(x2)-1)):y(3)];

x3=[x(3):x(4)];
y3=[y(3):((y(4)-y(3))./(length(x3)-1)):y(4)];

y4=[y(4):-1:y(5)];
%x4=[ones(1,length(y4))*x(4)];
x4=[x(4):((x(5)-x(4))./(length(y4)-1)):x(5)];

x5=[x(4):-1:x(6)];
y5=[y(5):((y(6)-y(5))./(length(x5)-1)):y(6)];

x6=[x(6):-1:x(7)];
%y6=[ones(1,length(x6))*y(6)];
y6=[y(6):((y(7)-y(6))./(length(x6)-1)):y(7)];

x7=[x(7):-1:x(8)];
y7=[y(7):((y(8)-y(7))./(length(x7)-1)):y(8)];

y8=[y(8):y(1)];
x8=[x(8):((x(1)-x(8))./(length(y8)-1)):x(1)];
%x8=[ones(1,length(y8))*x(1)];

GI(1,:)=[x1 x2 x3 x4 x5 x6 x7 x8];
GI(2,:)=[y1 y2 y3 y4 y5 y6 y7 y8];


hold on;
plot(GI(1,:),GI(2,:),'r');

LTraj=dsearchn(GI',Traj(:,3:4));%
LTrajS=dsearchn(GI',Traj(:,5:6));%

pause(2);
close(h);
return;