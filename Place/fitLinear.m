function [LTraj,LTrajS,TMpoints,GI]=fitLinear(Traj)
h=figure;
plot(Traj(:,3),Traj(:,4));
fprintf(['Check the 10 corners of the rectangule with the mouse\n' ...
         'in a clockwise direction\n']);
fprintf('Please start at the left reward zone\n');
[x,y]=ginput(10);

x1=[x(1):x(2)];
y1=[y(1):((y(2)-y(1))./(length(x1)-1)):y(2)];

x2=[x(2):x(3)];
y2=[y(2):((y(3)-y(2))./(length(x2)-1)):y(3)];

x3=[x(3):x(4)];
y3=[y(3):((y(4)-y(3))./(length(x3)-1)):y(4)];

x4=[x(4):x(5)];
y4=[y(4):((y(5)-y(4))./(length(x4)-1)):y(5)];

x5=[x(5):x(6)];
y5=[y(5):((y(6)-y(5))./(length(x5)-1)):y(6)];

x6=[x(6):-1:x(7)];
y6=[y(6):((y(7)-y(6))./(length(x6)-1)):y(7)];

x7=[x(8):-1:x(7)];
y7=[y(7):((y(8)-y(7))./(length(x7)-1)):y(8)];

x8=[x(8):-1:x(9)];
y8=[y(8):((y(9)-y(8))./(length(x8)-1)):y(9)];

y9=[y(9):y(10)];
x9=[x(9):((x(10)-x(9))./(length(y9)-1)):x(10)];

y10=[y(10):y(1)];
x10=[x(10):((x(1)-x(10))./(length(y10)-1)):x(1)];


GI(1,:)=[x1 x2 x3 x4 x5 x6 x7 x8 x9 x10];
GI(2,:)=[y1 y2 y3 y4 y5 y6 y7 y8 y9 y10];


hold on;
plot(GI(1,:),GI(2,:),'r');

LTraj=dsearchn(GI',Traj(:,3:4));%
if size(Traj,2)>=5
    LTrajS=dsearchn(GI',Traj(:,5:6));%
else
    LTrajS=[];
end

pause(2);

fprintf('Check a middle point of the 2 treadmills with the mouse, respectively\n');
[x,y]=ginput(2);


GI2(1,:)=[x(1) x(2)];
GI2(2,:)=[y(1) y(2)];
plot(GI2(1,:),GI2(2,:),'go');
TMpoints=dsearchn(GI',GI2');%
pause(2);
close(h);
return;