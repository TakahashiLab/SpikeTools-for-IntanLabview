function [LTraj,GI]=fitLinear(Traj)
h=figure;
plot(Traj(:,3),Traj(:,4));
fprintf(['Check the four corners of the rectangule with the mouse\n' ...
         'in a clockwise direction\n']);
fprintf('Please start at the upper left corner');
[x,y]=ginput(4)


x1=[x(1):x(2)];
y1=[ones(1,length(x1))*y(1)];

y2=[y(1):-1:y(3)];
x2=[ones(1,length(y2))*x(2)];

x3=[x(2):-1:x(4)];
y3=[ones(1,length(x3))*y(3)];

y4=[y(3):y(1)];
x4=[ones(1,length(y4))*x(1)];
GI(1,:)=[x1 x2 x3 x4];
GI(2,:)=[y1 y2 y3 y4];

hold on;
plot(GI(1,:),GI(2,:),'r');

LTraj=dsearchn(GI',Traj(:,3:4));

close(h);
return;