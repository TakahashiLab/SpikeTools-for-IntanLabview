function dispEps(EpsEdge)
Hz=25000;
cl={'r','g','k'};%ambulation,immobility,fine movement

hold on;
for Eps=[1:3]
  loop=size(EpsEdge{Eps},1);
  EpsEdge{Eps}=floor(EpsEdge{Eps});
  for i=1:loop
    x=[EpsEdge{Eps}(i,1) EpsEdge{Eps}(i,2) EpsEdge{Eps}(i,2) ...
	  EpsEdge{Eps}(i,1)];
    y=[1 1 5 5];
    

    f=fill(x,y,cl{Eps});
    set(f,'EdgeColor','none');
  end
end



return;