function printFig(FH,fn)

fn=sprintf('%s.pdf',fn);
%FH.WindowState='maximized';
%set(FH,'PaperOrientation','landscape');
%print(FH,fn,'-dpdf','-painters','-r450','-bestfit');
%print(FH,fn,'-dpdf');
export_fig(fn,FH);
%close(FH);
return;