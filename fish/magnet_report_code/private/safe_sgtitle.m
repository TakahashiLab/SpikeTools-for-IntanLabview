function safe_sgtitle(parentOrFig, titleStr)
titleStr = safe_text(titleStr);
try
    sgtitle(parentOrFig, titleStr, 'Interpreter','none');
catch
    try
        sgtitle(titleStr, 'Interpreter','none');
    catch
        if isempty(get(0,'CurrentFigure')), figure('Color','w'); end
        annotation(gcf,'textbox',[0 0.95 1 0.05], ...
            'String', titleStr, 'EdgeColor','none', ...
            'HorizontalAlignment','center', 'Interpreter','none');
    end
end
end
