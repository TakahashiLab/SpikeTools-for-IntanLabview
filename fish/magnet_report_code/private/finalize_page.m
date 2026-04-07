function finalize_page(figH, pdfPath, pngPath, pngRes, appendPDF, writePNG)
try
    exportgraphics(figH, pdfPath, 'ContentType','image', 'Append', appendPDF);
catch
    print(figH, pdfPath, '-dpdf', '-painters');
end
if writePNG
    try
        exportgraphics(figH, pngPath, 'Resolution', pngRes, 'Append', appendPDF);
    catch
        print(figH, pngPath, sprintf('-r%d', pngRes), '-dpng');
    end
end
end
