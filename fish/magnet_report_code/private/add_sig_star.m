function add_sig_star(ax, nStars)
if nargin<2, nStars = 1; end
stars = repmat('*',1,nStars);
text(ax, 0.95, 0.95, stars, 'Units','normalized', 'HorizontalAlignment','right', ...
    'VerticalAlignment','top', 'FontWeight','bold', 'FontSize', 14);
end
