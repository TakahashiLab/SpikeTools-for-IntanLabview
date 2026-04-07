function r = inset_rect(ax, x, y, w, h)
op = get(ax,'Position'); r = [op(1)+op(3)*x, op(2)+op(4)*y, op(3)*w, op(4)*h];
end
