function p_adj = holm_adjust(p)
[p_sorted, idx] = sort(p(:));
m = numel(p_sorted);
p_holm = zeros(m,1);
for i = 1:m
    p_holm(i) = min(1, max(p_sorted(1:i) .* (m - (0:i-1))'));
end
for i = m-1:-1:1
    p_holm(i) = max(p_holm(i), p_holm(i+1));
end
p_adj = zeros(size(p));
p_adj(idx) = p_holm;
end
