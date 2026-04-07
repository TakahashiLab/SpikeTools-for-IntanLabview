function [A,B] = align_pairs(A,B)
[~, ia, ib] = intersect(A.Pair, B.Pair, 'stable');
A = A(ia,:); B = B(ib,:);
end
