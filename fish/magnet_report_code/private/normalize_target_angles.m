function theta = normalize_target_angles(targetDirDeg, N)
if isscalar(targetDirDeg)
    theta = repmat(targetDirDeg, 1, N);
else
    theta = targetDirDeg(:).';
    if numel(theta) < N
        theta = [theta, nan(1, N-numel(theta))];
    elseif numel(theta) > N
        theta = theta(1:N);
    end
end
end
