function [pval, z, Rbar, mu_deg, n] = rayleigh_on_angles(alpha_rad)
alpha_rad = alpha_rad(:); alpha_rad = alpha_rad(isfinite(alpha_rad));
n = numel(alpha_rad);
if n==0, pval=NaN; z=NaN; Rbar=NaN; mu_deg=NaN; return; end
Cbar = mean(cos(alpha_rad)); Sbar = mean(sin(alpha_rad));
Rbar = hypot(Cbar, Sbar);
mu_deg = rad2deg(atan2(Sbar, Cbar));
if exist('circ_rtest','file') == 2
    [pval, z] = circ_rtest(alpha_rad);
else
    z = n * Rbar^2;
    pval = exp(sqrt(1+4*n+4*(n^2 - (n*Rbar)^2)) - (1+2*n));
    pval = min(max(pval,0),1);
end
end
