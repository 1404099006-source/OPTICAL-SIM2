function [Fz, Mx, My, Keff, g_geo_min, g_signed_min] = force_model(x, P)
% x = [u; v; z; thx; thy], but u,v don't affect normal contact here (plane is infinite)
z   = x(3);
thx = x(4);
thy = x(5);

% signed gap field on disk points
% g(x,y) = z + thx*y - thy*x
g = -z + thx .* P.disk_y - thy .* P.disk_x;   % z-up positive: higher z => smaller gap

g_signed_min = min(g);
g_geo_min = max(g_signed_min, 0);       % pure geometric gap, non-negative
w = max(-g, 0);                         % indentation field (um), >=0

% Winkler foundation: p = k_w * w (N/um^2)
p = P.k_w .* w;

% Integrals
dA = P.disk_dA;
Fz = sum(p) * dA;

% Moments about disk center (N*um)
Mx = sum(P.disk_y .* p) * dA;     % moment about x-axis from y lever arm
My = -sum(P.disk_x .* p) * dA;    % sign chosen so +thy (tilt about y) produces restoring My

% Effective stiffness proxy: dF/dz ≈ ∑ k_w * I(w>0) dA
Keff = sum(P.k_w .* (w>0)) * dA;
end
