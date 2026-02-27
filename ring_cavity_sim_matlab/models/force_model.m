function [Fz, Mx, My, Keff, g_geo_min, g_signed_min, Ac] = force_model(x, P)
%FORCE_MODEL Annulus contact model (Winkler) with optional weak uv->tilt coupling.
% x = [u; v; z; thx; thy] (um,um,um,rad,rad)

u   = x(1);
v   = x(2);
z   = x(3);
thx = x(4);
thy = x(5);

% --- contact grid: prefer annulus, fallback to legacy disk ---
if isfield(P,'contact_x') && isfield(P,'contact_y')
    Xc = P.contact_x;
    Yc = P.contact_y;
    dA = P.contact_dA;
else
    Xc = P.disk_x;
    Yc = P.disk_y;
    dA = P.disk_dA;
end

% --- weak coupling: robot uv through gripper flexibility induces tiny tilt ---
if ~isfield(P,'K_uv_to_th')
    P.K_uv_to_th = zeros(2,2); % [rad/um]
end
th_uv = P.K_uv_to_th * [u; v];
thx_eff = thx + th_uv(1);
thy_eff = thy + th_uv(2);

% signed gap field on annulus points
% g(x,y) = g0 - z + thx*y - thy*x  (small-angle plane approximation)
if ~isfield(P,'g0_um'); P.g0_um = 0; end
g = P.g0_um - z + thx_eff .* Yc - thy_eff .* Xc;

g_signed_min = min(g);
g_geo_min = max(g_signed_min, 0);   % geometric gap, non-negative
w = max(-g, 0);                      % indentation field (um)

% Winkler foundation: p = k_w * w  [N/um^2]
p = P.k_w .* w;

% Integrals
Fz = sum(p) * dA;
Mx = sum(Yc .* p) * dA;       % N*um
My = -sum(Xc .* p) * dA;      % N*um

Ac = nnz(w>0) * dA;           % um^2 contact area
Keff = sum(P.k_w .* (w>0)) * dA;   % dF/dz proxy (N/um)
end
