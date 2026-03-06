function Lppm = loss_true(x_true, P)
%LOSS_TRUE True cavity loss in ppm.
% Default: physically interpretable model
%   Lppm = 1e6 * (ell0 + ell_clip(rho))
% where rho is the aperture-plane spot offset from optics_dx_from_pose.
% Fallback: legacy quadratic model on normalized position/angle errors.

if ~isfield(P,'use_physical_loss')
    P.use_physical_loss = false;
end

if P.use_physical_loss
    Lppm = physical_loss_ppm(x_true, P);
else
    Lppm = legacy_quadratic_loss_ppm(x_true, P);
end

% optional tiny floor noise for "true" (usually keep off)
if isfield(P,'Ltrue_sigma_ppm') && P.Ltrue_sigma_ppm > 0
    Lppm = Lppm + P.Ltrue_sigma_ppm * randn();
end

Lppm = max(0, Lppm);
end

function Lppm = physical_loss_ppm(x_true, P)
% Physics-based round-trip loss:
% ell_rt = ell0 + ell_clip(rho),    Lppm = 1e6 * ell_rt

dx = optics_dx_from_pose(x_true, P);  % [dy; dz; dty; dtz], mm/rad
rho_mm = hypot(dx(1), dx(2));

if ~isfield(P,'L0_ppm'); P.L0_ppm = 500; end
if ~isfield(P,'aperture_phys_radius_mm'); P.aperture_phys_radius_mm = 0.5; end
if ~isfield(P,'w_eff_mm')
    wt = getfield_def(P,'w_t_mm',0.395);
    ws = getfield_def(P,'w_s_mm',0.333);
    P.w_eff_mm = sqrt(max(1e-12, wt*ws));
end

a_mm = P.aperture_phys_radius_mm;
w_mm = P.w_eff_mm;

eta = aperture_overlap_eta(rho_mm, a_mm, w_mm, P);  % transmitted power fraction in [0,1]
ell_clip = max(0, 1 - eta);
ell0 = max(0, P.L0_ppm * 1e-6);

Lppm = 1e6 * (ell0 + ell_clip);
end

function eta = aperture_overlap_eta(rho_mm, a_mm, w_mm, P)
% Compute Gaussian power fraction through a circular aperture.
% rho_mm: beam center offset magnitude at aperture plane
% a_mm: aperture radius
% w_mm: 1/e^2 intensity radius

if ~isfield(P,'loss_clip_method'); P.loss_clip_method = "approx_radial"; end
method = string(P.loss_clip_method);

rho_mm = max(0, rho_mm);
a_mm = max(1e-12, a_mm);
w_mm = max(1e-12, w_mm);

switch method
    case "numeric"
        if ~isfield(P,'loss_clip_grid_n'); P.loss_clip_grid_n = 81; end
        n = max(31, round(P.loss_clip_grid_n));
        x = linspace(-a_mm, a_mm, n);
        y = linspace(-a_mm, a_mm, n);
        [X, Y] = meshgrid(x, y);
        mask = (X.^2 + Y.^2) <= a_mm^2;
        dA = (x(2)-x(1)) * (y(2)-y(1));
        I = exp(-2 * (((X - rho_mm).^2 + Y.^2) / (w_mm^2)));
        P_in = sum(I(mask)) * dA;
        P_tot = pi * w_mm^2 / 2;  % integral of exp(-2 r^2 / w^2) over R^2
        eta = P_in / max(1e-18, P_tot);
    otherwise
        % Fast approximation: equivalent centered aperture radius a-rho.
        % Exact at rho=0; smooth/monotone for control and cheap in inner loops.
        a_eff = max(0, a_mm - rho_mm);
        eta = 1 - exp(-2 * (a_eff / w_mm)^2);
end

eta = min(1, max(0, eta));
end

function Lppm = legacy_quadratic_loss_ppm(x_true, P)
% Legacy empirical model kept for compatibility.

dx = optics_dx_from_pose(x_true, P);  % [mm; mm; rad; rad]
dy  = dx(1);  dz  = dx(2);
dty = dx(3);  dtz = dx(4);

if ~isfield(P,'aperture_ry_mm'); P.aperture_ry_mm = 1.0; end
if ~isfield(P,'aperture_rz_mm'); P.aperture_rz_mm = 1.0; end
ep = [dy / P.aperture_ry_mm;
      dz / P.aperture_rz_mm];

if ~isfield(P,'ang_ref_rad'); P.ang_ref_rad = deg2rad(0.02); end
ea = [dty / P.ang_ref_rad;
      dtz / P.ang_ref_rad];

if ~isfield(P,'L_min');  P.L_min  = 200;  end
if ~isfield(P,'Gp_ppm'); P.Gp_ppm = 2500; end
if ~isfield(P,'Ga_ppm'); P.Ga_ppm = 1200; end
if ~isfield(P,'Gc_ppm'); P.Gc_ppm = 200;  end

pos_term = P.Gp_ppm * (ep.'*ep);
ang_term = P.Ga_ppm * (ea.'*ea);
cross = P.Gc_ppm * (ep(1)*ea(1) + ep(2)*ea(2));

Lppm = P.L_min + pos_term + ang_term + cross;
end

function v = getfield_def(S, name, default_val)
if isfield(S, name)
    v = S.(name);
else
    v = default_val;
end
end
