function Lppm = loss_true(x_true, P)
%LOSS_TRUE  true cavity loss in ppm (optics-dominant), consistent with sensor_vision.
% Uses dx=[dy; dz; dty; dtz] from optics_dx_from_pose(x_true,P).

dx = optics_dx_from_pose(x_true, P);  % [mm; mm; rad; rad]
dy  = dx(1);  dz  = dx(2);
dty = dx(3);  dtz = dx(4);

% ---- normalize to dimensionless errors ----
if ~isfield(P,'aperture_ry_mm'); P.aperture_ry_mm = 1.0; end
if ~isfield(P,'aperture_rz_mm'); P.aperture_rz_mm = 1.0; end
ep = [dy / P.aperture_ry_mm;
      dz / P.aperture_rz_mm];

if ~isfield(P,'ang_ref_rad'); P.ang_ref_rad = deg2rad(0.02); end
ea = [dty / P.ang_ref_rad;
      dtz / P.ang_ref_rad];

% ---- quadratic loss model (ppm) ----
if ~isfield(P,'L_min');  P.L_min  = 200;  end
if ~isfield(P,'Gp_ppm'); P.Gp_ppm = 2500; end
if ~isfield(P,'Ga_ppm'); P.Ga_ppm = 1200; end
if ~isfield(P,'Gc_ppm'); P.Gc_ppm = 200;  end

pos_term = P.Gp_ppm * (ep.'*ep);
ang_term = P.Ga_ppm * (ea.'*ea);

% mild coupling (keep smooth; do NOT use abs if you want nice gradients)
cross = P.Gc_ppm * (ep(1)*ea(1) + ep(2)*ea(2));

Lppm = P.L_min + pos_term + ang_term + cross;

% optional tiny floor noise for "true" (usually keep off)
if isfield(P,'Ltrue_sigma_ppm') && P.Ltrue_sigma_ppm > 0
    Lppm = Lppm + P.Ltrue_sigma_ppm * randn();
end
end
