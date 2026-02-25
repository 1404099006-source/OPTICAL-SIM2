function dx = optics_dx_from_pose(x_true, P)
%OPTICS_DX_FROM_POSE  Unified optics output from pose with {C}->{O} bias removal.
% x_true: [u v z thx thy] lens pose relative cavity {C} (um,um,um,rad,rad)
% dx: [dy; dz; dty; dtz] at aperture plane (mm,mm,rad,rad) from linear optics model.

% ---- 0) defaults ----
if ~isfield(P,'optics_mirror_id'); P.optics_mirror_id = 4; end
if ~isfield(P,'xC_to_opt_bias');   P.xC_to_opt_bias   = zeros(5,1); end
if ~isfield(P,'optics_q0');        P.optics_q0        = zeros(20,1); end
if ~isfield(P,'optics') || ~isfield(P.optics,'J')
    error('P.optics.J not found. Ensure P.optics = optics_init_linear(P) in make_params.');
end

% ---- 1) remove {C}->{O} bias (THIS IS方案B的核心) ----
% After calibration, x_true == bias => optics sees zero.
x_opt = x_true(:) - P.xC_to_opt_bias(:);

u_um  = x_opt(1);
v_um  = x_opt(2);
z_um  = x_opt(3);
thx   = x_opt(4);
thy   = x_opt(5);

% ---- 2) unit conversion: um -> mm for translations ----
u_mm = u_um * 1e-3;
v_mm = v_um * 1e-3;
z_mm = z_um * 1e-3;

% ---- 3) build dq (20x1), only fill M4 block ----
dq = zeros(20,1);
m = P.optics_mirror_id;     % mirror index (4 for M4)
base = (m-1)*5;

% Convention per mirror: [dx, dy, dz, ax, ay]
dq(base + 1) = u_mm;
dq(base + 2) = v_mm;
dq(base + 3) = z_mm;
dq(base + 4) = thx;
dq(base + 5) = thy;

% ---- 4) linear optics ----
% If your J is linearized at q0: dx = J*(q - q0)
dq_rel = dq - P.optics_q0(:);
dx = P.optics.J * dq_rel;     % dx = [dy; dz; dty; dtz]
end
