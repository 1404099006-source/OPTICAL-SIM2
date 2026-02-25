function e = sensor_vision(x_true, P)
%SENSOR_VISION  normalized vision error e from optics model.
% e = [dy/ry; dz/rz] (dimensionless)

dx = optics_dx_from_pose(x_true, P);
dy = dx(1);   % mm
dz = dx(2);   % mm

% ---- normalization radii (mm) ----
if ~isfield(P,'aperture_ry_mm'); P.aperture_ry_mm = 1.0; end
if ~isfield(P,'aperture_rz_mm'); P.aperture_rz_mm = 1.0; end

eu = dy / P.aperture_ry_mm;
ev = dz / P.aperture_rz_mm;

e = [eu; ev];

% ---- optional noise on normalized measurement ----
if isfield(P,'vision_sigma') && P.vision_sigma > 0
    e = e + P.vision_sigma * randn(2,1);
end
end
