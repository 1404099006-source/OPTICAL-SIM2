function e = sensor_vision(x_true, P)
%SENSOR_VISION Vision error from optics model.
% Default (recommended): e = [dy; dz] in mm at aperture plane (physical metric).
% Optional legacy mode: e = [dy/ry; dz/rz] (dimensionless).

dx = optics_dx_from_pose(x_true, P);
dy = dx(1);   % mm
dz = dx(2);   % mm

if ~isfield(P,'vision_use_normalized'); P.vision_use_normalized = false; end
if P.vision_use_normalized
    % legacy normalized output (dimensionless)
    if ~isfield(P,'aperture_ry_mm'); P.aperture_ry_mm = 1.0; end
    if ~isfield(P,'aperture_rz_mm'); P.aperture_rz_mm = 1.0; end
    e = [dy / P.aperture_ry_mm;
         dz / P.aperture_rz_mm];

    if isfield(P,'vision_sigma') && P.vision_sigma > 0
        e = e + P.vision_sigma * randn(2,1);
    end
else
    % physical output (mm at aperture plane)
    e = [dy; dz];

    if isfield(P,'vision_sigma_mm') && P.vision_sigma_mm > 0
        e = e + P.vision_sigma_mm * randn(2,1);
    elseif isfield(P,'vision_sigma') && P.vision_sigma > 0
        % backward compatibility: treat legacy field as mm noise when using physical mode
        e = e + P.vision_sigma * randn(2,1);
    end
end
end
