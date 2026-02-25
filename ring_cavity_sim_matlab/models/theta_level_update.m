function th_new = theta_level_update(th, Fz, seated, P)
%THETA_LEVEL_UPDATE Ball-joint self-leveling: tilt decays toward 0 under contact force.
% th: [thx; thy] in rad
% Fz: normal force in N
% seated: boolean
% Returns updated th_new.

th_new = th;

if ~isfield(P, "enable_leveling") || ~P.enable_leveling
    return;
end

% activate only after touching
if Fz <= P.F_touch
    return;
end

% force-dependent leveling rate rho(F) in (0, rho_theta0)
rho = P.rho_theta0 * (1 - exp(-P.beta_theta * (Fz - P.F_touch)));
rho = min(max(rho, 0), P.rho_theta0);

% desired update: th <- (1-rho)*th  (exponential decay to 0)
th_target = (1 - rho) * th;

% cap per-step change for stability
dth = th_target - th;
dth = max(min(dth, P.dtheta_max), -P.dtheta_max);
th_new = th + dth;

% optional micro jitter after seating (models tiny stick-slip)
if seated && isfield(P, "theta_jitter") && P.theta_jitter > 0
    th_new = th_new + P.theta_jitter * randn(2,1);
end
end
