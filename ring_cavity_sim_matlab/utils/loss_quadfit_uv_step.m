function [duv_cmd, info] = loss_quadfit_uv_step(Xtrue, state, P, h_um, duv_max, e_guard)
%LOSS_QUADFIT_UV_STEP  2D quadratic fit step for loss minimization in (u,v).
% Samples loss at a 3x3 stencil around current (u,v), fits
%   L(u,v) = a*u^2 + b*v^2 + c*u*v + d*u + e*v + f
% then takes one damped Newton/LS step toward the minimizer.
%
% Inputs:
%   Xtrue.x : current pose [u v z thx thy] (um, um, um, rad, rad)
%   state   : struct with Keff, seated (used by sensor_lossfit)
%   P       : params (needs P.ls_lambda, P.ls_gain, optionally P.fit_always_ok)
%   h_um    : probe step in um (e.g. 1~5)
%   duv_max : clamp per-axis command (um)
%   e_guard : max allowed ||e|| during loss step (keep beam roughly in FOV)
%
% Outputs:
%   duv_cmd : [du; dv] in um
%   info    : fitted coeffs etc.

if nargin < 4 || isempty(h_um);    h_um = 2.0; end
if nargin < 5 || isempty(duv_max); duv_max = 1.0; end
if nargin < 6 || isempty(e_guard); e_guard = inf; end

if ~isfield(P,'ls_lambda'); P.ls_lambda = 1e-3; end
if ~isfield(P,'ls_gain');   P.ls_gain   = 1.0; end

x0 = Xtrue.x(:);
u0 = x0(1);
v0 = x0(2);

% -------- 1) 3x3 stencil points (du,dv) --------
D = [  0   0;
       h_um 0;
      -h_um 0;
       0   h_um;
       0  -h_um;
       h_um h_um;
       h_um -h_um;
      -h_um h_um;
      -h_um -h_um ];

n = size(D,1);
L = nan(n,1);
Eok = false(n,1);

% -------- 2) sample loss at each point --------
for i = 1:n
    du = D(i,1); dv = D(i,2);
    xt = x0;
    xt(1) = u0 + du;
    xt(2) = v0 + dv;

    % guard: keep vision error not too large
    et = sensor_vision(xt, P);
    if norm(et) <= e_guard
        [Lfit, ok, ~] = sensor_lossfit(xt, et, state.Keff, state.seated, P);
        if ok && isfinite(Lfit)
            L(i) = Lfit;
            Eok(i) = true;
        end
    end
end

% If too few valid points, fallback to e-based LS step (safe)
if nnz(Eok) < 6
    duv_cmd = fallback_uv_from_e(x0, P, duv_max);
    info.used_fallback = true;
    info.n_valid = nnz(Eok);
    return;
end

% Keep only valid samples
D2 = D(Eok,:);
L2 = L(Eok);

% -------- 3) quadratic fit around (0,0) in local coords (du,dv) --------
% L(du,dv)= a*du^2 + b*dv^2 + c*du*dv + d*du + e*dv + f
A = [ D2(:,1).^2, D2(:,2).^2, D2(:,1).*D2(:,2), D2(:,1), D2(:,2), ones(size(D2,1),1) ];
coef = A \ L2;   % least squares
a = coef(1); b = coef(2); c = coef(3); d = coef(4); e = coef(5); f = coef(6);

% -------- 4) solve for minimizer of quadratic (damped) --------
% grad = [2a*du + c*dv + d; c*du + 2b*dv + e] = 0
H = [2*a, c; c, 2*b];
g = [d; e];

% Damping to avoid singular / ill-conditioned
lam = P.ls_lambda;
H_d = H + lam * eye(2);

duv_star = - (H_d \ g);          % local optimal step in (du,dv)

% Apply gain and clamp
duv = P.ls_gain * duv_star;
duv_cmd = duv;

duv_cmd(1) = max(-duv_max, min(duv_max, duv_cmd(1)));
duv_cmd(2) = max(-duv_max, min(duv_max, duv_cmd(2)));

% If NaN/Inf, fallback
if any(~isfinite(duv_cmd))
    duv_cmd = fallback_uv_from_e(x0, P, duv_max);
    info.used_fallback = true;
else
    info.used_fallback = false;
end

info.coef = coef;
info.n_valid = nnz(Eok);
info.L_center = L(1);
info.D = D;
info.L = L;
info.Eok = Eok;

end


% ===== fallback: use one-step LS on vision error e =====
function duv_cmd = fallback_uv_from_e(x, P, duv_max)
if ~isfield(P,'Juv_step_um'); P.Juv_step_um = 1.0; end
if ~isfield(P,'Juv_lambda');  P.Juv_lambda  = 1e-3; end
if ~isfield(P,'k_uv');        P.k_uv        = 1.0; end

e0 = sensor_vision(x, P);
du = P.Juv_step_um;
dv = P.Juv_step_um;

x1 = x; x1(1) = x1(1) + du; e1 = sensor_vision(x1, P);
x2 = x; x2(2) = x2(2) + dv; e2 = sensor_vision(x2, P);

J = [ (e1-e0)/du, (e2-e0)/dv ];
lam = P.Juv_lambda;
duv = - ( (J.'*J + lam*eye(2)) \ (J.'*e0) );
duv = P.k_uv * duv;

duv_cmd = duv;
duv_cmd(1) = max(-duv_max, min(duv_max, duv_cmd(1)));
duv_cmd(2) = max(-duv_max, min(duv_max, duv_cmd(2)));

if any(~isfinite(duv_cmd)) || norm(J,'fro') < 1e-12
    duv_cmd = [-duv_max*sign(e0(1)); -duv_max*sign(e0(2))];
end
end
