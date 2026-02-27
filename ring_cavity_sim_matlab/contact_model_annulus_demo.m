function log = contact_model_annulus_demo()
%CONTACT_MODEL_ANNULUS_DEMO
% Standalone annulus contact simulation demo (quasi-static + stick-slip friction).
% This file is independent from main_sim.m, used to validate contact model behavior.
%
% Units:
%   position: um, angle: rad, force: N, moment: N*um

%% 1) Parameters
P = default_params();

% Build annulus grid
P = build_annulus_grid(P);

%% 2) Robot command profile (standalone test trajectory)
N = P.N;
[zr_cmd, ur_cmd, vr_cmd] = build_command_profile(P);

%% 3) State initialization
x = P.x0(:); % true [u,v,z,thx,thy]
s_t = [0;0];                                     % tangential stick displacement (um)

% logs
log.x = nan(5,N);
log.zr = zr_cmd;
log.ur = ur_cmd;
log.vr = vr_cmd;
log.Fz = nan(1,N);
log.Fz_meas = nan(1,N);
log.Mx = nan(1,N);
log.My = nan(1,N);
log.Fx = nan(1,N);
log.Fy = nan(1,N);
log.Ac = nan(1,N);
log.gmin = nan(1,N);
log.stick = false(1,N);
log.residual_norm = nan(1,N);
log.duv_cmd = nan(2,N);
log.dx_true = nan(5,N);
log.dx_cmd = nan(5,N);

% for contact map snapshots
snap_idx = unique(round([0.20,0.45,0.70,0.92] * N));
snapshots = struct('k',{},'p',{},'g',{},'x',{});

%% 4) Main loop (quasi-static per step)
for k = 1:N
    xr = [ur_cmd(k); vr_cmd(k); zr_cmd(k); 0; 0];

    % weak robot uv -> nominal tilt coupling
    th_uv = P.K_uv_to_th * xr(1:2);
    th_ref = xr(4:5) + th_uv;

    % solve static equilibrium for [z,thx,thy]
    q0 = x(3:5);
    [q, out] = solve_equilibrium_qs(q0, xr, th_ref, P);

    % update full state
    x(1:2) = xr(1:2);
    x(3:5) = q;

    % evaluate contact and friction
    [Fz, Mx, My, Keff, g_geo, gmin, Ac, p, g] = annulus_contact_eval(x, P); %#ok<ASGLU>

    if k == 1
        duv = [0;0];
        dz_cmd = 0;
        dx_cmd = zeros(5,1);
        dx_true = zeros(5,1);
    else
        duv = [ur_cmd(k)-ur_cmd(k-1); vr_cmd(k)-vr_cmd(k-1)];
        dz_cmd = zr_cmd(k)-zr_cmd(k-1);
        dx_cmd = [duv; dz_cmd; 0; 0];
        dx_true = x - log.x(:,k-1);
    end

    [Fx, Fy, s_t, is_stick] = tangential_friction_step(Fz, Ac, duv, s_t, P);

    Fz_meas = Fz + P.alpha_fx2fz * Fx + P.alpha_fy2fz * Fy;

    % logs
    log.x(:,k) = x;
    log.Fz(k) = Fz;
    log.Fz_meas(k) = Fz_meas;
    log.Mx(k) = Mx;
    log.My(k) = My;
    log.Fx(k) = Fx;
    log.Fy(k) = Fy;
    log.Ac(k) = Ac;
    log.gmin(k) = gmin;
    log.stick(k) = is_stick;
    log.residual_norm(k) = out.residual_norm;
    log.duv_cmd(:,k) = duv;
    log.dx_cmd(:,k) = dx_cmd;
    log.dx_true(:,k) = dx_true;

    if any(k == snap_idx)
        snapshots(end+1).k = k; %#ok<AGROW>
        snapshots(end).p = p;
        snapshots(end).g = g;
        snapshots(end).x = x;
    end
end

log.snapshots = snapshots;
log.P = P;

%% 5) Print summary
k_contact = find(log.Fz > 1e-6, 1, 'first');
if isempty(k_contact), k_contact = NaN; end
fprintf('\n[Annulus Contact Demo]\n');
fprintf('First contact step: %g\n', k_contact);
fprintf('Final Fz=%.4f N, Fz_meas=%.4f N\n', log.Fz(end), log.Fz_meas(end));
fprintf('Final theta=[%.4e, %.4e] rad\n', log.x(4,end), log.x(5,end));
fprintf('Max contact area ratio Ac/A_ring = %.3f\n', max(log.Ac)/P.A_ring_um2); 
fprintf('Stick ratio = %.3f\n', mean(log.stick));
fprintf('Mean |dx_true- dx_cmd| in [u,v,z] = [%.3e %.3e %.3e]\n', ...
    mean(abs(log.dx_true(1,:)-log.dx_cmd(1,:))), ...
    mean(abs(log.dx_true(2,:)-log.dx_cmd(2,:))), ...
    mean(abs(log.dx_true(3,:)-log.dx_cmd(3,:))));

%% 6) Visualization
plot_main_results(log);
plot_contact_snapshots(log);
plot_step_motion(log);

end


function P = default_params()
P.N = 420;

% geometry (um)
P.Ro_um = 11000;
P.Ri_um = 9000;
P.grid_step_um = 300;

% normal contact (Winkler)
P.k_w = 1.0e-9;       % N/um^3
P.g0_um = 0;

% compliance / equilibrium (robot->true)
P.k_z = 0.55;         % N/um
P.k_thx = 1.8e7;      % N*um/rad
P.k_thy = 1.8e7;
P.eq_max_iter = 25;
P.eq_relax = [0.7;0.7;0.7];
P.eq_tol = [2e-4;2e-7;2e-7];
P.dtheta_max = deg2rad(0.0035); % per iteration cap

% weak uv->theta coupling [rad/um]
P.K_uv_to_th = [0, 1.8e-9; -1.6e-9, 0];

% friction (lumped)
P.mu = 0.20;
P.k_t = 0.30 * P.k_w;      % N/um^3
P.alpha_fx2fz = 0.006;      % sensor axis misalignment projection
P.alpha_fy2fz = -0.004;

% initial state
P.x0 = [0; 0; -25; deg2rad(0.006); -deg2rad(0.005)];

% command profile settings (tunable)
P.k1 = 140;  % approach phase end
P.k2 = 300;  % build-up phase end
P.zr_start = -25;
P.zr_contact_hint = 0.8;
P.zr_end = 5.0;
P.uv_start = 120;
P.uv_amp_u = 6.0;
P.uv_amp_v = 4.5;
P.uv_Tu = 160;
P.uv_Tv = 200;
end


function P = build_annulus_grid(P)
xs = -P.Ro_um:P.grid_step_um:P.Ro_um;
ys = xs;
[X, Y] = meshgrid(xs, ys);
r2 = X.^2 + Y.^2;
mask = (r2 <= P.Ro_um^2) & (r2 >= P.Ri_um^2);
P.contact_x = X(mask);
P.contact_y = Y(mask);
P.contact_dA = P.grid_step_um^2;
P.A_ring_um2 = pi*(P.Ro_um^2 - P.Ri_um^2);
end




function [zr_cmd, ur_cmd, vr_cmd] = build_command_profile(P)
N = P.N;
zr_cmd = zeros(1,N);
ur_cmd = zeros(1,N);
vr_cmd = zeros(1,N);

for k = 1:N
    if k <= P.k1
        a = (k-1) / max(P.k1-1,1);
        zr_cmd(k) = (1-a) * P.zr_start + a * P.zr_contact_hint;
    elseif k <= P.k2
        a = (k-P.k1) / max(P.k2-P.k1,1);
        zr_cmd(k) = (1-a) * P.zr_contact_hint + a * (0.7*P.zr_end);
    else
        a = (k-P.k2) / max(P.N-P.k2,1);
        zr_cmd(k) = (1-a) * (0.7*P.zr_end) + a * P.zr_end;
    end

    if k > P.uv_start
        ur_cmd(k) = P.uv_amp_u * sin(2*pi*(k-P.uv_start)/P.uv_Tu);
        vr_cmd(k) = P.uv_amp_v * cos(2*pi*(k-P.uv_start)/P.uv_Tv);
    end
end
end

function [q, out] = solve_equilibrium_qs(q0, xr, th_ref, P)
% Solve for q=[z;thx;thy]
q = q0;
res = [inf;inf;inf];
for it = 1:P.eq_max_iter
    x_eval = [xr(1); xr(2); q(1); q(2); q(3)];
    [Fz, Mx, My] = annulus_contact_eval(x_eval, P);

    res = [ Fz - P.k_z   * (q(1) - xr(3));
            Mx - P.k_thx * (q(2) - th_ref(1));
            My - P.k_thy * (q(3) - th_ref(2)) ];

    if all(abs(res) <= P.eq_tol)
        break;
    end

    % Residual form: r = F_contact - K*(q-q_ref).
    % Use +r/K update so in free-space (F_contact=0) state moves toward q_ref
    % instead of diverging away.
    dq = [res(1)/max(P.k_z,1e-12);
          res(2)/max(P.k_thx,1e-12);
          res(3)/max(P.k_thy,1e-12)];
    dq = P.eq_relax .* dq;
    dq(2:3) = max(min(dq(2:3), P.dtheta_max), -P.dtheta_max);
    q = q + dq;
end

out.n_iter = it;
out.residual_norm = norm(res);
end


function [Fz, Mx, My, Keff, g_geo_min, g_signed_min, Ac, p, g] = annulus_contact_eval(x, P)
u = x(1); v = x(2); z = x(3); thx = x(4); thy = x(5);

th_uv = P.K_uv_to_th * [u; v];
thx_eff = thx + th_uv(1);
thy_eff = thy + th_uv(2);

X = P.contact_x; Y = P.contact_y; dA = P.contact_dA;

g = P.g0_um - z + thx_eff .* Y - thy_eff .* X;
g_signed_min = min(g);
g_geo_min = max(g_signed_min, 0);
w = max(-g, 0);

p = P.k_w * w;
Fz = sum(p) * dA;
Mx = sum(Y .* p) * dA;
My = -sum(X .* p) * dA;
Keff = sum(P.k_w .* (w>0)) * dA;
Ac = nnz(w>0) * dA;
end


function [Fx, Fy, s_t, is_stick] = tangential_friction_step(Fz, Ac, duv, s_t, P)
if Ac <= 0 || Fz <= 0
    Fx = 0; Fy = 0; s_t = [0;0]; is_stick = true; return;
end

s_trial = s_t + duv;
Ft_trial = (P.k_t * Ac) * s_trial;
Ft_lim = P.mu * Fz;
mag = norm(Ft_trial);

if mag <= Ft_lim || mag < 1e-15
    Ft = Ft_trial;
    s_t = s_trial;
    is_stick = true;
else
    Ft = Ft_lim * (Ft_trial / mag);
    s_t = Ft / max(P.k_t * Ac, 1e-15);
    is_stick = false;
end

Fx = Ft(1);
Fy = Ft(2);
end


function plot_main_results(log)
k = 1:numel(log.Fz);
A_ratio = log.Ac / log.P.A_ring_um2;

figure('Name','Annulus contact standalone demo','Color','w');
tiledlayout(3,2,'Padding','compact','TileSpacing','compact');

nexttile;
plot(k, log.Fz, 'LineWidth',1.2); hold on; grid on;
plot(k, log.Fz_meas, '--', 'LineWidth',1.1);
xlabel('step'); ylabel('Force (N)');
legend({'F_z true','F_{z,meas}'}, 'Location','best');
title('Normal force and measured projection');

nexttile;
plot(k, A_ratio, 'LineWidth',1.2); grid on;
xlabel('step'); ylabel('A_c / A_{ring}');
title('Contact area ratio');

nexttile;
plot(k, log.x(4,:), 'LineWidth',1.2); hold on; grid on;
plot(k, log.x(5,:), 'LineWidth',1.2);
xlabel('step'); ylabel('\theta (rad)');
legend({'\theta_x','\theta_y'}, 'Location','best');
title('Tilt evolution');

nexttile;
plot(k, log.Mx, 'LineWidth',1.2); hold on; grid on;
plot(k, log.My, 'LineWidth',1.2);
xlabel('step'); ylabel('Moment (N\cdotum)');
legend({'M_x','M_y'}, 'Location','best');
title('Contact moments');

nexttile;
plot(k, log.Fx, 'LineWidth',1.2); hold on; grid on;
plot(k, log.Fy, 'LineWidth',1.2);
xlabel('step'); ylabel('Tangential force (N)');
legend({'F_x','F_y'}, 'Location','best');
title('Stick-slip tangential force');

nexttile;
stairs(k, ~log.stick, 'LineWidth',1.2); hold on; grid on;
plot(k, log.residual_norm, 'LineWidth',1.1);
xlabel('step'); ylabel('flag / norm');
legend({'slip flag (1=slip)','equilibrium residual norm'}, 'Location','best');
title('Solver and friction state');
end


function plot_contact_snapshots(log)
if isempty(log.snapshots)
    return;
end

figure('Name','Contact pressure snapshots','Color','w');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
for i = 1:min(4,numel(log.snapshots))
    s = log.snapshots(i);
    nexttile;
    scatter(log.P.contact_x, log.P.contact_y, 9, s.p, 'filled');
    axis equal; grid on;
    cb = colorbar; ylabel(cb,'p (N/um^2)');
    xlabel('x (um)'); ylabel('y (um)');
    title(sprintf('k=%d, Fz=%.3fN, Ac/A=%.2f%%', s.k, ...
        log.Fz(s.k), 100*log.Ac(s.k)/log.P.A_ring_um2));
end
end

function plot_step_motion(log)
k = 1:numel(log.Fz);
figure('Name','Per-step commanded vs true motion','Color','w');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

nexttile;
plot(k, log.dx_cmd(1,:), 'LineWidth',1.1); hold on; grid on;
plot(k, log.dx_true(1,:), '--', 'LineWidth',1.1);
plot(k, log.dx_cmd(2,:), 'LineWidth',1.1);
plot(k, log.dx_true(2,:), '--', 'LineWidth',1.1);
xlabel('step'); ylabel('um/step');
legend({'du_{cmd}','du_{true}','dv_{cmd}','dv_{true}'}, 'Location','best');
title('In-plane step motion');

nexttile;
plot(k, log.dx_cmd(3,:), 'LineWidth',1.1); hold on; grid on;
plot(k, log.dx_true(3,:), '--', 'LineWidth',1.1);
xlabel('step'); ylabel('um/step');
legend({'dz_{cmd}','dz_{true}'}, 'Location','best');
title('Normal step motion');

nexttile;
plot(k, log.dx_true(4,:), 'LineWidth',1.1); hold on; grid on;
plot(k, log.dx_true(5,:), 'LineWidth',1.1);
xlabel('step'); ylabel('rad/step');
legend({'d\theta_x','d\theta_y'}, 'Location','best');
title('True angular increment per step');

nexttile;
plot(k, log.gmin, 'LineWidth',1.1); hold on; grid on;
yline(0,'--');
xlabel('step'); ylabel('g_{min} (um)');
title('Contact indicator per step');
end

end
