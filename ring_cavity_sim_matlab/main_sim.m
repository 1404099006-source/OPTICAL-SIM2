
function main_sim()
clear; clc;

% ---- add project subfolders to MATLAB path ----
rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir);
addpath(fullfile(rootDir,'models'));
addpath(fullfile(rootDir,'sensors'));
addpath(fullfile(rootDir,'controllers'));
addpath(fullfile(rootDir,'utils'));

P = make_params();
rng(P.seed);

% =============================
% 0) Safety defaults
% =============================
if ~isfield(P,'dz_approach');         P.dz_approach = +2.0; end    % um/step, z-up
if ~isfield(P,'duv_coarse');          P.duv_coarse  =  2.0; end
if ~isfield(P,'duv_fine');            P.duv_fine    =  0.5; end
if ~isfield(P,'force_tol');           P.force_tol   =  0.12; end
if ~isfield(P,'force_settle_steps');  P.force_settle_steps = 8; end
if ~isfield(P,'fine_iter_per_level'); P.fine_iter_per_level = 20; end
if ~isfield(P,'force_targets');       P.force_targets = [0.8 2 4 6 8]; end
if ~isfield(P,'L_thresh_ppm');        P.L_thresh_ppm = 1200; end
if ~isfield(P,'e_enter_cont');        P.e_enter_cont = 0.55; end   % 你用光路e时通常要放宽
if ~isfield(P,'Juv_step_um');         P.Juv_step_um  = 1.0; end
if ~isfield(P,'Juv_lambda');          P.Juv_lambda   = 1e-3; end
if ~isfield(P,'k_uv');                P.k_uv         = 1.0; end

% Loss quad-fit params (you said you've added them)
if ~isfield(P,'loss_probe_h_um');     P.loss_probe_h_um = 2.0; end
if ~isfield(P,'duv_ls_max');          P.duv_ls_max      = 1.0; end
if ~isfield(P,'e_guard');             P.e_guard         = 1.0; end
if ~isfield(P,'ls_lambda');           P.ls_lambda       = 1e-3; end
if ~isfield(P,'ls_gain');             P.ls_gain         = 1.0; end

% Force-hold gains (z-up)
if ~isfield(P,'dz_max');  P.dz_max = 0.5; end
if ~isfield(P,'kF_p');    P.kF_p   = 0.30; end
if ~isfield(P,'kF_i');    P.kF_i   = 0.00; end
if ~isfield(P,'dz_lpf');  P.dz_lpf = 0.90; end

if P.dz_approach <= 0
    error("P.dz_approach must be positive (z-up increases z, pushing toward cavity).");
end

% =============================
% 1) Initialize true world
% x_true is lens pose relative cavity {C}: [u v z thx thy]
% z starts negative: lens below cavity; moving up increases z.
% =============================
if ~isfield(P,'xC_init')
    P.xC_init = [80; -60; -200; deg2rad(0.02); deg2rad(-0.015)];
end
if ~isfield(P,'xR_to_C_bias')
    P.xR_to_C_bias = zeros(5,1);
end
Xtrue.xr = P.xC_init + P.xR_to_C_bias;
Xtrue.Dc = zeros(5,1);
Xtrue.x  = Xtrue.xr + Xtrue.Dc;

state.F_prev = 0;
state.Keff   = 0;
state.gmin   = NaN;     % signed gap proxy
state.g_geo  = NaN;     % geometric gap
state.Fz     = 0;
state.Mx     = 0;
state.My     = 0;
state.seated = false;
state.Fx     = 0;
state.Fy     = 0;
state.Mz     = 0;
state.Ac     = 0;
state.Fz_meas = 0;
state.stick_ratio = 1;
state.slip_ratio  = 0;

% =============================
% 2) Logs
% =============================
log_x     = nan(5, P.N);
log_xr    = nan(5, P.N);
log_duv   = nan(2, P.N);
log_dz    = nan(1, P.N);
log_F     = nan(1, P.N);
log_K     = nan(1, P.N);
log_g     = nan(1, P.N);
log_ggeo  = nan(1, P.N);
log_e     = nan(2, P.N);
log_Ltrue = nan(1, P.N);
log_Lfit  = nan(1, P.N);
log_ok    = false(1, P.N);
log_psucc = nan(1, P.N);
log_mode  = strings(1,P.N);
log_Fmeas = nan(1, P.N);
log_Fx    = nan(1, P.N);
log_Fy    = nan(1, P.N);
log_Ac    = nan(1, P.N);
log_stick = nan(1, P.N);
log_slip  = nan(1, P.N);

% =============================
% 3) Mode / counters
% =============================
mode = "approach";    % approach -> seat -> coarse -> cont_settle -> cont_fine
N_fit = 0;

level = 1;
F_target = P.force_targets(level);
settle_cnt = 0;
fine_iter  = 0;

% force controller internal
force_ctl.dz_bias = 0;
force_ctl.intF    = 0;

% =============================
% 4) Main loop
% =============================
for k = 1:P.N

    action.duv = [0;0];
    action.dz  = 0;
    action.do_seat = false;

    log_mode(k) = mode;

    % --- sensors ---
    e = sensor_vision(Xtrue.x, P);
    log_e(:,k) = e;

    log_Ltrue(k) = loss_true(Xtrue.x, P);

    % --- contact decision ---
    in_contact = false;
    if isfield(state,'gmin') && isfinite(state.gmin)
        in_contact = (state.gmin <= 0);
    else
        in_contact = (state.Fz > 1e-6);
    end

    % --- debug print ---
    if mod(k,200)==0
        fprintf("k=%d, mode=%s, z=%.2f, F=%.3f, gmin=%.2f, ||e||=%.3f\n", ...
            k, mode, Xtrue.x(3), state.Fz, state.gmin, norm(e));
    end

    % =========================
    % Control by mode
    % =========================
    if mode == "approach"
        if ~in_contact
            action.dz = P.dz_approach;      % move up to contact
        else
            mode = "seat";
        end

    elseif mode == "seat"
        action.do_seat = true;
        mode = "coarse";

    elseif mode == "coarse"

        % --- force hold (only meaningful after contact) ---
        if ~in_contact
            action.dz = P.dz_approach;
        else
            [dz_hold, force_ctl] = force_hold_dz_smooth_zup(state.Fz, F_target, P, force_ctl);
            action.dz = dz_hold;
        end

        % --- uv coarse step: one-step least squares on e (numeric Jacobian) ---
        action.duv = uv_step_from_J_numeric(Xtrue.x, e, P, P.duv_coarse);

        % --- gate into continuation ---
        if in_contact && abs(state.Fz - F_target) <= P.force_tol && norm(e) <= P.e_enter_cont
            fprintf("Coarse done at step %d, ||e||=%.3f, entering continuation\n", k, norm(e));
            mode = "cont_settle";

            % log a discrete fitted loss at transition (optional)
            [Lfit0, ok0, info0] = sensor_lossfit(Xtrue.x, e, state.Keff, state.seated, P);
            N_fit = N_fit + 1;
            log_Lfit(k)  = Lfit0;
            log_ok(k)    = ok0;
            log_psucc(k) = info0.p_succ;

            settle_cnt = 0;
            fine_iter  = 0;

            fprintf("Level %d/%d: target F=%.2fN\n", level, numel(P.force_targets), F_target);
        end

    elseif mode == "cont_settle"

        [dz_hold, force_ctl] = force_hold_dz_smooth_zup(state.Fz, F_target, P, force_ctl);
        action.dz = dz_hold;

        if abs(state.Fz - F_target) <= P.force_tol
            settle_cnt = settle_cnt + 1;
        else
            settle_cnt = 0;
        end

        % record discrete fitted loss during settle (optional)
        [LfitS, okS, infoS] = sensor_lossfit(Xtrue.x, e, state.Keff, state.seated, P);
        N_fit = N_fit + 1;
        log_Lfit(k)  = LfitS;
        log_ok(k)    = okS;
        log_psucc(k) = infoS.p_succ;

        if settle_cnt >= P.force_settle_steps
            mode = "cont_fine";
            fine_iter = 0;
        end

    elseif mode == "cont_fine"

        % ---- force hold ----
        [dz_hold, force_ctl] = force_hold_dz_smooth_zup(state.Fz, F_target, P, force_ctl);
        action.dz = dz_hold;

        fine_iter = fine_iter + 1;

        % =====================================================
        % 关键：用你写好的 2D 二次拟合 loss_quadfit_uv_step 一步更新 uv
        % =====================================================
        [duv_cmd, ~] = loss_quadfit_uv_step(Xtrue, state, P, ...
            P.loss_probe_h_um, P.duv_ls_max, P.e_guard);
        action.duv = duv_cmd;

        % ---- record discrete fitted loss point (真实扫频点) ----
        [Lnow, oknow, infonow] = sensor_lossfit(Xtrue.x, e, state.Keff, state.seated, P);
        N_fit = N_fit + 1;
        log_Lfit(k)  = Lnow;
        log_ok(k)    = oknow;
        log_psucc(k) = infonow.p_succ;

        if isfinite(Lnow) && oknow && (Lnow < P.L_thresh_ppm)
            fprintf("Reached target loss < %.0f ppm at step %d (level %d), Lfit=%.1f ppm, F=%.2fN, N_fit=%d\n", ...
                P.L_thresh_ppm, k, level, Lnow, state.Fz, N_fit);
            break;
        end

        % ---- force level budget ----
        if fine_iter >= P.fine_iter_per_level
            fprintf("Level %d done at step %d. F=%.2fN, N_fit=%d\n", level, k, state.Fz, N_fit);

            level = level + 1;
            if level > numel(P.force_targets)
                fprintf("Continuation finished. Final N_fit=%d\n", N_fit);
                break;
            else
                F_target = P.force_targets(level);
                fprintf("Level %d/%d: target F=%.2fN\n", level, numel(P.force_targets), F_target);

                mode = "cont_settle";
                settle_cnt = 0;
                fine_iter  = 0;
            end
        end
    end

    % =========================
    % Update truth world
    % =========================
    [Xtrue, state] = truth_update(Xtrue, state, action, P);

    % logs
    log_x(:,k) = Xtrue.x;
    log_xr(:,k) = Xtrue.xr;
    log_duv(:,k) = action.duv(:);
    log_dz(k) = action.dz;
    log_F(k)   = state.Fz;
    log_K(k)   = state.Keff;
    log_g(k)   = state.gmin;
    if isfield(state,'g_geo'), log_ggeo(k) = state.g_geo; end
    if isfield(state,'Fz_meas'), log_Fmeas(k) = state.Fz_meas; end
    if isfield(state,'Fx'), log_Fx(k) = state.Fx; end
    if isfield(state,'Fy'), log_Fy(k) = state.Fy; end
    if isfield(state,'Ac'), log_Ac(k) = state.Ac; end
    if isfield(state,'stick_ratio'), log_stick(k) = state.stick_ratio; end
    if isfield(state,'slip_ratio'), log_slip(k) = state.slip_ratio; end

end

Kend = find(isfinite(log_F), 1, 'last');
if isempty(Kend), Kend = 1; end
fprintf("Simulation ended at step %d, total fits N_fit=%d\n", Kend, N_fit);

% =============================
% (A) Reference poses
% =============================
x_geo  = zeros(5,1);
e_norm_all = vecnorm(log_e(:,1:Kend), 2, 1);
[~, k_opt] = min(e_norm_all);
x_opt = log_x(:, k_opt);
Ltrue_all = log_Ltrue(1:Kend);
[~, k_loss] = min(Ltrue_all);
x_loss = log_x(:, k_loss);

fprintf("\n[Ref poses]\n");
fprintf("k_opt (min ||e||)   = %d, ||e||=%.4f\n", k_opt, e_norm_all(k_opt));
fprintf("k_loss(min L_true)  = %d, L_true=%.2f ppm\n", k_loss, Ltrue_all(k_loss));
fprintf("x_opt  = [%.3f %.3f %.3f %.6g %.6g]\n", x_opt(1),x_opt(2),x_opt(3),x_opt(4),x_opt(5));
fprintf("x_loss = [%.3f %.3f %.3f %.6g %.6g]\n\n", x_loss(1),x_loss(2),x_loss(3),x_loss(4),x_loss(5));

% =============================
% (B) Build deviation signals
% =============================
X = log_x(:,1:Kend);
Xr = log_xr(:,1:Kend);
Dev_geo  = X - x_geo;
Dev_cmd  = X - Xr;

dev5_geo = vecnorm(Dev_geo, 2, 1);
devuv_geo  = vecnorm(Dev_geo(1:2,:), 2, 1);

% =============================
% (C) Plots
% =============================
close all;

dec = max(1, floor(Kend/600));
idx = 1:dec:Kend;

u = X(1,:);  v = X(2,:);  z = X(3,:);  F = log_F(1:Kend);
ur = Xr(1,:); vr = Xr(2,:); zr = Xr(3,:);

% ---- Figure 1: Force + gap ----
figure('Name','Force & Contact','Color','w');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

nexttile; plot(idx, log_F(idx), 'LineWidth',1.2); hold on; grid on;
if any(isfinite(log_Fmeas(1:Kend)))
    plot(idx, log_Fmeas(idx), '--', 'LineWidth',1.1);
    legend({'Fz true','Fz meas'}, 'Location','best');
end
xlabel('step'); ylabel('Fz (N)'); title('Normal force');

nexttile; plot(idx, log_K(idx), 'LineWidth',1.2); grid on;
xlabel('step'); ylabel('K_{eff}'); title('Effective stiffness');

nexttile;
h1 = plot(idx, log_g(idx), 'LineWidth',1.2); hold on;
if any(isfinite(log_ggeo(1:Kend)))
    h2 = plot(idx, log_ggeo(idx), 'LineWidth',1.2);
    legend([h1 h2], {'g\_signed','g\_geo'}, 'Location','best');
else
    legend(h1, 'g\_signed', 'Location','best');
end
xlabel('step'); ylabel('gap (\mum)'); title('Gap signals');

nexttile; plot(idx, z(idx), 'LineWidth',1.2); grid on;
xlabel('step'); ylabel('z (\mum)'); title('True z position');

% ---- Figure 2: Spot/aperture + loss ----
figure('Name','Optics','Color','w');
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

nexttile;
plot(idx, log_e(1,idx), 'LineWidth',1.2); hold on; grid on;
plot(idx, log_e(2,idx), 'LineWidth',1.2);
plot(idx, vecnorm(log_e(:,idx),2,1), 'LineWidth',1.2);
xlabel('step'); ylabel('e (normalized)');
title('Spot/aperture relative error');
legend({'e_y','e_z','||e||'},'Location','best');

nexttile;
plot(idx, log_Ltrue(idx), 'LineWidth',1.2); grid on; hold on;
fit_idx = find(isfinite(log_Lfit(1:Kend)));
if ~isempty(fit_idx)
    scatter(fit_idx, log_Lfit(fit_idx), 10, 'filled');
end
yline(P.L_thresh_ppm, '--', sprintf('%.0f ppm', P.L_thresh_ppm));
xlabel('step'); ylabel('Loss (ppm)');
title('Loss curve (true + fitted samples)');
legend({'L_{true}','L_{fit}','threshold'},'Location','best');

% ---- Figure 3: Robot trajectory vs true ----
figure('Name','Trajectory','Color','w');
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

nexttile; hold on; grid on; axis equal;
plot(ur, vr, 'LineWidth',1.2);
plot(u,  v,  'LineWidth',1.2);
xlabel('u (\mum)'); ylabel('v (\mum)');
title('In-plane trajectory');
legend({'commanded (x_r)','true (x)'},'Location','best');

nexttile;
scatter3(u, v, z, 14, F, 'filled'); grid on;
xlabel('u (\mum)'); ylabel('v (\mum)'); zlabel('z (\mum)');
title('True trajectory (color=Fz)'); cb=colorbar; ylabel(cb,'Fz (N)');
view(45,25); axis vis3d;

% ---- Figure 4: Pose error to cavity zero ----
figure('Name','Pose error to cavity zero','Color','w');
tiledlayout(3,2,'Padding','compact','TileSpacing','compact');

nexttile; plot(1:Kend, Dev_geo(1,:), 'LineWidth',1.1); grid on;
xlabel('step'); ylabel('u (\mum)'); title('\Delta u to geo(0)');

nexttile; plot(1:Kend, Dev_geo(2,:), 'LineWidth',1.1); grid on;
xlabel('step'); ylabel('v (\mum)'); title('\Delta v to geo(0)');

nexttile; plot(1:Kend, Dev_geo(3,:), 'LineWidth',1.1); grid on;
xlabel('step'); ylabel('z (\mum)'); title('\Delta z to geo(0)');

nexttile; plot(1:Kend, Dev_geo(4,:), 'LineWidth',1.1); grid on;
xlabel('step'); ylabel('\theta_x (rad)'); title('\Delta \theta_x to geo(0)');

nexttile; plot(1:Kend, Dev_geo(5,:), 'LineWidth',1.1); grid on;
xlabel('step'); ylabel('\theta_y (rad)'); title('\Delta \theta_y to geo(0)');

nexttile; plot(1:Kend, dev5_geo, 'LineWidth',1.1); grid on;
xlabel('step'); ylabel('||x - 0||'); title('5DoF norm to geo(0)');

% ---- Figure 5: Commanded vs true mismatch ----
figure('Name','Command vs true mismatch','Color','w');
tiledlayout(3,2,'Padding','compact','TileSpacing','compact');

nexttile; plot(1:Kend, Dev_cmd(1,:), 'LineWidth',1.1); grid on;
xlabel('step'); ylabel('\Delta u (\mum)'); title('x - x_r (u)');

nexttile; plot(1:Kend, Dev_cmd(2,:), 'LineWidth',1.1); grid on;
xlabel('step'); ylabel('\Delta v (\mum)'); title('x - x_r (v)');

nexttile; plot(1:Kend, Dev_cmd(3,:), 'LineWidth',1.1); grid on;
xlabel('step'); ylabel('\Delta z (\mum)'); title('x - x_r (z)');

nexttile; plot(1:Kend, Dev_cmd(4,:), 'LineWidth',1.1); grid on;
xlabel('step'); ylabel('\Delta \theta_x (rad)'); title('x - x_r (\theta_x)');

nexttile; plot(1:Kend, Dev_cmd(5,:), 'LineWidth',1.1); grid on;
xlabel('step'); ylabel('\Delta \theta_y (rad)'); title('x - x_r (\theta_y)');

nexttile; plot(1:Kend, devuv_geo, 'LineWidth',1.1); grid on;
xlabel('step'); ylabel('||[u v]|| (\mum)'); title('UV norm to geo(0)');


% ---- Figure 6: Convergence dashboard (key metrics) ----
figure('Name','Convergence dashboard','Color','w');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

x_star = zeros(5,1);
if isfield(P,'x_star') && numel(P.x_star)==5
    x_star = P.x_star(:);
end
Dev_star = X - x_star;
pose_uv_to_star = vecnorm(Dev_star(1:2,:),2,1);
pose_ang_to_star = vecnorm(Dev_star(4:5,:),2,1);

nexttile;
plot(1:Kend, F, 'LineWidth',1.2); grid on; hold on;
yline(P.force_targets(min(level,numel(P.force_targets))), '--');
xlabel('step'); ylabel('Fz (N)'); title('Force build-up and hold');

nexttile;
plot(1:Kend, vecnorm(log_e(:,1:Kend),2,1), 'LineWidth',1.2); grid on;
xlabel('step'); ylabel('||e||'); title('Spot/aperture convergence');

nexttile;
plot(1:Kend, log_Ltrue(1:Kend), 'LineWidth',1.2); grid on; hold on;
fit_idx2 = find(isfinite(log_Lfit(1:Kend)));
if ~isempty(fit_idx2)
    scatter(fit_idx2, log_Lfit(fit_idx2), 10, 'filled');
end
yline(P.L_thresh_ppm, '--', sprintf('%.0f ppm target', P.L_thresh_ppm));
xlabel('step'); ylabel('Loss (ppm)'); title('Loss convergence');

nexttile;
yyaxis left;
plot(1:Kend, pose_uv_to_star, 'LineWidth',1.2); ylabel('||[u,v]-[u*,v*]|| (\mum)');
yyaxis right;
plot(1:Kend, pose_ang_to_star, 'LineWidth',1.2); ylabel('||[\theta_x,\theta_y]-[\theta_x^*,\theta_y^*]|| (rad)');
grid on; xlabel('step');
title('Pose convergence to target');

% ---- Figure 7: Contact onset zoom (adhesion risk indicator) ----
k_contact = find(log_F(1:Kend) > 1e-6, 1, 'first');
if ~isempty(k_contact)
    w = 80;
    ks = max(1, k_contact-w):min(Kend, k_contact+w);
    figure('Name','Contact onset zoom','Color','w');
    tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

    nexttile;
    plot(ks, log_F(ks), 'LineWidth',1.2); grid on;
    xlabel('step'); ylabel('Fz (N)');
    title(sprintf('Force around first contact (k=%d)', k_contact));

    nexttile;
    plot(ks, log_g(ks), 'LineWidth',1.2); hold on; grid on;
    if any(isfinite(log_ggeo(ks)))
        plot(ks, log_ggeo(ks), 'LineWidth',1.2);
        legend({'g\_signed','g\_geo'},'Location','best');
    end
    yline(0,'--');
    xlabel('step'); ylabel('gap (\mum)');
    title('Gap transition near contact');
end


% ---- Figure 8: Friction & contact area ----
figure('Name','Friction and contact area','Color','w');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

nexttile;
plot(idx, log_Ac(idx), 'LineWidth',1.2); grid on;
if isfield(P,'Ro_um') && isfield(P,'Ri_um')
    A_ring = pi * (P.Ro_um^2 - P.Ri_um^2);
    hold on;
    plot(idx, log_Ac(idx)/A_ring, '--', 'LineWidth',1.1);
    legend({'A_c (um^2)','A_c/A_{ring}'}, 'Location','best');
end
xlabel('step'); ylabel('Area'); title('Contact area evolution');

nexttile;
plot(idx, log_Fx(idx), 'LineWidth',1.2); hold on; grid on;
plot(idx, log_Fy(idx), 'LineWidth',1.2);
xlabel('step'); ylabel('Force (N)'); title('Tangential friction force');
legend({'F_x','F_y'}, 'Location','best');

nexttile;
plot(idx, hypot(log_Fx(idx), log_Fy(idx)), 'LineWidth',1.2); grid on;
xlabel('step'); ylabel('||F_t|| (N)'); title('Tangential force magnitude');

nexttile;
plot(idx, log_stick(idx), 'LineWidth',1.2); hold on; grid on;
plot(idx, log_slip(idx), 'LineWidth',1.2);
xlabel('step'); ylabel('ratio'); title('Local stick/slip ratio');
legend({'stick ratio','slip ratio'}, 'Location','best');

end % ===== end main_sim =====


% ==========================================================
% Local helper: Force hold controller (z-up convention)
% dz>0 increases contact / force; dz<0 backs off.
% ==========================================================
function [dz, ctl] = force_hold_dz_smooth_zup(Fz, F_target, P, ctl)
eF = F_target - Fz;
ctl.intF = ctl.intF + eF;
raw = P.kF_p * eF + P.kF_i * ctl.intF;

raw_dz = +raw;  % z-up: +eF => push up
ctl.dz_bias = P.dz_lpf * ctl.dz_bias + (1-P.dz_lpf) * raw_dz;
dz = max(-P.dz_max, min(P.dz_max, ctl.dz_bias));
end


% ==========================================================
% Local helper: numeric Jacobian + damped LS for e -> uv
% ==========================================================
function duv_cmd = uv_step_from_J_numeric(x, e, P, duv_lim)
J = estimate_Juv_numeric(x, P);

lam = P.Juv_lambda;
A = (J.'*J + lam*eye(2));
b = (J.'*e);
duv = -(A \ b);                 % um

duv = P.k_uv * duv;

duv_cmd = duv;
duv_cmd(1) = max(-duv_lim, min(duv_lim, duv_cmd(1)));
duv_cmd(2) = max(-duv_lim, min(duv_lim, duv_cmd(2)));

if any(~isfinite(duv_cmd)) || norm(J,'fro') < 1e-12
    duv_cmd = [-duv_lim*sign(e(1)); -duv_lim*sign(e(2))];
end
end

function J = estimate_Juv_numeric(x, P)
du = P.Juv_step_um;
dv = P.Juv_step_um;

e0 = sensor_vision(x, P);

x1 = x; x1(1) = x1(1) + du;
e1 = sensor_vision(x1, P);

x2 = x; x2(2) = x2(2) + dv;
e2 = sensor_vision(x2, P);

J = [ (e1 - e0)/du, (e2 - e0)/dv ];   % 2x2
end
