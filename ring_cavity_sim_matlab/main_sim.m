
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
if ~isfield(P,'e_enter_cont');        P.e_enter_cont = 0.12; end   % mm, stricter coarse->fine gate
if ~isfield(P,'min_level_before_exit'); P.min_level_before_exit = numel(P.force_targets); end
if ~isfield(P,'final_attach_enable');   P.final_attach_enable = true; end
if ~isfield(P,'dz_attach_pulse');       P.dz_attach_pulse = 0.05; end   % um/step extra push
if ~isfield(P,'attach_push_steps');     P.attach_push_steps = 6; end
if ~isfield(P,'attach_hold_steps');     P.attach_hold_steps = 12; end
if ~isfield(P,'attach_force_tol');      P.attach_force_tol = max(P.force_tol, 0.15); end
if ~isfield(P,'Juv_step_um');         P.Juv_step_um  = 1.0; end
if ~isfield(P,'Juv_lambda');          P.Juv_lambda   = 1e-3; end
if ~isfield(P,'k_uv');                P.k_uv         = 1.0; end

% Loss quad-fit params (you said you've added them)
if ~isfield(P,'loss_probe_h_um');     P.loss_probe_h_um = 2.0; end
if ~isfield(P,'loss_probe_grid_n');    P.loss_probe_grid_n = 5; end
if ~isfield(P,'duv_ls_max');          P.duv_ls_max      = 1.0; end
if ~isfield(P,'e_guard');             P.e_guard         = 0.20; end  % mm
if ~isfield(P,'ls_lambda');           P.ls_lambda       = 1e-3; end
if ~isfield(P,'ls_gain');             P.ls_gain         = 1.0; end
if ~isfield(P,'tilt_level_enable');       P.tilt_level_enable = true; end
if ~isfield(P,'tilt_probe_dth');          P.tilt_probe_dth = deg2rad(2/3600); end
if ~isfield(P,'tilt_step_max');           P.tilt_step_max = deg2rad(4/3600); end
if ~isfield(P,'tilt_gain');               P.tilt_gain = 0.5; end
if ~isfield(P,'tilt_lambda');             P.tilt_lambda = 1e-6; end
if ~isfield(P,'tilt_rcond_min');          P.tilt_rcond_min = 1e-10; end
if ~isfield(P,'tilt_force_soft');         P.tilt_force_soft = 0.20; end
if ~isfield(P,'tilt_force_hard');         P.tilt_force_hard = 0.40; end
if ~isfield(P,'tilt_done_e');             P.tilt_done_e = 0.15; end
if ~isfield(P,'tilt_max_iter_per_level'); P.tilt_max_iter_per_level = 25; end
if ~isfield(P,'tilt_probe_settle_steps'); P.tilt_probe_settle_steps = 4; end
if ~isfield(P,'tilt_mode_step_budget');   P.tilt_mode_step_budget = 180; end
if ~isfield(P,'joint_coarse_enable');     P.joint_coarse_enable = true; end
if ~isfield(P,'joint_iter_per_level');    P.joint_iter_per_level = 20; end
if ~isfield(P,'joint_done_e');            P.joint_done_e = 0.08; end
if ~isfield(P,'joint_probe_uv');          P.joint_probe_uv = 1.5; end
if ~isfield(P,'joint_probe_dth');         P.joint_probe_dth = deg2rad(2/3600); end
if ~isfield(P,'joint_gain');              P.joint_gain = 0.6; end
if ~isfield(P,'joint_duv_max');           P.joint_duv_max = 1.0; end
if ~isfield(P,'joint_dth_max');           P.joint_dth_max = deg2rad(3/3600); end
if ~isfield(P,'joint_lam_u');             P.joint_lam_u = 1e-4; end
if ~isfield(P,'joint_lam_v');             P.joint_lam_v = 1e-4; end
if ~isfield(P,'joint_lam_thx');           P.joint_lam_thx = 2e-6; end
if ~isfield(P,'joint_lam_thy');           P.joint_lam_thy = 2e-6; end
if ~isfield(P,'joint_force_soft');        P.joint_force_soft = 0.20; end
if ~isfield(P,'joint_force_hard');        P.joint_force_hard = 0.40; end
if ~isfield(P,'joint_force_predict_scale_min'); P.joint_force_predict_scale_min = 0.10; end

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
if ~isfield(P,'Dc_init')
    P.Dc_init = zeros(5,1);
end
Xtrue.xr = P.xC_init + P.xR_to_C_bias;
Xtrue.Dc = P.Dc_init(:);
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
log_dth   = nan(2, P.N);
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
log_Ftarget = nan(1, P.N);

% fine-search probe cloud (u,v samples from loss_quadfit_uv_step)
probe_u = [];
probe_v = [];
probe_ok = [];

% =============================
% 3) Mode / counters
% =============================
mode = "approach";    % approach -> seat -> coarse -> cont_settle -> cont_fine -> final_attach
N_fit = 0;

level = 1;
F_target = P.force_targets(level);
settle_cnt = 0;
fine_iter  = 0;
tilt_iter  = 0;
joint_iter = 0;
attach_cnt = 0;
tilt_mode_steps = 0;
tilt_probe = init_tilt_probe_state();

% force controller internal
force_ctl.dz_bias = 0;
force_ctl.intF    = 0;

% =============================
% 4) Main loop
% =============================
for k = 1:P.N

    action.duv = [0;0];
    action.dth = [0;0];
    action.dz  = 0;
    action.do_seat = false;

    log_mode(k) = mode;
    log_Ftarget(k) = F_target;

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
            if P.tilt_level_enable
                mode = "tilt_level";
                tilt_iter = 0;
                tilt_mode_steps = 0;
                tilt_probe = init_tilt_probe_state();
            else
                mode = "cont_fine";
                fine_iter = 0;
            end
        end

    elseif mode == "tilt_level"

        % keep force hold active; lock uv and do ONLINE safe probe for theta
        [dz_hold, force_ctl] = force_hold_dz_smooth_zup(state.Fz, F_target, P, force_ctl);
        action.dz = dz_hold;
        action.duv = [0;0];
        tilt_mode_steps = tilt_mode_steps + 1;

        dF_now = abs(state.Fz - F_target);
        if dF_now > P.tilt_force_hard
            % global hard safety in this mode
            action.dth = [0;0];
            action.dz = max(-P.dz_max, min(P.dz_max, action.dz - 0.15));
            tilt_probe = init_tilt_probe_state();
        else
            ax = tilt_probe.axis;
            h = max(1e-12, P.tilt_probe_dth);

            if tilt_probe.phase == "probe_plus_cmd"
                action.dth(ax) = +h;
                tilt_probe.phase = "probe_plus_wait";
                tilt_probe.wait = 0;

            elseif tilt_probe.phase == "probe_plus_wait"
                tilt_probe.wait = tilt_probe.wait + 1;
                if tilt_probe.wait >= P.tilt_probe_settle_steps
                    dF = abs(state.Fz - F_target);
                    if dF <= P.tilt_force_soft
                        tilt_probe.e_plus(:,ax) = e;
                        tilt_probe.valid_plus(ax) = true;
                    elseif dF > P.tilt_force_hard
                        action.dz = max(-P.dz_max, min(P.dz_max, action.dz - 0.15));
                        tilt_probe = init_tilt_probe_state();
                    end
                    if tilt_probe.phase ~= "probe_plus_cmd"
                        tilt_probe.phase = "return_after_plus_cmd";
                    end
                end

            elseif tilt_probe.phase == "return_after_plus_cmd"
                action.dth(ax) = -h;
                tilt_probe.phase = "return_after_plus_wait";
                tilt_probe.wait = 0;

            elseif tilt_probe.phase == "return_after_plus_wait"
                tilt_probe.wait = tilt_probe.wait + 1;
                if tilt_probe.wait >= P.tilt_probe_settle_steps
                    tilt_probe.phase = "probe_minus_cmd";
                end

            elseif tilt_probe.phase == "probe_minus_cmd"
                action.dth(ax) = -h;
                tilt_probe.phase = "probe_minus_wait";
                tilt_probe.wait = 0;

            elseif tilt_probe.phase == "probe_minus_wait"
                tilt_probe.wait = tilt_probe.wait + 1;
                if tilt_probe.wait >= P.tilt_probe_settle_steps
                    dF = abs(state.Fz - F_target);
                    if dF <= P.tilt_force_soft
                        tilt_probe.e_minus(:,ax) = e;
                        tilt_probe.valid_minus(ax) = true;
                    elseif dF > P.tilt_force_hard
                        action.dz = max(-P.dz_max, min(P.dz_max, action.dz - 0.15));
                        tilt_probe = init_tilt_probe_state();
                    end
                    if tilt_probe.phase ~= "probe_plus_cmd"
                        tilt_probe.phase = "return_after_minus_cmd";
                    end
                end

            elseif tilt_probe.phase == "return_after_minus_cmd"
                action.dth(ax) = +h;
                tilt_probe.phase = "return_after_minus_wait";
                tilt_probe.wait = 0;

            elseif tilt_probe.phase == "return_after_minus_wait"
                tilt_probe.wait = tilt_probe.wait + 1;
                if tilt_probe.wait >= P.tilt_probe_settle_steps
                    if tilt_probe.axis < 2
                        tilt_probe.axis = tilt_probe.axis + 1;
                        tilt_probe.phase = "probe_plus_cmd";
                    else
                        tilt_probe.phase = "solve";
                    end
                end

            elseif tilt_probe.phase == "solve"
                J = nan(2,2);
                for jj = 1:2
                    if tilt_probe.valid_plus(jj) && tilt_probe.valid_minus(jj)
                        J(:,jj) = (tilt_probe.e_plus(:,jj) - tilt_probe.e_minus(:,jj)) / (2*h);
                    end
                end
                dth_cmd = tilt_level_solve_from_J(J, e, P);
                action.dth = dth_cmd;
                tilt_iter = tilt_iter + 1;
                tilt_probe = init_tilt_probe_state();
            end
        end

        % exit tilt-level when spot good enough or budget consumed
        if norm(e) <= P.tilt_done_e || tilt_iter >= P.tilt_max_iter_per_level
            if P.joint_coarse_enable
                mode = "joint_coarse";
                joint_iter = 0;
            else
                mode = "cont_fine";
                fine_iter = 0;
            end
        elseif tilt_mode_steps >= P.tilt_mode_step_budget
            if P.joint_coarse_enable
                mode = "joint_coarse";
                joint_iter = 0;
            else
                mode = "cont_fine";
                fine_iter = 0;
            end
        end

    elseif mode == "joint_coarse"

        [dz_hold, force_ctl] = force_hold_dz_smooth_zup(state.Fz, F_target, P, force_ctl);
        action.dz = dz_hold;

        dF_now = abs(state.Fz - F_target);
        if dF_now > P.joint_force_hard
            action.duv = [0;0];
            action.dth = [0;0];
            action.dz = max(-P.dz_max, min(P.dz_max, action.dz - 0.15));
        else
            [duv_cmd, dth_cmd] = joint_coarse_step_numeric(Xtrue.x, e, P);

            % predict-step force risk check (cheap safety layer)
            [duv_cmd, dth_cmd] = joint_force_safe_scale(Xtrue.x, duv_cmd, dth_cmd, F_target, P);

            if dF_now > P.joint_force_soft
                duv_cmd = 0.3 * duv_cmd;
                dth_cmd = 0.3 * dth_cmd;
            end
            action.duv = duv_cmd;
            action.dth = dth_cmd;
            if norm([duv_cmd; dth_cmd]) > 1e-12
                joint_iter = joint_iter + 1;
            end
        end

        if norm(e) <= P.joint_done_e || joint_iter >= P.joint_iter_per_level
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
        [duv_cmd, qfit_info] = loss_quadfit_uv_step(Xtrue, state, P, ...
            P.loss_probe_h_um, P.duv_ls_max, P.e_guard);
        action.duv = duv_cmd;

        % collect sampled probe points for trajectory visualization
        if isstruct(qfit_info) && isfield(qfit_info,'D') && ~isempty(qfit_info.D)
            Dk = qfit_info.D;
            probe_u = [probe_u; Xtrue.x(1) + Dk(:,1)]; %#ok<AGROW>
            probe_v = [probe_v; Xtrue.x(2) + Dk(:,2)]; %#ok<AGROW>
            if isfield(qfit_info,'Eok') && numel(qfit_info.Eok)==size(Dk,1)
                probe_ok = [probe_ok; logical(qfit_info.Eok(:))]; %#ok<AGROW>
            else
                probe_ok = [probe_ok; true(size(Dk,1),1)]; %#ok<AGROW>
            end
        end

        % ---- record discrete fitted loss point (真实扫频点) ----
        [Lnow, oknow, infonow] = sensor_lossfit(Xtrue.x, e, state.Keff, state.seated, P);
        N_fit = N_fit + 1;
        log_Lfit(k)  = Lnow;
        log_ok(k)    = oknow;
        log_psucc(k) = infonow.p_succ;

        reached_loss = isfinite(Lnow) && oknow && (Lnow < P.L_thresh_ppm);
        reached_process_level = (level >= P.min_level_before_exit);
        if reached_loss
            if reached_process_level
                if P.final_attach_enable
                    mode = "final_attach";
                    attach_cnt = 0;
                    fprintf("Loss met at level %d (step %d). Entering final_attach stage.\n", level, k);
                else
                    fprintf("Reached target loss < %.0f ppm at step %d (level %d), Lfit=%.1f ppm, F=%.2fN, N_fit=%d\n", ...
                        P.L_thresh_ppm, k, level, Lnow, state.Fz, N_fit);
                    break;
                end
            else
                fprintf("Loss met at level %d but hold exit until min level %d.\n", level, P.min_level_before_exit);
            end
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

    elseif mode == "final_attach"

        [dz_hold, force_ctl] = force_hold_dz_smooth_zup(state.Fz, F_target, P, force_ctl);
        dz_extra = 0;
        if attach_cnt < P.attach_push_steps
            dz_extra = P.dz_attach_pulse;
        end
        action.dz = max(-P.dz_max, min(P.dz_max, dz_hold + dz_extra));
        action.duv = [0;0];

        [Lnow, oknow, ~] = sensor_lossfit(Xtrue.x, e, state.Keff, state.seated, P);
        N_fit = N_fit + 1;
        log_Lfit(k) = Lnow;
        log_ok(k)   = oknow;

        attach_cnt = attach_cnt + 1;
        attach_done = (attach_cnt >= (P.attach_push_steps + P.attach_hold_steps));
        attach_force_ok = abs(state.Fz - F_target) <= P.attach_force_tol;
        attach_loss_ok = isfinite(Lnow) && oknow && (Lnow < P.L_thresh_ppm);
        if attach_done && attach_force_ok && attach_loss_ok
            fprintf("Final attach complete at step %d. Lfit=%.1f ppm, F=%.2fN, N_fit=%d\n", ...
                k, Lnow, state.Fz, N_fit);
            break;
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
    log_dth(:,k) = action.dth(:);
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
% (B) Build commonly used signals
% =============================
X  = log_x(:,1:Kend);
Xr = log_xr(:,1:Kend);

u   = X(1,:);   v   = X(2,:);
z   = X(3,:);   thx = X(4,:);  thy = X(5,:);

ur  = Xr(1,:);  vr  = Xr(2,:);

e_y = log_e(1,1:Kend);
e_z = log_e(2,1:Kend);
e_n = vecnorm(log_e(:,1:Kend),2,1);

th_n  = sqrt(thx.^2 + thy.^2);
dth_n = sqrt(log_dth(1,1:Kend).^2 + log_dth(2,1:Kend).^2);
duv_n = sqrt(log_duv(1,1:Kend).^2 + log_duv(2,1:Kend).^2);

Dev_cmd = X - Xr;
dev_uv_cmd = vecnorm(Dev_cmd(1:2,:),2,1);
dev_th_cmd = vecnorm(Dev_cmd(4:5,:),2,1);

if isfield(P,'Ro_um') && isfield(P,'Ri_um')
    A_ring = pi * (P.Ro_um^2 - P.Ri_um^2);
    Ac_ratio = log_Ac(1:Kend) / A_ring;
else
    Ac_ratio = nan(1,Kend);
end

F_target_trace = log_Ftarget(1:Kend);

dec = max(1, floor(Kend/800));
idx = 1:dec:Kend;

close all;

% =============================
% FIGURE 1: Force / target / contact
% =============================
figure('Name','Fig1_Force_Target_Contact','Color','w');
tiledlayout(4,1,'Padding','compact','TileSpacing','compact');

nexttile; hold on; grid on;
shade_modes(gca, log_mode(1:Kend));
plot(idx, log_F(idx), 'LineWidth',1.2);
plot(idx, F_target_trace(idx), '--', 'LineWidth',1.1);
ylabel('F_z (N)');
title('Force tracking');
legend({'F_z','F_{target}'}, 'Location','best');

nexttile; hold on; grid on;
shade_modes(gca, log_mode(1:Kend));
plot(idx, log_K(idx), 'LineWidth',1.2);
ylabel('K_{eff}');
title('Effective contact stiffness');

nexttile; hold on; grid on;
shade_modes(gca, log_mode(1:Kend));
if all(~isnan(Ac_ratio))
    plot(idx, Ac_ratio(idx), 'LineWidth',1.2);
    ylabel('A_c/A_{ring}');
else
    plot(idx, log_Ac(idx), 'LineWidth',1.2);
    ylabel('A_c');
end
title('Contact area ratio');

nexttile; hold on; grid on;
shade_modes(gca, log_mode(1:Kend));
plot(idx, log_g(idx), 'LineWidth',1.2);
if any(isfinite(log_ggeo(1:Kend)))
    plot(idx, log_ggeo(idx), 'LineWidth',1.1);
    legend({'g_{signed}','g_{geo}'}, 'Location','best');
end
yline(0,'--');
xlabel('step');
ylabel('gap (\mum)');
title('Gap evolution');

% =============================
% FIGURE 2: Spot + angle + action
% =============================
figure('Name','Fig2_Spot_Angle','Color','w');
tiledlayout(3,1,'Padding','compact','TileSpacing','compact');

nexttile; hold on; grid on;
shade_modes(gca, log_mode(1:Kend));
plot(idx, e_y(idx), 'LineWidth',1.2);
plot(idx, e_z(idx), 'LineWidth',1.2);
plot(idx, e_n(idx), 'LineWidth',1.3);
ylabel('e (mm)');
title('Spot centroid error');
legend({'e_y','e_z','||e||'}, 'Location','best');

nexttile; hold on; grid on;
shade_modes(gca, log_mode(1:Kend));
plot(idx, thx(idx), 'LineWidth',1.2);
plot(idx, thy(idx), 'LineWidth',1.2);
plot(idx, th_n(idx), 'LineWidth',1.3);
ylabel('\theta (rad)');
title('Angle evolution');
legend({'\theta_x','\theta_y','||\theta||'}, 'Location','best');

nexttile; hold on; grid on;
shade_modes(gca, log_mode(1:Kend));
yyaxis left; plot(idx, dth_n(idx), 'LineWidth',1.2); ylabel('||d\theta|| (rad)');
yyaxis right; plot(idx, duv_n(idx), 'LineWidth',1.2); ylabel('||duv|| (\mum)');
xlabel('step');
title('Action magnitudes');

% =============================
% FIGURE 3: Loss curve
% =============================
figure('Name','Fig3_Loss','Color','w');
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

nexttile; hold on; grid on;
shade_modes(gca, log_mode(1:Kend));
plot(idx, log_Ltrue(idx), 'LineWidth',1.2);
fit_idx = find(isfinite(log_Lfit(1:Kend)));
if ~isempty(fit_idx)
    scatter(fit_idx, log_Lfit(fit_idx), 10, 'filled');
end
yline(P.L_thresh_ppm, '--', sprintf('%.0f ppm', P.L_thresh_ppm));
ylabel('Loss (ppm)');
title('Loss evolution');
legend({'L_{true}','L_{fit}','threshold'}, 'Location','best');

nexttile; hold on; grid on;
shade_modes(gca, log_mode(1:Kend));
yyaxis left;
plot(idx, e_n(idx), 'LineWidth',1.2);
ylabel('||e|| (mm)');
yyaxis right;
plot(idx, th_n(idx), 'LineWidth',1.2);
ylabel('||\theta|| (rad)');
xlabel('step');
title('Loss bottleneck indicators');

% =============================
% FIGURE 4: Command vs true mismatch
% =============================
figure('Name','Fig4_Command_True_Mismatch','Color','w');
tiledlayout(3,1,'Padding','compact','TileSpacing','compact');

nexttile; hold on; grid on;
shade_modes(gca, log_mode(1:Kend));
plot(idx, dev_uv_cmd(idx), 'LineWidth',1.2);
ylabel('||\Delta uv|| (\mum)');
title('UV mismatch: true - commanded');

nexttile; hold on; grid on;
shade_modes(gca, log_mode(1:Kend));
plot(idx, dev_th_cmd(idx), 'LineWidth',1.2);
ylabel('||\Delta \theta|| (rad)');
title('Angle mismatch: true - commanded');

nexttile; hold on; grid on;
shade_modes(gca, log_mode(1:Kend));
plot(idx, Dev_cmd(3,idx), 'LineWidth',1.2);
xlabel('step');
ylabel('\Delta z (\mum)');
title('Z mismatch: true - commanded');

% =============================
% FIGURE 5: Contact quality
% =============================
figure('Name','Fig5_Contact_Quality','Color','w');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

nexttile; hold on; grid on;
shade_modes(gca, log_mode(1:Kend));
if all(~isnan(Ac_ratio))
    plot(idx, Ac_ratio(idx), 'LineWidth',1.2);
    ylabel('A_c/A_{ring}');
else
    plot(idx, log_Ac(idx), 'LineWidth',1.2);
    ylabel('A_c');
end
xlabel('step');
title('Contact area ratio');

nexttile; hold on; grid on;
shade_modes(gca, log_mode(1:Kend));
plot(idx, log_K(idx), 'LineWidth',1.2);
xlabel('step');
ylabel('K_{eff}');
title('Contact stiffness');

nexttile; hold on; grid on;
shade_modes(gca, log_mode(1:Kend));
plot(idx, hypot(log_Fx(idx), log_Fy(idx)), 'LineWidth',1.2);
xlabel('step');
ylabel('||F_t|| (N)');
title('Tangential friction magnitude');

nexttile; hold on; grid on;
shade_modes(gca, log_mode(1:Kend));
plot(idx, log_stick(idx), 'LineWidth',1.2);
plot(idx, log_slip(idx),  'LineWidth',1.2);
xlabel('step');
ylabel('ratio');
title('Stick / slip ratio');
legend({'stick','slip'}, 'Location','best');

% =============================
% FIGURE 6: UV trajectory + fine probe cloud
% =============================
figure('Name','Fig6_UV_Trajectory','Color','w');
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

nexttile; hold on; grid on; axis equal;
if ~isempty(probe_u)
    plot(probe_u, probe_v, '.', 'Color', [0.82 0.82 0.82], 'MarkerSize', 6);
    if any(probe_ok)
        plot(probe_u(probe_ok), probe_v(probe_ok), '.', 'Color', [0.50 0.70 1.00], 'MarkerSize', 7);
    end
end
plot(ur, vr, 'LineWidth',1.2);
plot(u,  v,  'LineWidth',1.2);
xlabel('u (\mum)'); ylabel('v (\mum)');
title('UV trajectory');
if ~isempty(probe_u)
    legend({'fine probe(all)','fine probe(valid)','commanded','true'}, 'Location','best');
else
    legend({'commanded','true'}, 'Location','best');
end

nexttile; hold on; grid on; axis equal;
scatter(u, v, 12, e_n, 'filled');
xlabel('u (\mum)'); ylabel('v (\mum)');
title('UV trajectory colored by ||e||');
cb = colorbar; ylabel(cb,'||e|| (mm)');

% =============================
% FIGURE 7: tilt_level zoom (contiguous segments)
% =============================
ranges_tilt = mode_contiguous_ranges(log_mode(1:Kend), "tilt_level");
if ~isempty(ranges_tilt)
    figure('Name','Fig7_TiltLevel_Zoom','Color','w');
    tiledlayout(4,1,'Padding','compact','TileSpacing','compact');
    for rr = 1:size(ranges_tilt,1)
        ks = ranges_tilt(rr,1):ranges_tilt(rr,2);
        nexttile(1); hold on; grid on; plot(ks, th_n(ks), 'LineWidth',1.2);
        nexttile(2); hold on; grid on; plot(ks, dth_n(ks), 'LineWidth',1.2);
        nexttile(3); hold on; grid on; plot(ks, e_n(ks), 'LineWidth',1.2);
        nexttile(4); hold on; grid on; plot(ks, log_F(ks)-F_target_trace(ks), 'LineWidth',1.2);
    end
    nexttile(1); ylabel('||\theta||'); title('tilt\_level: angle norm');
    nexttile(2); ylabel('||d\theta||'); title('tilt\_level: angle command norm');
    nexttile(3); ylabel('||e||'); title('tilt\_level: spot error norm');
    nexttile(4); yline(0,'--'); xlabel('step'); ylabel('F_z-F_{target}'); title('tilt\_level: force deviation');
end

% =============================
% FIGURE 8: joint_coarse zoom (contiguous segments)
% =============================
ranges_joint = mode_contiguous_ranges(log_mode(1:Kend), "joint_coarse");
if ~isempty(ranges_joint)
    figure('Name','Fig8_JointCoarse_Zoom','Color','w');
    tiledlayout(4,1,'Padding','compact','TileSpacing','compact');
    for rr = 1:size(ranges_joint,1)
        ks = ranges_joint(rr,1):ranges_joint(rr,2);
        nexttile(1); hold on; grid on; plot(ks, duv_n(ks), 'LineWidth',1.2);
        nexttile(2); hold on; grid on; plot(ks, dth_n(ks), 'LineWidth',1.2);
        nexttile(3); hold on; grid on; plot(ks, e_n(ks), 'LineWidth',1.2);
        nexttile(4); hold on; grid on; plot(ks, log_F(ks)-F_target_trace(ks), 'LineWidth',1.2);
    end
    nexttile(1); ylabel('||duv||'); title('joint\_coarse: uv command norm');
    nexttile(2); ylabel('||d\theta||'); title('joint\_coarse: angle command norm');
    nexttile(3); ylabel('||e||'); title('joint\_coarse: spot error norm');
    nexttile(4); yline(0,'--'); xlabel('step'); ylabel('F_z-F_{target}'); title('joint\_coarse: force deviation');
end

% =============================
% FIGURE 9: mode-by-mode improvement summary (segmented)
% =============================
figure('Name','Fig9_Mode_Improvement','Color','w');
mode_list = ["coarse","tilt_level","joint_coarse","cont_fine"];
[delta_e, delta_th, delta_L] = mode_improvement_stats_segmented(log_mode(1:Kend), e_n, th_n, log_Ltrue(1:Kend), mode_list);

tiledlayout(3,1,'Padding','compact','TileSpacing','compact');
nexttile; bar(categorical(cellstr(mode_list)), delta_e); grid on; ylabel('\Delta ||e||'); title('Mode-wise avg. spot improvement');
nexttile; bar(categorical(cellstr(mode_list)), delta_th); grid on; ylabel('\Delta ||\theta||'); title('Mode-wise avg. angle improvement');
nexttile; bar(categorical(cellstr(mode_list)), delta_L); grid on; ylabel('\Delta L_{true} (ppm)'); title('Mode-wise avg. loss improvement');

% =============================
% FIGURE 10: mechanism scatter plots
% =============================
figure('Name','Fig10_Mechanism_Scatters','Color','w');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

nexttile; hold on; grid on;
scatter(th_n, e_n, 12, 1:Kend, 'filled');
xlabel('||\theta|| (rad)'); ylabel('||e|| (mm)'); title('\theta \rightarrow e'); cb = colorbar; ylabel(cb,'step');

nexttile; hold on; grid on;
scatter(e_n, log_Ltrue(1:Kend), 12, 1:Kend, 'filled');
xlabel('||e|| (mm)'); ylabel('L_{true} (ppm)'); title('e \rightarrow loss'); cb = colorbar; ylabel(cb,'step');

nexttile; hold on; grid on;
scatter(th_n, log_Ltrue(1:Kend), 12, 1:Kend, 'filled');
xlabel('||\theta|| (rad)'); ylabel('L_{true} (ppm)'); title('\theta \rightarrow loss'); cb = colorbar; ylabel(cb,'step');

nexttile; hold on; grid on;
scatter(dth_n, e_n, 12, 1:Kend, 'filled');
xlabel('||d\theta|| (rad)'); ylabel('||e|| (mm)'); title('action vs residual spot'); cb = colorbar; ylabel(cb,'step');

end % ===== end main_sim =====


% ==========================================================
% Plot helpers
% ==========================================================
function shade_modes(ax, mode_vec)
hold(ax,'on');
K = numel(mode_vec);
if K < 1, return; end

mode_colors = containers.Map;
mode_colors('approach')     = [0.95 0.95 0.95];
mode_colors('seat')         = [0.90 0.95 1.00];
mode_colors('coarse')       = [1.00 0.96 0.88];
mode_colors('cont_settle')  = [0.92 1.00 0.92];
mode_colors('tilt_level')   = [0.90 0.95 1.00];
mode_colors('joint_coarse') = [0.96 0.90 1.00];
mode_colors('cont_fine')    = [1.00 0.92 0.92];
mode_colors('final_attach') = [0.95 0.90 0.90];

yl = ylim(ax);
s = 1;
while s <= K
    m = mode_vec(s);
    t = s;
    while t <= K && mode_vec(t) == m
        t = t + 1;
    end
    x0 = s;
    x1 = t - 1;
    key = char(m);
    if isKey(mode_colors, key)
        c = mode_colors(key);
    else
        c = [0.95 0.95 0.95];
    end
    patch(ax, [x0 x1 x1 x0], [yl(1) yl(1) yl(2) yl(2)], c, ...
        'FaceAlpha', 0.10, 'EdgeColor','none', 'HandleVisibility','off');
    s = t;
end
ylim(ax, yl);
end

function ranges = mode_contiguous_ranges(mode_vec, target_mode)
idx = find(mode_vec == target_mode);
if isempty(idx)
    ranges = zeros(0,2);
    return;
end
breaks = find(diff(idx) > 1);
starts = idx([1, breaks + 1]);
ends = idx([breaks, numel(idx)]);
ranges = [starts(:), ends(:)];
end

function [delta_e, delta_th, delta_L] = mode_improvement_stats_segmented(mode_vec, e_n, th_n, L_true, mode_list)
nm = numel(mode_list);
delta_e  = zeros(1,nm);
delta_th = zeros(1,nm);
delta_L  = zeros(1,nm);
for i = 1:nm
    ranges = mode_contiguous_ranges(mode_vec, mode_list(i));
    if isempty(ranges), continue; end
    de = nan(size(ranges,1),1);
    dth = nan(size(ranges,1),1);
    dL = nan(size(ranges,1),1);
    for r = 1:size(ranges,1)
        i0 = ranges(r,1); i1 = ranges(r,2);
        de(r) = e_n(i0) - e_n(i1);
        dth(r) = th_n(i0) - th_n(i1);
        dL(r) = L_true(i0) - L_true(i1);
    end
    delta_e(i) = mean(de,'omitnan');
    delta_th(i) = mean(dth,'omitnan');
    delta_L(i) = mean(dL,'omitnan');
end
end


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


% ==========================================================
% Local helpers: online tip/tilt safe probing (2x2)
% ==========================================================
function S = init_tilt_probe_state()
S.axis = 1;
S.phase = "probe_plus_cmd";
S.wait = 0;
S.e_plus = nan(2,2);
S.e_minus = nan(2,2);
S.valid_plus = false(1,2);
S.valid_minus = false(1,2);
end

function dth_cmd = tilt_level_solve_from_J(J, e, P)
if ~isfield(P,'tilt_step_max'); P.tilt_step_max = deg2rad(4/3600); end
if ~isfield(P,'tilt_gain');     P.tilt_gain = 0.5; end
if ~isfield(P,'tilt_lambda');   P.tilt_lambda = 1e-6; end
if ~isfield(P,'tilt_rcond_min'); P.tilt_rcond_min = 1e-10; end

if any(~isfinite(J(:)))
    dth_cmd = [0;0];
    return;
end

A = (J.'*J + P.tilt_lambda*eye(2));
if rcond(A) < P.tilt_rcond_min
    dth_cmd = [0;0];
    return;
end
b = (J.'*e);
dth = -(A \ b);
dth = P.tilt_gain * dth;

dth(1) = max(-P.tilt_step_max, min(P.tilt_step_max, dth(1)));
dth(2) = max(-P.tilt_step_max, min(P.tilt_step_max, dth(2)));

if any(~isfinite(dth))
    dth = [0;0];
end
dth_cmd = dth;
end


% ==========================================================
% Local helper: joint coarse correction (2x4 on [u v thx thy])
% ==========================================================
function [duv_cmd, dth_cmd] = joint_coarse_step_numeric(x, e, P)
if ~isfield(P,'joint_probe_uv');  P.joint_probe_uv = 1.5; end
if ~isfield(P,'joint_probe_dth'); P.joint_probe_dth = deg2rad(2/3600); end
if ~isfield(P,'joint_gain');      P.joint_gain = 0.6; end
if ~isfield(P,'joint_duv_max');   P.joint_duv_max = 1.0; end
if ~isfield(P,'joint_dth_max');   P.joint_dth_max = deg2rad(3/3600); end
if ~isfield(P,'joint_lam_u');     P.joint_lam_u = 1e-4; end
if ~isfield(P,'joint_lam_v');     P.joint_lam_v = 1e-4; end
if ~isfield(P,'joint_lam_thx');   P.joint_lam_thx = 2e-6; end
if ~isfield(P,'joint_lam_thy');   P.joint_lam_thy = 2e-6; end

du = max(1e-9, P.joint_probe_uv);
dv = max(1e-9, P.joint_probe_uv);
dtx = max(1e-12, P.joint_probe_dth);
dty = max(1e-12, P.joint_probe_dth);

e0 = sensor_vision(x, P);

xu = x; xu(1) = xu(1) + du;  eu = sensor_vision(xu, P);
xv = x; xv(2) = xv(2) + dv;  ev = sensor_vision(xv, P);
xtx = x; xtx(4) = xtx(4) + dtx; etx = sensor_vision(xtx, P);
xty = x; xty(5) = xty(5) + dty; ety = sensor_vision(xty, P);

J = [ (eu-e0)/du, (ev-e0)/dv, (etx-e0)/dtx, (ety-e0)/dty ];  % 2x4
R = diag([P.joint_lam_u, P.joint_lam_v, P.joint_lam_thx, P.joint_lam_thy]);

dq = -((J.'*J + R) \ (J.'*e));
dq = P.joint_gain * dq;

if any(~isfinite(dq)) || norm(J,'fro') < 1e-12
    dq = zeros(4,1);
end

dq(1) = max(-P.joint_duv_max, min(P.joint_duv_max, dq(1)));
dq(2) = max(-P.joint_duv_max, min(P.joint_duv_max, dq(2)));
dq(3) = max(-P.joint_dth_max, min(P.joint_dth_max, dq(3)));
dq(4) = max(-P.joint_dth_max, min(P.joint_dth_max, dq(4)));

duv_cmd = dq(1:2);
dth_cmd = dq(3:4);
end


function [duv_safe, dth_safe] = joint_force_safe_scale(x, duv_cmd, dth_cmd, F_target, P)
if ~isfield(P,'joint_force_hard'); P.joint_force_hard = 0.40; end
if ~isfield(P,'joint_force_predict_scale_min'); P.joint_force_predict_scale_min = 0.10; end

scales = [1.0, 0.6, 0.35, 0.2, P.joint_force_predict_scale_min];
duv_safe = [0;0];
dth_safe = [0;0];

for s = scales
    duv_try = s * duv_cmd;
    dth_try = s * dth_cmd;
    x_try = x;
    x_try(1:2) = x_try(1:2) + duv_try(:);
    x_try(4:5) = x_try(4:5) + dth_try(:);

    [F_try, ~, ~] = force_model(x_try, P);
    if abs(F_try - F_target) <= P.joint_force_hard
        duv_safe = duv_try;
        dth_safe = dth_try;
        return;
    end
end
end
