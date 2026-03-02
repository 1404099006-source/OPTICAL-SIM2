function [Xtrue, state] = truth_update(Xtrue, state, action, P)
%TRUTH_UPDATE Update true lens pose with quasi-static contact equilibrium.
%
% State x = [u v z thx thy] in cavity frame.
% Robot commands u,v,z through xr; z/theta are solved by static balance under contact.

    % ---------------------------
    % 0) Defaults / safety
    % ---------------------------
    if ~isfield(P,'alpha_def');   P.alpha_def = 0.25; end
    if ~isfield(P,'c_z');         P.c_z = 2.0; end
    if ~isfield(P,'c_th');        P.c_th = 5e-8; end
    if ~isfield(P,'theta_damp');  P.theta_damp = 0.2; end
    if ~isfield(P,'dtheta_max');  P.dtheta_max = inf; end

    % stiffness form (preferred); fallback from compliance
    if ~isfield(P,'k_z') || P.k_z <= 0
        P.k_z = 1/max(P.c_z, 1e-12); % N/um
    end
    if ~isfield(P,'k_thx') || P.k_thx <= 0
        P.k_thx = 1/max(P.c_th, 1e-18); % N*um/rad
    end
    if ~isfield(P,'k_thy') || P.k_thy <= 0
        P.k_thy = 1/max(P.c_th, 1e-18);
    end
    if ~isfield(P,'eq_max_iter'); P.eq_max_iter = 20; end
    if ~isfield(P,'eq_relax');    P.eq_relax = [0.6; 0.6; 0.6]; end
    if ~isfield(P,'eq_tol');      P.eq_tol = [1e-4; 1e-7; 1e-7]; end

    % friction defaults
    if ~isfield(P,'mu');              P.mu = 0.2; end
    if ~isfield(P,'k_t');             P.k_t = 0.3 * P.k_w; end
    if ~isfield(P,'g0_um');            P.g0_um = 0; end
    if ~isfield(P,'alpha_fx2fz');     P.alpha_fx2fz = 0.0; end
    if ~isfield(P,'alpha_fy2fz');     P.alpha_fy2fz = 0.0; end
    if ~isfield(P,'friction_mode');   P.friction_mode = "distributed"; end

    % ---------------------------
    % 1) Apply robot command to nominal pose xr (in cavity {C})
    % ---------------------------
    duv = action.duv(:);
    dz = action.dz;

    if isfield(P,'act_enable') && P.act_enable
        if ~isfield(P,'act_bias_uv');  P.act_bias_uv = [0; 0]; end
        if ~isfield(P,'act_bias_z');   P.act_bias_z = 0; end
        if ~isfield(P,'act_noise_uv'); P.act_noise_uv = 0; end
        if ~isfield(P,'act_noise_z');  P.act_noise_z = 0; end

        duv = duv + P.act_bias_uv(:) + P.act_noise_uv * randn(2,1);
        dz = dz + P.act_bias_z + P.act_noise_z * randn();
    end

    Xtrue.xr(1:2) = Xtrue.xr(1:2) + duv;
    Xtrue.xr(3)   = Xtrue.xr(3)   + dz;

    % ---------------------------
    % 2) Contact pre-check
    % ---------------------------
    Xtmp = Xtrue.xr + Xtrue.Dc;
    [Fz_tmp, ~, ~, Keff_c, ~, ~] = force_model(Xtmp, P);

    % ---------------------------
    % 3) Seating flag
    % ---------------------------
    if action.do_seat
        state.seated = true;
    end

    % ---------------------------
    % 4) Quasi-static equilibrium for [z thx thy] under contact
    % ---------------------------
    Dc_target = Xtrue.Dc;

    if state.seated && (Fz_tmp > 0)
        % robot uv induces weak extra nominal tilt through gripper compliance
        if ~isfield(P,'K_uv_to_th'); P.K_uv_to_th = zeros(2,2); end
        th_uv = P.K_uv_to_th * Xtrue.xr(1:2);

        z_r = Xtrue.xr(3);
        th_r = Xtrue.xr(4:5) + th_uv;

        q = [Xtmp(3); Xtmp(4); Xtmp(5)];
        for it = 1:P.eq_max_iter
            x_eval = [Xtmp(1); Xtmp(2); q(1); q(2); q(3)];
            [Fz_i, Mx_i, My_i] = force_model(x_eval, P);

            r = [ Fz_i - P.k_z   * (q(1)-z_r);
                  Mx_i - P.k_thx * (q(2)-th_r(1));
                  My_i - P.k_thy * (q(3)-th_r(2)) ];

            if all(abs(r) <= P.eq_tol)
                break;
            end

            % r = F_contact - K*(q-q_ref), so use +r/K update toward equilibrium
            dq = [ r(1)/max(P.k_z,1e-12);
                   r(2)/max(P.k_thx,1e-12);
                   r(3)/max(P.k_thy,1e-12) ];

            dq = P.eq_relax(:) .* dq;
            dq(2:3) = max(min(dq(2:3), P.dtheta_max), -P.dtheta_max);

            q = q + dq;
        end

        % convert equilibrium to deformation relative to robot nominal xr
        Dc_target(3)   = q(1) - Xtrue.xr(3);
        Dc_target(4:5) = q(2:3) - Xtrue.xr(4:5);
    else
        Dc_target(3)   = 0;
        Dc_target(4:5) = Xtrue.Dc(4:5);
    end

    % smooth deformation update
    alpha = min(max(P.alpha_def, 0.0), 1.0);
    Xtrue.Dc(3:5) = (1-alpha)*Xtrue.Dc(3:5) + alpha*Dc_target(3:5);

    % ---------------------------
    % 5) Drift (optional)
    % ---------------------------
    if exist('contact_drift_update','file') == 2
        Xtrue.Dc = contact_drift_update(Xtrue.Dc, Fz_tmp, Keff_c, state.seated, P);
    end

    % ---------------------------
    % 6) Final true pose and contact
    % ---------------------------
    Xtrue.x = Xtrue.xr + Xtrue.Dc;
    [Fz3, Mx3, My3, Keff3, g_geo3, g_signed3, Ac3] = force_model(Xtrue.x, P);

    % ---------------------------
    % 7) Friction (stick-slip)
    %   - distributed(default): per-node partial slip
    %   - lumped: legacy whole-contact stick/slip
    % ---------------------------
    Ft = [0;0];
    Mz = 0;
    stick_ratio = 1.0;
    slip_ratio = 0.0;

    if P.friction_mode == "distributed"
        if isfield(P,'contact_x') && isfield(P,'contact_y')
            Xc = P.contact_x(:).';
            Yc = P.contact_y(:).';
            dA = P.contact_dA;
        else
            Xc = P.disk_x(:).';
            Yc = P.disk_y(:).';
            dA = P.disk_dA;
        end

        Nn = numel(Xc);
        if ~isfield(state,'s_t') || ~isequal(size(state.s_t), [2, Nn])
            state.s_t = zeros(2, Nn); % per-node tangential stick displacement [um]
        end

        if Ac3 <= 0 || Fz3 <= 0 || Nn == 0
            state.s_t = zeros(2, Nn);
        else
            if ~isfield(P,'K_uv_to_th'); P.K_uv_to_th = zeros(2,2); end
            th_uv = P.K_uv_to_th * Xtrue.x(1:2);
            thx_eff = Xtrue.x(4) + th_uv(1);
            thy_eff = Xtrue.x(5) + th_uv(2);

            g = P.g0_um - Xtrue.x(3) + thx_eff .* Yc - thy_eff .* Xc;
            w = max(-g, 0);
            in_contact = (w > 0);
            p_i = P.k_w .* w; % [N/um^2]

            duv_mat = repmat(duv(:), 1, Nn);
            s_trial = state.s_t + duv_mat;
            t_trial = P.k_t .* s_trial; % [N/um^2]

            n_trial = sqrt(sum(t_trial.^2, 1));
            t_lim = P.mu .* p_i;

            stick_mask = (n_trial <= t_lim) | (n_trial < 1e-15) | (~in_contact);
            slip_mask = in_contact & ~stick_mask;

            t = zeros(2, Nn);
            t(:, stick_mask) = t_trial(:, stick_mask);
            if any(slip_mask)
                scale = t_lim(slip_mask) ./ max(n_trial(slip_mask), 1e-15);
                t(:, slip_mask) = t_trial(:, slip_mask) .* scale;
            end
            t(:, ~in_contact) = 0;

            Ft = sum(t, 2) .* dA;
            Mz = sum((Xc .* t(2,:) - Yc .* t(1,:)) .* dA);

            if P.k_t > 0
                state.s_t = t ./ P.k_t;
            else
                state.s_t = zeros(2, Nn);
            end

            n_contact = nnz(in_contact);
            if n_contact > 0
                stick_ratio = nnz(stick_mask & in_contact) / n_contact;
                slip_ratio = nnz(slip_mask) / n_contact;
            end
        end
    else
        % legacy lumped model
        if ~isfield(state,'s_t') || numel(state.s_t)~=2
            state.s_t = zeros(2,1); % um, equivalent tangential stick displacement
        end

        if Ac3 <= 0 || Fz3 <= 0
            state.s_t = zeros(2,1);
        else
            s_trial = state.s_t + duv;
            Ft_trial = (P.k_t * Ac3) * s_trial; % N

            Ft_lim = P.mu * Fz3;
            nFt = norm(Ft_trial);
            if nFt <= Ft_lim || nFt < 1e-15
                Ft = Ft_trial;               % stick
                state.s_t = s_trial;
                stick_ratio = 1.0;
                slip_ratio = 0.0;
            else
                Ft = Ft_lim * (Ft_trial / nFt); % slip saturation
                state.s_t = Ft / max(P.k_t * Ac3, 1e-15);
                stick_ratio = 0.0;
                slip_ratio = 1.0;
            end
        end
    end

    % ---------------------------
    % 8) Update state outputs
    % ---------------------------
    state.F_prev = Fz3;
    state.Keff   = Keff3;
    state.gmin   = g_signed3;
    state.g_geo  = g_geo3;

    state.Fz     = Fz3;
    state.Mx     = Mx3;
    state.My     = My3;
    state.Ac     = Ac3;

    state.Fx     = Ft(1);
    state.Fy     = Ft(2);
    state.Mz     = Mz;
    state.stick_ratio = stick_ratio;
    state.slip_ratio  = slip_ratio;
    state.Fz_meas = Fz3 + P.alpha_fx2fz * Ft(1) + P.alpha_fy2fz * Ft(2);
end
