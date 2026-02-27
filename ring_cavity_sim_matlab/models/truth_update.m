function [Xtrue, state] = truth_update(Xtrue, state, action, P)
%TRUTH_UPDATE Update "true" lens pose relative to cavity {C}
%
% Xtrue.xr : commanded/nominal pose xC (5x1) in cavity frame {C}
% Xtrue.Dc : deformation/drift (5x1)
% Xtrue.x  : true pose xC = xr + Dc
%
% action.duv: [du;dv] (um)  -- robot in-plane move (assume aligned with {C})
% action.dz : scalar (um)   -- robot z move (z-up positive)
% action.do_seat : bool

    % ---------------------------
    % 0) Defaults / safety
    % ---------------------------
    if ~isfield(P,'alpha_def'); P.alpha_def = 0.25; end
    if ~isfield(P,'c_z');       P.c_z = 0.0; end
    if ~isfield(P,'c_th');         P.c_th = 0.0; end
    if ~isfield(P,'theta_damp');   P.theta_damp = 0.0; end
    if ~isfield(P,'dtheta_max');   P.dtheta_max = inf; end

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

    Xtrue.xr(1:2) = Xtrue.xr(1:2) + duv; % u,v
    Xtrue.xr(3)   = Xtrue.xr(3)   + dz;  % z (z-up positive)
    % angles not directly actuated in this version

    % ---------------------------
    % 2) Temporary true pose
    % ---------------------------
    Xtmp = Xtrue.xr + Xtrue.Dc;

    % ---------------------------
    % 3) Contact result at temporary pose
    % ---------------------------
    [Fz, Mx, My, Keff_c, ~, ~] = force_model(Xtmp, P);

    % ---------------------------
    % 4) Seating flag
    % ---------------------------
    if action.do_seat
        state.seated = true;
    end

    % ---------------------------
    % 5) Force-induced deformation (z + tilt from contact moments)
    % ---------------------------
    Dc_target = Xtrue.Dc;

    if state.seated && (Fz > 0)
        % axial compliance under contact force
        Dc_target(3) = +P.c_z * Fz;

        % moment-driven tilt equilibrium (Scheme A)
        % Use restoring-moment direction directly; wrong sign here will diverge tilt.
        th_eq = [P.c_th * Mx; P.c_th * My];

        % optional damping toward equilibrium to avoid chattering
        theta_damp = min(max(P.theta_damp, 0.0), 1.0);
        th_prev = Xtrue.Dc(4:5);
        th_target = (1-theta_damp) * th_eq + theta_damp * th_prev;

        % per-step angle increment cap for numerical stability
        dth = th_target - th_prev;
        dth = max(min(dth, P.dtheta_max), -P.dtheta_max);
        Dc_target(4:5) = th_prev + dth;
    else
        Dc_target(3) = 0;
        Dc_target(4:5) = Xtrue.Dc(4:5);
    end

    alpha = min(max(P.alpha_def, 0.0), 1.0);
    Xtrue.Dc(3:5) = (1-alpha)*Xtrue.Dc(3:5) + alpha*Dc_target(3:5);

    % ---------------------------
    % 6) Drift (optional)
    % ---------------------------
    Keff = Keff_c;
    if exist('contact_drift_update','file') == 2
        Xtrue.Dc = contact_drift_update(Xtrue.Dc, Fz, Keff, state.seated, P);
    end

    % ---------------------------
    % 7) Final true pose
    % ---------------------------
    Xtrue.x = Xtrue.xr + Xtrue.Dc;

    % ---------------------------
    % 8) Recompute contact at final pose for consistent outputs
    % ---------------------------
    [Fz3, Mx3, My3, Keff3, g_geo3, g_signed3] = force_model(Xtrue.x, P);

    % ---------------------------
    % 9) Update state
    % ---------------------------
    state.F_prev = Fz3;
    state.Keff   = Keff3;

    state.gmin   = g_signed3;
    state.g_geo  = g_geo3;

    state.Fz     = Fz3;
    state.Mx     = Mx3;
    state.My     = My3;
end