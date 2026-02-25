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
    if ~isfield(P,'c_th');      P.c_th = 0.0; end
    if ~isfield(P,'enable_leveling'); P.enable_leveling = false; end

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
    [Fz, Mx, My, Keff_c, g_geo_min, g_signed_min] = force_model(Xtmp, P);

    % ---------------------------
    % 4) Seating flag
    % ---------------------------
    if action.do_seat
        state.seated = true;
    end

    % ---------------------------
    % 5) Force-induced deformation (z + tilt relaxation)
    % ---------------------------
    Dc_target = Xtrue.Dc;

    if state.seated && (Fz > 0)
        % axial compliance: under positive contact force, z tends to "move further into contact".
        % With your current force_model definition (w=max(-g,0)), you have been running z-up
        % convention successfully, so keep the same sign you are currently using:
        Dc_target(3) = +P.c_z * Fz;

        % tilt relaxation by restoring moments
        Dc_target(4) = -P.c_th * Mx;
        Dc_target(5) = -P.c_th * My;
    else
        Dc_target(3) = 0;
        Dc_target(4) = Xtrue.Dc(4);
        Dc_target(5) = Xtrue.Dc(5);
    end

    alpha = min(max(P.alpha_def, 0.0), 1.0);
    Xtrue.Dc(3:5) = (1-alpha)*Xtrue.Dc(3:5) + alpha*Dc_target(3:5);

    % ---------------------------
    % 6) Optional heuristic leveling (if you keep it)
    % ---------------------------
    if isfield(P,"enable_leveling") && P.enable_leveling
        Xtmp2 = Xtrue.xr + Xtrue.Dc;
        [Fz2, ~, ~, ~, ~, ~] = force_model(Xtmp2, P);

        th_old = Xtmp2(4:5);
        th_new = theta_level_update(th_old, Fz2, state.seated, P);

        Xtrue.Dc(4:5) = th_new - Xtrue.xr(4:5);
    end

    % ---------------------------
    % 7) Drift (optional)
    % ---------------------------
    Keff = Keff_c;
    if exist('contact_drift_update','file') == 2
        Xtrue.Dc = contact_drift_update(Xtrue.Dc, Fz, Keff, state.seated, P);
    end

    % ---------------------------
    % 8) Final true pose
    % ---------------------------
    Xtrue.x = Xtrue.xr + Xtrue.Dc;

    % ---------------------------
    % 9) Recompute contact at final pose for consistent outputs
    % ---------------------------
    [Fz3, Mx3, My3, Keff3, g_geo3, g_signed3] = force_model(Xtrue.x, P);

    % ---------------------------
    % 10) Update state
    % ---------------------------
    state.F_prev = Fz3;
    state.Keff   = Keff3;

    state.gmin   = g_signed3;
    state.g_geo  = g_geo3;

    state.Fz     = Fz3;
    state.Mx     = Mx3;
    state.My     = My3;
end