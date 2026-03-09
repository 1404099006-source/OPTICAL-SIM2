function Dc = contact_drift_update(Dc, Fz, Keff, seated, P)
%CONTACT_DRIFT_UPDATE Three-stage drift model:
% A) free-space: low random drift, near-zero mean
% B) critical contact: higher directional + random drift (peak)
% C) locked/pressed: random drift reduces, directional drift tends to lock value

if P.drift_mode == "smooth"
    % ---------- defaults ----------
    if ~isfield(P,'F_lock');            P.F_lock = 0.8; end
    if ~isfield(P,'drift_peak_scale');  P.drift_peak_scale = 1.0; end
    if ~isfield(P,'drift_lock_scale');  P.drift_lock_scale = 0.45; end
    if ~isfield(P,'kappa_lock');        P.kappa_lock = 2.0; end
    if ~isfield(P,'rw_scale_free');     P.rw_scale_free = 0.6; end
    if ~isfield(P,'rw_scale_peak');     P.rw_scale_peak = 2.0; end
    if ~isfield(P,'rw_scale_lock');     P.rw_scale_lock = 0.35; end
    if ~isfield(P,'creep_rate');        P.creep_rate = 0.03; end
    if ~isfield(P,'K_max');             P.K_max = inf; end

    Ft = P.F_touch;
    Fl = max(P.F_lock, Ft + 1e-6);
    D_peak = P.drift_max * P.drift_peak_scale;
    D_lock = P.drift_max * P.drift_lock_scale;

    % ---------- stage logic ----------
    if Fz <= Ft
        % Stage A: free-space / before meaningful contact
        D_target = zeros(5,1);
        rw_scale = P.rw_scale_free;

    elseif Fz < Fl
        % Stage B: critical contact (arc contact -> area expansion)
        s = (Fz - Ft) / (Fl - Ft); % 0..1
        D_target = D_peak * s;
        rw_scale = P.rw_scale_free + (P.rw_scale_peak - P.rw_scale_free) * s;

    else
        % Stage C: pressed/locked zone
        s_lock = 1 - exp(-P.kappa_lock * (Fz - Fl));
        D_target = D_peak + (D_lock - D_peak) * s_lock;

        % random term decays from peak to lock level
        rw_scale = P.rw_scale_lock + (P.rw_scale_peak - P.rw_scale_lock) * exp(-P.kappa_lock * (Fz - Fl));
        if seated
            rw_scale = min(rw_scale, P.rw_scale_lock);
        end
    end

    % If stiffness is abnormally high before lock, keep some extra randomness
    if (~seated) && (Fz > Ft) && (Fz < Fl) && (Keff > P.K_max)
        rw_scale = 1.25 * rw_scale;
    end

    % dynamic relaxation + random walk
    rho = P.rho_drift;
    rw = (P.sigma_drift_rw .* randn(5,1)) * rw_scale;
    Dc = (1-rho)*Dc + rho*D_target + rw;

    % Creep toward lock-direction target in stage C
    if Fz >= Fl
        Dc = Dc + P.creep_rate * (D_lock - Dc);
    end

    % Optional legacy restriction
    if isfield(P,'drift_rot_only') && P.drift_rot_only
        Dc(1:3) = 0;
    end

elseif P.drift_mode == "jump"
    % optional: implement event-driven jump later
end
end
