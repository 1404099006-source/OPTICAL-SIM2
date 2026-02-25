function Dc = contact_drift_update(Dc, Fz, Keff, seated, P)
%CONTACT_DRIFT_UPDATE Update hidden contact-induced drift Δx_c.
% Dc: current drift (5x1), Fz: normal force, Keff: effective stiffness,
% seated: boolean, P: params.

if P.drift_mode == "smooth"
    % target drift magnitude (force-driven saturation)
    if Fz <= P.F_touch
        D_target = zeros(5,1);
    else
        s = 1 - exp(-P.kappa_drift * (Fz - P.F_touch));
        D_target = P.drift_max * s; % elementwise scale
    end

    % dynamic relaxation + random walk (seating reduces randomness)
    rho = P.rho_drift;
    rw = P.sigma_drift_rw .* randn(5,1);
    if seated
        rw = 0.5 * rw; % seating stabilizes
    else
        % if contact is "hard", random effect could be larger
        if Keff > P.K_max
            rw = 2.0 * rw;
        end
    end

    Dc = (1-rho)*Dc + rho*D_target + rw;

    % Optional: restrict drift to rotation only
    if isfield(P,'drift_rot_only') && P.drift_rot_only
        Dc(1:3) = 0;
    end

elseif P.drift_mode == "jump"
    % optional: implement event-driven jump later (we can add in next step)
end
end