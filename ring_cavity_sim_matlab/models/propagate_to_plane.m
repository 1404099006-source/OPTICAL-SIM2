function [r_hit, k_out] = propagate_to_plane(r, k, A0, m)
    % 光线 (r,k) -> 平面 (r - A0)·m = 0
    denom = dot(k,m);
    if abs(denom) < 1e-12
        error('Ray parallel to plane.');
    end
    t = dot(A0 - r, m) / denom;
    if t <= 0
        error('Plane is behind ray start.');
    end
    r_hit = r + t*k;
    k_out = k;
end