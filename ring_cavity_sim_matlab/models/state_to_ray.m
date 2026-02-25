function [r, k] = state_to_ray(ap, x)
    y  = x(1);
    z  = x(2);
    ty = x(3);
    tz = x(4);

    % 位置
    r = ap.A0 + y*ap.ey + z*ap.ez;

    % 方向：k ≈ m + ty*ey + tz*ez
    k_approx = ap.m + ty*ap.ey + tz*ap.ez;
    k = k_approx / norm(k_approx);
end
