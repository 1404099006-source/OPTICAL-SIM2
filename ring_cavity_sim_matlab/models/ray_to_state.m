function x = ray_to_state(ap, r, k)
    % 位置投影
    dr = r - ap.A0;
    y  = dot(dr, ap.ey);
    z  = dot(dr, ap.ez);

    % 方向倾角
    km = dot(k, ap.m);
    ky = dot(k, ap.ey);
    kz = dot(k, ap.ez);
    ty = ky / km;
    tz = kz / km;

    x = [y; z; ty; tz];
end