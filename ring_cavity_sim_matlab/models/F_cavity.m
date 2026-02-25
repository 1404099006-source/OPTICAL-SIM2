function x_out = F_cavity(ap, mirrors, q, x_in)
    % F(q,x): 从光阑状态 x_in 出发，经过 4 面镜子，再回到光阑，得到 x_out

    % 1. 应用 DOF
    mirrors_eff = apply_DOFs_to_mirrors(mirrors, q);

    % 2. 光阑状态 -> 光线 (r,k)
    [r, k] = state_to_ray(ap, x_in);

    % 3. ap -> M1 (球面1)
    [r, k] = propagate_to_mirror(r, k, mirrors_eff(1));
    [r, k] = reflect_on_mirror(r, k, mirrors_eff(1));

    % 4. M1 -> M2 (球面2)
    [r, k] = propagate_to_mirror(r, k, mirrors_eff(2));
    [r, k] = reflect_on_mirror(r, k, mirrors_eff(2));

    % 5. M2 -> M3 (平面3)
    [r, k] = propagate_to_mirror(r, k, mirrors_eff(3));
    [r, k] = reflect_on_mirror(r, k, mirrors_eff(3));

    % 6. M3 -> M4 (平面4)
    [r, k] = propagate_to_mirror(r, k, mirrors_eff(4));
    [r, k] = reflect_on_mirror(r, k, mirrors_eff(4));

    % 7. M4 -> ap
    [r, k] = propagate_to_plane(r, k, ap.A0, ap.m);

    % 8. 回到光阑状态向量
    x_out = ray_to_state(ap, r, k);
end