function optics = optics_init_linear(P)
%OPTICS_INIT_LINEAR Build cavity geometry and linearize once: dx = J*dq

    % 这些函数来自你现有的光路模型代码（你那份 txt 里）
    [ap, mirrors] = build_cavity_geometry_adjacent_spheres();

    q0 = zeros(20,1);
    x_init = [0;0;0;0];
    x_star = find_eigen_ray(ap, mirrors, q0, x_init);

    % 线性化步长：按你的模型单位（通常 mm / rad）
    eps_x = 1e-6;
    eps_q = 1e-6;

    [~,~,J] = linearize_cavity(ap, mirrors, q0, x_star, eps_x, eps_q);

    optics.ap = ap;
    optics.mirrors = mirrors;
    optics.x_star = x_star;
    optics.J = J;
end
