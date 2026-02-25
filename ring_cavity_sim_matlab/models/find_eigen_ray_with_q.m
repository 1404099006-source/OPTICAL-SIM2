function x_star = find_eigen_ray_with_q(ap, mirrors, q, x_init)
    % 在给定 DOF q 的情况下，求闭环本征光路：
    %   x* 满足 x* = F_cavity(ap, mirrors, q, x*)
    maxIter = 200;
    x = x_init;
    for k = 1:maxIter
        x_new = F_cavity(ap, mirrors, q, x);
        if norm(x_new - x) < 1e-10
            break;
        end
        x = x_new;
    end
    x_star = x;
end
