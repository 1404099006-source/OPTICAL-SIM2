function x_star = find_eigen_ray(ap, mirrors, q0, x_init)
    % 简单迭代 F(x) 找闭环解： x* = F(q0, x*)
    maxIter = 100;
    x = x_init;
    for k = 1:maxIter
        x_new = F_cavity(ap, mirrors, q0, x);
        if norm(x_new - x) < 1e-9
            break;
        end
        x = x_new;
    end
    x_star = x;
end