function [A,B,J] = linearize_cavity(ap, mirrors, q0, x0, eps_x, eps_q)
    % 线性化单圈映射 F(q,x) 在 (q0,x0) 处的一阶导：
    %   A = dF/dx (4x4)
    %   B = dF/dq (4x20)
    % 然后闭环本征解的响应： (I - A) * dx = B * dq  =>  dx = (I-A)^(-1) B dq

    % ---- A: 对 x 的导数 ----
    A = zeros(4,4);
    f0 = F_cavity(ap, mirrors, q0, x0);
    for i = 1:4
        dx = zeros(4,1);
        dx(i) = eps_x;
        fi = F_cavity(ap, mirrors, q0, x0 + dx);
        A(:,i) = (fi - f0) / eps_x;
    end

    % ---- B: 对 q 的导数 ----
    B = zeros(4,20);
    for k = 1:20
        dq = zeros(20,1);
        dq(k) = eps_q;
        fk = F_cavity(ap, mirrors, q0 + dq, x0);
        B(:,k) = (fk - f0) / eps_q;
    end

    % ---- 闭环雅可比 ----
    J = (eye(4) - A) \ B;
end