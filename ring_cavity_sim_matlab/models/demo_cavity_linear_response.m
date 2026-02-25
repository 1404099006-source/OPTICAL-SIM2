function demo_cavity_linear_response()
set(0, 'defaultFigureColor', 'w');   % 所有新建 figure 默认白色背景
set(0, 'defaultAxesColor',   'w');   % 坐标区域默认白色

    % DEMO: 计算 4 镜环形腔的线性灵敏度矩阵 J，并用它把 20 个微动映射到光阑状态变化
    %
    % 结果：
    %   x_star : 4x1，本征光路在光阑面上的名义状态 (y,z,theta_y,theta_z)
    %   J      : 4x20，一阶雅可比矩阵
    %   dx     : 4x1，给定 dq 时的状态变化
    %
    % 你真正要用的就是： dx = J * dq;

    % 1. 建几何：光阑 + 4 面镜子
    [aperture, mirrors] = build_cavity_geometry_adjacent_spheres();

    % 2. 名义微动为 0
    q0 = zeros(20,1);   % 4 镜 * 每镜 5 DOF

    % 3. 先找一个本征光路的近似解 x_star
    x_init = [0;0;0;0]; % 初始猜测：光阑处位于轴心，角度 0
    x_star = find_eigen_ray(aperture, mirrors, q0, x_init);

    % 4. 数值线性化，求 A, B, J
    eps_x = 1e-6;   % 状态扰动步长（单位：mm / rad，对你的尺度足够小）
    eps_q = 1e-6;   % DOF 扰动步长（单位：mm / rad，同理）
    [A,B,J] = linearize_cavity(aperture, mirrors, q0, x_star, eps_x, eps_q);
    %整块灵敏度
    G = J(1:2, 16:20);   % 只看 M4 对 spot(y,z) 的 2D 灵敏度
disp('G = J(1:2,16:20) where u=[dx4 dy4 dz4 ax4 ay4]:');
disp(G);

s = vecnorm(G, 2, 1);  % 每列2范数
disp('Column norms s_i = ||G(:,i)|| : [dx4 dy4 dz4 ax4 ay4]');
disp(s);

% 看看2×5是否接近退化（对准会变得困难）
disp('singular values of G (2x5):');
disp(svd(G));


    fprintf('Nominal eigen state at aperture [y z ty tz]^T = [%.3e %.3e %.3e %.3e]\n', x_star);

%     镜片 1（右上球面镜）
% 
% dq(1) = Δx_1 —— 沿 局部 x₁ 方向平移（在腔平面内，差不多沿 M1→M2 那段光路方向）
% 
% dq(2) = Δy_1 —— 沿 局部 y₁ 方向平移（出平面方向，垂直于整个腔平面）
% 
% dq(3) = Δz_1 —— 沿 法向 z₁ 平移（piston，沿镜面法线前后推）
% 
% dq(4) = a_x1 —— 绕局部 x₁ 轴转角（tip/pitch，小角度，单位 rad）
% 
% dq(5) = a_y1 —— 绕局部 y₁ 轴转角（tilt/yaw，小角度，单位 rad）
% 
% 镜片 2（左上平面镜）
% 6. dq(6) = Δx_2
% 7. dq(7) = Δy_2
% 8. dq(8) = Δz_2
% 9. dq(9) = a_x2
% 10. dq(10) = a_y2
% 
% 镜片 3（左下平面镜）
% 11. dq(11) = Δx_3
% 12. dq(12) = Δy_3
% 13. dq(13) = Δz_3
% 14. dq(14) = a_x3
% 15. dq(15) = a_y3
% 
% 镜片 4（右下球面镜）
% 16. dq(16) = Δx_4
% 17. dq(17) = Δy_4
% 18. dq(18) = Δz_4
% 19. dq(19) = a_x4
% 20. dq(20) = a_y4
    dq = zeros(20,1);
    % 比如：第一块镜子的 Δy_1 = 5 um（沿出平面方向）
    dq(1) =0e-3;
    dq(2) = 0e-3;  % 单位这里用 mm：5μm = 5e-3 mm
    dq(3) =0;
    dq(4) =0;
    dq(5) =0;
    dq(11) = 0;
 dq(16)=5e-3;

    dx = J * dq;
    fprintf('dx from example dq: [dy dz dty dtz]^T = [%.3e %.3e %.3e %.3e]\n', dx);
  % ===== M4 平移自由度：3D + 切片 灵敏度图 =====
    plot_M4_translation_3D(J, x_star);

    % ===== M4 转动自由度：2D 灵敏度热力图 =====
plot_M4_rotation_3D(J, x_star);
demo_M4_move_and_hits();


end

       
function [ap, mirrors] = build_cavity_geometry_adjacent_spheres()
    % 标准 4 镜环形腔几何：
    % - 右边一竖排：上球面 M1，下球面 M4，中间是光阑
    % - 左边一竖排：上平面 M2，下平面 M3
    % - 光路顺序：ap -> M1 -> M2 -> M3 -> M4 -> ap
    % - 相邻镜面几何光程约为 L（这里取 L = 70 mm）
    % - 球面曲率半径 R = 6000 mm

    L = 70;        % 腔臂长度 ~ 70 mm
    R = 6000;      % 球面半径

    mirrors(4) = struct();

    % ========= 1. 定义名义光路上的 4 个反射点 =========
    % 坐标在 X-Z 平面内（Y 是出平面方向）
    % 右边为 X = +L/2，左边为 X = -L/2
    % 上下为 Z = ±L/2

    % M1：右上角球面镜上的名义打点
    P1 = [ +L/2; 0; +L/2 ];   % (x,y,z)

    % M2：左上角平面镜上的名义打点
    P2 = [ -L/2; 0; +L/2 ];

    % M3：左下角平面镜上的名义打点
    P3 = [ -L/2; 0; -L/2 ];

    % M4：右下角球面镜上的名义打点
    P4 = [ +L/2; 0; -L/2 ];

    % 光阑在右侧两球面之间的中点
    ap.A0 = [ +L/2; 0; 0 ];

    % ========= 2. 名义光路方向（保证能闭环） =========
    % 光路顺序：
    %   ap -> M1 -> M2 -> M3 -> M4 -> ap

    k_ap_to_M1 = [ 0; 0;  1];   % 从光阑往上到 M1
    k_M1_to_M2 = [-1; 0;  0];   % M1 -> 左上 M2
    k_M2_to_M3 = [ 0; 0; -1];   % M2 -> 左下 M3
    k_M3_to_M4 = [ 1; 0;  0];   % M3 -> 右下 M4
    k_M4_to_ap = [ 0; 0;  1];   % M4 -> 光阑（往上）

    % ========= 3. 每个镜面的法向（由入射/出射方向确定） =========
    % 镜面法向 n_j ∝ k_in - k_out

    n1 = k_ap_to_M1 - k_M1_to_M2;  n1 = n1 / norm(n1);   % 球面 1
    n2 = k_M1_to_M2 - k_M2_to_M3;  n2 = n2 / norm(n2);   % 平面 2
    n3 = k_M2_to_M3 - k_M3_to_M4;  n3 = n3 / norm(n3);   % 平面 3
    n4 = k_M3_to_M4 - k_M4_to_ap;  n4 = n4 / norm(n4);   % 球面 4

    % ========= 4. 定义 4 个镜子 =========

    % ----- Mirror 1: 球面，上右 -----
    mirrors(1).type      = 'sphere';
    mirrors(1).R         = R;
    mirrors(1).refPoint  = P1;             % 旋转中心先取打点
    mirrors(1).C0        = P1 - R * n1;    % 球心位置，使得 (P1 - C0)/R = n1
    mirrors(1).ez0       = n1;             % 局部 z 轴 = 法向
    mirrors(1).ex0       = k_M1_to_M2;     % 局部 x 轴先取沿光路切向
    mirrors               = normalize_local_axes(mirrors,1);

    % ----- Mirror 2: 平面，上左 -----
    mirrors(2).type      = 'plane';
    mirrors(2).R         = inf;
    mirrors(2).refPoint  = P2;
    mirrors(2).C0        = P2;             % 平面上一点即可
    mirrors(2).ez0       = n2;
    mirrors(2).ex0       = k_M2_to_M3;     % 切向
    mirrors               = normalize_local_axes(mirrors,2);

    % ----- Mirror 3: 平面，下左 -----
    mirrors(3).type      = 'plane';
    mirrors(3).R         = inf;
    mirrors(3).refPoint  = P3;
    mirrors(3).C0        = P3;
    mirrors(3).ez0       = n3;
    mirrors(3).ex0       = k_M3_to_M4;
    mirrors               = normalize_local_axes(mirrors,3);

    % ----- Mirror 4: 球面，下右 -----
    mirrors(4).type      = 'sphere';
    mirrors(4).R         = R;
    mirrors(4).refPoint  = P4;
    mirrors(4).C0        = P4 - R * n4;
    mirrors(4).ez0       = n4;
    mirrors(4).ex0       = k_M4_to_ap;
    mirrors               = normalize_local_axes(mirrors,4);

    % ========= 5. 光阑坐标系 =========
    % 光阑面法线沿着名义光路方向（从 ap 指向 M1 即 +z）
    ap.m  = k_ap_to_M1;              % 光阑法线 ≈ 光轴
    ap.m  = ap.m / norm(ap.m);

    % 光阑面内：
    %   ey = 出平面方向（取全局 Y）
    %   ez = 在腔平面内，和 m, ey 都正交（自动算出来，是沿着 -X）
    ap.ey = [0; 1; 0];
    ap.ey = ap.ey / norm(ap.ey);
    ap.ez = cross(ap.m, ap.ey);
    ap.ez = ap.ez / norm(ap.ez);

end


function mirrors = normalize_local_axes(mirrors,j)
    % 1) 先单位化 z 轴（法向）
    ez = mirrors(j).ez0;
    ez = ez / norm(ez);

    % 2) 想要 y 轴朝全局 +Y（往桌子外）
    ey_global = [0; 1; 0];

    % 对当前这套几何（法线都在 X-Z 平面），ey_global 本身就与 ez 正交，
    % 但写成通用一点：先减掉在 ez 方向上的分量，保证正交更鲁棒
    ey = ey_global - dot(ey_global, ez) * ez;
    if norm(ey) < 1e-12
        % 极端情况：如果 ez 恰好就是 ±Y（当前几何不会发生）
        % 就选一个备用方向，比如全局 Z，再正交化
        ey = [0; 0; 1];
        ey = ey - dot(ey, ez) * ez;
    end
    ey = ey / norm(ey);

    % 3) 用 y × z 得到 x 轴，自动与 y、z 都正交
    ex = cross(ey, ez);
    ex = ex / norm(ex);

    % 4) 回写到结构体
    mirrors(j).ez0 = ez;
    mirrors(j).ey0 = ey;
    mirrors(j).ex0 = ex;
end



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
function mirrors_eff = apply_DOFs_to_mirrors(mirrors, q)
    mirrors_eff = mirrors;
    for j = 1:4
        qj = q((j-1)*5 + (1:5));   % [Δx, Δy, Δz, ax, ay]^T

        ex0 = mirrors(j).ex0;
        ey0 = mirrors(j).ey0;
        ez0 = mirrors(j).ez0;

        % 平移（局部 -> 全局）
        tx = qj(1); ty = qj(2); tz = qj(3);
        t = tx*ex0 + ty*ey0 + tz*ez0;

        % 小转角矢量（局部 -> 全局）
        ax = qj(4); ay = qj(5);
        omega = ax*ex0 + ay*ey0;  % 绕局部 x,y 轴的小角度

        % 更新中心位置（假定绕 refPoint 旋转，这里用 C0/平面上一点）
        C0 = mirrors(j).C0;
        Om = mirrors(j).refPoint;
        dC_rot = cross(omega, (C0 - Om));   % 旋转引起的位移
        C  = C0 + t + dC_rot;

        % 更新法向和局部基底（线性近似）
        ez = ez0 + cross(omega, ez0);
        ez = ez / norm(ez);
        ex = ex0 + cross(omega, ex0);
        ex = ex - dot(ex,ez)*ez;
        ex = ex / norm(ex);
        ey = cross(ez, ex);
        ey = ey / norm(ey);

        mirrors_eff(j).C  = C;
        mirrors_eff(j).ez = ez;
        mirrors_eff(j).ex = ex;
        mirrors_eff(j).ey = ey;
    end
end
function [r_hit, k] = propagate_to_mirror(r, k, mirror)
    % 光线 (r,k) -> 与镜面交点 r_hit
    if strcmp(mirror.type,'sphere')
        % 球面：|r + t k - C|^2 = R^2
        C = mirror.C;
        R = mirror.R;
        oc = r - C;
        a = dot(k,k);
        b = 2*dot(oc,k);
        c = dot(oc,oc) - R^2;
        discriminant = b^2 - 4*a*c;
        if discriminant < 0
            error('No intersection with spherical mirror.');
        end
        t1 = (-b - sqrt(discriminant)) / (2*a);
        t2 = (-b + sqrt(discriminant)) / (2*a);
        t = min(t1,t2);
        if t < 0
            t = max(t1,t2);
        end
        if t <= 0
            error('Intersection behind ray start.');
        end
        r_hit = r + t*k;

    elseif strcmp(mirror.type,'plane')
        % 平面：(r - C)·n = 0，n = ez
        C = mirror.C;
        n = mirror.ez;
        denom = dot(k,n);
        if abs(denom) < 1e-12
            error('Ray parallel to plane.');
        end
        t = dot(C - r, n) / denom;
        if t <= 0
            error('Intersection behind ray start.');
        end
        r_hit = r + t*k;
    else
        error('Unknown mirror type.');
    end
end

function [r_out, k_out] = reflect_on_mirror(r_hit, k_in, mirror)
    % 在交点处计算法向并做理想镜面反射
    if strcmp(mirror.type,'sphere')
        C = mirror.C;
        R = mirror.R;
        n = (r_hit - C) / R;   % 球面法向
        n = n / norm(n);
    else
        n = mirror.ez;         % 平面镜，法向恒定
        n = n / norm(n);
    end

    k_out = k_in - 2*dot(k_in,n)*n;
    k_out = k_out / norm(k_out);
    r_out = r_hit;
end

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
function plot_M4_translation_3D(J, x_star)
    % M4 的 3 个平移 DOF: Δx_4, Δy_4, Δz_4
    % 用线性模型 δx = J δq，计算光阑中心距离 r = sqrt(y^2 + z^2)
    % 并在 3D 参数空间上画：
    %   1. 3D 散点“热力图” (Δx4, Δy4, Δz4) -> r
    %   2. 三个 2D 切片：Δx4=0, Δy4=0, Δz4=0

    % ---- 只取位置分量 y,z 的灵敏度（前两行） ----
    J_pos_M4 = J(1:2, 16:18);     % 2x3，对应 Δx4, Δy4, Δz4

    % 扫描范围（单位：mm）
    % 这里取 ±10 μm，可根据实际情况改
    amp_trans = 10e-3;           % 10 μm = 10e-3 mm
    N = 81;                       % 每个轴 21 个点

    dx4_vals = linspace(-amp_trans, amp_trans, N);
    dy4_vals = linspace(-amp_trans, amp_trans, N);
    dz4_vals = linspace(-amp_trans, amp_trans, N);

    [DX4, DY4, DZ4] = ndgrid(dx4_vals, dy4_vals, dz4_vals);   % N x N x N

    % 展平为向量，方便一次性矩阵乘
    dq_sub = [DX4(:)'; DY4(:)'; DZ4(:)'];   % 3 x (N^3)

    % y,z 的增量： 2 x (N^3)
    dYZ = J_pos_M4 * dq_sub;

    % 加上名义本征解的 y,z
    Y = x_star(1) + dYZ(1,:);    % 1 x (N^3)
    Z = x_star(2) + dYZ(2,:);

    % 到光阑中心的距离 r
    R = sqrt(Y.^2 + Z.^2);       % 1 x (N^3)

    % reshape 回 N x N x N
    R3 = reshape(R, N, N, N);

    % ================= 1. 3D 散点“热力图” =================
    figure; clf;
    scatter3(DX4(:)*1e3, DY4(:)*1e3, DZ4(:)*1e3, ...
             20, R(:), 'filled');    % 颜色代表 r
    xlabel('\Delta x_4 (\mum)', 'Interpreter','tex');
    ylabel('\Delta y_4 (\mum)', 'Interpreter','tex');
    zlabel('\Delta z_4 (\mum)', 'Interpreter','tex');
    title('M4 translation DOFs: distance from aperture center (3D)', ...
          'Interpreter','tex');
    grid on; box on; colorbar;
    colormap(jet);

    % ================= 2. 三个 2D 切片 =================
    % 选中间切片（0 附近），索引：
    mid = ceil(N/2);

    % (a) Δx4 = 0 切片：平面 (Δy4, Δz4)
    R_x0 = squeeze(R3(mid, :, :));  % N x N
    figure; clf;
    imagesc(dy4_vals*1e3, dz4_vals*1e3, R_x0');  % 注意转置，让坐标对应
    set(gca,'YDir','normal');
    xlabel('\Delta y_4 (\mum)', 'Interpreter','tex');
    ylabel('\Delta z_4 (\mum)', 'Interpreter','tex');
    title('\Delta x_4 = 0 slice: distance r(y,z)', 'Interpreter','tex');
    colorbar; colormap(jet);

    % (b) Δy4 = 0 切片：平面 (Δx4, Δz4)
    R_y0 = squeeze(R3(:, mid, :));  % N x N
    figure; clf;
    imagesc(dx4_vals*1e3, dz4_vals*1e3, R_y0');
    set(gca,'YDir','normal');
    xlabel('\Delta x_4 (\mum)', 'Interpreter','tex');
    ylabel('\Delta z_4 (\mum)', 'Interpreter','tex');
    title('\Delta y_4 = 0 slice: distance r(x,z)', 'Interpreter','tex');
    colorbar; colormap(jet);

    % (c) Δz4 = 0 切片：平面 (Δx4, Δy4)
    R_z0 = squeeze(R3(:, :, mid));  % N x N
    figure; clf;
    imagesc(dx4_vals*1e3, dy4_vals*1e3, R_z0');
    set(gca,'YDir','normal');
    xlabel('\Delta x_4 (\mum)', 'Interpreter','tex');
    ylabel('\Delta y_4 (\mum)', 'Interpreter','tex');
    title('\Delta z_4 = 0 slice: distance r(x,y)', 'Interpreter','tex');
    colorbar; colormap(jet);
end
function plot_M4_rotation_3D(J, x_star)
    % M4 的 2 个转动 DOF: a_x4, a_y4
    % 基于线性模型 δx = J δq，计算光阑中心距离 r = sqrt(y^2 + z^2)
    % 画：
    %   1) 3D surface: (a_x4, a_y4) -> r
    %   2) 俯视 2D 热力图（方便对比）

    % y,z 对 [a_x4, a_y4] 的灵敏度 (2x2)
    J_pos_rot_M4 = J(1:2, 19:20);

    % 扫描范围（单位：rad），这里取 ±100 μrad
    amp_rot = 100e-6;       % 100 μrad
    N = 81;                 % 网格分辨率

    ax4_vals = linspace(-amp_rot, amp_rot, N);
    ay4_vals = linspace(-amp_rot, amp_rot, N);

    [AX4, AY4] = ndgrid(ax4_vals, ay4_vals);   % N x N

    % 展成 2 x (N^2) 的小扰动向量
    dq_sub = [AX4(:)'; AY4(:)'];

    % 位置分量增量 (y,z)
    dYZ = J_pos_rot_M4 * dq_sub;   % 2 x (N^2)

    Y = x_star(1) + dYZ(1,:);
    Z = x_star(2) + dYZ(2,:);
    R = sqrt(Y.^2 + Z.^2);         % 1 x (N^2)

    R2 = reshape(R, N, N);         % N x N

    %% 1) 3D surface 图：高度和颜色都表示 r
    figure; clf;
    surf(ax4_vals*1e6, ay4_vals*1e6, R2', R2', 'EdgeColor','none');  % X=ax4, Y=ay4, Z=r
    xlabel('a_{x4} (\murad)', 'Interpreter','tex');
    ylabel('a_{y4} (\murad)', 'Interpreter','tex');
    zlabel('r (mm)');
    title('M4 rotation DOFs: distance from aperture center (3D surface)', ...
          'Interpreter','tex');
    colormap(jet);
    colorbar;
    grid on; box on;
    view(45, 30);          % 斜视角方便看
    shading interp;        % 颜色平滑一点

    %% 2) 俯视 2D 热力图（同一数据从上往下看）
    figure; clf;
    imagesc(ax4_vals*1e6, ay4_vals*1e6, R2');
    set(gca,'YDir','normal');
    xlabel('a_{x4} (\murad)', 'Interpreter','tex');
    ylabel('a_{y4} (\murad)', 'Interpreter','tex');
    title('M4 rotation DOFs: distance from aperture center (top view)', ...
          'Interpreter','tex');
    colormap(jet);
    colorbar;
end
function [Phits, x_out] = trace_cavity_hits(ap, mirrors, q, x_in)
    % TRACE_CAVITY_HITS
    % 从光阑状态 x_in 出发，带 DOF q，绕腔一圈，
    % 返回在 4 个镜子上的打点位置（全局坐标）和回到光阑的状态。
    %
    % 输出：
    %   Phits: 3x4 矩阵
    %          第 j 列 = [Xj; Yj; Zj] 为在镜子 Mj 上的打点
    %   x_out: 4x1，绕一圈回到光阑后的状态（应该接近 x_in）

    % 1. 应用 DOF
    mirrors_eff = apply_DOFs_to_mirrors(mirrors, q);

    % 2. 光阑状态 -> 光线 (r,k)
    [r, k] = state_to_ray(ap, x_in);

    Phits = nan(3,4);   % 预分配

    % ---- ap -> M1 ----
    [r, k] = propagate_to_mirror(r, k, mirrors_eff(1));
    Phits(:,1) = r;                      % 记录 M1 打点
    [r, k] = reflect_on_mirror(r, k, mirrors_eff(1));

    % ---- M1 -> M2 ----
    [r, k] = propagate_to_mirror(r, k, mirrors_eff(2));
    Phits(:,2) = r;                      % 记录 M2 打点
    [r, k] = reflect_on_mirror(r, k, mirrors_eff(2));

    % ---- M2 -> M3 ----
    [r, k] = propagate_to_mirror(r, k, mirrors_eff(3));
    Phits(:,3) = r;                      % 记录 M3 打点
    [r, k] = reflect_on_mirror(r, k, mirrors_eff(3));

    % ---- M3 -> M4 ----
    [r, k] = propagate_to_mirror(r, k, mirrors_eff(4));
    Phits(:,4) = r;                      % 记录 M4 打点
    [r, k] = reflect_on_mirror(r, k, mirrors_eff(4));

    % ---- M4 -> ap ----
    [r, k] = propagate_to_plane(r, k, ap.A0, ap.m);

    % 回到光阑状态向量
    x_out = ray_to_state(ap, r, k);
end
function demo_M4_move_and_hits()
    % 示例：只动 M4 的 5 个 DOF，
    % 求【新腔】的闭环本征光路，然后看 M1~M4 上的打点坐标

    % 1. 基本参数 & 几何
    L = 70;
    R = 6000;
    [ap, mirrors] = build_cavity_geometry_adjacent_spheres();

    % 2. 先求名义腔(q=0)的本征光路，作为后续迭代初始值
    q_zero = zeros(20,1);
    x_init = [0;0;0;0];
    x_star0 = find_eigen_ray(ap, mirrors, q_zero, x_init);

    % 3. 构造一个 q：只动 M4 的 5 个 DOF
    q = zeros(20,1);
    % M4 在 16~20：
    %   q(16) = Δx4
    %   q(17) = Δy4
    %   q(18) = Δz4
    %   q(19) = a_x4
    %   q(20) = a_y4

    % ------- 这里你自己设想演示用的位姿 -------
    % 提醒：为了画示意图好看，可以故意用“夸张一点”的数值
    % 比如 0.5 mm、几 mrad，而不是 μm/μrad
    dx4 = 0e-3;        % 0.5 mm，示意图用
    dy4 = 0.0e-3;
    dz4 = 0.0e-3;
    ax4 = 5e-3;        % 0.005 rad ~ 0.3°
    ay4 = 0e-3;       % 3 mrad，示意

    q(16) = dx4;
    q(17) = dy4;
    q(18) = dz4;
    q(19) = ax4;
    q(20) = ay4;

    % 4. 对这个“新腔”，再求一次闭环本征光路
    x_star_q = find_eigen_ray_with_q(ap, mirrors, q, x_star0);

    % 5. 用这个新的本征光路，追一圈打点
    [Phits, x_out] = trace_cavity_hits(ap, mirrors, q, x_star_q);

    % 6. 打印结果（用多一点小数位，看得更清楚）
    fprintf('=== M4 motion: [dx4 dy4 dz4 ax4 ay4] = [%.3e %.3e %.3e %.3e %.3e]\n', ...
        dx4, dy4, dz4, ax4, ay4);
    fprintf('New eigen state at aperture [y z ty tz]^T = [%.6f %.6f %.6e %.6e]\n', ...
        x_star_q(1), x_star_q(2), x_star_q(3), x_star_q(4));

    for j = 1:4
        P = Phits(:,j);
        fprintf('Hit on M%d: [X Y Z] = [%.6f  %.6f  %.6f] mm\n', ...
            j, P(1), P(2), P(3));
    end

    % 7. 看一下闭环误差（理论上应该很小）
    dx_loop = x_out - x_star_q;
    fprintf('Loop closure delta at aperture [dy dz dty dtz] = [%.3e %.3e %.3e %.3e]\n', ...
        dx_loop(1), dx_loop(2), dx_loop(3), dx_loop(4));
end

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


