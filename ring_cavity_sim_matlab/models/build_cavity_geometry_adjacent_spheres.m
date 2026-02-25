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