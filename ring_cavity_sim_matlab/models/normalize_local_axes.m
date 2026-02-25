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