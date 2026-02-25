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