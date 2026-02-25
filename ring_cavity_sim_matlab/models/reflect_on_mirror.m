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