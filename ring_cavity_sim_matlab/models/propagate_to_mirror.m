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