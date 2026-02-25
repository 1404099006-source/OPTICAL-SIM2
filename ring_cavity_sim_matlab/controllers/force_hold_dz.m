function dz = force_hold_dz(Fz, F_target, tol, dz_down, dz_up)
%FORCE_HOLD_DZ Simple bang-bang force holding via z motion.
% If force too low -> move down (dz_down negative)
% If force too high -> move up  (dz_up positive)
% If within band -> hold 0

if Fz < (F_target - tol)
    dz = dz_down;
elseif Fz > (F_target + tol)
    dz = dz_up;
else
    dz = 0;
end
end
