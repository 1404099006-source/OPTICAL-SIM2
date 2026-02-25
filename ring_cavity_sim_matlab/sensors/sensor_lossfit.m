function [Lfit, ok, info] = sensor_lossfit(x_true, e, Keff, seated, P)
%SENSOR_LOSSFIT Simulated curve-fitting output for cavity loss/linewidth/contrast metric.
% Inputs:
%   x_true: true pose [u v z thx thy]
%   e: vision error (2x1, normalized)
%   Keff: effective stiffness proxy (N/um)
%   seated: bool
% Outputs:
%   Lfit: measured fitted loss (NaN if fail)
%   ok: success flag
%   info: struct with fields p_succ, is_outlier, sigma_rel

% ---- true loss ----
Ltrue = loss_true(x_true, P);
% ---- ideal fit option: always succeed ----
if isfield(P,'fit_always_ok') && P.fit_always_ok
    ok = true;
    info.p_succ = 1.0;
    info.is_outlier = false;
    info.sigma_rel = 0.0;
    info.Ltrue = Ltrue;

    % 你可以选择：完全无噪声 or 仅保留一点点测量噪声
    if isfield(P,'loss_abs_sigma_ppm') && P.loss_abs_sigma_ppm > 0
        Lfit = Ltrue + P.loss_abs_sigma_ppm * randn();  % 只有加性噪声
    else
        Lfit = Ltrue;  % 完全理想
    end
    Lfit = max(0, Lfit);
    return;
end



% ---- success probability (soft gating) ----
e_norm = norm(e);

% (A) 基础成功率：保证拟合不是“永远失败”
% 建议在仿真阶段让它>=0.35，这样 Lfit 会持续出现，便于观察闭环效果
if ~isfield(P,'pfit_floor'); P.pfit_floor = 0.35; end

% (B) 视觉门控：e 越小越容易拟合
p_e = 1 ./ (1 + exp((e_norm - P.e_ok)/P.se_ok));   % 0~1, e小=>接近1

% (C) 刚度门控：你原来写的是 Keff>K_max 就几乎0，太狠了
% 改成“缓惩罚”：用一个更温和的 logistic，且允许 P.K_soft > P.K_max
% 推荐：P.K_soft 设成你实际 Keff 的 1.2~2 倍
if ~isfield(P,'K_soft'); P.K_soft = 1.5 * P.K_max; end
if ~isfield(P,'sK_soft'); P.sK_soft = 2.0 * P.sK; end

p_k = 1 ./ (1 + exp((Keff - P.K_soft)/P.sK_soft)); % Keff越大越低，但下降慢很多

% (D) seating boost：保持你原来的逻辑，但稍微加强一些（可选）
seat_boost = 0.20 * double(seated);

% (E) 合成成功率：用 floor + 加权，而不是“乘积直接压死”
% 你原来是 p_e*p_k，任何一个小都会把 p_succ 变成很小
% 这里用 convex blend 更稳：p = floor + (1-floor)*(w*pe + (1-w)*pk)
if ~isfield(P,'w_pe'); P.w_pe = 0.70; end  % 更相信视觉质量
p_core = P.w_pe * p_e + (1 - P.w_pe) * p_k;

p_succ = P.pfit_floor + (1 - P.pfit_floor) * p_core + seat_boost;
p_succ = min(0.98, max(0.02, p_succ));     % 夹住避免0/1极端

ok = (rand() < p_succ);

% ---- relative noise level (worse when misaligned/hard contact) ----
% 你原噪声项：*(1 + 0.6*max(0, Keff/P.K_max)) 也容易爆
% 改成更温和的线性：Keff 越大噪声稍涨，但不会把 Lfit 弄成离谱
if ~isfield(P,'k_noise'); P.k_noise = 0.25; end

sigma_rel = P.loss_rel_sigma0 * (1 + 0.8*(e_norm/P.e_ok)^2) * (1 + P.k_noise*(Keff/max(1e-9,P.K_soft)));

% ---- outlier probability ----
% outlier 不要太大，否则曲线看起来像随机脉冲
if ~isfield(P,'out_p0'); P.out_p0 = 0.01; end
if ~isfield(P,'out_p1'); P.out_p1 = 0.03; end
if ~isfield(P,'out_mag'); P.out_mag = 0.25; end

p_out = P.out_p0 + P.out_p1 * min(1, e_norm / P.e_ok);
is_out = (rand() < p_out);

if ~ok
    Lfit = NaN;
else
    % multiplicative + additive noise (ppm)
    if ~isfield(P,'loss_abs_sigma_ppm'); P.loss_abs_sigma_ppm = 30; end
    Lfit = Ltrue * (1 + sigma_rel * randn()) + P.loss_abs_sigma_ppm * randn();

    if is_out
        Lfit = Lfit * (1 + P.out_mag);
    end

    % 防止拟合出负值（数值上更稳）
    Lfit = max(0, Lfit);
end

info.p_succ = p_succ;
info.is_outlier = is_out;
info.sigma_rel = sigma_rel;
info.Ltrue = Ltrue;

end
