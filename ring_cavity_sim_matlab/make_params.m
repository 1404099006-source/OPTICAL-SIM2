function P = make_params()
%MAKE_PARAMS All simulation parameters in one place.

% ---------- Units ----------
P.unit_pos = "um";   % all positions in micrometers
P.unit_ang = "rad";  % internal angle in radians
% ===============================
% Coordinate / bias definitions
% ===============================
% xC: lens relative cavity small pose [u v z thx thy]
% optics model may not be perfectly aligned with cavity zero:
% x_opt = xC + xC_to_opt_bias
if ~isfield(P,'xC_to_opt_bias')
    P.xC_to_opt_bias = zeros(5,1);  % later you can tune this to align e(0)=0
end

% Initial misalignment relative cavity
if ~isfield(P,'xC_init')
    % Assume initial tilt within ~30 arcsec
    arcsec30 = 30/3600; % deg
    P.xC_init = [80; -60; -200; deg2rad(arcsec30); -deg2rad(arcsec30)];
end

% Initial robot-nominal vs lens-true mismatch (x - xr at step 0)
% [u; v; z; thx; thy] in [um; um; um; rad; rad]
if ~isfield(P,'Dc_init')
    % small initial robot-vs-lens mismatch: [u; v; z; thx; thy]
    P.Dc_init = [1.5; -1.0; 0.8; 2e-6; -2e-6];
end

% ===============================
% Vision normalization / stage switching (set later in block)
% ===============================

% During loss optimization, keep spot error within this guard
if ~isfield(P,'e_guard'); P.e_guard = 0.65; end

% ===============================
% LS uv controller params
% ===============================
if ~isfield(P,'ls_lambda');  P.ls_lambda  = 1e-3; end
if ~isfield(P,'ls_gain');    P.ls_gain    = 1.0; end
if ~isfield(P,'duv_ls_max'); P.duv_ls_max = 1.0; end

% ===============================
% Force-hold params (z-up)
% ===============================
if ~isfield(P,'dz_max'); P.dz_max = 0.5; end
if ~isfield(P,'kF_p');   P.kF_p   = 0.25; end
if ~isfield(P,'kF_i');   P.kF_i   = 0.00; end
if ~isfield(P,'dz_lpf'); P.dz_lpf = 0.90; end

% ===============================
% Debug printing (optional)
% ===============================
if ~isfield(P,'debug_print_every'); P.debug_print_every = 200; end


% ---------- Time / loop ----------
P.N = 2000;          % max steps
P.dt = 1.0;          % abstract step time

% ---------- Geometry (for vision normalization & contact) ----------
P.aperture_rx = 200;   % um, normalization radius in u direction
P.aperture_ry = 200;   % um, normalization radius in v direction
P.contact_R   = 5000;  % um, effective contact radius for tilt-to-gap (e.g., 5 mm)

% ---------- True loss landscape (L = f_L(x)) ----------

% weights for [u v z thx thy]
P.W = diag([1.0, 1.2, 0.8, 200.0, 200.0]);  % tune later
P.rug_amp = 0.02;    % small ruggedness amplitude (relative)
P.rug_k   = 2*pi/200; % spatial frequency (1/um)

% target optimum (ideal alignment)
P.x_star = [0; 0; 0; 0; 0];
P.loss_probe_h_um = 2.0;   % 2D扫频半径(um)，1~5都行
P.e_guard         = 1.0;   % 损耗阶段允许的最大 ||e||（放宽/收紧都可）

% ---- optics init ----
P.optics = optics_init_linear(P);
P.optics_mirror_id = 4;


% ---------- Contact-induced drift Δx_c (hidden) ----------
P.drift_mode = "smooth";  % "smooth" or "jump"
P.F_touch = 0.1;          % N, contact threshold
P.kappa_drift = 0.8;      % 1/N, drift saturation speed
P.drift_rot_only = false; % allow translational + rotational drift
P.drift_max = [0; 0; 0; 2e-5; -2e-5]; % [um;um;um;rad;rad], max drift

% Smooth drift dynamics (optional)
P.rho_drift = 0.15;       % update rate in dynamic version
% per-step random walk [u; v; z; thx; thy]
% translation terms are set to small non-zero values (um) to model free-space
% repeatability + contact micro-slip induced drift.
P.sigma_drift_rw = [0.02; 0.02; 0.01; 2e-6; 2e-6];

% Three-stage drift (A free-space / B critical-contact / C lock+creep)
P.F_lock = 0.8;             % N, enter lock/pressed stage threshold
P.drift_peak_scale = 1.0;   % stage-B directional drift peak scale
P.drift_lock_scale = 0.45;  % stage-C steady directional drift scale
P.kappa_lock = 2.0;         % lock-stage transition rate
P.rw_scale_free = 0.6;      % random drift scale in stage A
P.rw_scale_peak = 2.0;      % random drift peak scale in stage B
P.rw_scale_lock = 0.35;     % random drift scale in stage C
P.creep_rate = 0.03;        % stage-C slow creep toward lock target

% Jump drift parameters (if used)
P.Fc_mean = 3.0;          % N
P.Fc_std  = 0.6;          % N
P.J_mean  = [15; -10; 3; 1e-5; -1e-5];
P.J_std   = [4; 4; 1; 3e-6; 3e-6];

% ---------- Force model Fz = f_F(x) ----------
P.n_force = 1.2;       % exponent (1.0 linear, 1.5 Hertz-like)
P.k_force = 0.010;     % N/(um^n), tune to reach 5~10N at desired indentation
P.F_max   = 12.0;      % safety cap in sim (optional)

% ---------- Fiber pose sensor (z, thx, thy) ----------
P.fiber_sigma_z = 1.0;                % um (1σ)
P.fiber_sigma_th = deg2rad(0.005);    % rad (1σ)
P.fiber_bias_rw = [0.05; deg2rad(0.0002); deg2rad(0.0002)]; % per-step drift std

% contact-dependent fiber degradation (before seating)
P.fiber_contact_scale = 3.0;  % noise multiplier in contact before seating
P.fiber_seated_scale  = 1.0;  % noise multiplier after seating

% ---------- Vision sensor e ----------
P.vision_sigma = 0.0;  % normalized units (1σ) after normalization

% ---------- Loss fitting sensor (fit success + noise + outliers) ----------
P.e_ok = 0.8;           % if ||e|| below this, fitting likely succeeds
P.se_ok = 0.08;         % sigmoid softness for success probability
P.K_max = 0.08;         % N/um (example), threshold for "too hard"
P.sK    = 0.02;

P.loss_rel_sigma0 = 0.03; % base relative noise (3%)
P.out_p0 = 0.02;          % base outlier prob
P.out_p1 = 0.08;          % extra outlier prob when misaligned
P.out_mag = 0.20;         % outlier magnitude (relative, e.g. +20%)

% ---------- Gating windows ----------
P.F_pre_min = 0.5;   % N
P.F_pre_max = 1.2;   % N
% ---------- Force-continuation schedule ----------
P.force_targets = [0.4, 0.6, 0.8, 1.0];  % N 低力阶梯，贴合前控制
P.force_tol     = 0.12;   % N 目标力允许误差带（越小越严格，越大越稳）
P.force_settle_steps = 8; % 每个台阶先稳住力若干步再开始fine
P.min_level_before_exit = numel(P.force_targets); % 阶梯力优先：到该级后才允许loss触发退出
P.final_attach_enable = true;  % 最后一顶（贴合阶段）
P.dz_attach_pulse = 0.05;      % um/step, 贴合阶段额外z推进
P.attach_push_steps = 6;       % 额外推进步数
P.attach_hold_steps = 12;      % 推进后保持步数
P.attach_force_tol = 0.15;     % N, final_attach阶段力误差容忍
% ===== Loss model params =====
% Legacy empirical model (kept for fallback/ablation)
P.L_min  = 200;      % ppm best floor (legacy)
P.Gp_ppm = 2500;     % legacy position penalty
P.Ga_ppm = 388;      % legacy angle penalty
P.Gc_ppm = 200;      % legacy coupling
P.ang_ref_rad = deg2rad(0.02);

% Physical loss model (recommended)
% L_ppm = 1e6 * (ell0 + ell_clip(rho)), rho from optics aperture offset (dy,dz)
P.use_physical_loss = true;
P.L0_ppm = 500;                 % fixed round-trip baseline loss (ppm)
P.aperture_phys_radius_mm = 0.5; % real aperture radius a (mm)
% Gaussian mode radii at aperture plane (45 deg astigmatism equivalent)
P.w_t_mm = 0.395;
P.w_s_mm = 0.333;
P.w_eff_mm = sqrt(P.w_t_mm * P.w_s_mm);
% clip overlap evaluation: "approx_radial" (fast) or "numeric" (integral)
P.loss_clip_method = "approx_radial";
P.loss_clip_grid_n = 81;        % used when method="numeric"

% Fixed bias on optics output: [dy(mm); dz(mm); dty(rad); dtz(rad)]
% Reasonable: position bias ~0.02~0.05 mm if aperture radius ~0.5 mm
% Angle bias very small



% 每个力台阶的fine预算（迭代次数）
P.fine_iter_per_level = 25;   % 每级精调迭代数（一次迭代≈5次拟合）
% ---------- Moment-driven leveling (Scheme A) ----------
% Angle is driven by contact moments Mx/My (in truth_update), not by time-decay-to-zero.
P.enable_leveling = false;      % disable legacy exponential theta decay path
P.theta_damp = 0.20;            % 0: pure moment equilibrium, 1: hold previous angle
P.dtheta_max = deg2rad(0.02);   % rad/step, cap per-step angle increment

% ---------- Theoretical optical optimum (no micro-motion) ----------
% Use this as the "true" loss optimum. Controller doesn't know it.
P.z0_star = 0;
P.x_star = [0; 0; P.z0_star; 0; 0];  % [u v z thx thy]
% 如果你没有z0_star，就先用0或你初始化的z
% P.z0_star = 0;  % (um)
P.pfit_floor = 0.35;     % 拟合最小成功率（仿真阶段）
P.K_soft     = 0.45;     % 让它略高于你常见 Keff (比如0.38)
P.sK_soft    = 0.10;     % 刚度门控平滑程度
P.w_pe       = 0.70;     % 视觉 vs 刚度 权重
P.k_noise    = 0.25;     % 刚度对噪声的影响
P.loss_abs_sigma_ppm = 30;
P.out_p0     = 0.01;
P.out_p1     = 0.03;
P.out_mag    = 0.25;
P.fit_always_ok = true;   % 理想拟合：每次都成功
P.act_enable   = false;  % true 才启用
P.act_bias_uv  = [0;0];  % um/step
P.act_bias_z   = 0;      % um/step
P.act_noise_uv = 0;      % um std/step
P.act_noise_z  = 0;      % um std/step

% --- bias 1: equivalent optics q0 (20x1) ---
% 含义：当 x_C=0 时，在光学模型里等效为 q0 的偏置（C->O不重合、装夹偏差等）
P.optics_q0 = zeros(20,1);
% ===== Frame alignment (核心) =====
% xR: 机器人执行器指令坐标（你控制的五自由度）
% xC: 镜片相对腔体{C}的小位姿（你仿真“真实状态”）
% xO: 光学模型内部坐标（build cavity geometry / optics线性化坐标）

% 1) 夹爪/装夹导致的常值偏置：机器人坐标到镜片相对腔体坐标
%    若你希望“你控制的就是xC”，先设为0
P.xR_to_C_bias = zeros(5,1);      % [um;um;um;rad;rad]

% 2) 腔体几何坐标{C}到光学模型坐标{O}的偏置
%    你希望几何与光学重合 => 必须设为0
P.align_optical_geo = true;
if P.align_optical_geo
    P.xC_to_opt_bias = zeros(5,1);    % [um;um;um;rad;rad]
    % 3) loss_true 里如果还用了 dx_bias（mm/rad），要清零，否则最优点不会在光学零位
    P.dx_bias = zeros(4,1);           % [mm; mm; rad; rad]
    P.optics_q0 = zeros(20,1);
end



% ---- normalization radii on aperture plane (theoretical) ----
P.aperture_ry_mm = 2.0;
P.aperture_rz_mm = 2.0;
P.e_enter_cont   = 0.25;   % 更严格进入精调阈值
P.Juv_step_um    = 1.0;
P.Juv_lambda     = 1e-3;
% ===== Loss fine (quadratic fit) =====
P.loss_fit_step_um = 0.8;      % 采样半径 Δ（um），0.5~2um 之间调
P.loss_fit_maxstep_um = 1.0;   % 一步最多走多远（um）
P.e_loss_max = 0.8;            % 损耗阶段允许的 ||e|| 上限（无量纲）
P.loss_avg_n = 3;              % 每个点测几次取均值（噪声大就3~5）
P.qfit_ridge = 1e-6;           % 拟合/求解时的数值正则
P.qfit_min_improve = 0;        % 要求比当前点至少降低多少(ppm)，可先0




% ---------- Step sizes (initial; can tune) ----------
P.dz_approach = +0.5; % um/step, z-up positive, slower approach
P.dz_backoff  = +5.0;   % um
P.duv_coarse  = 2.0;    % um per step
P.duv_fine    = 0.5;    % um per step
% ---------- Loss model scales (set by "how sensitive" loss is) ----------
% ---------- Loss model in ppm ----------
% Threshold: 0.12% = 1200 ppm
P.L_thresh_ppm = 1200;

% Best achievable loss floor (ppm) near theoretical optimum
P.L_min = 200;            % ppm (你可以按实际希望的最好水平改，比如 100~300)

% Tilt-dominant improvement: tilt decreases -> loss decreases
P.L_theta_gain = 3000;    % ppm (倾角导致的最大增量，先给 3k)
P.theta_ref = deg2rad(0.02); % rad (0.02度为“明显变差”的尺度)

% Secondary effects: lateral / z errors
P.L_uv_gain = 800;        % ppm (先小一些，表示平移次要)
P.uv_ref = 30;            % um

P.L_z_gain = 350;         % ppm
P.z_ref = 30;             % um


P.loss_abs_sigma_ppm = 20;  % ppm, 拟合的绝对噪声底
% ===== Geometry (annulus contact domain) =====
P.Ro_um = 11000;   % outer radius, 22mm OD => 11mm
P.Ri_um = 9000;    % inner radius, 18mm ID => 9mm
P.R_disk_um = P.Ro_um; % keep legacy field for compatibility

% ===== Contact "foundation" stiffness (Winkler bed) =====
% pressure p = k_w * indentation (N/um^3)  [因为 p(N/um^2)=k_w(N/um^3)*w(um)]
P.k_w = 1e-9;
P.g0_um = 0;

% ===== Compliance / quasi-static equilibrium =====
% Fz = kz*(z-zr), Mx = kthx*(thx-thxr), My = kthy*(thy-thyr)
P.c_z = 2.0;            % um/N
P.c_th = 5e-8;          % rad/(N*um)
P.k_z = 1/max(P.c_z,1e-12);        % N/um
P.k_thx = 1/max(P.c_th,1e-18);     % N*um/rad
P.k_thy = P.k_thx;
P.alpha_def = 0.25;     % deformation 1st-order update rate (0~1)
P.eq_max_iter = 20;
P.eq_relax = [0.6; 0.6; 0.6];
P.eq_tol = [1e-4; 1e-7; 1e-7];

% weak uv -> tilt coupling through gripper flexibility [rad/um]
P.K_uv_to_th = [0, 2e-9; -2e-9, 0];

% ===== Friction (stick-slip + Coulomb limit) =====
P.friction_mode = "distributed"; % "distributed"(partial slip) or "lumped"(legacy)
P.mu = 0.20;
P.k_t = 0.30 * P.k_w;   % N/um^3
P.alpha_fx2fz = 0.0;    % sensor-axis projection coefficients
P.alpha_fy2fz = 0.0;

% ===== Numerical integration grid =====
P.grid_step_um = 300;   % 网格间距，越小越精细但更慢，300~500um够用




% ---------- Random seed ----------
P.seed = 42;
% Precompute annulus contact grid points and area weights
dx = P.grid_step_um;
xs = -P.Ro_um:dx:P.Ro_um;
ys = xs;
[Xg, Yg] = meshgrid(xs, ys);
r2 = Xg.^2 + Yg.^2;
mask_ann = (r2 <= P.Ro_um^2) & (r2 >= P.Ri_um^2);
P.contact_x = Xg(mask_ann);      % column vectors
P.contact_y = Yg(mask_ann);
P.contact_dA = dx*dx;            % each cell area (um^2)

% Legacy aliases (force_model fallback / compatibility)
P.disk_x = P.contact_x;
P.disk_y = P.contact_y;
P.disk_dA = P.contact_dA;

% ---------- Final alignment guard (keep optical & geometric zero aligned) ----------
if isfield(P,'align_optical_geo') && P.align_optical_geo
    P.xC_to_opt_bias = zeros(5,1);
    P.dx_bias = zeros(4,1);
    P.optics_q0 = zeros(20,1);
end



end
