function R = ardpf_run(N, dsep, taumax_s, nseed, traj, polyM, vsat, repform, phi_deg, mixrho, Katt, Kff, tau_v)
% traj: 'square'（论文原设）| 'circle'（航向平滑，用于检验拐角伪影）
% ARDPF_RUN  最小 ARD-PF + 均匀时延复现（对标 Drones 10(1):34）
%
%   R = ardpf_run(N, dsep, taumax_s, nseed)
%
%   严格按论文实现：
%     式(7)-(9)   逐链路均匀整数时延 tau_ij[k] ~ U{0,...,tau_max/dt}
%     式(11)(16)  吸引项：朝【延迟的】期望位置 p*_i = ptil_0 + R[k]d_i
%     式(12)(16)  排斥项：K_rep * sum_j (p_i - ptil_j)/(|ptil_j - p_i|^3 + eps)
%     式(17)(18)  速度积分 + 饱和：v <- sat_vmax(v + u*dt)，p <- p + v*dt
%     式(32)(35)  e_form = mean_i |(p_i-p_0) - R d_i|，>25 m 判失稳
%     式(36)      d_min < 0.5 m 判失稳
%
%   参数取自论文 Table 2 + p.11/p.12：
%     dt=0.1, v_max=1.0, K_att=K_rep=2.0, eps=1e-2
%     领航者：正方形，边长 30 m，匀速 1.0 m/s，5 圈

P.dt      = 0.1;
P.vmax    = 1.0;
P.Katt    = 2.0;
P.Krep    = 2.0;
P.eps     = 1e-2;
P.N       = N;
P.dsep    = dsep;
P.taustep = round(taumax_s / P.dt);
P.side    = 30;
P.laps    = 5;
P.Emax    = 25;      % 式(35)
P.dcoll   = 0.5;     % 式(36)
if nargin < 5 || isempty(traj),  traj  = 'square'; end
if nargin < 6 || isempty(polyM), polyM = 4;        end
if nargin < 7 || isempty(vsat),  vsat  = [];       end
P.traj  = traj;
P.polyM = polyM;
P.vsat  = vsat;      % [] -> 用 P.vmax（论文设定）；Inf -> 关掉速度饱和
if nargin < 8 || isempty(repform), repform = 'eq16'; end
P.repform = repform;
% P.repform: 'eq16'     -> K_rep/(r^3 + eps)      （论文式16 直接给的控制律）
%             'eq12grad' -> 2*K_rep/(r^2 + eps)^2 （式12 势能严格求梯度）
%             'eq16x10'  -> 10*K_rep/(r^3 + eps)  （【哨兵】用来验证开关真的生效）
if nargin < 11 || isempty(Katt), Katt = P.Katt; end
P.Katt = Katt;      % 吸引增益（默认 2.0）
if nargin < 12 || isempty(Kff), Kff = 0; end
P.Kff = Kff;        % 速度前馈增益（0=关；>0 把跟随者速度拉向领航者速度）
if nargin < 13 || isempty(tau_v), tau_v = 0; end
P.tau_v = tau_v;    % 一阶速度环时间常数（0=纯运动学；>0 实际速度以一阶惯性跟踪速度命令）

switch traj
    case 'square'
        Lper = 4*P.side;                 % 每圈 120 m
    case 'circle'
        P.Rc = 15*sqrt(2);               % 外接半径，与正方形同尺度
        Lper = 2*pi*P.Rc;
    case 'poly'
        % 正 M 边形，外接半径 P.Rc。M=4 且 Rc=15*sqrt(2) 时边长 30 m（= 论文的方形）
        P.Rc = 15*sqrt(2);
        P.M  = P.polyM;
        Lper = P.M * 2*P.Rc*sin(pi/P.M);
    case {'spiral','fig8','randpoly'}
        % 按弧长参数化的通用轨迹表（螺旋/八字/随机折线）
        P.Ttab = traj_table(traj, 15*sqrt(2), 1);
        Lper = P.Ttab.L;
    case 'mixpoly'
        % 混合轨迹：2 尖角 + 2 圆角（rho 由第 10 个参数给）
        if isempty(mixrho), mixrho = P.side*0.4; end
        P.Ttab = mixpoly_table(P.side, mixrho);
        P.mixrho = mixrho;
        Lper = P.Ttab.L;
    otherwise
        error('unknown traj');
end
P.Lper = Lper;
NT = round(P.laps * Lper / P.vmax / P.dt);
P.NT = NT;

if nargin < 9  || isempty(phi_deg), phi_deg = 45; end
if nargin < 10, mixrho = []; end
d0 = v_offsets(N, dsep, phi_deg);

R.unstable  = zeros(nseed,1);
R.tfail     = nan(nseed,1);
R.efrom_max = zeros(nseed,1);
R.dmin_min  = zeros(nseed,1);
R.N = N;  R.dsep = dsep;  R.taumax = taumax_s;

for s = 1:nseed
    rng(s, 'twister');
    [un, tf, em, dm, eft] = one_run(P, d0);
    R.unstable(s)  = un;
    R.tfail(s)     = tf;
    R.efrom_max(s) = em;
    R.dmin_min(s)  = dm;
    if s == 1, R.efrom_t = eft; end     % 仅存 seed 1 的时间序列，供诊断
end
R.stabRate = 1 - mean(R.unstable);
R.name = sprintf('N=%d, dsep=%g m, taumax=%g s', N, dsep, taumax_s);
end


%% ================= 单次仿真 =================
function [unstable, tfail, eform_max, dmin_min, ef_t] = one_run(P, d0)

N = P.N;  NT = P.NT;  dt = P.dt;
ef_t = zeros(NT+1,1);

Hx = zeros(NT+1, N);   Hy = zeros(NT+1, N);   % 历史缓冲

[p0, th0] = leader_state(0, P);
Rm0 = rotmat(th0);
Ppos = zeros(N,2);
Ppos(1,:) = p0;
for i = 2:N
    Ppos(i,:) = p0 + (Rm0 * d0(i,:)')';
end
Pvel = zeros(N,2);

Hx(1,:) = Ppos(:,1)';  Hy(1,:) = Ppos(:,2)';

unstable = 0;  tfail = NaN;  eform_max = 0;  dmin_min = inf;

for k = 1:NT
    t = k*dt;
    [p0k, th0k] = leader_state(t, P);
    Rm = rotmat(th0k);

    % ---- 式(7)-(9)：逐链路时延 ----
    tau = randi([0 P.taustep], N, N);
    tau(1:N+1:end) = 0;
    idx = max(k - tau, 1);
    col = repmat(1:N, N, 1);
    lin = sub2ind([NT+1, N], idx, col);
    PtilX = reshape(Hx(lin), N, N);
    PtilY = reshape(Hy(lin), N, N);

    % ---- 式(16) ARD-PF ----
    U = zeros(N,2);
    for i = 1:N
        plead = [PtilX(i,1), PtilY(i,1)];
        pstar = plead + (Rm * d0(i,:)')';
        uatt  = P.Katt * (pstar - Ppos(i,:));

        rel = [Ppos(i,1) - PtilX(i,:)', Ppos(i,2) - PtilY(i,:)'];
        dn  = sqrt(sum(rel.^2, 2));
        dn(i) = inf;
        switch P.repform
            case 'eq16'                        % 式(16) 直接给的控制律
                w = 1 ./ (dn.^3 + P.eps);
            case 'eq12grad'                    % 式(12) 势能严格求梯度
                w = 2 ./ (dn.^2 + P.eps).^2;
            case 'eq16x10'                     % 哨兵：排斥力放大 10 倍
                w = 10 ./ (dn.^3 + P.eps);
            otherwise
                error('unknown repform');
        end
        w(i) = 0;
        urep = P.Krep * sum(rel .* w, 1);

        U(i,:) = uatt + urep;
    end

    % ---- 式(17)(18) + 一阶速度环 + 前馈 ----
    Vcmd = Pvel + U*dt;                     % 速度命令 = 加速度积分
    if P.Kff > 0
        vlead = [cos(th0k), sin(th0k)];     % 领航者速度 v=1 沿航向
        Vcmd(2:end,:) = Vcmd(2:end,:) + P.Kff * (vlead - Vcmd(2:end,:)) * dt;
    end
    vs   = P.vmax;
    if isfield(P,'vsat') && ~isempty(P.vsat), vs = P.vsat; end
    if isfinite(vs)
        sp = sqrt(sum(Vcmd.^2, 2));
        ov = sp > vs;
        Vcmd(ov,:) = Vcmd(ov,:) .* (vs ./ sp(ov));   % 饱和作用在速度命令上
    end
    if P.tau_v > 0
        Pvel = Pvel + (Vcmd - Pvel) * (P.dt / P.tau_v);   % 一阶速度环（旋翼动力学）
    else
        Pvel = Vcmd;                                        % 纯运动学
    end
    Ppos = Ppos + Pvel*dt;
    Ppos(1,:) = p0k;

    Hx(k+1,:) = Ppos(:,1)';  Hy(k+1,:) = Ppos(:,2)';

    % ---- 式(32)(35)(36) ----
    ef = 0;
    for i = 2:N
        ef = ef + norm((Ppos(i,:) - Ppos(1,:)) - (Rm*d0(i,:)')');
    end
    ef = ef/(N-1);
    D = squareform(pdist(Ppos));
    D(1:N+1:end) = inf;
    dm = min(D(:));

    eform_max = max(eform_max, ef);
    dmin_min  = min(dmin_min, dm);
    ef_t(k+1) = ef;

    if (ef > P.Emax || dm < P.dcoll) && ~unstable
        unstable = 1;  tfail = t;
    end
end
end


%% ================= 辅助 =================
function [p, th] = leader_state(t, P)
L = P.Lper;  v = P.vmax;
s = mod(v*t, L);
switch P.traj
    case 'square'
        S = P.side;  c = [-S/2, -S/2];
        if s < S
            p = c + [s, 0];            th = 0;
        elseif s < 2*S
            p = c + [S, s-S];          th = pi/2;
        elseif s < 3*S
            p = c + [S-(s-2*S), S];    th = pi;
        else
            p = c + [0, S-(s-3*S)];    th = -pi/2;
        end
    case 'circle'
        ph = s / P.Rc;
        p  = P.Rc * [cos(ph), sin(ph)];
        th = ph + pi/2;
    case 'poly'
        M = P.M;  Rc = P.Rc;
        side = 2*Rc*sin(pi/M);
        kk = floor(s/side);
        uu = s - kk*side;
        a1 = 2*pi*kk/M;   a2 = 2*pi*(kk+1)/M;
        v1 = Rc*[cos(a1), sin(a1)];
        v2 = Rc*[cos(a2), sin(a2)];
        p  = v1 + (v2 - v1)*(uu/side);
        d  = v2 - v1;
        th = atan2(d(2), d(1));
    case {'spiral','fig8','randpoly','mixpoly'}
        T = P.Ttab;
        p  = [interp1(T.s, T.x,  s, 'linear', 'extrap'), ...
              interp1(T.s, T.y,  s, 'linear', 'extrap')];
        th = interp1(T.s, T.th, s, 'linear', 'extrap');
end
end


function T = traj_table(kind, Rc, seed)
% 构造按弧长参数化的闭合轨迹表
    switch kind
        case 'spiral'
            % 平滑花瓣状闭合曲线（曲率连续变化，无拐角）
            th = linspace(0, 2*pi, 6000);
            r  = Rc*(1 + 0.45*cos(3*th));
            x  = r.*cos(th);   y = r.*sin(th);
        case 'fig8'
            % 双纽线（Gerono）：有自交、曲率连续
            u  = linspace(0, 2*pi, 6000);
            x  = Rc*sin(u);
            y  = Rc*sin(u).*cos(u);
        case 'randpoly'
            % 随机闭合折线：顶点处航向【跳变】
            rng(seed,'twister');
            nv  = 7;
            ang = sort(2*pi*rand(1,nv));
            rad = Rc*(0.55 + 0.45*rand(1,nv));
            vx  = rad.*cos(ang);   vy = rad.*sin(ang);
            x = [];  y = [];
            for k = 1:nv
                k2 = mod(k,nv)+1;
                x = [x, linspace(vx(k), vx(k2), 400)];
                y = [y, linspace(vy(k), vy(k2), 400)];
            end
            x(end+1) = x(1);   y(end+1) = y(1);
    end
    x = x(:)';  y = y(:)';
    % 去掉重复点，否则弧长网格非唯一，interp1 会报错
    keep = [true, hypot(diff(x), diff(y)) > 1e-9];
    x = x(keep);  y = y(keep);
    ds = hypot(diff(x), diff(y));
    s  = [0, cumsum(ds)];
    T.s = s;  T.x = x;  T.y = y;
    % 切线方向（中心差分），unwrap 后插值
    dx = gradient(x, s);  dy = gradient(y, s);
    T.th = unwrap(atan2(dy, dx));
    T.L  = s(end);
end

function Rm = rotmat(th)
Rm = [cos(th), -sin(th); sin(th), cos(th)];
end

function d = v_offsets(N, dsep, phi_deg)
% V 形偏移。phi_deg = 半张角（度）。
%   d_i = (-r*dsep*cos(phi), ±r*dsep*sin(phi))
%   注意 |d_i| = r*dsep 【与 phi 无关】—— phi 只改变"沿航向 vs 侧向"的分配
if nargin < 3 || isempty(phi_deg), phi_deg = 45; end
phi = phi_deg*pi/180;
d = zeros(N,2);
for i = 1:N-1
    r = ceil(i/2);
    s = 2*mod(i,2) - 1;
    d(i+1,:) = [-r*dsep*cos(phi), s*r*dsep*sin(phi)];
end
end
