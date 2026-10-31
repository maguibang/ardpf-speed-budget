function R = sim_instr(N, dsep, traj, polyM, vsat)
%% sim_instr.m —— 带仪表的 ARD-PF 复现（仅 tau=0 路径）
%
%  与 2026-09-22\ardpf_run.m 的动力学逐行一致（tau=0 时该仿真完全确定，
%  故可安全去掉时延分支），但额外输出【逐跟随者误差矢量】与【速度/饱和】信息，
%  这是验证相量累积模型 A = 1/|1 - q e^{i(psi-dth)}| 所必需的。
%
%  输出：
%    R.t       (NT+1)x1        时间
%    R.efrom   (NT+1)x1        标量 e_form（= 原仿真器定义，用于交叉校验）
%    R.E       (NT+1)x(N-1)x2  逐跟随者误差矢量 e_i = (p_i-p_0) - R d_i（世界系）
%    R.Vmag    (NT+1)xN        速度模
%    R.sat     (NT+1)xN        饱和标志
%    R.p0, R.th0               领航者状态

if nargin < 3 || isempty(traj),   traj  = 'poly'; end
if nargin < 4 || isempty(polyM),  polyM = 4;      end
if nargin < 5 || isempty(vsat),   vsat  = 1.0;    end

P.dt=0.1; P.vmax=1.0; P.Katt=2.0; P.Krep=2.0; P.eps=1e-2;
P.N=N; P.side=30; P.laps=5; P.Rc=15*sqrt(2); P.M=polyM; P.traj=traj;

switch traj
    case 'square', Lper = 4*P.side;
    case 'circle', Lper = 2*pi*P.Rc;
    case {'poly','polycw'}, Lper = P.M*2*P.Rc*sin(pi/P.M);
    otherwise,     error('sim_instr: 目前只支持 square/circle/poly');
end
P.Lper = Lper;
NT = round(P.laps*Lper/P.vmax/P.dt);
P.NT = NT;

d0 = v_offsets(N, dsep);
dt = P.dt;

t    = (0:NT)'*dt;
efrom = zeros(NT+1,1);
E     = zeros(NT+1, N-1, 2);
Vmag  = zeros(NT+1, N);
SAT   = false(NT+1, N);
P0    = zeros(NT+1,2);
TH0   = zeros(NT+1,1);

[p0, th0] = leader_state(0, P);
Rm0 = rotmat(th0);
Ppos = zeros(N,2);  Ppos(1,:) = p0;
for i = 2:N, Ppos(i,:) = p0 + (Rm0*d0(i,:)')'; end
Pvel = zeros(N,2);

P0(1,:) = p0;  TH0(1) = th0;
[efrom(1), E(1,:,:)] = form_err(Ppos, p0, Rm0, d0, N);

p0prev = p0;   % 关键：ardpf_run 的历史缓冲在更新后才写入，
               % 因此 attract 用的是【上一步】的领航者位置 p0[k-1]，
               % 而旋转阵 R[k] 却是当前的。这个不一致必须复刻。
for k = 1:NT
    tt = k*dt;
    [p0k, th0k] = leader_state(tt, P);
    Rm = rotmat(th0k);

    U = zeros(N,2);
    for i = 1:N
        pstar = p0prev + (Rm*d0(i,:)')';
        uatt  = P.Katt*(pstar - Ppos(i,:));
        rel   = Ppos(i,:) - Ppos;              % p_i - p_j（排斥方向，勿反号）
        dn    = sqrt(sum(rel.^2,2));
        dn(i) = inf;
        w = 1./(dn.^3 + P.eps);   w(i) = 0;
        urep = P.Krep*sum(rel.*w, 1);
        U(i,:) = uatt + urep;
    end

    Pvel = Pvel + U*dt;
    sp = sqrt(sum(Pvel.^2,2));
    if isfinite(vsat)
        ov = sp > vsat;
        Pvel(ov,:) = Pvel(ov,:).*(vsat./sp(ov));
    else
        ov = false(N,1);
    end
    Ppos = Ppos + Pvel*dt;
    Ppos(1,:) = p0k;
    p0prev = p0k;

    P0(k+1,:) = p0k;  TH0(k+1) = th0k;
    Vmag(k+1,:) = sqrt(sum(Pvel.^2,2));
    SAT(k+1,:)  = ov';
    [efrom(k+1), E(k+1,:,:)] = form_err(Ppos, p0k, Rm, d0, N);
end

R = struct('t',t,'efrom',efrom,'E',E,'Vmag',Vmag,'sat',SAT, ...
           'p0',P0,'th0',TH0,'d0',d0,'dt',dt,'NT',NT,'Lper',Lper);
end


function [ef, Ev] = form_err(Ppos, p0, Rm, d0, N)
Ev = zeros(1,N-1,2);  ef = 0;
for i = 2:N
    e = (Ppos(i,:) - p0) - (Rm*d0(i,:)')';
    Ev(1,i-1,:) = reshape(e,1,1,2);
    ef = ef + norm(e);
end
ef = ef/(N-1);
end

function [p, th] = leader_state(t, P)
L = P.Lper;  v = P.vmax;
s = mod(v*t, L);
switch P.traj
    case 'square'
        S = P.side;  c = [-S/2, -S/2];
        if s < S,            p = c + [s, 0];          th = 0;
        elseif s < 2*S,      p = c + [S, s-S];        th = pi/2;
        elseif s < 3*S,      p = c + [S-(s-2*S), S];  th = pi;
        else,                p = c + [0, S-(s-3*S)];  th = -pi/2;
        end
    case 'circle'
        ph = s/P.Rc;  p = P.Rc*[cos(ph), sin(ph)];  th = ph + pi/2;
    case 'poly'
        M = P.M;  Rc = P.Rc;
        side = 2*Rc*sin(pi/M);
        kk = floor(s/side);  uu = s - kk*side;
        a1 = 2*pi*kk/M;  a2 = 2*pi*(kk+1)/M;
        v1 = Rc*[cos(a1), sin(a1)];  v2 = Rc*[cos(a2), sin(a2)];
        p  = v1 + (v2-v1)*(uu/side);
        dd = v2 - v1;  th = atan2(dd(2), dd(1));
    case 'polycw'                       % 顺时针（转向相反）——用于判据归因
        M = P.M;  Rc = P.Rc;
        side = 2*Rc*sin(pi/M);
        kk = floor(s/side);  uu = s - kk*side;
        a1 = -2*pi*kk/M;  a2 = -2*pi*(kk+1)/M;
        v1 = Rc*[cos(a1), sin(a1)];  v2 = Rc*[cos(a2), sin(a2)];
        p  = v1 + (v2-v1)*(uu/side);
        dd = v2 - v1;  th = atan2(dd(2), dd(1));
end
end

function Rm = rotmat(th)
Rm = [cos(th), -sin(th); sin(th), cos(th)];
end

function d = v_offsets(N, dsep)
phi = pi/4;  d = zeros(N,2);
for i = 1:N-1
    r = ceil(i/2);  s = 2*mod(i,2)-1;
    d(i+1,:) = [-r*dsep*cos(phi), s*r*dsep*sin(phi)];
end
end
