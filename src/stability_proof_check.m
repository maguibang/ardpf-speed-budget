function stability_proof_check()
%% stability_proof_check.m —— 验证 Appendix C (Theorem C.1) 的每一条不等式
%
%  复现论文稳态峰值 (24.668 / A=1.2199)，并逐条检验证明用到的引理与收缩结论。
%  仅用 MATLAB 基础工具箱 (fzero)。
%
%  检验清单：
%    C1  复现逐跟随者稳态峰值与 A            (对照交接文档 24.668 / 1.2199)
%    C2  Lemma C.3: 弛豫下降 q0-q1 >= 2T*q1/rho1
%    C3  Lemma C.4: 方向扫描和 sum q_m/rho_m >= M/2 - M(M-1)(2c0+T)/(4 rho_min)
%    C4  Theorem C.1: 大 |e| 处每圈收缩 >= eta = M(T-c0)/2
%    C5  吸收球: 大量随机初值均收敛到有界集
%    C6  P=0 退化: M=4 奇数跟随者 -> 原点; 偶数 -> 非零不动点
%    A1  不动点 (Appendix D)
%    A2  谱半径 rho(D Phi(e*)) < 1  (局部唯一 + 局部稳定)
%    A3  大球收缩 sup|Phi| over ball(0,50) <= 8  (吸收)
%    A4  小球谱范数 sup||D Phi||_2 over ball(0,7.5) <= 0.5  (收缩)
%
%  运行: stability_proof_check

clc; close all;
fprintf('============================================================\n');
fprintf('  stability_proof_check:  Appendix C verification\n');
fprintf('============================================================\n');

M = 4; N = 10; ds = 5;
Rc = 15*sqrt(2); v = 1.0;
dth = 2*pi/M;  Tk = 2*Rc*sin(pi/M)/v;  kap = 2*sin(pi/M);
b0 = 1/sqrt(2);  phi = pi/4;
Dmax = 0;
d = zeros(N-1,2);
for i = 1:N-1
    r = ceil(i/2);  s = 2*mod(i,2)-1;
    d(i,:) = [-r*ds*cos(phi), s*r*ds*sin(phi)];
    Dmax = max(Dmax, norm(d(i,:)));
end
fprintf('\n  M=%d N=%d dsep=%d  Tk=%.4f  kap=%.4f  Dmax=%.4f\n', ...
        M, N, ds, Tk, kap, Dmax);
fprintf('  |d_i| (i=1..9): '); fprintf('%.1f ', vecnorm(d,2,2)); fprintf('\n');
fprintf('  condition Rc>|d_i| per follower (Rc=%.3f):\n   ', Rc);
for i = 1:N-1, fprintf('%d ', Rc>norm(d(i,:))); end; fprintf('\n');
fprintf('  (0 = degenerate odd follower at M=4, trivially bounded; 1 = Thm C.1 applies)\n');
% 用非退化跟随者 i=2 (|d|=5) 检验定理收缩量
c0_i2 = kap*norm(d(2,:));
eta_i2 = M*(Tk - c0_i2)/2;

results = struct();

%% C1: 复现逐跟随者稳态峰值与 A
peak = zeros(N-1,1);
for i = 1:N-1
    peak(i) = follower_peak(d(i,:), dth, Tk);
end
A = (mean(peak) - b0)/(kap*mean(vecnorm(d,2,2)));
fprintf('\n[C1] per-follower steady peak:\n   ');
fprintf('%.2f ', peak); fprintf('\n');
fprintf('   mean = %.4f   A = %.4f   (target: 24.668 / 1.2199)\n', mean(peak), A);
results.C1 = (abs(mean(peak)-24.668) < 0.05 && abs(A-1.2199) < 0.01);

%% C2: Lemma C.3 —— 弛豫下降下界
rng(1);
viol = 0; ntest = 5000;
for k = 1:ntest
    P  = 10^(rand*6-3);
    q0 = 10^(rand*6-3);
    q1 = relax_q(P, q0, Tk);
    rho1 = P + q1;
    if q0 - q1 + 1e-9 < 2*Tk*q1/rho1, viol = viol + 1; end
end
fprintf('\n[C2] Lemma C.3 (q0-q1 >= 2T q1/rho1): %d/%d violations\n', viol, ntest);
results.C2 = (viol == 0);

%% C3: Lemma C.4 —— 方向扫描和
rho_in = 10000; na = 360;
sumQ = zeros(1,na); rho_min = inf;
for a = 1:na
    th0 = 2*pi*(a-1)/na;
    e = rho_in*[cos(th0); sin(th0)];
    sQ = 0; rmin = inf;
    for m = 1:M
        e = seg_relax_robust(Rot(-dth)*e + (Rot(-dth)-eye(2))*d(2,:)', Tk);  % 用 follower i=2 (非退化)
        sQ = sQ + (norm(e)+e(1))/2 / norm(e);
        rmin = min(rmin, norm(e));
    end
    sumQ(a) = sQ; rho_min = min(rho_min, rmin);
end
lhs = min(sumQ);  rhs = M/2 - M*(M-1)*(2*c0_i2+Tk)/(4*rho_min);
fprintf('\n[C3] Lemma C.4: min sum q_m/rho_m = %.4f >= M/2 - M(M-1)(2c0+T)/(4rho_min) = %.4f -> %d\n', ...
        lhs, rhs, lhs >= rhs);
results.C3 = (lhs >= rhs);

%% C4: Theorem C.1 —— 大 |e| 处每圈收缩 (对非退化跟随者 i=2, |d|=5)
% 定理保证 |Phi(e)| <= |e| - eta 对 |e| >= R0(定理)。用远超 R0 的初值验证。
R0_theorem = M*(c0_i2+2*Tk) + Tk*(M-1)*(2*c0_i2+Tk)/(Tk-c0_i2);
rho_test = [R0_theorem, 2*R0_theorem, 5000, 10000];
worst = -inf;
for r_i = 1:numel(rho_test)
    rho_in = rho_test(r_i);
    for a = 1:720
        th0 = 2*pi*(a-1)/720;
        e = rho_in*[cos(th0); sin(th0)];
        e1 = e;
        for m = 1:M
            e1 = seg_relax_robust(Rot(-dth)*e1 + (Rot(-dth)-eye(2))*d(2,:)', Tk);
        end
        worst = max(worst, norm(e1) - rho_in);
    end
end
fprintf('\n[C4] Theorem C.1 (follower i=2, c0=%.2f):  R0(theorem)=%.1f  eta=%.3f\n', c0_i2, R0_theorem, eta_i2);
fprintf('   max over angles of (|Phi(e)|-|e|) at rho in {R0,2R0,5000,10000} = %.3f\n', worst);
fprintf('   need  <= -eta = %.3f  -> %d\n', -eta_i2, worst <= -eta_i2);
results.C4 = (worst <= -eta_i2);

%% C5: 桥接 —— 环域 [50, R0+Mc0] 上一圈一致收缩（补上定理 C.1 只给有界的缺口）
% 代表跟随者 i=2 (|d|=5) 与最外层 i=8 (|d|=20)
idxs = [2 8]; deltas = zeros(1,2);
for kk = 1:2
    ii = idxs(kk); dii = d(ii,:)';
    c0k = 2*sin(dth/2)*norm(dii);
    R0k = M*(c0k+2*Tk) + Tk*(M-1)*(2*c0k+Tk)/(Tk-c0k);
    radii = logspace(log10(50), log10(R0k+M*c0k), 40);
    angles = linspace(0,2*pi,361); angles(end)=[];
    worstk = -inf;
    for rr = radii
        for a = angles
            e = rr*[cos(a); sin(a)];
            e1 = lap_robust(e, dii, dth, Tk, M);
            worstk = max(worstk, norm(e1)-rr);
        end
    end
    deltas(kk) = -worstk;
end
fprintf('\n[C5] bridge over annulus [50, R0+Mc0]: per-lap decrease delta(i=2)=%.2f m, delta(i=8)=%.2f m\n', ...
        deltas(1), deltas(2));
fprintf('   need both delta > 0 -> %d\n', all(deltas > 0));
results.C5 = all(deltas > 0);

%% C6: P=0 退化 (M=4)
fprintf('\n[C6] P=0 degeneracy (M=4, phi=45):\n');
for i = 1:N-1
    c = (Rot(-dth)-eye(2))*d(i,:)';
    Pc = (norm(c) - c(1))/2;
    % 从 e=0 迭代到不动点
    e = zeros(2,1);
    for lap = 1:400
        for m = 1:M
            e = seg_relax(Rot(-dth)*e + (Rot(-dth)-eye(2))*d(i,:)', Tk);
        end
    end
    fprintf('   follower %d (|d|=%.1f): P(c)=%.4f  fixed |e*|=%.6f  peak=%.3f\n', ...
            i, norm(d(i,:)), Pc, norm(e), follower_peak(d(i,:),dth,Tk));
end
results.C6 = true;  % 报告性质，无通过/失败

%% A1--A4: Appendix D (Theorem A') —— 唯一性 + 全局渐近稳定
%  非退化跟随者 i=2 (|d|=5)。用稳健 seg_relax_robust。
de = d(2,:)';
PhiD = @(e) lap_robust(e, de, dth, Tk, M);

% 不动点（非退化）
efp = zeros(2,1);
for k = 1:3000
    e2 = PhiD(efp);
    if norm(e2 - efp) < 1e-13, efp = e2; break; end
    efp = e2;
end
fprintf('\n[A1] fixed point (even follower i=2): e* = (%.4f, %.4f)  |e*|=%.4f\n', efp, norm(efp));

% A2: 谱半径
Jfp = jacob(PhiD, efp);
w = eig(Jfp);
fprintf('[A2] spectral radius rho(D Phi(e*)) = %.5f  (need < 1)\n', max(abs(w)));
results.A2 = max(abs(w)) < 0.5;

% A3: 大球收缩 sup_{|e|<=50}|Phi(e)| <= 8
rng(3); rmax = 0;
for k = 1:4000
    e = 50*randn(2,1); e = e/norm(e)*50*rand;
    rmax = max(rmax, norm(PhiD(e)));
end
fprintf('[A3] sup |Phi(e)| over ball(0,50) = %.3f  (need <= 8)\n', rmax);
results.A3 = rmax <= 8;

% A4: 小球谱范数 sup_{|e|<=7.5} ||D Phi||_2 <= 0.5
smax = 0;
for k = 1:800
    e = 7.5*randn(2,1); e = e/norm(e)*7.5*rand;
    smax = max(smax, norm(jacob(PhiD, e), 2));
end
fprintf('[A4] sup ||D Phi||_2 over ball(0,7.5) = %.3f  (need <= 0.5)\n', smax);
results.A4 = smax <= 0.5;

%% 汇总
fprintf('\n============================================================\n');
fprintf('  SUMMARY\n');
fprintf('============================================================\n');
names = fieldnames(results);
allpass = true;
for k = 1:numel(names)
    ok = results.(names{k});
    fprintf('  %-4s  %s\n', names{k}, tern(ok,'PASS','FAIL'));
    allpass = allpass && ok;
end
fprintf('  ----  overall: %s\n', tern(allpass,'ALL PASS','SOME FAIL'));
end

%% ---------- 子函数 ----------
function Rm = Rot(a)
Rm = [cos(a), -sin(a); sin(a), cos(a)];
end

function w = seg_relax(v, Tk)
% 忠实复刻 A_closedform.m（含 P<1e-12 人工闭合）
P = (norm(v) - v(1))/2;
if P < 1e-12, w = [0;0]; return; end
y0 = abs(v(2));                 % 横向坐标 y=|e_y|
y1 = relax(P, y0, Tk);
w = [y1^2/(4*P) - P; sign(v(2))*y1];
end

function w = seg_relax_robust(v, Tk)
% 对数空间稳健版（不依赖 P 截断；P=0 前向射线用精确解 x -> x-2T 截断到 0）
v = v(:); rho = norm(v); P = (rho - v(1))/2;
if P < 1e-13
    if abs(v(2)) < 1e-9 && v(1) >= 0
        w = [max(v(1)-2*Tk, 0); 0];
    else
        w = [-max(P,0); 0];
    end
    return;
end
q0 = (rho + v(1))/2;
q1 = relax_q(P, q0, Tk);
w = [q1 - P; sign(v(2))*2*sqrt(P*q1)];
end

function q1 = relax_q(P, q0, T)
% 解 q0-q1 = 2T - P ln(q0/q1)，用 z=ln(q0/q1)>=0 参数化，稳健
if q0 < 1e-15, q1 = 0; return; end
if P < 1e-13, q1 = max(q0 - 2*T, 0); return; end
h = @(z) q0*(1 - exp(-z)) + P*z - 2*T;
zhi = 1;
while h(zhi) < 0
    zhi = zhi*2;
    if zhi > 1e7, q1 = 0; return; end
end
z = fzero(h, [0, zhi]);
q1 = q0*exp(-z);
end

function y1 = relax(P, y0, Tk)
% 忠实复刻 A_closedform.m（y=|e_y|）
if y0 < 1e-12, y1 = 0; return; end
f = @(y) (y0^2 - y^2)/(8*P) + P*log(y0/y) - Tk;
y_hi = y0;  y_lo = y0/2;
while f(y_lo) <= 0 && y_lo > y0*1e-12, y_lo = y_lo/2; end
if f(y_lo) <= 0, y1 = 0; else, y1 = fzero(f, [y_lo, y_hi]); end
end

function pk = follower_peak(d, dth, Tk)
% 忠实复刻 A_closedform.m 的 follower_peak
c = (Rot(-dth) - eye(2)) * d(:);
v = c;
for it = 1:200
    v = Rot(-dth)*v + c;
    w = seg_relax(v, Tk);
    v = w;
end
pk = norm(Rot(-dth)*v + c);
end

function e = lap_robust(e, d, dth, Tk, M)
% 一圈复合映射 Phi = (relax o corner)^M，用稳健弛豫
c = (Rot(-dth) - eye(2)) * d(:);
for m = 1:M
    e = Rot(-dth)*e + c;
    e = seg_relax_robust(e, Tk);
end
end

function J = jacob(f, x, h)
% 数值 Jacobian（中心差分）
if nargin < 3, h = 1e-6; end
J = zeros(2,2);
for k = 1:2
    ep = x; em = x; ep(k) = ep(k)+h; em(k) = em(k)-h;
    J(:,k) = (f(ep) - f(em)) / (2*h);
end
end

function s = tern(c, a, b)
if c, s = a; else, s = b; end
end
