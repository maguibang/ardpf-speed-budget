function S = generalized_A()
%% generalized_A.m —— 把 A_closedform 从"正多边形"推广到"变段长/变转角多边形"
%
%  正多边形：Δθ 与 T_k 都被 M 锁死，故 A 可写成 A(M)。
%  变段长多边形：逐拐角的 Δθ_m、逐段 T_k,m 各自不同，没有单一 A；
%                正确量是"一整圈复合映射"的定点峰值。
%
%  复合映射（一圈 = nv 个拐角 + nv 段直线）：
%      e <- Rot(-Δθ_m) e + c_m       拐角注入 c_m = (Rot(-Δθ_m)-I)d
%      e <- seg_relax(e, T_k,m)      段内抛物线弛豫（解析，与 A_closedform 相同）
%  稳态 = 该复合映射的定点（迭代多圈收敛）。
%
%  本脚本在 randpoly（7 顶点、随机转角/半径，段长不一）上验证：
%  复合映射定点峰值应复现全仿真 max e_form，从而给 randpoly 一个【独立】的 C。

Rc = 15*sqrt(2); v = 1.0; N = 10; ds = 5; dt = 0.1;
phi = pi/4;  d = zeros(N-1,2);
for i = 1:N-1
    r = ceil(i/2);  s = 2*mod(i,2)-1;
    d(i,:) = [-r*ds*cos(phi), s*r*ds*sin(phi)];
end
b0 = 1/sqrt(2);  % baseline b0 = v_max/sqrt(K_att), v_max = 1

% ---- randpoly 的 7 个顶点（复刻 ardpf_run/traj_table 的 randpoly）----
rng(1,'twister');  nv = 7;
ang = sort(2*pi*rand(1,nv));  rad = Rc*(0.55 + 0.45*rand(1,nv));
vx = rad.*cos(ang);  vy = rad.*sin(ang);

% ---- 逐边方向、逐拐角转角、逐段长 ----
edge_dir = zeros(1,nv);
for k = 1:nv
    k2 = mod(k,nv)+1;
    edge_dir(k) = atan2(vy(k2)-vy(k), vx(k2)-vx(k));
end
dth_m = zeros(1,nv);   % 拐角 m 的航向变化（从边 m-1 到边 m）
for m = 1:nv
    dth_m(m) = atan2(sin(edge_dir(m)-edge_dir(mod(m-2,nv)+1)), ...
                     cos(edge_dir(m)-edge_dir(mod(m-2,nv)+1)));
end
Tk_m = zeros(1,nv);    % 段 m 的长度/速度
for m = 1:nv
    m2 = mod(m,nv)+1;
    Tk_m(m) = hypot(vx(m2)-vx(m), vy(m2)-vy(m))/v;
end

fprintf('\n====================================================================\n');
fprintf('  randpoly 顶点: nv=%d, 各段长 T_k (m/s):\n', nv);
fprintf('    '); fprintf('%.2f ', Tk_m); fprintf('\n');
fprintf('  各拐角 Δθ (deg):\n    '); fprintf('%.1f ', dth_m*180/pi); fprintf('\n');
fprintf('====================================================================\n');

% ---- 逐跟随者：复合映射定点 ----
peak = zeros(N-1,1);
for i = 1:N-1
    peak(i) = follower_peak_irreg(d(i,:), dth_m, Tk_m);
end
C = mean(peak) - b0;      % 拐角贡献（独立计算，非余项）

% ---- 对照：全仿真 ----
R = ardpf_run(N, ds, 0, 1, 'randpoly', 0, 1.0);
t = (0:numel(R.efrom_t)-1)'*dt;
full_steady = max(R.efrom_t(t >= 0.8*t(end)));

fprintf('\n  复合映射定点: 逐跟随者峰值 = '); fprintf('%.2f ', peak); fprintf('\n');
fprintf('  聚合 C = mean(peak) - b0 = %.3f - %.3f = %.3f\n', mean(peak), b0, C);
fprintf('  复合映射预测 max e_form = b0 + C = %.3f\n', b0 + C);
fprintf('  全仿真稳态  max e_form            = %.3f\n', full_steady);
fprintf('  偏差 = %.2f%%\n', 100*(b0+C - full_steady)/full_steady);
fprintf('\n  => 变段长多边形无需"单一 A"，复合映射定点即给独立拐角贡献。\n');
fprintf('====================================================================\n\n');
S = struct('C',C,'pred',b0+C,'full',full_steady,'peak',peak);
end


function pk = follower_peak_irreg(d, dth_m, Tk_m)
% 一整圈复合映射的定点：e_{m+1} = seg_relax(Rot(-Δθ_m) e_m + c_m, T_k,m)
nv = numel(dth_m);
e = [0;0];
for lap = 1:300
    for m = 1:nv
        c = (Rot(-dth_m(m)) - eye(2)) * d(:);
        e = Rot(-dth_m(m))*e + c;      % 拐角注入
        e = seg_relax(e, Tk_m(m));     % 段内弛豫
    end
end
% 峰值 = 一圈内拐角后的最大 |e|
pk = 0;
for m = 1:nv
    c = (Rot(-dth_m(m)) - eye(2)) * d(:);
    pk = max(pk, norm(Rot(-dth_m(m))*e + c));
end
end

function w = seg_relax(v, Tk)
P = (norm(v) - v(1))/2;
if P < 1e-12, w = [0;0]; return; end
y0 = abs(v(2));                 % 横向坐标 y=|e_y|
y1 = relax(P, y0, Tk);
w = [y1^2/(4*P) - P; sign(v(2))*y1];
end

function y1 = relax(P, y0, Tk)
if y0 < 1e-12, y1 = 0; return; end
f = @(y) (y0^2 - y^2)/(8*P) + P*log(y0/y) - Tk;
y_hi = y0;  y_lo = y0/2;
while f(y_lo) <= 0 && y_lo > y0*1e-12, y_lo = y_lo/2; end
if f(y_lo) <= 0, y1 = 0; else, y1 = fzero(f, [y_lo, y_hi]); end
end

function Rm = Rot(a)
Rm = [cos(a), -sin(a); sin(a), cos(a)];
end
