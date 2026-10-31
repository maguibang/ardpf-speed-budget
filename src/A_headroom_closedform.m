function [A, e_form] = A_headroom_closedform(M, N, dsep, vmax, Katt)
%% A_headroom_closedform.m —— 速度余量下累积因子 A(delta) 的闭式（附录 D 交付物）
%
%  【机理】v_max=v 时（论文 §6.2 / A_closedform），段内守恒 P=(|e|-e_x)/2，
%        误差沿抛物线 y²=4P(x+P) 弛豫，P 不衰减，拐角注入逐圈累积，A>1。
%
%        v_max>v 时（余量 delta=v_max-v>0），段内 P 不再守恒：
%           dP/dt = -delta·P/rho ,  dq/dt = -(2+delta)·q/rho ,  rho=P+q
%        守恒量推广为  q/P^{1+2/delta} = const（delta->0 退化为 P=const）。
%        弛豫时间推广为  T = (P0-P1)/delta + (q0-q1)/(2+delta)（式 D.3）。
%        且误差在有限时间 T_origin = P0/delta + q0/(2+delta) 内到原点，
%        故注入可在拐角间被清空，A 降到 1 以下。
%
%  【输入】M 边数, N 无人机数, dsep 机间距, vmax 跟随者速度上限, Katt 吸引增益
%  【输出】A = 聚合累积因子；e_form = b0 + A·kappa·mean|d|（方形轨迹 R=0）
%
%  【验证】M=4, N=10, dsep=5, Katt=2 时，vmax/v=1.0~2.0 与 Table 13 残差 <=0.05 m
%         （<=0.2%），A=1 交叉于 vmax/v=1.3977。

if nargin < 5, Katt = 2; end
v = 1.0;
Rc = 15*sqrt(2);
dth = 2*pi/M;  Tk = 2*Rc*sin(pi/M)/v;  kap = 2*sin(pi/M);
delta = vmax - v;
b0 = vmax/sqrt(Katt);
% mean|d_i| over the N-1 followers（奇偶分支，同 A_closedform）
if mod(N,2) == 0
    R = N/2;  md = dsep * R^2 / (2*R - 1);
else
    md = dsep * (N+1)/4;
end
phi = pi/4;  d = zeros(N-1,2);
for i = 1:N-1
    r = ceil(i/2);  s = 2*mod(i,2)-1;
    d(i,:) = [-r*dsep*cos(phi), s*r*dsep*sin(phi)];
end
kinj = kap*md;
peak = zeros(N-1,1);
for i = 1:N-1
    peak(i) = follower_peak_delta(d(i,:), dth, Tk, delta);
end
A = (mean(peak) - b0)/kinj;
e_form = b0 + A*kinj;      % 方形轨迹 R=0
end


function pk = follower_peak_delta(d, dth, Tk, delta)
% 单跟随者稳态拐角峰值 = 定点 |v|，段内弛豫用 delta 推广的闭式
c = (Rot(-dth) - eye(2)) * d(:);
v = c;
for it = 1:400
    v = Rot(-dth)*v + c;               % 拐角映射
    v = seg_relax_delta(v, Tk, delta); % 广义段内弛豫
end
pk = norm(Rot(-dth)*v + c);
end


function w = seg_relax_delta(v, Tk, delta)
% 段内弛豫 v->w，沿 q/P^{1+2/delta}=const 走 Tk 时间（delta>0）。
% delta=0 时退化为 A_closedform 的 seg_relax（P 守恒，抛物线）。
rho = norm(v);
if rho < 1e-12, w = [0;0]; return; end
P0 = (rho - v(1))/2;  q0 = (rho + v(1))/2;
if delta < 1e-9
    P = P0;  if P < 1e-12, w = [0;0]; return; end
    y0 = abs(v(2));                     % 横向坐标 y=|e_y|
    y1 = relax0(P, y0, Tk);
    w = [y1^2/(4*P) - P; sign(v(2))*y1];
    return;
end
C = q0 / max(P0,1e-300)^(1+2/delta);    % 守恒量（式 D.2）
T_origin = P0/delta + q0/(2+delta);     % 有限清空时间（式 D.4）
if Tk >= T_origin
    w = [0;0];                            % 注入被完全清空
    return;
end
% 解 T(P1)=Tk（式 D.3），T 关于 P1 单调减，唯一根 P1∈(0,P0)
f = @(P1) (P0-P1)/delta + (q0 - C*P1^(1+2/delta))/(2+delta) - Tk;
P1 = fzero(f, [P0*1e-12, P0]);
q1 = C*P1^(1+2/delta);
w = [q1 - P1; sign(v(2))*2*sqrt(P1*q1)];
end


function y1 = relax0(P, y0, Tk)
% delta=0 的弛豫（同 A_closedform 式 8，y=|e_y|）
if y0 < 1e-12, y1 = 0; return; end
f = @(y) (y0^2 - y^2)/(8*P) + P*log(y0/y) - Tk;
y_hi = y0;  y_lo = y0/2;
while f(y_lo) <= 0 && y_lo > y0*1e-12, y_lo = y_lo/2; end
if f(y_lo) <= 0, y1 = 0; else, y1 = fzero(f, [y_lo, y_hi]); end
end


function Rm = Rot(a)
Rm = [cos(a), -sin(a); sin(a), cos(a)];
end
