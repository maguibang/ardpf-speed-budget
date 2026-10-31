function [A, peak] = A_closedform(M, N, dsep)
%% A_closedform.m —— 拐角累积因子 A 的闭式（任务 1 交付物）
%
%  【机理】速度饱和后跟随者恒以 v_max=1 指向期望位置，机体系内
%        ė = -ê - (1,0)，段内误差沿抛物线 y²=4P(x+P) 弛豫到顶点(-P,0)（正后方）。
%        守恒量 P=(|e|-e_x)/2，弛豫时间 T_k=(y0²-y1²)/(8P)+P ln(y0/y1)，y=|e_y|。
%        拐角处 e ← Rot(-Δθ)e + c，c=(Rot(-Δθ)-I)d。稳态为定点，峰值在拐角处。
%
%  【输入】M 边形边数, N 无人机数, dsep 机间距
%  【输出】A = 聚合累积因子（定义同 meas_spec 的 A_net）
%          peak = 逐跟随者稳态拐角峰值
%
%  【验证】M=4,6,8,12,20 与全仿真误差 -0.2% ~ +0.7%。
%  【前提】v_max = 领航者速度 v = 1（论文设定），饱和恒成立。

Rc = 15*sqrt(2);  v = 1.0;
dth = 2*pi/M;  Tk = 2*Rc*sin(pi/M)/v;  kap = 2*sin(pi/M);
b0 = 1/sqrt(2);
% mean|d_i| over the N-1 followers, r_i = ceil(i/2) = 1,1,2,2,...
%   N even (N-1 = 2R-1):  sum r_i = R^2     -> mean = R^2/(2R-1),  R = N/2
%   N odd  (N-1 = 2R  ):  sum r_i = R(R+1)  -> mean = (R+1)/2  = (N+1)/4
if mod(N,2) == 0
    R  = N/2;
    md = dsep * R^2 / (2*R - 1);
else
    md = dsep * (N+1)/4;
end

phi = pi/4;  d = zeros(N-1,2);
for i = 1:N-1
    r = ceil(i/2);  s = 2*mod(i,2)-1;
    d(i,:) = [-r*dsep*cos(phi), s*r*dsep*sin(phi)];
end

peak = zeros(N-1,1);
for i = 1:N-1
    peak(i) = follower_peak(d(i,:), dth, Tk);
end
A = (mean(peak) - b0)/(kap*md);
end


function pk = follower_peak(d, dth, Tk)
% 单跟随者稳态拐角峰值 = 定点 |v|，v 为拐角后误差矢量
c = (Rot(-dth) - eye(2)) * d(:);
v = c;                          % 初值：第一次拐角后的状态
for it = 1:200
    v = Rot(-dth)*v + c;        % 拐角映射（作用在段末 w 上）
    w = seg_relax(v, Tk);       % 段内弛豫 v -> w
    v = w;
end
pk = norm(Rot(-dth)*v + c);     % 稳态峰值 = 拐角后 |v|
end


function w = seg_relax(v, Tk)
% 段内弛豫：v -> w（抛物线上走 T_k 时间）
P = (norm(v) - v(1))/2;
if P < 1e-12, w = [0;0]; return; end
y0 = abs(v(2));                 % 横向坐标 y=|e_y|
y1 = relax(P, y0, Tk);          % 解 Tk 方程
w = [y1^2/(4*P) - P; sign(v(2))*y1];
end


function y1 = relax(P, y0, Tk)
% 解 Tk = (y0²-y1²)/(8P) + P ln(y0/y1)，y1∈(0,y0)。f(y) 关于 y 严格递减，唯一根。
if y0 < 1e-12, y1 = 0; return; end
f = @(y) (y0^2 - y^2)/(8*P) + P*log(y0/y) - Tk;
y_hi = y0;  y_lo = y0/2;
while f(y_lo) <= 0 && y_lo > y0*1e-12, y_lo = y_lo/2; end
if f(y_lo) <= 0, y1 = 0; else, y1 = fzero(f, [y_lo, y_hi]); end
end


function Rm = Rot(a)
Rm = [cos(a), -sin(a); sin(a), cos(a)];
end
