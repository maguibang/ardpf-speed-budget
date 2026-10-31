function S = mixed_check(N, ds)
%% mixed_check.m —— 混合轨迹（2 尖角 + 2 圆角）验收
%
%  独立重算 §9 "Where the decomposition does not apply" 的混合轨迹测试
%  （两个 90° 尖角 + 两段半径 12 m 的 90° 圆弧，S=30，rho=12）。
%
%  【结论两半】
%    ① 降阶模型（统一离散映射）在混合轨迹上存活：resid ≈ -0.6%；
%    ② 混合误差超过任一纯端点（方形 25.1 m、圆 11.2 m）—— 耦合而非相加。
%
%  【输出】
%    full       全仿真 max_t e_form (m)
%    unified    统一离散映射 max_t e_form (m)
%    resid      (unified - full)/full * 100  (%)
%    full_sq    纯拐角方形（M=4）全仿真 max_t e_form (m)
%    full_circ  圆轨迹全仿真 max_t e_form (m)

dt = 0.1; v = 1.0; vmax = 1.0; laps = 5; S0 = 30; rho = 12;
phi = pi/4; d = zeros(N-1,2);
for i = 1:N-1
    r = ceil(i/2); s = 2*mod(i,2)-1;
    d(i,:) = [-r*ds*cos(phi), s*r*ds*sin(phi)];
end

% 全仿真
R = ardpf_run(N, ds, 0, 1, 'mixpoly', 0, 1.0, 'eq16', 45, rho);
full = R.efrom_max(1);

% 统一离散映射（降阶模型）
T = mixpoly_table(S0, rho, 40);
Lper = T.L;
t = 0:dt:laps*Lper/v;
s = mod(v*t, Lper);
thseq = interp1(T.s, T.th, s, 'linear', 'extrap');
NT = numel(thseq)-1;
e = zeros(N-1,2); emax = 0;
for j = 1:NT
    dth = thseq(j+1)-thseq(j);
    for i = 1:N-1
        ei = e(i,:)';
        ei = Rot(-dth)*ei + (Rot(-dth)-eye(2))*d(i,:)';
        ei = segrelax_rk4(ei, dt, vmax);
        e(i,:) = ei';
    end
    emax = max(emax, mean(vecnorm(e,2,2)));
end
unified = emax;
resid = 100*(unified - full)/full;

% 端点（同口径：全程 max）
full_sq   = ardpf_run(N, ds, 0, 1, 'square', 4, 1.0).efrom_max(1);
full_circ = ardpf_run(N, ds, 0, 1, 'circle', 0, 1.0).efrom_max(1);

S = struct('full', full, 'unified', unified, 'resid', resid, ...
           'full_sq', full_sq, 'full_circ', full_circ);
end


function e1 = segrelax_rk4(e0, dt, vmax)
h = min(0.005, dt); t = 0; e = e0;
while t < dt
    h = min(h, dt-t);
    k1 = f(e,vmax); k2 = f(e+0.5*h*k1,vmax); k3 = f(e+0.5*h*k2,vmax); k4 = f(e+h*k3,vmax);
    e = e + h*(k1+2*k2+2*k3+k4)/6; t = t+h;
    if norm(e) < 1e-6, break; end
end
e1 = e;
end


function dd = f(e, vmax)
rho = norm(e); if rho < 1e-9, dd = [-1;0]; else, dd = -vmax*e/rho - [1;0]; end
end


function Rm = Rot(a)
Rm = [cos(a), -sin(a); sin(a), cos(a)];
end
