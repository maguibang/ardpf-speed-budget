function S = unified_discrete()
%% unified_discrete.m —— 统一离散映射：一条模型覆盖全部 8 类轨迹
%
%  每步 dt：领航者航向变化 Δθ；机体系内
%     e ← Rot(-Δθ) e + (Rot(-Δθ)-I) d      （拐角=大 Δθ；平滑旋转=小 Δθ）
%     然后段内饱和弛豫 dt（RK4）
%  拐角(任务1)与平滑旋转(任务2)是同一映射的两个极限。
%
%  返回 S：.name(8x1 cell) .full(8x1) .model(8x1) .ratio(8x1)

Rc = 15*sqrt(2); v = 1.0; vmax = 1.0; N = 10; ds = 5; dt = 0.1;
phi = pi/4;  d = zeros(N-1,2);
for i = 1:N-1
    r = ceil(i/2);  s = 2*mod(i,2)-1;
    d(i,:) = [-r*ds*cos(phi), s*r*ds*sin(phi)];
end

specs = {'poly',4; 'poly',6; 'poly',8; 'poly',12; 'poly',20; 'spiral',0; 'fig8',0; 'randpoly',0};
docE  = [25.089 20.982 18.668 16.182 14.106 15.385 19.721 28.785];

fprintf('\n=====================================================================================\n');
fprintf('  统一离散映射 vs 全仿真   (8 类轨迹, N=%d, d_sep=%g, v_max=1)\n', N, ds);
fprintf('=====================================================================================\n');
fprintf('%-12s %12s %12s %9s\n','traj','full maxE','unified','ratio');
fprintf('%s\n', repmat('-',1,50));

S = struct('name',{{}}, 'full',[], 'model',[], 'ratio',[]);
for k = 1:size(specs,1)
    kind = specs{k,1};  M = specs{k,2};
    [thseq, Lper] = heading_seq(kind, M, Rc, v, dt, 5);
    NT = numel(thseq) - 1;

    e = zeros(N-1,2);   eform_max = 0;
    for j = 1:NT
        dth = thseq(j+1) - thseq(j);
        for i = 1:N-1
            ei = e(i,:)';
            ei = Rot(-dth)*ei + (Rot(-dth)-eye(2))*d(i,:)';
            ei = seg_relax(ei, dt, vmax);
            e(i,:) = ei';
        end
        eform_max = max(eform_max, mean(vecnorm(e,2,2)));
    end

    if M > 0, nm = sprintf('poly M=%d', M); else, nm = kind; end
    fprintf('%-12s %12.3f %12.3f %9.3f\n', nm, docE(k), eform_max, eform_max/docE(k));
    S.name{end+1,1} = nm;  S.full(end+1,1) = docE(k);
    S.model(end+1,1) = eform_max;  S.ratio(end+1,1) = eform_max/docE(k);
end
fprintf('=====================================================================================\n\n');
end


function w = seg_relax(e0, dt, vmax)
w = rk4(e0, dt, vmax);
end

function e1 = rk4(e0, dt, vmax)
h = min(0.005, dt);  t = 0;  e = e0;
while t < dt
    h = min(h, dt - t);
    k1 = f(e,vmax);  k2 = f(e+0.5*h*k1,vmax);  k3 = f(e+0.5*h*k2,vmax);  k4 = f(e+h*k3,vmax);
    e = e + h*(k1+2*k2+2*k3+k4)/6;  t = t + h;
    if norm(e) < 1e-6, break; end
end
e1 = e;
end

function dd = f(e, vmax)
rho = norm(e);
if rho < 1e-9, dd = [-1;0]; else, dd = -vmax*e/rho - [1;0]; end
end


function [thseq, Lper] = heading_seq(kind, M, Rc, v, dt, laps)
switch kind
    case 'poly'
        Lper = M*2*Rc*sin(pi/M);
        t = 0:dt:laps*Lper/v;
        thseq = zeros(size(t));
        for j = 1:numel(t)
            s = mod(v*t(j), Lper);
            side = 2*Rc*sin(pi/M);
            kk = floor(s/side);
            a1 = 2*pi*kk/M;  a2 = 2*pi*(kk+1)/M;
            v1 = Rc*[cos(a1), sin(a1)];  v2 = Rc*[cos(a2), sin(a2)];
            dd = v2 - v1;  thseq(j) = atan2(dd(2), dd(1));
        end
    case {'spiral','fig8','randpoly'}
        T = traj_table(kind, Rc, 1);
        Lper = T.L;
        t = 0:dt:laps*Lper/v;
        s = mod(v*t, Lper);
        thseq = interp1(T.s, T.th, s, 'linear', 'extrap');
end
end


function T = traj_table(kind, Rc, seed)
switch kind
    case 'spiral'
        th = linspace(0, 2*pi, 6000);
        r  = Rc*(1 + 0.45*cos(3*th));
        x  = r.*cos(th);  y = r.*sin(th);
    case 'fig8'
        u  = linspace(0, 2*pi, 6000);
        x  = Rc*sin(u);   y = Rc*sin(u).*cos(u);
    case 'randpoly'
        rng(seed,'twister');
        nv = 7; ang = sort(2*pi*rand(1,nv)); rad = Rc*(0.55 + 0.45*rand(1,nv));
        vx = rad.*cos(ang); vy = rad.*sin(ang);
        x = []; y = [];
        for k = 1:nv
            k2 = mod(k,nv)+1;
            x = [x, linspace(vx(k), vx(k2), 400)]; %#ok<AGROW>
            y = [y, linspace(vy(k), vy(k2), 400)]; %#ok<AGROW>
        end
        x(end+1) = x(1);  y(end+1) = y(1);
end
x = x(:)';  y = y(:)';
keep = [true, hypot(diff(x),diff(y)) > 1e-9];
x = x(keep);  y = y(keep);
ds_ = hypot(diff(x),diff(y));  s = [0, cumsum(ds_)];
dx = gradient(x, s);  dy = gradient(y, s);
T.s = s;  T.x = x;  T.y = y;  T.th = unwrap(atan2(dy, dx));  T.L = s(end);
end


function Rm = Rot(a)
Rm = [cos(a), -sin(a); sin(a), cos(a)];
end
