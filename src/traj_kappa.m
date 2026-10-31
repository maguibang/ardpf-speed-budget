function [x,y,s,kappa,L] = traj_kappa(kind, Rc)
%% traj_kappa.m —— 轨迹与曲率的【唯一】实现（供全部脚本共用）
%
%  此前 acceptance_test.m 有 bt/bt9 两份、decompose_independent.m 有第三份 bt，
%  导致口径修正只落到一份、兄弟拷贝漏改（T6/T8 的同类 bug 即源于此）。
%  现统一到本文件。
%
%  kind: 'circle' | 'spiral' | 'fig8' | 'randpoly'
%  返回：x,y 轨迹点；s 弧长；kappa 曲率 dθ/ds；L 周长

switch kind
    case 'circle'
        th = linspace(0, 2*pi, 12000)';
        x = Rc*cos(th);  y = Rc*sin(th);
    case 'spiral'
        th = linspace(0, 2*pi, 12000);
        r  = Rc*(1 + 0.45*cos(3*th));
        x  = r.*cos(th);  y = r.*sin(th);
    case 'fig8'
        u = linspace(0, 2*pi, 12000);
        x = Rc*sin(u);  y = Rc*sin(u).*cos(u);
    case 'randpoly'
        rng(1,'twister');
        nv = 7;  ang = sort(2*pi*rand(1,nv));  rad = Rc*(0.55 + 0.45*rand(1,nv));
        vx = rad.*cos(ang);  vy = rad.*sin(ang);
        x = [];  y = [];
        for k = 1:nv
            k2 = mod(k,nv)+1;
            x = [x, linspace(vx(k), vx(k2), 400)]; %#ok<AGROW>
            y = [y, linspace(vy(k), vy(k2), 400)]; %#ok<AGROW>
        end
        x(end+1) = x(1);  y(end+1) = y(1);
    otherwise
        error('traj_kappa: unknown kind %s', kind);
end

x = x(:)';  y = y(:)';
keep = [true, hypot(diff(x),diff(y)) > 1e-9];   % 去重复点
x = x(keep);  y = y(keep);

s = [0, cumsum(hypot(diff(x),diff(y)))];
dx = gradient(x, s);  dy = gradient(y, s);
kappa = gradient(unwrap(atan2(dy, dx)), s);
L = s(end);
end
