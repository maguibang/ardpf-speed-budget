function [emax_all, emax_last] = reduced_ode(kind, N, d, v, vmax, Rc, h, nseg)
%% reduced_ode.m —— 降阶单跟随者 ODE 的【唯一】积分实现
%
%   ė = -v_max ê - (v,0) - ω(t) J (e + d_i)
%
%  同时返回两种口径，供调用方按需取用（避免"模型用全程、参照用末圈"这类不一致）：
%     emax_all  : 5 圈【全程】max —— 含启动瞬态
%     emax_last : 最后【一圈】max  —— 稳态，与全仿真的末圈口径可比
%
%  nseg = 0 : 用真实曲率 ω(t) = v·κ(t)          （T3/T8/T9 用）
%  nseg > 0 : κ 用 nseg 段【分段常数】近似        （T6 分段圆弧用）
%
%  注：nseg=0 仅适用于曲率连续的轨迹（circle/spiral/fig8）；含拐角的轨迹
%      用复合映射定点（generalized_A.m），不走本函数。

[x,y,s,kappa,L] = traj_kappa(kind, Rc); %#ok<ASGLU>
Ttot = 5*L/v;
tt = 0:h:Ttot;
ss = mod(v*tt, L);

if nseg > 0
    edges = linspace(0, L, nseg+1);
    om_pc = zeros(1, nseg);
    for k = 1:nseg
        m = s >= edges(k) & s <= edges(k+1);
        om_pc(k) = v * mean(kappa(m));
    end
    segidx = min(floor(ss/L*nseg)+1, nseg);
    om = om_pc(segidx);
else
    om = interp1(s, kappa, ss, 'linear', 'extrap') * v;
end

J = [0 -1; 1 0];
e = zeros(N-1, 2);
emax_all = 0;  emax_last = 0;
lastT = 0.8 * tt(end);          % 末圈起点（5 圈 -> 后 20%）

for j = 1:numel(tt)-1
    wj = om(j);
    for i = 1:N-1
        ei = e(i,:)';  di = d(i,:)';
        f = @(x) vmax*(-x/max(norm(x),1e-12)) - [v;0] - wj*J*(x+di);
        k1 = f(ei);  k2 = f(ei+0.5*h*k1);  k3 = f(ei+0.5*h*k2);  k4 = f(ei+h*k3);
        e(i,:) = (ei + h*(k1+2*k2+2*k3+k4)/6)';
    end
    em = mean(vecnorm(e,2,2));
    emax_all = max(emax_all, em);
    if tt(j) >= lastT, emax_last = max(emax_last, em); end
end
end
