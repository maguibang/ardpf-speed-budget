%% dynamics_sweep.m —— 一阶速度环（旋翼动力学）对门的影响
%  在标称构型(N=10, square)上扫 tau_v（速度环时间常数），看临界 mean|d| 是否不变。
%  机理：一阶惯性只延迟速度响应，不改变"v_max=v 时沿航向追不上"的稳态性质。

clear; clc;
addpath(genpath(fileparts(mfilename('fullpath'))));

N = 10;  tau = 0;  nseed = 1;  traj = 'square';
D = [3.0 3.5 4.0 4.5 5.0 5.5 6.0];
c10 = 2.7777778;
tau_vs = [0 0.1 0.3 0.5 1.0];   % 速度环时间常数（s），0=纯运动学

fprintf('\n====================================================================\n');
fprintf('  旋翼动力学验证：临界 mean|d| 随速度环时间常数 tau_v 的变化\n');
fprintf('  (N=10, square, tau=0, c(10)=%.4f)\n', c10);
fprintf('====================================================================\n');

crit = zeros(size(tau_vs));
for ti = 1:numel(tau_vs)
    tv = tau_vs(ti);
    if tv == 0
        name = 'pure kinematics';
    else
        name = sprintf('1st-order speed loop, tau_v=%.1f s', tv);
    end
    fprintf('\n--- %s ---\n', name);
    fprintf('%-8s %-10s %-12s %-10s\n','d_sep','mean|d|','eform_max','稳定率');
    dstar = NaN;
    for d = D
        R = ardpf_run(N,d,0,1,'square',4,1.0,'eq16',45,[],2.0,0,tv);
        fprintf('%-8.1f %-10.2f %-12.3f %-8.0f%%\n', d, c10*d, R.efrom_max(1), 100*R.stabRate);
        if R.stabRate < 0.8 && isnan(dstar)
            dstar = d;
        end
    end
    if ~isnan(dstar)
        crit(ti) = c10*dstar;
        fprintf('  => 临界 mean|d|* = %.2f m\n', crit(ti));
    else
        crit(ti) = NaN;
        fprintf('  => 全扫描稳定（门 > %.2f m）\n', c10*D(end));
    end
end

fprintf('\n==== 汇总 ====\n');
fprintf('%-28s %-14s %-10s\n','tau_v (s)','临界 mean|d|* (m)','相对0');
for ti = 1:numel(tau_vs)
    if isnan(crit(ti))
        fprintf('%-28s %-14s %s\n', sprintf('%.1f', tau_vs(ti)), '>16.7', '—');
    else
        fprintf('%-28s %-14.2f %+.0f%%\n', sprintf('%.1f', tau_vs(ti)), crit(ti), ...
                100*(crit(ti)-crit(1))/crit(1));
    end
end
fprintf('\n');
