%% delay_sweep.m —— 时延区间：几何门随 τ 的退化
%  扫 τ ∈ {0..15} s，对每个 τ 扫 d_sep 找零延迟/非零延迟失稳的临界 mean|d|。
%  验证论文论断：延迟只是极限环的次级调制，几何门在 τ>0 时温和退化而非崩溃。

clear; clc;
addpath(genpath(fileparts(mfilename('fullpath'))));

N = 10;  traj = 'square';  nseed = 8;
tau_list = [0 2 5 10 15];
D = [3.5 4.0 4.5 5.0 5.5];
c10 = 2.7777778;

fprintf('\n====================================================================\n');
fprintf('  时延区间：临界 mean|d| 随 τ 的退化  (N=10, square, 8 seeds)\n');
fprintf('  c(10)=%.4f,  稳定阈值 >=80%%\n', c10);
fprintf('====================================================================\n');

crit = zeros(size(tau_list));
for ti = 1:numel(tau_list)
    tau = tau_list(ti);
    fprintf('\n--- tau = %g s ---\n', tau);
    fprintf('%-8s %-10s %-12s %-10s\n','d_sep','mean|d|','eform_max','稳定率');
    dstar = NaN;
    for d = D
        R = ardpf_run(N, d, tau, nseed, traj);
        fprintf('%-8.1f %-10.2f %-12.3f %-8.0f%%\n', d, c10*d, R.efrom_max(1), 100*R.stabRate);
        if R.stabRate < 0.8 && isnan(dstar)
            dstar = d;   % 首个失稳（<80%）
        end
    end
    if ~isnan(dstar)
        crit(ti) = c10*dstar;
        fprintf('  => 临界 mean|d|* = %.2f m\n', crit(ti));
    else
        crit(ti) = NaN;
        fprintf('  => 全扫描稳定（门 < %.2f m 或未触及）\n', c10*D(1));
    end
end

fprintf('\n==== 汇总：门随 τ 的退化 ====\n');
fprintf('%-10s %-16s %-12s\n','tau (s)','临界 mean|d|* (m)','相对 τ=0');
for ti = 1:numel(tau_list)
    if isnan(crit(ti))
        fprintf('%-10g %-16s %s\n', tau_list(ti), '<9.7', '—');
    else
        fprintf('%-10g %-16.2f %+.0f%%\n', tau_list(ti), crit(ti), 100*(crit(ti)-crit(1))/crit(1));
    end
end
fprintf('\n');
