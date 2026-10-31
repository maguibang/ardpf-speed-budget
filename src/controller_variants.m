%% controller_variants.m —— 多控制器验证：判据对势场变体的普适性
%  在标称构型(N=10)上扫 d_sep，找零延迟失稳的临界 mean|d|，对比四个势场变体：
%    1. eq16      ARD-PF 原控制律（基线）
%    2. eq12grad  排斥项取严格梯度（2*K_rep/r^4）
%    3. Katt=4    吸引增益加倍
%    4. Kff=1     速度前馈（跟随者速度被拉向领航者速度）
%    5. Kff=2     更强的前馈
%  预期：排斥形式与温和前馈不改变门（~13.5 m）；吸引增益按 b0=v_max/sqrt(Katt) 平移。

clear; clc;
addpath(genpath(fileparts(mfilename('fullpath'))));

N = 10;  tau = 0;  nseed = 1;  traj = 'square';
D = [3.0 3.5 4.0 4.5 5.0 5.5 6.0];
c10 = 2.7777778;   % mean|d| = c10 * d_sep

variants = {
  'ARD-PF eq16 (baseline)', @(d) ardpf_run(N,d,0,1,'square',4,1.0,'eq16',45,[],2.0,0);
  'repform eq12grad',       @(d) ardpf_run(N,d,0,1,'square',4,1.0,'eq12grad',45,[],2.0,0);
  'Katt = 4',               @(d) ardpf_run(N,d,0,1,'square',4,1.0,'eq16',45,[],4.0,0);
  'feedforward Kff = 1',    @(d) ardpf_run(N,d,0,1,'square',4,1.0,'eq16',45,[],2.0,1.0);
  'feedforward Kff = 2',    @(d) ardpf_run(N,d,0,1,'square',4,1.0,'eq16',45,[],2.0,2.0);
};

fprintf('\n====================================================================\n');
fprintf('  多控制器验证：临界 mean|d| 对比 (N=10, square, tau=0)\n');
fprintf('  c(10)=%.4f,  mean|d| = c(10)*d_sep\n', c10);
fprintf('====================================================================\n');

for vi = 1:size(variants,1)
    name = variants{vi,1};  fn = variants{vi,2};
    fprintf('\n--- %s ---\n', name);
    fprintf('%-8s %-10s %-12s %-10s\n','d_sep','mean|d|','eform_max','稳定率');
    dstar = NaN;  last_stable = NaN;
    for d = D
        R = fn(d);
        md = c10 * d;
        fprintf('%-8.1f %-10.2f %-12.3f %-8.0f%%\n', d, md, R.efrom_max(1), 100*R.stabRate);
        if R.stabRate < 0.8 && isnan(dstar)
            dstar = d;   % 首个失稳
        end
    end
    if ~isnan(dstar)
        fprintf('  => 临界 mean|d|* = c10 * %.1f = %.2f m\n', dstar, c10*dstar);
    else
        fprintf('  => 全扫描范围均稳定（门 > %.2f m）\n', c10*D(end));
    end
end
fprintf('\n');
