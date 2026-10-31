%% gate_vs_angle.m —— 几何门随转角 Δθ 的变化（60° vs 90° vs 45°）
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'code_release')));
N=10; c10=2+7/9; b0=1/sqrt(2); Emax=25;

% 论文 Table 6 的 A（仿真值）
Asim = [1.2226 1.4239 1.6278 2.0390 2.8767];
Mv   = [4 6 8 12 20];

fprintf('\n%-4s %-8s %-8s %-10s %-12s %-14s\n','M','dtheta','kappa','A(sim)','A*kappa','|d|* 预测');
fprintf('%s\n', repmat('-',1,62));
pred = zeros(size(Mv));
for k=1:numel(Mv)
    kap = 2*sin(pi/Mv(k));
    pred(k) = (Emax-b0)/(Asim(k)*kap);
    fprintf('%-4d %-8.0f %-8.4f %-10.4f %-12.4f %-14.2f\n', ...
        Mv(k), 360/Mv(k), kap, Asim(k), Asim(k)*kap, pred(k));
end

% 实测验证：M=6（60°），扫 d_sep 找首个 e_form>25
fprintf('\n=== 实测验证 M=6 (60°) ===\n');
fprintf('%-8s %-10s %-12s %-10s\n','d_sep','mean|d|','e_form_max','超25?');
for d = [5.4 5.6 5.8 6.0 6.2 6.4]
    R = ardpf_run(N, d, 0, 1, 'poly', 6);
    fprintf('%-8.1f %-10.2f %-12.3f %-10s\n', d, c10*d, R.efrom_max(1), ...
        tern_str(R.efrom_max(1) > 25));
end
fprintf('  预测门 = %.2f m\n', pred(Mv==6));

function s = tern_str(c)
if c, s='是'; else, s='否'; end
end
