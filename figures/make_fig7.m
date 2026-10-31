function make_fig7()
%% make_fig7.m —— 图 7：速度余量的非单调性（U 形）+ 闭式 A(delta) 下降段
%  散点 = Table 13（末圈稳态，N=10, d_sep=5, square）全仿真。
%  绿线 = 闭式 A(delta)（附录 D，A_headroom_closedform），只画饱和分支
%         v_max/v <= 2.0，即闭式有效区间。
%  标记 A=1 交叉点（v_max/v = 1.3977）。

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src')));

OUT = fullfile(fileparts(mfilename('fullpath')), '..', 'figs');
if ~exist(OUT,'dir'), mkdir(OUT); end

K.blue   = [  0 114 178]/255;
K.orange = [230 159   0]/255;
K.red    = [213  94   0]/255;
K.grey   = [ 90  90  90]/255;
K.green  = [  0 158 115]/255;
K.purple = [ 86  62 148]/255;

% Table 13 数据（末圈稳态 e_form）
vm  = [1.0 1.1 1.2 1.5 2.0 5.0 10.0];
ef  = [24.72 23.48 22.38 20.14 19.59 19.31 20.17];
stab = [0 1 1 1 1 1 0];   % 0=失稳, 1=稳定

f = figure('Position',[100 100 760 400]);
hold on;

% 阈值线
yline(25, '--', 'Color', K.grey, 'LineWidth', 1.1);
text(1.25, 25.25, 'E_{max} = 25 m', 'Color', K.grey, 'FontSize', 8.5);

% 闭式 A(delta)：饱和分支 v_max/v in [1.0, 2.0]
vcf = 1.0:0.05:2.0;
ecf = zeros(size(vcf));
for k = 1:numel(vcf)
    [~, ecf(k)] = A_headroom_closedform(4, 10, 5, vcf(k), 2);
end
pcf = plot(vcf, ecf, '-', 'Color', K.green, 'LineWidth', 2.0, ...
    'DisplayName', 'closed form  $A(\delta)$');

% A=1 交叉点
g = @(vm_) A_at_one(vm_);
vx = fzero(g, [1.05, 3.0]);
[~, ex] = A_headroom_closedform(4, 10, 5, vx, 2);
pcross = plot(vx, ex, 'd', 'MarkerSize', 9, 'MarkerFaceColor', K.purple, ...
    'MarkerEdgeColor', 'w', 'LineWidth', 0.7, 'DisplayName', 'A = 1');

% 稳定段（蓝）+ 失稳点（红）
hsim = gobjects(1,0);
for k = 1:numel(vm)
    if stab(k), c = K.blue; else, c = K.red; end
    hsim(end+1) = plot(vm(k), ef(k), 'o', 'MarkerSize', 7, ...
        'MarkerFaceColor', c, 'MarkerEdgeColor', 'w', 'LineWidth', 0.7);
end
% 仿真连线
plot(vm, ef, '-', 'Color', K.grey, 'LineWidth', 1.1);

% 谷底标注
[~, imin] = min(ef);
plot(vm(imin), ef(imin), 's', 'MarkerSize', 11, 'MarkerFaceColor', K.orange, ...
    'MarkerEdgeColor', 'w', 'LineWidth', 0.7);
text(vm(imin)+0.28, ef(imin), sprintf('optimum %s %.1f m', char(8776), ef(imin)), ...
    'Color', K.orange, 'FontSize', 8.5, 'FontWeight', 'bold');

% 失稳点标注
text(1.02, ef(1)-1.0, 'error > E_{max}', 'Color', K.red, 'FontSize', 8.5);
text(vm(end)-0.1, ef(end)-1.1, 'collision', 'Color', K.red, 'FontSize', 8.5, ...
    'HorizontalAlignment', 'right');

% 交叉点标注
text(vx+0.06, ex+0.35, sprintf('A=1 at %.2f', vx), 'Color', K.purple, ...
    'FontSize', 8.5, 'FontWeight', 'bold');

% 图例
legend([hsim(1), pcf, pcross], {'full simulation', 'closed form  A(\delta)', 'A = 1'}, ...
    'Location', 'northeast', 'FontSize', 8, 'Box', 'off');

set(gca, 'XScale', 'log');
xlim([0.9 12]); ylim([18 27]);
xticks([1 1.5 2 5 10]); xticklabels({'1.0','1.5','2.0','5.0','10.0'});
xlabel('speed limit ratio  v_{max}/v');
ylabel('last-lap  e_{form}  (m)');
title('Speed headroom is non-monotone: closed form on the falling branch', 'FontSize', 9);
box on;

exportfig(f, fullfile(OUT,'Fig7_headroom'), K);
fprintf('  已写 figs/Fig7_headroom.png / .eps\n');
end


function a1 = A_at_one(vm_)
% A(vmax) - 1，用于找 A=1 交叉点
[A, ~] = A_headroom_closedform(4, 10, 5, vm_, 2);
a1 = A - 1;
end


function exportfig(f, name, K) %#ok<INUSD>
set(f, 'ToolBar', 'none', 'MenuBar', 'none');
ax = findall(f, 'Type', 'axes');
for k = 1:numel(ax)
    try, ax(k).Toolbar = []; catch, end %#ok<CTCH>
end
print(f, [name '.png'], '-dpng', '-r400');
print(f, [name '.eps'], '-depsc', '-painters');
end
