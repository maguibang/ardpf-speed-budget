function make_figures()
%% make_figures.m —— 生成论文 9 张图中的 8 张（彩色，Okabe-Ito 色盲友好配色）
%
%  第 1 段：重算需要原始数据的图（图3 转向反转、图9 临界间距外推），存 fig_data.mat
%  第 2 段：由论文数据（表 2/4/10/11/13）与 fig_data 出这 8 张图
%  输出：figs/Fig1..Fig6, Fig8, Fig9 .png(400dpi) + .eps(矢量)
%        （图 7 由 make_fig7.m 生成）
%
%  配色：Okabe & Ito (2008) 色盲友好八色板。转灰度后仍可分辨。
%  图号 = 正文出现顺序：1 机理(§5.1) 2 降阶模型(§5.2) 3 转向反转(§5.4)
%                        4 零延迟(§5.5) 5 分段圆弧(§7.3) 6 三项分解(§8)
%                        7 速度余量(§9, make_fig7.m) 8 幂律对比(§10.1) 9 外推(§11)

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src')));

OUT = fullfile(fileparts(mfilename('fullpath')), '..', 'figs');
if ~exist(OUT,'dir'), mkdir(OUT); end

% ---------- 统一配色 ----------
K.blue   = [  0 114 178]/255;   % 主色：机理 / 全仿真
K.orange = [230 159   0]/255;   % 次色：模型 / 预测
K.green  = [  0 158 115]/255;   % 第三：解析闭式
K.red    = [213  94   0]/255;   % 警示：失稳 / 误差
K.purple = [204 121 167]/255;
K.sky    = [ 86 180 233]/255;
K.yellow = [240 228  66]/255;
K.grey   = [ 90  90  90]/255;

LW = 1.4;  MS = 5;  FS = 9;
set(0,'defaultAxesFontName','Arial','defaultAxesFontSize',FS, ...
      'defaultLineLineWidth',LW,'defaultAxesLineWidth',0.8, ...
      'defaultAxesTickDir','out','defaultFigureColor','w');

%% ==================== 第 1 段：重算原始数据 ====================
if exist('fig_data.mat','file')
    load('fig_data.mat');
    fprintf('已载入 fig_data.mat（删除该文件可强制重算）\n');
else
    N = 10; ds = 5; dt = 0.1; v = 1.0; Rc = 15*sqrt(2);
    fprintf('重算原始数据……\n');

    % ---- 图3：转向反转（M=8，末圈逐跟随者峰值）----
    M = 8;
    A = sim_instr(N, ds, 'poly',   M, 1.0);
    B = sim_instr(N, ds, 'polycw', M, 1.0);
    Tk = 2*Rc*sin(pi/M)/v;  Tlap = M*Tk;
    last = A.t >= 4*Tlap;
    pe1 = zeros(N-1,1);  pe2 = zeros(N-1,1);
    for i = 1:N-1
        pe1(i) = max(vecnorm(squeeze(A.E(last,i,:)),2,2));
        pe2(i) = max(vecnorm(squeeze(B.E(last,i,:)),2,2));
    end
    b_off = A.d0(2:end,2);
    turn = struct('b',b_off,'ccw',pe1,'cw',pe2, ...
                  'r_ccw',corr(pe1,b_off),'r_cw',corr(pe2,b_off));
    fprintf('  转向反转：corr = %+.3f (CCW) / %+.3f (CW)\n', turn.r_ccw, turn.r_cw);

    % ---- 图5：临界 d_sep* 随 N 外推（每格只跑一次）----
    cfgs = {10,[4 5]; 12,[3.5 4.5]; 15,[3 3.5]; 20,[2 2.5 2.6]};
    nrow = 0;
    for ci = 1:size(cfgs,1), nrow = nrow + numel(cfgs{ci,2}); end
    Nv = zeros(nrow,1); Dv = zeros(nrow,1); Sv = zeros(nrow,1); k = 0;
    for ci = 1:size(cfgs,1)
        n = cfgs{ci,1};  mdb = mean_d(n);
        for d = cfgs{ci,2}
            R = ardpf_run(n, d, 0, 1, 'square');
            k = k + 1;
            Nv(k) = n;  Dv(k) = d;  Sv(k) = R.stabRate;
            fprintf('  N=%-3d d_sep=%-4.1f  mean|d|=%5.1f  k*mean|d|=%5.1f  stab=%.0f%%\n', ...
                    n, d, mdb*d, sqrt(2)*mdb*d, 100*R.stabRate);
        end
    end
    % 每个 N 的临界点 = 最小的使 stab<0.8 的 d_sep
    Nu = unique(Nv);  Dstar = nan(size(Nu));
    for i = 1:numel(Nu)
        m = (Nv == Nu(i)) & (Sv < 0.8);
        if any(m), Dstar(i) = min(Dv(m)); end
    end
    ok = ~isnan(Dstar);
    prod_ = Nu(ok).*Dstar(ok);
    scal = struct('N',Nu,'dstar',Dstar,'prod',prod_, ...
                  'mean_prod',mean(prod_),'cv',100*std(prod_)/mean(prod_));
    fprintf('  d_sep*·N = %.1f, CV = %.2f%%\n', scal.mean_prod, scal.cv);

    save('fig_data.mat','turn','Nv','Dv','Sv','scal');
end

%% ==================== 图 1：机理示意 ====================
f1 = figure('Position',[100 100 900 340]);

% (a) 机体系与零预算方向
ax1 = subplot(1,2,1); hold(ax1,'on');
fill([-3 3 3 -3],[-1.6 -1.6 1.6 1.6], K.sky,'FaceAlpha',0.16,'EdgeColor','none');
text(0,1.38,'un-restored along-track direction','HorizontalAlignment','center', ...
     'Color',K.grey,'FontSize',8.5,'FontAngle','italic');
quiver(0,0,2.1,0,0,'Color',K.blue,'LineWidth',2.2,'MaxHeadSize',0.5);
text(1.05,0.22,'leader  v','Color',K.blue,'FontWeight','bold','FontSize',9);
plot(0,0,'o','MarkerSize',7,'MarkerFaceColor',K.blue,'MarkerEdgeColor','none');
plot(-1.6,-0.95,'s','MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','none');
text(-1.55,-1.28,'slot  d_i','FontSize',8.5);
quiver(-1.6,-0.95,1.05,0.62,0,'Color',K.red,'LineWidth',1.8,'MaxHeadSize',0.6);
text(-0.45,-0.15,'e','Color',K.red,'FontWeight','bold','FontSize',10);
quiver(-2.55,-1.62,0,0.62,0,'Color',K.green,'LineWidth',2.0,'MaxHeadSize',0.8);
quiver(2.45,-1.62,0,0.62,0,'Color',K.orange,'LineWidth',2.0,'MaxHeadSize',0.8);
text(-2.5,-2.02,'corner impulse','Color',K.green,'FontSize',8.5,'HorizontalAlignment','center');
text(2.45,-2.02,'curvature forcing','Color',K.orange,'FontSize',8.5,'HorizontalAlignment','center');
plot([-3 3],[0 0],':','Color',K.grey);  plot([0 0],[-2.2 2.2],':','Color',K.grey);
xlim([-3.2 3.2]); ylim([-2.3 2.2]); axis(ax1,'equal'); axis(ax1,'off');
title(ax1,'(a)  Body frame:  v_{max}=v  exhausts the along-track budget','FontSize',FS);

% (b) 增长率恒等式
ax2 = subplot(1,2,2); hold(ax2,'on');
th = linspace(0,2*pi,400);
rate = -1 - cos(th);
plot(th*180/pi, rate, 'Color',K.blue, 'LineWidth',2.0);
yline(0,'--','Color',K.grey);
plot(180,0,'o','MarkerSize',7,'MarkerFaceColor',K.red,'MarkerEdgeColor','none');
text(180,0.20,'decay rate = 0','Color',K.red,'FontSize',8.5, ...
     'HorizontalAlignment','center');
xlim([0 360]); ylim([-2.35 0.45]);
xticks(0:90:360); yticks(-2:0.5:0);
xlabel('error direction \theta_e  (deg)');
ylabel('d|e|/dt   (v_{max}=v=1)');
title(ax2,'(b)  Marginal along-track direction','FontSize',FS);
box(ax2,'on');
exportfig(f1, fullfile(OUT,'Fig1_mechanism'), K);

%% ==================== 图 2：零延迟失稳证据（表 4）====================
cfgname = {'5 m, N=10','8 m, N=7','8 m, N=10','3 m, N=5','3 m, N=7','5 m, N=5'};
stab = [  0  0  0  0;
          0  0  0  0;
          0  0  0  0;
        100 100 100  0;
        100 100 100  0;
        100 100 100 80];
tau_lab = {'\tau=0','\tau=0.1','\tau=2','\tau=15'};

f2 = figure('Position',[100 100 820 360]);
hb = bar(stab, 'grouped'); hold on;
cols = {K.red,K.orange,K.sky,K.blue};
for i = 1:4, hb(i).FaceColor = cols{i}; hb(i).EdgeColor = 'none'; end
set(gca,'XTickLabel',cfgname,'YLim',[0 108]);
ylabel('stability rate  (%)');
plot([0.5 3.5],[80 80],'--','Color',K.grey,'LineWidth',1.1);
text(2.0,84,'80% stability threshold','Color',K.grey,'FontSize',8, ...
     'HorizontalAlignment','center');
text(2.0,40,sprintf('unstable at ZERO delay'),'Color',K.red,'FontWeight','bold', ...
     'FontSize',10,'HorizontalAlignment','center');
title('Zero-delay instability: three configurations are unstable with no delay at all', ...
      'FontSize',FS);
legend(hb, tau_lab, 'Location','northoutside','Orientation','horizontal','Box','off');
box on;
exportfig(f2, fullfile(OUT,'Fig4_zerodelay'), K);

%% ==================== 图 3：转向反转 ====================
f3 = figure('Position',[100 100 820 350]);
for s = 1:2
    ax = subplot(1,2,s); hold(ax,'on');
    if s == 1, y = turn.ccw; r = turn.r_ccw; nm = 'counter-clockwise'; c = K.blue;
    else,      y = turn.cw;  r = turn.r_cw;  nm = 'clockwise';         c = K.orange;
    end
    yy = linspace(min(turn.b)*1.1, max(turn.b)*1.1, 50);
    pf = polyfit(turn.b, y, 1);
    plot(yy, polyval(pf,yy), '-', 'Color', c, 'LineWidth',1.6);
    scatter(turn.b, y, 46, 'MarkerFaceColor',c,'MarkerEdgeColor','w', ...
            'LineWidth',0.7,'MarkerFaceAlpha',0.9);
    xlabel('cross-track offset  b_i  (m)');
    if s == 1, ylabel('peak formation error  max_t |e_i|  (m)'); end
    title(sprintf('(%c)  %s      corr = %+.3f', 'a'+s-1, nm, r), 'FontSize',FS);
    xlim([min(yy) max(yy)]); ylim([0 47]); box(ax,'on');
end
sgtitle('Turn reversal mirrors the error distribution across followers', ...
        'FontSize',FS+1);
exportfig(f3, fullfile(OUT,'Fig3_turnreversal'), K);

%% ==================== 图 4：降阶模型 vs 全仿真（表 2）====================
tname = {'poly M=4','poly M=6','poly M=8','poly M=12','poly M=20', ...
         'spiral','figure-eight','random polygon'};
full_ = [25.089 20.982 18.668 16.182 14.106 15.385 19.721 28.785];
red_  = [24.847 20.916 18.657 16.141 14.071 15.150 19.552 28.554];

f4 = figure('Position',[100 100 780 400]);
hb = barh([full_; red_]', 'grouped'); hold on;
hb(1).FaceColor = K.blue;   hb(1).EdgeColor = 'none';
hb(2).FaceColor = K.orange; hb(2).EdgeColor = 'none';
set(gca,'YTickLabel',tname,'YDir','reverse');
xlabel('max_t e_{form}  (m)');
legend({'full N-body simulation','reduced single-follower model'}, ...
       'Location','northoutside','Orientation','horizontal','Box','off');
for i = 1:8
    dev = 100*(red_(i)-full_(i))/full_(i);
    text(max(full_(i),red_(i))+0.55, i, sprintf('%+.1f%%', dev), ...
         'FontSize',7.6,'Color',K.grey,'VerticalAlignment','middle');
end
xlim([0 33]);
title('Reduced model reproduces the full simulation within 1.53%', 'FontSize',FS);
box on;
exportfig(f4, fullfile(OUT,'Fig2_reducedmodel'), K);

%% ==================== 图 5：临界间距外推（E7）====================
f5 = figure('Position',[100 100 700 400]);
hold on;
nn = linspace(8, 21, 200);
plot(nn, scal.mean_prod./nn, '-', 'Color', K.grey, 'LineWidth',1.4);
ok = ~isnan(scal.dstar);
scatter(scal.N(ok), scal.dstar(ok), 70, 'MarkerFaceColor',K.blue, ...
        'MarkerEdgeColor','w','LineWidth',0.8);
scatter(scal.N(~ok), 2*ones(sum(~ok),1), 70, 'MarkerFaceColor','w', ...
        'MarkerEdgeColor',K.grey,'LineWidth',1.4);
text(20, 2.25, 'no violation found', 'Color',K.grey,'FontSize',8, ...
     'HorizontalAlignment','right');
text(13.2, 4.6, sprintf('d_{sep}^{\\ast} \\times N = %.1f   (CV = %.2f%%)', ...
     scal.mean_prod, scal.cv), 'Color',K.grey,'FontSize',9);
xlabel('swarm size  N'); ylabel('critical spacing  d_{sep}^{\ast}  (m)');
xlim([8 21]); ylim([1.6 5.6]);
title('Critical spacing scales as 1/N  —  the gate extrapolates beyond N=10', ...
      'FontSize',FS);
box on;
exportfig(f5, fullfile(OUT,'Fig9_extrapolation'), K);

%% ==================== 图 6：三项分解（表 11）====================
dn = {'poly M=4','poly M=6','poly M=8','poly M=12','poly M=20', ...
      'random poly','circle','spiral','figure-eight'};
b0 = [0.707 0.707 0.707 0.707 0.707 0.707 0 0 0];
Cc = [23.961 19.855 17.391 14.751 12.591 26.287 0 0 0];
Rr = [0 0 0 0 0 0 10.092 15.098 19.602];
fm = [24.722 20.483 18.011 15.367 13.208 27.989 10.227 15.243 19.719];

f6 = figure('Position',[100 100 820 380]);
hb = bar([b0; Cc; Rr]', 'stacked'); hold on;
hb(1).FaceColor = K.sky;    hb(1).EdgeColor = 'none';
hb(2).FaceColor = K.orange; hb(2).EdgeColor = 'none';
hb(3).FaceColor = K.green;  hb(3).EdgeColor = 'none';
plot(0.55:9.45, fm, 'd', 'Color',K.grey,'MarkerFaceColor','w', ...
     'MarkerEdgeColor',K.grey,'LineWidth',1.2,'MarkerSize',6);
set(gca,'XTickLabel',dn,'YLim',[0 33]);
ylabel('steady formation error  (m)');
legend({'b_0  baseline','C  corner','R  curvature','full simulation'}, ...
       'Location','northoutside','Orientation','horizontal','Box','off');
title('Three-term decomposition: every term from an independent model', 'FontSize',FS);
box on;
exportfig(f6, fullfile(OUT,'Fig6_decomposition'), K);

%% ==================== 图 7：幂律 vs 判据（表 13 + 断崖）====================
f7 = figure('Position',[100 100 850 350]);

% (a) 幂律预测 vs 实测边界（9 个格子）
ax = subplot(1,2,1); hold(ax,'on');
c_ = 17.65; al = 0.204; be = 0.249;
gN  = [5  5  5  7  7  7 10 10 10];
gd  = [3  5  8  3  5  8  3  5  8];
pred = c_./(gN.^al .* gd.^be);
% 实测 tau_max = 通过 80% 阈值的最大的被测时延
meas = [10 15  5  5 10  0 10  0  0];
xx = 1:9;
hb = bar(xx, [pred(:) meas(:)], 'grouped'); hold(ax,'on');
hb(1).FaceColor = K.red;  hb(1).EdgeColor = 'none';
hb(2).FaceColor = K.blue; hb(2).EdgeColor = 'none';
plot(xx(meas==0), zeros(1,sum(meas==0)), 'v', 'MarkerSize',7, ...
     'MarkerFaceColor',K.green,'MarkerEdgeColor','none');
text(5.5, 13.6, '\nabla  measured \tau_{max} = 0   (unstable with no delay)', ...
     'Color',K.green,'FontSize',8,'HorizontalAlignment','center');
set(ax,'XTick',xx,'XTickLabel', ...
    {'N5\newline3m','N5\newline5m','N5\newline8m','N7\newline3m','N7\newline5m', ...
     'N7\newline8m','N10\newline3m','N10\newline5m','N10\newline8m'}, 'FontSize',7.5);
ylim([0 18]);
ylabel('\tau_{max}  (s)');
legend({'power law (1)','measured'}, 'Location','northwest','Box','off','FontSize',8);
title(ax,'(a)  The law predicts 6.6–9.7 s everywhere','FontSize',FS); box(ax,'on');

% (b) 两项指标对比
ax = subplot(1,2,2);
% 显式给 x 坐标，避免 bar() 对矩阵分组方向的歧义
pl  = [36 67];    % 幂律：平均相对误差%, 二分准确率%
cr  = [ 0 100];   % 判据
dx  = 0.16;
hold(ax,'on');
bar((1:2)-dx, pl, 2*dx, 'FaceColor',K.red,  'EdgeColor','none');
bar((1:2)+dx, cr, 2*dx, 'FaceColor',K.green,'EdgeColor','none');
set(ax,'XTick',1:2, 'XTickLabel',{'mean relative error (%)','binary accuracy (%)'}, ...
    'YLim',[0 118], 'XLim',[0.45 2.55]);
labs = {'36%','67%'; '0%','100%'};
xs   = [1-dx, 2-dx; 1+dx, 2+dx];
ys   = [pl; cr];
for i = 1:2
    for j = 1:2
        if ys(i,j) < 6
            text(xs(i,j), ys(i,j)+3, labs{i,j}, 'HorizontalAlignment','center', ...
                 'FontSize',8,'Color',K.grey);
        else
            text(xs(i,j), ys(i,j)+3.5, labs{i,j}, 'HorizontalAlignment','center', ...
                 'FontSize',8,'FontWeight','bold');
        end
    end
end
legend({'fitted power law','mechanism-based criterion'},'Location','northoutside', ...
       'Orientation','horizontal','Box','off');
title(ax,'(b)  Criterion vs. fitted law','FontSize',FS); box(ax,'on');
exportfig(f7, fullfile(OUT,'Fig8_comparison'), K);

%% ==================== 图 8：分段圆弧收敛（表 10）====================
nseg = [4 8 16 32 64];
esp  = [-9.7  2.6  3.2 -1.3 -0.9];
ef8  = [-16.5 -3.5 -2.3 -2.5 -0.8];

f8 = figure('Position',[100 100 700 380]);
hold on;
hsp = semilogx(nseg, esp, '-o', 'Color',K.blue,   'MarkerFaceColor',K.blue,   'MarkerSize',MS);
hfg = semilogx(nseg, ef8, '-s', 'Color',K.orange, 'MarkerFaceColor',K.orange, 'MarkerSize',MS);
plot([4 64],[0 0], '--', 'Color',K.grey);
% 8–16 段的工作区用竖向色带标出（不要压在曲线上）
hband = fill([8 16 16 8], [-20 -20 20 20], K.red, ...
           'FaceAlpha',0.07, 'EdgeColor','none');
uistack(hband, 'bottom');
text(11.3, 17.6, '8–16 segments: |error| \leq 3.5%', 'Color',K.red, ...
     'FontSize',8.5,'HorizontalAlignment','center');
yline(15,'--','Color',K.grey); yline(-15,'--','Color',K.grey);
xlabel('number of constant-curvature segments');
ylabel('error vs. full simulation  (%)');
legend([hsp hfg],{'spiral','figure-eight'},'Location','southeast','Box','off');
xlim([3.6 72]); ylim([-20 20]); xticks(nseg); xticklabels(nseg);
title('Segmented-arc approximation: convergence is non-monotone', 'FontSize',FS);
box on;
exportfig(f8, fullfile(OUT,'Fig5_segmentedarc'), K);

fprintf('\n全部 8 张图已输出到 %s/\n', OUT);
end


%% ==================== 辅助函数 ====================
function exportfig(f, name, K) %#ok<INUSD>
% 同时导出 400 dpi PNG（投稿）与矢量 EPS（排版）
% 先剥掉坐标区工具栏，否则它会出现在导出的图里
set(f, 'ToolBar', 'none', 'MenuBar', 'none');
ax = findall(f, 'Type', 'axes');
for k = 1:numel(ax)
    try, ax(k).Toolbar = []; catch, end %#ok<CTCH>
end
print(f, '-dpng', '-r400', [name '.png']);
print(f, '-depsc', '-painters', [name '.eps']);
fprintf('  已写 %s.png / .eps\n', name);
end

function m = mean_d(N)
% 45° 半角 V 形编队的平均偏移（单位 d_sep），r_i = ceil(i/2) = 1,1,2,2,...
%   N even (N-1 = 2R-1):  sum r_i = R^2     -> m = R^2/(2R-1),  R = N/2
%   N odd  (N-1 = 2R  ):  sum r_i = R(R+1)  -> m = (R+1)/2  = (N+1)/4
if mod(N,2) == 0
    R = N/2;
    m = R^2 / (2*R - 1);
else
    m = (N+1)/4;
end
end
