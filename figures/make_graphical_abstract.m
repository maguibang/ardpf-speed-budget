function make_graphical_abstract()
%% make_graphical_abstract.m —— MDPI Graphical Abstract（横向三格：机理 → 模型 → 判据）
%  输出 graphical_abstract.png 2200x1120（MDPI 下限 1100x560）/ _hi.png 4400x2240 / .eps

K.blue   = [  0 114 178]/255;
K.orange = [230 159   0]/255;
K.green  = [  0 158 115]/255;
K.red    = [213  94   0]/255;
K.sky    = [ 86 180 233]/255;
K.grey   = [ 90  90  90]/255;

set(0,'defaultAxesFontName','Arial','defaultAxesFontSize',11, ...
      'defaultLineLineWidth',1.4,'defaultAxesTickDir','out','defaultFigureColor','w');

f = figure('Position',[80 80 1100 560],'Color','w','ToolBar','none','MenuBar','none');
set(f,'PaperPositionMode','auto');

% ---------------- 页眉 ----------------
annotation(f,'textbox',[0 0.925 1 0.075],'String', ...
  {'Along-Track Speed Budget — Not Delay — Sets Formation Error in Saturated UAV Swarms'}, ...
  'HorizontalAlignment','center','VerticalAlignment','middle','EdgeColor','none', ...
  'FontSize',15,'FontWeight','bold');
annotation(f,'textbox',[0 0.878 1 0.045],'String', ...
  {'ARD-PF leader–follower swarms   |   one reduced model, two closed forms, one design criterion'}, ...
  'HorizontalAlignment','center','VerticalAlignment','middle','EdgeColor','none', ...
  'FontSize',10,'Color',K.grey);

% ================= (a) 机理 =================
ax = axes('Position',[0.020 0.115 0.305 0.660]); hold(ax,'on'); axis(ax,'off');
fill([-2.3 2.4 2.4 -2.3],[-1.45 -1.45 1.45 1.45], K.sky,'FaceAlpha',0.16,'EdgeColor','none');
text(-0.05,1.22,'no restoring authority','HorizontalAlignment','center', ...
     'Color',K.grey,'FontSize',9.5,'FontAngle','italic');
% 跟随者（V 形，位于领航者之后）
fv = [-0.45 -0.58; -0.45 0.58; -1.30 -1.06; -1.30 1.06];
for i = 1:size(fv,1)
    plot(fv(i,1),fv(i,2),'^','MarkerSize',8.5,'MarkerFaceColor','k','MarkerEdgeColor','none');
end
% 领航者
quiver(0.30,0,1.55,0,0,'Color',K.blue,'LineWidth',2.6,'MaxHeadSize',0.55);
plot(0.30,0,'o','MarkerSize',8,'MarkerFaceColor',K.blue,'MarkerEdgeColor','none');
text(1.05,0.26,'leader  v','Color',K.blue,'FontWeight','bold','FontSize',10);
text(0.30,-0.32,'v_{max}=v','Color',K.blue,'FontSize',9.5,'HorizontalAlignment','center');
% 两条注入路径
quiver(-0.95,-1.98,0,0.50,0,'Color',K.green, 'LineWidth',2.4,'MaxHeadSize',0.95);
quiver( 0.75,-1.98,0,0.50,0,'Color',K.orange,'LineWidth',2.4,'MaxHeadSize',0.95);
text(-0.95,-2.32,'corner','Color',K.green, 'FontSize',9.5,'HorizontalAlignment','center');
text( 0.75,-2.32,'curvature','Color',K.orange,'FontSize',9.5,'HorizontalAlignment','center');
plot([-2.3 2.4],[0 0],':','Color',K.grey);
xlim([-2.5 2.6]); ylim([-2.55 1.65]);
title('(a)  saturation exhausts the along-track budget','FontSize',10.5,'FontWeight','bold');

% ================= (b) 统一模型 =================
ax = axes('Position',[0.370 0.200 0.225 0.560]); hold(ax,'on');
full_ = [25.089 20.982 18.668 16.182 14.106 15.385 19.721 28.785];
red_  = [24.847 20.916 18.657 16.141 14.071 15.150 19.552 28.554];
lim = [8 31];
plot(lim,lim,'--','Color',K.grey,'LineWidth',1.1);
scatter(full_, red_, 42,'MarkerFaceColor',K.orange,'MarkerEdgeColor','w','LineWidth',0.7);
xlim(lim); ylim(lim); axis(ax,'square');
xlabel('full N-body sim  (m)','FontSize',9.5);
ylabel('reduced model  (m)','FontSize',9.5);
text(9.3,30.2,sprintf('8 trajectories\nmax deviation 1.5%%'),'FontSize',9.5, ...
     'Color',K.red,'FontWeight','bold','VerticalAlignment','top');
title('(b)  one follower suffices','FontSize',10.5,'FontWeight','bold');

% ================= (c) 三项分解 =================
ax = axes('Position',[0.685 0.200 0.295 0.560]); hold(ax,'on');
b0 = [0.707 0.707 0     0];
Cc = [23.961 26.287 0   0];
Rr = [0     0     10.092 15.098];
fm = [24.722 27.989 10.227 15.243];
hb = bar(1:4, [b0; Cc; Rr]','stacked');
hb(1).FaceColor = K.sky;    hb(1).EdgeColor = 'none';
hb(2).FaceColor = K.orange; hb(2).EdgeColor = 'none';
hb(3).FaceColor = K.green;  hb(3).EdgeColor = 'none';
plot(1:4, fm,'d','Color',K.grey,'MarkerFaceColor','w','MarkerEdgeColor',K.grey, ...
     'LineWidth',1.1,'MarkerSize',5.5);
% 标签放在柱顶以上的空白区，避免横跨柱间空隙
text(1.35, 30.6,'corner  C','Color',K.orange,'FontSize',9.5,'FontWeight','bold', ...
     'HorizontalAlignment','center');
text(3.45, 30.6,'curvature  R','Color',K.green,'FontSize',9.5,'FontWeight','bold', ...
     'HorizontalAlignment','center');
text(2.4, 33.4,'max residual 3.6%','Color',K.red,'FontSize',9.5,'FontWeight','bold', ...
     'HorizontalAlignment','center');
set(ax,'XTick',1:4,'XTickLabel',{'M=4','random','circle','spiral'}, ...
    'FontSize',9.5,'YLim',[0 35.5],'XLim',[0.45 4.55]);
ylabel('steady e_{form}  (m)','FontSize',9.5);
title('(c)  closed-form error budget','FontSize',10.5,'FontWeight','bold');

% ---------------- 页脚 ----------------
annotation(f,'textbox',[0 0.018 1 0.060],'String', ...
  {'Critical mean formation offset 13.5 m — independent of swarm size   |   three of nine fitted configurations are already unstable at ZERO delay'}, ...
  'HorizontalAlignment','center','VerticalAlignment','middle','EdgeColor','none', ...
  'FontSize',10.5,'Color',K.red,'FontWeight','bold');

% ---------------- 导出 ----------------
for a = findall(f,'Type','axes')', try, a.Toolbar = []; catch, end, end
print(f,'-dpng','-r200','graphical_abstract.png');
print(f,'-dpng','-r400','graphical_abstract_hi.png');
print(f,'-depsc','-painters','graphical_abstract.eps');
fprintf('已写 graphical_abstract.png (2200x1120) / _hi.png / .eps\n');
end
