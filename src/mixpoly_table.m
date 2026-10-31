function T = mixpoly_table(S, rho, n_per_m)
%% mixpoly_table.m —— 混合轨迹：2 个尖角 + 2 个圆角（闭合）
%
%  补充"拐角+旋转"混合型轨迹，用于检验三项加法分解
%  （此前 8 条轨迹各自被单一主项主导，从未测过混合型）
%
%  路径（逆时针，起点 (rho,0)，航向 0）：
%    (rho,0) ──线──> (S,0)  【尖角 +90°】
%            ──线──> (S,S-rho) ──弧(中心(S-rho,S-rho), 90°)──> (S-rho,S)
%            ──线──> (0,S)  【尖角 +90°】
%            ──线──> (0,rho) ──弧(中心(rho,rho), 90°)──> (rho,0)   闭合
%
%  特殊情形：
%    rho = 0    -> 退化为普通正方形（4 个尖角）      = 纯拐角
%    rho = S/2  -> 圆角占满，尖角边消失              = 接近纯旋转
%    中间 rho   -> 尖角与圆弧并存                    = 混合
%
%  输入：S 边长，rho 圆角半径，n_per_m 每米采样点数
%  输出：T .s(弧长) .x .y .th(航向) .L

if nargin < 3 || isempty(n_per_m), n_per_m = 20; end
ds = 1/n_per_m;

X = rho;  Y = 0;  TH = 0;        % 起点
cur = [rho, 0];

L1 = S - rho;                    % 直线段长度（四段相同）
if L1 < 0, error('rho 必须 <= S/2'); end
Larc = rho * (pi/2);

segs = { 'line', L1, 0;
         'jump', pi/2,  [];
         'line', L1, pi/2;
         'arc',  pi/2, [S-rho, S-rho];
         'line', L1, pi;
         'jump', 3*pi/2, [];
         'line', L1, 3*pi/2;
         'arc',  pi/2, [rho, rho] };

for k = 1:size(segs,1)
    switch segs{k,1}
        case 'line'
            L = segs{k,2};  hh = segs{k,3};
            if L <= 0, continue; end
            n = max(2, round(L/ds));
            t = linspace(0, L, n+1);  t(1) = [];
            X = [X, cur(1) + t*cos(hh)];   %#ok<AGROW>
            Y = [Y, cur(2) + t*sin(hh)];   %#ok<AGROW>
            TH = [TH, hh*ones(1,n)];       %#ok<AGROW>
            cur = cur + [L*cos(hh), L*sin(hh)];

        case 'jump'
            hh = segs{k,2};
            X = [X, cur(1)];  Y = [Y, cur(2)];  TH = [TH, hh];  %#ok<AGROW>

        case 'arc'
            ang = segs{k,2};  c = segs{k,3};
            if Larc < ds/2        % rho~0：退化为尖角
                hh = atan2(cur(2)-c(2), cur(1)-c(1)) + pi/2 + ang;
                X = [X, cur(1)];  Y = [Y, cur(2)];  TH = [TH, hh];  %#ok<AGROW>
                continue;
            end
            a0 = atan2(cur(2)-c(2), cur(1)-c(1));
            n  = max(2, round(Larc/ds));
            aa = linspace(a0, a0+ang, n+1);  aa(1) = [];
            X = [X, c(1)+rho*cos(aa)];        %#ok<AGROW>
            Y = [Y, c(2)+rho*sin(aa)];        %#ok<AGROW>
            TH = [TH, aa + pi/2];             %#ok<AGROW>  切向航向 = 极角 + 90°
            cur = [c(1)+rho*cos(a0+ang), c(2)+rho*sin(a0+ang)];
    end
end

X = X(:)';  Y = Y(:)';  TH = TH(:)';
% 去重（避免弧长网格非唯一）
keep = [true, hypot(diff(X), diff(Y)) > 1e-9];
X = X(keep);  Y = Y(keep);  TH = TH(keep);

s = [0, cumsum(hypot(diff(X), diff(Y)))];
T.s = s;  T.x = X;  T.y = Y;  T.th = unwrap(TH);  T.L = s(end);
end
