function S = decompose_independent()
%% decompose_independent.m —— 三项分解：每一项都由【独立模型】给出，残差 = 模型误差
%
%  独立来源（均不调用全仿真、无拟合常数）：
%    poly        : C 由 A_closedform 的定点迭代（标量超越方程）
%    randpoly    : C 由复合映射定点 generalized_A
%    circle      : R 由恒定 ω 闭式 g = sqrt(S²-v_max²)/ω
%    spiral/fig8 : R 由降阶单跟随者 ODE ė=-v_max ê-(v,0)-ω(t)J(e+d)
%
%  记账：拐角主导 pred = b0 + C；旋转主导 pred = R（b0 被 R 内含）。
%  返回 S：.name .full .b0 .C .R .pred .resid .source

Rc = 15*sqrt(2); v = 1.0; vmax = 1.0; N = 10; ds = 5; dt = 0.1;
md = ds*ceil((N-1)/2)^2/(2*ceil((N-1)/2)-1);  b0 = 1/sqrt(2);  % N=10 (even): closed form for mean|d_i|
phi = pi/4;  d = zeros(N-1,2);
for i = 1:N-1
    r = ceil(i/2);  s = 2*mod(i,2)-1;
    d(i,:) = [-r*ds*cos(phi), s*r*ds*sin(phi)];
end

R_ = {};   % 每行: {name, full, b0, C, R, source}

%% ---- poly: C 由 A_closedform 独立给出，R=0 ----
for M = [4 6 8 12 20]
    C = A_closedform(M, N, ds) * 2*sin(pi/M) * md;
    Rs = sim_instr(N, ds, 'poly', M, 1.0);
    Tk = 2*Rc*sin(pi/M)/v;  Tlap = M*Tk;
    full = max(Rs.efrom(Rs.t >= 4*Tlap));
    R_(end+1,:) = {sprintf('poly M=%d',M), full, b0, C, 0, 'A_closedform'}; %#ok<SAGROW>
end

%% ---- randpoly: C 由复合映射定点独立给出，R=0 ----
g = generalized_A();
R_(end+1,:) = {'randpoly', g.full, b0, g.C, 0, 'composite fixed point'};

%% ---- circle: R 由恒定 ω 闭式独立给出 ----
Rcl = circle_closed(N, d, v, vmax, Rc);
Rf  = ardpf_run(N, ds, 0, 1, 'circle', 0, 1.0);
tf  = (0:numel(Rf.efrom_t)-1)'*dt;
fullC = max(Rf.efrom_t(tf >= 0.8*tf(end)));
R_(end+1,:) = {'circle', fullC, 0, 0, Rcl, 'closed form g (const omega)'};

%% ---- spiral / fig8: R 由降阶 ODE 独立给出 ----
for kind = {'spiral','fig8'}
    Rm = reduced_rotation(kind{1}, N, d, v, vmax, Rc);
    Rf = ardpf_run(N, ds, 0, 1, kind{1}, 0, 1.0);
    tf = (0:numel(Rf.efrom_t)-1)'*dt;
    full = max(Rf.efrom_t(tf >= 0.8*tf(end)));
    R_(end+1,:) = {kind{1}, full, 0, 0, Rm, 'reduced ODE (independent)'}; %#ok<SAGROW>
end

%% ---- 组装 struct 并打印 ----
n = size(R_,1);
S = struct('name',{R_(:,1)}, 'full',[], 'b0',[], 'C',[], 'R',[], ...
           'pred',[], 'resid',[], 'source',{R_(:,6)});
S.full = cell2mat(R_(:,2));  S.b0 = cell2mat(R_(:,3));
S.C    = cell2mat(R_(:,4));  S.R  = cell2mat(R_(:,5));
S.pred = S.b0 + S.C + S.R;
S.resid = 100*(S.pred - S.full)./S.full;

fprintf('\n================================================================================================\n');
fprintf('  三项分解（每项由独立模型给出；残差 = (预测-全仿真)/全仿真）\n');
fprintf('================================================================================================\n');
fprintf('%-12s %10s %8s %10s %10s %10s %9s  %s\n', ...
        'traj','full_sim','b0','C','R','pred','resid_%','independent_source');
fprintf('%s\n', repmat('-',1,104));
for k = 1:n
    fprintf('%-12s %10.3f %8.3f %10.3f %10.3f %10.3f %9.2f  %s\n', ...
        S.name{k}, S.full(k), S.b0(k), S.C(k), S.R(k), S.pred(k), S.resid(k), S.source{k});
end
fprintf('%s\n', repmat('-',1,104));
fprintf('注：拐角主导 pred = b0 + C；旋转主导 pred = R（b0 被 R 内含，内侧跟随者近零误差）。\n');
fprintf('    口径统一为【最后一圈 max】：模型与全仿真均同（修正 spiral 此前用全程 max 的不一致）。\n');
fprintf('    9 条轨迹各自由其主项给出、无余项，残差 <= %.2f%%。\n', max(abs(S.resid)));
fprintf('    ※ 9 行中没有任何一行三项同时非零：加式只在"单主项"退化区被检验。\n');
fprintf('================================================================================================\n\n');

%% ---- 纯 ASCII 输出 ----
fid = fopen('decompose_independent.txt','w');
cl = onCleanup(@() fclose(fid));
p = @(varargin) fprintf(fid, varargin{:});
p('THREE-TERM DECOMPOSITION - every term from an INDEPENDENT model\n');
p('residual = (pred - full_sim)/full_sim ;  pred = b0 + C + R ; steady state (last lap)\n\n');
p('%-12s %10s %8s %10s %10s %10s %9s  %s\n', ...
  'traj','full_sim','b0','C','R','pred','resid_%','independent_source');
for k = 1:n
    p('%-12s %10.3f %8.3f %10.3f %10.3f %10.3f %9.2f  %s\n', ...
      S.name{k}, S.full(k), S.b0(k), S.C(k), S.R(k), S.pred(k), S.resid(k), S.source{k});
end
p('\nNOTE: corner-dominated -> pred = b0 + C ; rotation-dominated -> pred = R\n');
p('      (b0 absorbed inside R: inner followers have ~0 error).\n');
p('      Every row uses an independent model; max |residual| = %.2f %%.\n', max(abs(S.resid)));
p('\nSCALAR GROWTH-RATE IDENTITY (body frame, saturated):\n');
p('  d|e|/dt = -v_max - v*cos(theta_e) - omega*(e_hat . J d)\n');
p('  corners inject via impulse; rotation drives via -omega*(e_hat . J d).\n');
end


%% ================= 旋转项：降阶 ODE（独立） =================
function emax = reduced_rotation(kind, N, d, v, vmax, Rc)
% 口径：与 full 一致 —— 取【最后一圈】max（用公共函数 reduced_ode，
%       它同时返回全程/末圈两种口径，这里取末圈）。
[~, emax] = reduced_ode(kind, N, d, v, vmax, Rc, 0.005, 0);
end


%% ================= 圆：恒定 ω 闭式 R =================
function R = circle_closed(N, d, v, vmax, Rc)
% 恒定 ω 闭式，对全部 N-1 个跟随者取均值（内侧 g_i = 0）
om = v/Rc;  R = 0;
for i = 1:N-1
    a = -d(i,1);  b = d(i,2);
    if b < 0
        S = sqrt((v-om*b)^2 + (om*a)^2);
        R = R + sqrt(S^2 - vmax^2)/om;
    end
end
R = R/(N-1);
end
