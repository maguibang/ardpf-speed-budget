%% circle_headroom.m —— 圆轨迹：v_max 超过 max S_i 后误差是否归零？
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'code_release')));
N=10; ds=5; v=1; Rc=15*sqrt(2); phi=pi/4;
r = ceil((1:N-1)/2);
a = r*ds*cos(phi);  b = r*ds*sin(phi);
om = v/Rc;
S_out = sqrt((v + om*b).^2 + (om*a).^2);   % 外侧跟随者 (b_i<0)
fprintf('\n  Rc=%.3f  omega=%.4f rad/s\n', Rc, om);
fprintf('  |d_i| : '); fprintf('%.1f ', r*ds); fprintf('\n');
fprintf('  S_out : '); fprintf('%.3f ', S_out); fprintf('\n');
fprintf('  max S_i = %.4f  -> 理论所需 v_max 上界\n\n', max(S_out));

fprintf('%-10s %-14s %-14s %-12s\n','v_max','e_form_max','末圈','是否 S<v_max');
fprintf('%s\n', repmat('-',1,56));
for vs = [1.0 1.5 2.0 2.1 2.5 3.0 5.0]
    R = ardpf_run(N, ds, 0, 1, 'circle', 0, vs);
    dt=0.1; t=(0:numel(R.efrom_t)-1)*dt;
    last = max(R.efrom_t(t >= 0.8*t(end)));
    fprintf('%-10.2f %-14.4f %-14.4f %-12s\n', vs, R.efrom_max(1), last, ...
        tern_str(vs >= max(S_out)));
end

function s = tern_str(c)
if c, s='是（无饱和）'; else, s='否（仍有饱和）'; end
end
