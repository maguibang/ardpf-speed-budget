%% mixpoly_check.m —— 混合轨迹（尖角+圆角）上验证降阶模型与三项分解
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'code_release')));
N=10; ds=5; v=1; vmax=1; dt=0.1; laps=5;
phi=pi/4; d=zeros(N-1,2);
for i=1:N-1, r=ceil(i/2); s=2*mod(i,2)-1; d(i,:)=[-r*ds*cos(phi), s*r*ds*sin(phi)]; end
S=30;

rhos = [6 12 15];   % 圆角半径：6 偏拐角，12 混合，15 接近纯旋转
fprintf('\n=== 全仿真 mixpoly（2尖角+2圆角, S=30）===\n');
fprintf('%-6s %-14s %-14s %-10s\n','rho','full max e_form','last-lap','稳定性');
fprintf('%s\n', repmat('-',1,48));
fulls = zeros(size(rhos));
for k=1:numel(rhos)
    R = ardpf_run(N, ds, 0, 1, 'mixpoly', 0, 1.0, 'eq16', 45, rhos(k));
    t=(0:numel(R.efrom_t)-1)*dt;
    fulls(k)=R.efrom_max(1); last=max(R.efrom_t(t>=0.8*t(end)));
    fprintf('%-6g %-14.3f %-14.3f %-10.0f%%\n', rhos(k), fulls(k), last, 100*R.stabRate);
end

fprintf('\n=== 统一离散映射（降阶模型）vs 全仿真 ===\n');
fprintf('%-6s %-14s %-14s %-9s\n','rho','full maxE','unified','ratio');
fprintf('%s\n', repmat('-',1,48));
for k=1:numel(rhos)
    rho = rhos(k);
    T = mixpoly_table(S, rho, 40);
    Lper = T.L;
    t = 0:dt:laps*Lper/v;
    s = mod(v*t, Lper);
    thseq = interp1(T.s, T.th, s, 'linear', 'extrap');
    % 尖角处：航向跳变，interp1 会抹圆；用差分检测跳变并强制一步跳
    NT = numel(thseq)-1;
    e = zeros(N-1,2); emax=0;
    for j=1:NT
        dth = thseq(j+1)-thseq(j);
        for i=1:N-1
            ei = e(i,:)';
            ei = Rot(-dth)*ei + (Rot(-dth)-eye(2))*d(i,:)';
            ei = segrelax_rk4(ei, dt, vmax);
            e(i,:) = ei';
        end
        emax = max(emax, mean(vecnorm(e,2,2)));
    end
    fprintf('%-6g %-14.3f %-14.3f %-9.3f\n', rho, fulls(k), emax, emax/fulls(k));
end

function e1 = segrelax_rk4(e0, dt, vmax)
h=min(0.005,dt); t=0; e=e0;
while t<dt
    h=min(h,dt-t);
    k1=f(e,vmax); k2=f(e+0.5*h*k1,vmax); k3=f(e+0.5*h*k2,vmax); k4=f(e+h*k3,vmax);
    e=e+h*(k1+2*k2+2*k3+k4)/6; t=t+h;
    if norm(e)<1e-6, break; end
end
e1=e;
end
function dd=f(e,vmax)
rho=norm(e); if rho<1e-9, dd=[-1;0]; else, dd=-vmax*e/rho-[1;0]; end
end
function Rm=Rot(a), Rm=[cos(a),-sin(a);sin(a),cos(a)]; end
