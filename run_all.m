function run_all()
%% run_all.m —— one-shot reproduction: acceptance suite + all figures
%
%  用法 / Usage:
%      cd <this folder>
%      run_all
%
%  输出 / Outputs (written to the current folder):
%      acceptance_result.txt          acceptance suite log
%      decompose_independent.txt      three-term decomposition table
%      figs/Fig1..Fig9 .png / .eps    manuscript figures
%      graphical_abstract*.png/.eps   graphical abstract
%
%  预计耗时 15–25 分钟（大部分在轨迹仿真）。只需验收套件时改跑 acceptance_test。

here = fileparts(mfilename('fullpath'));
addpath(genpath(here));
cd(here);

fprintf('\n########## 1/3   acceptance suite ##########\n');
acceptance_test();

fprintf('\n########## 2/3   manuscript figures ##########\n');
make_figures();
make_fig7();

fprintf('\n########## 3/3   graphical abstract ##########\n');
make_graphical_abstract();

fprintf('\n=========================================================\n');
fprintf(' Done. See acceptance_result.txt, decompose_independent.txt,\n');
fprintf(' figs/ and graphical_abstract*.png in %s\n', here);
fprintf('=========================================================\n');
end
