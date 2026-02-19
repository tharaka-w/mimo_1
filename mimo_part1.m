
clear; clc; close all;


% 1) Nominal parameters (from proposal)
C1  = 3.0e6;      % [J/°C]
C2  = 2.5e6;      % [J/°C]
R1o = 0.010;      % [°C/W]
R2o = 0.012;      % [°C/W]
R12 = 0.020;      % [°C/W]
Ebat = 13.5 * 3.6e6;   % 13.5 kWh -> Joules [J]

% 2) Build state-space model for u -> y
% States:  x = [T1; T2; s]
% Inputs:  u = [Q1; Q2; Pbat]
% Outputs: y = [T1; T2; s]
A = [ -(1/(C1*R1o) + 1/(C1*R12)),    1/(C1*R12),               0;
       1/(C2*R12),                 -(1/(C2*R2o) + 1/(C2*R12)), 0;
       0,                           0,                         0 ];

B2 = [ 1/C1,   0,      0;
       0,     1/C2,    0;
       0,      0,   -1/Ebat ];

C2 = eye(3);
D22 = zeros(3,3);

Gss = ss(A, B2, C2, D22);
Gss.InputName  = {'Q1','Q2','Pbat'};
Gss.OutputName = {'T1','T2','s'};


% 3) Print state-space matrices
disp('--- State-space (u -> y) ---');
disp('A =');   disp(A);
disp('B2 =');  disp(B2);
disp('C2 =');  disp(C2);
disp('D22 ='); disp(D22);


% 4) Print FULL transfer matrix G(s) entries
Gtf = tf(Gss);

disp(' ');
disp('--- FULL Transfer matrix G_yu(s) entries ---');
disp('Rows: y = [T1; T2; s], Cols: u = [Q1; Q2; Pbat]');
for i = 1:size(Gtf,1)
    for j = 1:size(Gtf,2)
        fprintf('\nG(%d,%d)(s):  %s <- %s\n', i, j, Gss.OutputName{i}, Gss.InputName{j});
        disp(Gtf(i,j));
    end
end
disp(' ');


% 5) Poles and transmission zeros (built-in)
polesG = pole(Gss);
zerosG = tzero(Gss);

disp('--- Poles of G_yu(s) ---');
disp(polesG);

disp('--- Transmission zeros from tzero(G_yu) ---');
if isempty(zerosG)
    disp('No finite transmission zeros (empty).');
else
    disp(zerosG);
end

if all(real(polesG) < 0)
    disp('Stability check: All poles have negative real part (stable).');
else
    disp('Stability check: Not all poles have negative real part (NOT asymptotically stable).');
end

% 6) Rosenbrock rank-drop verification (complex grid + random samples)
% Rosenbrock matrix:
%   R(s) = [ sI - A   -B2
%            C2        D22 ]
% Transmission zeros = finite s where rank(R(s)) drops below normal rank.
disp(' ');
disp('--- Rosenbrock rank-drop verification ---');

n = size(A,1);
R = @(s) [s*eye(n) - A, -B2;
          C2,        D22];

% Estimate normal rank by random complex sampling near the dynamics scale
rng(1);
Nsamp = 80;
rank_samp = zeros(Nsamp,1);
for k = 1:Nsamp
    re = (2*rand-1) * 5e-3;    % [-5e-3, 5e-3]
    im = (2*rand-1) * 5e-3;
    s  = re + 1j*im;
    rank_samp(k) = rank(R(s), 1e-10);
end
normal_rank = max(rank_samp);
fprintf('Estimated normal rank of R(s): %d (max over random samples)\n', normal_rank);

% Verify at each tzero, if any
if isempty(zerosG)
    disp('tzero returned empty => no finite transmission zeros to verify.');
else
    disp('Rank check at each transmission zero candidate:');
    for k = 1:numel(zerosG)
        s0 = zerosG(k);
        rk = rank(R(s0), 1e-10);
        fprintf('  s0 = %+0.6e %+0.6ej : rank(R(s0)) = %d (normal rank = %d)\n', ...
            real(s0), imag(s0), rk, normal_rank);
    end
end

% Grid search for any rank drop (stronger than real-axis only)
re_grid = linspace(-5e-3, 5e-3, 61);
im_grid = linspace(-5e-3, 5e-3, 61);
min_rank = normal_rank;
min_s    = 0;

for a = 1:numel(re_grid)
    for b = 1:numel(im_grid)
        s = re_grid(a) + 1j*im_grid(b);
        rk = rank(R(s), 1e-10);
        if rk < min_rank
            min_rank = rk;
            min_s = s;
        end
    end
end

fprintf('Minimum rank found on complex grid: %d at s = %+0.3e %+0.3ej\n', ...
    min_rank, real(min_s), imag(min_s));

if min_rank < normal_rank
    disp('=> Rank drop detected somewhere on the grid: possible transmission zero exists.');
else
    disp('=> No rank drop detected on the grid: consistent with no finite transmission zeros.');
end



%%Rosenbrock rank test 
% For (A,B,C,D), transmission zeros are s where rank(R(s)) drops below normal rank.
% Rosenbrock matrix:
%   R(s) = [ sI - A   -B
%            C       D  ]

A = A; B = B2; C = C2; D = D22;   % (u -> y) plant
n = size(A,1);

R = @(s) [s*eye(n) - A, -B;
          C,        D];

% Estimate normal rank using random complex samples
rng(1);
Ns = 80;
ranks = zeros(Ns,1);
for k = 1:Ns
    s = ((2*rand-1) + 1j*(2*rand-1)) * 5e-3;   % sample around 0 (thermal poles ~1e-4)
    ranks(k) = rank(R(s), 1e-10);
end
normal_rank = max(ranks);

fprintf('\n===== Rosenbrock Transmission-Zero Proof =====\n');
fprintf('R(s) size = %d x %d\n', size(R(0),1), size(R(0),2));
fprintf('Estimated normal rank of R(s): %d\n', normal_rank);

% Grid scan in complex plane to try to find any rank drop
re_grid = linspace(-5e-3, 5e-3, 61);
im_grid = linspace(-5e-3, 5e-3, 61);

min_rank = normal_rank;
min_s    = 0;

for a = 1:numel(re_grid)
    for b = 1:numel(im_grid)
        s = re_grid(a) + 1j*im_grid(b);
        rk = rank(R(s), 1e-10);
        if rk < min_rank
            min_rank = rk;
            min_s = s;
        end
    end
end

fprintf('Minimum rank found on grid: %d at s = %+0.3e %+0.3ej\n', ...
    min_rank, real(min_s), imag(min_s));

if min_rank < normal_rank
    fprintf('=> Rank drop detected: finite transmission zero exists near that s.\n');
else
    fprintf('=> No rank drop detected on grid => no finite transmission zeros (consistent).\n');
end

% Cross-check with MATLAB's tzero (optional but nice)
tz = tzero(ss(A,B,C,D));
fprintf('tzero() returned %d transmission zeros.\n', numel(tz));
if isempty(tz)
    fprintf('=> Transmission zeros set is empty.\n');
else
    disp('Transmission zeros from tzero:'); disp(tz);
end
fprintf('=============================================\n');

