%% ==========================================================
% Assignment 3 — Part 3 (FIXED + COPY/PASTE READY)
% Balanced realization + Hankel singular values (thermal subsystem)
%
% Thermal subsystem:
%   states  x_T = [T1; T2]
%   inputs  u_T = [Q1; Q2]
%   outputs y_T = [T1; T2]
%
% Prints:
%   - A_th, B_th, C_th, D_th
%   - Hankel singular values (HSV)
%   - Balanced realization matrices (Ab,Bb,Cb,Db)
%   - Most/least important balanced state (largest/smallest HSV)
%   - Optional r=1 balanced truncation + error bound
%
% Requires: Control System Toolbox
%% ==========================================================
clear; clc; close all;

fprintf('==========================================================\n');
fprintf('Assignment 3 — Part 3: Balanced Realization + Hankel SVs\n');
fprintf('Using stable thermal subsystem (T1,T2)\n');
fprintf('==========================================================\n\n');

%% --------------------------
% 1) Nominal parameters (from proposal)
C1  = 3.0e6;      % [J/°C]
C2  = 2.5e6;      % [J/°C]
R1o = 0.010;      % [°C/W]
R2o = 0.012;      % [°C/W]
R12 = 0.020;      % [°C/W]

%% --------------------------
% 2) Thermal subsystem matrices
A_th = [ -(1/(C1*R1o) + 1/(C1*R12)),    1/(C1*R12);
          1/(C2*R12),                 -(1/(C2*R2o) + 1/(C2*R12)) ];

B_th = [ 1/C1, 0;
         0,   1/C2 ];

C_th = eye(2);
D_th = zeros(2,2);

sys = ss(A_th, B_th, C_th, D_th);
sys.InputName  = {'Q1','Q2'};
sys.OutputName = {'T1','T2'};

disp('--- Thermal subsystem matrices ---');
disp('A_th ='); disp(A_th);
disp('B_th ='); disp(B_th);
disp('C_th ='); disp(C_th);
disp('D_th ='); disp(D_th);

%% --------------------------
% 3) Hankel singular values (HSV)
hsv = hsvd(sys);

disp('--- Hankel singular values (descending) ---');
disp(hsv);

fprintf('Largest HSV  = %.6e\n', hsv(1));
fprintf('Smallest HSV = %.6e\n', hsv(end));

%% --------------------------
% 4) Balanced realization
[sysb, g] = balreal(sys);    % g are HSVs in balanced coordinates

[Ab, Bb, Cb, Db] = ssdata(sysb);

disp('--- Balanced realization (sysb) ---');
disp('Ab ='); disp(Ab);
disp('Bb ='); disp(Bb);
disp('Cb ='); disp(Cb);
disp('Db ='); disp(Db);

disp('--- HSV returned by balreal (g) ---');
disp(g);

% Most/least important balanced state indices
[~, idx_max] = max(g);
[~, idx_min] = min(g);

fprintf('\nBalanced-state importance (by HSV):\n');
fprintf('Most important balanced state index  = %d (largest HSV)\n', idx_max);
fprintf('Least important balanced state index = %d (smallest HSV)\n', idx_min);

%% --------------------------
% 5) Optional: Balanced truncation to order r=1 and error bound
r = 1;

% balred requires stable sys (true here)
sysr = balred(sys, r);

% BT error bound: ||G - G_r||_inf <= 2 * sum_{i=r+1}^n hsv_i
err_bound = 2 * sum(hsv(r+1:end));

disp('--- Optional: Balanced truncation to order r=1 ---');
disp('Reduced model sysr (tf form):');
disp(tf(sysr));

fprintf('BT H-infinity error bound: ||G - G_r||_inf <= %.6e\n', err_bound);

disp(' ');
disp('==========================================================');
disp('DONE: Part 3 outputs printed above.');
disp('==========================================================');
