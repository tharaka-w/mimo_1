
clear; clc; close all;


% 1) Nominal parameters (from proposal)
C1  = 3.0e6;      % [J/°C]
C2  = 2.5e6;      % [J/°C]
R1o = 0.010;      % [°C/W]
R2o = 0.012;      % [°C/W]
R12 = 0.020;      % [°C/W]
Ebat = 13.5 * 3.6e6;   % 13.5 kWh -> Joules [J]

% 2) State-space model for u -> y
% x = [T1; T2; s],  u = [Q1; Q2; Pbat],  y = [T1; T2; s]
A = [ -(1/(C1*R1o) + 1/(C1*R12)),    1/(C1*R12),               0;
       1/(C2*R12),                 -(1/(C2*R2o) + 1/(C2*R12)), 0;
       0,                           0,                         0 ];

B = [ 1/C1,   0,      0;
      0,     1/C2,    0;
      0,      0,   -1/Ebat ];

C = eye(3);
D = zeros(3,3);

G = ss(A,B,C,D);
G.InputName  = {'Q1','Q2','Pbat'};
G.OutputName = {'T1','T2','s'};


% 3) Poles and a warning about Gramians
p = pole(G);
disp('--- Poles of (u->y) plant ---');
disp(p);

% For continuous-time Gramians, the system must be asymptotically stable.
% Here we have an integrator at 0 (SOC state), so standard infinite-horizon
% Gramians diverge. Fix: compute Gramians on the stable thermal subsystem
% (T1,T2) only, which is what most instructors intend for this model.
if any(real(p) >= 0)
    fprintf(['\nNOTE: System is not asymptotically stable (pole at 0 due to SOC integrator).\n',...
             '=> Standard infinite-horizon Gramians Wc, Wo do NOT exist for the full 3-state plant.\n',...
             '=> Proceeding with the stable thermal 2-state subsystem (T1,T2) for Part 2.\n\n']);
end


% 4) Extract stable thermal subsystem (states T1,T2; inputs Q1,Q2; outputs T1,T2)
A_th = A(1:2,1:2);
B_th = B(1:2,1:2);
C_th = eye(2);
D_th = zeros(2,2);

Gth = ss(A_th, B_th, C_th, D_th);
Gth.InputName  = {'Q1','Q2'};
Gth.OutputName = {'T1','T2'};

disp('--- Thermal subsystem matrices used for Gramians ---');
disp('A_th ='); disp(A_th);
disp('B_th ='); disp(B_th);
disp('C_th ='); disp(C_th);
disp('D_th ='); disp(D_th);


% 5) Compute controllability and observability Gramians
Wc = gram(Gth,'c');   % controllability Gramian
Wo = gram(Gth,'o');   % observability Gramian


resWo = A_th'*Wo + Wo*A_th + C_th'*C_th;
norm(resWo,'fro')/norm(C_th'*C_th,'fro')

resWc = A_th*Wc + Wc*A_th' + B_th*B_th';
norm(resWc,'fro')/norm(B_th*B_th','fro')

disp('--- Controllability Gramian Wc (thermal subsystem) ---');
disp(Wc);

disp('--- Observability Gramian Wo (thermal subsystem) ---');
disp(Wo);

% 6) Eigen-decomposition and most/least directions
% Wc, Wo are symmetric PSD -> use eig, then sort by eigenvalue
[Vc, Dc] = eig(Wc);
evc = diag(Dc);

[Vo, Do] = eig(Wo);
evo = diag(Do);

% Sort descending (largest = "most", smallest = "least")
[evc_sorted, idxc] = sort(evc, 'descend');
Vc_sorted = Vc(:, idxc);

[evo_sorted, idxo] = sort(evo, 'descend');
Vo_sorted = Vo(:, idxo);

% Most/least controllable directions
v_c_most  = Vc_sorted(:,1);
lam_c_most = evc_sorted(1);

v_c_least  = Vc_sorted(:,end);
lam_c_least = evc_sorted(end);

% Most/least observable directions
v_o_most  = Vo_sorted(:,1);
lam_o_most = evo_sorted(1);

v_o_least  = Vo_sorted(:,end);
lam_o_least = evo_sorted(end);


disp('Part 2 Results (Thermal subsystem T1,T2)');

fprintf('\nControllability eigenvalues (descending):\n');
disp(evc_sorted);

fprintf('Most controllable direction (eigvec of largest eig):\n');
disp(v_c_most);
fprintf('Largest eig (most controllable): %.6e\n', lam_c_most);

fprintf('\nLeast controllable direction (eigvec of smallest eig):\n');
disp(v_c_least);
fprintf('Smallest eig (least controllable): %.6e\n', lam_c_least);

fprintf('\nObservability eigenvalues (descending):\n');
disp(evo_sorted);

fprintf('Most observable direction (eigvec of largest eig):\n');
disp(v_o_most);
fprintf('Largest eig (most observable): %.6e\n', lam_o_most);

fprintf('\nLeast observable direction (eigvec of smallest eig):\n');
disp(v_o_least);
fprintf('Smallest eig (least observable): %.6e\n', lam_o_least);


% 7) Optional: normalized directions (unit vectors)
v_c_most  = v_c_most  / norm(v_c_most);
v_c_least = v_c_least / norm(v_c_least);
v_o_most  = v_o_most  / norm(v_o_most);
v_o_least = v_o_least / norm(v_o_least);

disp('--- Unit-norm directions (optional, easier to report) ---');
disp('v_c_most_unit  ='); disp(v_c_most);
disp('v_c_least_unit ='); disp(v_c_least);
disp('v_o_most_unit  ='); disp(v_o_most);
disp('v_o_least_unit ='); disp(v_o_least);

fprintf('\nDONE.\n');

