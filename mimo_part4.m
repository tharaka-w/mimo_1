
clear; clc; close all;


% 1) Nominal parameters (from proposal)
C1  = 3.0e6;
C2  = 2.5e6;
R1o = 0.010;
R2o = 0.012;
R12 = 0.020;

Ebat = 13.5 * 3.6e6;     % 13.5 kWh -> Joules
COP  = 3.0;
k1   = 1/COP;
k2   = 1/COP;

% z-penalty weights (set as needed)
lam1 = 1;
lam2 = 1;
lamb = 1;


% 2) Build u->y subsystem (thermal only), then BALANCE it
%    states xT = [T1; T2], inputs uT = [Q1; Q2], outputs yT = [T1; T2]

A_T = [ -(1/(C1*R1o) + 1/(C1*R12)),    1/(C1*R12);
         1/(C2*R12),                 -(1/(C2*R2o) + 1/(C2*R12)) ];

B_T = [ 1/C1, 0;
        0,   1/C2 ];

C_T = eye(2);
D_T = zeros(2,2);

sysT = ss(A_T,B_T,C_T,D_T);
sysT.InputName  = {'Q1','Q2'};
sysT.OutputName = {'T1','T2'};

% Balanced realization + transformation
% balreal gives balanced system, and returns T such that:
%   Ab = inv(T)*A*T , Bb = inv(T)*B , Cb = C*T
[sysTb, hsv, Tbal] = balreal(sysT);

Ab = sysTb.A;  Bb = sysTb.B;  Cb = sysTb.C;  Db = sysTb.D; %#ok<NASGU>

fprintf('--- u->y thermal subsystem is balanced ---\n');
disp('Hankel singular values (descending) ='); disp(hsv);
disp('Balancing transform Tbal (x = Tbal * x_bal) ='); disp(Tbal);


% 3) Build FULL generalized plant matrices in ORIGINAL coordinates
%    but with only the THERMAL states x = [T1; T2] (stable)
%    Inputs:
%      u = [Q1; Q2; Pbat]   (3)
%      w = [Tout; qint1; qint2; Ppv; rT1; rT2; nT1; nT2]  (8)
%    Outputs:
%      y = [T1_meas; T2_meas] (2)
%      z = [e1; e2; Pgrid; sqrt(lam1)Q1; sqrt(lam2)Q2; sqrt(lamb)Pbat] (6)
%
% NOTE: SOC state s is excluded here because balancing/Gramians require stability.
%       If instructor demands s included, we must do finite-horizon balancing.

% State: [T1; T2]
A = A_T;

% Control inputs u = [Q1; Q2; Pbat]
B2 = zeros(2,3);
B2(:,1:2) = B_T;       % Q1,Q2 affect thermal states
% Pbat does not affect temperatures => column 3 = 0

% Exogenous inputs w = [Tout; qint1; qint2; Ppv; rT1; rT2; nT1; nT2]
nw = 8;
B1 = zeros(2,nw);
B1(1,1) = 1/(C1*R1o);  % Tout -> T1
B1(1,2) = 1/C1;        % qint1 -> T1
B1(2,1) = 1/(C2*R2o);  % Tout -> T2
B1(2,3) = 1/C2;        % qint2 -> T2
% others do not enter state dynamics

% Measured outputs y = [T1_meas; T2_meas]
C2 = eye(2);
D22 = zeros(2,3);

D21 = zeros(2,nw);
D21(1,7) = 1;          % nT1 added to T1 measurement
D21(2,8) = 1;          % nT2 added to T2 measurement

% Performance outputs z
% z = [e1; e2; Pgrid; sqrt(lam1)Q1; sqrt(lam2)Q2; sqrt(lamb)Pbat]
C1 = zeros(6,2);
C1(1,1) = 1;           % e1 includes +T1
C1(2,2) = 1;           % e2 includes +T2
% row3 Pgrid has no x dependence here

D11 = zeros(6,nw);
D11(1,5) = -1;         % e1 = T1 - rT1
D11(2,6) = -1;         % e2 = T2 - rT2
D11(3,4) = -1;         % Pgrid includes -Ppv

D12 = zeros(6,3);
D12(3,:) = [k1, k2, -1];                 % Pgrid = k1 Q1 + k2 Q2 - Pbat - Ppv
D12(4,:) = [sqrt(lam1), 0, 0];           % sqrt(lam1)*Q1
D12(5,:) = [0, sqrt(lam2), 0];           % sqrt(lam2)*Q2
D12(6,:) = [0, 0, sqrt(lamb)];           % sqrt(lamb)*Pbat


% 4) Transform the ENTIRE four-block plant to BALANCED coordinates
%    x = Tbal * x_bal
%    => A_bal = inv(Tbal)*A*Tbal
%       B1_bal = inv(Tbal)*B1
%       B2_bal = inv(Tbal)*B2
%       C1_bal = C1*Tbal
%       C2_bal = C2*Tbal
%    D blocks unchanged

Ti = inv(Tbal);

A_bal  = Ti*A*Tbal;
B1_bal = Ti*B1;
B2_bal = Ti*B2;
C1_bal = C1*Tbal;
C2_bal = C2*Tbal;


% 5) Print the balanced four-block realization matrices
fprintf('\n==========================================================\n');
fprintf('FOUR-BLOCK REALIZATION IN BALANCED COORDINATES (thermal)\n');
fprintf('Inputs:  [w; u],  Outputs: [z; y]\n');
fprintf('w = [Tout qint1 qint2 Ppv rT1 rT2 nT1 nT2]^T\n');
fprintf('u = [Q1 Q2 Pbat]^T\n');
fprintf('z = [e1 e2 Pgrid sqrt(lam1)Q1 sqrt(lam2)Q2 sqrt(lamb)Pbat]^T\n');
fprintf('y = [T1_meas T2_meas]^T\n');
fprintf('==========================================================\n\n');

disp('A_bal =');  disp(A_bal);
disp('B1_bal ='); disp(B1_bal);
disp('B2_bal ='); disp(B2_bal);

disp('C1_bal ='); disp(C1_bal);
disp('D11 =');    disp(D11);
disp('D12 =');    disp(D12);

disp('C2_bal ='); disp(C2_bal);
disp('D21 =');    disp(D21);
disp('D22 =');    disp(D22);

% 6) Build generalized plant P_bal : [w;u] -> [z;y]
B = [B1_bal, B2_bal];
C = [C1_bal; C2_bal];
D = [D11, D12;
     D21, D22];

P_bal = ss(A_bal, B, C, D);

% Name signals
w_names = {'Tout','qint1','qint2','Ppv','rT1','rT2','nT1','nT2'};
u_names = {'Q1','Q2','Pbat'};
z_names = {'e1','e2','Pgrid','sqrt(lam1)Q1','sqrt(lam2)Q2','sqrt(lamb)Pbat'};
y_names = {'T1_meas','T2_meas'};

P_bal.InputName  = [w_names, u_names];
P_bal.OutputName = [z_names, y_names];

fprintf('\n--- Generalized plant P_bal built: [w;u] -> [z;y] ---\n');
fprintf('Size(P_bal): outputs = %d, inputs = %d\n', size(P_bal,1), size(P_bal,2));


% 7) Extract the four blocks (for reporting)
Pzw = P_bal(1:6, 1:nw);
Pzu = P_bal(1:6, nw+1:nw+3);
Pyw = P_bal(7:8, 1:nw);
Pyu = P_bal(7:8, nw+1:nw+3);

fprintf('\n--- Four-block partitions (balanced coordinates) ---\n');
fprintf('Pzw (w->z): %dx%d\n', size(Pzw,1), size(Pzw,2));
fprintf('Pzu (u->z): %dx%d\n', size(Pzu,1), size(Pzu,2));
fprintf('Pyw (w->y): %dx%d\n', size(Pyw,1), size(Pyw,2));
fprintf('Pyu (u->y): %dx%d\n', size(Pyu,1), size(Pyu,2));




