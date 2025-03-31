
%Discrete system model
A1 = eye(4) + delta_t*[[0 1 0 0];[0 0 -K_2 0]; [0 0 0 1];[0 0 -K_1*K_pp -K_1*K_pd]];
B1 = delta_t*[0 0 0 K_1*K_pp]';


% Initial values
x1_0 = pi;                               % Lambda
x2_0 = 0;                               % r
x3_0 = 0;                               % p
x4_0 = 0;                               % p_dot
x0 = [x1_0 x2_0 x3_0 x4_0]';           % Initial values

% Time horizon and initialization
N  = 100;                                  % Time horizon for states

% Bounds
ul 	    = -(60*pi)/360;                   % Lower bound on control
uu 	    = (60*pi)/360;                   % Upper bound on control

% Generate constraints on measurements and inputs
[vlb,vub]       = gen_constraints(N,M, xl, xu, ul, uu); 

% Generate the matrix Q and the vector c (objecitve function weights in the QP problem) 
Q1 = zeros(mx,mx);
Q1(1,1) = 1;                            % Weight on state x1
Q1(2,2) = 0;                            % Weight on state x2
Q1(3,3) = 0;                            % Weight on state x3
Q1(4,4) = 0;                            % Weight on state x4
P1 = 12;                                % P1 - Weight on inputs                        
Q = gen_q(Q1, P1, N,M);                 % Generate Q, hint: gen_q
c = 0 ;                                 % Generate c

%% Generate system matrixes for linear model

Aeq = gen_aeq(A1, B1, N, mx,mu);             % Generate A, hint: gen_aeq
beq = zeros(N*mx,1);             % Generate b
beq(1:mx) = A1*x0;

%% Solve QP problem with linear model
[z,lambda] = quadprog(Q,[],[],[],Aeq,beq,vlb,vub) ; 

%To Simulink 
ts_lambda_ref = timeseries(u,t);
