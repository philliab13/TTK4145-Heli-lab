% TTK4135 - Helicopter lab
% Hints/template for problem 2.
% Updated spring 2018, Andreas L. Flåten

%% Initialization and model definition
init01; % Change this to the init file corresponding to your helicopter
Day3;
% Discrete time system model. x = [lambda r p p_dot]'
delta_t	= 0.25; % sampling time
A1 = eye(4) + delta_t*[[0 1 0 0];[0 0 -K_2 0]; [0 0 0 1];[0 0 -K_1*K_pp -K_1*K_pd]];
B1 = delta_t*[0 0 0 K_1*K_pp]';

% Number of states and inputs
mx = size(A1,2); % Number of states (number of columns in A)
mu = size(B1,2); % Number of inputs(number of columns in B)

% Initial values
x1_0 = pi;                               % Lambda
x2_0 = 0;                               % r
x3_0 = 0;                               % p
x4_0 = 0;                               % p_dot
x0 = [x1_0 x2_0 x3_0 x4_0]';           % Initial values

% Time horizon and initialization
N  = 100;                                  % Time horizon for states
M  = N;                                 % Time horizon for inputs
z  = zeros(N*mx+M*mu,1);                % Initialize z for the whole horizon
z0 = z;                                 % Initial value for optimization

% Bounds
ul 	    = -(60*pi)/360;                   % Lower bound on control
uu 	    = (60*pi)/360;                   % Upper bound on control

xl      = -Inf*ones(mx,1);              % Lower bound on states (no bound)
xu      = Inf*ones(mx,1);               % Upper bound on states (no bound)
xl(3)   = ul;                           % Lower bound on state x3
xu(3)   = uu;                           % Upper bound on state x3

% Generate constraints on measurements and inputs
[vlb,vub]       = gen_constraints(N,M, xl, xu, ul, uu); % hint: gen_constraints
vlb(N*mx+M*mu)  = 0;                    % We want the last input to be zero
vub(N*mx+M*mu)  = 0;                    % We want the last input to be zero

% Generate the matrix Q and the vector c (objecitve function weights in the QP problem) 
Q1 = zeros(mx,mx);
Q1(1,1) = 1;                            % Weight on state x1
Q1(2,2) = 0;                            % Weight on state x2
Q1(3,3) = 0;                            % Weight on state x3
Q1(4,4) = 0;                            % Weight on state x4
P1 = 1; 
% Q1 - Weight on states (mx*mx matrix)                       
% P1 - Weight on inputs (mu*mu matrix)                          
% N  - Time horizon for states
% M  - Time horizon for inputs% Weight on input
Q = gen_q(Q1, P1, N,M);                                  % Generate Q, hint: gen_q
c = 0 ;                                  % Generate c, this is the linear constant term in the QP

%% Generate system matrixes for linear model

Aeq = gen_aeq(A1, B1, N, mx,mu);             % Generate A, hint: gen_aeq
beq = zeros(N*mx,1);             % Generate b
beq(1:mx) = A1*x0;


%% Solve QP problem with linear model
tic
[z,lambda] = quadprog(Q,[],[],[],Aeq,beq,vlb,vub) ; % hint: quadprog. Type 'doc quadprog' for more info 
t1=toc;





% Calculate objective value
phi1 = 0.0;
PhiOut = zeros(N*mx+M*mu,1);
for i=1:N*mx+M*mu
  phi1=phi1+Q(i,i)*z(i)*z(i);
  PhiOut(i) = phi1;
end

%% Extract control inputs and states
u  = [z(N*mx+1:N*mx+M*mu);z(N*mx+M*mu)]; % Control input from solution
x1 = [x0(1);z(1:mx:N*mx)];              % State x1 from solution
x2 = [x0(2);z(2:mx:N*mx)];              % State x2 from solution
x3 = [x0(3);z(3:mx:N*mx)];              % State x3 from solution
x4 = [x0(4);z(4:mx:N*mx)];              % State x4 from solution

num_variables = 5/delta_t;
zero_padding = zeros(num_variables,1);
unit_padding  = ones(num_variables,1);

u   = [zero_padding; u; zero_padding];

x1  = [pi*unit_padding; x1; zero_padding];
x2  = [zero_padding; x2; zero_padding];
x3  = [zero_padding; x3; zero_padding];
x4  = [zero_padding; x4; zero_padding];

x_vector=[x1 x2 x3 x4];

%% Plotting
t = 0:delta_t:delta_t*(length(u)-1);
ts_x_vector = timeseries(x_vector, t);
ts_u = timeseries(u, t);


figure(2)
subplot(511)
stairs(t,u),grid
ylabel('u')
subplot(512)
plot(t,x1,'m',t,x1,'mo'),grid
ylabel('lambda')
subplot(513)
plot(t,x2,'m',t,x2','mo'),grid
ylabel('r')
subplot(514)
plot(t,x3,'m',t,x3,'mo'),grid
ylabel('p')
subplot(515)
plot(t,x4,'m',t,x4','mo'),grid
xlabel('tid (s)'),ylabel('pdot')

ts_lambda_ref = timeseries(u,t);


A= [
0 1 0 0 0 0;   
0 0 -K_2 0 0 0;
0 0 0 1 0 0;
0 0 K_1*K_pp -K_1*K_pd 0 0;
0 0 0 0 0 1;
0 0 0 0 -K_3*K_ep -K_3*K_ed;
];

B = [
0 0;
0 0;
0 0;
K_1*K_pp 0;
0 0;
0 K_3*K_ed
];

x0=[pi; 0; 0; 0; 0.1; 0];
mx = size(A,2); % Number of states (number of columns in A)
mu = size(B,2);
nz = mx + mu; 
%Vi er her sjekk at elevetaion cons er good
ul 	    = -(60*pi)/360;               % Lower bound on control
uu 	    = (60*pi)/360;                   % Upper bound on control

xl      = -Inf*ones(mx,1);              % Lower bound on states (no bound)
xu      = Inf*ones(mx,1);               % Upper bound on states (no bound)
xl(3)   = ul;                           % Lower bound on state x3
xu(3)   = uu;
N=40;
beta=20;
alpha=0.2;
lambda_t=2*pi/3;
Aeq = gen_aeq(A, B, N, mx,mu);
beq = zeros(N*mx,1);             % Generate b
beq(1:mx) = A*x0;
dt=delta_t;

[vlb,vub]  = gen_constraints(N,N, xl, xu, ul, uu);



q1=1;
q2=1;
lambda_f=0;

z0_2=zeros(N*8,1);




cost_function = @(z) cost_func(z, N, lambda_f, q1, q2);
nonlincon     = @(z) nonlinear_constraints(z, N, alpha, beta, lambda_t);

%% Initial Guess for Decision Variable z
% z has dimension 8*N. A simple guess is zero.
z0 = zeros(nz*N, 1);

%% Solve the Optimization Problem with fmincon
options = optimoptions('fmincon','Algorithm','sqp','Display','iter');
[z_opt, fval] = fmincon(cost_function, z0, [], [], Aeq, beq, vlb, vub, nonlincon, options);

%% Extract Optimized States and Inputs
%% Extract Optimized States and Inputs (if not already done)
% Reshape z_opt into an 8-by-N matrix, where each column corresponds to one time step.
Z = reshape(z_opt, nz, N);
X = Z(1:mx, :);      % states for k = 1,...,N
U = Z(mx+1:end, :);  % inputs for k = 1,...,N

% Prepend the initial state (x0) to build a full trajectory.
X_full = [x0, X];  % Now X_full is 6 x (N+1)

%% Plot All States
t = 0:dt:dt*N;  % time vector for states, length N+1

figure;
stateNames = {'\lambda', 'r', 'p', '\dot{p}', 'e', '\dot{e}'};
for i = 1:mx
    subplot(mx, 1, i);
    plot(t, X_full(i,:), 'o-','LineWidth',1.5);
    ylabel(stateNames{i});
    grid on;
    if i == 1
        title('State Trajectories');
    end
    if i == mx
        xlabel('Time (s)');
    end
end


%% --- Helper Functions ---

function cost = cost_func(z, N, lambda_f, q1, q2)
    % Reshape decision vector z into an 8-by-N matrix.
    Z = reshape(z, 8, N);
    % Extract state variables:
    % - lambda is row 1,
    % - p (pitch) is row 3,
    % - elevation e is row 5.
    lambda_vals = Z(1, :);
    p_vals      = Z(3, :);
    elev_vals   = Z(5, :);
    
    % Compute the cost as a sum over the horizon.
    cost = sum((lambda_vals - lambda_f).^2) + q1*sum(p_vals.^2) + q2*sum(elev_vals.^2);
end

function [c, ceq] = nonlinear_constraints(z, N, alpha, beta, lambda_t)
    % Reshape decision vector z into an 8-by-N matrix.
    Z = reshape(z, 8, N);
    % Extract lambda (row 1) and elevation (row 5)
    lambda_vals = Z(1, :);
    elev_vals   = Z(5, :);
    
    % Inequality constraint at each time step:
    %   alpha*exp(-beta*(lambda - lambda_t)^2) - elevation <= 0.
    c = alpha*exp(-beta*(lambda_vals - lambda_t).^2) - elev_vals;
    c = c';  % ensure column vector
    ceq = [];
end



