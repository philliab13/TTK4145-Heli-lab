% TTK4135 - Helicopter lab
% Hints/template for problem 2.
% Updated spring 2018, Andreas L. Flåten

%% Initialization and model definition
init01; % Change this to the init file corresponding to your helicopter
day_4_K_calc;

%% --- 4-State Model (for QP solution) ---
% Discrete time system model. x = [lambda r p p_dot]'
delta_t = 0.25; % Sampling time
A1 = eye(4) + delta_t * [0 1 0 0; 0 0 -K_2 0; 0 0 0 1; 0 0 -K_1*K_pp -K_1*K_pd];
B1 = delta_t * [0 0 0 K_1*K_pp]';

mx = size(A1,2);  % number of states (4)
mu = size(B1,2);  % number of inputs (1)

% Initial state for 4-state model
x1_0 = pi;   % lambda
x2_0 = 0;    % r
x3_0 = 0;    % p
x4_0 = 0;    % p_dot
x0 = [x1_0; x2_0; x3_0; x4_0];

% Time horizon for QP solution
N = 100;      % States horizon
M = N;        % Inputs horizon
z = zeros(N*mx + M*mu, 1);  % Decision variable vector
z0 = z;  % initial guess

% Bounds for 4-state model
ul_bound = -(60*pi)/360;  % lower bound on control
uu_bound = (60*pi)/360;   % upper bound on control
xl = -Inf*ones(mx,1); 
xu = Inf*ones(mx,1);
% For example, bound state x3 (p)
xl(3) = ul_bound;
xu(3) = uu_bound;

[vlb, vub] = gen_constraints(N, M, xl, xu, ul_bound, uu_bound);
% Force the last input to be zero:
vlb(N*mx+M*mu) = 0;
vub(N*mx+M*mu) = 0;

% Generate cost weighting matrices
Q1 = zeros(mx,mx);
Q1(1,1) = 1;  % Weight on lambda
Q1(2,2) = 0;
Q1(3,3) = 0;
Q1(4,4) = 0;
P1 = 1;
Q = gen_q(Q1, P1, N, M);  % Using the provided gen_q
c = 0;  % constant term in QP cost

% Generate equality constraints for the 4-state model
Aeq = gen_aeq(A1, B1, N, mx, mu);
beq = zeros(N*mx, 1);
beq(1:mx) = A1 * x0;

%% Solve QP problem for the 4-state model
tic
[z, lambda_qp] = quadprog(Q, [], [], [], Aeq, beq, vlb, vub);
t_qp = toc;

% (Optional) Compute QP objective value
phi1 = 0.0;
PhiOut = zeros(length(z), 1);
for i = 1:length(z)
    phi1 = phi1 + Q(i,i) * z(i)^2;
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


%% --- Augmented 6-State Model with Elevation Dynamics ---
% Continuous model:
%   x = [lambda; r; p; p_dot; e; e_dot]
%   u = [p_c; e_c]

A = [0     1       0         0     0         0;   
     0     0      -K_2       0     0         0;
     0     0       0         1     0         0;
     0     0   K_1*K_pp  -K_1*K_pd 0         0;
     0     0       0         0     0         1;
     0     0       0         0 -K_3*K_ep -K_3*K_ed];
 
B = [0      0;
     0      0;
     0      0;
  K_1*K_pp  0;
     0      0;
     0  K_3*K_ed];
dt=0.25;

x0_aug = [pi; 0; 0; 0; 0.1; 0];  % initial state for augmented model
mx_aug = size(A,2);  % 6 states
mu_aug = size(B,2);  % 2 inputs
nz = mx_aug + mu_aug;  % 8 variables per time step

% Discretize the augmented system (as you did with the 4-state system)
A_d_aug = eye(mx_aug) + dt * A;
B_d_aug = dt * B;

% Update bounds for the augmented model.
% For example, we bound the elevation (state x5).
ul_bound = -(60*pi)/360;
uu_bound = (60*pi)/360;
xl_aug = -Inf*ones(mx_aug,1);
xu_aug = Inf*ones(mx_aug,1);
xl_aug(3) = ul_bound;
xu_aug(3) = uu_bound;

N_aug = 40;  % shorter horizon for augmented model

% Use the discretized matrices in generating equality constraints.
Aeq_aug = gen_aeq(A_d_aug, B_d_aug, N_aug, mx_aug, mu_aug);
beq_aug = zeros(N_aug*mx_aug, 1);
beq_aug(1:mx_aug) = A_d_aug * x0_aug;

[vlb_aug, vub_aug] = gen_constraints(N_aug, N_aug*2, xl_aug, xu_aug, ul_bound, uu_bound);

% Cost function and nonlinear constraint parameters for augmented model
q1 = 1;
q2 = 1;
lambda_f = 0;
alpha = 0.2;
beta = 20;
lambda_t = 2*pi/3;

% Cost function: penalize deviation in lambda (state 1), pitch p (state 3),
% and elevation e (state 5)
cost_function = @(z) cost_func(z, N_aug, lambda_f, q1, q2);

% Initial guess for augmented model optimization
z0_aug = zeros(nz * N_aug, 1);
z_opt=z0_aug;
for i=1:N_aug
    nonlincon = @(z) nonlinear_constraints(z_opt, N_aug, alpha, beta, lambda_t, i);
    [z_opt, fval] = fmincon(cost_function, z0_aug, [], [], Aeq_aug, beq_aug, vlb_aug, vub_aug, nonlincon);

    
end


%% Solve the Optimization Problem with fmincon for the augmented model

%% Extract Optimized States and Inputs for the Augmented Model
Z = reshape(z_opt, nz, N_aug);
X_aug = Z(1:mx_aug, :);      % states for k = 1,...,N_aug
U_aug = Z(mx_aug+1:end, :);  % inputs for k = 1,...,N_aug

% Prepend the initial state to form a full state trajectory.
X_full_aug = [x0_aug, X_aug];  % size is 6 x (N_aug+1)
x_aug_timeseries=timeseries(X_aug);
U_aug_timeseries=timeseries(U_aug);

%% Plot All States for the Augmented Model
t_aug = 0:dt:dt*N_aug;  % time vector (length N_aug+1)
figure(2)
stateNames = {'\lambda', 'r', 'p', 'p\_dot', 'e', 'e\_dot'};
for i = 1:mx_aug
    subplot(mx_aug, 1, i);
    plot(t_aug, X_full_aug(i,:), 'o-', 'LineWidth', 1.5);
    ylabel(stateNames{i});
    grid on;
    if i == 1
        title('State Trajectories for Augmented Model');
    end
    if i == mx_aug
        xlabel('Time (s)');
    end
end


%% --- Helper Functions ---



function cost = cost_func(z, N, lambda_f, q1, q2)
% Cost function for the augmented model.
% Assumes the decision variable z is organized in blocks of 8:
%   [lambda; r; p; p_dot; e; e_dot; p_c; e_c]
Z = reshape(z, [], N);
lambda_vals = Z(1, :);
p_vals      = Z(7, :);
elev_vals   = Z(8, :);
cost = sum((lambda_vals - lambda_f).^2 + q1 * (p_vals.^2) + q2 * (elev_vals.^2));
end

function [c, ceq] = nonlinear_constraints(z, N, alpha, beta, lambda_t, i)
% Nonlinear inequality constraint:
%   alpha * exp(-beta*(lambda - lambda_t)^2) - elevation <= 0
% for each time step.
Z = reshape(z, [], N);
lambda_vals = Z(1, :);
lambda_vals = lambda_vals(1,i);
elev_vals   = Z(5, :);
elev_vals =elev_vals(1,i);
c = alpha * exp(-beta * (lambda_vals - lambda_t).^2) - elev_vals;
c = c';  
ceq = [];
end
