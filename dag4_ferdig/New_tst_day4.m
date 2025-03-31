% TTK4135 - Helicopter lab
% Hints/template for problem 2.
% Updated spring 2018, Andreas L. Flåten
clear all;
%% Initialization and model definition
init01; % Change this to the init file corresponding to your helicopter
day_4_K_calc;

% Initial state for 4-state model
x1_0 = pi;   % lambda
x2_0 = 0;    % r
x3_0 = 0;    % p
x4_0 = 0;    % p_dot
x0 = [x1_0; x2_0; x3_0; x4_0];


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
     0  K_3*K_ep];

%some values
dt=0.25;
N_aug = 40;  % shorter horizon for augmented model
lambda_f = 0;
alpha = 0.2;
beta = 20;
lambda_t = 2*pi/3;

%define area for zero padding
num_variables = 2/dt;
zero_padding = zeros(num_variables,1);
unit_padding  = ones(num_variables,1);


x0_aug = [pi; 0; 0; 0; 0; 0];  % initial state for augmented model
x=zeros(N_aug*6,1);
x(1:6)=x0_aug;
mx_aug = size(A,2);  % 6 states
mu_aug = size(B,2);  % 2 inputs
nz = mx_aug + mu_aug;  % 8 variables per time step



% Discretize the augmented system (as you did with the 4-state system)
A_d_aug = eye(mx_aug) + dt * A;
B_d_aug = dt * B;


% Use the discretized matrices in generating equality constraints.
Aeq_aug = gen_aeq(A_d_aug, B_d_aug, N_aug, mx_aug, mu_aug);
beq_aug = zeros(N_aug*mx_aug, 1);
beq_aug(1:mx_aug) = A_d_aug * x0_aug;

%create lower bound
vlb_aug=-inf*ones(N_aug*nz,1);
vlb_aug(3:mx_aug:mx_aug*N_aug)=-(60*pi)/360;

vlb_aug(mx_aug*N_aug+1:mu_aug:nz*N_aug)=-(60*pi)/360;

%create upper bound
vub_aug=inf*ones(N_aug*nz,1);
vub_aug(3:mx_aug:mx_aug*N_aug)=(60*pi)/360;

vub_aug(mx_aug*N_aug+1:mu_aug:nz*N_aug)=(60*pi)/360;


q1 = 1;
q2 = 1;

% Create Q
Q_aug = zeros(mx_aug,mx_aug);
Q_aug(1,1) = 1;  % Weight on lambda
Q_aug(2,2) = 0;
Q_aug(3,3) = 0;
Q_aug(4,4) = 0;
Q_aug(5,5) =0;
Q_aug(6,6)=0;

%create G
P_aug = [q1 0;
          0 q2];
G = gen_q(Q_aug, P_aug, N_aug, N_aug);

%create z
z_aug = zeros(nz * N_aug, 1);
z0_aug=z_aug;

% Cost function:
cost_function = @(z) 1/2*z_aug'*G*z_aug;

%non linear constraints
nonlincon = @(z_aug) nonlinear_constraints(z_aug, N_aug, alpha, beta, lambda_t, mx_aug);

%run fmincon
opt=optimoptions('fmincon','Algorithm','sqp', 'MaxFunEvals', 40000);
[z_opt, fval] = fmincon(cost_function, z0_aug, [], [], Aeq_aug, beq_aug, vlb_aug, vub_aug, nonlincon,opt);

    
%% Extract Optimized States and Inputs for the Augmented Model

u1=[z_opt(N_aug*mx_aug+1:mu_aug:nz*N_aug);z_opt(nz*N_aug-1)];
u2=[z_opt(N_aug*mx_aug+2:mu_aug:nz*N_aug);z_opt(nz*N_aug)];
x1=[x0_aug(1);z_opt(1:mx_aug:N_aug*mx_aug)];
x2=[x0_aug(2);z_opt(2:mx_aug:N_aug*mx_aug)];
x3=[x0_aug(3);z_opt(3:mx_aug:N_aug*mx_aug)];
x4=[x0_aug(4);z_opt(4:mx_aug:N_aug*mx_aug)];
x5=[x0_aug(5);z_opt(5:mx_aug:N_aug*mx_aug)];
x6=[x0_aug(6);z_opt(6:mx_aug:N_aug*mx_aug)];

u1= [zero_padding; u1; zero_padding];
u2=[zero_padding; u2; zero_padding];

x1  = [pi*unit_padding; x1; zero_padding];
x2  = [zero_padding; x2; zero_padding];
x3  = [zero_padding; x3; zero_padding];
x4  = [zero_padding; x4; zero_padding];
x5  = [zero_padding; x5; zero_padding];
x6  = [zero_padding; x6; zero_padding];


t=0:size(x1)-1;

x_vector=[x1 x2 x3 x4 x5 x6];
x_vector_trans = x_vector';
u_vec=[u1 u2];

x_aug_timeseries=timeseries(x_vector);
U_aug_timeseries=timeseries(u_vec);


%% Plot All States for the Augmented Model
figure(2)
subplot(811)
stairs(t,u1),grid
ylabel('u1')
subplot(812)
plot(t,x1,'m',t,x1,'mo'),grid
ylabel('lambda')
subplot(813)
plot(t,x2,'m',t,x2','mo'),grid
ylabel('r')
subplot(814)
plot(t,x3,'m',t,x3,'mo'),grid
ylabel('p')
subplot(815)
plot(t,x4,'m',t,x4','mo'),grid
xlabel('tid (s)'),ylabel('pdot')
subplot(816)
plot(t,x5,'m',t,x5','mo'),grid
xlabel('tid (s)'),ylabel('e')
subplot(817)
plot(t,x6,'m',t,x6','mo'),grid
xlabel('tid (s)'),ylabel('edot')
subplot(818)
stairs(t,u2),grid
ylabel('u2')



%% --- Helper Functions ---
function [c,ceq] = nonlinear_constraints(z,N,alpha,beta,lambda_t,nx_aug)
    c=zeros(N,1);
    for k=1:N
        c(k)=alpha*exp(-beta*(z(1+(k-1)*nx_aug)-lambda_t)^2)-z(5+(k-1)*nx_aug);
    end
    save('non_lin_con.mat', 'c');
    plot(c)
    ceq=[];
end