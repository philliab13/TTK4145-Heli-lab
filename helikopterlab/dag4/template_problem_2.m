% TTK4135 - Helicopter lab
% Hints/template for problem 2.
% Updated spring 2018, Andreas L. Fl�ten

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

A_k_4=eye(6)-delta_t*A;
B_k_4=delta_t*B;

x0=[pi; 0; 0; 0; 0.1; 0];
mx = size(A_k_4,2); % Number of states (number of columns in A)
mu = size(B_k_4,2);

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
Aeq = gen_aeq(A_k_4, B_k_4, N, mx,mu);
beq = zeros(N*mx,1);             % Generate b
beq(1:mx) = A_k_4*x0;

[vlb,vub]  = gen_constraints(N,N, xl, xu, ul, uu);


q1=1;
q2=1;
lambda_f=0;

z0_2=zeros(N*8,1);

[k_1_1_1, klm]=nonlinear_constraints(z,N,alpha,beta,lambda_t);



cost_function= @(z) sum((z(7:8:end) - lambda_f).^2) + q1*sum(z(8:8:end).^2) +q2*sum(z(9:8:end).^2);


[z_opt, fval]=fmincon(cost_function, z0_2, [], [], Aeq, beq, vlb, vub, @(z) nonlinear_constraints(z,N,alpha,beta,lambda_t));










function [c,ceq] = nonlinear_constraints(z,N,alpha,beta, lambda_t)
    c=zeros(N,1);
    ceq=[];
    hele_lambda=z(1:6:N*6);
    hele_elev=z(5:6:N*6);

    

    for k=1:N
        lambda_k=hele_lambda(k);
        e_k=hele_elev(k);
        c(k)=alpha * exp(-beta * (lambda_k-lambda_t).^2) - e_k;
    end
end



