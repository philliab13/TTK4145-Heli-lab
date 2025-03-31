

%Oppgave 10.3.1 - 1
delta_t	= 0.25; % sampling time

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
mx_aug = size(A,2);  % 6 states
mu_aug = size(B,2);  % 2 inputs
nz = mx_aug + mu_aug;  % 8 variables per time step

A_d_aug = eye(mx_aug) + dt * A;
B_d_aug = dt * B;

q=[100 1 1 1 1 1];
r=[1 1];

Q_raaa=diag(q);
R_raaa=diag(r);


K=dlqr(A_d_aug,B_d_aug,Q_raaa,R_raaa);

%Oppgave 10.3.1 - 2
%Simulink   