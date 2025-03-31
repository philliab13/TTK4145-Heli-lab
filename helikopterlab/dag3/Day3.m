

%Oppgave 10.3.1 - 1
delta_t	= 0.25; % sampling time

A1 = eye(4) + delta_t*[[0 1 0 0];[0 0 -K_2 0]; [0 0 0 1];[0 0 -K_1*K_pp -K_1*K_pd]];
B1 = delta_t*[0 0 0 K_1*K_pp]';

q=[10 10 100 100];
r=[1];

Q_raaa=diag(q);
R_raaa=diag(r);


K=dlqr(A1,B1,Q_raaa,R_raaa);

