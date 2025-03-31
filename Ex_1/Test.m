syms L2 q1 q2 L1 dq1 dq2    

Ry = [cos(q2) 0 sin(q2); 0 1 0; -sin(q2) 0 cos(q2)];
Rz = [cos(q1) -sin(q1) 0; sin(q1) cos(q1) 0; 0 0 1];
r=[L2; 0; 0];

r_1_2_3=Rz*Ry*r;


r_B_A=[L1*cos(q1); L1*sin(q1);0];

r_1_2_3=Rz*Ry*r;
r_C_A=r_B_A+r_1_2_3;

q_BC=[0; dq2; 0];

v_A_B=Rz*q_BC;

omega_A_B=[0;0; dq1];
omega_B_C=[-dq2*sin(q1);dq2*cos(q1);dq1];

v_B_0=cross(omega_A_B,r_B_A);

v_C_0=v_B_0+cross(omega_B_C,Rz*Ry*r);