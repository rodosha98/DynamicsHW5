clear; clc;
% 3 Dofs manipularot
%Given:
syms theta1 theta2 d3 q1 q2 q3 l1 l2 t m1 m2 m3 real;    
% Initial constants
l1 = 2;
l2 = 1;
m1 = 5;
m2 = 2;
m3 = 1;
g0 = 9.82;
q1 = 1*sin(t);
q2 = 3*cos(2*t);
q3 = 2*sin(3*t);
q = [0, q1, q2, q3];

dq = [0, diff(q1,t), diff(q2,t), diff(q3,t)];
ddq = [0, diff(q1,t, 2), diff(q2,t, 2), diff(q3,t, 2)];



%% Forward Kinematics
O1 = trotz(q1)*transl(0, 0, l1)*trotx(pi/2);
O2 = O1*trotz(q2)*transl(l2, 0, 0)*troty(pi/2);
O3 = O2*transl(0, 0, q3);

O21 = trotz(q2)*transl(l2, 0, 0)*troty(pi/2);
O32 = transl(0, 0, q3);
%Rotation matrices
R_10 = rotz(q1)*rotx(pi/2);
R_20 = R_10*rotz(q2)*roty(pi/2);
R_30 = R_20;


R_21 = rotz(q2)*roty(pi/2);
R_32 = rotz(0);
R = [R_10, R_21, R_32];
R_0 = [R_10, R_20, R_30];


d_10 = [O1(1,4); O1(2,4); O1(3,4)];
d_20 = [O2(1,4); O2(2,4); O2(3,4)];
d_30 = [O3(1,4); O3(2,4); O3(3,4)];
d_21 = [O21(1,4); O21(2,4); O21(3,4)];
d_32 = [O32(1,4); O32(2,4); O32(3,4)];

r = [0, d_10, d_21, d_32];
rc = r/2;
%% Initial parameters
g = 9.8;
lc1 = l1/2;
lc2 = l2/2;
lc3 = d3/2;



% Unit vectors
k = [0 0 1]';
%unit z vector 
z0 = k;


% Initial vectors
w = zeros(3, 4, 'sym');
dw = zeros(3, 4, 'sym');
ddpc = zeros(3, 4, 'sym'); %Linear velocies of the centre 
ddp = zeros(3, 4, 'sym'); %Linear acceleration of the end 
Centr = zeros(3, 4, 'sym');
Coriol = zeros(3, 4, 'sym');


dpc = zeros(3, 4, 'sym'); %Linear velocies of the centre 
dp= zeros(3, 4, 'sym');


z = zeros(3, 4, 'sym');
z(:, 1)=[0; 0; 1];
B = zeros(3, 4, 'sym');

g0= -9.82*k;
G = zeros(3, 4, 'sym');
G(:, 1) = g0;

%%
for i = 2:4
    i
    %bounds
    al = 1+3*(i-2);
    bh = 3*(i-1);
    z(:, i) = R_0(:, bh);
    z0 = [0; 0; 1];
    B(:,i) = R_0(:, al:bh)*z(:, i-1);
    
    % Velocities
    %Prismatic
    if i==4
         w(:, i) = R(:, al:bh)'*w(:, i-1);
         dp(:, i) = R(:, al:bh)'*(dp(:, i-1) + dq(i)*z0) + cross(w(:, i), r(:, i));
         dpc(:, i) = dp(:, i) + cross(w(:, i), -rc(:, i)); 
    else
        %Revolute
        w(:, i) = R(:, al:bh)'*(w(:, i-1) + dq(i)*z0);
        dp(:, i) = R(:, al:bh)'*dp(:, i-1) + cross(w(:, i), r(:, i));
        dpc(:, i) = dp(:, i) + cross(w(:, i), -rc(:, i)); 
    end
   

    

    
    %Accelerations
    dw(:, i) = diff(w(:, i), t);
    %Prismatic
    if i==4
         dw(:, i) = R(:, al:bh)'*dw(:, i-1);
         Centr(:,i) = cross(w(:, i), cross(w(:, i), r(: ,i)));
         Coriol(:,i) =2*cross(dq(i)*w(:, i), R(:, al:bh)'*z0);
         ddp(:, i) = R(:, al:bh)'*(ddp(:, i-1) + ddq(i)*z0) + Coriol(:, i)+ cross(dw(:, i), r(:, i)) + Centr(:,i);
         ddpc(:, i) = ddp(:, i) + cross(dw(:, i), -rc(:, i)) +cross(w(:, i), cross(w(:, i), -rc(: ,i)));
    else 
         dw(:, i) = R(:, al:bh)'*(dw(:, i-1) + ddq(i)*z0 + cross(dq(i)*w(:, i-1), z0));
         Centr(:,i) = cross(w(:, i), cross(w(:, i), r(: ,i)));
         ddp(:, i) = R(:, al:bh)'*ddp(:, i-1) +  cross(dw(:, i), r(:, i)) + Centr(:,i);
         ddpc(:, i) = ddp(:, i) + cross(dw(:, i), -rc(:, i)) +cross(w(:, i), cross(w(:, i), -rc(: ,i)));
    end
    
    G(:, i) = R(:, al:bh)'*G(:, i-1); 
    
    
end



%%


w = simplify(w);
dw = simplify(dw);
dp = simplify(dp);
dpc = simplify(dpc);
ddp = simplify(ddp);
ddpc = simplify(ddpc);
G = simplify(G);


%% Backward
f = zeros(3, 4, 'sym');
mu = zeros(3, 4, 'sym');
tau = zeros(1, 4, 'sym');
m = [0,  m1, m2, m3];
R = [R_10, R_21, R_32];

I1 = [(1/12)*m1*l1^2 0 0;
       0 (1/12)*m1*l1^2 0;
       0 0 0];

I2 = [0 0 0;
      0 (1/12)*m2*l2^2 0;
       
      0 0 (1/12)*m2*l2^2];
  
I3 = [(1/12)*m3*q3^2 0 0;
       0 (1/12)*m3*q3^2 0;
       0 0 0];
I= [I1, I2, I3];
%%
for i = 4:-1:2
    if i == 4
        al = 1+3*(i-2);
        bh = 3*(i-1);
        f(:, i) = m(i)*ddpc(:,i) - m(i)*G(:, i);
        mu(:, i) = -cross(f(:, i), rc(:, i)) + I(:, al:bh)*dw(:, i) + cross(w(:,i), I(:, al:bh)*w(:,i));
        tau(:,i) = f(:, i)'*(R(:, al:bh)'*z0);
    else 
        
        al = 1+3*(i-1);
        bh = 3*(i);
        ali = 1+3*(i-2);
        bhi = 3*(i-1);
        f(:, i) = R(:, al:bh)*f(:, i+1) + m(i)*ddpc(:,i) - m(i)*G(:, i);
        mu(:, i) = -cross(f(:, i), rc(:, i)) + R(:, al:bh)*mu(:, i+1) + cross(R(:, al:bh)*f(:, i+1),-rc(:, i)) + I(:, ali:bhi)*dw(:, i) + cross(w(:,i), I(:, ali:bhi)*w(:,i));
        tau(:,i) = mu(:, i)'*(R(:, ali:bhi)'*z0);
    end
    
    
end
f = simplify(f);
mu = simplify(mu);
tau = simplify(tau);
%% Generalized forces
%% plots
% Assignment of variables
time = 0: 0.01 :10;
%% Link 1
%Velocity plot Link 1
Plotme(t, time, dp(:, 2), 'v_{ex}', 'v_{ey}', 'v_{ez}', 'time', 'Velocity link1', 'Velocity of first link', '101')
Plotme(t, time, dpc(:, 2), 'v_{cx}', 'v_{cy}', 'v_{cz}', 'time', 'Velocity link1', 'Velocity of first link', '102')
Plotme(t, time, w(:, 2), '\omega_x', '\omega_y', '\omega_z', 'time', 'Angular velocity of link1', 'Angular velocity of first link', '103')

%Acceleration Link 1
Plotme(t, time, ddp(:, 2), 'a_{ex}', 'a_{ey}', 'a_{ez}', 'time', 'Acceleration link1', 'Acceleration  of first link', '104')
Plotme(t, time, ddpc(:, 2), 'a_{cx}', 'a_{cy}', 'a_{cz}', 'time', 'Acceleration link1', 'Acceleration of first link', '105')
Plotme(t, time, dw(:, 2), '\alpha_x', '\alpha_y', '\alpha_z', 'time', 'Angular acceleration of link1', 'Angular acceleration of first link', '106')

%Coriolis, Centrifugal, Gravity forces of link 1
Plotme(t, time, Centr(:, 2), 'acen_{x}', 'acen_{y}', 'acen_{z}', 'time', 'Centrifugal force of link 1', 'Centrifugal force of link 1', '107')
Plotme(t, time, G(:, 2), 'g_{x}', 'g_{y}', 'g_{z}', 'time', 'Gravity term link2', 'Gravity term link1', '108')

% Link 1 Force and torque
Plotme(t, time, f(:, 2), 'f_{x}', 'f_{y}', 'f_{z}', 'time', 'Force link1', 'Force of first link', '109')
Plotme(t, time, mu(:, 2), '\mu_{x}', '\mu_{y}', '\mu_{z}', 'time', 'moment link1', 'Moment of first link', '110')

%% Link 2
%Velocity plot Link 2
Plotme(t, time, dp(:, 3), 'v_{ex}', 'v_{ey}', 'v_{ez}', 'time', 'Velocity link2', 'Velocity of second link', '201')
Plotme(t, time, dpc(:, 3), 'v_{cx}', 'v_{cy}', 'v_{cz}', 'time', 'Velocity link2', 'Velocity of second link', '202')
Plotme(t, time, w(:, 3), '\omega_x', '\omega_y', '\omega_z', 'time', 'Angular velocity of link2', 'Angular velocity of second link', '203')

%Acceleration Link 2
Plotme(t, time, ddp(:, 3), 'a_{ex}', 'a_{ey}', 'a_{ez}', 'time', 'Acceleration link2', 'Acceleration  of second link', '204')
Plotme(t, time, ddpc(:, 3), 'a_{cx}', 'a_{cy}', 'a_{cz}', 'time', 'Acceleration link2', 'Acceleration of second link', '205')
Plotme(t, time, dw(:, 3), '\alpha_x', '\alpha_y', '\alpha_z', 'time', 'Angular acceleration of link2', 'Angular acceleration of second link', '206')

%Coriolis, Centrifugal, Centrifugal forces of link 2
Plotme(t, time, Centr(:, 3), 'acen_{x}', 'acen_{y}', 'acen_{z}', 'time', 'Centrifugal force of link 2', 'Centrifugal force of link 2', '207')
Plotme(t, time, G(:, 3), 'g_{x}', 'g_{y}', 'g_{z}', 'time', 'Gravity term link2', 'Gravity term link2', '208')

% Link 2 Force and torque
Plotme(t, time, f(:, 3), 'f_{x}', 'f_{y}', 'f_{z}', 'time', 'Force link2', 'Force of second link', '209')
Plotme(t, time, mu(:, 3), '\mu_{x}', '\mu_{y}', '\mu_{z}', 'time', 'moment link2', 'Moment of second link', '210')

%% Link 3
%Velocity plot
Plotme(t, time, dp(:, 4), 'v_{ex}', 'v_{ey}', 'v_{ez}', 'time', 'Velocity link3', 'Velocity of end effector', '301')
Plotme(t, time, dpc(:, 4), 'v_{cx}', 'v_{cy}', 'v_{cz}', 'time', 'Velocity link3', 'Velocity of centre of mass effector', '302')
Plotme(t, time, w(:, 4), '\omega_x', '\omega_y', '\omega_z', 'time', 'Angular velocity of link3', 'Angular velocity of end effector', '303')

%Acceleration
Plotme(t, time, ddp(:, 4), 'a_{ex}', 'a_{ey}', 'a_{ez}', 'time', 'Acceleration link3', 'Acceleration of end effector', '304')
Plotme(t, time, ddpc(:, 4), 'a_{cx}', 'a_{cy}', 'a_{cz}', 'time', 'Acceleration link3', 'Acceleration of centre of mass effector', '305')
Plotme(t, time, dw(:, 4), '\alpha_x', '\alpha_y', '\alpha_z', 'time', 'Angular acceleration of link3', 'Angular acceleration of end effector', '306')

%Coriolis, Centrifugal, Centrifugal forces of link 3
Plotme(t, time, Coriol(:, 4), 'acor_{x}', 'acor_{y}', 'acor_{z}', 'time', 'Coriolis acceleration link3', 'Coriolis acceleration link3', '307')
Plotme(t, time, Centr(:, 4), 'acen_{x}', 'acen_{y}', 'acen_{z}', 'time', 'Centrifugal force of link 3', 'Centrifugal force of link 3', '308')
Plotme(t, time, G(:, 4), 'g_{x}', 'g_{y}', 'g_{z}', 'time', 'Gravity term link3', 'Gravity term link3', '309')

% Link 3 Force and torque
Plotme(t, time, f(:, 4), 'f_{x}', 'f_{y}', 'f_{z}', 'time', 'Force link3', 'Force of end effector', '310')
Plotme(t, time, mu(:, 4), '\mu_{x}', '\mu_{y}', '\mu_{z}', 'time', 'moment link3', 'Moment of end effector ', '311')


%% Generalized force
taut = [tau(:,2); tau(:, 3); tau(:, 4)];
Plotme(t, time, taut, '\tau_1', '\tau_2', '\tau_3', 'time', 'torques', 'Generalized force', 'Torque')



