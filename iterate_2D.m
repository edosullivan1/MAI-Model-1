% Script to compute internall forces in joints using OpenPose 2D position
% co-ordinates.

pos = csvread('ciaran_2D_position.csv',4,1);

t = pos(:, 1);

Fpx = pos(:, 5); Fpy = pos(:, 6);    %Neck    2
Gpx = pos(:, 53); Gpy = pos(:, 54);     %Head   18

Hpx = pos(:, 68); Hpy = pos(:, 69);    %Right Toe   23
Ipx = pos(:, 35); Ipy = pos(:, 36);    %Right Ankle 12
Jpx = pos(:, 32); Jpy = pos(:, 33);    %Right Knee  11
Kpx = pos(:, 29); Kpy = pos(:, 30);    %Right Hip   10
Lpx = pos(:, 8); Lpy = pos(:, 9);    %Right Shoulder  3
Mpx = pos(:, 26); Mpy = pos(:, 27);    %Pelvis  

Opx = pos(:, 14); Opy = pos(:, 15);    %Right Wrist

length_scale = sqrt( (Jpx(30)-Ipx(30))^2 + (Jpy(30)-Ipy(30))^2 );
m_per_pixel = 460/length_scale; %assumes were 0.46m apart (length of lower leg)

Fx = Fpx; Fy = -Fpy;    %Neck    2
Gx = Gpx; Gy = -Gpy;     %Head   18

Hx = Hpx; Hy = -Hpy;    %Right Toe   23
Ix = Ipx; Iy = -Ipy;    %Right Ankle 12
Jx = Jpx; Jy = -Jpy;    %Right Knee  11
Kx = Kpx; Ky = -Kpy;    %Right Hip   10
Lx = Lpx; Ly = -Lpy;    %Right Shoulder  3
Mx = Mpx; My = -Mpy;    %Pelvis  

Ox = Opx; Oy = -Opy;    %Right Wrist

m_total = 80;
L_total = 1.86;

density = 1000;
g = 9.81;
g_z = -9.81;

head_scale = 1.1;

% Estimates for trunk are from Plagenhoef et al., 1983 
% https://exrx.net/Kinesiology/Segments

m_foot_L  = 0.5*0.0143*m_total;
m_leg_L   = 0.5*0.0475*m_total;
m_thigh_L = 0.5*0.2416*m_total;    %And pelvis
m_arm_L = 0.5*0.057*m_total;

m_torso = 0.551*m_total;     %not including arms

m_foot_R  = 0.5*0.0143*m_total;
m_leg_R   = 0.5*0.0475*m_total;
m_thigh_R = 0.5*0.2416*m_total;    %And pelvis
m_arm_R = 0.5*0.057*m_total;

m_head_neck  = 0.0826*m_total;

m_check = m_foot_L + m_leg_L + m_thigh_L + m_arm_L + m_foot_R + m_leg_R + m_thigh_R + m_torso + m_arm_R + m_head_neck;

% mass of body above hip
m_above_hip = (m_arm_L + m_arm_R + m_torso +  m_head_neck);

% mass of body above knee
m_above_knee_L = m_arm_L + m_thigh_L + 0.5 * (m_torso +  m_head_neck);

% mass of body above knee
m_above_knee_R = m_arm_R + m_thigh_R + 0.5 * (m_torso +  m_head_neck);

% mass of body above ankle
m_above_ankle_L = m_arm_L + m_leg_L + m_thigh_L + 0.5 * (m_torso +  m_head_neck);

% mass of body above ankle
m_above_ankle_R = m_arm_R + m_leg_R + m_thigh_R + 0.5 * (m_torso +  m_head_neck);

%vertical force through hip
F_matrix_hip_R_y = 0.5 * m_above_hip * g;

%vertical force through knee
F_matrix_knee_R_y = m_above_knee_R * g;

%vertical force through ankle
F_matrix_ankle_R_y = m_above_ankle_R * g;

num_rows = height(Hx);


for i = 1:num_rows


L_foot_R(i) = sqrt( (Hx(i) - Ix(i))^2 + (Hy(i) - Iy(i))^2);
L_leg_R(i) = sqrt( (Jx(i) - Ix(i))^2 + (Jy(i) - Iy(i))^2);
L_thigh_R(i) = sqrt( (Kx(i) - Jx(i))^2 + (Ky(i) - Jy(i))^2);
L_torso_R(i) = sqrt( (Fx(i) - Mx(i))^2 + (Fy(i) - My(i))^2);
L_arm_R(i) = sqrt( (Ox(i) - Lx(i))^2 + (Oy(i) - Ly(i))^2);

L_head_neck(i)  = head_scale*sqrt( (Gx(i) - Fx(i))^2 + (Gy(i) - Fy(i))^2);

L_check(i) = L_leg_R(i) + L_thigh_R(i) + L_torso_R(i) + L_head_neck(i);

theta_torso(i) = atand((Fy(i) - My(i))/sqrt(((Fx(i) - Mx(i))^2)));

theta_foot_R(i)  = atand((Iy(i) - Hy(i))/sqrt(((Ix(i) - Hx(i))^2)));
theta_leg_R(i) = atand((Jy(i) - Iy(i))/sqrt(((Jx(i) - Ix(i))^2)));
theta_thigh_R(i) = atand((Ky(i) - Jy(i))/sqrt(((Kx(i) - Jx(i))^2)));
theta_arm_R(i) = atand((Ly(i) - Oy(i))/sqrt(((Lx(i) - Ox(i))^2)));

theta_head(i)  = atand((Gy(i) - Fy(i))/sqrt(((Gx(i) - Fx(i))^2)));

torso_cg_x(i) =  (Fx(i) + Mx(i))/2;
torso_cg_y(i) =  (Fy(i) + My(i))/2;

head_neck_cg_x(i) =  (Fx(i) + Gx(i))/2;
head_neck_cg_y(i) =  (Fy(i) + Gy(i))/2;

foot_R_cg_x(i) = (Ix(i) + Hx(i))/2;
foot_R_cg_y(i) = (Iy(i) + Hy(i))/2;

leg_R_cg_x(i) =  (Jx(i) + Ix(i))/2;
leg_R_cg_y(i) =  (Jy(i) + Iy(i))/2;

thigh_R_cg_x(i) =  (Kx(i) + Jx(i))/2;
thigh_R_cg_y(i) =  (Ky(i) + Jy(i))/2;

arm_R_cg_x(i) =  (Ox(i) + Lx(i))/2;
arm_R_cg_y(i) =  (Oy(i) + Ly(i))/2;


whole_body_cg_x(i) =  (m_foot_R*foot_R_cg_x(i) + m_leg_R*leg_R_cg_x(i) ...
    + m_thigh_R*thigh_R_cg_x(i) + m_foot_R*foot_R_cg_x(i) ...
    + m_leg_R*leg_R_cg_x(i) + m_thigh_R*thigh_R_cg_x(i) ...
    + m_arm_R*arm_R_cg_x(i) + m_arm_R*arm_R_cg_x(i) ...
    + m_torso*torso_cg_x(i) +  m_head_neck*head_neck_cg_x(i))/m_total;

whole_body_cg_y(i) =  (m_foot_R*foot_R_cg_y(i) + m_leg_R*leg_R_cg_y(i) + m_thigh_R*thigh_R_cg_y(i) ...
    + m_foot_R*foot_R_cg_y(i) + m_leg_R*leg_R_cg_y(i) + m_thigh_R*thigh_R_cg_y(i) ...
    + m_arm_R*arm_R_cg_y(i) + m_arm_R*arm_R_cg_y(i) ...
    + m_torso*torso_cg_y(i) + m_head_neck*head_neck_cg_y(i))/m_total;

%Right Hip
% cg of body above hip in pixels
cg_above_hip_x(i) =  (m_arm_R*arm_R_cg_x(i) + m_arm_R*arm_R_cg_x(i) ...
        + m_torso*torso_cg_x(i) ...
    +  m_head_neck*head_neck_cg_x(i))/m_above_hip;

cg_above_hip_y(i) =  (m_arm_R*arm_R_cg_y(i) + m_arm_R*arm_R_cg_y(i) ...
+ m_torso*torso_cg_y(i) ...
    +  m_head_neck*head_neck_cg_y(i))/m_above_hip;

%d is horizontal moment arm (meters) about the hip of the weight above hip 
d_hip_R_x(i) = m_per_pixel * (cg_above_hip_x(i) - Kx(i));
d_hip_R_y(i) = m_per_pixel * (cg_above_hip_y(i) - Ky(i));

%moment about hip (Nm)
M_hip_R(i) = d_hip_R_x(i) .* F_matrix_hip_R_y;

%M_hip_mag_R(i) = M_hip_R(i);



% cg of body above knee in pixels
cg_above_knee_R_x(i) =  (m_thigh_R*thigh_R_cg_x(i) + 0.5 * m_torso*torso_cg_x(i) ...
    + m_arm_R*arm_R_cg_x(i)  ...
    +  0.5 * m_head_neck*head_neck_cg_x(i))/m_above_knee_R;

cg_above_knee_R_y(i) =  (m_thigh_R*thigh_R_cg_y(i) + 0.5 * m_torso*torso_cg_y(i) ...
    + m_arm_R*arm_R_cg_y(i) ...
    +  0.5 * m_head_neck*head_neck_cg_y(i))/m_above_knee_R;

%d is horizontal moment arm (meters) about the knee
d_knee_R_x(i) = m_per_pixel * (cg_above_knee_R_x(i) - Jx(i));
d_knee_R_y(i) = m_per_pixel * (cg_above_knee_R_y(i) - Jy(i));

%d_matrix_knee_R(i,:) = [d_knee_R_x(i), d_knee_R_y(i), 0];
d_knee_R_mag = transpose(d_knee_R_x);


%moment about knee (Nm)
M_knee_R(i) = -d_knee_R_x(i) .* F_matrix_knee_R_y; 

%M_knee_mag_R(i) = sqrt((M_knee_R(i, 1).^2) + (M_knee_R(i, 2).^2));


%Right Ankle

% cg of body above ankle in pixels
cg_above_ankle_R_x(i) =  (m_leg_R*leg_R_cg_x(i) + m_thigh_R*thigh_R_cg_x(i) + 0.5 * m_torso*torso_cg_x(i) ...
    + m_arm_R*arm_R_cg_x(i) ...
    +  0.5 * m_head_neck*head_neck_cg_x(i))/m_above_ankle_R;

cg_above_ankle_R_y(i) =  (m_leg_R*leg_R_cg_y(i) +m_thigh_R*thigh_R_cg_y(i) + 0.5 * m_torso*torso_cg_y(i) ...
    + m_arm_R*arm_R_cg_y(i) ...
    +  0.5 * m_head_neck*head_neck_cg_y(i))/m_above_ankle_R;

%d is horizontal moment arm (meters) about the ankle
d_ankle_R_x(i) = m_per_pixel * (cg_above_ankle_R_x(i) - Ix(i));
d_ankle_R_y(i) = m_per_pixel * (cg_above_ankle_R_y(i) - Iy(i));

%d_matrix_ankle_R(i,:) = [d_ankle_R_x(i), d_ankle_R_y(i), 0];


%moment about ankle (Nm)
M_ankle_R(i) = d_ankle_R_x(i) .* F_matrix_ankle_R_y;

%M_ankle_mag_R(i) = sqrt((M_ankle_R(i, 1).^2) + (M_ankle_R(i, 2).^2));



end
% uncomment to write knee moment arm to csv file
%writematrix(d_knee_R_mag,'d_knee_2D.csv') 

M_hip_R_transpose = transpose(M_hip_R);
M_knee_R_transpose = transpose(M_knee_R);
M_ankle_R_transpose = transpose(M_ankle_R);

figure (1)
plot(t,M_knee_R)
title('Torque-Time plot of Right Knee')
xlabel('Time (s)')
ylabel('Torque (Nmm)')


figure (2)
plot(t,M_hip_R)
title('Torque-Time plot of Right Hip')
xlabel('Time (s)')
ylabel('Torque (Nmm)')


figure (3)
plot(t,M_ankle_R)
title('Torque-Time plot of Right Ankle')
xlabel('Time (s)')
ylabel('Torque (Nmm)')






Check_x = [ Hx(40), Ix(40), Jx(40), Kx(40), Lx(40), Ox(40), Fx(40), Gx(40)];

Check_y = [ Hy(40), Iy(40), Jy(40), Ky(40), Ly(40), Oy(40), Fy(40), Gy(40)];


Check_px = [ Hpx(40), Ipx(40), Jpx(40), Kpx(40), Lpx(40), Opx(40), Fpx(40), Gpx(40)];

Check_py = [ Hpy(40), Ipy(40), Jpy(40), Kpy(40), Lpy(40), Opy(40), Fpy(40), Gpy(40)];


Above_CG_x = [ cg_above_hip_x; cg_above_knee_R_x;   cg_above_ankle_R_x];
Above_CG_y = [ cg_above_hip_y; cg_above_knee_R_y;   cg_above_ankle_R_y];




%figure(5)
%plot(Check_px,Check_py,'-o', 'Color', 'b', 'MarkerSize',10,'MarkerFaceColor','#D9FFFF');
%axis equal

cg_check_x = [foot_R_cg_x(40), leg_R_cg_x(40), thigh_R_cg_x(40), torso_cg_x(40), head_neck_cg_x(40), arm_R_cg_x(40)];

cg_check_y = [foot_R_cg_y(40), leg_R_cg_y(40), thigh_R_cg_y(40), torso_cg_y(40), head_neck_cg_y(40), arm_R_cg_y(40)];

figure(4)
plot(Check_x,Check_y,'-o', 'Color', 'b', 'MarkerSize',10,'MarkerFaceColor','#D9FFFF');
axis equal
hold on
plot(cg_check_x, cg_check_y, '-o', 'Color', 'g', 'MarkerSize',10,'MarkerFaceColor','#D9FFFF');
hold on
plot(Above_CG_x(:,1), Above_CG_y(:,1), 'o', 'Color', 'r', 'MarkerSize',10,'MarkerFaceColor','#D9FFFF')