close all; clear

% pica = imread('IMG_4017a.jpg');
% 
% figure(1); subplot(1,2,1)
% imagesc(pica); hold on; 
% 
% picb = imread('IMG_4017b.jpg');
% 
% figure(1); subplot(1,2,2)
% imagesc(picb); hold on; 
% 

%[x,y] = ginput(7) click to get 7 key points of body joints... pose estimation gives this
%[x1,y1] = ginput(2) click to get 2 key points for scale... 
%{
001 Keypoint - frame 1
[x] = [592.376, 435.351, 466.749, 445.878, 419.57, 440.551, 503.444];
[y] = [1671.28, 1613.86, 1220.93, 890.989, 440.484, 445.762, 309.59];

[x1] = [435.351, 466.749];
[y1] = [1613.86, 1220.93];

%016 Keypoint - frame 16
[x] = [592.251, 440.35, 477.057, 440.494, 435.251, 456.123, 513.666];
[y] = [1676.52, 1613.79, 1226.18, 885.734, 440.623, 445.896, 314.745];

[x1] = [440.35, 477.057];
[y1] = [1613.79, 1226.18];
%}
%{
%031 Keypoint - frame 31
[x] = [597.576, 440.637, 576.612, 325.372, 466.717, 487.683, 555.659];
[y] = [1671.33, 1619.05, 1257.58, 1006.29, 602.764, 618.486, 487.69];

[x1] = [440.637, 576.612];
[y1] = [1619.05, 1257.58];

%046 Keypoint - frame 46
[x] = [592.418, 445.863, 608.129, 283.253, 518.991, 513.8, 587.253];
[y] = [1671.31, 1618.88, 1273.32, 1137.22, 775.824, 780.903, 686.455];

[x1] = [445.863, 608.129];
[y1] = [1618.88, 1273.32];

%061 Keypoint - frame 61
[x] = [592.422, 451.085, 629.036, 267.729, 534.656, 524.405, 592.277];
[y] = [1676.54, 1618.86, 1289.05, 1273.15, 927.497, 927.632, 828.262];

[x1] = [451.085, 629.036];
[y1] = [1618.86, 1289.05];


%076 Keypoint - frame 76
[x] = [592.446, 450.963, 634.358, 272.847, 524.274, 508.544, 571.31];
[y] = [1676.46, 1619.03, 1299.38, 1351.78, 985.192, 990.32, 890.891];

[x1] = [450.963, 634.358];
[y1] = [1619.03, 1299.38];
%}
%{
%091 Keypoint - frame 91
[x] = [592.361, 450.972, 608.3, 272.977, 513.652, 503.262, 566.058];
[y] = [1671.45, 1619.03, 1278.68, 1205.17, 849.158, 854.096, 749.449];

[x1] = [450.972, 608.3];
[y1] = [1619.03, 1278.68];
%}
%106 Keypoint - frame 106
[x] = [597.639, 440.626, 519.01, 346.286, 456.123, 482.389, 519.128];
[y] = [1676.55, 1613.93, 1236.72, 938.141, 508.521, 524.117, 382.979];

[x1] = [440.626, 519.01];
[y1] = [1613.93, 1236.72];
%}


length_scale = sqrt( (x1(2)-x1(1))^2 + (y1(2)-y1(1))^2 );

m_per_pixel = 0.46/length_scale; %assumes x1 y1 were 0.28m apart (head height)

%m_total = input('total mass (kg) ')
%L_total = input('total height (m) ')
m_total = 120;
L_total = 1.86;

density = 1000;
g = 9.81;

%lengths in m
L_foot  = m_per_pixel*sqrt( (x(2)-x(1))^2 + (y(2)-y(1))^2 );
L_leg   = m_per_pixel*sqrt( (x(3)-x(2))^2 + (y(3)-y(2))^2 );
L_thigh = m_per_pixel*sqrt( (x(4)-x(3))^2 + (y(4)-y(3))^2 );
L_torso = m_per_pixel*sqrt( (x(5)-x(4))^2 + (y(5)-y(4))^2 );
%L_neck  = m_per_pixel*sqrt( (x(6)-x(5))^2 + (y(6)-y(5))^2 );
L_head_neck  = 1.5*m_per_pixel*sqrt( (x(7)-x(6))^2 + (y(7)-y(6))^2 );

L_check = L_leg + L_thigh + L_torso + L_head_neck;

area_factor = m_total/(density*L_total); 

m_foot  = L_foot*density*area_factor;
m_leg   = L_leg*density*area_factor;
m_thigh = L_thigh*density*area_factor;
m_torso = L_torso*density*area_factor;
%m_neck  = L_neck*density*area_factor;
m_head_neck  = L_head_neck*density*area_factor;

m_check = m_foot + m_leg + m_thigh + m_torso + m_head_neck

Ax = x(1); Ay = y(1);
Bx = x(2); By = y(2);
Cx = x(3); Cy = y(3);
Dx = x(4); Dy = y(4);
Ex = x(5); Ey = y(5);
Fx = x(6); Fy = y(6);
Gx = x(7); Gy = y(7);

theta_foot  = atand((By - Ay)/(Bx - Ax));
theta_leg   = atand((Cy - By)/(Cx - Bx));
theta_thigh = atand((Dy - Cy)/(Dx - Cx));
theta_torso = atand((Ey - Dy)/(Ex - Dx));
%theta_neck  = atand((Fy - Ey)/(Fx - Ex));
theta_head  = atand((Gy - Fy)/(Gx - Fx));

foot_cg_x = (Ax + Bx)/2;
foot_cg_y = (Ay + By)/2;

leg_cg_x =  (Bx + Cx)/2;
leg_cg_y =  (By + Cy)/2;

thigh_cg_x =  (Cx + Dx)/2;
thigh_cg_y =  (Cy + Dy)/2;

torso_cg_x =  (Dx + Ex)/2;
torso_cg_y =  (Dy + Ey)/2;

%neck_cg_x =  (Ex + Fx)/2;
%neck_cg_y =  (Ey + Fy)/2;

head_neck_cg_x =  (Fx + Gx)/2;
head_neck_cg_y =  (Fy + Gy)/2;

%get whole body cg
whole_body_cg_x =  (m_foot*foot_cg_x + m_leg*leg_cg_x + m_thigh*thigh_cg_x ...
    + m_torso*torso_cg_x +  m_head_neck*head_neck_cg_x)/m_total;

whole_body_cg_y =  (m_foot*foot_cg_y + m_leg*leg_cg_y + m_thigh*thigh_cg_y ...
    + m_torso*torso_cg_y + m_head_neck*head_neck_cg_y)/m_total;

%mass of body above hip
m_above_hip = m_torso +  m_head_neck;

%cg of body above hip in pixels
cg_above_hip_x =  (m_torso*torso_cg_x ...
    +  m_head_neck*head_neck_cg_x)/m_above_hip;

cg_above_hip_y =  (m_torso*torso_cg_y ...
    +  m_head_neck*head_neck_cg_y)/m_above_hip;

%d is horizontal moment arm (meters) about the hip of the weight above hip 
d_hip = m_per_pixel*(cg_above_hip_x - Dx);

%moment about hip (Nm)
M_hip = m_above_hip*g*d_hip

%vertical force through hip
F_hip_y = (m_torso + m_head_neck)*g


%mass of body above knee 
m_above_knee = m_thigh + m_torso +  m_head_neck;

%cg of body above knee in pixels
cg_above_knee_x =  (m_thigh*thigh_cg_x + m_torso*torso_cg_x ...
    +  m_head_neck*head_neck_cg_x)/m_above_knee;

cg_above_knee_y =  (m_thigh*thigh_cg_y + m_torso*torso_cg_y ...
    +  m_head_neck*head_neck_cg_y)/m_above_knee;

%d is horizontal moment arm (meters) about the knee of the weight above knee 
d_knee = m_per_pixel*(cg_above_knee_x - Cx);

%moment about knee (Nm)
M_knee = m_above_knee*g*d_knee

%vertical force through knee
F_knee_y = (m_thigh + m_torso + m_head_neck)*g


%mass of body above ankle
m_above_ankle = m_leg + m_thigh + m_torso +  m_head_neck;

%cg of body above ankle in pixels
cg_above_ankle_x =  (m_leg*leg_cg_x + m_thigh*thigh_cg_x + m_torso*torso_cg_x ...
    +  m_head_neck*head_neck_cg_x)/m_above_ankle;

cg_above_knee_y =  (m_leg*leg_cg_y + m_thigh*thigh_cg_y + m_torso*torso_cg_y ...
    +  m_head_neck*head_neck_cg_y)/m_above_ankle;

%d is horizontal moment arm (meters) about the ankle of the weight above ankle 
d_ankle = m_per_pixel*(cg_above_ankle_x - Bx);

%moment about ankle (Nm)
M_ankle = m_above_ankle*g*d_ankle

%vertical force through ankle
F_ankle_y = (m_leg + m_thigh + m_torso + m_head_neck)*g


%plots joints and segments mass centres
h1 = plot([Ax Bx Cx Dx Ex Fx Gx],[Ay By Cy Dy Ey Fy Gy],'ro');
set(h1,'markersize',10,'markeredgecolor','r','markerfacecolor','r')

h2 = plot([Ax Bx Cx Dx Ex Fx Gx],[Ay By Cy Dy Ey Fy Gy],'g-');
set(h2,'linewidth',2)

h3 = plot([foot_cg_x leg_cg_x thigh_cg_x torso_cg_x head_neck_cg_x],[foot_cg_y leg_cg_y thigh_cg_y torso_cg_y head_neck_cg_y],'ks');
set(h3,'markersize',10,'markeredgecolor','g','markerfacecolor','g')

h4 = plot([whole_body_cg_x],[whole_body_cg_x],'m*');
set(h3,'markersize',10,'markeredgecolor','m','markerfacecolor','m')





