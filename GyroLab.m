syms t real

%alpha__, beta_, gamma_ and delta_ time functions
syms alpha_(t) beta_(t) gamma_(t) delta_(t)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up measurements (symbolic in part 1)
% Constants
syms g                                                      % gravity [ms^-2]
% Frame Measurements
syms L                                                      % Distance from O to G [m]
syms h_rod                                                  % Height of vertical rod [m]
syms R_min_H_tor                                            % Horizontal torus minor radius [m]
syms R_maj_H_tor                                            % Horizontal torus marjor radius[m]
syms R_min_V_tor                                            % Vertical torus minor radius [m]
syms R_maj_V_tor                                            % Vertical torus major radius [m]

% Rotor measurements
syms m_frame                                                % Mass of frame [kg] (Measured)
syms m_rotor                                                % Mass of rotor [kg] (Measured)
syms H_rot                                                  % Rotor rod height [m]
syms r_rot                                                  % Rotor rod radius [m] (Vertical rod)
syms R_min_rotor                                            % Minor radius of rotor torus [m]
syms R_maj_rotor                                            % Major radius of rotor torus [m]

% Top and bottom spheres (For animation only, DO NOT MODEL THEM AS A PART OF THE FRAME)
syms R_sph                                                  %radius of the top and bottom spheres [m]

%Calculated constants
%%%%%% Frame %%%%%

%Total frame density [kg/m^3] (Test).
% Variable Name: rho_T
% g = 9.81;
% m_frame = 0.023;                                            % Mass of frame [kg] (Measured)
% m_rotor = 0.045;                                            % Mass of rotor [kg] (Measured)
% H_rot = 60*10^-3;
% r_rot = 1.555*10^-3;
% R_min_rotor = 3.305*10^-3;
% R_maj_rotor = 28*10^-3;
% h_rod = 88*10^-3;
% R_min_H_tor = 1.155*10^-3;;
% R_maj_H_tor = 36*10^-3;
% R_sph = 2.365*10^-3;
% R_min_V_tor = 1.715*10^-3
% R_maj_V_tor = 31.5*10^-3;
% L = h_rod/2;           

Vol_V_tor = (2*pi*R_maj_V_tor)*(pi*R_min_V_tor^2);                   % Correct      
Vol_H_tor =(2*pi*R_maj_H_tor)*(pi*R_min_H_tor^2);                    % Correct      
Vol_rod_ends = (h_rod-H_rot)*pi*r_rot^2;                             % Correct      

rho_T = m_frame/(Vol_V_tor+Vol_H_tor+Vol_rod_ends);                  % Correct              

                          
%%%%%% Rotor %%%%%

% Desity of the rotor [kg/m^3] (Pretest).
% Variable Name: rho_rotor

Vol_central_rod = H_rot*pi*r_rot^2;                                                              % Correct      
Vol_rot_tor = (2*pi*R_maj_rotor)*(pi*R_min_rotor^2); % Assume mass is only on the torus          % Correct      

rho_rotor = m_rotor/(Vol_rot_tor+Vol_central_rod);                                               % Correct      

% Central rod considered as a cylinder belonging to the frame

% Recalculated mass of frame [kg] (Test).
% Variable Name: m_frame 
m_frame = m_frame + (Vol_central_rod * rho_rotor);                              % Correct 

% Central rod considered as a cylinder belonging to the frame

% Recalculated mass of frame [kg] (Test).
% Variable Name: m_frame 
             

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inertia tensor
% Inerta tensor of the 3 frame components in frame 3 

% Inertia tensor of horizontal torus about its center of mass expressed in frame 3 (Pretest). 
% Variable Name: IGtorH_3
m_H_tor = rho_T*Vol_H_tor;                                                             % Correct      

IGtorH_3 = [m_H_tor/8*(5*R_min_H_tor^2 + 4*R_maj_H_tor^2) 0 0; 
            0 m_H_tor/8*(5*R_min_H_tor^2 + 4*R_maj_H_tor^2) 0;                         % Correct      
            0 0 m_H_tor/4*(4*R_maj_H_tor^2 + 3*R_min_H_tor^2)]

% Inertia tensor of vertical torus about its center of mass expressed in frame 3 (Test)
% Variable Name: IGtorV_3
m_V_tor = rho_T*Vol_V_tor;                     % Correct 

IGtorV_3 = [m_V_tor/8*(5*R_min_V_tor^2 + 4*R_maj_V_tor^2) 0 0;                      % Correct 
            0 m_V_tor/4*(4*R_maj_V_tor^2 + 3*R_min_V_tor^2) 0;
            0 0 m_V_tor/8*(5*R_min_V_tor^2 + 4*R_maj_V_tor^2)]

% Inertia tensor of the rod about its center of mass expressed in frame 3 
% Variable Name: IGrod_3
m_rod = rho_T*Vol_rod_ends + rho_rotor*Vol_central_rod;                     % Correct 

IGrod_3 = [m_rod/12*h_rod^2 + 1/4*m_rod*r_rot^2 0 0;                     % Correct 
           0 m_rod/12*h_rod^2 + 1/4*m_rod*r_rot^2 0;
           0 0 1/2*m_rod*r_rot^2]

% Inertia tensor of the frame about its center of mass expressed in frame 3 (Test)
% Variable Name: IGframe_3
IGframe_3 = IGrod_3 + IGtorV_3 + IGtorH_3                     % Correct 

% Inertia tensor of rotor about its center of mass in frame 4 
% Variable Name: IGrotor_4
m_rot = rho_rotor*pi*(R_maj_rotor)*2*pi*R_min_rotor^2;                         % Correct

IGrotor_4 = [m_rot/8*(5*R_min_rotor^2 + 4*R_maj_rotor^2) 0 0;                      % Correct 
            0 m_rot/8*(5*R_min_rotor^2 + 4*R_maj_rotor^2) 0;
            0 0 m_rot/4*(4*R_maj_rotor^2 + 3*R_min_rotor^2)]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up rotational matrices

%Rotation matrix from frame 4 to frame 0 (Pretest).
%Variable Name: R04
R01 = [cos(alpha_) -sin(alpha_) 0;
       sin(alpha_) cos(alpha_) 0;               % Correct      
       0 0 1];
R12 = [1 0 0;
       0 cos(beta_) -sin(beta_);
       0 sin(beta_) cos(beta_)];
R23 = [cos(gamma_) -sin(gamma_) 0;              % Correct      
       sin(gamma_) cos(gamma_) 0;
       0 0 1];
R34 = [cos(delta_) -sin(delta_) 0;
       sin(delta_) cos(delta_) 0;               % Correct      
       0 0 1];
R04 = R01*R12*R23*R34                           % Correct      


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Kinematics
% Centers of mass positions from origin
% Center of mass position of the frame. Variable Name: rOG_3
rOG_3 = [0;0;L]                             % Correct

% Center of mass position of the rotor. Variable Name: rOG_4
rOG_4 = inv(R34)*rOG_3                           % Correct

% Absolute Angular velocity of frame 1 represented in frame 1 
% Variable Name: w1_1
w1_1 = [0;0;diff(alpha_,t)]                            % Correct


% Relative angular velocity of frame 2 relative to 1 represented in frame 2 
% Variable Name: w21_2
w21_2 = [diff(beta_,t);0;0]                           % Correct


% Relative angular velocity of frame 3 relative to 2 represented in frame 3
% Variable Name: w32_3
w32_3 = [0;0;diff(gamma_,t)]                            % Correct


% Relative angular velocity of frame 4 relative to 3 represented in frame 4
% Variable Name: w43_4
w43_4 = [0;0;diff(delta_,t)]                            % Correct


% Absolute angular velocity of frame 3 represented in frame 3 (Test)
% Variable Name: w3_3
w3_3 = w32_3 + inv(R23)*w21_2 + inv(R12*R23)*w1_1                           % Correct




% Velocity of center of mass of the frame represented in frame 3 (Test)
% Variable Name: rOG_3_dot
rOG_3_dot = diff(rOG_3,t) + cross(w3_3, rOG_3)                       % Correct 



% Acceleration of center of mass of the rotor represented in frame 4 (Test)
% Variable Name: rOG_4_dotdot
w4_4 = w43_4 + inv(R34)*w3_3;                                           % Correct
rOG_4_dot = diff(rOG_4,t) + cross(w4_4, rOG_4);                         % Correct
rOG_4_dotdot = diff(rOG_4_dot,t) + cross(w4_4, rOG_4_dot)               % Correct

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Newton-Euler Equation Setup
% Force due to gravity on the torus in frame zero
% Variable Name: FGrotor_0
FGrotor_0 = [0;0;-m_rot*g]                           % Correct


% Force due to gravity on the frame in frame zero
% Variable Name: FGframe_0
FGframe_0 = [0;0; -m_frame*g]                      % Correct 


% rotor angular momentum and its time derivative
% Variable Name: hGrotor_4
hGrotor_4 = IGrotor_4 * w4_4                     % Correct 


% Variable Name: hGrotor_4_dot (Test)
hGrotor_4_dot = diff(hGrotor_4,t) + cross(w4_4, hGrotor_4)                     % Correct 


% frame angular momentum and its time derivative
% Variable Name: hGframe_3
hGframe_3 = IGframe_3 * w3_3                     % Correct 


% Variable Name: hGframe_3_dot (Test)
hGframe_3_dot = diff(hGframe_3,t) + cross(w3_3, hGframe_3)                     % Correct 


%Symbolic variables for reaction Forces and Moments
syms F_Ox F_Oy F_Oz real
syms F_Gx F_Gy F_Gz real
syms M_Gx M_Gy M_Gz real
syms M_Ox M_Oy M_Oz real

% reaction forces
% Reaction force to the rotor at Point G
Frotor_4 = [F_Gx; F_Gy; F_Gz];
% Reaction moment to the rotor
Mrotor_4 = [M_Gx; M_Gy; M_Gz];
% Reaction force to the frame at Point O
Fframe_3 = [F_Ox; F_Oy; F_Oz];
% Reaction moment to the frame
Mframe_3 = [M_Ox; M_Oy; M_Oz];

%zero equations for the reaction forces/moments that don't exist (Prestest)
% variable name: zero_reaction
zero_reaction = [M_Gz; M_Ox; M_Oy; M_Oz] == 0                           % Correct      


%Linear NE Equations for the Rotor (Pretest)
% variable name: lin_NE_rotor
%lin_NE_rotor = Frotor_4  ==  m_rot*rOG_4_dotdot - inv(R04)*FGrotor_0                            % Correct

% %Linear NE Equations for the Frame (Test)
% variable name: lin_NE_frame
%rOG_3_dotdot = diff(rOG_3_dot,t) + cross(w3_3, rOG_3_dot)
%lin_NE_frame = Fframe_3  == m_frame*rOG_3_dotdot + R34*Frotor_4 - inv(R01*R12*R23)*FGframe_0

% %Angular NE Equations for the Rotor (Test)
% variable name: ang_NE_rotor
%ang_NE_rotor = Mrotor_4 == hGrotor_4_dot + cross(rOG_4, Frotor_4)

% %Angular NE Equations for the Frame (Test)
% variable name: ang_NE_frame
%ang_NE_frame = Mframe_3 == R34*Mrotor_4 + hGframe_3_dot + cross(rOG_3, Fframe_3)


lin_NE_rotor = [F_Gx == (2*R_maj_rotor*R_min_rotor^2*g*m_rotor*pi^2*sin(beta_(t))*(cos(delta_(t))*sin(gamma_(t)) + cos(gamma_(t))*sin(delta_(t))))/(2*R_maj_rotor*pi^2*R_min_rotor^2 + H_rot*pi*r_rot^2) - (R_maj_rotor*R_min_rotor^2*m_rotor*pi^2*(L*(cos(delta_(t))*(sin(gamma_(t))*diff(beta_(t), t, 2) + cos(gamma_(t))*diff(gamma_(t), t)*diff(beta_(t), t) - cos(gamma_(t))*sin(beta_(t))*diff(alpha_(t), t, 2) - cos(beta_(t))*cos(gamma_(t))*diff(beta_(t), t)*diff(alpha_(t), t) + sin(beta_(t))*sin(gamma_(t))*diff(gamma_(t), t)*diff(alpha_(t), t)) + sin(delta_(t))*(cos(gamma_(t))*diff(beta_(t), t, 2) - sin(gamma_(t))*diff(gamma_(t), t)*diff(beta_(t), t) + sin(beta_(t))*sin(gamma_(t))*diff(alpha_(t), t, 2) + cos(beta_(t))*sin(gamma_(t))*diff(beta_(t), t)*diff(alpha_(t), t) + cos(gamma_(t))*sin(beta_(t))*diff(gamma_(t), t)*diff(alpha_(t), t)) + cos(delta_(t))*(cos(gamma_(t))*diff(beta_(t), t) + sin(beta_(t))*sin(gamma_(t))*diff(alpha_(t), t))*diff(delta_(t), t) - sin(delta_(t))*(sin(gamma_(t))*diff(beta_(t), t) - cos(gamma_(t))*sin(beta_(t))*diff(alpha_(t), t))*diff(delta_(t), t)) - L*(cos(delta_(t))*(cos(gamma_(t))*diff(beta_(t), t) + sin(beta_(t))*sin(gamma_(t))*diff(alpha_(t), t)) - sin(delta_(t))*(sin(gamma_(t))*diff(beta_(t), t) - cos(gamma_(t))*sin(beta_(t))*diff(alpha_(t), t)))*(cos(beta_(t))*diff(alpha_(t), t) + diff(delta_(t), t) + diff(gamma_(t), t)))*2)/(2*R_maj_rotor*pi^2*R_min_rotor^2 + H_rot*pi*r_rot^2); F_Gy == (R_maj_rotor*R_min_rotor^2*m_rotor*pi^2*(L*(- cos(delta_(t))*(cos(gamma_(t))*diff(beta_(t), t, 2) - sin(gamma_(t))*diff(gamma_(t), t)*diff(beta_(t), t) + sin(beta_(t))*sin(gamma_(t))*diff(alpha_(t), t, 2) + cos(beta_(t))*sin(gamma_(t))*diff(beta_(t), t)*diff(alpha_(t), t) + cos(gamma_(t))*sin(beta_(t))*diff(gamma_(t), t)*diff(alpha_(t), t)) + sin(delta_(t))*(sin(gamma_(t))*diff(beta_(t), t, 2) + cos(gamma_(t))*diff(gamma_(t), t)*diff(beta_(t), t) - cos(gamma_(t))*sin(beta_(t))*diff(alpha_(t), t, 2) - cos(beta_(t))*cos(gamma_(t))*diff(beta_(t), t)*diff(alpha_(t), t) + sin(beta_(t))*sin(gamma_(t))*diff(gamma_(t), t)*diff(alpha_(t), t)) + cos(delta_(t))*(sin(gamma_(t))*diff(beta_(t), t) - cos(gamma_(t))*sin(beta_(t))*diff(alpha_(t), t))*diff(delta_(t), t) + sin(delta_(t))*(cos(gamma_(t))*diff(beta_(t), t) + sin(beta_(t))*sin(gamma_(t))*diff(alpha_(t), t))*diff(delta_(t), t)) - L*(cos(delta_(t))*(sin(gamma_(t))*diff(beta_(t), t) - cos(gamma_(t))*sin(beta_(t))*diff(alpha_(t), t)) + sin(delta_(t))*(cos(gamma_(t))*diff(beta_(t), t) + sin(beta_(t))*sin(gamma_(t))*diff(alpha_(t), t)))*(cos(beta_(t))*diff(alpha_(t), t) + diff(delta_(t), t) + diff(gamma_(t), t)))*2)/(2*R_maj_rotor*pi^2*R_min_rotor^2 + H_rot*pi*r_rot^2) + (2*R_maj_rotor*R_min_rotor^2*g*m_rotor*pi^2*sin(beta_(t))*(cos(delta_(t))*cos(gamma_(t)) - sin(delta_(t))*sin(gamma_(t))))/(2*R_maj_rotor*pi^2*R_min_rotor^2 + H_rot*pi*r_rot^2); F_Gz == (2*R_maj_rotor*R_min_rotor^2*g*m_rotor*pi^2*cos(beta_(t)))/(2*R_maj_rotor*pi^2*R_min_rotor^2 + H_rot*pi*r_rot^2) - (R_maj_rotor*R_min_rotor^2*m_rotor*pi^2*(L*(cos(delta_(t))*(sin(gamma_(t))*diff(beta_(t), t) - cos(gamma_(t))*sin(beta_(t))*diff(alpha_(t), t)) + sin(delta_(t))*(cos(gamma_(t))*diff(beta_(t), t) + sin(beta_(t))*sin(gamma_(t))*diff(alpha_(t), t)))^2 + L*(cos(delta_(t))*(cos(gamma_(t))*diff(beta_(t), t) + sin(beta_(t))*sin(gamma_(t))*diff(alpha_(t), t)) - sin(delta_(t))*(sin(gamma_(t))*diff(beta_(t), t) - cos(gamma_(t))*sin(beta_(t))*diff(alpha_(t), t)))^2)*2)/(2*R_maj_rotor*pi^2*R_min_rotor^2 + H_rot*pi*r_rot^2)]
lin_NE_frame = [F_Ox == F_Gx*cos(delta_(t)) - (L*(sin(gamma_(t))*diff(beta_(t), t, 2) + cos(gamma_(t))*diff(gamma_(t), t)*diff(beta_(t), t) - cos(gamma_(t))*sin(beta_(t))*diff(alpha_(t), t, 2) - cos(beta_(t))*cos(gamma_(t))*diff(beta_(t), t)*diff(alpha_(t), t) + sin(beta_(t))*sin(gamma_(t))*diff(gamma_(t), t)*diff(alpha_(t), t)) - L*(cos(gamma_(t))*diff(beta_(t), t) + sin(beta_(t))*sin(gamma_(t))*diff(alpha_(t), t))*(cos(beta_(t))*diff(alpha_(t), t) + diff(gamma_(t), t)))*((2*R_maj_H_tor*R_min_H_tor^2*m_frame*pi^2)/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2) - (pi*m_frame*r_rot^2*(H_rot - h_rod))/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2) + (2*R_maj_V_tor*R_min_V_tor^2*m_frame*pi^2)/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2) + (pi*H_rot*m_rotor*r_rot^2)/(2*R_maj_rotor*pi^2*R_min_rotor^2 + H_rot*pi*r_rot^2)) - F_Gy*sin(delta_(t)) + g*sin(beta_(t))*sin(gamma_(t))*((2*R_maj_H_tor*R_min_H_tor^2*m_frame*pi^2)/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2) - (pi*m_frame*r_rot^2*(H_rot - h_rod))/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2) + (2*R_maj_V_tor*R_min_V_tor^2*m_frame*pi^2)/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2) + (pi*H_rot*m_rotor*r_rot^2)/(2*R_maj_rotor*pi^2*R_min_rotor^2 + H_rot*pi*r_rot^2)); F_Oy == F_Gy*cos(delta_(t)) - (L*(cos(gamma_(t))*diff(beta_(t), t, 2) - sin(gamma_(t))*diff(gamma_(t), t)*diff(beta_(t), t) + sin(beta_(t))*sin(gamma_(t))*diff(alpha_(t), t, 2) + cos(beta_(t))*sin(gamma_(t))*diff(beta_(t), t)*diff(alpha_(t), t) + cos(gamma_(t))*sin(beta_(t))*diff(gamma_(t), t)*diff(alpha_(t), t)) + L*(sin(gamma_(t))*diff(beta_(t), t) - cos(gamma_(t))*sin(beta_(t))*diff(alpha_(t), t))*(cos(beta_(t))*diff(alpha_(t), t) + diff(gamma_(t), t)))*((2*R_maj_H_tor*R_min_H_tor^2*m_frame*pi^2)/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2) - (pi*m_frame*r_rot^2*(H_rot - h_rod))/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2) + (2*R_maj_V_tor*R_min_V_tor^2*m_frame*pi^2)/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2) + (pi*H_rot*m_rotor*r_rot^2)/(2*R_maj_rotor*pi^2*R_min_rotor^2 + H_rot*pi*r_rot^2)) + F_Gx*sin(delta_(t)) + g*cos(gamma_(t))*sin(beta_(t))*((2*R_maj_H_tor*R_min_H_tor^2*m_frame*pi^2)/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2) - (pi*m_frame*r_rot^2*(H_rot - h_rod))/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2) + (2*R_maj_V_tor*R_min_V_tor^2*m_frame*pi^2)/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2) + (pi*H_rot*m_rotor*r_rot^2)/(2*R_maj_rotor*pi^2*R_min_rotor^2 + H_rot*pi*r_rot^2)); F_Oz == F_Gz - (L*(cos(gamma_(t))*diff(beta_(t), t) + sin(beta_(t))*sin(gamma_(t))*diff(alpha_(t), t))^2 + L*(sin(gamma_(t))*diff(beta_(t), t) - cos(gamma_(t))*sin(beta_(t))*diff(alpha_(t), t))^2)*((2*R_maj_H_tor*R_min_H_tor^2*m_frame*pi^2)/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2) - (pi*m_frame*r_rot^2*(H_rot - h_rod))/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2) + (2*R_maj_V_tor*R_min_V_tor^2*m_frame*pi^2)/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2) + (pi*H_rot*m_rotor*r_rot^2)/(2*R_maj_rotor*pi^2*R_min_rotor^2 + H_rot*pi*r_rot^2)) + g*cos(beta_(t))*((2*R_maj_H_tor*R_min_H_tor^2*m_frame*pi^2)/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2) - (pi*m_frame*r_rot^2*(H_rot - h_rod))/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2) + (2*R_maj_V_tor*R_min_V_tor^2*m_frame*pi^2)/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2) + (pi*H_rot*m_rotor*r_rot^2)/(2*R_maj_rotor*pi^2*R_min_rotor^2 + H_rot*pi*r_rot^2))]
ang_NE_rotor = [M_Gx == (R_maj_rotor*R_min_rotor^2*m_rotor*pi^2*(cos(delta_(t))*(sin(gamma_(t))*diff(beta_(t), t) - cos(gamma_(t))*sin(beta_(t))*diff(alpha_(t), t)) + sin(delta_(t))*(cos(gamma_(t))*diff(beta_(t), t) + sin(beta_(t))*sin(gamma_(t))*diff(alpha_(t), t)))*(R_maj_rotor^2/2 + (5*R_min_rotor^2)/8)*(cos(beta_(t))*diff(alpha_(t), t) + diff(delta_(t), t) + diff(gamma_(t), t))*2)/(2*R_maj_rotor*pi^2*R_min_rotor^2 + H_rot*pi*r_rot^2) - (R_maj_rotor*R_min_rotor^2*m_rotor*pi^2*(R_maj_rotor^2 + (3*R_min_rotor^2)/4)*(cos(delta_(t))*(sin(gamma_(t))*diff(beta_(t), t) - cos(gamma_(t))*sin(beta_(t))*diff(alpha_(t), t)) + sin(delta_(t))*(cos(gamma_(t))*diff(beta_(t), t) + sin(beta_(t))*sin(gamma_(t))*diff(alpha_(t), t)))*(cos(beta_(t))*diff(alpha_(t), t) + diff(delta_(t), t) + diff(gamma_(t), t))*2)/(2*R_maj_rotor*pi^2*R_min_rotor^2 + H_rot*pi*r_rot^2) - (R_maj_rotor*R_min_rotor^2*m_rotor*pi^2*(R_maj_rotor^2/2 + (5*R_min_rotor^2)/8)*(- cos(delta_(t))*(cos(gamma_(t))*diff(beta_(t), t, 2) - sin(gamma_(t))*diff(gamma_(t), t)*diff(beta_(t), t) + sin(beta_(t))*sin(gamma_(t))*diff(alpha_(t), t, 2) + cos(beta_(t))*sin(gamma_(t))*diff(beta_(t), t)*diff(alpha_(t), t) + cos(gamma_(t))*sin(beta_(t))*diff(gamma_(t), t)*diff(alpha_(t), t)) + sin(delta_(t))*(sin(gamma_(t))*diff(beta_(t), t, 2) + cos(gamma_(t))*diff(gamma_(t), t)*diff(beta_(t), t) - cos(gamma_(t))*sin(beta_(t))*diff(alpha_(t), t, 2) - cos(beta_(t))*cos(gamma_(t))*diff(beta_(t), t)*diff(alpha_(t), t) + sin(beta_(t))*sin(gamma_(t))*diff(gamma_(t), t)*diff(alpha_(t), t)) + cos(delta_(t))*(sin(gamma_(t))*diff(beta_(t), t) - cos(gamma_(t))*sin(beta_(t))*diff(alpha_(t), t))*diff(delta_(t), t) + sin(delta_(t))*(cos(gamma_(t))*diff(beta_(t), t) + sin(beta_(t))*sin(gamma_(t))*diff(alpha_(t), t))*diff(delta_(t), t))*2)/(2*R_maj_rotor*pi^2*R_min_rotor^2 + H_rot*pi*r_rot^2); M_Gy == (R_maj_rotor*R_min_rotor^2*m_rotor*pi^2*(cos(delta_(t))*(cos(gamma_(t))*diff(beta_(t), t) + sin(beta_(t))*sin(gamma_(t))*diff(alpha_(t), t)) - sin(delta_(t))*(sin(gamma_(t))*diff(beta_(t), t) - cos(gamma_(t))*sin(beta_(t))*diff(alpha_(t), t)))*(R_maj_rotor^2/2 + (5*R_min_rotor^2)/8)*(cos(beta_(t))*diff(alpha_(t), t) + diff(delta_(t), t) + diff(gamma_(t), t))*2)/(2*R_maj_rotor*pi^2*R_min_rotor^2 + H_rot*pi*r_rot^2) - (R_maj_rotor*R_min_rotor^2*m_rotor*pi^2*(R_maj_rotor^2 + (3*R_min_rotor^2)/4)*(cos(delta_(t))*(cos(gamma_(t))*diff(beta_(t), t) + sin(beta_(t))*sin(gamma_(t))*diff(alpha_(t), t)) - sin(delta_(t))*(sin(gamma_(t))*diff(beta_(t), t) - cos(gamma_(t))*sin(beta_(t))*diff(alpha_(t), t)))*(cos(beta_(t))*diff(alpha_(t), t) + diff(delta_(t), t) + diff(gamma_(t), t))*2)/(2*R_maj_rotor*pi^2*R_min_rotor^2 + H_rot*pi*r_rot^2) - (R_maj_rotor*R_min_rotor^2*m_rotor*pi^2*(R_maj_rotor^2/2 + (5*R_min_rotor^2)/8)*(cos(delta_(t))*(sin(gamma_(t))*diff(beta_(t), t, 2) + cos(gamma_(t))*diff(gamma_(t), t)*diff(beta_(t), t) - cos(gamma_(t))*sin(beta_(t))*diff(alpha_(t), t, 2) - cos(beta_(t))*cos(gamma_(t))*diff(beta_(t), t)*diff(alpha_(t), t) + sin(beta_(t))*sin(gamma_(t))*diff(gamma_(t), t)*diff(alpha_(t), t)) + sin(delta_(t))*(cos(gamma_(t))*diff(beta_(t), t, 2) - sin(gamma_(t))*diff(gamma_(t), t)*diff(beta_(t), t) + sin(beta_(t))*sin(gamma_(t))*diff(alpha_(t), t, 2) + cos(beta_(t))*sin(gamma_(t))*diff(beta_(t), t)*diff(alpha_(t), t) + cos(gamma_(t))*sin(beta_(t))*diff(gamma_(t), t)*diff(alpha_(t), t)) + cos(delta_(t))*(cos(gamma_(t))*diff(beta_(t), t) + sin(beta_(t))*sin(gamma_(t))*diff(alpha_(t), t))*diff(delta_(t), t) - sin(delta_(t))*(sin(gamma_(t))*diff(beta_(t), t) - cos(gamma_(t))*sin(beta_(t))*diff(alpha_(t), t))*diff(delta_(t), t))*2)/(2*R_maj_rotor*pi^2*R_min_rotor^2 + H_rot*pi*r_rot^2); M_Gz == (2*R_maj_rotor*R_min_rotor^2*m_rotor*pi^2*(R_maj_rotor^2 + (3*R_min_rotor^2)/4)*(cos(beta_(t))*diff(alpha_(t), t, 2) - sin(beta_(t))*diff(beta_(t), t)*diff(alpha_(t), t) + diff(delta_(t), t, 2) + diff(gamma_(t), t, 2)))/(2*R_maj_rotor*pi^2*R_min_rotor^2 + H_rot*pi*r_rot^2)]
ang_NE_frame = [M_Ox == ((2*R_maj_H_tor*R_min_H_tor^2*m_frame*pi^2*(R_maj_H_tor^2/2 + (5*R_min_H_tor^2)/8))/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2) - ((pi*m_frame*r_rot^2*(H_rot - h_rod))/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2) - (pi*H_rot*m_rotor*r_rot^2)/(2*R_maj_rotor*pi^2*R_min_rotor^2 + H_rot*pi*r_rot^2))*(h_rod^2/12 + r_rot^2/4) + (2*R_maj_V_tor*R_min_V_tor^2*m_frame*pi^2*(R_maj_V_tor^2/2 + (5*R_min_V_tor^2)/8))/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2))*(cos(gamma_(t))*diff(beta_(t), t, 2) - sin(gamma_(t))*diff(gamma_(t), t)*diff(beta_(t), t) + sin(beta_(t))*sin(gamma_(t))*diff(alpha_(t), t, 2) + cos(beta_(t))*sin(gamma_(t))*diff(beta_(t), t)*diff(alpha_(t), t) + cos(gamma_(t))*sin(beta_(t))*diff(gamma_(t), t)*diff(alpha_(t), t)) + M_Gx*cos(delta_(t)) - F_Oy*L - M_Gy*sin(delta_(t)) + (sin(gamma_(t))*diff(beta_(t), t) - cos(gamma_(t))*sin(beta_(t))*diff(alpha_(t), t))*(cos(beta_(t))*diff(alpha_(t), t) + diff(gamma_(t), t))*((2*R_maj_V_tor*R_min_V_tor^2*m_frame*pi^2*(R_maj_V_tor^2 + (3*R_min_V_tor^2)/4))/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2) - ((pi*m_frame*r_rot^2*(H_rot - h_rod))/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2) - (pi*H_rot*m_rotor*r_rot^2)/(2*R_maj_rotor*pi^2*R_min_rotor^2 + H_rot*pi*r_rot^2))*(h_rod^2/12 + r_rot^2/4) + (2*R_maj_H_tor*R_min_H_tor^2*m_frame*pi^2*(R_maj_H_tor^2/2 + (5*R_min_H_tor^2)/8))/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2)) - (sin(gamma_(t))*diff(beta_(t), t) - cos(gamma_(t))*sin(beta_(t))*diff(alpha_(t), t))*(cos(beta_(t))*diff(alpha_(t), t) + diff(gamma_(t), t))*((2*R_maj_H_tor*R_min_H_tor^2*m_frame*pi^2*(R_maj_H_tor^2 + (3*R_min_H_tor^2)/4))/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2) - (r_rot^2*((pi*m_frame*r_rot^2*(H_rot - h_rod))/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2) - (pi*H_rot*m_rotor*r_rot^2)/(2*R_maj_rotor*pi^2*R_min_rotor^2 + H_rot*pi*r_rot^2)))/2 + (2*R_maj_V_tor*R_min_V_tor^2*m_frame*pi^2*(R_maj_V_tor^2/2 + (5*R_min_V_tor^2)/8))/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2)); M_Oy == M_Gy*cos(delta_(t)) - ((2*R_maj_V_tor*R_min_V_tor^2*m_frame*pi^2*(R_maj_V_tor^2 + (3*R_min_V_tor^2)/4))/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2) - ((pi*m_frame*r_rot^2*(H_rot - h_rod))/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2) - (pi*H_rot*m_rotor*r_rot^2)/(2*R_maj_rotor*pi^2*R_min_rotor^2 + H_rot*pi*r_rot^2))*(h_rod^2/12 + r_rot^2/4) + (2*R_maj_H_tor*R_min_H_tor^2*m_frame*pi^2*(R_maj_H_tor^2/2 + (5*R_min_H_tor^2)/8))/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2))*(sin(gamma_(t))*diff(beta_(t), t, 2) + cos(gamma_(t))*diff(gamma_(t), t)*diff(beta_(t), t) - cos(gamma_(t))*sin(beta_(t))*diff(alpha_(t), t, 2) - cos(beta_(t))*cos(gamma_(t))*diff(beta_(t), t)*diff(alpha_(t), t) + sin(beta_(t))*sin(gamma_(t))*diff(gamma_(t), t)*diff(alpha_(t), t)) + F_Ox*L + M_Gx*sin(delta_(t)) - (cos(gamma_(t))*diff(beta_(t), t) + sin(beta_(t))*sin(gamma_(t))*diff(alpha_(t), t))*(cos(beta_(t))*diff(alpha_(t), t) + diff(gamma_(t), t))*((2*R_maj_H_tor*R_min_H_tor^2*m_frame*pi^2*(R_maj_H_tor^2 + (3*R_min_H_tor^2)/4))/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2) - (r_rot^2*((pi*m_frame*r_rot^2*(H_rot - h_rod))/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2) - (pi*H_rot*m_rotor*r_rot^2)/(2*R_maj_rotor*pi^2*R_min_rotor^2 + H_rot*pi*r_rot^2)))/2 + (2*R_maj_V_tor*R_min_V_tor^2*m_frame*pi^2*(R_maj_V_tor^2/2 + (5*R_min_V_tor^2)/8))/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2)) + (cos(gamma_(t))*diff(beta_(t), t) + sin(beta_(t))*sin(gamma_(t))*diff(alpha_(t), t))*(cos(beta_(t))*diff(alpha_(t), t) + diff(gamma_(t), t))*((2*R_maj_H_tor*R_min_H_tor^2*m_frame*pi^2*(R_maj_H_tor^2/2 + (5*R_min_H_tor^2)/8))/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2) - ((pi*m_frame*r_rot^2*(H_rot - h_rod))/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2) - (pi*H_rot*m_rotor*r_rot^2)/(2*R_maj_rotor*pi^2*R_min_rotor^2 + H_rot*pi*r_rot^2))*(h_rod^2/12 + r_rot^2/4) + (2*R_maj_V_tor*R_min_V_tor^2*m_frame*pi^2*(R_maj_V_tor^2/2 + (5*R_min_V_tor^2)/8))/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2)); M_Oz == M_Gz + (cos(beta_(t))*diff(alpha_(t), t, 2) - sin(beta_(t))*diff(beta_(t), t)*diff(alpha_(t), t) + diff(gamma_(t), t, 2))*((2*R_maj_H_tor*R_min_H_tor^2*m_frame*pi^2*(R_maj_H_tor^2 + (3*R_min_H_tor^2)/4))/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2) - (r_rot^2*((pi*m_frame*r_rot^2*(H_rot - h_rod))/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2) - (pi*H_rot*m_rotor*r_rot^2)/(2*R_maj_rotor*pi^2*R_min_rotor^2 + H_rot*pi*r_rot^2)))/2 + (2*R_maj_V_tor*R_min_V_tor^2*m_frame*pi^2*(R_maj_V_tor^2/2 + (5*R_min_V_tor^2)/8))/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2)) - (cos(gamma_(t))*diff(beta_(t), t) + sin(beta_(t))*sin(gamma_(t))*diff(alpha_(t), t))*(sin(gamma_(t))*diff(beta_(t), t) - cos(gamma_(t))*sin(beta_(t))*diff(alpha_(t), t))*((2*R_maj_V_tor*R_min_V_tor^2*m_frame*pi^2*(R_maj_V_tor^2 + (3*R_min_V_tor^2)/4))/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2) - ((pi*m_frame*r_rot^2*(H_rot - h_rod))/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2) - (pi*H_rot*m_rotor*r_rot^2)/(2*R_maj_rotor*pi^2*R_min_rotor^2 + H_rot*pi*r_rot^2))*(h_rod^2/12 + r_rot^2/4) + (2*R_maj_H_tor*R_min_H_tor^2*m_frame*pi^2*(R_maj_H_tor^2/2 + (5*R_min_H_tor^2)/8))/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2)) + (cos(gamma_(t))*diff(beta_(t), t) + sin(beta_(t))*sin(gamma_(t))*diff(alpha_(t), t))*(sin(gamma_(t))*diff(beta_(t), t) - cos(gamma_(t))*sin(beta_(t))*diff(alpha_(t), t))*((2*R_maj_H_tor*R_min_H_tor^2*m_frame*pi^2*(R_maj_H_tor^2/2 + (5*R_min_H_tor^2)/8))/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2) - ((pi*m_frame*r_rot^2*(H_rot - h_rod))/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2) - (pi*H_rot*m_rotor*r_rot^2)/(2*R_maj_rotor*pi^2*R_min_rotor^2 + H_rot*pi*r_rot^2))*(h_rod^2/12 + r_rot^2/4) + (2*R_maj_V_tor*R_min_V_tor^2*m_frame*pi^2*(R_maj_V_tor^2/2 + (5*R_min_V_tor^2)/8))/(2*R_maj_H_tor*pi^2*R_min_H_tor^2 + 2*R_maj_V_tor*pi^2*R_min_V_tor^2 - pi*(H_rot - h_rod)*r_rot^2))]



% By this point, we have 16 unknown.
% Make sure you have 16 equation in symbolic variable equations. From top to bottom
% in sequence of zero_reaction, rotor linear equation, frame linear equation
% rotor angular equation, frame angular equation
% The following variable is not tested, but it will help you with the rest of your modelling. 
% It is initially set to zeros for technical reasons. Change it accordingly.
% variable name: equations
% Remove comment of the next line once you have all the above equations and delete the last line with the zeros
equations = [zero_reaction; lin_NE_rotor; lin_NE_frame; ang_NE_rotor; ang_NE_frame];

eqs = [subs(ang_NE_rotor(3), M_Gz, 0); subs(ang_NE_frame, [M_Oz, M_Ox, M_Oy], [0,0,0])]

%%
syms alpha_d(t) alpha_dd(t) delta_d(t) delta_dd(t) beta_d(t) beta_dd(t) gamma_d(t) gamma_dd(t)
eqs = subs(eqs, [diff(alpha_,t,t), diff(beta_,t,t), diff(gamma_,t,t), diff(delta_,t,t)], [alpha_dd beta_dd delta_dd gamma_dd])
[a,b] = equationsToMatrix(eqs, [alpha_dd beta_dd delta_dd gamma_dd])

function xdot = eom(t,x)

    tau1 = 10*sin(t);
    tau2 = 5*cos(t);

    theta = x(1);
    alpha = x(2);
    theta_d = x(3);
    alpha_d = x(4);
    
    g = 9.8;
       
    l2 = 0.9;
    m2 = 1;
    
    Ix2 = 0.5;
    Iy2 = Ix2;
    Iz2 = 1;
    
    Iz1 = 1;
        
    f1_num = tau1 + alpha_d*theta_d*sin(alpha)*cos(alpha)*(2*(Iz2-Ix2)-m2*l2^2/4);
    f1_den = Iz1 + Iz2*cos(alpha)^2+(Ix2+m2*l2^2/4)*sin(alpha)^2;
    
    f1 = f1_num/f1_den;
    
    f2_num = tau2+theta_d^2*sin(alpha)*cos(alpha)*(Ix2-Iz2+m2*l2^2/4)-m2*l2*g*sin(alpha)/2;
    f2_den = Iy2 + m2*l2^2/4;
    
    f2 = f2_num/f2_den;
    
    xdot = [theta_d;
            alpha_d;
            f1
            f2]; 

    %%% Gyro Simulation
    alpha_ = x(1);
    beta_ = x(2);
    gamma_ = x(3);
    delta_ = x(4);
    alpha_d = x(5);
    beta_d = x(6);
    gamma_d = x(7);
    delta_d = x(8);





end



