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
                             


                          
%%%%%% Rotor %%%%%

% Desity of the rotor [kg/m^3] (Pretest).
% Variable Name: rho_rotor
 



      

% Central rod considered as a cylinder belonging to the frame

% Recalculated mass of frame [kg] (Test).
% Variable Name: m_frame 
             

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inertia tensor
% Inerta tensor of the 3 frame components in frame 3 

% Inertia tensor of horizontal torus about its center of mass expressed in frame 3 (Pretest). 
% Variable Name: IGtorH_3


% Inertia tensor of vertical torus about its center of mass expressed in frame 3 (Test)
% Variable Name: IGtorV_3


% Inertia tensor of the rod about its center of mass expressed in frame 3 
% Variable Name: IGrod_3


% Inertia tensor of the frame about its center of mass expressed in frame 3 (Test)
% Variable Name: IGframe_3


% Inertia tensor of rotor about its center of mass in frame 4 
% Variable Name: IGrotor_4


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up rotational matrices



%Rotation matrix from frame 4 to frame 0 (Pretest).
%Variable Name: R04


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Kinematics
% Centers of mass positions from origin
% Center of mass position of the frame. Variable Name: rOG_3

% Center of mass position of the rotor. Variable Name: rOG_4


% Absolute Angular velocity of frame 1 represented in frame 1 
% Variable Name: w1_1


% Relative angular velocity of frame 2 relative to 1 represented in frame 2 
% Variable Name: w21_2


% Relative angular velocity of frame 3 relative to 2 represented in frame 3
% Variable Name: w32_3


% Relative angular velocity of frame 4 relative to 3 represented in frame 4
% Variable Name: w43_4




% Absolute angular velocity of frame 3 represented in frame 3 (Test)
% Variable Name: w3_3




% Velocity of center of mass of the frame represented in frame 3 (Test)
% Variable Name: rOG_3_dot





% Acceleration of center of mass of the rotor represented in frame 4 (Test)
% Variable Name: rOG_4_dotdot


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Newton-Euler Equation Setup
% Force due to gravity on the torus in frame zero
% Variable Name: FGrotor_0


% Force due to gravity on the frame in frame zero
% Variable Name: FGframe_0


% rotor angular momentum and its time derivative
% Variable Name: hGrotor_4


% Variable Name: hGrotor_4_dot (Test)


% frame angular momentum and its time derivative
% Variable Name: hGframe_3


% Variable Name: hGframe_3_dot (Test)


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


%Linear NE Equations for the Rotor (Pretest)
% variable name: lin_NE_rotor


% %Linear NE Equations for the Frame (Test)
% variable name: lin_NE_frame


% %Angular NE Equations for the Rotor (Test)
% variable name: ang_NE_rotor


% %Angular NE Equations for the Frame (Test)
% variable name: ang_NE_frame


% By this point, we have 16 unknown.
% Make sure you have 16 equation in symbolic variable equations. From top to bottom
% in sequence of zero_reaction, rotor linear equation, frame linear equation
% rotor angular equation, frame angular equation
% The following variable is not tested, but it will help you with the rest of your modelling. 
% It is initially set to zeros for technical reasons. Change it accordingly.
% variable name: equations
% Remove comment of the next line once you have all the above equations and delete the last line with the zeros
%equations = [zero_reaction; lin_NE_rotor; lin_NE_frame; ang_NE_rotor; ang_NE_frame];

equations = zeros(4,1);







