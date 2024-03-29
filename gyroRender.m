close all
clear all
clc

m_frame = 0.023;                                            % Mass of frame [kg] (Measured)
m_rotor = 0.045;                                            % Mass of rotor [kg] (Measured)
H_rot = 60*10^-3;
r_rot = 1.555*10^-3;
R_min_rotor = 3.305*10^-3;
R_maj_rotor = 28*10^-3;
h_rod = 88*10^-3;
R_min_H_tor = 1.155*10^-3;
R_maj_H_tor = 36*10^-3;
R_sph = 2.365*10^-3;
R_min_V_tor = 1.715*10^-3;
R_maj_V_tor = 31.5*10^-3;
L = h_rod/2; 
R_shp = 2.365*10^-3;

[theta,phi] = meshgrid(linspace(0,2*pi,50));
r_v = R_min_V_tor;
R_v = R_maj_V_tor;
z_V_tor = (R_v + r_v*cos(theta)).*cos(phi) + L;
x_V_tor = (R_v + r_v*cos(theta)).*sin(phi);
y_V_tor = r_v*sin(theta);
surf(x_V_tor,y_V_tor,z_V_tor)
%shading interp
axis equal
hold on

[theta,phi] = meshgrid(linspace(0,2*pi,50));
r_h = R_min_H_tor;
R_h = R_maj_H_tor;
x_H_tor = (R_h + r_h*cos(theta)).*cos(phi);
y_H_tor = (R_h + r_h*cos(theta)).*sin(phi);
z_H_tor = r_h*sin(theta) + L;
surf(x_H_tor,y_H_tor,z_H_tor)


[x_rod,y_rod,z_rod] = cylinder(r_rot, 100);
z_rod = z_rod * h_rod;
surf(x_rod,y_rod,z_rod)

[x_rotor,y_rotor,z_rotor] = cylinder(R_maj_rotor, 100);
z_rotor = z_rotor * R_min_rotor;
z_rotor = z_rotor - R_min_rotor/2 + L;
surf(x_rotor,y_rotor,z_rotor)
fill3(x_rotor(1,:),y_rotor(1,:),z_rotor(1,:), 'b');
fill3(x_rotor(2,:),y_rotor(2,:),z_rotor(2,:),'b');

[x_T_sph,y_T_sph,z_T_sph] = sphere;
x_T_sph = x_T_sph* R_shp;
y_T_sph = y_T_sph* R_shp;
z_T_sph = z_T_sph* R_shp;

[x_B_sph,y_B_sph,z_B_sph] = sphere;
x_B_sph = x_B_sph* R_shp;
y_B_sph = y_B_sph* R_shp;
z_B_sph = z_B_sph* R_shp;

z_T_sph = z_T_sph + h_rod;

surf(x_T_sph,y_T_sph,z_T_sph);
surf(x_B_sph,y_B_sph,z_B_sph);

