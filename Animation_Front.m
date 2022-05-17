function two_dof_animation()
    clear all
    close all
    clc
    %% DO I WANT TO RECORD THE VIDEO
    VIDEO = 1;
    
    %% SETUP THE PROBLEM
    X_init = [0;0.2;0;0;0;0;0;2000];                         % initial conditions
    %Very Important I.c, try to understand what variales means, important
    %for gyroscope    
    tspan = [0 9];                                 % start and finish times
    options = odeset('RelTol',1e-7,'AbsTol',1e-7'); % solver options,optional
    sol = ode45(@eom,tspan,X_init,options);         % SOLVE the eoms
    

    %% EVAULATE THE SOLUTION
    %dt is equally time steps, better to have constant time step
    dt = 0.03;                                  % set time step                        
    t = tspan(1):dt:tspan(2);                   % creat time vector
    X = deval(sol,t);                           % deval
    %% PLOT THE STATES
    plot(t,X)  %%X is the solution that we obtained before
    xlabel('time')
    ylabel('states')

    %% Create Gyroscope
%Parameters
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

%Torus Vertical
[theta,phi] = meshgrid(linspace(0,2*pi,50));
r_v = R_min_V_tor;
R_v = R_maj_V_tor;
z_V_tor = (R_v + r_v*cos(theta)).*cos(phi) + L;
x_V_tor = (R_v + r_v*cos(theta)).*sin(phi);
y_V_tor = r_v*sin(theta);
surf(x_V_tor,y_V_tor,z_V_tor)
shading interp
axis equal
hold on

%Torus Horizontal
[theta,phi] = meshgrid(linspace(0,2*pi,50));
r_h = R_min_H_tor;
R_h = R_maj_H_tor;
x_H_tor = (R_h + r_h*cos(theta)).*cos(phi);
y_H_tor = (R_h + r_h*cos(theta)).*sin(phi);
z_H_tor = r_h*sin(theta) + L;
surf(x_H_tor,y_H_tor,z_H_tor)

%Rod
[x_rod,y_rod,z_rod] = cylinder(r_rot, 100);
z_rod = z_rod * h_rod;
surf(x_rod,y_rod,z_rod)

%Rotor
[x_rotor,y_rotor,z_rotor] = cylinder(R_maj_rotor, 100);
z_rotor = z_rotor * R_min_rotor;
z_rotor = z_rotor - R_min_rotor/2 + L;
surf(x_rotor,y_rotor,z_rotor)
fill3(x_rotor(1,:),y_rotor(1,:),z_rotor(1,:), 'b');
fill3(x_rotor(2,:),y_rotor(2,:),z_rotor(2,:),'b');

%Sphere
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

    %% SETUP VIDEO IF REQUIRED
    if VIDEO
        fps = 1/dt;
        MyVideo = VideoWriter('two_dof_animation','MPEG-4');
        MyVideo.FrameRate = fps;
        open(MyVideo);
    end

    %% CREATE ANIMATION
    handle = figure;
    hold on % ; grid on
    for i = 1:length(t)
        cla

        alpha = X(1,i);
        beta = X(2,i);
        gamma = X(3,i);
        delta = X(4,i);

       %Rotation frame
        R01 = [cos(alpha) -sin(alpha) 0;
               sin(alpha) cos(alpha) 0;                    
               0 0 1];
        R12 = [1 0 0;
               0 cos(beta) -sin(beta);
               0 sin(beta) cos(beta)];
        R23 = [cos(gamma) -sin(gamma) 0;                  
               sin(gamma) cos(gamma) 0;
               0 0 1];
        R34 = [cos(delta) -sin(delta) 0;
               sin(delta) cos(delta) 0;                 
               0 0 1];
        R02 = R01*R12;
        R03 = R01*R12*R23;
        R04 = R01*R12*R23*R34; 

        %Vertical Torus
        [x_V_tor_rotated,y_V_tor_rotated,z_V_tor_rotated] = rotation(x_V_tor,y_V_tor,z_V_tor,R03);
        %Horizontal Torus
        [x_H_tor_rotated,y_H_tor_rotated,z_H_tor_rotated] = rotation(x_H_tor,y_H_tor,z_H_tor,R03);
        %Rotor
        [x_rotor_rotated,y_rotor_rotated,z_rotor_rotated] = rotation(x_rotor,y_rotor,z_rotor,R04);
        %Rod
        [x_rod_rotated,y_rod_rotated,z_rod_rotated] = rotation(x_rod,y_rod,z_rod,R02);
        %Sphere_Top
        [x_T_sph_rotated,y_T_sph_rotated,z_T_sph_rotated] = rotation(x_T_sph,y_T_sph,z_T_sph,R02);
        %Sphere_Bottom
        [x_B_sph_rotated,y_B_sph_rotated,z_B_sph_rotated] = rotation(x_B_sph,y_B_sph,z_B_sph,R02);
        
        surf(x_V_tor_rotated,y_V_tor_rotated,z_V_tor_rotated,x_V_tor, 'EdgeColor','none')
        surf(x_H_tor_rotated,y_H_tor_rotated,z_H_tor_rotated,x_H_tor, 'EdgeColor','none')
        surf(x_rotor_rotated,y_rotor_rotated,z_rotor_rotated,x_rotor, 'EdgeColor','none')
        surf(x_rod_rotated,y_rod_rotated,z_rod_rotated,x_rod, 'EdgeColor','none')
        surf(x_T_sph_rotated,y_T_sph_rotated,z_T_sph_rotated,x_T_sph, 'EdgeColor','none')
        surf(x_B_sph_rotated,y_B_sph_rotated,z_B_sph_rotated,x_B_sph, 'EdgeColor','none')
        
        axis square
        view(90,0)
        axis([-1 1 -1 1 -1 1]*0.1)
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        if VIDEO
            writeVideo(MyVideo,getframe(handle));
        else
            pause(dt)
        end
    end

    if VIDEO
    close(MyVideo)
    end

end

%%
function [Xf,Yf,Zf]=rotation(Xi,Yi,Zi,R)

    I=size(Xi,1);
    J=size(Xi,2);

    Xf=zeros(I,J);
    Yf=zeros(I,J);
    Zf=zeros(I,J);

    for ii=1:I
        for jj=1:J
            vector=[Xi(ii,jj);Yi(ii,jj);Zi(ii,jj)];
            vector=R*vector;
                Xf(ii,jj)=vector(1);
                Yf(ii,jj)=vector(2);
                Zf(ii,jj)=vector(3);
        end
    end

end

function xdot = eom(t,x)
    alpha = x(1);
    beta = x(2);
    gamma = x(3);
    delta = x(4);
    alpha_d = x(5);
    beta_d = x(6);
    gamma_d = x(7);
    delta_d = x(8);


    Xdd = [(76423946057877712394330035967879099009524405524722382819492096050750000*beta_d*gamma_d*cos(gamma)^2 + 34251616451986631058281102331118006535803460925322758651373190328875000*beta_d*gamma_d*sin(gamma)^2 + 3926504415204791525763742537222856920440791396518143483995412096000000000*cos(gamma)*sin(beta)*sin(gamma) + 5603743290379891353298224402569617329729943735611563139766681600000000000*cos(delta)^2*cos(gamma)*sin(beta)*sin(gamma) + 5603743290379891353298224402569617329729943735611563139766681600000000000*cos(gamma)*sin(beta)*sin(delta)^2*sin(gamma) - 389349273897089330679102905535573103957075922287450521541745378047582896*alpha_d*beta_d*cos(beta)*cos(gamma)^2 - 431521603502980412015151839172334196430796866886850145709864283769457896*alpha_d*beta_d*cos(beta)*sin(gamma)^2 + 94950864807426754678937732083392106038328296936928220323750644314906250*beta_d*delta_d*cos(delta)^2*cos(gamma)^2 + 126456334998878665904778462366289547144880391480736929124695670784000000*beta_d*delta_d*cos(delta)^4*cos(gamma)^2 + 196732898538575103481763771790756260204605292022023999634732780122906250*beta_d*gamma_d*cos(delta)^2*cos(gamma)^2 + 126456334998878665904778462366289547144880391480736929124695670784000000*beta_d*gamma_d*cos(delta)^4*cos(gamma)^2 + 105235458706737837857068841181373116707212523169672996291361587918421875*beta_d*delta_d*cos(delta)^2*sin(gamma)^2 + 94950864807426754678937732083392106038328296936928220323750644314906250*beta_d*delta_d*cos(gamma)^2*sin(delta)^2 + 126456334998878665904778462366289547144880391480736929124695670784000000*beta_d*delta_d*cos(delta)^4*sin(gamma)^2 + 126456334998878665904778462366289547144880391480736929124695670784000000*beta_d*delta_d*cos(gamma)^2*sin(delta)^4 + 146393961150061211819617563484744388266235439347784841100286968526421875*beta_d*gamma_d*cos(delta)^2*sin(gamma)^2 + 196732898538575103481763771790756260204605292022023999634732780122906250*beta_d*gamma_d*cos(gamma)^2*sin(delta)^2 + 126456334998878665904778462366289547144880391480736929124695670784000000*beta_d*gamma_d*cos(delta)^4*sin(gamma)^2 + 126456334998878665904778462366289547144880391480736929124695670784000000*beta_d*gamma_d*cos(gamma)^2*sin(delta)^4 + 105235458706737837857068841181373116707212523169672996291361587918421875*beta_d*delta_d*sin(delta)^2*sin(gamma)^2 + 126456334998878665904778462366289547144880391480736929124695670784000000*beta_d*delta_d*sin(delta)^4*sin(gamma)^2 + 146393961150061211819617563484744388266235439347784841100286968526421875*beta_d*gamma_d*sin(delta)^2*sin(gamma)^2 + 126456334998878665904778462366289547144880391480736929124695670784000000*beta_d*gamma_d*sin(delta)^4*sin(gamma)^2 - 983284705050862066298482222806307362745187186508912157495488445855719798*alpha_d*beta_d*cos(beta)*cos(delta)^2*cos(gamma)^2 - 618952735733775353148350650853440531620292924742921993000815935363219456*alpha_d*beta_d*cos(beta)*cos(delta)^4*cos(gamma)^2 - 1033623642439375957960628431112319234683557039183151316029934257452204173*alpha_d*beta_d*cos(beta)*cos(delta)^2*sin(gamma)^2 - 983284705050862066298482222806307362745187186508912157495488445855719798*alpha_d*beta_d*cos(beta)*cos(gamma)^2*sin(delta)^2 - 618952735733775353148350650853440531620292924742921993000815935363219456*alpha_d*beta_d*cos(beta)*cos(delta)^4*sin(gamma)^2 - 618952735733775353148350650853440531620292924742921993000815935363219456*alpha_d*beta_d*cos(beta)*cos(gamma)^2*sin(delta)^4 - 1033623642439375957960628431112319234683557039183151316029934257452204173*alpha_d*beta_d*cos(beta)*sin(delta)^2*sin(gamma)^2 - 618952735733775353148350650853440531620292924742921993000815935363219456*alpha_d*beta_d*cos(beta)*sin(delta)^4*sin(gamma)^2 + 42172329605891081336048933636761092473720944599399624168118905721875000*alpha_d^2*cos(beta)*cos(gamma)*sin(beta)*sin(gamma) + 42172329605891081336048933636761092473720944599399624168118905721875000*alpha_d*gamma_d*cos(gamma)*sin(beta)*sin(gamma) + 252912669997757331809556924732579094289760782961473858249391341568000000*beta_d*delta_d*cos(delta)^2*cos(gamma)^2*sin(delta)^2 + 252912669997757331809556924732579094289760782961473858249391341568000000*beta_d*gamma_d*cos(delta)^2*cos(gamma)^2*sin(delta)^2 + 252912669997757331809556924732579094289760782961473858249391341568000000*beta_d*delta_d*cos(delta)^2*sin(delta)^2*sin(gamma)^2 + 252912669997757331809556924732579094289760782961473858249391341568000000*beta_d*gamma_d*cos(delta)^2*sin(delta)^2*sin(gamma)^2 - 1237905471467550706296701301706881063240585849485843986001631870726438912*alpha_d*beta_d*cos(beta)*cos(delta)^2*cos(gamma)^2*sin(delta)^2 - 1237905471467550706296701301706881063240585849485843986001631870726438912*alpha_d*beta_d*cos(beta)*cos(delta)^2*sin(delta)^2*sin(gamma)^2 + 50338937388513891662146208306011871938369852674239158534445811596484375*alpha_d^2*cos(beta)*cos(delta)^2*cos(gamma)*sin(beta)*sin(gamma) + 50338937388513891662146208306011871938369852674239158534445811596484375*alpha_d^2*cos(beta)*cos(gamma)*sin(beta)*sin(delta)^2*sin(gamma) - 10284593899311083178131109097981010668884226232744775967610943603515625*alpha_d*delta_d*cos(delta)^2*cos(gamma)*sin(beta)*sin(gamma) + 50338937388513891662146208306011871938369852674239158534445811596484375*alpha_d*gamma_d*cos(delta)^2*cos(gamma)*sin(beta)*sin(gamma) - 10284593899311083178131109097981010668884226232744775967610943603515625*alpha_d*delta_d*cos(gamma)*sin(beta)*sin(delta)^2*sin(gamma) + 50338937388513891662146208306011871938369852674239158534445811596484375*alpha_d*gamma_d*cos(gamma)*sin(beta)*sin(delta)^2*sin(gamma))/(16*(107921345109286508779897069447938048*cos(delta)^2 + 107921345109286508779897069447938048*sin(delta)^2 + 81033702656323036153327122834079417)*(179621720919313878356379449073236959*cos(gamma)^2*sin(beta) + 179621720919313878356379449073236959*sin(beta)*sin(gamma)^2 + 215842690218573017559794138895876096*cos(delta)^2*cos(gamma)^2*sin(beta) + 215842690218573017559794138895876096*cos(delta)^2*sin(beta)*sin(gamma)^2 + 215842690218573017559794138895876096*cos(gamma)^2*sin(beta)*sin(delta)^2 + 215842690218573017559794138895876096*sin(beta)*sin(delta)^2*sin(gamma)^2));
        -(34251616451986631058281102331118006535803460925322758651373190328875000*alpha_d*gamma_d*cos(gamma)^2*sin(beta) - 36250822691098013220624550051113113475693986403247587930539612869005803520*sin(beta)*sin(gamma)^2 - 105618545950988687849115948941837238703687924377436348081564432647197818880*cos(delta)^2*cos(gamma)^2*sin(beta) - 68901975684568924003213441800005129256665468370073916060536797394894848000*cos(delta)^4*cos(gamma)^2*sin(beta) - 100014802660608796495817724539267621373957980641824784941797751047197818880*cos(delta)^2*sin(beta)*sin(gamma)^2 - 105618545950988687849115948941837238703687924377436348081564432647197818880*cos(gamma)^2*sin(beta)*sin(delta)^2 - 68901975684568924003213441800005129256665468370073916060536797394894848000*cos(delta)^4*sin(beta)*sin(gamma)^2 - 68901975684568924003213441800005129256665468370073916060536797394894848000*cos(gamma)^2*sin(beta)*sin(delta)^4 - 100014802660608796495817724539267621373957980641824784941797751047197818880*sin(beta)*sin(delta)^2*sin(gamma)^2 - 68901975684568924003213441800005129256665468370073916060536797394894848000*sin(beta)*sin(delta)^4*sin(gamma)^2 - 40177327106302804746388292588335970396134777799765731414535024965005803520*cos(gamma)^2*sin(beta) + 76423946057877712394330035967879099009524405524722382819492096050750000*alpha_d*gamma_d*sin(beta)*sin(gamma)^2 - 137803951369137848006426883600010258513330936740147832121073594789789696000*cos(delta)^2*cos(gamma)^2*sin(beta)*sin(delta)^2 - 137803951369137848006426883600010258513330936740147832121073594789789696000*cos(delta)^2*sin(beta)*sin(delta)^2*sin(gamma)^2 + 42172329605891081336048933636761092473720944599399624168118905721875000*beta_d*gamma_d*cos(gamma)*sin(gamma) - 198634993525496890478435368420608094947496702980763693529245546720291448*alpha_d^2*cos(beta)*cos(gamma)^2*sin(beta) - 156462663919605809142386434783847002473775758381364069361126640998416448*alpha_d^2*cos(beta)*sin(beta)*sin(gamma)^2 - 10284593899311083178131109097981010668884226232744775967610943603515625*beta_d*delta_d*cos(gamma)*sin(delta)^2*sin(gamma) + 50338937388513891662146208306011871938369852674239158534445811596484375*beta_d*gamma_d*cos(gamma)*sin(delta)^2*sin(gamma) - 443614840644657373070505433813787423208660799917683237464823644462891149*alpha_d^2*cos(beta)*cos(delta)^2*cos(gamma)^2*sin(beta) - 246248200367448343621786094243575492237706266631092531938060132289609728*alpha_d^2*cos(beta)*cos(delta)^4*cos(gamma)^2*sin(beta) - 393275903256143481408359225507775551270290947243444078930377832866406774*alpha_d^2*cos(beta)*cos(delta)^2*sin(beta)*sin(gamma)^2 - 443614840644657373070505433813787423208660799917683237464823644462891149*alpha_d^2*cos(beta)*cos(gamma)^2*sin(beta)*sin(delta)^2 - 246248200367448343621786094243575492237706266631092531938060132289609728*alpha_d^2*cos(beta)*cos(delta)^4*sin(beta)*sin(gamma)^2 - 246248200367448343621786094243575492237706266631092531938060132289609728*alpha_d^2*cos(beta)*cos(gamma)^2*sin(beta)*sin(delta)^4 - 393275903256143481408359225507775551270290947243444078930377832866406774*alpha_d^2*cos(beta)*sin(beta)*sin(delta)^2*sin(gamma)^2 - 246248200367448343621786094243575492237706266631092531938060132289609728*alpha_d^2*cos(beta)*sin(beta)*sin(delta)^4*sin(gamma)^2 + 105235458706737837857068841181373116707212523169672996291361587918421875*alpha_d*delta_d*cos(delta)^2*cos(gamma)^2*sin(beta) + 126456334998878665904778462366289547144880391480736929124695670784000000*alpha_d*delta_d*cos(delta)^4*cos(gamma)^2*sin(beta) + 146393961150061211819617563484744388266235439347784841100286968526421875*alpha_d*gamma_d*cos(delta)^2*cos(gamma)^2*sin(beta) + 126456334998878665904778462366289547144880391480736929124695670784000000*alpha_d*gamma_d*cos(delta)^4*cos(gamma)^2*sin(beta) + 94950864807426754678937732083392106038328296936928220323750644314906250*alpha_d*delta_d*cos(delta)^2*sin(beta)*sin(gamma)^2 + 105235458706737837857068841181373116707212523169672996291361587918421875*alpha_d*delta_d*cos(gamma)^2*sin(beta)*sin(delta)^2 + 126456334998878665904778462366289547144880391480736929124695670784000000*alpha_d*delta_d*cos(delta)^4*sin(beta)*sin(gamma)^2 + 126456334998878665904778462366289547144880391480736929124695670784000000*alpha_d*delta_d*cos(gamma)^2*sin(beta)*sin(delta)^4 + 196732898538575103481763771790756260204605292022023999634732780122906250*alpha_d*gamma_d*cos(delta)^2*sin(beta)*sin(gamma)^2 + 146393961150061211819617563484744388266235439347784841100286968526421875*alpha_d*gamma_d*cos(gamma)^2*sin(beta)*sin(delta)^2 + 126456334998878665904778462366289547144880391480736929124695670784000000*alpha_d*gamma_d*cos(delta)^4*sin(beta)*sin(gamma)^2 + 126456334998878665904778462366289547144880391480736929124695670784000000*alpha_d*gamma_d*cos(gamma)^2*sin(beta)*sin(delta)^4 + 94950864807426754678937732083392106038328296936928220323750644314906250*alpha_d*delta_d*sin(beta)*sin(delta)^2*sin(gamma)^2 + 126456334998878665904778462366289547144880391480736929124695670784000000*alpha_d*delta_d*sin(beta)*sin(delta)^4*sin(gamma)^2 + 196732898538575103481763771790756260204605292022023999634732780122906250*alpha_d*gamma_d*sin(beta)*sin(delta)^2*sin(gamma)^2 + 126456334998878665904778462366289547144880391480736929124695670784000000*alpha_d*gamma_d*sin(beta)*sin(delta)^4*sin(gamma)^2 + 42172329605891081336048933636761092473720944599399624168118905721875000*alpha_d*beta_d*cos(beta)*cos(gamma)*sin(gamma) - 10284593899311083178131109097981010668884226232744775967610943603515625*beta_d*delta_d*cos(delta)^2*cos(gamma)*sin(gamma) + 50338937388513891662146208306011871938369852674239158534445811596484375*beta_d*gamma_d*cos(delta)^2*cos(gamma)*sin(gamma) - 492496400734896687243572188487150984475412533262185063876120264579219456*alpha_d^2*cos(beta)*cos(delta)^2*cos(gamma)^2*sin(beta)*sin(delta)^2 - 492496400734896687243572188487150984475412533262185063876120264579219456*alpha_d^2*cos(beta)*cos(delta)^2*sin(beta)*sin(delta)^2*sin(gamma)^2 + 252912669997757331809556924732579094289760782961473858249391341568000000*alpha_d*delta_d*cos(delta)^2*cos(gamma)^2*sin(beta)*sin(delta)^2 + 252912669997757331809556924732579094289760782961473858249391341568000000*alpha_d*gamma_d*cos(delta)^2*cos(gamma)^2*sin(beta)*sin(delta)^2 + 252912669997757331809556924732579094289760782961473858249391341568000000*alpha_d*delta_d*cos(delta)^2*sin(beta)*sin(delta)^2*sin(gamma)^2 + 252912669997757331809556924732579094289760782961473858249391341568000000*alpha_d*gamma_d*cos(delta)^2*sin(beta)*sin(delta)^2*sin(gamma)^2 + 50338937388513891662146208306011871938369852674239158534445811596484375*alpha_d*beta_d*cos(beta)*cos(delta)^2*cos(gamma)*sin(gamma) + 50338937388513891662146208306011871938369852674239158534445811596484375*alpha_d*beta_d*cos(beta)*cos(gamma)*sin(delta)^2*sin(gamma))/(16*(107921345109286508779897069447938048*cos(delta)^2 + 107921345109286508779897069447938048*sin(delta)^2 + 81033702656323036153327122834079417)*(179621720919313878356379449073236959*cos(gamma)^2 + 179621720919313878356379449073236959*sin(gamma)^2 + 215842690218573017559794138895876096*cos(delta)^2*cos(gamma)^2 + 215842690218573017559794138895876096*cos(delta)^2*sin(gamma)^2 + 215842690218573017559794138895876096*cos(gamma)^2*sin(delta)^2 + 215842690218573017559794138895876096*sin(delta)^2*sin(gamma)^2));
        -(1172503524893525*(beta_d*cos(gamma) + alpha_d*sin(beta)*sin(gamma))*(beta_d*sin(gamma) - alpha_d*cos(gamma)*sin(beta)))/2764574727905982;
        -(105639854939229361219673820582119727779756472841096806600042942774864523779750293250000*beta_d*gamma_d*cos(beta)*cos(gamma)^2 + 47345576616545498571000475772876628087827011255299812010776344678528514531329915125000*beta_d*gamma_d*cos(beta)*sin(gamma)^2 + 136530185549551498916930460463168183800193567078764874592309425437632189282209601224600*alpha_d^2*cos(gamma)*sin(beta)^3*sin(gamma)^3 + 136530185549551498916930460463168183800193567078764874592309425437632189282209601224600*alpha_d^2*cos(gamma)^3*sin(beta)^3*sin(gamma) - 538192581472218698160313212831614405611802689722670649671553673013420919406687319641936*alpha_d*beta_d*cos(beta)^2*cos(gamma)^2 - 596486859794902560808986557640857505303732151308467644260820271109756928655107697766936*alpha_d*beta_d*cos(beta)^2*sin(gamma)^2 - 321916218205724029689993516706867066695779581281883728135798307894142721593218806445968*alpha_d*beta_d*cos(gamma)^2*sin(beta)^2 + 136530185549551498916930460463168183800193567078764874592309425437632189282209601224600*alpha_d*beta_d*cos(gamma)^4*sin(beta)^2 - 321916218205724029689993516706867066695779581281883728135798307894142721593218806445968*alpha_d*beta_d*sin(beta)^2*sin(gamma)^2 - 136530185549551498916930460463168183800193567078764874592309425437632189282209601224600*alpha_d*beta_d*sin(beta)^2*sin(gamma)^4 + 5427557437643211752263139713717827000172949976593942171906872967234414616779136000000000*cos(beta)*cos(gamma)*sin(beta)*sin(gamma) - 136530185549551498916930460463168183800193567078764874592309425437632189282209601224600*beta_d^2*cos(gamma)*sin(beta)*sin(gamma)^3 - 136530185549551498916930460463168183800193567078764874592309425437632189282209601224600*beta_d^2*cos(gamma)^3*sin(beta)*sin(gamma) - 345893699911256339997794145344959090070795057460607577908389211039830077725186247884800*beta_d^2*cos(delta)^2*cos(gamma)*sin(beta)*sin(gamma)^3 - 345893699911256339997794145344959090070795057460607577908389211039830077725186247884800*beta_d^2*cos(delta)^2*cos(gamma)^3*sin(beta)*sin(gamma) - 218498690730410934784737927712893923642655650092914486436237191725678158211014302105600*beta_d^2*cos(delta)^4*cos(gamma)*sin(beta)*sin(gamma)^3 - 218498690730410934784737927712893923642655650092914486436237191725678158211014302105600*beta_d^2*cos(delta)^4*cos(gamma)^3*sin(beta)*sin(gamma) - 345893699911256339997794145344959090070795057460607577908389211039830077725186247884800*beta_d^2*cos(gamma)*sin(beta)*sin(delta)^2*sin(gamma)^3 - 345893699911256339997794145344959090070795057460607577908389211039830077725186247884800*beta_d^2*cos(gamma)^3*sin(beta)*sin(delta)^2*sin(gamma) - 218498690730410934784737927712893923642655650092914486436237191725678158211014302105600*beta_d^2*cos(gamma)*sin(beta)*sin(delta)^4*sin(gamma)^3 - 218498690730410934784737927712893923642655650092914486436237191725678158211014302105600*beta_d^2*cos(gamma)^3*sin(beta)*sin(delta)^4*sin(gamma) + 131249380619714751144491179690988118586621125026800009017820446014154152181025572093750*beta_d*delta_d*cos(beta)*cos(delta)^2*cos(gamma)^2 + 174798993960756348197390608542753849846273537113462224214112443769004354588114944000000*beta_d*delta_d*cos(beta)*cos(delta)^4*cos(gamma)^2 + 271941399723718215280434589139767797336610288733615494227168411281654523070827300093750*beta_d*gamma_d*cos(beta)*cos(delta)^2*cos(gamma)^2 + 174798993960756348197390608542753849846273537113462224214112443769004354588114944000000*beta_d*gamma_d*cos(beta)*cos(delta)^4*cos(gamma)^2 + 145465644810120481252162058223311942238581989596706857283884827451949003747230406078125*beta_d*delta_d*cos(beta)*cos(delta)^2*sin(gamma)^2 + 131249380619714751144491179690988118586621125026800009017820446014154152181025572093750*beta_d*delta_d*cos(beta)*cos(gamma)^2*sin(delta)^2 + 174798993960756348197390608542753849846273537113462224214112443769004354588114944000000*beta_d*delta_d*cos(beta)*cos(delta)^4*sin(gamma)^2 + 174798993960756348197390608542753849846273537113462224214112443769004354588114944000000*beta_d*delta_d*cos(beta)*cos(gamma)^2*sin(delta)^4 + 202358522656754687205081577811001572628472766711423925475374515772655267842228934078125*beta_d*gamma_d*cos(beta)*cos(delta)^2*sin(gamma)^2 + 271941399723718215280434589139767797336610288733615494227168411281654523070827300093750*beta_d*gamma_d*cos(beta)*cos(gamma)^2*sin(delta)^2 + 174798993960756348197390608542753849846273537113462224214112443769004354588114944000000*beta_d*gamma_d*cos(beta)*cos(delta)^4*sin(gamma)^2 + 174798993960756348197390608542753849846273537113462224214112443769004354588114944000000*beta_d*gamma_d*cos(beta)*cos(gamma)^2*sin(delta)^4 + 145465644810120481252162058223311942238581989596706857283884827451949003747230406078125*beta_d*delta_d*cos(beta)*sin(delta)^2*sin(gamma)^2 + 174798993960756348197390608542753849846273537113462224214112443769004354588114944000000*beta_d*delta_d*cos(beta)*sin(delta)^4*sin(gamma)^2 + 202358522656754687205081577811001572628472766711423925475374515772655267842228934078125*beta_d*gamma_d*cos(beta)*sin(delta)^2*sin(gamma)^2 + 174798993960756348197390608542753849846273537113462224214112443769004354588114944000000*beta_d*gamma_d*cos(beta)*sin(delta)^4*sin(gamma)^2 + 345893699911256339997794145344959090070795057460607577908389211039830077725186247884800*alpha_d^2*cos(delta)^2*cos(gamma)*sin(beta)^3*sin(gamma)^3 + 345893699911256339997794145344959090070795057460607577908389211039830077725186247884800*alpha_d^2*cos(delta)^2*cos(gamma)^3*sin(beta)^3*sin(gamma) + 218498690730410934784737927712893923642655650092914486436237191725678158211014302105600*alpha_d^2*cos(delta)^4*cos(gamma)*sin(beta)^3*sin(gamma)^3 + 218498690730410934784737927712893923642655650092914486436237191725678158211014302105600*alpha_d^2*cos(delta)^4*cos(gamma)^3*sin(beta)^3*sin(gamma) + 345893699911256339997794145344959090070795057460607577908389211039830077725186247884800*alpha_d^2*cos(gamma)*sin(beta)^3*sin(delta)^2*sin(gamma)^3 + 345893699911256339997794145344959090070795057460607577908389211039830077725186247884800*alpha_d^2*cos(gamma)^3*sin(beta)^3*sin(delta)^2*sin(gamma) + 218498690730410934784737927712893923642655650092914486436237191725678158211014302105600*alpha_d^2*cos(gamma)*sin(beta)^3*sin(delta)^4*sin(gamma)^3 + 218498690730410934784737927712893923642655650092914486436237191725678158211014302105600*alpha_d^2*cos(gamma)^3*sin(beta)^3*sin(delta)^4*sin(gamma) - 1359182022960050380851586254089165936210934896910692044241394114102654227548800940015818*alpha_d*beta_d*cos(beta)^2*cos(delta)^2*cos(gamma)^2 - 855570545488932589517493676539735847038203969899805669263484289173443648082264400592896*alpha_d*beta_d*cos(beta)^2*cos(delta)^4*cos(gamma)^2 - 1428764900027013908926939265417932160919072418932883612993188009611653482777399306031443*alpha_d*beta_d*cos(beta)^2*cos(delta)^2*sin(gamma)^2 - 1359182022960050380851586254089165936210934896910692044241394114102654227548800940015818*alpha_d*beta_d*cos(beta)^2*cos(gamma)^2*sin(delta)^2 - 815561711341884298066010421614466866773772592822153769234281262692154375309814120054784*alpha_d*beta_d*cos(delta)^2*cos(gamma)^2*sin(beta)^2 - 855570545488932589517493676539735847038203969899805669263484289173443648082264400592896*alpha_d*beta_d*cos(beta)^2*cos(delta)^4*sin(gamma)^2 - 855570545488932589517493676539735847038203969899805669263484289173443648082264400592896*alpha_d*beta_d*cos(beta)^2*cos(gamma)^2*sin(delta)^4 + 345893699911256339997794145344959090070795057460607577908389211039830077725186247884800*alpha_d*beta_d*cos(delta)^2*cos(gamma)^4*sin(beta)^2 - 515184769724844468857442142541244848442238753506633946738798366471224001335189672296448*alpha_d*beta_d*cos(delta)^4*cos(gamma)^2*sin(beta)^2 + 218498690730410934784737927712893923642655650092914486436237191725678158211014302105600*alpha_d*beta_d*cos(delta)^4*cos(gamma)^4*sin(beta)^2 - 1428764900027013908926939265417932160919072418932883612993188009611653482777399306031443*alpha_d*beta_d*cos(beta)^2*sin(delta)^2*sin(gamma)^2 - 815561711341884298066010421614466866773772592822153769234281262692154375309814120054784*alpha_d*beta_d*cos(delta)^2*sin(beta)^2*sin(gamma)^2 - 815561711341884298066010421614466866773772592822153769234281262692154375309814120054784*alpha_d*beta_d*cos(gamma)^2*sin(beta)^2*sin(delta)^2 - 855570545488932589517493676539735847038203969899805669263484289173443648082264400592896*alpha_d*beta_d*cos(beta)^2*sin(delta)^4*sin(gamma)^2 - 345893699911256339997794145344959090070795057460607577908389211039830077725186247884800*alpha_d*beta_d*cos(delta)^2*sin(beta)^2*sin(gamma)^4 - 515184769724844468857442142541244848442238753506633946738798366471224001335189672296448*alpha_d*beta_d*cos(delta)^4*sin(beta)^2*sin(gamma)^2 - 515184769724844468857442142541244848442238753506633946738798366471224001335189672296448*alpha_d*beta_d*cos(gamma)^2*sin(beta)^2*sin(delta)^4 + 345893699911256339997794145344959090070795057460607577908389211039830077725186247884800*alpha_d*beta_d*cos(gamma)^4*sin(beta)^2*sin(delta)^2 - 218498690730410934784737927712893923642655650092914486436237191725678158211014302105600*alpha_d*beta_d*cos(delta)^4*sin(beta)^2*sin(gamma)^4 + 218498690730410934784737927712893923642655650092914486436237191725678158211014302105600*alpha_d*beta_d*cos(gamma)^4*sin(beta)^2*sin(delta)^4 + 58294278322683862648673344809243099691929461585796994589266598096336009248420378125000*alpha_d^2*cos(beta)^2*cos(gamma)*sin(beta)*sin(gamma) + 7745983541128480209019527002681234934260727617113545615320285912254010464665600000000000*cos(beta)*cos(delta)^2*cos(gamma)*sin(beta)*sin(gamma) - 815561711341884298066010421614466866773772592822153769234281262692154375309814120054784*alpha_d*beta_d*sin(beta)^2*sin(delta)^2*sin(gamma)^2 - 345893699911256339997794145344959090070795057460607577908389211039830077725186247884800*alpha_d*beta_d*sin(beta)^2*sin(delta)^2*sin(gamma)^4 - 515184769724844468857442142541244848442238753506633946738798366471224001335189672296448*alpha_d*beta_d*sin(beta)^2*sin(delta)^4*sin(gamma)^2 - 218498690730410934784737927712893923642655650092914486436237191725678158211014302105600*alpha_d*beta_d*sin(beta)^2*sin(delta)^4*sin(gamma)^4 + 7745983541128480209019527002681234934260727617113545615320285912254010464665600000000000*cos(beta)*cos(gamma)*sin(beta)*sin(delta)^2*sin(gamma) - 436997381460821869569475855425787847285311300185828972872474383451356316422028604211200*beta_d^2*cos(delta)^2*cos(gamma)*sin(beta)*sin(delta)^2*sin(gamma)^3 - 436997381460821869569475855425787847285311300185828972872474383451356316422028604211200*beta_d^2*cos(delta)^2*cos(gamma)^3*sin(beta)*sin(delta)^2*sin(gamma) + 58294278322683862648673344809243099691929461585796994589266598096336009248420378125000*alpha_d*gamma_d*cos(beta)*cos(gamma)*sin(beta)*sin(gamma) + 349597987921512696394781217085507699692547074226924448428224887538008709176229888000000*beta_d*delta_d*cos(beta)*cos(delta)^2*cos(gamma)^2*sin(delta)^2 + 349597987921512696394781217085507699692547074226924448428224887538008709176229888000000*beta_d*gamma_d*cos(beta)*cos(delta)^2*cos(gamma)^2*sin(delta)^2 + 349597987921512696394781217085507699692547074226924448428224887538008709176229888000000*beta_d*delta_d*cos(beta)*cos(delta)^2*sin(delta)^2*sin(gamma)^2 + 349597987921512696394781217085507699692547074226924448428224887538008709176229888000000*beta_d*gamma_d*cos(beta)*cos(delta)^2*sin(delta)^2*sin(gamma)^2 + 436997381460821869569475855425787847285311300185828972872474383451356316422028604211200*alpha_d^2*cos(delta)^2*cos(gamma)*sin(beta)^3*sin(delta)^2*sin(gamma)^3 + 436997381460821869569475855425787847285311300185828972872474383451356316422028604211200*alpha_d^2*cos(delta)^2*cos(gamma)^3*sin(beta)^3*sin(delta)^2*sin(gamma) - 1711141090977865179034987353079471694076407939799611338526968578346887296164528801185792*alpha_d*beta_d*cos(beta)^2*cos(delta)^2*cos(gamma)^2*sin(delta)^2 - 1711141090977865179034987353079471694076407939799611338526968578346887296164528801185792*alpha_d*beta_d*cos(beta)^2*cos(delta)^2*sin(delta)^2*sin(gamma)^2 - 1030369539449688937714884285082489696884477507013267893477596732942448002670379344592896*alpha_d*beta_d*cos(delta)^2*cos(gamma)^2*sin(beta)^2*sin(delta)^2 + 436997381460821869569475855425787847285311300185828972872474383451356316422028604211200*alpha_d*beta_d*cos(delta)^2*cos(gamma)^4*sin(beta)^2*sin(delta)^2 + 69582877066963528075353011328766224708137522022191568751793895508999255228598366015625*alpha_d^2*cos(beta)^2*cos(delta)^2*cos(gamma)*sin(beta)*sin(gamma) - 1030369539449688937714884285082489696884477507013267893477596732942448002670379344592896*alpha_d*beta_d*cos(delta)^2*sin(beta)^2*sin(delta)^2*sin(gamma)^2 - 436997381460821869569475855425787847285311300185828972872474383451356316422028604211200*alpha_d*beta_d*cos(delta)^2*sin(beta)^2*sin(delta)^2*sin(gamma)^4 + 69582877066963528075353011328766224708137522022191568751793895508999255228598366015625*alpha_d^2*cos(beta)^2*cos(gamma)*sin(beta)*sin(delta)^2*sin(gamma) - 14216264190405730107670878532323823651960864569906848266064381437794851566204833984375*alpha_d*delta_d*cos(beta)*cos(delta)^2*cos(gamma)*sin(beta)*sin(gamma) + 69582877066963528075353011328766224708137522022191568751793895508999255228598366015625*alpha_d*gamma_d*cos(beta)*cos(delta)^2*cos(gamma)*sin(beta)*sin(gamma) - 14216264190405730107670878532323823651960864569906848266064381437794851566204833984375*alpha_d*delta_d*cos(beta)*cos(gamma)*sin(beta)*sin(delta)^2*sin(gamma) + 69582877066963528075353011328766224708137522022191568751793895508999255228598366015625*alpha_d*gamma_d*cos(beta)*cos(gamma)*sin(beta)*sin(delta)^2*sin(gamma))/(22116597823247856*(107921345109286508779897069447938048*cos(delta)^2 + 107921345109286508779897069447938048*sin(delta)^2 + 81033702656323036153327122834079417)*(179621720919313878356379449073236959*cos(gamma)^2*sin(beta) + 179621720919313878356379449073236959*sin(beta)*sin(gamma)^2 + 215842690218573017559794138895876096*cos(delta)^2*cos(gamma)^2*sin(beta) + 215842690218573017559794138895876096*cos(delta)^2*sin(beta)*sin(gamma)^2 + 215842690218573017559794138895876096*cos(gamma)^2*sin(beta)*sin(delta)^2 + 215842690218573017559794138895876096*sin(beta)*sin(delta)^2*sin(gamma)^2))];

    
    xdot = [x(5);
            x(6);
            x(7);
            x(8);
            Xdd]; 

end