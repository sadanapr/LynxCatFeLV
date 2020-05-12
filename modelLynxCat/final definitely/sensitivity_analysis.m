function sensitivity_analysis()

% %BETA_CL
%     %Parameters
%     param.L0 = 30; %OK 30
%     param.S_L0 = param.L0;
%     param.I_L0 = 0;
%     init_infect = 0.03; %initial infected in Cats
%     param.C0 = 500; %ESTIMATION 500
%     param.S_C0 = (1 - init_infect)*param.C0;
%     param.I_C0 = init_infect * param.C0;
%     param.R_C0 = 0;
% 
%     param.R0 = 800;
%     
%     param.e_LR = 0.127; %OK  0.12
%     param.e_CR = 0.24; %ESTIMATION 0.24
%     param.e_LC = 0.0001; %ESTIMATION 0.01
%     param.beta_L = (26/250000) * 25 * 0.9; %ESTIMATION 0.07 
%     param.beta_C = (500/250000) * 25 * 0.3; %OK 0.027;
%     param.beta_CL = 0.95 * (26/250000) * 25 * 0.7; %ESTIMATION 
%     param.gamma_L = 0.00011; %OK 0.00046
%     param.rho_C = 0.0018; %OK 
%     param.delta_L = 0.0028; %OK 0.0028
%     param.delta_C = 0.0007; %OK 0.0007
%     param.q_L = 0.1; %ESTIMATION 0.036
%     param.p_C = 0.5; %OK 
%     param.p_L = 0.5; %OK 
%     param.x_L = 0.12; %OK 0.14
%     param.x_C = 0.08; %BOF 0.002; 0.1
%     param.y_R = 6.1; %OK 6.1
%     param.K_C = 600;
%     param.K_R = 800;
%     
%     range = linspace (0.0005,0.008,300);
%     S = [];
%     for r = range
%         param.beta_CL = r;
%         %Model
%         X0 = [param.S_L0 param.I_L0 param.S_C0 param.I_C0 param.R_C0 param.R0];
%         % simulation parameters
%             tf = 365 * 5; %Time
%             options = odeset('NonNegative',1:6);
%         [t,X] = ode45( @(t,x) Mod_rhs(t,x,param), [0,tf], X0(:), options);
% 
%         X = X';
%         t = t';
%         i = 1;
%         s = t(1);
%         keep = true;
%         while (s ~= tf) && (keep)
%             s = t(i);
%             if (X(1,i) + X(2,i) < 15)
%                 keep = false;
%             end
%             i = i + 1;
%         end
%         S = [S s];
%     end
%     p = polyfit(range,S, 3)
%     figure
%     plot(range, S,'b')
%     hold on;
%     plot(range, polyval(p,range),'r')
%     hold off;
%     xlabel('beta_CL');
%     ylabel('t_extinction')
%     
%    
%         param.beta_CL = 0.95 * (26/250000) * 25 * 0.7;
%         %Model
%         X0 = [param.S_L0 param.I_L0 param.S_C0 param.I_C0 param.R_C0 param.R0];
%         % simulation parameters
%             tf = 365 * 5; %Time
%             options = odeset('NonNegative',1:6);
%         [t,X] = ode45( @(t,x) Mod_rhs(t,x,param), [0,tf], X0(:), options);
% 
%         X = X';
%         t = t';
%         i = 1;
%         s = t(1);
%         keep = true;
%         while (s ~= tf) && (keep)
%             s = t(i);
%             if (X(1,i) + X(2,i) < 15)
%                 keep = false;
%             end
%             i = i + 1;
%         end
%     extinction = s
%     k = polyder(p);
%     sensitivity = polyval(k,param.beta_CL)*param.beta_CL/s
    
% %BETA_CL (centered thing)
%Parameters
    param.L0 = 30; %OK 30
    param.S_L0 = param.L0;
    param.I_L0 = 0;
    init_infect = 0.03; %initial infected in Cats
    param.C0 = 500; %ESTIMATION 500 = 470 + 30
    param.S_C0 = (1 - init_infect)*param.C0;
    param.I_C0 = init_infect * param.C0;
    param.R_C0 = 0;

    param.R0 = 800;
    
    param.e_LR = 0.127; %OK  0.12
    param.e_CR = 0.24; %ESTIMATION 0.24
    param.e_LC = 0.0001; %ESTIMATION 0.01
    param.beta_L = (30/250000) * 25 * 0.9; %ESTIMATION 0.07 
    param.beta_C = (500/250000) * 25 * 0.3; %OK 0.027;
    param.beta_CL = 0.95 * (30/250000) * 25 * 0.7; %ESTIMATION 
    param.gamma_L = 0.00011; %OK 0.00046
    param.rho_C = 0.00046; %OK 
    param.delta_L = 0.0032; %OK 0.0032
    param.delta_C = 0.0007; %OK 0.0007
    param.q_L = 0.1; %ESTIMATION 0.036
    param.p_C = 0.5; %OK 
    param.p_L = 0.5; %OK 
    param.x_L = 0.12; %OK 0.14
    param.x_C = 0.08; %BOF 0.002; 0.1
    param.y_R = 6.1; %OK 6.1
    param.K_C = 600;
    param.K_R = 2000;
    param.K_L = 800;
    
    
    h = 0.01 * param.beta_CL;
    param.beta_CL = param.beta_CL + h;
    X0 = [param.S_L0 param.I_L0 param.S_C0 param.I_C0 param.R_C0 param.R0];
    % simulation parameters
        tf = 365 * 15; %Time
        options = odeset('NonNegative',1:6);
    [t,X] = ode45( @(t,x) Mod_rhs(t,x,param), [0,tf], X0(:), options);
    X = X';
    t = t';
    i = 1;
    s1 = t(1);
    keep = true;
    while (s1 ~= tf) && (keep)
        s1 = t(i);
        if (X(1,i) + X(2,i) < 8)
            keep = false;
        end
        i = i + 1;
    end
    
    param.beta_CL = param.beta_CL - 2*h;
    X0 = [param.S_L0 param.I_L0 param.S_C0 param.I_C0 param.R_C0 param.R0];
    % simulation parameters
        tf = 365 * 15; %Time
        options = odeset('NonNegative',1:6);
    [t,X] = ode45( @(t,x) Mod_rhs(t,x,param), [0,tf], X0(:), options);
    X = X';
    t = t';
    i = 1;
    s2 = t(1);
    keep = true;
    while (s2 ~= tf) && (keep)
        s2 = t(i);
        if (X(1,i) + X(2,i) < 8)
            keep = false;
        end
        i = i + 1;
    end
    
    param.beta_CL = param.beta_CL + h;
    X0 = [param.S_L0 param.I_L0 param.S_C0 param.I_C0 param.R_C0 param.R0];
    % simulation parameters
        tf = 365 * 15; %Time
        options = odeset('NonNegative',1:6);
    [t,X] = ode45( @(t,x) Mod_rhs(t,x,param), [0,tf], X0(:), options);
    X = X';
    t = t';
    i = 1;
    s = t(1);
    keep = true;
    while (s ~= tf) && (keep)
        s = t(i);
        if (X(1,i) + X(2,i) < 8)
            keep = false;
        end
        i = i + 1;
    end
    sensitivityCL = ((s1 - s2)/(2*h)) * param.beta_CL/s

%BETA_L (centered thing)
    %Parameters
    param.L0 = 30; %OK 30
    param.S_L0 = param.L0;
    param.I_L0 = 0;
    init_infect = 0.03; %initial infected in Cats
    param.C0 = 500; %ESTIMATION 500 = 470 + 30
    param.S_C0 = (1 - init_infect)*param.C0;
    param.I_C0 = init_infect * param.C0;
    param.R_C0 = 0;

    param.R0 = 800;
    
    param.e_LR = 0.127; %OK  0.12
    param.e_CR = 0.24; %ESTIMATION 0.24
    param.e_LC = 0.0001; %ESTIMATION 0.01
    param.beta_L = (30/250000) * 25 * 0.9; %ESTIMATION 0.07 
    param.beta_C = (500/250000) * 25 * 0.3; %OK 0.027;
    param.beta_CL = 0.95 * (30/250000) * 25 * 0.7; %ESTIMATION 
    param.gamma_L = 0.00011; %OK 0.00046
    param.rho_C = 0.00046; %OK 
    param.delta_L = 0.0032; %OK 0.0032
    param.delta_C = 0.0007; %OK 0.0007
    param.q_L = 0.1; %ESTIMATION 0.036
    param.p_C = 0.5; %OK 
    param.p_L = 0.5; %OK 
    param.x_L = 0.12; %OK 0.14
    param.x_C = 0.08; %BOF 0.002; 0.1
    param.y_R = 6.1; %OK 6.1
    param.K_C = 600;
    param.K_R = 2000;
    param.K_L = 800;
    
    
    
    h = 0.1 * param.beta_L;
    param.beta_L = param.beta_L + h;
    X0 = [param.S_L0 param.I_L0 param.S_C0 param.I_C0 param.R_C0 param.R0];
    % simulation parameters
        tf = 365 * 15; %Time
        options = odeset('NonNegative',1:6);
    [t,X] = ode45( @(t,x) Mod_rhs(t,x,param), [0,tf], X0(:), options);
    X = X';
    t = t';
    i = 1;
    s1 = t(1);
    keep = true;
    while (s1 ~= tf) && (keep)
        s1 = t(i);
        if (X(1,i) + X(2,i) < 8)
            keep = false;
        end
        i = i + 1;
    end
    
    param.beta_L = param.beta_L - 2*h;
    X0 = [param.S_L0 param.I_L0 param.S_C0 param.I_C0 param.R_C0 param.R0];
    % simulation parameters
        tf = 365 * 15; %Time
        options = odeset('NonNegative',1:6);
    [t,X] = ode45( @(t,x) Mod_rhs(t,x,param), [0,tf], X0(:), options);
    X = X';
    t = t';
    i = 1;
    s2 = t(1);
    keep = true;
    while (s2 ~= tf) && (keep)
        s2 = t(i);
        if (X(1,i) + X(2,i) < 8)
            keep = false;
        end
        i = i + 1;
    end
    
    param.beta_L = param.beta_L + h;
    X0 = [param.S_L0 param.I_L0 param.S_C0 param.I_C0 param.R_C0 param.R0];
    % simulation parameters
        tf = 365 * 15; %Time
        options = odeset('NonNegative',1:6);
    [t,X] = ode45( @(t,x) Mod_rhs(t,x,param), [0,tf], X0(:), options);
    X = X';
    t = t';
    i = 1;
    s = t(1);
    keep = true;
    while (s ~= tf) && (keep)
        s = t(i);
        if (X(1,i) + X(2,i) < 8)
            keep = false;
        end
        i = i + 1;
    end
    sensitivityL = ((s1 - s2)/(2*h)) * param.beta_L/s
    
%BETA_L
 %Parameters
%     param.L0 = 30; %OK 30
%     param.S_L0 = param.L0;
%     param.I_L0 = 0;
%     init_infect = 0.03; %initial infected in Cats
%     param.C0 = 500; %ESTIMATION 500
%     param.S_C0 = (1 - init_infect)*param.C0;
%     param.I_C0 = init_infect * param.C0;
%     param.R_C0 = 0;
% 
%     param.R0 = 800;
%     
%     param.e_LR = 0.127; %OK  0.12
%     param.e_CR = 0.24; %ESTIMATION 0.24
%     param.e_LC = 0.0001; %ESTIMATION 0.01
%     param.beta_L = (26/250000) * 25 * 0.9; %ESTIMATION 0.07 
%     param.beta_C = (500/250000) * 25 * 0.3; %OK 0.027;
%     param.beta_CL = 0.95 * (26/250000) * 25 * 0.7; %ESTIMATION 
%     param.gamma_L = 0.00011; %OK 0.00046
%     param.rho_C = 0.0018; %OK 
%     param.delta_L = 0.0028; %OK 0.0028
%     param.delta_C = 0.0007; %OK 0.0007
%     param.q_L = 0.1; %ESTIMATION 0.036
%     param.p_C = 0.5; %OK 
%     param.p_L = 0.5; %OK 
%     param.x_L = 0.12; %OK 0.14
%     param.x_C = 0.08; %BOF 0.002; 0.1
%     param.y_R = 6.1; %OK 6.1
%     param.K_C = 600;
%     param.K_R = 800;
%     
%     range = linspace (0.001,0.01,100);
%     S = [];
%     for r = range
%         param.beta_L = r;
%         %Model
%         X0 = [param.S_L0 param.I_L0 param.S_C0 param.I_C0 param.R_C0 param.R0];
%         % simulation parameters
%             tf = 365 * 5; %Time
%             options = odeset('NonNegative',1:6);
%         [t,X] = ode45( @(t,x) Mod_rhs(t,x,param), [0,tf], X0(:), options);
% 
%         X = X';
%         t = t';
%         i = 1;
%         s = t(1);
%         keep = true;
%         while (s ~= tf) && (keep)
%             s = t(i);
%             if (X(1,i) + X(2,i) < 2)
%                 keep = false;
%             end
%             i = i + 1;
%         end
%         S = [S s];
%     end
%     p = polyfit(range,S, 3)
%     figure
%     plot(range, S,'b')
%     hold on;
%     plot(range, polyval(p,range),'r')
%     hold off;
%     xlabel('beta_L');
%     ylabel('t_extinction')
%     
%    
%         
%         %Model
%         param.beta_L = (26/250000) * 25 * 0.9; %ESTIMATION 0.07 
%         X0 = [param.S_L0 param.I_L0 param.S_C0 param.I_C0 param.R_C0 param.R0];
%         % simulation parameters
%             tf = 365 * 5; %Time
%             options = odeset('NonNegative',1:6);
%         [t,X] = ode45( @(t,x) Mod_rhs(t,x,param), [0,tf], X0(:), options);
% 
%         X = X';
%         t = t';
%         i = 1;
%         s = t(1);
%         keep = true;
%         while (s ~= tf) && (keep)
%             s = t(i);
%             if (X(1,i) + X(2,i) < 2)
%                 keep = false;
%             end
%             i = i + 1;
%         end
%     extinction = s
%     k = polyder(p);
%     sensitivity = polyval(k,param.beta_CL)*param.beta_L/s
    
end    

%FUNCTION (mix of SIR and Rosenzweig-MacArthur model)
function dxdt = Mod_rhs(t,x,param)
% function that computes the right-hand side of the differential equation
    S_L = x(1);
    I_L = x(2);
    L = S_L + I_L;
    S_C = x(3);
    I_C = x(4);
    R_C = x(5);
    C = S_C + I_C + R_C;
    R = x(6);
    
    dxdt = [  (param.e_LC * C * S_L/(C + 60)*(1 - L/param.K_L)) +  (param.e_LR * R * S_L/(R + 48)*(1 - L/param.K_L)) + (param.gamma_L * I_L) - (param.x_L * S_L) - (param.beta_L * I_L * S_L/L) - (param.beta_CL * I_C * S_L/L), ... dS_L/dt
              (param.e_LC * C * I_L/(C + 60)*(1 - L/param.K_L)) + (param.e_LR * R * I_L/(R + 48)*(1 - L/param.K_L)) + (param.beta_L * I_L * S_L/L) + (param.beta_CL * I_C * S_L/L) - (param.x_L * I_L) - (param.gamma_L * I_L) - (param.delta_L * I_L), ... dI_L/dt
              (param.e_CR * R * S_C/(R + 48) * (1 - C/param.K_C)) - (param.x_C * S_C) - (param.beta_C * I_C * S_C/C) - (param.q_L * L * S_C/(C + 60))  , ... dS_C/dt
              (param.e_CR * R * I_C/(R + 48) * (1 - C/param.K_C)) + (param.beta_C * I_C * S_C/C) - (param.x_C * I_C) - (param.q_L * L * I_C/(C + 60)) - (param.rho_C * I_C) - (param.delta_C * I_C)  , ... dI_C/dt
              (param.e_CR * R * R_C/(R + 48) * (1 - C/param.K_C)) + (param.rho_C * I_C) - (param.x_C * R_C) - (param.q_L * L * R_C/(C + 60))  , ... dR_C/dt
              (param.y_R * R*(1 - R/param.K_R)) - (param.p_C * R * C/(R + 48)) - (param.p_L * R * L/(R + 48))    ... dR/dt
           ];
	dxdt = dxdt (:);
end