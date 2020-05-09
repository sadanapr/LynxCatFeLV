function isolation()
    %Parameters
    param.L0 = 30; %OK 30
    param.S_L0 = param.L0;
    param.IS_L0 = 0;
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
    param.K_C = 800;
    param.K_R = 2000;
    param.K_L = 300;
    param.isol_L = 0.011; %isolated infected lynx each day 0.011
    param.cure_L = 0.001; %curing rate in isolation (0.001 normal, 0.1 with great cure)
    
    %Model
    X0 = [param.S_L0 param.I_L0 param.S_C0 param.I_C0 param.R_C0 param.R0 param.IS_L0];
    % simulation parameters
        tf = 365 * 13; %Time
        options = odeset('NonNegative',1:7);
    [t,X] = ode45( @(t,x) Mod_rhs(t,x,param), [0,tf], X0(:), options);

    X = X';
    t = t';
    
    %plot all
    figure
    plot(t,X(1,:),'k','DisplayName','SL');
    hold on;
    plot(t,X(2,:),'r','DisplayName','IL');
    plot(t,X(3,:),'b','DisplayName','SC');
    plot(t,X(4,:),'c','DisplayName','IC');
    plot(t,X(5,:),'m','DisplayName','RC');
    plot(t,X(6,:),'g','DisplayName','R');
    plot(t,X(7,:),'g','DisplayName','QL');    
    hold off;
    legend show;
    
    %plot R
    figure
    plot(t,X(6,:),'g','DisplayName','R');
    legend show;
    
    %plot C 
    figure
    plot(t,X(3,:),'k','DisplayName','SC');
    hold on;
    plot(t,X(4,:),'r','DisplayName','IC');
    plot(t,X(5,:),'m','DisplayName','RC');
    plot(t,X(3,:) + X(4,:) + X(5,:),'b','DisplayName','C');
    hold off;
    legend show;
    
    %plot L
    figure
    plot(t,X(1,:),'k','DisplayName','Susceptible Lynx');
    hold on;
    plot(t,X(2,:),'r','DisplayName','Infected Lynx');
    plot(t,X(7,:),'g','DisplayName','Quarantined Lynx');    
    plot(t,X(1,:) + X(2,:) + X(7,:),'b','DisplayName','Total Lynx');
    xlim([0 365*5])
    xlabel('Time (days)')
    ylabel('Population')
    hold off;
    legend show;
    
    %plot LCR
    figure
    plot(t,X(1,:) + X(2,:),'r','DisplayName','L');
    hold on;
    plot(t,X(3,:)+X(4,:)+X(5,:),'b','DisplayName','C');
    plot(t,X(6,:),'g','DisplayName','R');
    xlim([0 365*13])
    xlabel('Time (days)')
    ylabel('Population')
    hold off;
    legend show;
    
    %plot 3D
    figure
    scatter3(X(1,:) + X(2,:) + X(7,:),X(3,:)+X(4,:)+X(5,:),X(6,:));
    xlabel('L')
    ylabel('C')
    zlabel('R')
end    

% FUNCTION (mix of SIR and Rosenzweig-MacArthur model)
function dxdt = Mod_rhs(t,x,param)
% function that computes the right-hand side of the differential equation
    S_L = x(1);
    I_L = x(2);
    S_C = x(3);
    I_C = x(4);
    R_C = x(5);
    C = S_C + I_C + R_C;
    R = x(6);
    IS_L = x(7);
    L = S_L + I_L + IS_L;
   dxdt = [  (param.e_LC * C * S_L/(C + 60)*(1 - L/param.K_L)) +  (param.e_LR * R * S_L/(R + 48)*(1 - L/param.K_L)) + (param.gamma_L * I_L) + (param.gamma_L * IS_L) - (param.x_L * S_L) - (param.beta_L * I_L * S_L/L) - (param.beta_CL * I_C * S_L/L), ... dS_L/dt
              (param.e_LC * C * I_L/(C + 60)*(1 - L/param.K_L)) + (param.e_LR * R * I_L/(R + 48)*(1 - L/param.K_L)) + (param.beta_L * I_L  *  S_L/L) + (param.beta_CL * I_C * S_L/L) - (param.x_L * I_L) - (param.gamma_L * I_L) - (param.delta_L * I_L) - (param.isol_L* IS_L), ... dI_L/dt
              (param.e_CR * R * S_C/(R + 48) * (1 - C/param.K_C)) - (param.x_C * S_C) - (param.beta_C * I_C * S_C/C) - (param.q_L * L * S_C/(C + 60))  , ... dS_C/dt
              (param.e_CR * R * I_C/(R + 48) * (1 - C/param.K_C)) + (param.beta_C * I_C * S_C/C) - (param.x_C * I_C) - (param.q_L * L * I_C/(C + 60)) - (param.rho_C * I_C) - (param.delta_C * I_C) - param.isol_L*I_C , ... dI_C/dt
              (param.e_CR * R * R_C/(R + 48) * (1 - C/param.K_C)) + (param.rho_C * I_C) - (param.x_C * R_C) - (param.q_L * L * R_C/(C + 60))  , ... dR_C/dt
              (param.y_R * R*(1 - R/param.K_R)) - (param.p_C * R * C/(R + 48)) - (param.p_L * R * L/(R + 48))  ,  ... dR/dt
              (param.isol_L*I_L) - (param.cure_L * IS_L) - (param.delta_L * IS_L)%quarantined
           ];
	dxdt = dxdt (:);
end