% Heat Transfer - Spring 2023 Computational Assignment 1
% Michael Allen, Cullen Hirstius, Hayden Payne
% 2/8/2022

clc; clear;  close all;
format short 

for WINDOW_HEIGHT = 1:3:4 %loops twice, once for .5m^2 and once for 2m^2 (math for window area shown on line 20)
    %constants
    k_glass = 1.4;              %W/(m*k)
    t_glass = 0.007;            %m
    k_air = .0245;              %W/(m*k)
    t_air = 0.007;              %m
    
    T_inside = 21+273;          %°K
    T_outside = -2+273;         %°K
    
    %WINDOW_HEIGHT = 1;          %m
    WINDOW_WIDTH = 2;           %m
    windowArea = WINDOW_WIDTH * WINDOW_HEIGHT / 4;
    
    h_i = 10;                   %W/(m^2*K)
    h_0 = [10:10:100];          %W/(m^2*K) range
    EMISSIVITY = 0.05;  
    ABSORPTIVITY = 0.2; 
    SOLAR_IRRADIATION = 1094.4; %W/m^2
    
    
    
    PANE_PRICE_USD = 223; %$223 per pane installed on window
    HEAT_LOSS_COST = .12; % USD/kWh
    WINTER_DAYS = 30 * 4; % 4 months * 30 days/mo
    YEARS = [1:10];         
    OPERATING_HOURS = YEARS * WINTER_DAYS * 24;
    watts_to_kWh = @(q) q/(1000) * OPERATING_HOURS; %convert W to kW/h
    
    
    
    %Part 1 (testing math)
    
    %LHS
    R_convInside = 1/(h_i*windowArea);
    R_condGlassPane = t_glass / (k_glass*windowArea);
    q_lhs = @(T_int) (T_inside-T_int)/(R_convInside + R_condGlassPane); %W
    
    %RHS
    p1_h0 = h_0(1);
    sigma = 5.67*10^(-8);
    Irradiation_OUTSIDE = ABSORPTIVITY*SOLAR_IRRADIATION*windowArea; 
    R_convOutside = 1/(p1_h0*windowArea);
    q_rhs = @(T_int) -(T_int - T_outside)/(R_convOutside) + (Irradiation_OUTSIDE) - (sigma*EMISSIVITY)*(T_int^4)*windowArea; %W
    
    %Solution p1
    Q = @(T_int) q_lhs(T_int) + q_rhs(T_int) ; %W = 0
    T_int = fsolve(Q, 0);
    %disp(T_int)            %Test case soln
    %disp(q_lhs(T_int))     %Test case soln q
    
    
    %Given N panes and Varying h (nested for loops)
    
    %Constant resistance values across variation
    Irradiation_OUTSIDE = ABSORPTIVITY*SOLAR_IRRADIATION*windowArea; 
    R_convInside = 1/(h_i*windowArea);
    R_condAir = t_air / (k_air*windowArea);
    R_condGlassPane = t_glass / (k_glass*windowArea);
    
    
    %Initalize results vector, which stores definitvie metricies such as T_int,
    %q_total, q_radiation, and q_conv
    resultsVector = zeros(3, length(h_0), 4);
    
    for N = 1:3
        %Calculate resistances of glass and air in series, depending on the
        %number of panes (1-3)
        R_glass = N*R_condGlassPane;
        R_air = (N-1)*R_condAir;
    
    
        for h_outside = 1:length(h_0)
            %Calculate convection resistance given h_outside
            R_convOutside = 1/(h_0(h_outside)*windowArea);
            
            %Define q equations
            q_lhs = @(T_int) (T_inside-T_int)/(R_convInside + R_glass + R_air); %W
            q_rhs = @(T_int) -(T_int - T_outside)/(R_convOutside) + (Irradiation_OUTSIDE) - (sigma*EMISSIVITY)*(T_int^4)*windowArea; %W
            Q = @(T_int) q_lhs(T_int) + q_rhs(T_int) ; %W = 0
    
            %Results
            T_int = fsolve(Q, 0);
            q_total = q_lhs(T_int);
            q_convOutside = @(T_int) -(T_int - T_outside)/R_convOutside;
            q_radLoss = @(T_int) (sigma*EMISSIVITY)*(T_int^4)*windowArea;
           
            %Store results in matrix
            resultsVector(N, h_outside,1) = T_int;
            resultsVector(N, h_outside,2) = q_total;
            resultsVector(N, h_outside,3) = q_convOutside(T_int);
            resultsVector(N, h_outside,4) = q_radLoss(T_int);
        end
    end
    
    %disp(resultsVector)
    
    
    %Cost Analysis
    heatLost = resultsVector(:,4,2);                         %W
    heatLost_cost = watts_to_kWh(heatLost) * HEAT_LOSS_COST; %USD formatted as (NumPanes, Year)
    overallCost = heatLost_cost;                            
    
    %add window pane price to each respective row
    for i = 1:3
        for j = 1:10
            overallCost(i,j) = heatLost_cost(i,j) + PANE_PRICE_USD*i;
        end
    end
    
    figure('Name', 'Cost over 10 years')  %initalize first figure
    hold on
    title(sprintf('Cost of Heat Loss for %.1f m^2 area', windowArea));
    xlabel('Num Years')
    ylabel('Cost ($USD)')
    plot(YEARS, overallCost);
    legend('1 Pane', '2 Panes', '3 Panes');
    hold off
    
    
    %Heat Transfer Rate (q) / convection coefficient
    heatTransferRate = resultsVector(:,:,2);
    
    figure('Name', 'Heat Transfer Rate (q) / Convection Coefficient (h)')  %initalize first figure
    hold on
    title(sprintf('Heat Transfer Rate (q) / Convection Coefficient (h) for %.1f m^2 area', windowArea));
    xlabel('Convection Coeff. (W/m^2*K)');
    ylabel('Rate of Heat Transfer (W)');
    plot(h_0, heatTransferRate);
    legend('1 Pane', '2 Panes', '3 Panes');
    hold off
    
    
    %Interface Temp / convection Coefficient 
    interfaceTemp = resultsVector(:,:,1);
    
    figure('Name', 'Interface Temp (K) / Convection Coefficient (h)')  %initalize first figure
    hold on
    title(sprintf('Interface Temp (K)  / Convection Coefficient (h) for %.1f m^2 area', windowArea));
    xlabel('Convection Coefficient (h)');
    ylabel('Interface Temp (K)');
    plot(h_0, interfaceTemp);
    legend('1 Pane', '2 Panes', '3 Panes');
    hold off
    
    
    
    %Heat Transfer / Num panes
    heatTransferRate = resultsVector(:,:,2);
    
    figure('Name', 'Heat Transfer Rate / Num Panes')  %initalize first figure
    hold on
    title(sprintf('Heat Transfer Rate / Num Panes for %.1f m^2 area', windowArea));
    xlabel('Num Panes');
    xticks([1 2 3]);
    ylabel('Heat Transfer Rate (q)');
    bar([1:3], heatTransferRate);
    legendStrings = "h = " + string(h_0);
    legend(legendStrings);
    hold off
    
    
    %Convection / Radiation Cost
    % ratio related to each entry in results 
    q_convOutside = resultsVector(:,:,3);
    q_radLoss = resultsVector(:,:,4);
    
    convection_to_radiation_ratio = zeros(3, length(h_0), 1);
    
    for N = 1:3
        for h_outside = 1:length(h_0)
            convection_to_radiation_ratio(N, h_outside,1) = abs(q_convOutside(N, h_outside) / q_radLoss(N, h_outside));
        end
    end
    
    figure('Name', 'Ratio of Convective (q_h) to Radiation (q_r) Losses')  %initalize first figure
    hold on
    bar3([1 2 3], convection_to_radiation_ratio);
    xlabel('Convection Coefficient (*10 W/m^2*K)');
    yticks([1 2 3]);
    ylabel('Num Panes');
    zlabel('Ratio (q_h / q_r)')
    view([225 40]);
    title(sprintf('Ratio of Convective (q_h) to Radiation (q_r) Losses for %.1f m^2 area', windowArea));
    %s = pcolor(h_0, [1:3], convection_to_radiation_ratio);
    hold off
end









