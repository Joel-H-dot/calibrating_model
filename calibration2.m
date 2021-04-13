function [] = calibration()



for kk=1:3
extension = ['_data_' num2str(kk) '.xlsx'];
table_95K = xlsread(['C:\Users\EEE Admin\Desktop\measured_data\95k' extension ]);
frequencies(1,:) = table_95K(:,2);
real_95K = table_95K(:,3);
imag_95K = table_95K(:,4);

table_45K = xlsread(['C:\Users\EEE Admin\Desktop\measured_data\45k' extension]);
real_45K = table_45K(:,3);
imag_45K = table_45K(:,4);

table_10K = xlsread(['C:\Users\EEE Admin\Desktop\measured_data\10k' extension]);
real_10K = table_10K(:,3);
imag_10K = table_10K(:,4);

table_test = xlsread(['C:\Users\EEE Admin\Desktop\measured_data\test' extension]);
real_test = table_test(:,3);
imag_test = table_test(:,4);

table_air = xlsread(['C:\Users\EEE Admin\Desktop\measured_data\air' extension]);
real_air  = table_air (:,3);
imag_air  = table_air (:,4);

inductance_95k_meas_all(kk,:)  =  (real_95K+1i*imag_95K)./(1i*2*pi*frequencies');
inductance_45k_meas_all(kk,:)  =  (real_45K+1i*imag_45K)./(1i*2*pi*frequencies');
inductance_10k_meas_all(kk,:)  =  (real_10K+1i*imag_10K)./(1i*2*pi*frequencies');
inductance_test_meas_all(kk,:) =  (real_test+1i*imag_test)./(1i*2*pi*frequencies');
inductance_air_meas_all(kk,:)  =  (real_air+1i*imag_air)./(1i*2*pi*frequencies');

end


inductance_95k_meas =  mean(inductance_95k_meas_all);
inductance_45k_meas =  mean(inductance_45k_meas_all);
inductance_10k_meas =  mean(inductance_10k_meas_all);
inductance_test_meas=  mean(inductance_test_meas_all);
inductance_air_meas =  mean(inductance_air_meas_all);




number_of_layers = 10;
number_of_frequencies = length(frequencies); % at least equal to number of layers.

V=2.7e-3
I=11.92
s=7.5e-3
rho=2*pi*s*(V/I)
sigma=1/rho


conductivity_low = 1.10E+04; % measured with 4 point probe
conductivity_med = 3.85E+04;
conductivity_high = 8.6E+04;


search=true;
if search
    number_search = 10;
else 
    number_search=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% COIL TUNING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


conductivity_search = false;
if conductivity_search
    parameters = [59.5115000000000,61.6205000000000,64.5493000000000e-3,35.8267000000000e-3,0.00864190000000000,0.00436450000000000,conductivity_high,conductivity_med,conductivity_low];
    conductivity_high_initial = parameters(7);
    conductivity_med_initial = parameters(8);
    conductivity_low_initial = parameters(9);
else
    parameters = [59.5115000000000,61.6205000000000,64.5493000000000e-3,35.8267000000000e-3,0.00864190000000000,0.00436450000000000];
end



frequencies=[];
frequencies(1,:) = table_95K(:,2);
frequency_selection = round(linspace(2,11,length(parameters)));
frequencies =frequencies(frequency_selection);

inductance_differential_95K_meas = (inductance_95k_meas(frequency_selection)-inductance_air_meas(frequency_selection));
inductance_differential_45K_meas = (inductance_45k_meas(frequency_selection)-inductance_air_meas(frequency_selection));
inductance_differential_10K_meas = (inductance_10k_meas(frequency_selection)-inductance_air_meas(frequency_selection));
inductance_differential_test_meas = (inductance_test_meas(frequency_selection)-inductance_air_meas(frequency_selection));

% Bisection Conductivity

conductivity_high_upper = 87e3;
conductivity_high_lower = 80e3;

conductivity_med_upper = 42e3;
conductivity_med_lower = 36e3;

conductivity_low_upper = 13e3;
conductivity_low_lower = 10e3;

coil_thickness= 10e-3;
sigma_air=0.01;

turns_Tx = parameters(1);
turns_Rx = parameters(2);
inner_diameter_Tx= parameters(3);
inner_diameter_Rx= parameters(4);
lift_off_Tx = parameters(5);
lift_off_Rx = parameters(6);


save('C:\Users\EEE Admin\Desktop\ALL CODE\Applied Actual\Comsol\comsol_models.mat', 'frequencies', 'number_of_layers','lift_off_Tx','lift_off_Rx','inner_diameter_Tx','inner_diameter_Rx','coil_thickness','turns_Tx','turns_Rx','sigma_air')




inductance_in_air = inductance_air();
    
for kk=1:10
    
    
        
   
    s= ['Computing original'];
    disp(s);
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%    Highh conductivity block
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    conductivity_high_middle = (conductivity_high_lower+conductivity_high_upper)/2;
    calibration_conductivity = conductivity_high_middle*ones(10,1);
    inductance_unknown_brick = inductance_brick(calibration_conductivity);
    inductance_differential_95K_sim_0 = (inductance_unknown_brick - inductance_in_air);
    DMI_sim = [real(inductance_differential_95K_sim_0) imag(inductance_differential_95K_sim_0)];
    DMI_meas_95K = [real(inductance_differential_95K_meas) imag(inductance_differential_95K_meas)];
    
    res = (DMI_sim - DMI_meas_95K)./DMI_meas_95K;
    if sum(res)>0 % then the conductivity is too large
        conductivity_high_upper = conductivity_high_middle;
    else
        conductivity_high_lower = conductivity_high_middle;
    end
    
    
    if kk ==1
        inductance_differential_95K_sim_0_iteration_1 = inductance_differential_95K_sim_0;
    end
    calibration_95K = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2,1,1)
    semilogx(frequencies,real(inductance_differential_95K_meas),'k','linewidth',3)
    hold on
    scatter(frequencies,real(inductance_differential_95K_sim_0),100,'k','linewidth',3)
    hold on
    plot(frequencies,real(inductance_differential_95K_sim_0_iteration_1),'k:','linewidth',3)
    grid on
    set(gca,'fontsize',20,'yscale','log')
    ylabel('Real\{DMI\}','fontweight','bold','interpreter','latex','FontSize',35)
    xlabel('Frequency(Hz)','fontweight','bold','interpreter','latex','FontSize',35)
    legend({'Measured','Simulated','Start'},'fontweight','bold','interpreter','latex','FontSize',20,'location','southeast')
    
    subplot(2,1,2)
    semilogx(frequencies,imag(inductance_differential_95K_meas),'r','linewidth',3)
    hold on
    scatter(frequencies,imag(inductance_differential_95K_sim_0 ),100,'r','linewidth',3)
    hold on
    plot(frequencies,imag(inductance_differential_95K_sim_0_iteration_1),'r:','linewidth',3)
    grid on
    set(gca,'fontsize',20,'yscale','log')
    ylabel('Imag\{DMI\}','fontweight','bold','interpreter','latex','FontSize',35)
    xlabel('Frequency(Hz)','fontweight','bold','interpreter','latex','FontSize',35)
    legend({'Measured','Simulated','Start'},'fontweight','bold','interpreter','latex','FontSize',20,'location','southeast')
    
    saveas(calibration_95K,'Pictures\calibration_95K_tuning', 'fig')
    saveas(calibration_95K,'Pictures\calibration_95K_tuning', 'png')
    close(calibration_95K)
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%    med conductivity block
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    conductivity_med_middle = (conductivity_med_lower+conductivity_med_upper)/2;
    calibration_conductivity = conductivity_med_middle*ones(10,1);
    inductance_unknown_brick = inductance_brick(calibration_conductivity);
    inductance_differential_45K_sim_0 = (inductance_unknown_brick - inductance_in_air);
    DMI_sim = [real(inductance_differential_45K_sim_0) imag(inductance_differential_45K_sim_0)];
    DMI_meas_45K = [real(inductance_differential_45K_meas) imag(inductance_differential_45K_meas)];
    
    res = (DMI_sim - DMI_meas_45K)./DMI_meas_45K;
    if sum(res)>0 % then the conductivity is too large
        conductivity_med_upper = conductivity_med_middle;
    else
        conductivity_med_lower = conductivity_med_middle;
    end
    
    
    if kk ==1
        inductance_differential_45K_sim_0_iteration_1 = inductance_differential_45K_sim_0;
    end
    calibration_95K = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2,1,1)
    semilogx(frequencies,real(inductance_differential_45K_meas),'k','linewidth',3)
    hold on
    scatter(frequencies,real(inductance_differential_45K_sim_0),100,'k','linewidth',3)
    hold on
    plot(frequencies,real(inductance_differential_45K_sim_0_iteration_1),'k:','linewidth',3)
    grid on
    set(gca,'fontsize',20,'yscale','log')
    ylabel('Real\{DMI\}','fontweight','bold','interpreter','latex','FontSize',35)
    xlabel('Frequency(Hz)','fontweight','bold','interpreter','latex','FontSize',35)
    legend({'Measured','Simulated','Start'},'fontweight','bold','interpreter','latex','FontSize',20,'location','southeast')
    
    subplot(2,1,2)
    semilogx(frequencies,imag(inductance_differential_45K_meas),'r','linewidth',3)
    hold on
    scatter(frequencies,imag(inductance_differential_45K_sim_0 ),100,'r','linewidth',3)
    hold on
    plot(frequencies,imag(inductance_differential_45K_sim_0_iteration_1),'r:','linewidth',3)
    grid on
    set(gca,'fontsize',20,'yscale','log')
    ylabel('Imag\{DMI\}','fontweight','bold','interpreter','latex','FontSize',35)
    xlabel('Frequency(Hz)','fontweight','bold','interpreter','latex','FontSize',35)
    legend({'Measured','Simulated','Start'},'fontweight','bold','interpreter','latex','FontSize',20,'location','southeast')
    
    saveas(calibration_95K,'Pictures\calibration_45K_tuning', 'fig')
    saveas(calibration_95K,'Pictures\calibration_45K_tuning', 'png')
    close(calibration_95K)
    
     % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%    low conductivity block
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    conductivity_low_middle = (conductivity_low_lower+conductivity_low_upper)/2;
    calibration_conductivity = conductivity_low_middle*ones(10,1);
    inductance_unknown_brick = inductance_brick(calibration_conductivity);
    inductance_differential_10K_sim_0 = (inductance_unknown_brick - inductance_in_air);
    DMI_sim = [real(inductance_differential_10K_sim_0) imag(inductance_differential_10K_sim_0)];
    DMI_meas_10K = [real(inductance_differential_10K_meas) imag(inductance_differential_10K_meas)];
    
    res = (DMI_sim - DMI_meas_10K)./DMI_meas_10K;
    if sum(res)>0 % then the conductivity is too large
        conductivity_low_upper = conductivity_low_middle;
    else
        conductivity_low_lower = conductivity_low_middle;
    end
    
    
   if kk ==1
        inductance_differential_10K_sim_0_iteration_1 = inductance_differential_10K_sim_0;
    end 
    calibration_95K = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2,1,1)
    semilogx(frequencies,real(inductance_differential_10K_meas),'k','linewidth',3)
    hold on
    scatter(frequencies,real(inductance_differential_10K_sim_0),100,'k','linewidth',3)
    hold on
    plot(frequencies,real(inductance_differential_10K_sim_0_iteration_1),'k:','linewidth',3)
    grid on
    set(gca,'fontsize',20,'yscale','log')
    ylabel('Real\{DMI\}','fontweight','bold','interpreter','latex','FontSize',35)
    xlabel('Frequency(Hz)','fontweight','bold','interpreter','latex','FontSize',35)
    legend({'Measured','Simulated','Start'},'fontweight','bold','interpreter','latex','FontSize',20,'location','southeast')
    
    subplot(2,1,2)
    semilogx(frequencies,imag(inductance_differential_10K_meas),'r','linewidth',3)
    hold on
    scatter(frequencies,imag(inductance_differential_10K_sim_0 ),100,'r','linewidth',3)
    hold on
    plot(frequencies,imag(inductance_differential_10K_sim_0_iteration_1),'r:','linewidth',3)
    grid on
    set(gca,'fontsize',20,'yscale','log')
    ylabel('Imag\{DMI\}','fontweight','bold','interpreter','latex','FontSize',35)
    xlabel('Frequency(Hz)','fontweight','bold','interpreter','latex','FontSize',35)
    legend({'Measured','Simulated','Start'},'fontweight','bold','interpreter','latex','FontSize',20,'location','southeast')
    
    saveas(calibration_95K,'Pictures\calibration_10K_tuning', 'fig')
    saveas(calibration_95K,'Pictures\calibration_10K_tuning', 'png')
    close(calibration_95K)
end

conductivity_low = 1.0926e+04; % measured with 4 point probe
conductivity_med = 3.8086e+04;
conductivity_high = 8.6973e+04;

grad_last=NaN;
step_last = NaN(length(parameters),1);
B_last=NaN;
J_last=NaN(2*length(frequencies),length(parameters),3);
    
coil_thickness= 10e-3;
sigma_air=0.01;

normalise = false;
if normalise
    gamma=2.0922;
else
    gamma = 1e-8;
end
for kk=1:number_search
    
    start = tic;
    disp('Parameters New')
    s= ['Turns Tx = ' num2str(parameters(1)), ' || Turns Rx = ' num2str(parameters(2)), ' || Inner D Tx = ' num2str(parameters(3)*1e3) ' mm ' , ' || Inner D Rx = ' num2str(parameters(4)*1e3) ' mm ' , 'LO Tx= ' num2str(parameters(5)*1e3) ' mm' , ' || LO Rx= ' num2str(parameters(6)*1e3) ' mm ',' || High = ' num2str(conductivity_high/1e3) ' kS/m ',' || Med = ' num2str(conductivity_med/1e3) ' kS/m ',' || Low = ' num2str(conductivity_low/1e3) ' kS/m '];
    disp(s)

    if search
        
        turns_Tx = parameters(1);
        turns_Rx = parameters(2);
        inner_diameter_Tx= parameters(3);
        inner_diameter_Rx= parameters(4);
        lift_off_Tx = parameters(5);
        lift_off_Rx = parameters(6);
        
        if conductivity_search
            conductivity_high = parameters(7);
            conductivity_med = parameters(8);
            conductivity_low = parameters(9);
            
            conductivity_high_original = parameters(7);
            conductivity_med_original = parameters(8);
            conductivity_low_original = parameters(9);
            
        else 
            conductivity_high_original = conductivity_high;
            conductivity_med_original = conductivity_med;
            conductivity_low_original = conductivity_low;
        end
            
           
        
        save('C:\Users\EEE Admin\Desktop\ALL CODE\Applied Actual\Comsol\comsol_models.mat', 'frequencies', 'number_of_layers','lift_off_Tx','lift_off_Rx','inner_diameter_Tx','inner_diameter_Rx','coil_thickness','turns_Tx','turns_Rx','sigma_air')
    
    
    else
        offset= 0.0018 ; %offset_arr(kk);0.0018
        turns_Tx = 50;
        turns_Rx = 50;
        save('C:\Users\EEE Admin\Desktop\ALL CODE\Applied Actual\Comsol\comsol_models.mat', 'frequencies', 'number_of_layers','lift_off','offset','turns_Tx','turns_Rx')
    end
    
    
    
    inductance_in_air = inductance_air();
    s= ['Computing original'];
    disp(s);
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%    Highh conductivity block
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    calibration_conductivity = conductivity_high*ones(10,1);
    inductance_unknown_brick = inductance_brick(calibration_conductivity);
    inductance_differential_95K_sim_0 = (inductance_unknown_brick - inductance_in_air);
    DMI_sim = [real(inductance_differential_95K_sim_0) imag(inductance_differential_95K_sim_0)];
    DMI_meas_95K = [real(inductance_differential_95K_meas) imag(inductance_differential_95K_meas)];
    if normalise
        NORM = DMI_meas_95K;
    else
        NORM=ones(1,length(DMI_meas_95K));
    end
    r(:,1,1)= (DMI_sim-DMI_meas_95K)./NORM;
    
    
    if kk ==1
        inductance_differential_95K_sim_0_iteration_1 = inductance_differential_95K_sim_0;
    end 
    calibration_95K = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2,1,1)
    semilogx(frequencies,real(inductance_differential_95K_meas),'k','linewidth',3)
    hold on
    scatter(frequencies,real(inductance_differential_95K_sim_0),100,'k','linewidth',3)
    hold on
    plot(frequencies,real(inductance_differential_95K_sim_0_iteration_1),'k:','linewidth',3)
    grid on
    set(gca,'fontsize',20,'yscale','log')
    ylabel('Real\{DMI\}','fontweight','bold','interpreter','latex','FontSize',35)
    xlabel('Frequency(Hz)','fontweight','bold','interpreter','latex','FontSize',35)
    legend({'Measured','Simulated','Start'},'fontweight','bold','interpreter','latex','FontSize',20,'location','southeast')
    
    subplot(2,1,2)
    semilogx(frequencies,imag(inductance_differential_95K_meas),'r','linewidth',3)
    hold on
    scatter(frequencies,imag(inductance_differential_95K_sim_0 ),100,'r','linewidth',3)
    hold on
    plot(frequencies,imag(inductance_differential_95K_sim_0_iteration_1),'r:','linewidth',3)
    grid on
    set(gca,'fontsize',20,'yscale','log')
    ylabel('Imag\{DMI\}','fontweight','bold','interpreter','latex','FontSize',35)
    xlabel('Frequency(Hz)','fontweight','bold','interpreter','latex','FontSize',35)
    legend({'Measured','Simulated','Start'},'fontweight','bold','interpreter','latex','FontSize',20,'location','southeast')
    
    saveas(calibration_95K,'Pictures\calibration_95K_tuning', 'fig')
    saveas(calibration_95K,'Pictures\calibration_95K_tuning', 'png')
    close(calibration_95K)
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%    med conductivity block
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    calibration_conductivity = conductivity_med*ones(10,1);
    inductance_unknown_brick = inductance_brick(calibration_conductivity);
    inductance_differential_45K_sim_0 = (inductance_unknown_brick - inductance_in_air);
    DMI_sim = [real(inductance_differential_45K_sim_0) imag(inductance_differential_45K_sim_0)];
    DMI_meas_45K = [real(inductance_differential_45K_meas) imag(inductance_differential_45K_meas)];
    if normalise
        NORM= DMI_meas_45K;
    else
        NORM=ones(1,length(DMI_meas_45K));
    end
    r(:,1,2)= (DMI_sim-DMI_meas_45K)./NORM;
    
    
    if kk ==1
        inductance_differential_45K_sim_0_iteration_1 = inductance_differential_45K_sim_0;
    end 
    calibration_95K = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2,1,1)
    semilogx(frequencies,real(inductance_differential_45K_meas),'k','linewidth',3)
    hold on
    scatter(frequencies,real(inductance_differential_45K_sim_0),100,'k','linewidth',3)
    hold on
    plot(frequencies,real(inductance_differential_45K_sim_0_iteration_1),'k:','linewidth',3)
    grid on
    set(gca,'fontsize',20,'yscale','log')
    ylabel('Real\{DMI\}','fontweight','bold','interpreter','latex','FontSize',35)
    xlabel('Frequency(Hz)','fontweight','bold','interpreter','latex','FontSize',35)
    legend({'Measured','Simulated','Start'},'fontweight','bold','interpreter','latex','FontSize',20,'location','southeast')
    
    subplot(2,1,2)
    semilogx(frequencies,imag(inductance_differential_45K_meas),'r','linewidth',3)
    hold on
    scatter(frequencies,imag(inductance_differential_45K_sim_0 ),100,'r','linewidth',3)
    hold on
    plot(frequencies,imag(inductance_differential_45K_sim_0_iteration_1),'r:','linewidth',3)
    grid on
    set(gca,'fontsize',20,'yscale','log')
    ylabel('Imag\{DMI\}','fontweight','bold','interpreter','latex','FontSize',35)
    xlabel('Frequency(Hz)','fontweight','bold','interpreter','latex','FontSize',35)
    legend({'Measured','Simulated','Start'},'fontweight','bold','interpreter','latex','FontSize',20,'location','southeast')
    
    saveas(calibration_95K,'Pictures\calibration_45K_tuning', 'fig')
    saveas(calibration_95K,'Pictures\calibration_45K_tuning', 'png')
    close(calibration_95K)
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%    low conductivity block
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    calibration_conductivity = conductivity_low*ones(10,1);
    inductance_unknown_brick = inductance_brick(calibration_conductivity);
    inductance_differential_10K_sim_0 = (inductance_unknown_brick - inductance_in_air);
    
    DMI_sim = [real(inductance_differential_10K_sim_0) imag(inductance_differential_10K_sim_0)];
    DMI_meas_10K = [real(inductance_differential_10K_meas) imag(inductance_differential_10K_meas)];
    if normalise
        NORM= DMI_meas_10K;
    else
        NORM=ones(1,length(DMI_meas_10K));
    end
    r(:,1,3)= (DMI_sim-DMI_meas_10K)./NORM;
    
    
   if kk ==1
        inductance_differential_10K_sim_0_iteration_1 = inductance_differential_10K_sim_0;
    end 
    calibration_95K = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2,1,1)
    semilogx(frequencies,real(inductance_differential_10K_meas),'k','linewidth',3)
    hold on
    scatter(frequencies,real(inductance_differential_10K_sim_0),100,'k','linewidth',3)
    hold on
    plot(frequencies,real(inductance_differential_10K_sim_0_iteration_1),'k:','linewidth',3)
    grid on
    set(gca,'fontsize',20,'yscale','log')
    ylabel('Real\{DMI\}','fontweight','bold','interpreter','latex','FontSize',35)
    xlabel('Frequency(Hz)','fontweight','bold','interpreter','latex','FontSize',35)
    legend({'Measured','Simulated','Start'},'fontweight','bold','interpreter','latex','FontSize',20,'location','southeast')
    
    subplot(2,1,2)
    semilogx(frequencies,imag(inductance_differential_10K_meas),'r','linewidth',3)
    hold on
    scatter(frequencies,imag(inductance_differential_10K_sim_0 ),100,'r','linewidth',3)
    hold on
    plot(frequencies,imag(inductance_differential_10K_sim_0_iteration_1),'r:','linewidth',3)
    grid on
    set(gca,'fontsize',20,'yscale','log')
    ylabel('Imag\{DMI\}','fontweight','bold','interpreter','latex','FontSize',35)
    xlabel('Frequency(Hz)','fontweight','bold','interpreter','latex','FontSize',35)
    legend({'Measured','Simulated','Start'},'fontweight','bold','interpreter','latex','FontSize',20,'location','southeast')
    
    saveas(calibration_95K,'Pictures\calibration_10K_tuning', 'fig')
    saveas(calibration_95K,'Pictures\calibration_10K_tuning', 'png')
    close(calibration_95K)
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%    Sensitiity
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    f1(kk) = norm(r(:,:,1))^2;
    f2(kk) = norm(r(:,:,2))^2;
    f3(kk) = norm(r(:,:,3))^2;
    
    resid_curves = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2,2,1)
    title('High')
    scatter(1:kk,f1(1:kk),100,'linewidth',2)
    grid on
    set(gca,'yscale','log')
    subplot(2,2,2)
    title('Med')
    scatter(1:kk,f2(1:kk),100,'linewidth',2)
    grid on
    set(gca,'yscale','log')
    subplot(2,2,[3,4])
    title('Low')
    scatter(1:kk,f3(1:kk),100,'linewidth',2)
    grid on
    set(gca,'yscale','log')
    saveas(resid_curves,'Pictures\callibration_residual_curves', 'fig')
    saveas(resid_curves,'Pictures\callibration_residual_curves', 'png')
    close(resid_curves)
    
    purturbations = eye(length(parameters),length(parameters))*0.01*diag((parameters));
        
    s= ['Computing Perturbed'];
    disp(s);
    for mm=1:length(parameters)
        s= ['     Parameter = ' num2str(mm)];
        disp(s);
        
        turns_Tx = parameters(1)+purturbations(mm,1);
        turns_Rx = parameters(2)+purturbations(mm,2);
        inner_diameter_Tx= parameters(3)+purturbations(mm,3);
        inner_diameter_Rx= parameters(4)+purturbations(mm,4);
        lift_off_Tx = parameters(5)+purturbations(mm,5);
        lift_off_Rx = parameters(6)+purturbations(mm,6);
        
        
        if conductivity_search
            conductivity_high = parameters(7)+purturbations(mm,7);
            conductivity_med = parameters(8)+purturbations(mm,8);
            conductivity_low = parameters(9)+purturbations(mm,9);
        end
        
        save('C:\Users\EEE Admin\Desktop\ALL CODE\Applied Actual\Comsol\comsol_models.mat', 'frequencies', 'number_of_layers','lift_off_Tx','lift_off_Rx','inner_diameter_Tx','inner_diameter_Rx','coil_thickness','turns_Tx','turns_Rx','sigma_air')
    
        
        s= ['           air'];
        disp(s);
        inductance_in_air = inductance_air();
        
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % %%%%%%%%%%%%    Highh conductivity block
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        s= ['           High Cond'];
        disp(s);
        calibration_conductivity = conductivity_high*ones(10,1);
        inductance_unknown_brick = inductance_brick(calibration_conductivity);
        inductance_differential_95K_sim = (inductance_unknown_brick - inductance_in_air);
        
        DMI_pert = [real(inductance_differential_95K_sim) imag(inductance_differential_95K_sim)];
        DMI_orig = [real(inductance_differential_95K_sim_0) imag(inductance_differential_95K_sim_0)];
        if normalise
            NORM= DMI_meas_95K;
        else
            NORM=ones(1,length(DMI_meas_95K));
        end

        J(:,mm,1) = ((DMI_pert-DMI_orig)./NORM)./purturbations(mm,mm);
        
        
        
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % %%%%%%%%%%%%    med conductivity block
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        s= ['           Med Cond'];
        disp(s);
        
        calibration_conductivity = conductivity_med*ones(10,1);
        inductance_unknown_brick = inductance_brick(calibration_conductivity);
        inductance_differential_45K_sim = (inductance_unknown_brick - inductance_in_air);
        
        DMI_pert = [real(inductance_differential_45K_sim) imag(inductance_differential_45K_sim)];
        DMI_orig = [real(inductance_differential_45K_sim_0) imag(inductance_differential_45K_sim_0)];
        
        if normalise
            NORM= DMI_meas_45K;
        else
            NORM=ones(1,length(DMI_meas_45K));
        end

        J(:,mm,2) = ((DMI_pert-DMI_orig)./NORM)./purturbations(mm,mm);
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % %%%%%%%%%%%%    low conductivity block
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        s= ['           Low Cond'];
        disp(s);
        
        calibration_conductivity = conductivity_low*ones(10,1);
        inductance_unknown_brick = inductance_brick(calibration_conductivity);
        inductance_differential_10K_sim = (inductance_unknown_brick - inductance_in_air);
        
        DMI_pert = [real(inductance_differential_10K_sim) imag(inductance_differential_10K_sim)];
        DMI_orig = [real(inductance_differential_10K_sim_0) imag(inductance_differential_10K_sim_0)];
        
        if normalise
            NORM= DMI_meas_10K;
        else
            NORM=ones(1,length(DMI_meas_10K));
        end

        J(:,mm,3) = ((DMI_pert-DMI_orig)./NORM)./purturbations(mm,mm);
        
    end
    
    
    %     
    %     gamma = 1e-50;
    %
    %     step = -inv(H+D*gamma)*grad;
    turns_upper = 90;
    turns_lower = 40;
    lift_off_upper = 20e-3;
    lift_off_lower =  3e-3;
    
    for mm=1:3
        Gradients(:,mm) = J(:,:,mm)'*r(:,1,mm);
        Hessian(:,:,mm) = J(:,:,mm)'*J(:,:,mm);
        y_s(:,:,mm) = J(:,:,mm)'*r(:,1,mm)-J_last(:,:,mm)'*r(:,1,mm);
    end

    
    grad = sum(Gradients,2); % don't average if just summing in objective function
    y_sharp =sum(y_s,3);
    
  
    
    y = grad-grad_last;
    
    s = step_last;
    
    
    if (y'*s)==0
        rho =0;
    else
        rho = 1/(y'*s);
    end
    
    if isnan(rho)
        B = zeros(length(parameters));
    elseif ~isnan(rho) && ~isinf(norm(rho))
        tau = min([1 norm(s'*y_sharp)/norm(s'*B_last*s)]);
        B_last = tau*B_last;
        term = (y_sharp-B_last*s)*rho;
        B = B_last+term*y'+y*term' -(term'*rho)*s*y*y'; 
    end
    
    
    H = sum(Hessian,3)+B;
    
    D = diag(diag(H)/max(diag(H)));
    
    ratio=0;
    
    upper_factor=4;
    lower_factor=3;
    count=1;
    while ratio <0.8 || ratio>1.2 || ratio<(1+(1-0.99)) && ratio > 1 || ratio<1 && ratio > 0.99
        
        if count == 10
            break
        end
        count=count+1;
        
        s= ['Computing damping parameter'];
        disp(s);
    
        step = -inv(H+D*gamma)*grad;
        
        turns_Tx = parameters(1)+step(1);
        turns_Rx = parameters(2)+step(2);
        inner_diameter_Tx= parameters(3)+step(3);
        inner_diameter_Rx= parameters(4)+step(4);
        lift_off_Tx = parameters(5)+step(5);
        lift_off_Rx = parameters(6)+step(6);
        
        diameter_of_wire = 2*sqrt(2.827431e-7/pi);
        area_wire_takes = diameter_of_wire*diameter_of_wire;
        
        area =area_wire_takes *turns_Tx;
        width = area/5e-3;
        
        outer_diameter_Tx = inner_diameter_Tx+2*width;
        
        area = area_wire_takes*turns_Rx;
        width = area/5e-3;
        
        outer_diameter_Rx = inner_diameter_Rx+2*width;

        diameter_difference = inner_diameter_Tx-(outer_diameter_Rx);
        
        Turns_Tx_condition = turns_Tx <turns_lower || turns_Tx >turns_upper;
        Turns_Rx_condition = turns_Rx <turns_lower || turns_Rx >turns_upper;
        diameter_condition = diameter_difference < 1e-3 || inner_diameter_Rx <1e-3;
        
        if conductivity_search
            conductivity_high=parameters(7)+step(7);
            conductivity_med=parameters(8)+step(8);
            conductivity_low=parameters(9)+step(9);
            bound = 4e3;
            cond_step_condition_upper = (parameters(7)+step(7))>(conductivity_high_initial+bound) || (parameters(8)+step(8))>(conductivity_med_initial+bound) || (parameters(9)+step(9))>(conductivity_low_initial+bound);
            cond_step_condition_lower = (parameters(7)+step(7))>(conductivity_high_initial-bound) || (parameters(8)+step(8))>(conductivity_med_initial-bound) || (parameters(9)+step(9))>(conductivity_low_initial-bound);
            cond_step_condition = cond_step_condition_upper || cond_step_condition_lower;
        else
            cond_step_condition = false;
        end
        
        
        
        lift_off_condition = lift_off_Tx>lift_off_upper || lift_off_Tx<lift_off_lower || lift_off_Rx>lift_off_upper || lift_off_Rx<lift_off_lower;
        
        
        
        if Turns_Tx_condition || Turns_Rx_condition || lift_off_condition || diameter_condition || cond_step_condition
            low = log10(gamma);
            upper =5;
            for ii=1:100
                mid = (low + upper)/2;
                
                step = -inv(H+D*10^mid)*grad;
                
                
                
                turns_Tx = parameters(1)+step(1);
                turns_Rx = parameters(2)+step(2);
                inner_diameter_Tx= parameters(3)+step(3);
                inner_diameter_Rx= parameters(4)+step(4);
                lift_off_Tx = parameters(5)+step(5);
                lift_off_Rx = parameters(6)+step(6);
                
                diameter_of_wire = 2*sqrt(2.827431e-7/pi);
                area_wire_takes = diameter_of_wire*diameter_of_wire;
                
                area =area_wire_takes *turns_Tx;
                width = area/5e-3;
                
                outer_diameter_Tx = inner_diameter_Tx+2*width;
                
                area = area_wire_takes*turns_Rx;
                width = area/5e-3;
                
                outer_diameter_Rx = inner_diameter_Rx+2*width;
                
                diameter_difference = inner_diameter_Tx-(outer_diameter_Rx);
                
                Turns_Tx_condition = turns_Tx <turns_lower || turns_Tx >turns_upper;
                Turns_Rx_condition = turns_Rx <turns_lower || turns_Rx >turns_upper;
                diameter_condition = diameter_difference < 1e-3 || inner_diameter_Rx <1e-3;
                
                lift_off_condition = lift_off_Tx>lift_off_upper || lift_off_Tx<lift_off_lower || lift_off_Rx>lift_off_upper || lift_off_Rx<lift_off_lower;
                
                
                if conductivity_search
                    conductivity_high=parameters(7)+step(7);
                    conductivity_med=parameters(8)+step(8);
                    conductivity_low=parameters(9)+step(9);
                    bound = 4e3;
                    cond_step_condition_upper = (parameters(7)+step(7))>(conductivity_high_initial+bound) || (parameters(8)+step(8))>(conductivity_med_initial+bound) || (parameters(9)+step(9))>(conductivity_low_initial+bound);
                    cond_step_condition_lower = (parameters(7)+step(7))>(conductivity_high_initial-bound) || (parameters(8)+step(8))>(conductivity_med_initial-bound) || (parameters(9)+step(9))>(conductivity_low_initial-bound);
                    cond_step_condition = cond_step_condition_upper || cond_step_condition_lower;
                else
                    cond_step_condition = false;
                end
                
                if   Turns_Tx_condition || Turns_Rx_condition || lift_off_condition || diameter_condition || cond_step_condition
                    low =  mid;
                else
                    upper =mid;
                end
                
            end
            gamma=10^mid;
            
            step = -inv(H+D*gamma)*grad;
            
        end
        
        
          
        turns_Tx = parameters(1)+step(1);
        turns_Rx = parameters(2)+step(2);
        inner_diameter_Tx= parameters(3)+step(3);
        inner_diameter_Rx= parameters(4)+step(4);
        lift_off_Tx = parameters(5)+step(5);
        lift_off_Rx = parameters(6)+step(6);
        

        if conductivity_search
            conductivity_high = parameters(7)+step(7);
            conductivity_med = parameters(8)+step(8);
            conductivity_low = parameters(9)+step(9);
        end
         save('C:\Users\EEE Admin\Desktop\ALL CODE\Applied Actual\Comsol\comsol_models.mat', 'frequencies', 'number_of_layers','lift_off_Tx','lift_off_Rx','inner_diameter_Tx','inner_diameter_Rx','coil_thickness','turns_Tx','turns_Rx','sigma_air')
        
        
        s= ['           air'];
        disp(s);
        inductance_in_air = inductance_air();
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%    Highh conductivity block
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        s= ['           High'];
        disp(s);
        calibration_conductivity = conductivity_high*ones(10,1);
        inductance_unknown_brick = inductance_brick(calibration_conductivity);
        inductance_differential_95K_sim = (inductance_unknown_brick - inductance_in_air);
        
        DMI_pred = [real(inductance_differential_95K_sim) imag(inductance_differential_95K_sim)];
        DMI_meas = [real(inductance_differential_95K_meas) imag(inductance_differential_95K_meas)];
        
        
        if normalise
            NORM= DMI_meas_95K;
        else
            NORM=ones(1,length(DMI_meas_95K));
        end
        r_new(:,1,1)= (DMI_pred-DMI_meas)./NORM;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%    med conductivity block
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        s= ['           Med'];
        disp(s);
        calibration_conductivity = conductivity_med*ones(10,1);
        inductance_unknown_brick = inductance_brick(calibration_conductivity);
        inductance_differential_45K_sim = (inductance_unknown_brick - inductance_in_air);
        
        DMI_pred = [real(inductance_differential_45K_sim) imag(inductance_differential_45K_sim)];
        DMI_meas = [real(inductance_differential_45K_meas) imag(inductance_differential_45K_meas)];
 
        
        if normalise
            NORM= DMI_meas_45K;
        else
            NORM=ones(1,length(DMI_meas_45K));
        end
        r_new(:,1,2)= (DMI_pred-DMI_meas)./NORM;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%    low conductivity block
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        s= ['           Low'];
        disp(s);
        calibration_conductivity = conductivity_low*ones(10,1);
        inductance_unknown_brick = inductance_brick(calibration_conductivity);
        inductance_differential_10K_sim = (inductance_unknown_brick - inductance_in_air);
        
        DMI_pred = [real(inductance_differential_10K_sim) imag(inductance_differential_10K_sim)];
        DMI_meas = [real(inductance_differential_10K_meas) imag(inductance_differential_10K_meas)];
        
        if normalise
            NORM= DMI_meas_10K;
        else
            NORM=ones(1,length(DMI_meas_10K));
        end
        r_new(:,1,3)= (DMI_pred-DMI_meas)./NORM;
        
        m_0 = 0.5*norm(r(:,:,1))^2+0.5*norm(r(:,:,2))^2+0.5*norm(r(:,:,3))^2;
        m_p = m_0+step'*grad+0.5*step'*H*step;
        
        f_0 = m_0;
        f = 0.5*norm(r_new(:,:,1))^2+0.5*norm(r_new(:,:,2))^2+0.5*norm(r_new(:,:,3))^2;
        
        delta_f= f_0 - f;
        delta_m= m_0 - m_p;
        
        ratio = delta_f/delta_m;
        
        s = ['Delta F = ' num2str(delta_f) ' || Delta M = ' num2str(delta_m) ' || ratio = ' num2str(ratio) ' || Gamma = ' num2str(gamma)];
        disp(s)
        if ratio <0.8 || ratio>1.2
            gamma=gamma*upper_factor;
%             if upper_factor>1
%                 upper_factor = upper_factor*0.9;
%             end
        elseif  ratio<(1+(1-0.99)) && ratio > 1 || ratio<1 && ratio > 0.99
            gamma=gamma/lower_factor;
%             if lower_factor>1
%                 lower_factor = lower_factor*0.9;
%             end
        end
    end
    disp('Parameters Old')
    s= ['Turns Tx = ' num2str(parameters(1)), ' || Turns Rx = ' num2str(parameters(2)), ' || Inner D Tx = ' num2str(parameters(3)*1e3) ' mm ' , ' || Inner D Rx = ' num2str(parameters(4)*1e3) ' mm ' , 'LO Tx= ' num2str(parameters(5)*1e3) ' mm' , ' || LO Rx= ' num2str(parameters(6)*1e3) ' mm ',' || High = ' num2str(conductivity_high_original/1e3) ' kS/m ',' || Med = ' num2str(conductivity_med_original/1e3) ' kS/m ',' || Low = ' num2str(conductivity_low_original/1e3) ' kS/m '];
    disp(s)
    
    if count == 10
        
    else
        parameters = parameters+step';
        grad_last=grad;
        step_last = step;
        B_last=B;
        J_last=J;
        
    end
    
    
    
    
end

% frequencies=[];
% extension = ['_data_' num2str(1) '.xlsx'];
% table_95K = xlsread(['C:\Users\EEE Admin\Desktop\measured_data\95k' extension ]);
% frequencies(1,:) = table_95K(:,2);
    
% turns_Tx = parameters(1);
% turns_Rx = parameters(1);
% offset= parameters(2);
% lift_off = parameters(3);
% save('Comsol\comsol_models.mat', 'frequencies', 'number_of_layers','lift_off','offset','turns_Tx','turns_Rx')

frequencies=[];
frequencies(1,:) = table_95K(:,2);
frequency_selection = 2:11;
frequencies =frequencies(frequency_selection);

inductance_differential_95K_meas = (inductance_95k_meas(frequency_selection)-inductance_air_meas(frequency_selection));
inductance_differential_45K_meas = (inductance_45k_meas(frequency_selection)-inductance_air_meas(frequency_selection));
inductance_differential_10K_meas = (inductance_10k_meas(frequency_selection)-inductance_air_meas(frequency_selection));
inductance_differential_test_meas = (inductance_test_meas(frequency_selection)-inductance_air_meas(frequency_selection));


turns_Tx = parameters(1);
turns_Rx = parameters(2);
inner_diameter_Tx= parameters(3);
inner_diameter_Rx= parameters(4);
lift_off_Tx = parameters(5);
lift_off_Rx = parameters(6);
        

        
save('C:\Users\EEE Admin\Desktop\ALL CODE\Applied Actual\Comsol\comsol_models.mat', 'frequencies', 'number_of_layers','lift_off_Tx','lift_off_Rx','inner_diameter_Tx','inner_diameter_Rx','coil_thickness','turns_Tx','turns_Rx','sigma_air')

        
inductance_in_air = inductance_air();


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%    Highh conductivity block
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

calibration_conductivity = conductivity_high*ones(10,1);
inductance_unknown_brick = inductance_brick(calibration_conductivity);
inductance_differential_95K_sim = (inductance_unknown_brick - inductance_in_air);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%    med conductivity block
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


calibration_conductivity = conductivity_med*ones(10,1);
inductance_unknown_brick = inductance_brick(calibration_conductivity);
inductance_differential_45K_sim = (inductance_unknown_brick - inductance_in_air);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%    low conductivity block
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


calibration_conductivity = conductivity_low*ones(10,1);
inductance_unknown_brick = inductance_brick(calibration_conductivity);
inductance_differential_10K_sim = (inductance_unknown_brick - inductance_in_air);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% CALIBRATION FACTORS
%%%%%%%%%%%%%%%%%%%%%%%%%%%

DATA_abs = [abs(inductance_differential_10K_meas) ; abs(inductance_differential_45K_meas) ; abs(inductance_differential_95K_meas)];
std_DATA_abs  = std(DATA_abs);
DATA_abs = DATA_abs./std_DATA_abs ;

DATA_abs_targets = [abs(inductance_differential_10K_sim) ; abs(inductance_differential_45K_sim) ; abs(inductance_differential_95K_sim )];
DATA_abs_targets = DATA_abs_targets./std_DATA_abs ;

DATA_angle = [unwrap(angle(inductance_differential_10K_meas.'))  unwrap(angle(inductance_differential_45K_meas.'))  unwrap(angle(inductance_differential_95K_meas.'))  ]';
DATA_angle_targets = [unwrap(angle(inductance_differential_10K_sim.'))  unwrap(angle(inductance_differential_45K_sim.')) unwrap(angle(inductance_differential_95K_sim.')) ]';


coefficients_abs=[];
coefficients_angle=[];

X_abs = [DATA_abs(1,:) DATA_abs(2,:) DATA_abs(3,:)]';%

X_abs(:,2) = ones(length(X_abs(:,1)),1);

Y_abs = [DATA_abs_targets(1,:) DATA_abs_targets(2,:) DATA_abs_targets(3,:)]';

coefficients_abs = pinv(X_abs'*X_abs)*X_abs'*Y_abs;

%%%% Don't need to normalise for angle as all of the same magnitude.

X_angle = [DATA_angle(1,:) DATA_angle(2,:) DATA_angle(3,:)]';%  ones(length(inductance_differential_10K_meas),1)];
X_angle(:,2) = ones(length(X_angle(:,1)),1);
Y_angle= [DATA_angle_targets(1,:) DATA_angle_targets(2,:) DATA_angle_targets(3,:)]';
coefficients_angle = pinv(X_angle'*X_angle)*X_angle'*Y_angle;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%    Callibrate
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inductance_calibrated_amplitude_95K_meas_cal = ((abs(inductance_differential_95K_meas))./std_DATA_abs)*coefficients_abs(1)+coefficients_abs(2);
inductance_calibrated_amplitude_95K_meas_cal = inductance_calibrated_amplitude_95K_meas_cal.*std_DATA_abs;
inductance_calibrated_angle_95K_meas_cal = unwrap(angle(inductance_differential_95K_meas))*coefficients_angle(1)+coefficients_angle(2);



inductance_calibrated_amplitude_45K_meas_cal = ((abs(inductance_differential_45K_meas))./std_DATA_abs )*coefficients_abs(1)+coefficients_abs(2);
inductance_calibrated_amplitude_45K_meas_cal = inductance_calibrated_amplitude_45K_meas_cal.*std_DATA_abs ;
inductance_calibrated_angle_45K_meas_cal = unwrap(angle(inductance_differential_45K_meas))*coefficients_angle(1)+coefficients_angle(2);


inductance_calibrated_amplitude_10K_meas_cal = ((abs(inductance_differential_10K_meas))./std_DATA_abs )*coefficients_abs(1)+coefficients_abs(2);
inductance_calibrated_amplitude_10K_meas_cal = inductance_calibrated_amplitude_10K_meas_cal.*std_DATA_abs ;
inductance_calibrated_angle_10K_meas_cal = unwrap(angle(inductance_differential_10K_meas))*coefficients_angle(1)+coefficients_angle(2);


inductance_differential_95K_meas_cal=  inductance_calibrated_amplitude_95K_meas_cal.*(cos(inductance_calibrated_angle_95K_meas_cal )+1i*sin(inductance_calibrated_angle_95K_meas_cal ));

inductance_differential_45K_meas_cal=  inductance_calibrated_amplitude_45K_meas_cal.*(cos(inductance_calibrated_angle_45K_meas_cal)+1i*sin(inductance_calibrated_angle_45K_meas_cal));

inductance_differential_10K_meas_cal=  inductance_calibrated_amplitude_10K_meas_cal.*(cos(inductance_calibrated_angle_10K_meas_cal)+1i*sin(inductance_calibrated_angle_10K_meas_cal));


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%    test
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 

% %inductance_differential_test_sim = inductance_brick_test() - inductance_air_test();
% calibration_conductivity = [41000,41000,41000,41000,81600,93700,93700,93700,93700,93700];
% inductance_unknown_brick = inductance_brick(calibration_conductivity);
% inductance_differential_test_sim = (inductance_unknown_brick - inductance_in_air);


inductance_calibrated_amplitude = abs((inductance_differential_test_meas)./std_DATA_abs)*coefficients_abs(1)+coefficients_abs(2);
inductance_calibrated_amplitude = inductance_calibrated_amplitude.*std_DATA_abs;
inductance_calibrated_angle = unwrap(angle(inductance_differential_test_meas))*coefficients_angle(1)+coefficients_angle(2);

inductance_differential_test_cal =  inductance_calibrated_amplitude.*(cos(inductance_calibrated_angle)+1i*sin(inductance_calibrated_angle));


inductance_differential_test = inductance_differential_test_cal;


inductance_in_air = inductance_air();

save('DATA\unknown_test_data.mat','inductance_differential_test','inductance_in_air')


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%    Discrepancy
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


discrepancy_complex_all(1,:)=inductance_differential_10K_sim - inductance_differential_10K_meas_cal;
discrepancy_complex_all(2,:)=inductance_differential_45K_sim - inductance_differential_45K_meas_cal;
discrepancy_complex_all(3,:)=inductance_differential_95K_sim - inductance_differential_95K_meas_cal;

discrepancy_sep = [real(discrepancy_complex_all) imag(discrepancy_complex_all)];
envelope = max(abs(discrepancy_sep));

discrepancy_full = envelope;

largest_SNR = find(max(abs(discrepancy_full))==abs(discrepancy_full));
discrepancy_full_floor(1,:) = discrepancy_full(largest_SNR(1))*ones(length(discrepancy_full),1);

save('DATA\discrepancy.mat','discrepancy_full','discrepancy_full_floor')


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%    Plot discrepancy
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dis = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,1,1)
plot(frequencies/1e3,abs(real(discrepancy_complex_all(1,:)))*1e9,'k--','linewidth',3)
hold on 
plot(frequencies/1e3,abs(real(discrepancy_complex_all(2,:)))*1e9,'k:','linewidth',3)
hold on 
plot(frequencies/1e3,abs(real(discrepancy_complex_all(3,:)))*1e9,'k','linewidth',3)
hold on 
scatter(frequencies/1e3,envelope(1:length(envelope)/2)*1e9,200,'k','linewidth',3)
grid on
set(gca,'fontsize',35)
% s = ['$ \frac{ \mid \mid \vec{d}_t \mid \mid - \mid \mid \vec{d}_p \mid \mid }{\mid \mid \vec{d}_t \mid \mid} * 100  $ = ' num2str(error) '\%'];
 title('(a)','fontweight','bold','interpreter','latex','FontSize',40)
xlabel('Frequency(kHz)','fontweight','bold','interpreter','latex','FontSize',40)
ylabel('$|| d_r ||$ (nH)','fontweight','bold','interpreter','latex','FontSize',60)
legend({'10K','45K','95K','pred.'},'fontweight','bold','interpreter','latex','FontSize',40)

subplot(2,1,2)
plot(frequencies/1e3,abs(imag(discrepancy_complex_all(1,:)))*1e9,'r--','linewidth',3)
hold on 
plot(frequencies/1e3,abs(imag(discrepancy_complex_all(2,:)))*1e9,'r:','linewidth',3)
hold on 
plot(frequencies/1e3,abs(imag(discrepancy_complex_all(3,:)))*1e9,'r','linewidth',3)
hold on 
scatter(frequencies/1e3,envelope(length(envelope)/2+1:end)*1e9,200,'r','linewidth',3)
grid on
set(gca,'fontsize',35)
% s = ['$ \frac{ \mid \mid \vec{d}_t \mid \mid - \mid \mid \vec{d}_p \mid \mid }{\mid \mid \vec{d}_t \mid \mid} * 100  $ = ' num2str(error) '\%'];
 title('(b)','fontweight','bold','interpreter','latex','FontSize',40)
xlabel('Frequency(kHz)','fontweight','bold','interpreter','latex','FontSize',40)
ylabel('$|| d_i ||$ (nH)','fontweight','bold','interpreter','latex','FontSize',60)
legend({'10K','45K','95K','pred.'},'fontweight','bold','interpreter','latex','FontSize',40)
saveas(dis,'Pictures\discrepancy_complex', 'fig')
saveas(dis,'Pictures\discrepancy_complex', 'png')
close(dis)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%    Plot SNR
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



SNR_95 = 20*log10(abs(inductance_differential_95K_meas_cal)./abs(inductance_differential_95K_meas_cal-inductance_differential_95K_sim));
SNR_45 = 20*log10(abs(inductance_differential_45K_meas_cal)./abs(inductance_differential_45K_meas_cal-inductance_differential_45K_sim));
SNR_10 = 20*log10(abs(inductance_differential_10K_meas_cal)./abs(inductance_differential_10K_meas_cal-inductance_differential_10K_sim));


dis = figure('units','normalized','outerposition',[0 0 1 1]);
scatter(frequencies/1e3,SNR_95 ,400,'r+','linewidth',4)
hold on 
scatter(frequencies/1e3,SNR_45 ,400,'b+','linewidth',4)
hold on 
scatter(frequencies/1e3,SNR_10 ,400,'g+','linewidth',4)
grid on
ylim([0 60])
set(gca,'xscale','log','fontsize',35)
xlabel('Frequency(kHz)','fontweight','bold','interpreter','latex','FontSize',50)
ylabel('SNR (dB)','fontweight','bold','interpreter','latex','FontSize',50)
legend({'95K','45K','10K'},'fontweight','bold','interpreter','latex','FontSize',50)
saveas(dis,'Pictures\SNR', 'fig')
saveas(dis,'Pictures\SNR', 'png')
close(dis)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
calibration_95K = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,1,1)
semilogx(frequencies,real(inductance_differential_95K_meas),'k','linewidth',3)
hold on 
scatter(frequencies,real(inductance_differential_95K_sim),100,'k','linewidth',3)
hold on
semilogx(frequencies,real(inductance_differential_95K_meas_cal),'k-.','linewidth',3)
grid on
set(gca,'fontsize',20,'yscale','log')
ylabel('Real\{DMI\}','fontweight','bold','interpreter','latex','FontSize',35)
xlabel('Frequency(Hz)','fontweight','bold','interpreter','latex','FontSize',35)
legend({'Measured','Simulated','Calibrated'},'fontweight','bold','interpreter','latex','FontSize',20,'location','southeast')

subplot(2,1,2)
semilogx(frequencies,imag(inductance_differential_95K_meas),'r','linewidth',3)
hold on 
scatter(frequencies,imag(inductance_differential_95K_sim),100,'r','linewidth',3)
hold on
semilogx(frequencies,imag(inductance_differential_95K_meas_cal),'r-.','linewidth',3)
grid on
set(gca,'fontsize',20,'yscale','log')
ylabel('Imag\{DMI\}','fontweight','bold','interpreter','latex','FontSize',35)
xlabel('Frequency(Hz)','fontweight','bold','interpreter','latex','FontSize',35)
legend({'Measured','Simulated','Calibrated'},'fontweight','bold','interpreter','latex','FontSize',20,'location','southeast')
saveas(calibration_95K,'Pictures\calibration_95K', 'fig')
saveas(calibration_95K,'Pictures\calibration_95K', 'png')
close(calibration_95K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

calibration_45K = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,1,1)
semilogx(frequencies,real(inductance_differential_45K_meas),'k','linewidth',3)
hold on 
scatter(frequencies,real(inductance_differential_45K_sim),100,'k','linewidth',3)
hold on
semilogx(frequencies,real(inductance_differential_45K_meas_cal),'k-.','linewidth',3)
grid on
set(gca,'fontsize',20,'yscale','log')
ylabel('Real\{DMI\}','fontweight','bold','interpreter','latex','FontSize',35)
xlabel('Frequency(Hz)','fontweight','bold','interpreter','latex','FontSize',35)
legend({'Measured','Simulated','Calibrated'},'fontweight','bold','interpreter','latex','FontSize',20,'location','southeast')

subplot(2,1,2)
semilogx(frequencies,imag(inductance_differential_45K_meas),'r','linewidth',3)
hold on 
scatter(frequencies,imag(inductance_differential_45K_sim),100,'r','linewidth',3)
hold on
semilogx(frequencies,imag(inductance_differential_45K_meas_cal),'r-.','linewidth',3)
grid on
set(gca,'fontsize',20,'yscale','log')
ylabel('Imag\{DMI\}','fontweight','bold','interpreter','latex','FontSize',35)
xlabel('Frequency(Hz)','fontweight','bold','interpreter','latex','FontSize',35)
legend({'Measured','Simulated','Calibrated'},'fontweight','bold','interpreter','latex','FontSize',20,'location','southeast')
saveas(calibration_45K,'Pictures\calibration_45K', 'fig')
saveas(calibration_45K,'Pictures\calibration_45K', 'png')
close(calibration_45K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
calibration_10K = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,1,1)
semilogx(frequencies,real(inductance_differential_10K_meas),'k','linewidth',3)
hold on 
scatter(frequencies,real(inductance_differential_10K_sim),100,'k','linewidth',3)
hold on
semilogx(frequencies,real(inductance_differential_10K_meas_cal),'k-.','linewidth',3)
grid on
set(gca,'fontsize',20,'yscale','log')
ylabel('Real\{DMI\}','fontweight','bold','interpreter','latex','FontSize',35)
xlabel('Frequency(Hz)','fontweight','bold','interpreter','latex','FontSize',35)
legend({'Measured','Simulated','Calibrated'},'fontweight','bold','interpreter','latex','FontSize',20,'location','southeast')

subplot(2,1,2)
semilogx(frequencies,imag(inductance_differential_10K_meas),'r','linewidth',3)
hold on 
scatter(frequencies,imag(inductance_differential_10K_sim),100,'r','linewidth',3)
hold on
semilogx(frequencies,imag(inductance_differential_10K_meas_cal),'r-.','linewidth',3)
grid on
set(gca,'fontsize',20,'yscale','log')
ylabel('Imag\{DMI\}','fontweight','bold','interpreter','latex','FontSize',35)
xlabel('Frequency(Hz)','fontweight','bold','interpreter','latex','FontSize',35)
legend({'Measured','Simulated','Calibrated'},'fontweight','bold','interpreter','latex','FontSize',20,'location','southeast')
saveas(calibration_10K,'Pictures\calibration_10K', 'fig')
saveas(calibration_10K,'Pictures\calibration_10K', 'png')
close(calibration_10K)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%    Plot Test
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calibration_test = figure('units','normalized','outerposition',[0 0 1 1]);
% subplot(2,1,1)
% semilogx(frequencies,real(inductance_differential_test_meas),'k','linewidth',3)
% hold on 
% scatter(frequencies,real(inductance_differential_test_sim),100,'k','linewidth',3)
% hold on
% semilogx(frequencies,real(inductance_differential_test_cal),'k-.','linewidth',3)
% grid on
% set(gca,'fontsize',20,'yscale','log')
% ylabel('Real\{DMI\}','fontweight','bold','interpreter','latex','FontSize',35)
% xlabel('Frequency(Hz)','fontweight','bold','interpreter','latex','FontSize',35)
% legend({'Measured','sim','Calibrated'},'fontweight','bold','interpreter','latex','FontSize',20)
% subplot(2,1,2)
% semilogx(frequencies,imag(inductance_differential_test_meas),'r','linewidth',3)
% hold on 
% scatter(frequencies,imag(inductance_differential_test_sim),100,'r','linewidth',3)
% hold on
% semilogx(frequencies,imag(inductance_differential_test_cal),'r-.','linewidth',3)
% grid on
% set(gca,'fontsize',20,'yscale','log')
% ylabel('Imag\{DMI\}','fontweight','bold','interpreter','latex','FontSize',35)
% xlabel('Frequency(Hz)','fontweight','bold','interpreter','latex','FontSize',35)
% legend({'Measured','sim','Calibrated'},'fontweight','bold','interpreter','latex','FontSize',20)
% saveas(calibration_test,'Pictures\calibration_test', 'fig')
% saveas(calibration_test,'Pictures\calibration_test', 'png')
% close(calibration_test)




end
