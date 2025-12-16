
%%% PROBLEM 1   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Analytical solution (calculated exactly exactly)

analytical_Change_in_internal_energy   =  700*(820-320) + ((0.35/2)*(820.^2 - 320.^2)) - (((2*10.^-4)/3)*(820.^3 - 320.^3));
% energy_change = 4.151766666666667e+05;



% solution using integeral() function

cv = @(T)700 + 0.35*T - 2*10.^-4*T.^2;
integeral_energy_change = integral(cv,320,820);
% eneryg_change = 4.151766666666667e+05


 
 % solution using trapz() function

T = linspace(320,820,1000);
cv_new = cv(T);
trapz_energy_change = trapz(T,cv_new);
% energy_change = 4.151766624916542e+05


energy_change = 4.151766666666667e+05;

for i=2:1000
    
    T_new = linspace(320,820,i);
    cv_min = cv(T_new);
    change_in_energy = trapz(T_new,cv_min);
    error = abs(change_in_energy - energy_change)/energy_change;
    if error < 0.001
        minimum_grid = i;
       break; 
    end
    
    
end

compare1 = analytical_Change_in_internal_energy-integeral_energy_change;
disp(compare1)
compare2 = analytical_Change_in_internal_energy-trapz_energy_change;
disp(compare2)
disp(minimum_grid)
disp(analytical_Change_in_internal_energy)
disp(integeral_energy_change)
disp(trapz_energy_change)







%%% PROBLEM 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% All the equations
internal_energy = @(T)450 + 1.1*T + 0.0012*T.^2;   % internal energy equation
q_ext = @(t)5000*exp(-0.002*t);   % external heat addition
r_heat = @(T)1500*(1-exp(-0.01*T));  % reaction heat addition
der_int_energy = @(T)1.1 + 0.0024*T;  % derivative of internal enegy equation with respect to T (temperature)
energy_balance = @(t,T)q_ext(t) + r_heat(T);  %  energy balance equation
ode = @(t,T) q_ext(t) + r_heat(T)./der_int_energy(T);   % ODE as a function of (t,T)

tspan = [0 4000];     % interval in which we have to integrate
[t,T] = ode45(ode,tspan,300);  %  solving the ODE using ode45
internal_energy_value = internal_energy(T);   %  value of internal energy at all values of T(temperature) obtained from ODE
surpass_index = find((q_ext(t)< r_heat(T)) ,1, 'first');  % finding the index at which reaction heat contribution is surpassing external heating

q_values = q_ext(t);     % external heating contribution at all values of t(time)
r_values = r_heat(T);    % reaction heating contribution at all values of T (temperature)
surpass_t = t(surpass_index);  %  time when reaction heat contribution is surpassing 
surpass_T = T(surpass_index);  %  Temperature when reaction heta contribution is surpassing

% Plotting the figures;

%figure;
subplot(3,1,1);   % this figure recovers the Temperature history 
plot(t,T,'LineWidth',1.2)
xlabel('t(sec)'); 
ylabel('T(kelvin)');
title('Temperature(history)')
grid on;

subplot(3,1,2);    % this figure gives the evolution of the internal energy
plot(T,internal_energy_value,'LineWidth',1.2);
xlabel('T(Kelvin)');
ylabel('Internal_energy(kJ/kg)');
title('internal_energy evolution');
grid on;

subplot(3,1,3);     % tis figures clearly shows when the reaction heat is surpassing the external heat contribution
plot(t,q_values,'b-',t,r_values,'r--','LineWidth',1.2);
xlabel('t(sec)');
ylabel('KJ/sec');
title('reaction heat surpassing external heating');
grid on;







%%%  PROBLEM 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cp = 1.05;   % KJ/Kg.K
R = 0.287;   
T1 = 300;  % kelvin       %% initial guss for temperature

T = @(P)T1*P.^0.21875;   % temperature as a function of pressure

Z = @(P)1 + 0.0008*P - (120./T(P));  % analytic model for real gas
p1 = 1; p2 = 20;     % pressure range for integration


ds_real = @(P)Cp*(0.21875./P) - R*(Z(P)./P);   % infinitesimal entropy for real gas
entropy_change_real = integral(ds_real,p1,p2);  % entropy chage for real gas

ds_ideal = @(P)Cp*(0.21875./P) - (R./P);      % infinitesimal entropy change for ideal gas
entropy_change_ideal = integral(ds_ideal,p1,p2);   % entropy change for ideal gas

% percentage deviation from ideal gas
percentage_deviation = ((entropy_change_real-entropy_change_ideal)/abs(entropy_change_ideal))*100;   

% Real gas influence
%%%  Because in real gas Z<1 it decrease the second negative term in the
%%%  integeral and the positive term temperature dominates which gives us
%%%  net positive entropy change






%%% PROBLEM 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


T1 = 310; T2 = 670; Tb = 300;  % temperatures are in kelvin


enthalpy = @(Te)300 + 2.5*Te + 0.0007*Te.^2;
entropy = @(Te)2.0*log(Te) + 0.001*Te;

%energy_range = mass_flow*(enthalpy(T2)-enthalpy(T1));
%S_gen = mass_flow*(entropy(T2)-entropy(T1)) - energy/Tb;


energy_range =linspace(20e3,100e3,1000);   % rate of energy in watt

mass_flow = energy_range/(Tb*(entropy(T2)-entropy(T1)));   % mass flow rate

figure;
plot(mass_flow,energy_range,'r--',LineWidth=1.2);
xlabel('mass flow (Kg/sec)')
ylabel('energy(watt)')
title('feasible energy and mass_flow region')
grid on;







%%% PROBLEM 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T0 = 298; % in kelvin
cp = @(T)1200 + 0.4*T - (1.2e-4)*T.^2;  % temperature dependent heat capacity
equation = @(T)cp(T)./T;    % integrand of the integral for entropy change
entropy_change = integral(equation,350,900);  % entropy change using integral function
exergy_destroyed1 = T0*0.02*entropy_change;   % exergy destruction with 2% irreversibility
exergy_destroyed2 = T0*0.10*entropy_change;   % exergy destruction with 10% irreversibility

irre = linspace(0,0.20,100);    % Spanning irreversibility from 0% to 20%  with 100 equal gaps
exergy_destroyed = T0*irre*entropy_change; % exrgy destruction between 0% to 20% 

figure;
plot(irre,exergy_destroyed,'g--')   % pllotting of figure for irreversibility and exegy destruction
xlabel('irreversibility');
ylabel('exergy_destroyed');
title('exergy destruction vs irreversibility level');
grid on;








%%% PROBLEM 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Th = @(t)900 - 300*exp(-0.0008*t);   % higher boundary teperature function of time
Tc = @(t)300 + 40*sin(0.002*t);    % lower boundary temperature function of time

effic = @(t)1 - (Tc(t)./Th(t));    % efficiency function of time

Qin = @(t)20000*(1 + 0.3*sin(0.003*t));   % heat input function of time in kW

power = @(t)effic(t).*Qin(t);  % power in kW 
work = integral(power,0,100);  % work in kJ
time = linspace(0,100,1000);   % time frpm 0s to 100s divided into 1000 equal points
power1 = power(time);        % power at all points of time
work1 = trapz(time,power1);  % finding work using trapz function  and all time grid with power at that points

entropy_gen = @(t)(Qin(t)./Th(t) - Qin(t)./Tc(t));  % rate of entropy generation as a function of time
time1 = linspace(0,100,1000);   % time grid spanning from 0s to 100s with 1000 equally sapced points
entropy_gen = entropy_gen(time1); % entropy gen rate at all ponts of time between 0s to 100s 


figure;      % plotting time vs rate of entropy gen
plot(time1,entropy_gen,'b-')
xlabel('time(sec)');
ylabel('entropy generation (KW/K')
title('entropy generation vs time')
grid on;







%%% PROBLEM 7 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = 1.25;
R = 0.287;
Pa = 1; Pb = 10;   %pressure in bar
TA = 300;    % in kelvin
TB = (10.^0.2)*TA;  % final temperature using polytropic relation

 U = @(T)500 + 0.8*T + (1.5e-3)*T.^2;  % internal energy

 change_in_U = U(TB)-U(TA);  %change in internal energy(it is state function therefore final - initial)
 W = (R./(1-m))*(TB-TA);   % work using polytropic relation and then integrating it
 H = change_in_U + W;    % heat in the process
 






 %%% PROBLEM 8 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 C = 2.5; N = 1000;  % in kJ/K
 q_total = 5e5;  % in kJ
 Tb = 300;   % in Kelvin
 tf = 2000;  % time length in second
 dt = tf/N;  % differential temperature
 Tf = linspace(0,tf,N);  % time grid from 0s to 2000s with N equal points in between 

 q_uni = q_total/tf * ones(1,N);  % uniform heating
 q_opt = linspace(0,2*q_total/tf,N); % optimal heating descretizing
 q_opt = q_opt * (q_total/sum(q_opt*dt)); % optimal heating profile 


 T_uni = T0 + cumsum(q_uni)*dt/C;  % uniform temperature evolution
 T_opt = T0 + cumsum(q_opt)*dt/C;  % optimal temperature evolution

 S_uni = sum((q_uni./T_uni - q_uni/Tb)*dt); % uniform entropy 
 S_opt = sum((q_opt./T_opt - q_opt/Tb)*dt); % optimal entropy

 figure; % plotting between time length vs (uniform heat and optimal heat)
 plot(Tf,q_uni,'b--',Tf,q_opt,'LineWidth',1.2);
 xlabel('Time (s)'); ylabel('q(t)');
 legend('Uniform','Optimal'); grid on;

 figure;  % plotting between time length vs (uniform temperature evolution and optimal temperature evolution)
 plot(Tf,T_uni,'r--',Tf,T_opt,'LineWidth',1.2);
 xlabel('Time (s)'); ylabel('Temperature (K)');
 legend('Uniform','Optimal'); grid on;





 

 %%% SOLUTION 9 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 A = 1e5;
 E = 45;  % kJ/mol
 R = 8.314e-3; % gas constant
 c = 1.8; % kj/K

 t_guess = 300;  % guess temperature 
 tspan = linspace(0,1000,1000); % time length from 0s to 1000s with 1000 equal points in between

 r = @(Temp)A*exp(-E./(R*Temp)); % modeled heat reaction
 differ_T = @(tim,Temp)(r(Temp)+ 2000*exp(-0.001*tim))./c; % differential equation as with respect to time(t)

 [tim,Temp] = ode45(differ_T,tspan,t_guess);  %solving differential equation using ode45 function

 noise_T = Temp.*(1 + 0.01*randn(size(Temp))); % addding noise to the temperature-time data

 figure;  % plotting time vs (noise temperature and temperature from ode)
 plot(tim,noise_T,'o',tim,Temp,'LineWidth',1.2);
 xlabel('Time (s)'); ylabel('Temperature (K)');
 legend('Noisy data','True');
 grid on;






 %%% SOLUTION 10 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Tc = 270;          % lower sink temperature in kelvin 
Th = 320;          %  higher source temperature in kelvin
alpha = 0.02;  % given
k = 50;    % given


COP = @(rp) 5.4*(1 - alpha*(rp-1).^2./rp); % real coefficient of performance
Wc  = @(rp) k*(sqrt(rp) - 1);    % compressor work


obj = @(rp) -COP(rp);  %objective function to maximize COP


lb = 1.1;  %lower bound
ub = 6;    % upper bound


rp0 = 2;   %guess pressure ratio


options = optimoptions('fmincon','Display','iter');  % optimizations options
[rp_opt, fval] = fmincon(obj,rp0,[],[],[],[],lb,ub,@entropy_constraint, options); % running the optimizations

COP_opt = COP(rp_opt);  % optimized COP at all rp optimized

rp = linspace(1.1,6,400); % pressure ratio length with 400 points in between 
COP_vals = COP(rp);  % values of COP at all pressure ratio values
Wc_vals  = Wc(rp);   % compressor work at all pressure ratio values

Sgen = Wc_vals .*(COP_vals.*(1/Th - 1/Tc) + 1/Tc);  % entropy generation 

feasible = Sgen <= 0.05;   % feasible region

% Plot 1: COP vs pressure ratio
figure;
plot(rp, COP_vals,'LineWidth',2); grid on;
xlabel('Pressure ratio r_p');
ylabel('COP');
title('COP vs Pressure Ratio');

% Plot 2: Entropy generation vs pressure ratio
figure;
plot(rp, Sgen,'LineWidth',2); hold on;
yline(0.05,'r--','Entropy limit');
grid on;
xlabel('Pressure ratio r_p');
ylabel('Entropy generation (kJ/K)');
title('Entropy Generation vs Pressure Ratio');

% Plot 3: Feasible region
figure;
plot(rp, COP_vals,'b','LineWidth',2); hold on;
plot(rp(~feasible), COP_vals(~feasible),'r','LineWidth',2);
grid on;
xlabel('Pressure ratio r_p');
ylabel('COP');
legend('Feasible','Infeasible');
title('Feasible Region under Entropy Constraint');


% function for energy constraint
function [c_new, ceq] = entropy_constraint(rp)

Tc = 270;
Th = 320;
alpha = 0.02;
k = 50;

COP = 5.4*(1 - alpha*(rp-1).^2./rp);
Wc  = k*(sqrt(rp) - 1);

Sgen = Wc*(COP*(1/Th - 1/Tc) + 1/Tc);

c_new = Sgen - 0.05;   
ceq = [];

end





























