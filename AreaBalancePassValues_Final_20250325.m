function [maxEpsValue, plotArray, t] = AreaBalancePassValues_Final_20250325(molarity, t_d, k5_max, initial_v_T, fEis, eps_Ca, Caout, exoOn, endoOn, plotting)
global y CaOverTime count

maxTime = 60*10;

timestep = 0.01;
t = 0.0:timestep:maxTime;
t = t';

% height of gap between slides in m
gap_height = 4.5e-6;
small_r = gap_height*0.5;

% del molarity in mol/m^3
molarity = molarity*1000;

Berro_ratio = 0.5E-3 / 1.08E-5;
% stretching modulus (N/m)
kappa = 2E-3*Berro_ratio;

fOsm = 0.7;

eps_threshold = 0.01;

Pslow = 4.879E-9;
Pfast = 2.311E-7;

% kT at 300 K
kT = 300*1.3806e-23;
Na = 6.022e23;

% molar volume of water (m^3/mol)
vs = 1.803e-5;


initial_R = sqrt((pi^2-32/3)*small_r^4 + 8*small_r*initial_v_T/pi)/4/small_r - pi*small_r/4;

initial_Area = 2*pi^2*small_r*initial_R +4*pi*small_r^2 + 2*pi*initial_R^2;
initial_Volume = pi^2*small_r^2*initial_R + 4/3*pi*small_r^3 + 2*pi*small_r*initial_R^2;

initial_sigma = 4.973782828322014e-04;
initial_eps = initial_sigma / kappa;

% inital concentration drop across cell at steady state before
% concentration perturbation [mol/m^3]
delC_0 = initial_sigma / (kT * Na / (1/small_r + 1/(initial_R+small_r)));

initial_A_0 = initial_Area / (initial_eps + 1.0);

% initial values of
% [1] protoplast radius (m)
% [2] del molarity (mol/m^3)
% [3] Aves (m^3)
% [4] A0 (m^3)
% [5] Aeis (m^3)
% [6] Asterol (m^3)
% [7] Ca (mo/m^3)

initial_y=[initial_R molarity 0.0 0.0 0.0 9E-12 0.0];


eps_c = 0.01;
eps_s = eps_c*3;


beta = 3.08E-10;
gamma = 1.4E-15;


% Ca leakage mechanosensitive
betaCa = 7E-11;
gammaCa = 7E-12;


%um^2/min
maxExoRatePerMinute = 8.0;

V_T = pi^2*small_r^2*initial_y(1) + 4/3*pi*small_r^3 + 2*pi*small_r*initial_y(1)^2;

fAves = 0.3;

initial_y(3) = fAves*initial_Area;

initial_y(4) = initial_A_0;

initial_y(5) = fEis*initial_Area;

initial_y(6) = 0.09*initial_Area;

format long

% initial rates of exo and endocytosis, respectively (1/s)
% new values of exo/endocytosis initial are taken from Gerganova et al. Sci
% Adv 2021
r_exo_WT = 1.284e-12 / 60.0;
r_endo_WT = 0.942e-12 / 60.0;

k10 = r_exo_WT / initial_y(3);


if  initial_y(3) <= 0.0
    k10 = 0.0;
end

k20 = r_endo_WT / initial_y(4);
if  initial_y(4) <= 0.0
    k20 = 0.0;
end

%BFA
if exoOn == false
    k10 = 0.0*k10;
    k20 = 0.0*k20;
end

%LatA
if endoOn == false
    k20 = 0.0*k20;
    k10 = 0.5*k10;
end

% initial rates of eisosome removal and reassembly (1/s)
k30 = 0.05E-12 /60 / initial_y(5);
if  initial_y(5) <= 0.0
    k30 = 0.0;
end

k40 = k30 * initial_y(5) / initial_y(4);


% Ca dependence of exo
lambda1 = 5E3*1.0;
lambda1a = 0.5*lambda1*0.3;

lambda1b = 100;

% stress dependence of endocytosis
lambda2 = 5.0E1*0.3;

% stress dependence of eisosome removal and formation
lambda3 = 200.0*0.3;
lambda4 = 50.0*0.3;


initialCin = 1.2*1000;
Cout = initialCin - molarity;


burstTime = 0.0;

Ca_thres = 1.0;

CaOverTime = [];
plotArray = [];

num_timesteps = length(t);
num_entries_plotArray = 25;

plotArray = zeros(num_timesteps, num_entries_plotArray);
count = 1;

% solve ode using fixed step solver
[y] = ode1(@rhs, t, initial_y);



maxEpsValue = max(plotArray(:,12));



final_delC = plotArray(end,3)*1000;
final_internal_C = Cout + final_delC;

final_V_T = pi^2*small_r^2*(plotArray(end,1)/1E6) + 4/3*pi*small_r^3 + 2*pi*small_r*(plotArray(end,1)/1E6)^2;

totalMolRemoved = initialCin*initial_Volume - final_internal_C*final_V_T;

totalMolarityRemoved = totalMolRemoved / final_V_T;


%plotting
if plotting == true
figure(1)
clf()
plot(t/60, (plotArray(:,1) + small_r*1E6)*2,'b')


hold on;

plot(x_data_WT, y_data_WT, 'ko');


hold off;
xlabel('t (min)')
ylabel('Medial D [{\mu}m]')


set(findall(gca, 'Type', 'Line'),'LineWidth',2);

figure(2)
clf()
plot(t/60, plotArray(:,3))
xlabel('t (min)')
ylabel('\Delta c')
set(findall(gca, 'Type', 'Line'),'LineWidth',2);

figure(3)
clf()
plot(t/60, plotArray(:,4))
xlabel('t (min)')
ylabel('A_{Cyto} [{\mu}m^2]')

set(findall(gca, 'Type', 'Line'),'LineWidth',2);

figure(4)
clf()
plot(t/60, plotArray(:,5))
xlabel('t (min)')
ylabel('A_{Eis} [{\mu}m^2]')

set(findall(gca, 'Type', 'Line'),'LineWidth',2);

figure(5)
clf()

plot(t/60, plotArray(:,6))
xlabel('t (min)')
ylabel('\sigma')

set(findall(gca, 'Type', 'Line'),'LineWidth',2);

figure(6)
clf()
plot(t/60, plotArray(:,8), t/60, plotArray(:,9))
legend('Exocytosis, k1*Aves', 'Endocytosis, k2*A0')
xlabel('t (min)')
ylabel('Rate [{\mu}m^2/min]')
set(findall(gca, 'Type', 'Line'),'LineWidth',2);

figure(7)
clf()
plot(t/60, plotArray(:,10), t/60, plotArray(:,11))
legend('Eisosome Destruction, k3*Aeis', 'Eisosome Creation, k4*A0')
xlabel('t (min)')
ylabel('Rate [{\mu}m^2/min]')
set(findall(gca, 'Type', 'Line'),'LineWidth',2);


figure(8)
plot(t/60, plotArray(:,12))
xlabel('t (min)')
ylabel('\epsilon')

set(findall(gca, 'Type', 'Line'),'LineWidth',2);

figure(9)
plot(t/60, plotArray(:,13))
xlabel('t (min)')
ylabel('c_{in}')


set(findall(gca, 'Type', 'Line'),'LineWidth',2);

figure(10)
plot(t/60, plotArray(:,14))
xlabel('t (min)')
%ylabel('Permeability / P_{0}')
ylabel('Permeability [m/s]')

set(findall(gca, 'Type', 'Line'),'LineWidth',2);

figure(11)
plot(t/60, plotArray(:,15), t/60,plotArray(:,16))
xlabel('t (min)')
ylabel('Rates')
legend('Exocytosis, k1', 'Endocytosis, k2')

set(findall(gca, 'Type', 'Line'),'LineWidth',2);

figure(12)
plot(t/60, plotArray(:,17), t/60,plotArray(:,18))
xlabel('t (min)')
ylabel('Rates')
legend('Eisosome Destruction, k3', 'Eisosome Creation, k4')

set(findall(gca, 'Type', 'Line'),'LineWidth',2);

figure(13)
plot(t/60, plotArray(:,19), t/60,plotArray(:,2))
xlabel('t (min)')
ylabel('Surface Area [{\mu}m^2]')
legend('A', 'A0')

set(findall(gca, 'Type', 'Line'),'LineWidth',2);

figure(14)
plot(t/60, plotArray(:,7))
xlabel('t (min)')
ylabel('Net Area Addition Rate [{\mu}m^2/min]')

set(findall(gca, 'Type', 'Line'),'LineWidth',2);

figure(15)
plot(t/60, plotArray(:,20))
xlabel('t (min)')
ylabel('Ca')

set(findall(gca, 'Type', 'Line'),'LineWidth',2);


figure(16)
plot(t/60, plotArray(:,21))
xlabel('t (min)')
ylabel('k5 [{\mu}m^2/min]')

set(findall(gca, 'Type', 'Line'),'LineWidth',2);


figure(17)
plot(t/60, plotArray(:,22))
xlabel('t (min)')
ylabel('j_mech [mol / m^2/ s]')

set(findall(gca, 'Type', 'Line'),'LineWidth',2);

dvdt = pi^2*small_r^2.*plotArray(:,23) + 0 + 4*pi*small_r*plotArray(:,1).*plotArray(:,23);
figure(18)
plot(plotArray(:,3),  dvdt)
xlabel('\Delta c')
ylabel('dV/dt [m^3/s]')


set(findall(gca, 'Type', 'Line'),'LineWidth',2);

figure(19)
plot(t/60, plotArray(:,24))
xlabel('t (min)')
ylabel('j_mech_Ca [mol / um^2 s]')

set(findall(gca, 'Type', 'Line'),'LineWidth',2);

figure(20)
plot(t/60, plotArray(:,25))
xlabel('t (min)')
ylabel('j_act_Ca [mol / um^2 s]')

set(findall(gca, 'Type', 'Line'),'LineWidth',2);


figure(21)
Area_values = 2.0*pi^2*small_r*1E6.*plotArray(:,1)  + 4*pi*(small_r*1E6)^2 + 2*pi.*plotArray(:,1) .* plotArray(:,1);
Area_values_normalized = Area_values/Area_values(1);
plot(t, Area_values_normalized, 'b');

end


    % differential equation definition
    function dydt = rhs(t,y)
        
        dydt = zeros(7,1);
        
        A_T = 2*pi^2*small_r*y(1) + 4*pi*small_r^2 + 2*pi*y(1)^2;
        V_T = pi^2*small_r^2*y(1) + 4/3*pi*small_r^3 + 2*pi*small_r*y(1)^2;
        
        eps = (A_T - y(4)) / y(4);
        
        
        sigma = kappa*eps;


        P = Pslow;
        
        if(eps > eps_threshold)
            P = Pfast;
        end

        
        if(eps < 1000.0)
            constant = P*vs/kT/Na * A_T / (pi^2*small_r^2 + 4.0*pi*small_r*y(1))/fOsm;

            firstTerm = (kT*Na*y(2)) * constant;
            secondTerm = (sigma * (1/small_r + 1/(small_r + y(1)))) * constant;

            A_growth = (r_exo_WT - r_endo_WT);
            V_growth = (A_T^0.5) * A_growth / (4*pi^0.5);

            % differential equation for proplast radius
            dydt(1) = (firstTerm - secondTerm) + V_growth / (pi^2*small_r^2 + 4.0*pi*small_r*y(1));

            fsigma = 0;            
            
            if eps > eps_c
                if eps > eps_s
                    fsigma = (eps_s - eps_c);
                else
                    fsigma = (eps - eps_c);
                end   
            end
            
            
            j_mech = -beta*kT*Na*y(2)*fsigma / fOsm;
            j_act = gamma*kT*Na * (delC_0 - y(2)) / fOsm;

            %differential equation for delta_c
            dydt(2) = (-(pi^2*small_r^2+4*pi*y(1)*small_r)*dydt(1)*(y(2)+Cout) + A_T*(j_mech+j_act)) / V_T;
            
            %balance on Ca (7)
            fsigma3 = 0;            
            
            
            if eps > eps_Ca
                fsigma3 = 1.0;
            end
            

            
            j_mech_Ca = -betaCa*kT*Na*(y(7)-Caout)*fsigma3;
            j_act_Ca = -gammaCa*kT*Na * (y(7));
            
            
            %differential equation for Ca
            dydt(7) = (-(pi^2*small_r^2+4*pi*y(1)*small_r)*dydt(1)*y(7) + A_T*(j_mech_Ca+j_act_Ca)) / V_T / fOsm;
            

            
            CaOverTime = [CaOverTime, y(7)];
            

            k1 = k10;
            

            if(length(CaOverTime) >= t_d / timestep)
                k1 = 1 + (lambda1a*eps)*(CaOverTime(1));
                
                k1 = k1 + lambda1b*(eps-eps_threshold);

                
                k1 = k10 * k1;
                
            
                CaOverTime = CaOverTime(2:end);
            end

            
            if(k1*y(3)/1E-12*60 > maxExoRatePerMinute)
                k1 = maxExoRatePerMinute/(y(3)/1E-12*60);
            end
            
            if(k1 < 0)
                k1 = 0.0;
            end
            
            
            k2 = k20 * (1 - lambda2*(eps - initial_eps));
            
            if(k2 < 0)
                k2 = 0.0;
            end

            k3 = k30 * (1 + lambda3 *(eps - initial_eps));

            k4 = k40 * (1 - lambda4 *(eps - initial_eps));
            if(k4 < 0)
                k4 = 0.0;
            end

            % differential equation for Aves
            dydt(3) = -k1*y(3) + k2*y(4);
            
          
            k5 = 0.0;
            
            if(y(7) > Ca_thres)
                k5 = k5_max;
            end
            
            % differential equation for A0
            dydt(4) = k1*y(3) - k2*y(4) + k3*y(5) - k4*y(4) + (k5*1E-12/60.0);
            
            % differential equation for Aeis
            dydt(5) = -k3*y(5) + k4*y(4);
            
            dydt(6) = 0.0;
           
            % medial diameter (um), A0 (um^2), 
            tmpArray = [y(1)*1E6, y(4)/1e-12, y(2) / 1000.0, y(3)/(10^-12), y(5) / (10^-12), sigma, ...
                60*(k1*y(3)/(10^-12) - k2*y(4)/1e-12 + k3*y(5) / (10^-12) - k4*y(4)/1e-12), k1*y(3)/(10^-12) * 60,  ...
                k2*y(4)/1e-12*60, k3*y(5) / (10^-12)*60, k4*y(4)/1e-12*60, eps, Cout + y(2), P, k1, k2, k3, k4, ...
                (2.0*pi^2*small_r*y(1)+4*pi*small_r^2+2*pi*y(1)^2) *1E12, y(7), k5, j_mech, dydt(1), j_mech_Ca, j_act_Ca];
            plotArray(count,:) = tmpArray;

            count = count + 1;
        else
            if(burstTime == 0.0)
                burstTime = t;
            end
        end
    end
end