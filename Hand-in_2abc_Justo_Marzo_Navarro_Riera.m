%% Hand-in 2 - Section a) b) c)

% Authors: Alexandre Justo, Núria Marzo, Pablo Navarro, Pau Riera.

% In this script we'll use the ODEs given by the HH model to predict the
% behaviour of the system given some initial conditions and a constant
% intensity, which will be 0 for sections a) and b).

%% a)

% We want to prove that the system arrives to a steady state (shown in
% the graphs).

% We set the initial conditions:
Vo=0.5; no=0.5; mo=0.5; ho=0.5; Tmax=10000; Io=0;
yo=[Vo no mo ho Io];

% odefun_V encodes the ODEs of the HH model. ode23 is used to solve the
% system from 0 to Tmax, where time is in ms:
[T,Y] = ode23(@odefun_V,[0 Tmax],yo);

% Now we show the results, we plot only the first 200 components of the
% solution vector for clarity (A). However, the equilibrium value
% displayed is the last one to improve accuracy (B).

% Equilibrium state of the transmembrane potential in mV (B).
Vf = Y(end,1)

% Equilibrium value of the probability of having the sub-unit n open.
nf = Y(end,2)

% Equilibrium value of the probability of having the sub-unit m open.
mf = Y(end,3)

% Equilibrium value of the probability of having the sub-unit h open.
hf = Y(end,4)

%% b)

f=100; % Factor we use to properly see the probabilities (0-100)
figure(1)
plot(T(1:200),Y(1:200,1))
hold on
plot(T(1:200),f*Y(1:200,2))
hold on
plot(T(1:200),f*Y(1:200,3))
hold on
plot(T(1:200),f*Y(1:200,4))
grid on
xlabel('Time [ms]')
legend('Membrane potential','Probability(subunit_n=open)[%]',...
    'Probability(subunit_m=open)[%]','Probability(subunit_h=open)[%]')
hold off

%% c)

% We show how for different values of the intensity we aproach the exact
% value at which we have periodic oscillations of the membrane potential.
I = [1.054 1.0541 1.05418];

for j=1:3
    
    Io = I(j); 
    yo = [Vf nf mf hf Io]; % Using the steady state as initial conditions
    [T,Y] = ode23(@odefun_V,[0 Tmax],yo);
    V_t = Y(:,1);
    
    figure(1+j)

    plot(T,V_t)
    xlabel('Time [ms]')
    ylabel('Membrane potential [mV]')
    title(['Membrane potential for I =',num2str(Io)])
    grid on
end

%% Functions:

% The function that encodes the HH model:

function dy = odefun_V(t,y)

% Constants:

g1 = 120; g2 = 34; g3 = 0.33; V1 = 51; V2 = -75; V3 = -55;

% Variables:
 
V = y(1); n = y(2); m = y(3); h = y(4); u = V+65; 

alpha_n = 0.01*(10-u)/(exp(1-(0.1*u))-1);
beta_n = 0.125*exp(-u/80);

alpha_m = 0.1*(24-u)/(exp(2.4-(0.1*u))-1);
beta_m = 4*exp(-u/17);

alpha_h = 0.07*exp(-u/20);
beta_h = 1/(exp(3-(0.1*u))+1);

% I(t):

I = y(5);

% ODEs:

dV = I-g1*(m^3)*h*(V-V1)-g2*(n^4)*(V-V2)-g3*(V-V3);
dn = alpha_n*(1-n)-beta_n*n;
dm = alpha_m*(1-m)-beta_m*m;
dh = alpha_h*(1-h)-beta_h*h;
dI = 0;

dy = [dV; dn; dm; dh; dI];
end