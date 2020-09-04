%% Hand-in 2 - Section d)

% Authors: Alexandre Justo, Núria Marzo, Pablo Navarro, Pau Riera.

% Our intensity is not constant anymore. Actually we have a time-dependent
% intensity that varies exponentially.

%% d)

% Vectors for DNI numbers and time. ii serves as a reference for the
% future values of the intensity. z is defined to select our peaks later. 

Ic = 0.5; m = [2 23 238 2384]; % Pau's ID
t = linspace(0,1000,100000); ii = 0; z = [0 0 2.5 0];


for j = 1:length(m) % Loop for the different values of n.
    n = m(j); 
    Ta = 1000*rand(1,n); % We pick random times between 0 and 1000 ms.
    
    for l = 1:n
        Is = Ic*heaviside(t-Ta(l)).*exp(-(t-Ta(l))/10); 
        % Heaviside function prevents from getting values out of range.
        ii = ii+Is; % ii becomes I(t)
    end
    
    figure('units','normalized','outerposition',[0 0 1 1]) % Full screen.
    subplot(2,1,1)
    plot(t,ii)
    grid minor
    axis([0 1000 -25 25])
    title(['I(t) for n=' num2str(m(j))])
    xlabel('Time (ms)')
    ylabel('Intensity (mA)')
    
    % After computing I(t), let's integrate the system with a simple
    % Euler method. Initial values & parameters:
    
    k = 1e-2; N = 1000/k;
    g1 = 120; g2 = 34; g3 = 0.33; V1 = 51; V2 = -75; V3 = -55; 
    Vo = 20; no = 0.5; mo = 0.5; ho = 0.5; yo = [Vo no mo ho];
    u = yo(1)+65; v = zeros(length(yo), N);
    
    % The derivatives:
    
    dV = ii(1)-g1*yo(3)^3*yo(4)*(yo(1)-V1)-g2*yo(2)^4*(yo(1)-V2)...
        -g3*(yo(1)-V3);
    dn = 0.01*(10-u)/(exp(1-0.1*u)-1)*(1-yo(2))-yo(2)*0.125*exp(-u/80);
    dm = 0.1*(24-u)/(exp(2.4-0.1*u)-1)*(1-yo(3))-yo(3)*4*exp(-u/17);
    dh = 0.07*exp(-u/20)*(1-yo(4))-yo(4)/(1+exp(3-0.1*u));
    f0 = [dV dn dm dh]; a = yo'; % Saving in a the initial values.
   
    for q=1:(N-1)
        v(:,q) = yo' + k*f0'; 
        yo = v(:,q); yo = yo'; u = yo(1)+65; % Redefining initial values.
        
        dV = ii(q+1)-g1*yo(3)^3*yo(4)*(yo(1)-V1)-g2*yo(2)^4*(yo(1)-V2)...
            -g3*(yo(1)-V3);
        dn = 0.01*(10-u)/(exp(1-0.1*u)-1)*(1-yo(2))-yo(2)*0.125*exp(-u/80);
        dm = 0.1*(24-u)/(exp(2.4-0.1*u)-1)*(1-yo(3))-yo(3)*4*exp(-u/17);
        dh = 0.07*exp(-u/20)*(1-yo(4))-yo(4)/(1+exp(3-0.1*u));
        f0 = [dV dn dm dh];
    end
    
    v = [a v]; % Adding the original initial values to the solution.
    % Selecting our peaks that surpass -63 mV except for the 3rd case:
    [pks,locs] = findpeaks(v(1,1:(end-1)),'MinPeakHeight',-63+z(j));
    subplot(2,1,2)
    % Adjusting their lengths so that both can be represented:
    plot(t,v(1,1:(end-1)),t(locs),pks,'ok')
    grid minor
    axis([0 1000 -80 40])
    title(['V(t) for n=' num2str(m(j))])
    xlabel('Time (ms)')
    ylabel('Voltage (mV)')
    legend('Voltage Profile','Firing times')
end