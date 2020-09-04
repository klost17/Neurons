%% Hand-in 2 - Section e)

% Authors: Alexandre Justo, Núria Marzo, Pablo Navarro, Pau Riera.

% For this section, it is required the same script as the one used in
% section d), just modifying that now the number of pulses per second is
% fixed at 49 (Alexandre's ID), and the strength of the pulse is not a
% fixed value anymore, but it varies. Simulation is run 100 times.

% WARNING: when pressing Run, you are likely to wait for about 1h 30min.

%% e)

Icmin = 0.1; Icmax = 12.5; Icvec = Icmin:0.1:Icmax;
m = 49;
reps = 100;
numspikes = zeros(reps,numel(Icvec));

for iii = 1:numel(Icvec)
    Ic = Icvec(iii);
    
    for num = 1:reps
        t = linspace(0,1000,100000); ii = 0; 
        Ta = 1000*rand(1,m);
    
        for kk = 1:m
            Is = Ic*heaviside(t-Ta(kk)).*exp(-(t-Ta(kk))/10); 
            ii = ii+Is;
        end
    
        k = 1e-2; N = 1000/k;
        g1 = 120; g2 = 34; g3 = 0.33; V1 = 51; V2 = -75; V3 = -55; 
        Vo = 20; no = 0.5; mo = 0.5; ho = 0.5; yo = [Vo no mo ho];
        u = yo(1)+65; v = zeros(length(yo), N);
    
        dV = ii(1)-g1*yo(3)^3*yo(4)*(yo(1)-V1)-g2*yo(2)^4*(yo(1)-V2)...
            -g3*(yo(1)-V3);
        dn = 0.01*(10-u)/(exp(1-0.1*u)-1)*(1-yo(2))-yo(2)*0.125*exp(-u/80);
        dm = 0.1*(24-u)/(exp(2.4-0.1*u)-1)*(1-yo(3))-yo(3)*4*exp(-u/17);
        dh = 0.07*exp(-u/20)*(1-yo(4))-yo(4)/(1+exp(3-0.1*u));
        f0 = [dV dn dm dh]; a = yo';
   
        for q=1:(N-1)
            v(:,q) = yo' + k*f0'; 
            yo = v(:,q); yo = yo'; u = yo(1)+65;
        
            dV = ii(q+1)-g1*yo(3)^3*yo(4)*(yo(1)-V1)...
                -g2*yo(2)^4*(yo(1)-V2)-g3*(yo(1)-V3);
            dn = 0.01*(10-u)/(exp(1-0.1*u)-1)*(1-yo(2))...
                -yo(2)*0.125*exp(-u/80);
            dm = 0.1*(24-u)/(exp(2.4-0.1*u)-1)*(1-yo(3))...
                -yo(3)*4*exp(-u/17);
            dh = 0.07*exp(-u/20)*(1-yo(4))-yo(4)/(1+exp(3-0.1*u));
            f0 = [dV dn dm dh];
        end
    
        v = [a v];

        spikes = findpeaks(v(1,1:(end-1)));
        numspikes(num,iii) = numel(find(spikes>-40));
    end
    
end

figure('units','normalized','outerposition',[0 0 1 1])
errorbar(Icvec,mean(numspikes),var(numspikes),var(numspikes))
grid minor
axis([Icmin Icmax 0 100])
title(['Average number of spikes per second, with error bars indicating variance (n=' num2str(m) ')'])
xlabel('I_c (mA)')
ylabel('# Spikes')