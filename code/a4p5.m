%Jay Maini 101037537
set(0, 'DefaultFigureWindowStyle','docked')
clear all
close all
global G C F

%Circuit parameters
R1 = 1;
R2 = 2;
R3 = 418;
R4 = 0.1;
Ro = 1000;

L1 = 0.2;
a = 100;
C1 = 0.25;
G1 = 1/R1;
G2 = 1/R2;
G3 = 1/R3;
G4 = 1/R4;
Go = 1/Ro;

%Bandwidth limiting capacitor
Cn = 0.00001;
%Cn = 0;

C = [0 0 0 0 0 0 0;
    -C1 C1 0 0 0 0 0;
    0 0 -L1 0 0 0 0;
    0 0 0 Cn 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 Cn 0 0 0;
    0 0 0 0 0 0 0];

G = [1 0 0 0 0 0 0;
    -G2 G1+G2 -1 0 0 0 0;
    0 1 0 -1 0 0 0;
    0 0 -1 G3 0 0 0;
    0 0 0 0 -a 1 0;
    0 0 0 G3 -1 0 0;
    0 0 0 0 0 -G4 G4+Go];

F = [1;0;0;0;0;0;0];

%[V1 V2 IL V3 I3 V4 Vo]

%Time step
time = 1000;

%Initialize three input signals 
Vin = zeros(1,time);
f = 1/(0.03*time);

V = zeros(7,1);

In = 0.01*rand(time,1);
figure
%Gaussian Input
for i=1:time
    %gaussian parameters
    a = 1;
    b = 0.06*time;
    c = 0.03*time;
    Vin(i) = a*exp(-(i-b)^2/(2*c^2));
    
    %recalculate F vector
    F = [Vin(i);0;0;0;In(i);0;0];

    A = C/0.001 + G;
   V_p = V;
   V = A\((C*V_p/0.001) + F);
   Vo(i) = V(7);

    if (i >= 2)
        plot([i i-1], [Vin(i) Vin(i-1)], 'r')
        hold on
        plot([i i-1], [Vo(i) Vo(i-1)], 'b')
    end
end
title('Gaussian Input with std. dev 30ms and delay 60ms with Noise')
ylabel('Voltage')
xlabel('Time')

figure
%Gaussian Input - Zoomed in result
for i=1:time
    %gaussian parameters
    a = 1;
    b = 0.06*time;
    c = 0.03*time;
    Vin(i) = a*exp(-(i-b)^2/(2*c^2));
    
    %recalculate F vector
    F = [Vin(i);0;0;0;In(i);0;0];

    A = C/0.001 + G;
   V_p = V;
   V = A\((C*V_p/0.001) + F);
   Vo(i) = V(7);

    if (i >= 2)
        plot([i i-1], [Vin(i) Vin(i-1)], 'r')
        hold on
        plot([i i-1], [Vo(i) Vo(i-1)], 'b')
    end
end
title('Close-up of previous graph, showing Noise')
ylabel('Voltage')
xlabel('Time')

xlim([142.8 244.3])
ylim([-0.068 0.054])

%Gaussian with noise - FFT plot
In = 0.1*rand(time,1);
for i=1:time
    %gaussian parameters
    a = 1;
    b = 0.06*time;
    c = 0.03*time;
    Vin(i) = a*exp(-(i-b)^2/(2*c^2));
    
    %recalculate F vector
    F = [Vin(i);0;0;0;In(i);0;0];

    A = C/0.001 + G;
   V_p = V;
   V = A\((C*V_p/0.001) + F);
   Vo(i) = V(7);

%     if (i >= 2)
%         plot([i i-1], [Vin(i) Vin(i-1)], 'r')
%         hold on
%         plot([i i-1], [Vo(i) Vo(i-1)], 'b')
%     end
end

figure
shrek = fftshift(fft(Vin));
plot(abs(shrek), 'r')

hold on
shrek2 = fftshift(fft(Vo));
plot(abs(shrek2), 'b')
xlim([460 540])
title('Freq response of Gaussian with Noise')
ylabel('Voltage')
xlabel('Time')

%recalculate Cn with different values

Cn = 0.0001;

C = [0 0 0 0 0 0 0;
    -C1 C1 0 0 0 0 0;
    0 0 -L1 0 0 0 0;
    0 0 0 Cn 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 Cn 0 0 0;
    0 0 0 0 0 0 0];

figure
%Gaussian Input - Zoomed in result, Cn LARGE ------------------
for i=1:time
    %gaussian parameters
    a = 1;
    b = 0.06*time;
    c = 0.03*time;
    Vin(i) = a*exp(-(i-b)^2/(2*c^2));
    
    %recalculate F vector
    F = [Vin(i);0;0;0;In(i);0;0];

    A = C/0.001 + G;
   V_p = V;
   V = A\((C*V_p/0.001) + F);
   Vo(i) = V(7);

    if (i >= 2)
        plot([i i-1], [Vin(i) Vin(i-1)], 'r')
        hold on
        plot([i i-1], [Vo(i) Vo(i-1)], 'b')
    end
end
title('Gaussian Input with larger Cn=0.0001')
ylabel('Voltage')
xlabel('Time')

Cn = 0.0000001;

C = [0 0 0 0 0 0 0;
    -C1 C1 0 0 0 0 0;
    0 0 -L1 0 0 0 0;
    0 0 0 Cn 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 Cn 0 0 0;
    0 0 0 0 0 0 0];

figure
%Gaussian Input - Zoomed in result, Cn SMALL ------------------
for i=1:time
    %gaussian parameters
    a = 1;
    b = 0.06*time;
    c = 0.03*time;
    Vin(i) = a*exp(-(i-b)^2/(2*c^2));
    
    %recalculate F vector
    F = [Vin(i);0;0;0;In(i);0;0];

    A = C/0.001 + G;
   V_p = V;
   V = A\((C*V_p/0.001) + F);
   Vo(i) = V(7);

    if (i >= 2)
        plot([i i-1], [Vin(i) Vin(i-1)], 'r')
        hold on
        plot([i i-1], [Vo(i) Vo(i-1)], 'b')
    end
end
title('Gaussian Input with smaller Cn=0.0000001')
ylabel('Voltage')
xlabel('Time')


figure
shrek = fftshift(fft(Vin));
plot(abs(shrek), 'r')

hold on
shrek2 = fftshift(fft(Vo));
plot(abs(shrek2), 'b')
xlim([460 540])
title('Freq response of Gaussian with Noise')
ylabel('Voltage')
xlabel('Time')



%Different TIME STEPS ------------------------------------------

%Time step
time = 500;

%Initialize three input signals 
Vin = zeros(1,time);
f = 1/(0.03*time);

V = zeros(7,1);

In = 0.001*rand(time,1);
figure
subplot(1,2,1)
%Gaussian Input
for i=1:time
    %gaussian parameters
    a = 1;
    b = 0.06*time;
    c = 0.03*time;
    Vin(i) = a*exp(-(i-b)^2/(2*c^2));
    
    %recalculate F vector
    F = [Vin(i);0;0;0;In(i);0;0];

    A = C/0.001 + G;
   V_p = V;
   V = A\((C*V_p/0.001) + F);
   Vo(i) = V(7);

    if (i >= 2)
        plot([i i-1], [Vin(i) Vin(i-1)], 'r')
        hold on
        plot([i i-1], [Vo(i) Vo(i-1)], 'b')
    end
end
title('Gaussian Input with 2ms time step')
ylabel('Voltage')
xlabel('Time')


%Time step
time = 2000;

%Initialize three input signals 
Vin = zeros(1,time);
f = 1/(0.03*time);

V = zeros(7,1);

In = 0.001*rand(time,1);
%Gaussian Input
subplot(1,2,2)
for i=1:time
    %gaussian parameters
    a = 1;
    b = 0.06*time;
    c = 0.03*time;
    Vin(i) = a*exp(-(i-b)^2/(2*c^2));
    
    %recalculate F vector
    F = [Vin(i);0;0;0;In(i);0;0];

    A = C/0.001 + G;
   V_p = V;
   V = A\((C*V_p/0.001) + F);
   Vo(i) = V(7);

    if (i >= 2)
        plot([i i-1], [Vin(i) Vin(i-1)], 'r')
        hold on
        plot([i i-1], [Vo(i) Vo(i-1)], 'b')
    end
end
title('Gaussian Input with 0.5ms time step')
ylabel('Voltage')
xlabel('Time')

