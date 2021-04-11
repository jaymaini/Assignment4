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

C = [0 0 0 0 0 0 0;
    -C1 C1 0 0 0 0 0;
    0 0 -L1 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
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
time = 500;

%Initialize three input signals 
Vin = zeros(1,time);

f = 1/(0.03*time);
%v2 = sin(2*pi*f*i);
%v3 = 0;
V = zeros(7,1);
%Vout = zeros(1,1000);
figure
subplot(1,2,1)
%Step input -----------------------------------------------
for i=1:time
   if (i >= 0.03*time)
    Vin(i) = 1;
   end
   
   A = C/0.001 + G;
   V_p = V;
   V = A\((C*V_p/0.001) + F*Vin(i));
  % Vout(i) = V(7);
   Vo(i) = V(7);
   
   if (i >= 2)
        plot([i i-1], [Vin(i) Vin(i-1)], 'r')
        hold on
        plot([i i-1], [Vo(i) Vo(i-1)], 'b')
        
   end
    
end
title('Step Input with 0 -> 1 at 30ms')
xlabel('time (ms)')
ylabel('Voltage (V)')
legend('Vin','Vout')

subplot(1,2,2)
shrek = fftshift(fft(Vin));
plot(abs(shrek), 'r')

hold on
shrek2 = fftshift(fft(Vo));
plot(abs(shrek2), 'b')
%xlim([460 540])

title('Fourier frequency response for Step Input')
xlabel('Frequency (Hz)')
ylabel('Magnitude')
legend('Vin','Vout')

%Sinusodal input -------------------------------------
figure
subplot(2,2,1)
for i=1:time
    Vin(i) = sin(2*pi*f*i);
    
   A = C/0.001 + G;
   V_p = V;
   V = A\((C*V_p/0.001) + F*Vin(i));
   Vo(i) = V(7);
   
   if (i >= 2)
        plot([i i-1], [Vin(i) Vin(i-1)], 'r')
        hold on
        plot([i i-1], [Vo(i) Vo(i-1)], 'b')
   end
    
end
title('Sinusoid Input with period of 30ms')
xlabel('time (ms)')
ylabel('Voltage (V)')
legend('Vin','Vout')

subplot(2,2,2)
shrek = fftshift(fft(Vin));
plot(abs(shrek), 'r')

hold on
shrek2 = fftshift(fft(Vo));
plot(abs(shrek2), 'b')
%xlim([400 600])

title('Fourier frequency response for Sinusoid Input')
xlabel('Frequency (Hz)')
ylabel('Magnitude')
legend('Vin','Vout')

%Sinusodal input - 2nd frequency -------------------------
subplot(2,2,3)
f = 1/(4*0.03*time);
for i=1:time
    Vin(i) = sin(2*pi*f*i);
    
   A = C/0.001 + G;
   V_p = V;
   V = A\((C*V_p/0.001) + F*Vin(i));
   Vo(i) = V(7);
   
   if (i >= 2)
        plot([i i-1], [Vin(i) Vin(i-1)], 'r')
        hold on
        plot([i i-1], [Vo(i) Vo(i-1)], 'b')
   end
    
end
title('Sinusoid Input with period of 120ms')
xlabel('time (ms)')
ylabel('Voltage (V)')
legend('Vin','Vout')

subplot(2,2,4)
shrek = fftshift(fft(Vin));
plot(abs(shrek), 'r')

hold on
shrek2 = fftshift(fft(Vo));
plot(abs(shrek2), 'b')
%xlim([400 600])

title('Fourier frequency response for slower Sinusoid Input')
xlabel('Frequency (Hz)')
ylabel('Magnitude')
legend('Vin','Vout')

%Gaussian Input ------------------------------------
figure
subplot(1,2,1)
for i=1:time
    %gaussian parameters
    a = 1;
    b = 60;
    c = 30;
    Vin(i) = a*exp(-(i-b)^2/(2*c^2));

    A = C/0.001 + G;
   V_p = V;
   V = A\((C*V_p/0.001) + F*Vin(i));
   Vo(i) = V(7);

    if (i >= 2)
        plot([i i-1], [Vin(i) Vin(i-1)], 'r')
        hold on
        plot([i i-1], [Vo(i) Vo(i-1)], 'b')
    end
end
title('Gaussian Input with std. dev 30ms and delay 60ms')
xlabel('time (ms)')
ylabel('Voltage (V)')
legend('Vin','Vout')

subplot(1,2,2)
shrek = fftshift(fft(Vin));
plot(abs(shrek), 'r')

hold on
shrek2 = fftshift(fft(Vo));
plot(abs(shrek2), 'b')
%xlim([460 540])

title('Fourier frequency response for Gaussian Input')
xlabel('Frequency (Hz)')
ylabel('Magnitude')
legend('Vin','Vout')
