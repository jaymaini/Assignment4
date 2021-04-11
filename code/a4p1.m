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

%DC Case

s = 0; %no capacitive/inductive effects
Vin = linspace(-10, 10, 100);
V3 = zeros(size(Vin));
Vo = zeros(size(Vin));
for i=1:length(Vin)
   V = (G+s.*C)\(F.*Vin(i)); 
   V3(i) = V(4)*G3;
   Vo(i) = V(7);
end


figure(1)
%subplot(2,3,1)
plot(Vin,Vo)
hold on
plot(Vin,V3)
title('DC Case: Vin vs Vo & V3')
xlabel('Voltage')
ylabel('Gain (Vo/Vi)')
legend('Vo','V3')

f = linspace(0.01, 1000, 10000);
w = 2*pi*f;
Vo = zeros(size(f));
for i = 1:length(f)
    s = j*w(i);
    V = (G+s*C)\F;  % Solving
    Vo(i) = V(7);
end
%subplot(2,3,2)
figure(2)
plot(w,abs(Vo));
title('AC Sweep');
xlabel('\omega (rad/s)');
ylabel('Voltage');
ylim([0.1 0.3])
xlim([0 4000])

%subplot(2,3,3)
figure(3)
semilogx(f, 20*log(abs(Vo)),'LineWidth',3);
grid;
title('Frequency Response', 'FontSize',12);
xlabel('Frequency  (Hz)','FontSize',12);
ylabel('V_{out}/V_{1}  (dB)','FontSize',12);

%subplot(2,3,4)
figure(4)
subplot(2,1,1)
%normal distribution of Capacitance
C1_dist = C1 + 0.05 * randn(1000,1);
hist(C1_dist, 20)
title('C_in distribution')

%subplot(2,3,5)
subplot(2,1,2)
Vo = zeros(size(C1_dist));
s = 1i*pi;
%Solved distribution of gain 
for i = 1:length(C1_dist)
    C1 = C1_dist(i);
    %recalculate C matrix
    C = [0 0 0 0 0 0 0;
    -C1 C1 0 0 0 0 0;
    0 0 -L1 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0];
    
    V = (G+s.*C)\(F);
    Vo(i) = abs(V(7));
end
hist(Vo, 20)
title('Distribution of Gain')


