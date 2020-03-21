% Create Bessel function of order 0 and 1

close all
clear all

s = linspace(0,10,501);

J0 = besselj(0,s);
J1 = besselj(1,s);

figure; 
hold on
plot(s,J0,'k','LineWidth',2)
plot(s,J1,'r','LineWidth',2)

f = [s',J0',J1'];
save('besselj01.txt','f','-ascii')


s = linspace(0,40,501);

J0 = besselj(0,s);
J1 = besselj(1,s);

figure; 
hold on
plot(s,J0,'k','LineWidth',2)
plot(s,J1,'r','LineWidth',2)

f = [s',J0',J1'];
save('besselj01_0_to_40.txt','f','-ascii')

f = [s',J0'];
save('besselj0_0_to_40.txt','f','-ascii')

f = [s',J1'];
save('besselj1_0_to_40.txt','f','-ascii')