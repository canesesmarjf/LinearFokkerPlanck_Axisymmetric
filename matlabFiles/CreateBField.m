% Objective: create magnetic field profile

clear all
close all
clc

saveData = 1;
saveFig  = 1;

% Select coil geometry:
% =========================================================================
coilType = "c";        

switch coilType
    case {"a","b"}
        Lx = 6;         % Domain length [m]
        Lx_offset = 0;  % Offset [m]
        N = 501;        %  Number of elements on profile 
        rM       = 0.6;   % mirror coil radius in [m]
        IM       = 1.5e2; % Mirror coil current [kA]
        zM       = 1.5;   % Coil location along "z" [m]
        n  =     [+5  ,+5  ]*1e3;
        R  =     [+rM ,+rM ];
        z0 =     [-zM ,+zM ] + Lx_offset;
        I  =     [+IM ,+IM];
    case {"c"}
        Lx = 3;             % Domain length [m]
        Lx_offset = 0;      % Offset [m]
        N = 501;            %  Number of elements on profile 
        rM       = 0.22;    % mirror coil radius in [m]
        IM       = 15*15e3; % Mirror coil current [kA]
        zM       = 0.875;   % Coil location along "z" [m]
        n  =     [+25 ,+25 ];
        R  =     [+rM ,+rM ];
        z0 =     [-zM ,+zM ] + Lx_offset;
        I  =     [+IM ,+IM];        
end

% Calculate magnetic field:
% =========================================================================
% Vacuum magnetic field function:
% Function based on simple current loops:
mu0 = 4*pi*1e-7;
f = @(s) 0.5*mu0*sum((n.*I./R).*(1 + ((s-z0)./R).^2 ).^(-3/2));

% Axial domain:
z_B = linspace(-Lx/2,Lx/2,N)' + Lx_offset;

% Calculate vacuum magnetic field profile:
for ii = 1:numel(z_B)
    B(ii,1) = f(z_B(ii));
end

% Plot data:
% =========================================================================
figure('color','w'); 
plot(z_B,B,'LineWidth',2)
ylabel('B [T]','Interpreter','latex','FontSize',13)
xlabel('x [m]','Interpreter','latex','FontSize',13)
set(gcf,'color','w')
box on
grid on
ylim([0,1.2*max(B)])

% Save data to text file:
% =========================================================================
if saveData
    fileName = 'Bfield_c.txt';
    f = [z_B,B];
    save(fileName,'f','-ascii');
end

% Create figure of magnetic field:
% =========================================================================
if saveFig
    saveas(gcf,[fileName(1:end-4)],'tif');
end