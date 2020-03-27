clear all
% close all

homeAddress = cd;
addpath(homeAddress);
rootAddress = 'C:\Users\nfc\Documents\Programming\LinearFokkerPlanck\outputFiles';

cd(rootAddress);
folderName = 'Case13';
cd(folderName);

% Simulation conditions:
% =========================================================================
fileName = 'Metadata.out';
fileName = 'data.out';
metadata = GetMetadata(fileName,1);
dt  = metadata.DT;
Te0 = metadata.TE0;
Ti0 = metadata.TI0;
Ne0 = metadata.NE0;
try
    zmin = metadata.ZDUMP;
    zmax = metadata.ZTARGET;
catch
    zmin = metadata.ZMIN;
    zmax = metadata.ZMAX;
end

for ii = 1;
% Extract Binary data:
% =========================================================================
% Create file IDs:
% =========================================================================
fid{1}  = fopen('zp.out','r');
fid{2}  = fopen('kep.out','r');
fid{3}  = fopen('xip.out','r');
fid{4}  = fopen('tp.out','r');
fid{5}  = fopen('pcount1.out','r');
fid{6}  = fopen('pcount2.out','r');
fid{7}  = fopen('pcount3.out','r');
fid{8}  = fopen('pcount4.out','r');
fid{9}  = fopen('ecount1.out','r');
fid{10} = fopen('ecount2.out','r');
fid{11} = fopen('ecount3.out','r');
fid{12} = fopen('ecount4.out','r');

% Read binary data:
% =========================================================================
f1 = ExtractBinaryData(fid(1:4),8);
f2 = ExtractBinaryData(fid(5:12),8);

% Assign data to variables:
% =========================================================================
zp{ii}        = f1{1};
kep{ii}       = f1{2};
xip{ii}       = f1{3};
tp{ii}        = f1{4};
pcount1_0{ii} = f2{1};
pcount2_0{ii} = f2{2};
pcount3_0{ii} = f2{3};
pcount4_0{ii} = f2{4};
ecount1_0{ii} = f2{5};
ecount2_0{ii} = f2{6};
ecount3_0{ii} = f2{7};
ecount4_0{ii} = f2{8};
tc{ii} = 0:dt:((numel(pcount1_0{ii})-1)*dt);
clearvars f1 f2 fid
end
rmpath(homeAddress)

% Magnetic field data
% =========================================================================
BfieldAddress = 'C:\Users\nfc\Documents\Programming\LinearFokkerPlanck\BfieldData\';
[~,ii] = find(metadata.BFIELDFILE == 't');
BfieldFileName = metadata.BFIELDFILE(3:ii(end));
f = load([BfieldAddress,BfieldFileName]);
zb = f(:,1);
b  = f(:,2);
bmax = max(b);

cd(homeAddress);

%% Plot data:
% =========================================================================
if 0
    figure; 
   
    EmaxPlot = 45000;   
    EminPlot = 0;
    n_persist = 35;
    k = 1;
    for ii = 1:6:(size(zp{k},2)-n_persist)
        t_rng = ii:(ii+n_persist);
        plot(zp{k}(:,t_rng),kep{k}(:,t_rng)*1e-3,'k.','MarkerSize',4)
        hold on
        plot(zb,0.8*b*EmaxPlot*(1e-3)/bmax)
        hold off
        title(['t: ',num2str(tp{k}(ii)*1e6),' {\mu}s'])
        ylim([EminPlot,EmaxPlot*1e-3])
        ylabel('E$_\alpha$ [keV]','Interpreter','latex','FontSize',13)
        zlabel('z [m]','Interpreter','latex','FontSize',13)
        set(gca,'Yscale','lin')
        xlim([zmin,zmax])
        drawnow
        pause(0.01)
    end
end

%% Calculate the slowing down of ions

% Plasma conditions:
Tb = metadata.TE0;
Ea = metadata.KEP_INIT;
nb = metadata.NE0;
switch metadata.SPECIES_A
    case 1
        Ma = m_e;
    case 2
        Ma = m_p*metadata.AION;
end
Mb1 = m_e;
Mb2 = m_p*metadata.AION;
t0 = 0;

% RK4 solution
% -------------------------------------------------------------------------
                     %nb,Tb,Zb,Mb ,Za,Ma,y,yType,EqType
eqType = 3;
nu_E = @(E) nu_energy(nb,Tb,1 ,Mb1,1 ,Ma,E,'E'  ,eqType) + nu_energy(nb,Tb,1,Mb2,1,Ma,E,'E',eqType);
frk4 = @(E) -nu_E(E).*E;

% Energy relaxation time:
t_E = 1./nu_E(Ea);

tn = linspace(0,5*t_E,1e3) - t0;
dtn = mean(diff(tn));

E_RK4(1) = Ea;
for i = 1:(length(tn)-1)
    K1 = frk4(E_RK4(i) + 0.0   )*dtn;
    K2 = frk4(E_RK4(i) + 0.5*K1)*dtn;
    K3 = frk4(E_RK4(i) + 0.5*K2)*dtn;
    K4 = frk4(E_RK4(i) + 1.0*K3)*dtn;
    E_RK4(i+1) = E_RK4(i) + (K1 + 2*K2 + 2*K3 + K4)/6;
end

figure('color','w')
subplot(2,2,[1:2])
hold on
for ii = 1:metadata.NPARTS
    plot(tp{1}*1e3,kep{1}(ii,:),'k-','LineWidth',0.5)
end
plot(tn*1e3   ,E_RK4 ,'r-','LineWidth',2)
xlim([0,max(tp{1})]*1e3)
ylabel('E$_\alpha$ [keV]','Interpreter','latex','FontSize',13)
xlabel('t [ms]','Interpreter','latex','FontSize',13)
box on

subplot(2,2,[3]);
plot(tc{1}*1e3,ecount4_0{1},'k.')
ylabel('ecount4 [eV]','Interpreter','latex','FontSize',13)
xlabel('t [ms]','Interpreter','latex','FontSize',13)

subplot(2,2,[4]);
plot(tc{1}*1e3,pcount4_0{1},'r.')
ylabel('pcount4','Interpreter','latex','FontSize',13)
xlabel('t [ms]','Interpreter','latex','FontSize',13)

figure('color','w');
subplot(1,2,1)
hold on
plot(tc{1}*1e-3,cumsum(pcount1_0{1}),'k','LineWidth',2)
plot(tc{1}*1e-3,cumsum(pcount2_0{1}),'r','LineWidth',2)
box on
subplot(1,2,2)
hold on
plot(tc{1}*1e-3,cumsum(ecount1_0{1}),'k','LineWidth',2)
plot(tc{1}*1e-3,cumsum(ecount2_0{1}),'r','LineWidth',2)
box on

figure('color','w');
subplot(1,2,1)
hold on
plot(tc{1}*1e-3,cumsum(pcount3_0{1}),'k.','LineWidth',2)
box on
title('cumsum of pcount3')
subplot(1,2,2)
hold on
plot(tc{1}*1e-3,cumsum(ecount3_0{1}),'k.','LineWidth',2)
box on
title('cumsum of ecount3')

% What is the temperature with flow and it is it very different to the no
% flow case?
% We need to have an option to recycle particles into the source and/or to
% inject them like with a NBI
% 