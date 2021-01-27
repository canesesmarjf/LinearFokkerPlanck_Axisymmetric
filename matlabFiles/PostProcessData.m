% This script provides an example how to post-process the simulation data:

GetData

%% Animate data:
if 1
    figure; 
   
    EmaxPlot = 1000;   
    EminPlot = 0;
    n_persist = 0;
    k = 1;
    for ii = 1:6:(size(zp,2)-n_persist)
        t_rng = ii:(ii+n_persist);
        plot(zp(:,t_rng),kep(:,t_rng)*1e-3,'k.','MarkerSize',4)
        hold on
        plot(zb,0.8*b*EmaxPlot*(1e-3)/bmax)
        hold off
        title(['t: ',num2str(tp(ii)*1e6),' {\mu}s'])
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

if 0
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
        plot(tp*1e3,kep(ii,:),'k-','LineWidth',0.5)
    end
    plot(tn*1e3   ,E_RK4 ,'r-','LineWidth',2)
    xlim([0,max(tp)]*1e3)
    ylabel('E$_\alpha$ [keV]','Interpreter','latex','FontSize',13)
    xlabel('t [ms]','Interpreter','latex','FontSize',13)
    box on

    subplot(2,2,[3]);
    plot(tc*1e3,ecount4_0{1},'k.')
    ylabel('ecount4 [eV]','Interpreter','latex','FontSize',13)
    xlabel('t [ms]','Interpreter','latex','FontSize',13)

    subplot(2,2,[4]);
    plot(tc*1e3,pcount4_0{1},'r.')
    ylabel('pcount4','Interpreter','latex','FontSize',13)
    xlabel('t [ms]','Interpreter','latex','FontSize',13)

    figure('color','w');
    subplot(1,2,1)
    hold on
    plot(tc*1e-3,cumsum(pcount1_0),'k','LineWidth',2)
    plot(tc*1e-3,cumsum(pcount2_0),'r','LineWidth',2)
    box on
    subplot(1,2,2)
    hold on
    plot(tc*1e-3,cumsum(ecount1_0),'k','LineWidth',2)
    plot(tc*1e-3,cumsum(ecount2_0),'r','LineWidth',2)
    box on

    figure('color','w');
    subplot(1,2,1)
    hold on
    plot(tc*1e-3,cumsum(pcount3_0),'k.','LineWidth',2)
    box on
    title('cumsum of pcount3')
    subplot(1,2,2)
    hold on
    plot(tc*1e-3,cumsum(ecount3_0),'k.','LineWidth',2)
    box on
    title('cumsum of ecount3')
end