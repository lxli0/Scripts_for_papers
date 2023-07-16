function [pbc_out,p_undrained_out,p_diffusion_out]= oned_pore_pressure(D,z,t0,Q,out_length)
    %% parameter
    B=0.5;
    v=0.31;%泊松比
    a=B*(1+v)/(3*(1-v));
    dp=[Q(1),diff(Q)];
    %% calculation
    t=0.5+t0;
    pbc=zeros(1,length(t));
    p_undrained=pbc;
    p_diffusion=pbc;
    for j=1:length(t)
        t_c=t(j);
        for k=1:length(find(t0<t_c))
            t_load=t0(k);
            pbc(j)=pbc(j)+((1-a)*erfc(z/sqrt(4*D*(t_c-t_load)))+a)*dp(k);
            p_undrained(j)=p_undrained(j)+a*erf(z/sqrt(4*D*(t_c-t_load)))*dp(k);
            p_diffusion(j)=p_diffusion(j)+erfc(z/sqrt(4*D*(t_c-t_load)))*dp(k);
        end
    end
   
    %% Output
    t_out=t(end-out_length+1:end);
    pbc_out=pbc(end-out_length+1:end);
    p_undrained_out=p_undrained(end-out_length+1:end);
    p_diffusion_out=p_diffusion(end-out_length+1:end);
    
    %{
    plot(t/24/3600/365.25,pbc,'-','color',1/255*[46 89 167],'linewidth',2.5)
    hold on;
    plot(t/24/3600/365.25,p_undrained,'-','color',1/255*[242 200 103],'linewidth',2.5)
    hold on;
    plot(t/24/3600/365.25,p_diffusion,'-','color',1/255*[237 109 61],'linewidth',2.5)
    hold on;
    plot(t0/24/3600/365.25,Q,'-','color',1/255*[21 29 41],'linewidth',2.5)
    legend('Total','Undrained','Diffusion','Water Load','NumColumns',1);
    xlabel('Year')
    ylabel('Pore Pressure')
    grid on;
    box on;
    %xlim([0 10])
    set(gca,'fontsize',16)
    %}
end