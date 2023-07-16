clc,clear
dh=100/1e3;
depth=dh/2:dh:100-dh/2;
depth=depth*1e3;
eta=10.^(-5-depth/1e3);
eta(eta<1e-9)=1e-9;
Z1=Calculatiion_Z1(depth,eta);
diff=[0.1 1 10];

[dt0,year2dt0,t0,F]=generate_waterload(350,1,0);%generate TWS
tws2gws=[1/2 2/3];
gwsdelaytws=[0 1 2];
num2str=['1to2_0delay';'1to2_1delay';'1to2_2delay';'2to3_0delay';'2to3_1delay';'2to3_2delay'];
for m=1:length(tws2gws)
    for n=1:length(gwsdelaytws)
        [dt0_pp,year2dt0_pp,t0_pp,F_pp]=generate_waterload(tws2gws(m)*350,10,gwsdelaytws(n));%generate GWS
        %% elastic load
        for j=1:length(depth)
            s11_0=-3/7;s22_0=-3/7;s33_0=-1;s12_0=0;s13_0=0;s23_0=0;% compression negative, niu=0.3
            for i=1:length(t0)
                s11(i)=F(i)*s11_0;
                s22(i)=F(i)*s22_0;
                s33(i)=F(i)*s33_0;
                s12(i)=F(i)*s12_0;
                s13(i)=F(i)*s13_0;
                s23(i)=F(i)*s23_0;
                delta_stress_3d_Boussinesq=[
                   s11(i) s12(i) s13(i)
                   s12(i) s22(i) s23(i)
                   s13(i) s23(i) s33(i)];
                [dsigma1,dsigma2,dsigma3]=principal_stress_change(delta_stress_3d_Boussinesq);
                [dVp1,dVp2,dVp3,dVs1,dVs2,dVs3]=stress2velocity(-dsigma1,-dsigma2,-dsigma3,Z1(j));
                dVs2Vs_loading(i,j)=(sqrt((2*(2e3+dVs3)^2+(2e3+dVs1)^2)/3)-2e3)/2e3;                  
            end
        end

        %% pore pressure and effective stress
        for k=1:3
            for j=1:length(depth)
                [pp,p_undrained,p_diffusion]=oned_pore_pressure(diff(k),depth(j),t0_pp,F_pp,length(t0));
                for i=1:length(t0)
                    delta_stress_pp=[
                       pp(i) 0 0
                       0 pp(i) 0
                       0 0 pp(i)];
                    [dsigma1,dsigma2,dsigma3]=principal_stress_change(delta_stress_pp);
                    [dVp1,dVp2,dVp3,dVs1,dVs2,dVs3]=stress2velocity(-dsigma1,-dsigma2,-dsigma3,Z1(j));
                    dVs2Vs_pp(i,j,k)=(sqrt((2*(2e3+dVs3)^2+(2e3+dVs1)^2)/3)-2e3)/2e3;  
                end
                s11_0=-3/7;s22_0=-3/7;s33_0=-1;s12_0=0;s13_0=0;s23_0=0;% compression negative, niu=0.3
                for i=1:length(t0)
                        s11(i)=F(i)*s11_0;
                        s22(i)=F(i)*s22_0;
                        s33(i)=F(i)*s33_0;
                        s12(i)=F(i)*s12_0;
                        s13(i)=F(i)*s13_0;
                        s23(i)=F(i)*s23_0; 
                        delta_stress=[
                           s11(i)+pp(i) s12(i) s13(i)
                           s12(i) s22(i)+pp(i) s23(i)
                           s13(i) s23(i) s33(i)+pp(i)];
                        [dsigma1,dsigma2,dsigma3]=principal_stress_change(delta_stress);
                        [dVp1,dVp2,dVp3,dVs1,dVs2,dVs3]=stress2velocity(-dsigma1,-dsigma2,-dsigma3,Z1(j));
                        dVs2Vs_eff(i,j,k)=(sqrt((2*(2e3+dVs3)^2+(2e3+dVs1)^2)/3)-2e3)/2e3;
                end
            end
        end
        fileName = [num2str(length(gwsdelaytws)*(m-1)+n,1:11),'.mat'];
        save(fileName, 'dVs2Vs_loading','dVs2Vs_pp','dVs2Vs_eff');
    end
end