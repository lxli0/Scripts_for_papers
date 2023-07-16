function [dt0,year2dt0,t0,Water_Load]= generate_waterload(amplitude,span,delay)% mm year month
    dt0=24*3600*365.25/12/10; %per 3 day
    year2dt0=24*3600*365.25/dt0;
    t0=0:dt0:span*year2dt0*dt0;
    FBCD=zeros(1,length(t0));
    for i=1:length(t0)
        FBCD(i)=sin(2*pi*(t0(i)/year2dt0/dt0-delay/12-0.25));
    end
    Water_Load=FBCD*1e3*9.8*amplitude/1000;
end
