function dTr = EvoDyn(t,Tr,n,m,Food,Juve,Adul,Traitres,Ea,AMAX,THETAF,TAU,mu,sigma,bmort,delta,phi,gamma,omega,maintenance)

    %Evolutionary dynamics
    % Here the rate of change of the feeding trait is calculated using the
    % canonical equation of adaptive dynamics for structured populations
    % (eq. 4).

Traitchange = zeros(m,1);
ingestFresjuv = zeros(m,1);
ingestFresadu = zeros(m,1);
ingestFmutjuv = zeros(m,1);
ingestFmutadu = zeros(m,1);
birthres = zeros(m,1);
maturationres = zeros(m,1);
birthmut = zeros(m,1);
maturationmut = zeros(m,1);
FitnessGradR0 = zeros(m,1);
ARres   = zeros(m,n);
ARmut   = zeros(m,n);
popsize = zeros(m,1);
Ts      = zeros(m,1);
R0mut   = zeros(m,1);

for j=1:m
    %Attack rates
    ARres(j,:) = AMAX.*exp(-((Traitres(j,1)-THETAF).^2)./(2*(TAU^2))); %attack rate of resident
    ARmut(j,:) = AMAX.*exp(-((Tr(j,1)-THETAF).^2)./(2*(TAU^2))); %attack rate of mutant
end

for j=1:m
    %feeding on basal resources
    for i=1:n %number of resources
        ingestFresjuv(j,1) = ingestFresjuv(j,1) + gamma*ARres(j,i)* Food(i,1); 
        ingestFresadu(j,1) = ingestFresadu(j,1) + (2-gamma)*ARres(j,i)* Food(i,1); 
        ingestFmutjuv(j,1) = ingestFmutjuv(j,1) + gamma*ARmut(j,i)* Food(i,1);
        ingestFmutadu(j,1) = ingestFmutadu(j,1) + (2-gamma)*ARmut(j,i)* Food(i,1); 
    end
    birthres(j,1)=max(Ea*ingestFresadu(j,1)-(2-omega)*maintenance,0);     %birth rate resident
    maturationres(j,1)=max(Ea*ingestFresjuv(j,1)-omega*maintenance,0);    %maturation rate resident
    
    birthmut(j,1)=max(Ea*ingestFmutadu(j,1)-(2-omega)*maintenance,0);     %birth rate mutant
    maturationmut(j,1)=max(Ea*ingestFmutjuv(j,1)-omega*maintenance,0);    %maturation rate mutant
    
    popsize(j,1)=Juve(j,1)+Adul(j,1); %population size of the resident
    Ts(j,1)=1/(maturationres(j,1)+delta*phi+bmort)*(1+maturationres(j,1)/(delta*(2-phi)+bmort));
    R0mut(j,1) = maturationmut(j,1)./(maturationmut(j,1)+delta*phi+bmort).*birthmut(j,1)./((2-phi)*delta+bmort); %R0 of mutant
    
    %Fitness gradient expression
    
    C1=phi*delta+bmort;
    C2=delta*(2-phi)+bmort;
    Cb=Ea*(2-gamma);
    Cm=Ea*gamma;

    dadeta = AMAX*(THETAF-Tr(j,1))./TAU^2.*exp(-(THETAF-Tr(j,1)).^2./(2*TAU^2));
    dfideta= sum(dadeta.*Food');
    fi= sum(ARmut(j,:).*Food');
    dR0dfi = (Cb*Cm^2*fi^2+2*fi*Cb*Cm*(C1-omega*maintenance)+Cb*omega^2*maintenance^2+C1*omega*maintenance*(Cm-Cb)-2*Cm*C1*maintenance)/(C2*(Cm*fi-omega*maintenance+C1)^2);

    FitnessGradR0(j,1) = dfideta.*dR0dfi;
    Traitchange(j,1) = popsize(j,1)./Ts(j,1) .*mu*sigma.* 1./R0mut(j,1) .* FitnessGradR0(j,1);
end

dTr = Traitchange;
end