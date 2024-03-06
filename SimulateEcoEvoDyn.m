%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function definitions - Simulate evolutionary dynamics of an ancestral
% population colonizing an environment with n food resources
%
% Other m-files required: EcoEqui.m, EvoDyn.m
% Subfunctions: none
% MAT-files required: none
%
% Author: Catalina Chaparro
%
%   original version: 10.08.2022
%   last version: 30.03.2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------
%ANCESTRAL POPULATION

%Initial conditions
initialTrai = 0.4;                  %Initial feeding trait
m=1;                                %One ancestral population

%eco-evo dynamics

tmax1 = 10000;
Food = zeros(n,tmax1);
Juve = zeros(m,tmax1);
Adul = zeros(m,tmax1);
Trai = zeros(m,tmax1);

Trai(:,1) = initialTrai;

AR   = zeros(m,n);
traitchange = 1;
time=1;

while (time<tmax1-1 && traitchange>1E-8)
    for j=1:m
    %Attack rates
        AR(j,:) = AMAX.*exp(-((Trai(j,time)-THETA).^2)./(2*(TAU^2))); %attack rate of morph j on preexisting resources
    end
    
    %find the ecological equilibrium
    %Z* from analytical results
    Z = (deltaB*gamma - 8*deltaS - 4*deltaB + 2*deltaS*gamma - 2*gamma*maintenance + 4*deltaS.*phi + (- 3*deltaB^2*gamma^2 + 8*deltaB^2*gamma - 2*deltaB*deltaS*gamma^2.*phi - 4*deltaB*deltaS*gamma^2 + 16*deltaB*deltaS*gamma + 4*deltaB*gamma^2*maintenance - 4*deltaB*gamma*maintenance*omega + 5*deltaS^2*gamma^2.*phi.^2 - 12*deltaS^2*gamma^2.*phi + 4*deltaS^2*gamma^2 - 8*deltaS^2*gamma.*phi.^2 + 16*deltaS^2*gamma.*phi - 4*deltaS*gamma^2*maintenance.*phi + 8*deltaS*gamma^2*maintenance + 4*deltaS*gamma*maintenance*omega.*phi - 8*deltaS*gamma*maintenance*omega + 4*gamma^2*maintenance^2 - 8*gamma*maintenance^2*omega + 4*maintenance^2*omega^2)^(1/2) + 2*maintenance*omega - deltaS*gamma.*phi)/(2*(deltaB*gamma - 4*deltaS - 2*deltaB - 2*gamma*maintenance + 2*deltaS.*phi + 2*maintenance*omega + deltaS*gamma.*phi));

    %Define compound parameters
    dTm = deltaB+((phi.*Z)+(2-phi).*(1-Z)).*deltaS + (2-omega).*(1-Z).*maintenance;
    G   = gamma.*Z + (2-gamma).*(1-Z);
    E   = Ea.*(2-gamma).*(1-Z);
    B   = E.*FTmax/2.*G.*AR(1,1).*AR(1,2);

    %Calculate food density in equilibrium
    RFood(1,1) = (2.*AR(1,2).*rho.*FTmax/2.*dTm.*G)./((rho^2.*(4*B.^2+dTm.^2.*G.^2.*(AR(1,2)-AR(1,1)).^2)).^(1/2)+rho.*(2.*B+dTm.*G.*(AR(1,2)-AR(1,1))));
    RFood(1,2) = (2.*AR(1,1).*rho.*FTmax/2.*dTm.*G)./((rho^2.*(4*B.^2+dTm.^2.*G.^2.*(AR(1,2)-AR(1,1)).^2)).^(1/2)+rho.*(2.*B-dTm.*G.*(AR(1,2)-AR(1,1))));

    %Calculate population density
    aN = dTm.*G.^2.*AR(1,1).*AR(1,2);
    bN = rho*dTm.*G.*(AR(1,1)+AR(1,2)) - 2.*E.*G.*(rho*FTmax/2*AR(1,1).*AR(1,2));
    cN = rho^2*dTm-E.*(rho^2*FTmax/2*(AR(1,1)+AR(1,2)));
    Nequi = (-bN+(bN.^2-4.*aN.*cN).^(1/2))./(2.*aN);
    
    %Save ecological equilibrium in terms of juveniles and adults
    Food(:,time+1) = RFood;
    Juve(:,time+1) = Z*Nequi;
    Adul(:,time+1) = (1-Z)*Nequi;

    Traitres=Trai(:,time);
    fod   = @(t,Tr) EvoDyn(t,Tr,n,m,Food(:,time+1),Juve(:,time+1),Adul(:,time+1),Traitres,Ea,AMAX,THETA,TAU,mu,sigma,deltaB,deltaS,phi,gamma,omega,maintenance);
    [tevo,xevo] = ode23(fod,[0 1E6],Trai(:,time));
    Trai(:,time+1) = xevo(end,:);
    
    traitvar = zeros(m,1);
    for j=1:m
        traitvar(j,1) = abs(Trai(j,time+1)-Trai(j,time));
    end
    traitchange = max(traitvar);
    time = time+1;
end

TFood{1} = Food(:,1:time);
TJuve{1} = Juve(:,1:time);
TAdul{1} = Adul(:,1:time);
TTrai{1} = Trai(:,1:time);
tspeciation{1} = time;

%Calculate the 2nd derivative to determine if selection is disruptive

disruptive = zeros(1,m);

% defining some constants
C1=phi*deltaS+deltaB;
C2=deltaS*(2-phi)+deltaB;
    
for j=1:m
    ARres = AMAX.*exp(-((TTrai{1}(j,end)-THETA).^2)./(2*(TAU^2))); %attack rate of resident
    ingestFres = 0;
    for i=1:n %number of resources
        ingestFres = ingestFres + ARres* TFood{1}(i,end); 
    end
    birthres=max(Ea*ingestFres(j,1),0);     %birth rate resident
    maturationres=max(Ea*ingestFres(j,1),0);    %maturation rate resident
    
    dadeta = AMAX*(THETA-TTrai{1}(j,end))./TAU^2.*exp(-(THETA-TTrai{1}(j,end)).^2./(2*TAU^2));
    d2adeta2 = -(AMAX*exp(-(THETA-TTrai{1}(j,end)).^2./(2*TAU^2)).*(- TTrai{1}(j,end)^2 + 2*TTrai{1}(j,end).*THETA + TAU^2 - THETA.^2))/TAU^4;
    fi     = sum(ARres.*TFood{1}(:,end)'); %fi is landa in the appendix
    dfideta= sum(dadeta.*TFood{1}(:,end)');
    d2fideta2= sum(d2adeta2.*TFood{1}(:,end)');
    
    dR0dfi = (Ea^3*fi^2+2*Ea^2*C1*fi)/(C2*(Ea*fi+C1)^2);
    d2R0dfi2= 2*Ea^2*C1^2/(C2*(Ea*fi+C1)^3);   
    d2R0deta2 = d2R0dfi2*(dfideta)^2 + dR0dfi*d2fideta2; %second derivative of R0 with respect to eta mutant

    if d2R0deta2>1E-12
        disruptive(1,j)=1;
    else
        disruptive(1,j)=0;
    end
end

%--------------------------------------------------------------------
%SPECIATION
%If selection is disruptive split the population(s) with disruptive
%selection

timemaxrad=1E8;
timerad=0;
cont = 1;
numsp = 1;
numspnopred = 1;

while (sum(disruptive)>0 && m<n+10) && timerad<timemaxrad
    
    timerad = timerad + time;
    
    initialFood = TFood{cont}(:,end);
    initialJuve = [];
    initialAdul = [];
    initialTrai = [];
    
    %Diversification

    for j=1:m
        if disruptive(1,j)==1
            initialJuve = [initialJuve; TJuve{cont}(j,end)/2; TJuve{cont}(j,end)/2];
            initialAdul = [initialAdul; TAdul{cont}(j,end)/2; TAdul{cont}(j,end)/2];
            initialTrai = [initialTrai; TTrai{cont}(j,end)*1-1E-3; TTrai{cont}(j,end)*1+1E-3];
        else
            initialJuve = [initialJuve; TJuve{cont}(j,end)];
            initialAdul = [initialAdul; TAdul{cont}(j,end)];
            initialTrai = [initialTrai; TTrai{cont}(j,end)];
        end
    end
    
    m=m+sum(disruptive);
    
    %eco-evo dynamics

    tmax1 = 10000;
    Food = zeros(n,tmax1);
    Juve = zeros(m,tmax1);
    Adul = zeros(m,tmax1);
    Trai = zeros(m,tmax1);

    Food(:,1) = initialFood;
    Juve(:,1) = initialJuve;
    Adul(:,1) = initialAdul;
    Trai(:,1) = initialTrai;

    AR   = zeros(m,n);
    traitchange = 1;
    time=1;

    tic

    while (time<tmax1-1 && traitchange>1E-8)
        for j=1:m
        %Attack rates
            AR(j,:) = AMAX.*exp(-((Trai(j,time)-THETA).^2)./(2*(TAU^2))); %attack rate of morph j on preexisting resources
        end

        %find the ecological equilibrium using numerical simulations
        %because we do not have analytical expressions for a lineage with 
        %more than 1 population
        
        tmaxe = 100;
        dif = 1;
        rep = 0;
        x0    = [Food(:,time)' Juve(:,time)' Adul(:,time)'];
        while (dif>1E-8 && rep<20)
            fod   = @(t,x) EcoEqui(t,x,AR,n,m,rho,FTmax,Ea,deltaB,deltaS,phi);
            [t1,x1] = ode45(fod,[0 tmaxe],x0);
            dif = abs(x1(end-1,1)-x1(end,1));
            rep = rep + 1;
            x0  = x1(end,:);
        end

        Food(:,time+1)=x1(end,1:n);
        Juve(:,time+1)=x1(end,n+1:n+m);
        Adul(:,time+1)=x1(end,n+m+1:n+2*m);

        Traitres=Trai(:,time);
        fod   = @(t,Tr) EvoDyn(t,Tr,n,m,Food(:,time+1),Juve(:,time+1),Adul(:,time+1),Traitres,Ea,AMAX,THETA,TAU,mu,sigma,deltaB,deltaS,phi,gamma,omega,maintenance);
        [tevo,xevo] = ode23(fod,[0 1E6],Trai(:,time));
        Trai(:,time+1) = xevo(end,:);

        traitvar = zeros(m,1);
        for j=1:m
            traitvar(j,1) = abs(Trai(j,time+1)-Trai(j,time));
        end
        traitchange = max(traitvar);
        time = time+1;
    end

    TFood{cont+1} = Food(:,1:time);
    TJuve{cont+1} = Juve(:,1:time);
    TAdul{cont+1} = Adul(:,1:time);
    TTrai{cont+1} = Trai(:,1:time);
    tspeciation{cont+1} = time;

    %Calculate the 2nd derivative to determine if selection is disruptive

    disruptive = zeros(1,m);

    for j=1:m
        ARres = AMAX.*exp(-((TTrai{cont+1}(j,end)-THETA).^2)./(2*(TAU^2))); %attack rate of resident
        ingestFres = 0;
        for i=1:n %number of resources
            ingestFres = ingestFres + ARres(1,i)* TFood{cont+1}(i,end); 
        end
        birthres=max(Ea*ingestFres,0);     %birth rate resident
        maturationres=max(Ea*ingestFres,0);    %maturation rate resident

        dadeta = AMAX*(THETA-TTrai{cont+1}(j,end))./TAU^2.*exp(-(THETA-TTrai{cont+1}(j,end)).^2./(2*TAU^2));
        d2adeta2 = -(AMAX*exp(-(THETA-TTrai{cont+1}(j,end)).^2./(2*TAU^2)).*(- TTrai{cont+1}(j,end)^2 + 2*TTrai{cont+1}(j,end).*THETA + TAU^2 - THETA.^2))/TAU^4;
        fi     = sum(ARres.*TFood{cont+1}(:,end)'); %fi is landa in the appendix
        dfideta= sum(dadeta.*TFood{cont+1}(:,end)');
        d2fideta2= sum(d2adeta2.*TFood{cont+1}(:,end)');
        
        dR0dfi = (Ea^3*fi^2+2*Ea^2*C1*fi)/(C2*(Ea*fi+C1)^2);
        d2R0dfi2= 2*Ea^2*C1^2/(C2*(Ea*fi+C1)^3);
        d2R0deta2 = d2R0dfi2*(dfideta)^2 + dR0dfi*d2fideta2; %second derivative of R0 with respect to eta mutant
        
        if d2R0deta2>0
            disruptive(1,j)=1;
        else
            disruptive(1,j)=0;
        end
    end
    cont = cont+1;
end