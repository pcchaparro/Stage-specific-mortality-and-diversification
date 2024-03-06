%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function definitions - Calculate ecological equilibrium of an ancestral
% population colonizing an environment with 2 food resources
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Catalina Chaparro
%
%   original version: 10.08.2022,
%   last version: 30.03.2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate attack rates
AR = AMAX.*exp(-((Trait-THETA).^2)./(2*(TAU^2))); %attack rate on preexisting resources
dadeta = AMAX*(THETA-Trait)./TAU^2.*exp(-(THETA-Trait).^2./(2*TAU^2)); %da/dTrait

%Initialize vectors to save results
svPhi=size(Phivar);
minProd=zeros(svPhi(1,2),1);
Zequi=zeros(svPhi(1,2),1);
Nequi=zeros(svPhi(1,2),1);
mortS=zeros(svPhi(1,2),1);
SFood=zeros(svPhi(1,2),1);
dm=zeros(svPhi(1,2),1);
Ev=zeros(svPhi(1,2),1);
initialpoint=2;

for i=1:svPhi(1,2)
    phi   = Phivar(1,i);
    deltaS = stagespmort(i,1);

    %Z* from analytical results
    Zequi(i,1) = (deltaB*gamma - 8*deltaS - 4*deltaB + 2*deltaS*gamma - 2*gamma*maintenance + 4*deltaS.*phi + (- 3*deltaB^2*gamma^2 + 8*deltaB^2*gamma - 2*deltaB*deltaS*gamma^2.*phi - 4*deltaB*deltaS*gamma^2 + 16*deltaB*deltaS*gamma + 4*deltaB*gamma^2*maintenance - 4*deltaB*gamma*maintenance*omega + 5*deltaS^2*gamma^2.*phi.^2 - 12*deltaS^2*gamma^2.*phi + 4*deltaS^2*gamma^2 - 8*deltaS^2*gamma.*phi.^2 + 16*deltaS^2*gamma.*phi - 4*deltaS*gamma^2*maintenance.*phi + 8*deltaS*gamma^2*maintenance + 4*deltaS*gamma*maintenance*omega.*phi - 8*deltaS*gamma*maintenance*omega + 4*gamma^2*maintenance^2 - 8*gamma*maintenance^2*omega + 4*maintenance^2*omega^2)^(1/2) + 2*maintenance*omega - deltaS*gamma.*phi)/(2*(deltaB*gamma - 4*deltaS - 2*deltaB - 2*gamma*maintenance + 2*deltaS.*phi + 2*maintenance*omega + deltaS*gamma.*phi));

    Z = Zequi(i,1);

    %Define compound parameters
    dTm = deltaB+((phi.*Z)+(2-phi).*(1-Z)).*deltaS + (2-omega).*(1-Z).*maintenance;
    G   = gamma.*Z + (2-gamma).*(1-Z);
    E   = Ea.*(2-gamma).*(1-Z);
    B   = E.*FTmax/2.*G.*AR(1,1).*AR(1,2);
    dm(i,1) = dTm;
    Ev(i,1) = E;

    %Calculate food density in equilibrium

    FoodR(1,1) = (2.*AR(1,2).*rho.*FTmax/2.*dTm.*G)./((rho^2.*(4*B.^2+dTm.^2.*G.^2.*(AR(1,2)-AR(1,1)).^2)).^(1/2)+rho.*(2.*B+dTm.*G.*(AR(1,2)-AR(1,1))));
    FoodR(1,2) = (2.*AR(1,1).*rho.*FTmax/2.*dTm.*G)./((rho^2.*(4*B.^2+dTm.^2.*G.^2.*(AR(1,2)-AR(1,1)).^2)).^(1/2)+rho.*(2.*B-dTm.*G.*(AR(1,2)-AR(1,1))));
    SFood(i,1)= sum(FoodR);

    %Calculate population density

    aN = dTm.*G.^2.*AR(1,1).*AR(1,2);
    bN = rho*dTm.*G.*(AR(1,1)+AR(1,2)) - 2.*E.*G.*(rho*FTmax/2*AR(1,1).*AR(1,2));
    cN = rho^2*dTm-E.*(rho^2*FTmax/2*(AR(1,1)+AR(1,2)));

    Nequi(i,1) = (-bN+(bN.^2-4.*aN.*cN).^(1/2))./(2.*aN);

    %Calculate mortality due to stage-specific mortality

    mortS(i,1) = ((phi.*Z)+(2-phi).*(1-Z)).*deltaS;

    %Calculate minimum productivity

    minProd(i,1)=dTm*rho*(THETA(1,2)-THETA(1,1))^2/(4*E*TAU^2*AMAX)*exp((THETA(1,2)-THETA(1,1))^2/(8*TAU^2));
end

