%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function definitions - Fitness landscape of an ancestral
% population colonizing an environment with 2 food resources
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

% Initialize vectors to save information
svPhi = size(Phivar);
svtrait = size(Traitvar);
FitnessGrad = zeros(svPhi(1,2),svtrait(1,2)); %Fitness gradient
CurvFit = zeros(svPhi(1,2),svtrait(1,2)); %Curvature of fitness landscape

for i=1:svPhi(1,2)
    for j=1:svtrait(1,2)

        Trait = Traitvar(1,j);
        phi   = Phivar(1,i);
        deltaS = stagespmort(i,1);
        AR    = AMAX.*exp(-((Trait-THETA).^2)./(2*(TAU^2))); %attacked rate on each resource

        %Z* from analytical results
        Z = (deltaB*gamma - 8*deltaS - 4*deltaB + 2*deltaS*gamma - 2*gamma*maintenance + 4*deltaS.*phi + (- 3*deltaB^2*gamma^2 + 8*deltaB^2*gamma - 2*deltaB*deltaS*gamma^2.*phi - 4*deltaB*deltaS*gamma^2 + 16*deltaB*deltaS*gamma + 4*deltaB*gamma^2*maintenance - 4*deltaB*gamma*maintenance*omega + 5*deltaS^2*gamma^2.*phi.^2 - 12*deltaS^2*gamma^2.*phi + 4*deltaS^2*gamma^2 - 8*deltaS^2*gamma.*phi.^2 + 16*deltaS^2*gamma.*phi - 4*deltaS*gamma^2*maintenance.*phi + 8*deltaS*gamma^2*maintenance + 4*deltaS*gamma*maintenance*omega.*phi - 8*deltaS*gamma*maintenance*omega + 4*gamma^2*maintenance^2 - 8*gamma*maintenance^2*omega + 4*maintenance^2*omega^2)^(1/2) + 2*maintenance*omega - deltaS*gamma.*phi)/(2*(deltaB*gamma - 4*deltaS - 2*deltaB - 2*gamma*maintenance + 2*deltaS.*phi + 2*maintenance*omega + deltaS*gamma.*phi));

        %Define compound parameters
        dTm = deltaB+((phi.*Z)+(2-phi).*(1-Z)).*deltaS + (2-omega).*(1-Z).*maintenance;
        G   = gamma.*Z + (2-gamma).*(1-Z);
        E   = Ea.*(2-gamma).*(1-Z);
        B   = E.*FTmax/2.*G.*AR(1,1).*AR(1,2);

        %Calculate food density in equilibrium

        FoodR(1,1) = (2.*AR(1,2).*rho.*FTmax/2.*dTm.*G)./((rho^2.*(4*B.^2+dTm.^2.*G.^2.*(AR(1,2)-AR(1,1)).^2)).^(1/2)+rho.*(2.*B+dTm.*G.*(AR(1,2)-AR(1,1))));
        FoodR(1,2) = (2.*AR(1,1).*rho.*FTmax/2.*dTm.*G)./((rho^2.*(4*B.^2+dTm.^2.*G.^2.*(AR(1,2)-AR(1,1)).^2)).^(1/2)+rho.*(2.*B-dTm.*G.*(AR(1,2)-AR(1,1))));

        %Calculate fitness gradient

        C1=phi*deltaS+deltaB;
        C2=deltaS*(2-phi)+deltaB;

        dadeta  = AMAX*(THETA-Trait)./TAU^2.*exp(-(THETA-Trait).^2./(2*TAU^2));       
        dfideta= sum(dadeta.*FoodR);
        fi= sum(AR.*FoodR);
        dR0dfi = (Ea^3*fi^2+2*Ea^2*C1*fi)/(C2*(Ea*fi+C1)^2);

        FitnessGrad(i,j) = dfideta.*dR0dfi; %eq. S1.6
        
        %Calculate the curvature of the fitness landscape
        d2adeta2  = -(AMAX*exp(-(THETA-Trait).^2./(2*TAU^2)).*(- Trait^2 + 2*Trait.*THETA + TAU^2 - THETA.^2))/TAU^4;
        d2fideta2 = sum(d2adeta2.*FoodR);
        d2R0dfi2  = 2*Ea^2*C1^2/(C2*(Ea*fi+C1)^3);   
        
        CurvFit(i,j) = d2R0dfi2*(dfideta)^2 + dR0dfi*d2fideta2; %eq. S1.7
    end
end