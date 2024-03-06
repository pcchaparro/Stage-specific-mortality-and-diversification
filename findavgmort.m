function [avgstgmort] = findavgmort(deltaB,deltaS,phi,gamma,omega,maintenance,objdeltaS)

    %Z* from analytical results
    Z = (deltaB*gamma - 8*deltaS - 4*deltaB + 2*deltaS*gamma - 2*gamma*maintenance + 4*deltaS.*phi + (- 3*deltaB^2*gamma^2 + 8*deltaB^2*gamma - 2*deltaB*deltaS*gamma^2.*phi - 4*deltaB*deltaS*gamma^2 + 16*deltaB*deltaS*gamma + 4*deltaB*gamma^2*maintenance - 4*deltaB*gamma*maintenance*omega + 5*deltaS^2*gamma^2.*phi.^2 - 12*deltaS^2*gamma^2.*phi + 4*deltaS^2*gamma^2 - 8*deltaS^2*gamma.*phi.^2 + 16*deltaS^2*gamma.*phi - 4*deltaS*gamma^2*maintenance.*phi + 8*deltaS*gamma^2*maintenance + 4*deltaS*gamma*maintenance*omega.*phi - 8*deltaS*gamma*maintenance*omega + 4*gamma^2*maintenance^2 - 8*gamma*maintenance^2*omega + 4*maintenance^2*omega^2)^(1/2) + 2*maintenance*omega - deltaS*gamma.*phi)/(2*(deltaB*gamma - 4*deltaS - 2*deltaB - 2*gamma*maintenance + 2*deltaS.*phi + 2*maintenance*omega + deltaS*gamma.*phi));

    %Calculate de average mortality

    avgstgmort = ((phi*Z)+(2-phi)*(1-Z))*deltaS-objdeltaS;
end