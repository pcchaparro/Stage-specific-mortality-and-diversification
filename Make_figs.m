%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function definitions - Make figure 2-4 of manuscript "Differential
% stage-specific mortality as a mechanism for diversification"
%
% Other m-files required: SimulateEcoEvoDyn.m, EcoequivsPhi.m,
% findavgmort.m, FitnessLandscapePhi.m and files required by this rutines
% Subfunctions: none
% MAT-files required: none
%
% Author: Catalina Chaparro
%
%   original version: 13.10.2023
%   last version: 13.10.2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
%-------- PARAMETERS -----------

% Fixed Parameters

n=2;                      % Number of resources
rho=.01;                  % Resource growth rate
FTmax=1;                  % Total carrying capacity
AMAX=0.1;                 % Maximum attack rate
TAU=.4;                   % Standard deviation of the Gaussian that determines how attack rate varies with trait
THETA=[0 1];              % Optimal trait values to feed on the resources
Ea=.6;                    % Efficiency to transform ingested mass into biomass
deltaB=.002;              % Background mortality
deltaS=.01;               % Stage-dependent mortality
maintenance=0;            % Metabolic maintenance
omega=1;                  % Ontogenetic asymmetry in metabolic maintenance. If omega>1, then juveniles have higher metabolic maintenance cost.
phi=1;                    % Ontogenetic asymmetry in mortality. If phi>1, then juveniles have higher mortality rate.
gamma=1;                  % Ontogenetic asymmetry in feeding. If gamma>1, then juveniles have a higher attack rate

mu = 1E-3;                % mutation rate
sigma = .01;              % variance of mutational steps

%Parameters to be varied
Phivar = 0:.01:1.99;

%Calculate attack rates
Trait = 0.5;
AR = AMAX.*exp(-((Trait-THETA).^2)./(2*(TAU^2))); %attack rate on preexisting resources
dadeta = AMAX*(THETA-Trait)./TAU^2.*exp(-(THETA-Trait).^2./(2*TAU^2)); %da/dTrait

%Calculate the deltaS values for Phivar that make the average mortality constant 
svPhi=size(Phivar);
stagespmort = zeros(svPhi(1,2),1);
objdeltaS = deltaS; %target stage-specific mortality to be kept constant while varying phi
initialpoint = deltaS;
for i=1:svPhi(1,2)
    phi = Phivar(1,i);
    fun=@(deltaS) findavgmort(deltaB,deltaS,phi,gamma,omega,maintenance,objdeltaS);
    stagespmort(i,1) = fzero(fun,initialpoint);
    initialpoint = stagespmort(i,1);
end

%%%%% FIGURE 2 %%%%%%
figure
suptitle('Figure 2')
% -----Simulate the evolutionary dynamics of phi = 1 and phi = 1.5

%FIGURE 2A
%parameter
phi=1;                    % Differential mortality between life stages. If phi>1, then juveniles have higher mortality rate.
deltaS=stagespmort(Phivar==phi,1);
%Simulate
SimulateEcoEvoDyn 
%Plot
subplot(2,2,1)
tfig=0;
for i=1:cont
    plot(tfig+(1:tspeciation{i}),TTrai{i},'k')
    hold on
    tfig = tfig + tspeciation{i};
end
ylim([0 1])
title('phi = 1')
xlabel('Time')
ylabel('Trait')

% FIGURE 2B

%parameter
phi=1.6;                    % Differential mortality between life stages. If phi>1, then juveniles have higher mortality rate.
deltaS=stagespmort(Phivar==phi,1);
%Simulate
SimulateEcoEvoDyn 
%Plot
subplot(2,2,2)
tfig=0;
for i=1:cont
    plot(tfig+(1:tspeciation{i}),TTrai{i},'k')
    hold on
    tfig = tfig + tspeciation{i};
end
title('phi = 1.6')
ylim([0 1])
xlabel('Time')
ylabel('Trait')

% FIGURE 2C

% Parameters to be varied
Traitvar = 0:.01:1;       % Trait value

% Compute
FitnessLandscapePhi
%plot
numcolor=2;
cm=[52 77 126;222 222 166]./255;
subplot(2,2,[3 4])
X = Phivar;
Y = Traitvar;
C = contourf(Phivar,Traitvar,(CurvFit'),'LevelList',0);
Cbreak = max(max(C)) + 1;
x = [C(1,2:Cbreak), X(end), X(1), X(1), X(end), C(1,Cbreak+2:end)];
y = [C(2,2:Cbreak), Y(1), Y(1), Y(end), Y(end), C(2,Cbreak+2:end)];
contourf(Phivar,Traitvar,(FitnessGrad'),[-10, 0, 10]);
colormap(cm)
caxis([-3 3])
G = colorbar;
G.Label.String = 'Fitness gradient';
hold on
curv=fill(x,y,'m','FaceAlpha',0.3,'FaceColor',[0.85 0.325 .098]);
hold off
xlabel('Differential mortality')
ylabel('Feeding niche trait')
[~,h_legend] = legend(curv(1),'Curvature of fitness landscape is positive','Location','southoutside');
PatchInLegend = findobj(h_legend, 'type', 'patch');
set(PatchInLegend(1), 'FaceAlpha', 0.3);


%%%%% FIGURE 3 %%%%%%
figure
suptitle('Figure 3')

%Find the threshold of minimum productivity as a function of phi
%Parameter
Trait = 0.5;
%Compute
EcoequivsPhi
%FIGURE 3A
subplot(1,3,1)
plot(Phivar,1-Zequi,'k')
ylim([0 1])
xlabel('Differential mortality')
ylabel('Fraction of adults')
subplot(1,3,2)
plot(Phivar,Nequi,'k')
xlabel('Differential mortality')
ylabel('Population density')
subplot(1,3,3)
plot(Phivar,minProd,'k')
xlabel('Differential mortality')
ylabel({'Minimum productivity required'; 'for diversification'})


%%%%% FIGURE 4 %%%%%%
fig4=figure(3);
fig4.Renderer='Painters';
suptitle('Figure 4')
FTmax=10;                  % Total carrying capacity

%Calculate the deltaS values for Phivar that make the average mortality
%constant for varying differential metabolic cost between life stages

omegavar = 0:.01:2;
svomega = size(omegavar);
stagespmort2D = zeros(svPhi(1,2),svomega(1,2));
maintenance = 0.01;
initialpoint = objdeltaS;
for j=1:svomega(1,2)
    omega = omegavar(1,j);
    for i=1:svPhi(1,2)
        phi = Phivar(1,i);
        fun=@(deltaS) findavgmort(deltaB,deltaS,phi,gamma,omega,maintenance,objdeltaS);
        stagespmort2D(i,j) = fzero(fun,initialpoint);
        if i<svPhi(1,2)
            initialpoint = stagespmort2D(i,j);
        else
            initialpoint = stagespmort2D(1,j);
        end
    end
end

%Find the threshold of minimum productivity as a function of phi
%Parameter
Trait = 0.5;
%create vectors to save output
Zequi2D = zeros(svPhi(1,2),svomega(1,2));
Nequi2D = zeros(svPhi(1,2),svomega(1,2));
minProd2D = zeros(svPhi(1,2),svomega(1,2));
%Compute
for j=1:svomega(1,2)
    omega = omegavar(1,j);
    stagespmort=stagespmort2D(:,j);
    EcoequivsPhi
    Zequi2D(:,j)=Zequi;
    Nequi2D(:,j)=Nequi;
    minProd2D(:,j)=minProd;
end

%FIGURE 4A
sp1=subplot(2,3,1);
numcolor=10;
limcolor=[222 222 166;52 77 126]./255;
cmr=linspace(limcolor(1,1),limcolor(2,1),numcolor);
cmg=linspace(limcolor(1,2),limcolor(2,2),numcolor);
cmb=linspace(limcolor(1,3),limcolor(2,3),numcolor);
cm1=[cmr' cmg' cmb'];
contourf(Phivar,omegavar,(1-Zequi2D)')
colormap(sp1,flipud(cm1));
caxis([0 1])
colorbar
xlabel('Differential mortality')
ylabel('Differential metabolic cost')
title('Fraction of adults')

%FIGURE 4B
sp2=subplot(2,3,2);
numcolor=6;
limcolor=[222 222 166;52 77 126]./255;
cmr=linspace(limcolor(1,1),limcolor(2,1),numcolor);
cmg=linspace(limcolor(1,2),limcolor(2,2),numcolor);
cmb=linspace(limcolor(1,3),limcolor(2,3),numcolor);
cm1=[cmr' cmg' cmb'];
contourf(Phivar,omegavar,(Nequi2D)')
colormap(sp2,flipud(cm1));
caxis([1 2.5])
colorbar
xlabel('Differential mortality')
ylabel('Differential metabolic cost')
title('Population density')

%FIGURE 4C
sp3=subplot(2,3,3);
numcolor=6;
limcolor=[222 222 166;52 77 126]./255;
cmr=linspace(limcolor(1,1),limcolor(2,1),numcolor);
cmg=linspace(limcolor(1,2),limcolor(2,2),numcolor);
cmb=linspace(limcolor(1,3),limcolor(2,3),numcolor);
cm1=[cmr' cmg' cmb'];
contourf(Phivar,omegavar,(minProd2D)',numcolor-1)
colormap(sp3,flipud(cm1));
caxis([.01 .022])
colorbar
xlabel('Differential mortality')
ylabel('Differential metabolic cost')
title({'Minimum productivity required'; 'for diversification'})

%Calculate the deltaS values for Phivar that make the average mortality
%constant for varying differential feeding rate between life stages

maintenance = 0;
omega = 1;
gammavar = 0.1:.01:1.9;
svgamma = size(gammavar);
stagespmort2D = zeros(svPhi(1,2),svgamma(1,2));
initialpoint = objdeltaS*10;
for j=1:svgamma(1,2)
    gamma = gammavar(1,j);
    for i=1:svPhi(1,2)
        phi = Phivar(1,i);
        fun=@(deltaS) findavgmort(deltaB,deltaS,phi,gamma,omega,maintenance,objdeltaS);
        stagespmort2D(i,j) = fzero(fun,initialpoint);
        if i<svPhi(1,2)
            initialpoint = real(stagespmort2D(i,j));
        else
            initialpoint = real(stagespmort2D(1,j));
        end
    end
end

%Find the threshold of minimum productivity as a function of phi
%Parameter
Trait = 0.5;
%create vectors to save output
Zequi2D = zeros(svPhi(1,2),svgamma(1,2));
Nequi2D = zeros(svPhi(1,2),svgamma(1,2));
minProd2D = zeros(svPhi(1,2),svgamma(1,2));
%Compute
for j=1:svgamma(1,2)
    gamma = gammavar(1,j);
    stagespmort=stagespmort2D(:,j);
    EcoequivsPhi
    Zequi2D(:,j)=Zequi;
    Nequi2D(:,j)=Nequi;
    minProd2D(:,j)=minProd;
end

%FIGURE 4D
sp4=subplot(2,3,4);
numcolor=10;
limcolor=[222 222 166;52 77 126]./255;
cmr=linspace(limcolor(1,1),limcolor(2,1),numcolor);
cmg=linspace(limcolor(1,2),limcolor(2,2),numcolor);
cmb=linspace(limcolor(1,3),limcolor(2,3),numcolor);
cm1=[cmr' cmg' cmb'];
contourf(Phivar,gammavar,(1-Zequi2D)')
colormap(sp4,flipud(cm1));
caxis([0 1])
colorbar
xlabel('Differential mortality')
ylabel('Differential foraging capacity')
title('Fraction of adults')

%FIGURE 4E
sp5=subplot(2,3,5);
numcolor=6;
limcolor=[222 222 166;52 77 126]./255;
cmr=linspace(limcolor(1,1),limcolor(2,1),numcolor);
cmg=linspace(limcolor(1,2),limcolor(2,2),numcolor);
cmb=linspace(limcolor(1,3),limcolor(2,3),numcolor);
cm1=[cmr' cmg' cmb'];
contourf(Phivar,gammavar,(Nequi2D)')
colormap(sp5,flipud(cm1));
caxis([1.5 4.5])
colorbar
xlabel('Differential mortality')
ylabel('Differential foraging capacity')
title('Population density')

%FIGURE 4F
sp6=subplot(2,3,6);
numcolor=6;
limcolor=[222 222 166;52 77 126]./255;
cmr=linspace(limcolor(1,1),limcolor(2,1),numcolor);
cmg=linspace(limcolor(1,2),limcolor(2,2),numcolor);
cmb=linspace(limcolor(1,3),limcolor(2,3),numcolor);
cm1=[cmr' cmg' cmb'];
contourf(Phivar,gammavar,log10(minProd2D)',numcolor-1)
colormap(sp6,flipud(cm1));
caxis([-2.2 -1])
colorbar
xlabel('Differential mortality')
ylabel('Differential foraging capacity')
title({'log10(Minimum productivity required'; 'for diversification)'})
