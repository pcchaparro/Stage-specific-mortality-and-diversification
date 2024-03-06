function F = EcoEqui(t,x,AR,n,m,rho,FTmax,Ea,bmort,delta,phi)
    %Ecological dynamics
    % AR is the matrix of attack rates of all ecomorphs on preexisting
    % resources. Rows correspond to each consumer ecomorph and
    % columns to the attacked resource.
    % x is the vector of variables to solve, where first n entries are
    % food abundance, n+1 to n+m entries correspond to abundance of stage
    % 1, n+m+1 to n+2m entries correspond to abundance of stage 2,...
    
F = zeros(n+2*m,1);

Food=x(1:n,1);
stages=zeros(m,2);


for i=1:2
    stages(:,i) = x(n+1+m*(i-1):n+m*i,1);
end

juvstage=1;
adustage=2;

ingestF = zeros(m,1);
birth = zeros(m,1);
maturation = zeros(m,1);
Juvchange = zeros(m,1);
Aduchange = zeros(m,1);
for j=1:m %number of morphs
    %feeding on basal resources
    for i=1:n %number of resources
        ingestF(j,1) = ingestF(j,1) + AR(j,i)* Food(i,1); 
    end
    birth(j,1)=max(Ea*ingestF(j,1),0);     %birth rate
    maturation(j,1)=max(Ea*ingestF(j,1),0);    %maturation rate
    
    Juvchange(j,1) = birth(j,1)*stages(j,adustage)-maturation(j,1)*stages(j,juvstage)-(phi*delta+bmort)*stages(j,juvstage);
    Aduchange(j,1) = maturation(j,1)*stages(j,juvstage)-((2-phi)*delta+bmort)*stages(j,adustage);
end

Fmax=FTmax/n;
    
grazing= zeros(n,1);
Foodchange = zeros(n,1);
for i=1:n %number of food resources
    for j=1:m %number of morphs
        grazing(i,1) = grazing(i,1) + AR(j,i)*Food(i,1)*(stages(j,adustage)+stages(j,juvstage));
    end
    Foodchange(i,1) = rho*(Fmax-Food(i,1)) - grazing(i,1);
end

F(1:n,1)=Foodchange;
F(n+1:n+m,1)=Juvchange;
F(n+m+1:n+2*m,1)=Aduchange;
end