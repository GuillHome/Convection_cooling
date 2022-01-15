clear;
clc;

%DONNEES
sigma=5.67*10^(-12); %W/cm2.K4 
ray=0.93; 
cp= 460;	% capacité calorifique de l'acier (J/kg/°C)
Tamb=20; % Température ambiante (°C)
lambda=30*10^(-2);  % conductivité de l'acier (W/cm.K)
rho=7850*10^(-6); % masse volumique de l'acier (kg/cm3)
K=40;   % nombre de lignes
C=10;   % nombre de colonnes
N=160000;   % nombre de calculs sur le temps
Rext=14.6;  % rayon exterieur de la couronne (cm)
Rint=10;	% rayon interieur de la couronne (cm)
tetude=2*3600;  % temps de chauffage (s)
dt=tetude/(N-1);    % pas de temps  (s)
dr=(Rext-Rint)/(K-1);   % pas sur le rayon (cm)
dteta=15/(C-1);     % pas sur l'angle (degré)
r1=12.5-(3.75/2);   % bord inferieur du tube (cm)
r2=12.5+(3.75/2);   % bord superieur du tube (cm)
k1Ini=round(((Rext-r1)/dr)+1);  %indice pour le bord inferieur du tube (le + grand)
k2Ini=round(((Rext-r2)/dr)+1);  %indice pour le bord superieur du tube (le + petit)
Rc=1.875;   % rayon du tube (cm)
phio=9000*10^(-4);  % flux de chauffage (kW/m2)
mu=1.877*10^(-5);   % viscosité de l'air (Pa.s)
omega=500;  % vitesse de rotation de la couronne 
lambdaair=0.0242;   % conductivité de l'air  (W/m.K)
relax=0.001;  % coefficient de relaxation pour le calcul de h
epsilon=0.05;   % précision souhaitée sur h
pmax=60;   % nombre de boucles max à effectuer par le programme pour le calcul de h

%Calcul du he sur la surface extèrieure
he=(((0.073*lambdaair)/(2*Rext*10^(-2)))*((omega*(2*Rext*10^(-2))^2)/(2*mu))^0.7)*10^(-4);	

%CREATION DES MATRICES DE TRAVAIL
T=zeros(K,C,2);
deltaT=zeros(C);
%Tmes=file;    % matrice des températures mesurées par la caméra infrarouge
h=50*10^(-4);
deltah=0;

%CONDITIONS INITIALES DE T
for c=1:C
    for k=1:K
        T(k,c,1)=Tamb;
    end;
end;

%GEOMETRIE:CONSTRUCTION DU CERCLE
k1=k1Ini;
k2=k2Ini;
x=Rc+dr;

while k2-k1<=0
    x=x-dr;
    tetaC=asin(x/Rc);
    y=Rc*cos(tetaC);
    
    r=Rext-(k2-1)*dr;
    teta=y/r;
    c=round((teta*180/pi/dteta)+1);
    
    for i=1:c
       T(k2,i,1)=0;
    end

    r=Rext-(k1-1)*dr;
    teta=y/r;
    c=round((teta*180/pi/dteta)+1);
    for i=1:c
       T(k1,i,1)=0;
    end
    
    k1=k1-1;
    k2=k2+1;
end


%CALCUL DE LA TEMPERATURE ET DE h
r=zeros(1,K);
%suivi=zeros(1,N);
r(1)=Rext;

Tini=zeros(K,C,1);
for c=1:C
    for k=1:K
        Tini(k,c,1)=T(k,c,1);
    end;
end;

p=1;
while p<pmax  % boucle pour le calcul de h
    for c=1:C   % remise aux conditions initiales
        for k=1:K
            T(k,c,1)=Tini(k,c,1);
        end;
    end;
    for n=1:N-1           
        for k=2:K-1     % boucle sur les lignes   
            r(k)=Rext-(k-1)*dr;
            for cinv=2:C-1  %boucle sur les colonnes
                c=-cinv+C-1+2;  %inversion du sens de parcours pour faciliter la mise en oeuvre du programme
                if T(k,c,1)>0
            
                    if (T(k+1,c,1)>0)&&(T(k-1,c,1)>0)&&((T(k,c+1,1)>0)&&(T(k,c-1,1)>0))                                  
                        T(k,c,2)=T(k,c,1)+((dt*lambda)/(rho*cp))*((T(k+1,c,1)-T(k,c,1))/(dr*r(k))+((T(k+1,c,1)-2*T(k,c,1)+T(k-1,c,1))/dr^2)+((T(k,c+1,1)-2*T(k,c,1)+T(k,c-1,1))/(r(k)*dteta)^2));
                        T(K,c,2)=T(K-1,c,2);
                        T(1,c,2)=(1/(1+(he*dr/lambda)))*((T(2,c,2)+273)+(dr/lambda)*(phio+he*(Tamb+273)-sigma*ray*(((T(1,c,1)+273)^4)-(Tamb+273)^4)))-273; %condition de flux sur face ext
                    else
                        if T(k,c-1,1)==0
                            T(k,c,2)=(1/(h*(pi/2)+(2*lambda/dteta)))*((2*lambda/dteta)*T(k,c+1,2)+Tamb*h*(pi/2)); %condition de flux au niveau du tube
                        else
                            T(k,c,2)=T(k,c+1,2);
                        end;
                    end;                 
                end;
        
            end;
            T(k,C,2)=T(k,C-1,2);  % flux nul sur face latérale
        end;
    
        for k=1:k2Ini
            T(k,1,2)=T(k,2,2);  % flux nul au-dessus du tube 
        end;
    
        for k=k1Ini:K
            T(k,1,2)=T(k,2,2);  % flux nul en-dessous du tube
        end;
        
        %prend en charge les quatre bouts 
        T(K,C,2)=T(K-1,C,2);
        T(1,1,2)=T(1,2,2);
        T(1,C,2)=T(1,C-1,2);
        T(K,1,2)=T(K-1,1,2);
        
    %petite boucle permettant d'écraser les anciennes valeurs et éviter les
    %problèmes de mémoire vive
    for c=1:C
        for k=1:K
            T(k,c,1)=T(k,c,2);
        end;
    end;

    end;  %fin de la boucle sur le temps

    
    %suivi(n)=T(2,2,2); % permet le tracé de la température en fonction
    %du temps pour un point

    deltaT=T(1,:,2)'-suivi(:)

    if max(abs(deltaT))<epsilon %respect de la précision
        p=pmax-1;
    else  %calcul du nouveau h

    deltah=max(abs(deltaT))*relax
    h=h+sign(mean(T(1,:,2)'-suivi(:)))*deltah

    end;

    p=p+1;
end; %fin de la boucle pour le calcul de h

