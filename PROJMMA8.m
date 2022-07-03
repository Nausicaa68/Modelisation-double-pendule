%%authors
%DILLY, DUMAS, GRELAUD, HASSOUN, ROSSIGNOL

clc;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%% INITIALISATION DES CONSTANTES ET CONDITIONS INITIALES%%%%%%%%%%%%%%%%%%%%%%%%
T=input('fin intervalle ?');
N=input('nbr de pts à approximer');
%y10=input('init 1 ?');
y10=pi/2;
%y20=input('init 2 ?');
y20=pi/2;
%y30=input('init 3 ?');
y30=0;
%y40=input('init 4 ?');
y40=0;
%M1=input('M1 ?');
M1=1;
%M2=input('M2 ?');
M2=60;
g=1;
%L1=input('L1 ?');
L1=1;
%L2=input('L2 ?');
L2=1;
%%%%%%%%%%%%%%%%%%%%%%%% INITIALISATION DES CONSTANTES ET CONDITIONS INITIALES%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%% DEFINITION DE LA FONCTION F%%%%%%%%%%%%%%%%%%%%%%%%
%f1=@(y1,y2,y3,y4) y3;
%f2=@(y1,y2,y3,y4) y4;
f3=@(y1,y2,y3,y4) (1/((L1*L2)*(M1 + (M2*sin(y1-y2)*sin(y1-y2)))) )*( (-(M1+M2)*g*L2*sin(y1)) + (M2*g*L2*sin(y2)*cos(y1-y2)) - (M2*L1*L2*sin(y1-y2)*cos(y1-y2)*y3*y3) - (M2*L2*L2*sin(y1-y2)*y4*y4) );
f4=@(y1,y2,y3,y4) (1/((L1*L2)*(M1 + (M2*sin(y1-y2)*sin(y1-y2)))) )*( ((M1+M2)*g*L1*sin(y1)*cos(y1-y2)) - ((M1+M2)*g*L1*sin(y2)) + ((M1+M2)*L1*L1*sin(y1-y2)*y3*y3) + ((M2*L1*L2*sin(y1-y2))*(cos(y1-y2)*y4*y4))) ;
%K = (1/((L1*L2)*(M1 + (M2*sin(y1-y2)*sin(y1-y2)))) );
%%%%%%%%%%%%%%%%%%%%%%%% DEFINITION DE LA FONCTION F%%%%%%%%%%%%%%%%%%%%%%%%


h=T/(N-1);  %creation du pas
t=0:h:T;%initialisation pour la boucle for


%%%%%%%%%%%%%%%%%%%%%%%% INITIALISATION EULER%%%%%%%%%%%%%%%%%%%%%%%%
%les INDICES commencent tjrs à 1 dans matlab
%Initialisation, conditions initiales methode euler
y1(1)=y10;
y2(1)=y20;
y3(1)=y30;
y4(1)=y40;
M2x(1)= sin(y10)*L1 + sin(y20)*L2;
M2y(1)= L1+L2 - (cos(y10)*L1 + cos(y20)*L2);
%%%%%%%%%%%%%%%%%%%%%%%% INITIALISATION EULER%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% INITIALISATION POINT MILIEU %%%%%%%%%%%%%%%%%%%%%%%%
y1milieu(1)=y10;
y2milieu(1)=y20;
y3milieu(1)=y30;
y4milieu(1)=y40;

M2xmilieu(1)= sin(y10)*L1 + sin(y20)*L2;
M2ymilieu(1)= L1+L2 - (cos(y10)*L1 + cos(y20)*L2);
%%%%%%%%%%%%%%%%%%%%%%%% INITIALISATION POINT MILIEU %%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% INITIALISATION HEUN %%%%%%%%%%%%%%%%%%%%%%%%
%Initialisation, conditions initiales
y1heun(1)=y10;
y2heun(1)=y20;
y3heun(1)=y30;
y4heun(1)=y40;

M2xheun(1)= sin(y10)*L1 + sin(y20)*L2;
M2yheun(1)= L1+L2 - (cos(y10)*L1 + cos(y20)*L2);
%%%%%%%%%%%%%%%%%%%%%%%% INITIALISATION HEUN %%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% INITIALISATION RK4 %%%%%%%%%%%%%%%%%%%%%%%%
%Initialisation, conditions initiales
y1RK4(1)=y10;
y2RK4(1)=y20;
y3RK4(1)=y30;
y4RK4(1)=y40;

M2xRK4(1)= sin(y10)*L1 + sin(y20)*L2;
M2yRK4(1)= L1+L2 - (cos(y10)*L1 + cos(y20)*L2);
%%%%%%%%%%%%%%%%%%%%%%%% INITIALISATION RK4 %%%%%%%%%%%%%%%%%%%%%%%%



%boucle calculant les 4 schemas
for i=2:N
    
   %%%%%%%%%%%%%%%%%%%%%%%% EULER %%%%%%%%%%%%%%%%%%%%%%%%
   %formule d'Euler : y(n+1) = y(n) + h*F(y(n))
   
   y1(i) = y1(i-1) + h*y3(i-1);
   y2(i) = y2(i-1) + h*y4(i-1);
   y3(i) = y3(i-1) + h*f3(y1(i-1),y2(i-1),y3(i-1),y4(i-1));
   y4(i) = y4(i-1) + h*f4(y1(i-1),y2(i-1),y3(i-1),y4(i-1));
   
   %calcul des pos x et y du point M2 en fc du temps
   M2x(i)= sin(y1(i))*L1 + sin(y2(i))*L2;
   M2y(i)= L1+L2 - (cos(y1(i))*L1 + cos(y2(i))*L2);
   %%%%%%%%%%%%%%%%%%%%%%%% EULER %%%%%%%%%%%%%%%%%%%%%%%%
   
   
   %%%%%%%%%%%%%%%%%%%%%%%% POINT MILIEU %%%%%%%%%%%%%%%%%%%%%%%%
   %formule du point milieu : y(n+1) = y(n) + h*F( y(n) + h/2*F(y(n)) )
   
   %calculs intermédiaries
   % y(n)+ (h/2)*f(y(n))
   z1milieu = y1milieu(i-1) + (h/2)*y3milieu(i-1);
   z2milieu = y2milieu(i-1) + (h/2)*y4milieu(i-1);
   z3milieu = y3milieu(i-1) + (h/2)*f3(y1milieu(i-1),y2milieu(i-1),y3milieu(i-1),y4milieu(i-1));
   z4milieu = y4milieu(i-1) + (h/2)*f4(y1milieu(i-1),y2milieu(i-1),y3milieu(i-1),y4milieu(i-1));

   %calcul des y(n+1)
   %y1(i) = y1(i-1) + h*f1(z1,z2,z3,z4);   
   %y2(i) = y2(i-1) + h*f2(z1,z2,z3,z4);   
   %y3(i) = y3(i-1) + h*f3(z1,z2,z3,z4);   
   %y4(i) = y4(i-1) + h*f4(z1,z2,z3,z4);
   
   y1milieu(i) = y1milieu(i-1) + h*z3milieu;   
   y2milieu(i) = y2milieu(i-1) + h*z4milieu;   
   y3milieu(i) = y3milieu(i-1) + h*f3(z1milieu,z2milieu,z3milieu,z4milieu);   
   y4milieu(i) = y4milieu(i-1) + h*f4(z1milieu,z2milieu,z3milieu,z4milieu); 
   
   %calcul des pos x et y du point M2 en fc du temps
   M2xmilieu(i)= sin(y1milieu(i))*L1 + sin(y2milieu(i))*L2;
   M2ymilieu(i)= L1+L2 - (cos(y1milieu(i))*L1 + cos(y2milieu(i))*L2);
   %%%%%%%%%%%%%%%%%%%%%%%% POINT MILIEU %%%%%%%%%%%%%%%%%%%%%%%%
   
   %%%%%%%%%%%%%%%%%%%%%%%% HEUN %%%%%%%%%%%%%%%%%%%%%%%%
   %formule de Heun : y(n+1) = y(n) + h/2*(F(y(n)) + F( y(n)+ h*F(y(n)) ) )
   
   %calculs intermédiaries
   % y(n) + h*f(y(n))
   z1heun = y1heun(i-1) + (h)*y3heun(i-1);
   z2heun = y2heun(i-1) + (h)*y4heun(i-1);
   z3heun = y3heun(i-1) + (h)*f3(y1heun(i-1),y2heun(i-1),y3heun(i-1),y4heun(i-1));
   z4heun = y4heun(i-1) + (h)*f4(y1heun(i-1),y2heun(i-1),y3heun(i-1),y4heun(i-1));

    
   %calcul des y(n+1)
   %y1(i) = y1(i-1) + (h/2) * (f1(y1,y2,y3,y4) + f1(z1,z2,z3,z4));   
   %y2(i) = y2(i-1) + (h/2) * (f2(y1,y2,y3,y4) + f2(z1,z2,z3,z4));   
   %y3(i) = y3(i-1) + (h/2) * (f3(y1,y2,y3,y4) + f3(z1,z2,z3,z4));   
   %y4(i) = y4(i-1) + (h/2) * (f4(y1,y2,y3,y4) + f4(z1,z2,z3,z4));
   
   y1heun(i) = y1heun(i-1) + (h/2)*( y3heun(i-1) + z3heun);  
   y2heun(i) = y2heun(i-1) + (h/2)*( y4heun(i-1) + z4heun);   
   y3heun(i) = y3heun(i-1) + (h/2)*( f3(y1heun(i-1),y2heun(i-1),y3heun(i-1),y4heun(i-1)) + f3(z1heun,z2heun,z3heun,z4heun) );   
   y4heun(i) = y4heun(i-1) + (h/2)*( f4(y1heun(i-1),y2heun(i-1),y3heun(i-1),y4heun(i-1)) + f4(z1heun,z2heun,z3heun,z4heun) ); 
   
   
   %calcul des pos x et y du point M2 en fc du temps
   M2xheun(i)= sin(y1heun(i))*L1 + sin(y2heun(i))*L2;
   M2yheun(i)= L1+L2 - (cos(y1heun(i))*L1 + cos(y2heun(i))*L2);
   %%%%%%%%%%%%%%%%%%%%%%%% HEUN %%%%%%%%%%%%%%%%%%%%%%%%
   
   %%%%%%%%%%%%%%%%%%%%%%%% RK4 %%%%%%%%%%%%%%%%%%%%%%%%
   %calcul k1
   k11 = y3RK4(i-1);
   k12 = y4RK4(i-1);
   k13 = f3(y1RK4(i-1),y2RK4(i-1),y3RK4(i-1),y4RK4(i-1));
   k14 = f4(y1RK4(i-1),y2RK4(i-1),y3RK4(i-1),y4RK4(i-1));
    
   %calcul k2
   z21 = y1RK4(i-1) + ((h/2)*(k11));
   z22 = y2RK4(i-1) + ((h/2)*(k12));
   z23 = y3RK4(i-1) + ((h/2)*(k13));
   z24 = y4RK4(i-1) + ((h/2)*(k14));
   
   k21 = z23;
   k22 = z24;
   k23 = f3(z21,z22,z23,z24);
   k24 = f4(z21,z22,z23,z24);
   
   %calcul k3
   z31 = y1RK4(i-1) + ((h/2)*(k21));
   z32 = y2RK4(i-1) + ((h/2)*(k22));
   z33 = y3RK4(i-1) + ((h/2)*(k23));
   z34 = y4RK4(i-1) + ((h/2)*(k24));
   
   k31 = z33;
   k32 = z34;
   k33 = f3(z31,z32,z33,z34);
   k34 = f4(z31,z32,z33,z34);
   
   %calcul k4
   z41 = y1RK4(i-1) + ((h)*(k31));
   z42 = y2RK4(i-1) + ((h)*(k32));
   z43 = y3RK4(i-1) + ((h)*(k33));
   z44 = y4RK4(i-1) + ((h)*(k34));
   
   k41 = z43;
   k42 = z44;
   k43 = f3(z41,z42,z43,z44);
   k44 = f4(z41,z42,z43,z44);
    
   %calcul des y(n+1)
   y1RK4(i) = y1RK4(i-1) + ((h/6)*(k11+(2*k21)+(2*k31)+k41));
   y2RK4(i) = y2RK4(i-1) + ((h/6)*(k12+(2*k22)+(2*k32)+k42));
   y3RK4(i) = y3RK4(i-1) + ((h/6)*(k13+(2*k23)+(2*k33)+k43));
   y4RK4(i) = y4RK4(i-1) + ((h/6)*(k14+(2*k24)+(2*k34)+k44));
   
   
   %calcul des pos x et y du point M2 en fc du temps
   M2xRK4(i)= sin(y1RK4(i))*L1 + sin(y2RK4(i))*L2;
   M2yRK4(i)= L1+L2 - (cos(y1RK4(i))*L1 + cos(y2RK4(i))*L2);
   %%%%%%%%%%%%%%%%%%%%%%%% RK4 %%%%%%%%%%%%%%%%%%%%%%%%
   
end


% --- Tracés ---
% on peut améliorer l'approximation en augmentant le nbr de pts 
% c'est surtout necessaire dans le cas de la methode d'Euler

%plot des 4schemas pos angulaire y2 en fc du temps
plot(t,y2,t,y2milieu,t,y2heun,t,y2RK4);legend('Euler','Point milieu','Heun','RK4');

%plot pour la position angulaire y1
plot(t,y1RK4);legend('RK4');

%coords (x,y) de M2 (le trapéziste) en fc du temps
%les 4 schémas
%plot(M2x,M2y,M2xmilieu,M2ymilieu,M2xheun,M2yheun,M2xRK4,M2yRK4);legend('Euler','Point milieu','Heun','RK4');

%affichage des shcémas mais sans Euler
%plot(M2xmilieu,M2ymilieu,M2xheun,M2yheun,M2xRK4,M2yRK4);legend('Point milieu','Heun','RK4');

%plot pour le sous-espace des phases entre  y2 et y'2=y4 
%plot(y2RK4,y4RK4);legend('Espace des phases y2');

%plot pour les 2 sous-espaces des phases (y2,y'2) et (y1,y'1)
%plot(y2RK4,y4RK4,y1RK4,y3RK4);legend('Espace des phases y2','Espace des phases y1');
