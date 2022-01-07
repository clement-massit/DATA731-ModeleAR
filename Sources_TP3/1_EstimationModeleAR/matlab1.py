
clear;
close all;
clc;

n=1536; 
a1(1)=1;						
a1(2)=-0.1344;
a1(3)=0.9025;
% coefficients du second processus AR
a2(1)=1;						
a2(2)=-1.6674;
a2(3)=0.9025;
% coefficients du troisieme processus AR
a3(1)=1;						
a3(2)=1.7820;
a3(3)=0.8100;

t=-2:1:n-1;
y=0*t;
for k=3:n/3
    z=randn;
    y(k)=-a1(2)*y(k-1)-a1(3)*y(k-2)+z;
end
%
%  Generation du 2nd bloc de la serie (largeur n/3)
for k=n/3+1:2*n/3
    z=randn;
    y(k+2)=-a2(2)*y(k+1)-a2(3)*y(k)+z;
end
%
%  Generation du 3e bloc de la serie (largeur n/3)
for k=2*n/3+1:n
    z=randn;
    y(k+2)=-a3(2)*y(k+1)-a3(3)*y(k)+z;
end
echo on
% Trace- de la serie
t = t(3:end);
y = y(3:end);
figure
plot(t,y);
title('Serie = juxtapososition de 3 sous-series stationnaires');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Calcul et trace-s des spectres des trois sous-series a partir de freqz
%  Puisque l on connait les coefficients, le calcul est fait directement
%  a- partir des coefficients (ne depend donc pas du nombre d-echantillons) 
Np = 256; % nombre de points du spectre
[hh1,w]=freqz(1,a1,Np); % calcul du spectre de la sous-s�rie 1
mag11=abs(hh1);
[hh2,w]=freqz(1,a2,Np); % calcul du spectre de la sous-s�rie 2
mag12=abs(hh2);
[hh3,w]=freqz(1,a3,Np); % calcul du spectre de la sous-s�rie 3
mag13=abs(hh3);
f=w; %/(2*pi)*1000;		% �chantillonnage de l axe des fr�quences
figure
semilogy(f,mag11,'-g',f,mag12,':b',f,mag13,'--r'); % trac� du spectre
title('Spectres des sous-s�ries');
legend('Sous-s�rie 1', 'Sous-s�rie 2', 'Sous-s�rie 3');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pause
%% On choisit de decrire y par un modele AR d-ordre 3, puis d'ordre 4.
%   Estimation des coefficients des modeles AR d-ordres 3 et 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcul global sur n donn�es 
%
% Calcul des coefficients des mod�les AR d 'ordre 3 et 4
% On cr�e les versions translat�es de y (style autocorr�lation)
y1=[y 0 0 0 0]';
y2=[0 y 0 0 0]';
y3=[0 0 y 0 0]';
y4=[0 0 0 y 0]';
y5=[0 0 0 0 y]';
A=[y1 y2 y3 y4 y5];
D=cov(A);			% matrice de covariance / d'autocorr�lation 5*5
%
% Calcul des coefficients par Yule-Walker
E=-inv(D(1:4,1:4))*D(1,2:5)'; % ordre 4
H=-inv(D(1:3,1:3))*D(1,2:4)'; % ordre 3
E1=[1 E']'; % vecteur de coefficients incluant a0(ordre 4)
H1=[1 H']'; % vecteur de coefficients incluant a0(ordre 3)
%
% Trac� des spectres estim�s
[h11,w]=freqz(1,E1,Np); % ordre 4
[h12,w]=freqz(1,H1,Np); % ordre 3
mag1=abs(h11);
mag2=abs(h12);
figure
semilogy(f,mag1,f,mag2,':r', 'LineWidth',4);% trac� du spectre
title(['Spectre / Calcul sur l intervalle [1 ',num2str(n) ']']);
legend('ordre4', 'ordre3');


%% Calculs locaux sur m1 donn�es avec m1 << n 
pause;
%
m1=256; % largeur d une trame d analyse
m2=350; % param�tre utilis� pour d�finir le d�but d une trame
%
% Calcul des coefficients des mod�les AR d 'ordre 3 et 4
%   � partir de la donn�e 1 sur m1 donn�es
%
%calcul de l' autocorrelation
y1=[y(1:m1) 0 0 0 0]';
y2=[0 y(1:m1) 0 0 0]';
y3=[0 0 y(1:m1) 0 0]';
y4=[0 0 0 y(1:m1) 0]';
y5=[0 0 0 0 y(1:m1)]';
A=[y1 y2 y3 y4 y5];
D=cov(A);             % matrice de covariance / d'autocorr�lation    
%
% Calcul des coefficients par Yule-Walker
E=-inv(D(1:4,1:4))*D(1,2:5)'; % ordre 4
H=-inv(D(1:3,1:3))*D(1,2:4)'; % ordre 3
E1=[1 E']'; % vecteur de coefficients incluant a0 (ordre 4)
H1=[1 H']'; % vecteur de coefficients incluant a0 (ordre 3)
[h11,w]=freqz(1,E1,Np);
[h12,w]=freqz(1,H1,Np);
mag1=abs(h11);
mag2=abs(h12);
figure
plot(y(1:m1));				%trac� de la trame (sous-s�rie)
title(['\{y(1) ,..., y(' num2str(m1) ')\}']);
pause
figure
semilogy(f,mag11,f,mag1,'-r',f,mag2,'--k');%trace du spectre
title(['Spectre / Calcul sur l intervalle [1 ',num2str(m1) ']']);
legend('vrai spectre-y1', 'ordre4', 'ordre 3');
pause;

%calcul des coefficients des mod�les AR d 'ordre 3 et 4
%� partir de la donn�e m2 sur m1 donn�es
y1=[ y(m2:m1+m2) 0 0 0 0]';
y2=[0 y(m2:m1+m2) 0 0 0]';
y3=[0 0 y(m2:m1+m2) 0 0]';
y4=[0 0 0 y(m2:m1+m2) 0]';
y5=[0 0 0 0 y(m2:m1+m2)]';
A=[y1 y2 y3 y4 y5];
D=cov(A);             % matrice de covariance / d'autocorr�lation    
E=-inv(D(1:4,1:4))*D(1,2:5)'; % ordre 4
H=-inv(D(1:3,1:3))*D(1,2:4)'; % ordre 3
E1=[1 E']'; % vecteur de coefficients incluant a0 (ordre 4)
H1=[1 H']'; % vecteur de coefficients incluant a0 (ordre 3)
[h11,w]=freqz(1,E1,Np);
[h12,w]=freqz(1,H1,Np);
mag1=abs(h11);
mag2=abs(h12);
figure
plot(y(m2:m1+m2-1));			%trac�
title(['\{y(' num2str(m2) ') ,..., y(' num2str(m2+m1-1) ')\}']);
pause
figure
semilogy(f,mag11,'r',f,mag12,'g',f,mag1,'--r',f,mag2,'--g');% trac� du spectre
title(['Spectre / Calcul sur l intervalle [' num2str(m2) ',' num2str(m2+m1-1) ']']);
legend('vrai spectre-y1', 'vrai spectre-y2', 'ordre4', 'ordre 3');
pause;

%calcul des coefficients des mod�les AR d 'ordre 3 et 4
%� partir de la donn�e m2+m1 sur m1 donn�?es
y1=[ y(m2+m1:2*m1+m2) 0 0 0 0]';
y2=[0 y(m2+m1:2*m1+m2) 0 0 0]';
y3=[0 0 y(m2+m1:2*m1+m2) 0 0]';
y4=[0 0 0 y(m2+m1:2*m1+m2) 0]';
y5=[0 0 0 0 y(m2+m1:2*m1+m2)]';
A=[y1 y2 y3 y4 y5];
D=cov(A);             % matrice de covariance / d'autocorr�lation    
E=-inv(D(1:4,1:4))*D(1,2:5)'; % ordre 4
H=-inv(D(1:3,1:3))*D(1,2:4)'; % ordre 3
E1=[1 E']'; % vecteur de coefficients incluant a0 (ordre 4)
H1=[1 H']'; % vecteur de coefficients incluant a0 (ordre 3)
[h11,w]=freqz(1,E1,Np);
[h12,w]=freqz(1,H1,Np);
mag1=abs(h11);
mag2=abs(h12);
figure
plot(y(m2+m1:2*m1+m2-1));				%trac�
title(['\{y(' num2str(m2+m1) ') ,..., y(' num2str(m2+2*m1-1) ')\}']);
pause
figure
semilogy(f,mag12,f,mag1,'--r',f,mag2,'--g');%trac� du spectre
title(['Spectre / Calcul sur l intervalle [' num2str(m2+m1) ',' num2str(m2+2*m1-1) ']']);
legend('vrai spectre-y2', 'ordre4', 'ordre 3');
pause

%calcul des coefficients des mod�les AR d 'ordre 3 et 4
%� partir de la donn�e 2*m2 sur m1 donn�es
y1=[ y(2*m2:m1+2*m2) 0 0 0 0]';
y2=[0 y(2*m2:m1+2*m2) 0 0 0]';
y3=[0 0 y(2*m2:m1+2*m2) 0 0]';
y4=[0 0 0 y(2*m2:m1+2*m2) 0]';
y5=[0 0 0 0 y(2*m2:m1+2*m2)]';
A=[y1 y2 y3 y4 y5];
D=cov(A);             % matrice de covariance / d'autocorr�lation    
E=-inv(D(1:4,1:4))*D(1,2:5)'; % ordre 4
H=-inv(D(1:3,1:3))*D(1,2:4)'; % ordre 3
E1=[1 E']'; % vecteur de coefficients incluant a0 (ordre 4)
H1=[1 H']'; % vecteur de coefficients incluant a0 (ordre 3)
[h11,w]=freqz(1,E1,Np);
[h12,w]=freqz(1,H1,Np);
mag1=abs(h11);
mag2=abs(h12);
figure
plot(y(2*m2:m1+2*m2-1));			%trac�
title(['\{y(' num2str(2*m2) ') ,..., y(' num2str(2*m2+m1-1) ')\}']);
pause
figure
semilogy(f,mag12,'b',f,mag1,'--r',f,mag2,'--g');%trac� du spectre
title(['Spectre / Calcul sur l intervalle [' num2str(2*m2) ',' num2str(2*m2+m1-1) ']']);
legend('vrai spectre-y2', 'ordre4', 'ordre 3');
pause

%calcul des coefficients des mod�les AR d 'ordre 3 et 4
%� partir de la donn�e 2*m2+m1 sur m1 donn�es
y1=[ y(2*m2+m1:2*m1+2*m2) 0 0 0 0]';
y2=[0 y(2*m2+m1:2*m1+2*m2) 0 0 0]';
y3=[0 0 y(2*m2+m1:2*m1+2*m2) 0 0]';
y4=[0 0 0 y(2*m2+m1:2*m1+2*m2) 0]';
y5=[0 0 0 0 y(2*m2+m1:2*m1+2*m2)]';
A=[y1 y2 y3 y4 y5];
D=cov(A);             % matrice de covariance / d'autocorr�lation    
E=-inv(D(1:4,1:4))*D(1,2:5)'; % ordre 4
H=-inv(D(1:3,1:3))*D(1,2:4)'; % ordre 3
E1=[1 E']'; % vecteur de coefficients incluant a0 (ordre 4)
H1=[1 H']'; % vecteur de coefficients incluant a0 (ordre 3)
[h11,w]=freqz(1,E1,Np);
[h12,w]=freqz(1,H1,Np);
mag1=abs(h11);
mag2=abs(h12);
figure
plot(y(2*m2+m1:2*m1+2*m2-1));			%trac� 
title(['\{y(' num2str(2*m2+m1) ') ,..., y(' num2str(2*m2+2*m1-1) ')\}']);
pause
figure
semilogy(f,mag12,'r',f,mag13,'g',f,mag1,'--r',f,mag2,'--g');%trac� du spectre
title(['Spectre / Calcul sur l intervalle [' num2str(2*m2+m1) ',' num2str(2*m2+2*m1-1) ']']);
legend('vrai spectre-y2', 'vrai spectre-y3', 'ordre4', 'ordre 3');
pause

%calcul des coefficients des mod�les AR d 'ordre 3 et 4
%� partir de la donn�e 2*m2+2*m1 sur m1 donn�es
y1=[ y(2*m2+2*m1:3*m1+2*m2) 0 0 0 0]';
y2=[0 y(2*m2+2*m1:3*m1+2*m2) 0 0 0]';
y3=[0 0 y(2*m2+2*m1:3*m1+2*m2) 0 0]';
y4=[0 0 0 y(2*m2+2*m1:3*m1+2*m2) 0]';
y5=[0 0 0 0 y(2*m2+2*m1:3*m1+2*m2)]';
A=[y1 y2 y3 y4 y5];
D=cov(A);             % matrice de covariance / d'autocorr�lation    
E=-inv(D(1:4,1:4))*D(1,2:5)'; % ordre 4
H=-inv(D(1:3,1:3))*D(1,2:4)'; % ordre 3
E1=[1 E']'; % vecteur de coefficients incluant a0 (ordre 4)
H1=[1 H']'; % vecteur de coefficients incluant a0 (ordre 3)
[h11,w]=freqz(1,E1,Np);
[h12,w]=freqz(1,H1,Np);
mag1=abs(h11);
mag2=abs(h12);
figure
plot(y(2*m2+2*m1:3*m1+2*m2-1));		%trac�
title(['\{y(' num2str(2*m2+2*m1) ') ,..., y(' num2str(2*m2+3*m1-1) ')\}']);
pause
figure
semilogy(f,mag13,f,mag1,'--r',f,mag2,'--g');%trac� du spectre
title(['Spectre / Calcul sur l intervalle [' num2str(2*m2+2*m1) ',' num2str(2*m2+3*m1-1) ']']);
legend('vrai spectre-y3','ordre4', 'ordre 3');




