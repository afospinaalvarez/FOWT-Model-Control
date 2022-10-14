% Parámetros de entrada turbina eólica flotante
% Trabajo fin de máster
% Andrés Felipe Ospina

clear
clc

% Parámetros modelo de viento
v=11; %Velocidad media de viento

% Parámetros turbina
cuerda=3; % Longitud de cuerda de la pala [m]
densidad=1.25; % Densidad del aire en condiciones de operación [kg/m^3]
radio=35; % Radio de la pala dada en [m]
area=pi*radio^2;
angulo_pala=5; % Ángulo de pala en caso de que se simule aerogenerador de pala fija
betta=angulo_pala*pi/180; %Angulo en radianes de la pala
m=0.5*densidad*radio^3*pi;
N=83.5; % Relación de transimisión
Jt=4.4535e5; % Momento de inercia total 
Bt=2836.49; % Coeficiente de fricción total

% Parámetros de la torre y barcaza
Ibase=1.766e9; % Momento de inercia de la base de la turbina flotante [kg*m^2]
Itorre=3.3428e9; % Momento de inercia de la turbina de la turbina flotante [kg*m^2]
Itotal=Ibase+Itorre; % Momento de inercia total [kg*m^2]
altura=85; % Áltura de la torre [m]
kbase=1.888e9; % Coeficiente elástico de la base [kg*m^2*s^-2]
ktorre=1.2519e9; % Coeficiente elástico de la torre [kg*m^2*s^-2]
ktotal=kbase+ktorre; % Coeficiente elástico total [kg*m^2*s^-2]
dbase=5.123e7; % Coeficiente de fricción de la base [kg*m^2*s^-1]
dtorre=2.869e7; % Coeficiente de fricción de la torre[kg*m^2*s^-1]
dtotal=dbase+dtorre; % Coeficiente de fricción total [kg*m^2*s^-1]

% Parámetros Generador DFIG 1.5 MW

omega_s=2*pi*60; % Frecuencia angular de la red
Ls=1.61598e-3; % Inductancia del estator
Lr=1.608088e-3; %Inductancia del rotor
Lm=1.526e-3; % Inductancia de acoplamient
Le=Ls*Lr-Lm^2;
Rr=0.99187e-3; % Resistencia del rotor
Vs=460*(2/3)^0.5; %Voltaje del generador   
K=Ls*Rr/Le;
p=2; %Número de polos del generador

% Parámetros del filtro

%Spatial filtering
alfa=0.55; % [s^-1] 
decay_factor=1.3;
Vtprom=v; % Campo de velocidad de viento promedio [m/s]
betaf=decay_factor*radio/v;

%Lag filter
alfai=1.17;
thaui=9; % Constante de tiempo del filtro de retardo [s]

% Modelo viento de mediano-largo plazo

Time=1800; % Tiempo de simulación [s]
Ts1=180; % Tiempo de muestreo de la velocidad media [m/s]
Ts=1; % Tiempo de muestreo de la simulación de viento [m/s]
Vm=v; % Velocidad de viento promedio [m/s]
sigma=1.3; % Desviación estándar de la turbulencia [m/s]
L=75.30; % Longitud de la turbulencia [m]
ksigmav=0.15; % Parámetro experimental (Turbinas offshore varía entre 0.1-0.15)
fk=[];M=[]; 
for k=-3:2
   for i=1:9
       fi=i*10^k;
       fk(i)=fi;
   end
   M((k+4),:)=[fk];
end
f=horzcat(M(1,:),M(2,:),M(3,:),M(4,:),M(5,:),M(6,:));
f=f/3600; % Frecuencia discreta de la señal [Hz]
w=2*pi*f; % Frecuencia angular [rad/s]
Svv=[];A=[];Vmli=[];
a=0;b=2*pi; % Rango de la función aleatoria phi
Nv=30; % Número de iteraciones de la función de velocidad media
for k=1:Time/Ts1-1
    Vmli=Vm;
    Vml(1)=Vm;
    for i=1:Nv+1
        Svv(i)=(0.475*sigma^2*(L/Vm))/(1+(w(i)*L/Vm)^2)^(5/6);
    end
    for i=1:Nv
        A(i)=2/pi*sqrt(1/2*(Svv(i)+Svv(i+1))*(w(i+1)-w(i)));
        Vmli=A(i)*cos(w(i)*k+(a+(b-a)*rand(1)))+Vmli;
    end
    Vml(k+1)=Vmli;
end

Vmedia=Vml; 
m1=0.4;m2=0.25; % Constantes de la función de transferencia
cnt=[];wn=[];
for k=1:(length(Vml))
    Tf=L/Vmedia(k);
    Kf=sqrt((2*pi*Tf)/(beta(1/2,1/3)*Ts));
    s=tf('s');
    Hf=Kf*(m1*Tf*s+1)/((Tf*s+1)*(m2*Tf*s+1)); % Función de transferencia
    h=impulse(Hf,[0:Ts:Ts1-1]); % Respuesta del filtro ante una señal impulso
    wn(k,:)=wgn(length(h),1,0); % Ruido blanco
    cn=conv(h,wn(k,:)); % Convolución de la componente turbulenta de la velocidad 
    for i=1:length(wn(k,:))
        cnt(k,i)=cn(i);
    end
end
A=wn;B=A';wn=B(:)';
A=cnt;B=A';cnt=B(:)';
sigmav=[];Vt=[];
for k=0:length(Vml)-1
    for t=(k*Ts1+1):(k+1)*Ts1
        sigmav(k+1)=ksigmav*Vmedia(k+1);
        Vt(t)=sigmav(k+1).*cnt(t);
    end
end
V=[];
for k=0:length(Vml)-1
    for t=(k*Ts1+1):(k+1)*Ts1
        V(t)=Vml(k+1)+Vt(t); % Velocidad de viento
    end
end
t = 1*[0:Time-1]'; 
wave.time=[t]; % Señal de tiempo para simulimk
wave.signals.values=[V]'; % Señal de viento para simulink


% Componente turbulenta modelo de corto plazo > 600 s

Nt=50; % Número de iteraciones
V50=[];
Vwnt=[];
Dw=0.5; % Paso frecuencia
F=500; % longitud de turbulencia
mu=v;
Kn=0.004; % Coeficiente de arrastre de la superficie
Tsg=10; % Tiempo inicio ráfaga
Teg=20; % Tiempo final ráfaga
Tsr=30; % Tiempo inicio de la rampa
Ter=50; % Tiempo final rampa
Ag=4; % Amplitud ráfaga
Ar=4; % Amplitud rampa



for tt=0:179
    Vwn=0;
    at=0;bt=2*pi; % Rango de la variable aleatoria
for i=1:Nt
    wt(i)=(i-0.5)*Dw;
    Sv(i)=(2*Kn*F^2*abs(wt(i)))/(pi^2*(1+((F*wt(i))/(mu*pi))^2)^(4/3)); % Densidad espectral
    Vwn=2*(((Sv(i)*Dw)^(0.5))*cos(wt(i)*tt+(at+(bt-at)*rand(1))))+Vwn; % Componente turbulenta
    V50(i)=Vwn;
end
Vwnt(tt+1)=Vwn;
end

tt = 1*[0:179]';
turbulent.time=[tt];
turbulent.signals.values=[Vwnt]';
turbulent.signals.dimensions=1;
