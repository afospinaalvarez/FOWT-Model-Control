% Graficas modelo turbina

figure(1)
plot(out.cp.Time,out.cp.Data,'k')
grid on
title('Coeficiente de potencia')
xlabel('Time [s]')
ylabel('Coeficiente de potencia')
figure(2)
plot(out.velocidad.Time,out.velocidad.Data,'k')
grid on
title('Velocidad de giro')
xlabel('Time [s]')
ylabel('Velocidad de giro [rad/s]')
figure(3)
plot(out.potencia.Time,out.potencia.Data/1000000,'k')
grid on
title('Potencia Eléctrica')
xlabel('Time [s]')
ylabel('Potencia Eléctrica [MW]')
ax=gca;
%ax.YAxis.Exponent=6;
figure(4)
plot(wave.time,wave.signals.values,'k')
grid on
title('Perfil de velocidad de viento')
xlabel('Time [s]')
ylabel('Velocidad de viento [m/s]')
figure(5)
plot(out.angulopala.Time,out.angulopala.Data,'k')
grid on
title('Control Ángulo paso de pala')
xlabel('Time [s]')
ylabel('Ángulo pala [°]')
figure(6)
plot(out.angulotorre.Time,out.angulotorre.Data,'k')
grid on
title('Ángulo de la torre')
xlabel('Time [s]')
ylabel('Ángulo Torre [°]')
figure(7)
plot(out.ganancia.Time,out.ganancia.Data,'k')
grid on
title('Ganancia coeficiente de potencia')
xlabel('Time [s]')
ylabel('Ganancia neta [%]')
figure(8)
plot(out.velocidadrpm.Time,out.velocidadrpm.Data,'k')
grid on
title('Velocidad de giro')
xlabel('Time [s]')
ylabel('Velocidad de giro [rpm]')

