clearvars
close all

t = 1000;

k = 0.5:0.001:1;
y1 = exp(0.25*t);
y2 = -0.25./(0.25 - (0.5*k))*exp(0.25*t) + ((0.5-k)./(0.25-k)).*exp(0.5*(1-k)*t);

%figure 1 
figure
plot(k, y1./(y1+y2))
hold on
plot(k, y2./(y1+y2))

ax = gca;
h = findobj(gca, 'Type', 'Line');

x = h.XData;
y = h.YData;


title ('Optimal frequency of transmitters and utilisers depending on utiliser public good')
xlabel('Utiliser public good')
ylabel('Frequency')
legend('transmitters', 'utilisers')
legend('location', 'southeast')
set(gcf, 'PaperPositionMode', 'auto')
axis([0.5 1 0 1])



return