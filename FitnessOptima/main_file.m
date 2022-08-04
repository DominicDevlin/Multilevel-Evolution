%going to do some parameter checks
clearvars
close all 





alpha_val = 1;
beta_val = 1;
max_val = 2;
n = 2;
step = (max_val - alpha_val) / n;

k_val = 0.02;
d_val = 0.5;
len = 50;
d_step = (1 - d_val) / len;
k_step = (1 - k_val) / len;

C = {'k','b','r','g','c', 'm', [.5 .6 .7],[.8 .2 .6]}; 
tspan =  [0 2000];
init = [0.5 0.5];


sol = ode45(@(t,y) my_system(t, y, 1, 1, 1, 0.5), tspan, init);
baseline = log((deval(sol,1000,1)) + (deval(sol,1000,2)));

x_list = zeros(0, len * n);
y_list = zeros(0, len * n);
x1_list = zeros(0, len * n);
y1_list = zeros(0, len * n);

% f1 = figure;
% f2 = figure;

% One parameter iteration
% for i = 1:len
%     sol = ode45(@(t,y) my_system(t, y, alpha_val, beta_val, k_val, d_val), tspan, init);
%     growth = log((deval(sol,1000,1)) + (deval(sol,1000,2)));
%     set(0, 'CurrentFigure', f1)
%     scatter(k_val, growth / baseline, 'b', 'filled');
%     hold on
%     %% Second figure for transmitter frequency
%     set(0, 'CurrentFigure', f2)
%     scatter(k_val, (deval(sol,1000,1)/(deval(sol,1000,1) + deval(sol,1000,2))), 'g', 'filled');
%     hold on
%     k_val = k_val + k_step;
    % beta_val = beta_val + step;
% end


%%% This function does multiple iterations
for j = 1:n
    for i = 1:len
        sol = ode45(@(t,y) my_system(t, y, alpha_val, beta_val, k_val, d_val), tspan, init);
        growth = log((deval(sol,1000,1)) + (deval(sol,1000,2)));
%         set(0, 'CurrentFigure', f1)
%         scatter(k_val, growth / baseline, C{j}, 'filled');
        x_list(len*(j-1)+i) = k_val;
        y_list(len*(j-1)+i) = growth / baseline;
%         hold on
        %% Second figure for transmitter frequency
%         set(0, 'CurrentFigure', f2)
%         scatter(k_val, (deval(sol,1000,1)/(deval(sol,1000,1) + deval(sol,1000,2))), C{j}, 'filled');
%         hold on
        x1_list(len*(j-1)+i) = k_val;
        y1_list(len*(j-1)+i) = (deval(sol,1000,1)/(deval(sol,1000,1) + deval(sol,1000,2)));
        k_val = k_val + k_step;
        % beta_val = beta_val + step;
    end
    alpha_val = 1;
    beta_val = 1.1;
    k_val = 0.02;
end




%  title('Title here');
% set(0, 'CurrentFigure', f1)
% xlabel('Returns on reproduction investment (\beta)');
%xlabel('Returns on public good (\alpha)');
% xlabel('Transmitter switching rate')
% xlabel('Utiliser public good')
% ylabel('Growth rate');
% legend('y_1','y_2', 'y_3');
% xlim([-0.03 1.03])
% ylim([-0.03 1.6])


% set(0, 'CurrentFigure', f2)
% xlabel('Transmitter switching rate')
% xlabel('Utiliser public good')
% ylabel('Transmitter frequency');
% legend('y_1','y_2', 'y_3');
% ylim([0 1])
% xlim([-0.02 1.02])


f1 = figure;
set(0, 'CurrentFigure', f1)
xlabel('Utiliser public good');
ylabel('Growth rate');
plot(x_list(1:len),y_list(1:len))
hold on
plot(x_list(len+1:end),y_list(len+1:end))
ylim([-0.02 1.02])
xlim([0 1])


f2 = figure;
set(0, 'CurrentFigure', f2)
xlabel('Utiliser public good')
ylabel('Transmitter frequency')
plot(x1_list(1:len),y1_list(1:len))
hold on
plot(x1_list(len+1:end),y1_list(len+1:end))
ylim([-0.01 0.51])
xlim([0 1])

return

function dydt = my_system(t,y,a,b,k,d1)
    k2 = k;
    dydt = [   (((y(2).*k2.^a)./(y(1)+y(2))).*((1-d1).*y(1)));
               (((y(2).*k2.^a)./(y(1)+y(2))).*(d1.*y(1)+((1-k2).^b).*y(2)));
           ];
    
end