%going to do some parameter checks
clearvars
close all 

k_val = 0.5;

alpha_val = 1;
beta_val = 1;
max_val = 1;
n = 50;
step = (max_val - k_val) / n;

d_val = 0.5;
d_len = 5;
d_step = 0.1;

C = {'k','b','r','g','c',[.5 .6 .7],[.8 .2 .6]}; 
tspan =  [0 1500];
init = [0.5 0.5];


[t,y] = ode45(@(t,y) my_system(t, y, 1, 1, 1, 0.5), tspan, init);
baseline = y(500,1) + y(500,2);



% One parameter iteration
for i = 1:n
    [t,y] = ode45(@(t,y) my_system(t, y, k_val, alpha_val, beta_val, d_val), tspan, init);
    if (y(500,1)/(y(500,1) + y(500,2))) > 0.001
        scatter(k_val, y(500,1)/(y(500,1) + y(500,2)), 'b');
        hold on
    end
    k_val = k_val + step;
end

alpha_val = 1.1;
k_val = 0.5;

% for i = 1:n
%     [t,y] = ode45(@(t,y) my_system(t, y, k_val, alpha_val, beta_val, d_val), tspan, init);
%     if (y(500,1)/(y(500,1) + y(500,2))) > 0.001
%         scatter(k_val, y(500,1)/(y(500,1) + y(500,2)), 'g');
%         hold on
%     end
%     k_val = k_val + step;
% end





% Two parameter iteration
% for j = 1:d_len
%     for i = 1:n
%         [t,y] = ode45(@(t,y) my_system(t, y, k_val, alpha_val, beta_val, d_val), tspan, init);
%         if (y(500,1)/(y(500,1) + y(500,2))) > 0.001
%             scatter(alpha_val, y(500,1)/(y(500,1) + y(500,2)), C{j});
%             hold on
%         end
%         alpha_val = alpha_val + step;
%     end
%     k_val = k_val + d_step;
%     alpha_val = 0.5;
% end

%  title('Title here');
% xlabel('Returns on reproduction investment (\beta)');
% xlabel('Returns on public good production (\alpha)')
xlabel('Utiliser public good');
ylabel('Growth rate');
% ylabel('Transmitter frequency');
% legend('y_1','y_2', 'y_3');
% ylim([0 0.75])
% xlim([0.4 2.1])



return

function dydt = my_system(t,y,k,a,b,d1)
    k2 = k;
    dydt = [   (((y(2).*k2.^a)./(y(1)+y(2))).*((1-d1).*y(1)));
               (((y(2).*k2.^a)./(y(1)+y(2))).*(d1.*y(1)+((1-k2).^b).*y(2)));
           ];
    
end