%going to do some parameter checks
clearvars
close all 



alpha_val = 0.5;
beta_val = 1;
max_val = 2;
n = 3;
step = (max_val - alpha_val) / n;

k1_val = 0;
k2_val = 0;
d_val = 0;

p_step = 0.05;


C = {'k','b','r','g','c', 'm', [.5 .6 .7],[.8 .2 .6]}; 
tspan =  [0 1100];
init = [0.5 0.5];


% sol = ode45(@(t,y) my_system(t, y, 1, 1, 1, 0.5), tspan, init);
% baseline = log((deval(sol,1000,1)) + (deval(sol,1000,2)));

f1 = figure;

%%% This function does multiple iterations
for j = 1:n
    max_fitness = 0;
    max_k1 = 0;
    max_k2 = 0;
    max_d = 0;
    while d_val < 1.01
        while k2_val < 1.01
            while k1_val < 0.51
                sol = ode45(@(t,y) my_system(t, y, alpha_val, beta_val, k1_val, k2_val, d_val), tspan, init);
                fitness = log((deval(sol,1000,1)) + (deval(sol,1000,2)));
                if fitness > max_fitness
                    max_fitness = fitness;
                    max_k1 = k1_val;
                    max_k2 = k2_val;
                    max_d = d_val;
                end
                k1_val = k1_val + p_step;
            end
            k2_val = k2_val + p_step;
            k1_val = 0;
        end
        %% Second figure for transmitter frequency
        k2_val = 0;
        d_val = d_val + p_step;
    end
    % scatter(max_d, max_k, C{j}, 'filled');
    scatter3(max_k1, max_k2, max_d)
    hold on
    disp(alpha_val)
    disp(max_fitness)
    disp(max_k1)
    disp(max_k2)
    disp(max_d)
    alpha_val = alpha_val + step;
    d_val = 0;
    max_k1 = 0;
    max_k2 = 0;
    max_d = 0;
    max_fitness = 0;
end



%  title('Title here');
set(0, 'CurrentFigure', f1)
% xlabel('Returns on reproduction investment (\beta)');
%xlabel('Returns on public good (\alpha)');
zlabel('Transmitter switching rate');
% xlabel('Utiliser public good');
xlabel('Transmitter public good');
ylabel('Utiliser public good');
% legend('y_1','y_2', 'y_3');
xlim([-0.03 1.03])
ylim([-0.03 1.03])
zlim([-0.03 1.03])



return

function dydt = my_system(t,y,a,b,k1,k,d1)
    k2 = k;
    dydt = [   (((y(1).*k1.^a + y(2).*k2.^a)./(y(1)+y(2))).*((1-d1).*((1-k1).^b).*y(1)));
               (((y(1).*k1.^a + y(2).*k2.^a)./(y(1)+y(2))).*((d1.*((1-k1).^b)).*y(1)+((1-k2).^b).*y(2)));
           ];
    
end