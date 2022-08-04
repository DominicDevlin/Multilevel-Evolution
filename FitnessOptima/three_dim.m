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
d_val = 0.3;

k1_len = 20;

len = 20;
p_step = 0.05;

C = {'k','b','r','g','c', 'm', [.5 .6 .7],[.8 .2 .6]}; 
tspan =  [0 2000];
init = [0.5 0.5];


sol = ode45(@(t,y) my_system(t, y, 1, 1, 0, 1, 0.5), tspan, init);
baseline = log((deval(sol,1000,1)) + (deval(sol,1000,2)));

f1 = figure;
f2 = figure;

x = sort(rand(1, k1_len * len * n), 'ascend');
y = sort(rand(1, k1_len * len * n), 'ascend');
z = sort(rand(1, k1_len * len * n), 'ascend');

%%% This function does multiple iterations
for j = 1:k1_len
    for k = 1:len
        for i = 1:n
            sol = ode45(@(t,y) my_system(t, y, alpha_val, beta_val, k1_val, k2_val, d_val), tspan, init);
            growth = log((deval(sol,1000,1)) + (deval(sol,1000,2)));
            growth = growth / baseline;
            set(0, 'CurrentFigure', f1)
            scatter3(k1_val, k2_val, growth, C{i}, 'filled');
            hold on
            %% Second figure for transmitter frequency
            set(0, 'CurrentFigure', f2)
            scatter3(k1_val, k2_val, (deval(sol,1000,1)/(deval(sol,1000,1) + deval(sol,1000,2))), C{i}, 'filled');
            hold on
            
            x((j - 1)*n*len + (k-1)*n + i) = k1_val;
            y((j - 1)*n*len + (k-1)*n + i) = k2_val;
            z((j - 1)*n*len + (k-1)*n + i) = growth;
            
            alpha_val = alpha_val + step;
            % beta_val = beta_val + step;
        end
        k2_val = k2_val + p_step;
        alpha_val = 0.5;
        beta_val = 0.5;
    end
    k2_val = 0.0;
    k1_val = k1_val + p_step;
    disp(k1_val)
end


%%%% This function does a single iteration
% for i = 1:n
%     set(0, 'CurrentFigure', f1)
%     sol = ode45(@(t,y) my_system(t, y, alpha_val, beta_val, k_val, d_val), tspan, init);
%     growth = log(deval(sol,1000,1)) + log(deval(sol,1000,2));
%     scatter(alpha_val, growth / baseline);
%     hold on
%     % set(0, 'CurrentFigure', f2)
%     % scatter(d_val, (deval(sol,1000,1)/(deval(sol,1000,1) + deval(sol,1000,2))));
%     % hold on
%     alpha_val = alpha_val + step;
% end


%  title('Title here');
set(0, 'CurrentFigure', f1)
% xlabel('Returns on reproduction investment (\beta)');
%xlabel('Returns on public good (\alpha)');
xlabel('Transmitter public good')
ylabel('Utiliser public good')
% xlabel('Utiliser public good')
zlabel('Growth rate');
% legend('y_1','y_2', 'y_3');
xlim([-0.03 1.03])
ylim([-0.03 1.03])
zlim([-0.03 2])


set(0, 'CurrentFigure', f2)
%xlabel('Transmitter switching rate')
xlabel('Transmitter public good')
ylabel('Utiliser public good')
zlabel('Transmitter frequency')
% legend('y_1','y_2', 'y_3');
zlim([0 1])
xlim([-0.02 1.02])
ylim([-0.02 1.02])

% [X,Y,Z] = meshgrid(x,y,z);



f3 = figure;
set(0, 'CurrentFigure', f3)
tri = delaunay(x,y);
trisurf(tri,x,y,z)
axis vis3d
% l = light('Position',[-50 -15 29]);
% set(gca,'CameraPosition',[0.4 1 2])
zlim([0 2])
xlim([-0.02 0.4])
ylim([0.4 1.02])


% plot3(x, y, z) 
% xlim([-0.03 1.03])
% ylim([-0.03 1.03])
% zlim([-0.03 2])

return

function dydt = my_system(t,y,a,b,k1,k,d1)
    k2 = k;
    dydt = [   (((y(1).*k1.^a + y(2).*k2.^a)./(y(1)+y(2))).*((1-d1).*((1-k1).^b).*y(1)));
               (((y(1).*k1.^a + y(2).*k2.^a)./(y(1)+y(2))).*((d1.*((1-k1).^b)).*y(1)+((1-k2).^b).*y(2)));
           ];
    
end