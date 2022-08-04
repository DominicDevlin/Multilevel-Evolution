clearvars
close all 

alpha_val = 1.5;
beta_val = 1;
k_val = 0.5;
d_val = 1;


tspan =  [0 1100];
init = [0.5 0.5];

sol = ode45(@(t,y) my_system(t, y, alpha_val, beta_val, k_val, d_val), tspan, init);
fitness = log((deval(sol,1000,1)) + (deval(sol,1000,2)));
disp(fitness)

% figure; semilogy(t,y(:,1), '-')
% hold on; semilogy(t,y(:,2), '--')
% title('Title here');
% xlabel('time');
% ylabel('solution');
% legend('y_1','y_2');

% figure; plot(k_val, y1vals)
% hold on; plot(t, y2vals)
% hold on; plot(t, y3vals)
%figure; plot(t,y(:,1)./(y(:,1) + y(:,2)),y(:,3)./(y(:,3) + y(:,4)))
% y(:,2),y(:,3),y(:,4),y(:,5), '-');
% figure; plot(t,y(:,1)./(y(:,1) + y(:,2)), '-')
%hold on; plot(t,y(:,2), '--')
% title('Title here');
% xlabel('parameter');
% ylabel('transmitter frequency');
% legend('y_1','y_2', 'y_3');
% ylim([0 1])



return

function dydt = my_system(t,y,a,b,k,d1)
    k2 = k;
    dydt = [   (((y(2).*k2.^a)./(y(1)+y(2))).*((1-d1).*y(1)));
               (((y(2).*k2.^a)./(y(1)+y(2))).*(d1.*y(1)+((1-k2).^b).*y(2)));
           ];
    
end


%% logistic function
% function dydt = my_system(t,y)
%     k2 = 0.7;
%     d1 = 0.5;
%     a = 1;
%     dydt = [   (((y(2)*k2^a)/(y(1)+y(2)))*((1-d1)*y(1)))*(1-(y(1)+y(2)));
%                (((y(2)*k2^a)/(y(1)+y(2)))*(d1*y(1)+(1-k2)*y(2)))*(1-(y(1)+y(2)));
%            ];
% 
% %     dydt = [   d1*((1-d1)*y(1))*(1-(y(1)+y(2)));
% %                d1*(d1*y(1)+(1-k2)*y(2))*(1-(y(1)+y(2)));
% %            ];
%     
% end

%% birth-death function - can't seem to get it working
% function dydt = my_system(t,y)
%     k2 = 0.8;
%     d1 = 0.5;
%     a = 1;
%     dydt = [   (k2^a)*(y(2)/(y(1)+y(2)))*((1-d1)*y(1)) - (((y(2)*k2^a)/(y(1)+y(2)))*((1-d1)*y(1)) + ((y(2)*k2^a)/(y(1)+y(2)))*(d1*y(1)+(1-k2)*y(2)))/(y(1)+y(2));
%                (k2^a)*(y(2)/(y(1)+y(2)))*(d1*y(1)+(1-k2)*y(2)) - (((y(2)*k2^a)/(y(1)+y(2)))*((1-d1)*y(1))+ ((y(2)*k2^a)/(y(1)+y(2)))*(d1*y(1)+(1-k2)*y(2)))/(y(1)+y(2));
%            ];
%     
% end

%fractional






