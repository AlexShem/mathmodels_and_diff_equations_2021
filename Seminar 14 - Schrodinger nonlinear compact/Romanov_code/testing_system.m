
D = 0.01;
tau = 2*pi/256;
h = 2*pi/128;
n1 = ceil(2*pi/h);
m1 = ceil(2*pi/tau);
%9-26-36
%10-27-37

[U1,V1, ep1, ep2, integral10, integral27, integral37, anss, spectr2] = schredinger_MAK_eps(h,  D, tau, y1, x1);
% % % % 
% [U2,V2, ep11, ep22, integral1, integral2,integral3, test1] = schredinger_MAK(h, D, tau, y1, x1);

%[U1, FF] = Burgers_MAK(h, D, tau);
%[U2, FF] = Burgers_Mak_eps_sglas(h, D, tau);
%[U1, ep1, ep2, max_epsilon1, max_epsilon2] = Burgers_AB_eps(h, D, tau);
%[U1, ep1, ep2, max_epsilon1, max_epsilon2, max_epsilon3] = Burgers_Mak_eps(h, D, tau);
%[U2, ep1, ep2, max_epsilon1, max_epsilon2, max_epsilon3] = Burgers_Mak_eps(h, D, tau); 
%[U, FF]  = Burgers(h, D, tau);
%[U1, ep1, ep2] = Burgers_full(h, D, tau);
%[U1, FF]  = Burgers_full(h, D, tau);

% X = linspace(0,2.3, length(y1));
% X = X(1:length(y1)-1)';
% 
% for k = 1:1000
%     
%     if mod(k,200) == 0
%         
%         %u_etal = -2*D*cos(X)./(2*exp(D*tau*(k-1))+sin(X));
%         
%         plot(X,  U1(:, k)-U2(:, k));
%         hold on
% 
%     end
%     
% end

Y = linspace(1, 1025, 1025);
plot(Y, integral10);

% Y = linspace(1,125,125)*8*pi/125;
% plot(Y(2:125), integral26(2:125));
% hold on
% plot(Y(2:125), integral27(2:125));
% hold on
% plot(Y(2:125), integral2(2:125));
% hold on

% plot(Y(2:125), integral1(2:125));
% hold on
%plot(Y(1:2000), integral6(1:2000));
%hold on
% plot(Y(2:16000), spectr(2:16000));
% hold on
% plot(Y(2:16000), spectr1(2:16000));
% hold on
% plot(Y(2:16000), spectr2(2:16000));
% hold on
% plot(Y(1:16000), integral17(1:16000));
% hold on
% plot(Y(1:16000), integral15(1:16000));
% hold on
% plot(Y(1:4000), integral8(1:4000));
% hold on
% % 
% plot(Y(1:4000), integral10(1:4000));
% hold on
% plot(Y(1:4000), integral11(1:4000));
% hold on
% plot(Y(1:4000), integral12(1:4000));
% hold on
% plot(Y(1:4000), integral13(1:4000));
% hold on

% plot(Y(1:4000), integral3(1:4000));


% xlabel('x')
% ylabel('error')
% %legend('time = 1t/5', 'time = 2t/5', 'time = 3t/5', 'time = 4t/5')
% legend('Compact АБ', 'АБ', 'Compact МК', 'МК')

title(['D =',num2str(D),', tau =',num2str(tau),', h =',num2str(h), ', I(x)'])

% hold off;
%xlabel('Период до 2.3')
xlabel('x')
ylabel('I(x)')
% legend('Compact АБ', 'АБ')
legend( 't = 0.3', 't = 0.6', 't = 0.9', 't = 1.2', 't = 1.5', 't= 1.6')
legend('t = 1.25', 't = 2.5', 't = 3.75', 't = 5', 't = 6.25')
legend('Первый интеграл для компактной МК');
% legend('Первый интеграл компактн', 'Первый интеграл явн', 'Компактная, a = 1/4', 'a =  -1/4', 'Второе значение')
% legend('Первый интеграл явн', 'Первый интеграл компакт a = -1/4', 'Компактн a = 0', 'Компактн a = -1', 'Компактн a = -10')
% %legend('I(t), компактная, производная - компакт ', 'I(t), явная','I(t), компактная, производная -центральные ' )
% legend('Комп a = -1', 'Явная', 'Комп a = -5', 'Комп a = 1','Комп a = -0.53', 'Комп a = -0.52', 'Комп a = -0.51', 'a = -0.501'); 
% legend('компа a = -0.49', 'комп a = -0.51', 'Явн')
% legend('t = 9.8', 't = 19.63');
%legend('0.075', '0.15', '0.225', '0.3', '0.375', '0.45', '0.525', '0.6')
% Y= linspace(1, 4001, 4001);
% ep1 = ep1(1:4001);
% plot(Y, ep1);
% hold on
% plot(Y, max(ep2));
%legend('компактная')
