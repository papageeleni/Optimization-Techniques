%This is the code for the Levenberg-Marquardt method with gamma constand
%author: Eleni Papageorgiou 
%aem:9530%
%Optimiazation techniques lab2%

clc;
clear;
close all;

%%
% function visualization
syms f x y mk;
f(x,y) = (x^5).*exp(-(x^2)-(y^2));
% fsurf(x,y,f)
% xlabel('x')
% ylabel('y')
% zlabel('f')
%%
% initial point
x0 = [0.0 -1 1];
y0 = [0.0 1 -1];
gamma = 0.4;
e = 0.001;
I = eye(2);


grad = gradient(f,[x,y]); % gradient calculation
grad_2 = jacobian(grad); %jacobian of gradient
    
for i=1:length(x0)
    disp(['initial point: ', num2str(x0(i)),',',num2str(y0(i))]);

    X=[];
    Y=[];
    k=1;
    X(k) = x0(i);
    Y(k) = y0(i);
    
    figure
    fcontour(f)
    grid on
    title(['initial point: x_0=',num2str(x0(i)),' y_0=',num2str(y0(i))])
    hold on

    while vpa(norm(grad(X(k), Y(k)))) >= e

        plot(X,Y,"o")
        
        step = 0.1;
        f_grad_2 = vpa(grad_2(X(k), Y(k)));
        d = eig(f_grad_2);
        m = max(abs(d))+1;
        
        A = (vpa(grad_2(X(k), Y(k))) + m*I);
        B = vpa(-grad(X(k),Y(k)));
        dk = double(A^(-1)*(B));

        temp = [X(k); Y(k)] + gamma*dk;
        X(k+1) = temp(1);
        Y(k+1) = temp(2);
        k = k+1; 
    end
    disp(['iterations: ',num2str(k-1)]);

    hold off

end
%plot(X(end),Y(end),'o')

%%
