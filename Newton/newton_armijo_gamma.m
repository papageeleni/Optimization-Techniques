%This is the code for the newton method with gamma selected based on armijo rule 
%author: Eleni Papageorgiou 
%aem:9530%
%Optimiazation techniques lab2%

clc;
clear;
close all;

%%
% function visualization
syms f x y;
f(x,y) = (x^5).*exp(-(x^2)-(y^2));
% fsurf(x,y,f)
% xlabel('x')
% ylabel('y')
% zlabel('f')


%%
% initial point

x0 = [0.0 -1 1];
y0 = [0.0 1 -1];
mk = 0;
a = 10e-4;
b = 0.2;
s = 0.5;
e = 0.001;
grad = gradient(f,[x,y]); %gradient calculation
grad_2 = jacobian(grad);

%%

for i=1:length(x0)
    disp(['initial point: ', num2str(x0(i)),',',num2str(y0(i))]);
    X = [];
    Y = [];
    k=1;
    X(k) = x0(i);
    Y(k) = y0(i);
    
    figure
    fcontour(f)
    grid on
    title(['initial point: x_0=',num2str(x0(i)),' y_0=',num2str(y0(i))])
    hold on
    
    while 1    %finds the point and the grad

        f_grad = 0;
        temp = 0;
        mk=0;    
        plot(X,Y,"o")

        f_grad = double(grad(X(k), Y(k)));%gradient for xk,yk
        
        if (abs(vpa(norm(f_grad))) < e) 
            break;
        end
        f_grad_2 = double(grad_2(X(k), Y(k)));% 2nd gradient for xk,yk

        dk = -f_grad_2^(-1)*f_grad;
        %for constand gamma
        gamma = s*b^mk;
        
        temp = [X(k); Y(k)] + gamma*dk;
        X(k+1) = temp(1);
        Y(k+1) = temp(2);
        %%armijo
        while vpa(f(X(k),Y(k)) - f(X(k+1),Y(k+1))) < -a*b^mk * s * dk.' * f_grad 
            mk = mk+1;
            gamma = s*b^mk;
            temp = [X(k); Y(k)] + gamma*dk;
            X(k+1) = temp(1);
            Y(k+1) = temp(2);
            if mk>15
                gamma = 0.5; 
                temp = [X(k); Y(k)] + gamma*dk;
                X(k+1) = temp(1);
                Y(k+1) = temp(2);
                break;
            end

        end 
        %%

        k = k+1;
        
    
    end
    disp(['iterations: ',num2str(k-1)]);
    hold off
end
