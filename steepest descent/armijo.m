%This is the code for the steepest descent method with gamma selected based on armijo rule 
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
        plot(X,Y,"o")

        f_grad = 0;
        temp = 0;

        f_grad = double(grad(X(k), Y(k)));%gradient for xk,yk
        
        if (abs(vpa(norm(f_grad))) < e) 
            break;
        end
        dk = -f_grad;
        %for constand gamma
        mk = 0;

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
        end 
        %%

        k = k+1;
        
    
    end
    disp(['iterations: ',num2str(k-1)]);
    hold off
end
