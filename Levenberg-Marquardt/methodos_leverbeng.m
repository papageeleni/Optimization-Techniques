%Vellios Georgios Serafeim AEM:9471

function [X,Y,counter]= methodos_leverbeng(f, e, xarx, yarx)
syms x y;
gradf = gradient(f, [x,y]);
grad2f = hessian(f,[x,y]);
X=[xarx];
Y=[yarx];
counter = 0;


while norm(gradf(X(end), Y(end))) >= e
    %epilogi toy mk se eyros [0.01, 100]
    start = 0.01;
    step = 0.01;
    m = -2;
    for mk = start:step:100
        M = double(grad2f(X(end), Y(end))) + mk*[1 0; 0 1];
        [~,p] = chol(M);
        if p == 0
            m = mk;
            break;
        end
    end
    
    if m == -2
        disp("error");
        break;
    end
    
    %epilisi tiw e3isosis gia na bre8ei to dk
    A = double((grad2f(X(end), Y(end)) + m*[1 0; 0 1]));
    B = double((-gradf(X(end),Y(end))));
    dk = double(inv(A)*B);
    gk = 0.5;
    xnext = [X(end); Y(end)] + gk*dk
    X = [X xnext(1)];
    Y =[Y xnext(2)];
    counter = counter +1;

end



