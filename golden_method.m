function golden_method=golden_method(h,a,b,l)

%syms f x y g;

gamma=0.618;


        %%
        %initialization
        x1=[];
        x2=[];
        k=1;

        A(k)=a;
        B(k)=b;
        x1(k)=A(k)+(1-gamma)*(B(k)-A(k));
        x2(k)=A(k)+gamma*(B(k)-A(k));
        
    
        phi1 = vpa(h(x1(k)));
        phi2 = vpa(h(x2(k)));
    
       while b(k)-a(k)>=l
    
            if phi1>phi2 %f(x1)>f(x2)

                a(k+1)=x1(k);
                b(k+1)=b(k);
                x2(k+1)=a(k+1)+gamma*(b(k+1)-a(k+1));
                x1(k+1)=x2(k);
                phi1=phi2;
                phi2 = vpa(h(x2(k+1)));
                
            else  % f(x1)<f(x2)

                a(k+1)=a(k);
                b(k+1)=x2(k);
                x2(k+1)=x1(k);
                x1(k+1)= a(k+1)+(1-gamma)*(b(k+1)-a(k+1));
                phi2=phi1;
                phi1 = vpa(h(x1(k+1)));
        
            end
    
            k=k+1;
        
       end
    
       golden_method= (x1(end)+x2(end))/2;
   end