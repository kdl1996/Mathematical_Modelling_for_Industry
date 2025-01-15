function ret = nf( u,v ) 

for i = 1:length(u) 
    if (f(u(i)) >= 0 && f(v(i))>=0)
        ustar(i)=u(i);
    else if (f(u(i)) < 0 && f(v(i))<0)
            ustar(i)=v(i); 
    else if (f(u(i)) >= 0 && f(v(i))<0)
            s=(f(v(i))-f(u(i)))/(v(i)-u(i));
            if (s>=0)
                ustar(i)=u(i); 
            else
                ustar(i)=v(i); 
            end
        else if((f(u(i)) < 0 && f(v(i))<=0))
               ret= 0;
               
            end
            
        
                
        end
        end
    end
end
         
ret =f(ustar);


%{
for i = 1:length(u) 
    if (u(i) >= v(i))
        if ((u(i)+v(i))/2 > 0) 
            ustar(i)=u(i);
        else
            ustar(i)=v(i); 
        end
    else if (u(i)>0) 
            ustar(i)=u(i);
        elseif (v(i)<0)
            ustar(i)=v(i);
        else ustar(i)=0;
        end
    end
end
ret =f(ustar);
%}

