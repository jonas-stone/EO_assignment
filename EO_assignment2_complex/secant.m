function [f_opt, x_opt] = secant(fun, x1, x2, tol, max_iter, max_limit, low_limit)
    x(1) = x1;
    x(2) = x2;
    i = 2;
    f_prev = fun(x(1));
    f_curr = fun(x(2));
    while abs(f_curr) > tol && abs(x(i)-x(i-1)) > tol && i < max_iter
        denom = f_curr-f_prev;
        if abs(denom) < 1e-12
            break;
        end

        x(i+1)=x(i)-f_curr*((x(i)-x(i-1))/denom);
        
        if x(i+1) > max_limit
            x(i+1) = max_limit;
        elseif x(i+1) < low_limit
            x(i+1) = low_limit;
        end   

        f_prev = f_curr;            
        f_curr = fun(x(i+1));
        i = i+1;
    end
    f_opt = f_curr;
    x_opt = x(i);
end