function [c_root_best, c_tip_best] = simplex(fun, x0, ub, lb, max_iter, tol_f, tol_x)
    simplex_points = 3;
    F = zeros(simplex_points, 1);
    X = zeros(simplex_points, 2);
    
    X(1, :) = x0; 
    X(2, :) = [x0(1)*1.05, x0(2)];
    X(3, :) = [x0(1), x0(2)*1.05];
    
    for i = 1:3
        if X(i, 1) < lb(1) || X(i, 2) < lb(2) || X(i, 1) > ub(1) || X(i, 2) > ub(2)
            F(i)= 1e6;
        else
            F(i)= fun([X(i, 1), X(i, 2)]);
        end
    end
    
    for iter = 1:max_iter
        % Sorting
        [F_new, I] = sort(F);
        F = F_new;
        X = X(I, :);
        m = (X(1, :) + X(2, :)) / 2;

        f_dist = max(F) - min(F);
        x_dist = max([norm(X(1,:) - X(2,:)), norm(X(1,:) - X(3,:))]);
        
        if f_dist < tol_f && x_dist < tol_x
            fprintf('Convergence reached at iteration %d!\n', iter);
        break; % This exits the for-loop immediately
        end

        % Reflecting
        r = 2*m - X(3,:);
        if  r(1) < lb(1,1) || r(2) < lb(1,2) || r(1) > ub(1,1) || r(2) > ub(1,2)
            f_r = 1e6;
        else
            f_r = fun(r);
        end

        if f_r >= F(1) && f_r < F(2)
            % --- STANDARD REFLECTION ---
            % ACTION: Accept 'r'. 
            % TODO: Overwrite the worst point (X(3,:) and F(3)) with r and f_r.
            X(3,:) = r;
            F(3) = f_r;
            
        elseif f_r < F(1)
            % --- EXPANSION ---
            % 'r' is a new global best! We are heading in a great direction.
            % Let's push even further in that same direction to a point 's'.
            % ACTION: Calculate 's', evaluate 'f_s' (remember bounds!), 
            % then compare f_s and f_r to see which one to accept.
            % TODO: Write the Expansion math here.
            s = m + 2*(m-X(3,:));
            if  s(1) < lb(1,1) || s(2) < lb(1,2) || s(1) > ub(1,1) || s(2) > ub(1,2)
                f_s = 1e6;
            else
                f_s = fun(s);
            end
            
            if f_s < f_r
               X(3,:) = s;
               F(3) = f_s;

            else
                X(3,:) = r;
                F(3) = f_r;
            end         
                      
        elseif f_r >= F(2)
            % CONTRACTION
            do_shrink = false; % Flag to tell us if we need to shrink
            
            if f_r < F(3)
                % Contract OUTSIDE
                c = m + (r - m)/2;
                if c(1) < lb(1) || c(2) < lb(2) || c(1) > ub(1) || c(2) > ub(2)
                    f_c = 1e6;
                else
                    f_c = fun(c);
                end
                
                if f_c < f_r
                    X(3,:) = c;
                    F(3) = f_c;
                else
                    do_shrink = true;
                end
            else
                % Contract INSIDE
                cc = m + (X(3,:) - m)/2;
                if cc(1) < lb(1) || cc(2) < lb(2) || cc(1) > ub(1) || cc(2) > ub(2)
                    f_cc = 1e6;
                else
                    f_cc = fun(cc);
                end
                
                if f_cc < F(3)
                    X(3,:) = cc;
                    F(3) = f_cc;
                else
                    do_shrink = true;
                end
            end
            
            % SHRINK
            if do_shrink
                % Pull point 2 halfway to point 1
                X(2,:) = X(1,:) + (X(2,:) - X(1,:))/2;
                if X(2,1) < lb(1) || X(2,2) < lb(2) || X(2,1) > ub(1) || X(2,2) > ub(2)
                    F(2) = 1e6;
                else
                    F(2) = fun(X(2,:));
                end
                
                % Pull point 3 halfway to point 1
                X(3,:) = X(1,:) + (X(3,:) - X(1,:))/2;
                if X(3,1) < lb(1) || X(3,2) < lb(2) || X(3,1) > ub(1) || X(3,2) > ub(2)
                    F(3) = 1e6;
                else
                    F(3) = fun(X(3,:));
                end
            end
        end 

    end

    c_root_best = X(1,1);
    c_tip_best = X(1,2);
end
