function M = getM(n_seg, n_order, ts)
    M = [];
    n_one_poly = n_order + 1;

    for k = 1:n_seg
        M_k = zeros(8,n_one_poly);

        for der = 0:3
            % calculate  p v a j at t_start
            for i = 1:n_one_poly
                paper_i = i - 1;
                if paper_i - der >= 0
                    M_k(der + 1, i) = factorial(paper_i) / factorial(paper_i - der) * 0^(paper_i - der);
                else
                    M_k(der + 1, i) = 0;
                end
            end  
            % calculate  p v a j at t_end
            for i = 1:n_one_poly
                paper_i = i - 1;
                if paper_i - der >= 0
                    M_k(der + 1 + 4, i) = factorial(paper_i) / factorial(paper_i - der) * ts(k)^(paper_i - der);
                else
                    M_k(der + 1 + 4, i) = 0;
                end
            end  
        end
        %#####################################################
        % STEP 1.1: calculate M_k of the k-th segment 
        %
        %
        %
        %
        M = blkdiag(M, M_k);
    end
end