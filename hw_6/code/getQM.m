function [Q, M] = getQM(n_seg, n_order, ts)
    Q = [];
    M = [];
    M_k = getM(n_order);
    n_one_poly = n_order + 1;
    for k = 1:n_seg
        %#####################################################
        % STEP 2.1 calculate Q_k of the k-th segment 
        Q_k = zeros(n_one_poly,n_one_poly);
        % disp(Q_k(2,2));
        % disp(ts(k));
        % fprintf('K %d \n', k);
        der = 3;
        for j = 1:n_one_poly
            paper_i = j - 1;
            for h = 1:n_one_poly
                paper_l = h - 1;
                l_factor = 0;
                i_factor = 0;
                if paper_l >= der && paper_i >= der
                    i_factor = factorial(paper_i) / factorial(paper_i - der);
                    l_factor = factorial(paper_l) / factorial(paper_l - der);
                    root = (paper_i + paper_l - (der * 2 - 1));
                    Q_k(j,h) = i_factor * l_factor / root * ts(k)^(root);    
                    % fprintf('Root number: %d, paperi: %d, paperL:%d, i_factor:%d, l_factor:%d \n',root, paper_i, paper_l, i_factor, l_factor);
                end
            end
        end

        Q = blkdiag(Q, Q_k);
        M = blkdiag(M, M_k);
    end
end