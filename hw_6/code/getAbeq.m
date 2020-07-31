function [Aeq, beq] = getAbeq(n_seg, n_order, ts, start_cond, end_cond)
    n_all_poly = n_seg*(n_order+1);
    n_one_poly = n_order + 1;
    M = [];
    M_k = getM(n_order);
    for k = 1:n_seg
        M = blkdiag(M, M_k);
    end
    %#####################################################
    % STEP 2.1 p,v,a constraint in start 
    Aeq_start = zeros(3, n_all_poly);
    beq_start = zeros(3, 1);
    for j = 1:3
        der = j - 1;
        for i = 1:n_one_poly
            paper_i = i - 1;
            if paper_i - der >= 0
                Aeq_start(j, i) = factorial(paper_i) / factorial(paper_i - der) * 0^(paper_i - der);
            else
                Aeq_start(j, i) = 0;
            end
        end
        beq_start(j, 1) = start_cond(j);
    end
    % fprintf('Aeq_start \n');
    % disp(Aeq_start)
    % fprintf('beq_start \n');
    % disp(beq_start)
    %#####################################################
    % STEP 2.2 p,v,a constraint in end
    Aeq_end = zeros(3, n_all_poly);
    beq_end = zeros(3, 1);
    for j = 1:3
        der = j - 1;
        for i = 1:n_one_poly
            paper_i = i - 1;
            if paper_i - der >= 0
                Aeq_end(j,n_all_poly-(n_one_poly-i)) = factorial(paper_i) / factorial(paper_i - der) * ts(n_seg)^(paper_i - der);
            else
                Aeq_end(j,n_all_poly-(n_one_poly-i)) = 0;
            end
        end
        beq_end(j,1) = end_cond(j);
    end
    % fprintf('Aeq_end \n');
    % disp(Aeq_end)
    % fprintf('beq_end \n');
    % disp(beq_end)
    
    %#####################################################
    % STEP 2.3 position continuity constrain between 2 segments
    Aeq_con_p = zeros(n_seg-1, n_all_poly);
    beq_con_p = zeros(n_seg-1, 1);
    % disp(Aeq_con_p);
    for j = 1:n_seg-1
        der = 0;
        for i = 1:n_one_poly
            paper_i = i - 1;
            if paper_i - der >= 0
                Aeq_con_p(j, i + (j-1) * n_one_poly) = factorial(paper_i) / factorial(paper_i - der) * ts(j)^(paper_i - der);
            else
                Aeq_con_p(j, i + (j-1) * n_one_poly) = 0;
            end
        end  

        for l = 1:n_one_poly
            paper_l = l - 1;
            if paper_l - der >= 0
                Aeq_con_p(j, l + j * n_one_poly) = factorial(paper_l) / factorial(paper_l - der) * 0^(paper_l - der) * -1;
            else
                Aeq_con_p(j, l + j * n_one_poly) = 0;
            end
        end   
    end
    % fprintf('Aeq_con_p \n');
    % disp(Aeq_con_p)
    % fprintf('beq_con_p \n');
    % disp(beq_con_p)

    %#####################################################
    % STEP 2.4 velocity continuity constrain between 2 segments
    Aeq_con_v = zeros(n_seg-1, n_all_poly);
    beq_con_v = zeros(n_seg-1, 1);
    for j = 1:n_seg-1
        der = 1;
        for i = 1:n_one_poly
            paper_i = i - 1;
            if paper_i - der >= 0
                Aeq_con_v(j, i + (j-1) * n_one_poly) = factorial(paper_i) / factorial(paper_i - der) * ts(j)^(paper_i - der);
            else
                Aeq_con_v(j, i + (j-1) * n_one_poly) = 0;
            end
        end  

        for l = 1:n_one_poly
            paper_l = l - 1;
            if paper_l - der >= 0
                Aeq_con_v(j, l + j * n_one_poly) = factorial(paper_l) / factorial(paper_l - der) * 0^(paper_l - der) * -1;
            else
                Aeq_con_v(j, l + j * n_one_poly) = 0;
            end
        end   
    end
    % fprintf('Aeq_con_v \n');
    % disp(Aeq_con_v)
    % fprintf('beq_con_v \n');
    % disp(beq_con_v)

    %#####################################################
    % STEP 2.5 acceleration continuity constrain between 2 segments
    Aeq_con_a = zeros(n_seg-1, n_all_poly);
    beq_con_a = zeros(n_seg-1, 1);

    for j = 1:n_seg-1
        der = 2;
        for i = 1:n_one_poly
            paper_i = i - 1;
            if paper_i - der >= 0
                Aeq_con_a(j, i + (j-1) * n_one_poly) = factorial(paper_i) / factorial(paper_i - der) * ts(j)^(paper_i - der);
            else
                Aeq_con_a(j, i + (j-1) * n_one_poly) = 0;
            end
        end  

        for l = 1:n_one_poly
            paper_l = l - 1;
            if paper_l - der >= 0
                Aeq_con_a(j, l + j * n_one_poly) = factorial(paper_l) / factorial(paper_l - der) * 0^(paper_l - der) * -1;
            else
                Aeq_con_a(j, l + j * n_one_poly) = 0;
            end
        end   
    end

    %#####################################################
    % combine all components to form Aeq and beq   
    Aeq_con = [Aeq_con_p; Aeq_con_v; Aeq_con_a];
    beq_con = [beq_con_p; beq_con_v; beq_con_a];
    Aeq = [Aeq_start; Aeq_end; Aeq_con];
    Aeq = Aeq * M;
    beq = [beq_start; beq_end; beq_con];
end