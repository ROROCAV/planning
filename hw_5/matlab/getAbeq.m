function [Aeq beq]= getAbeq(n_seg, n_order, waypoints, ts, start_cond, end_cond)
    n_one_poly = n_order + 1;
    n_all_poly = n_seg*(n_one_poly);
    %#####################################################
    % p,v,a,j constraint in start, 此时约束的应该是第一段poly轨迹的p参数
    Aeq_start = zeros(4, n_all_poly);
    beq_start = zeros(4, 1);
    for j = 1:4
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
    fprintf('Aeq_start \n');
    disp(Aeq_start)
    fprintf('beq_start \n');
    disp(beq_start)
    % STEP 2.1: write expression of Aeq_start and beq_start
    %
    %
    %
    %
    
    %#####################################################
    % p,v,a j constraint in end 此时约束的应该是最后一段poly轨迹的p参数
    Aeq_end = zeros(4, n_all_poly);
    beq_end = zeros(4, 1);
    for j = 1:4
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
    fprintf('Aeq_end \n');
    disp(Aeq_end)
    fprintf('beq_end \n');
    disp(beq_end)
    % STEP 2.2: write expression of Aeq_end and beq_end
    %
    %
    %
    %
    
    % %#####################################################
    % % position constrain in all middle waypoints
    Aeq_wp = zeros(n_seg-1, n_all_poly);
    beq_wp = zeros(n_seg-1, 1);
    % STEP 2.3: write expression of Aeq_wp and beq_wp
    %
    %
    %
    %
    for i=1:n_seg-1
        start_idx_2 = n_one_poly * i;
        Aeq_wp(i,start_idx_2 + 1) = 1;
        beq_wp(i,1) = waypoints(i+1);
    end
    %#####################################################
    % position continuity constrain between each 2 segments
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
    % STEP 2.4: write expression of Aeq_con_p and beq_con_p
    %
    %
    %
    %
    
    %#####################################################
    % velocity continuity constrain between each 2 segments
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
    fprintf('Aeq_con_v \n');
    disp(Aeq_con_v)
    fprintf('beq_con_v \n');
    disp(beq_con_v)
    % STEP 2.5: write expression of Aeq_con_v and beq_con_v
    %
    %
    %
    %

    %#####################################################
    % acceleration continuity constrain between each 2 segments
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

    % fprintf('Aeq_con_a \n');
    % disp(Aeq_con_a)
    % fprintf('beq_con_a \n');
    % disp(beq_con_a)

    % STEP 2.6: write expression of Aeq_con_a and beq_con_a
    %
    %
    %
    %
    
    %#####################################################
    % jerk continuity constrain between each 2 segments
    Aeq_con_j = zeros(n_seg-1, n_all_poly);
    beq_con_j = zeros(n_seg-1, 1);
    for j = 1:n_seg-1
        der = 3;
        for i = 1:n_one_poly
            paper_i = i - 1;
            if paper_i - der >= 0
                Aeq_con_j(j, i + (j-1) * n_one_poly) = factorial(paper_i) / factorial(paper_i - der) * ts(j)^(paper_i - der);
            else
                Aeq_con_j(j, i + (j-1) * n_one_poly) = 0;
            end
        end  

        for l = 1:n_one_poly
            paper_l = l - 1;
            if paper_l - der >= 0
                Aeq_con_j(j, l + j * n_one_poly) = factorial(paper_l) / factorial(paper_l - der) * 0^(paper_l - der) * -1;
            else
                Aeq_con_j(j, l + j * n_one_poly) = 0;
            end
        end   
    end

    % fprintf('Aeq_con_j \n');
    % disp(Aeq_con_j)
    % fprintf('beq_con_j \n');
    % disp(beq_con_j)

    % STEP 2.7: write expression of Aeq_con_j and beq_con_j
    %
    %
    %
    %
    
    %#####################################################
    % combine all components to form Aeq and beq   
    Aeq_con = [Aeq_con_p; Aeq_con_v; Aeq_con_a; Aeq_con_j];
    beq_con = [beq_con_p; beq_con_v; beq_con_a; beq_con_j];
    Aeq = [Aeq_start; Aeq_end; Aeq_wp; Aeq_con];
    beq = [beq_start; beq_end; beq_wp; beq_con];
end