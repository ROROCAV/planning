function Ct = getCt(n_seg, n_order)
    n_all_con = 4*(n_seg+1);
    dp_start = n_seg + 7;
    %#####################################################
    % STEP 2.1: finish the expression of Ct
    %
    %
    %
    %
    %
    
    Ct_c = zeros(n_all_con,n_all_con);
    Ct = zeros(8*n_seg,n_all_con);
    for i = 1:n_all_con
        % pvaj begin
        if i <=4
            Ct_c(i,i) = 1;
            continue;
        end
        % pvaj end
        if i > n_all_con - 4
            Ct_c(i, dp_start - (n_all_con - i)) = 1;
            continue;
        end
        % is p
        if i > 4 && mod(i,4) == 1
            Ct_c(i, 4 + idivide(int32(i),int32(4),'floor')) = 1;
            continue;
        end
        % is v
        if i > 4 && mod(i,4) == 2
            Ct_c(i,dp_start + (idivide(int32(i),int32(4),'floor')-1)*3+1) = 1;
            continue;
        end
        % a
        if i > 4  && mod(i,4) == 3
            Ct_c(i,dp_start + (idivide(int32(i),int32(4),'floor')-1)*3+2) = 1;
            continue;
        end
        % j
        if i > 4 && mod(i,4) == 0
            Ct_c(i,dp_start + (idivide(int32(i-1),int32(4),'floor')-1)*3+3) = 1;
            continue;
        end
    end

    % append middle waypoint constrain
    for i = 1:n_all_con
        if i == 4 
            Ct(1:4,:) = Ct_c(1:4,:);
        end
        if i == n_all_con - 3
            Ct(8*n_seg-3:8*n_seg,:) = Ct_c(i:n_all_con,:);
        end
        if i > 4 && i < n_all_con - 2 && mod(i,4) == 0
            index = ((i-4)/4-1)*8+4;
            
            Ct(index+1:index+4,:) = Ct_c(i-3:i,:);
            Ct(index+5:index+8,:) = Ct_c(i-3:i,:);
        end
    end

end