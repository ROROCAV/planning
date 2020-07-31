function [Aieq, bieq] = getAbieq(n_seg, n_order, corridor_range, ts, v_max, a_max)
    n_all_poly = n_seg*(n_order+1);
    n_one_poly = n_order + 1;
    %#####################################################
    % STEP 3.2.1 p constraint
    % 2个不等式* n_seg段 * 系数
    Aieq_p = zeros(2*n_all_poly,n_all_poly);
    bieq_p = zeros(2*n_all_poly,1);

    for i = 1:n_all_poly
        Aieq_p(i,i) = 1;
        index = idivide(uint32(i),uint32(n_one_poly),'ceil');
        bieq_p(i,1) = corridor_range(1,index);

        Aieq_p(i+n_all_poly,i) = -1;
        index = idivide(uint32(i),uint32(n_one_poly),'ceil');
        bieq_p(i+n_all_poly,1) = -corridor_range(2,index);
    end
    % corridor_range
    % Aieq_p
    % bieq_p
    %#####################################################
    % STEP 3.2.2 v constraint   
    Aieq_v = zeros(2*n_all_poly,n_all_poly);
    bieq_v = zeros(2*n_all_poly,1);

    for i = 2:n_all_poly
        Aieq_v(i,i) = 1 * n_order;
        Aieq_v(i,i-1) = -1 * n_order;
        bieq_v(i,1) = v_max;

        Aieq_v(i,i) = -1 * n_order;
        Aieq_v(i,i-1) = 1 * n_order;
        bieq_v(i,1) = -v_max;
    end
    %#####################################################
    % STEP 3.2.3 a constraint   
    Aieq_a = zeros(2*n_all_poly,n_all_poly);
    bieq_a = zeros(n_all_poly,1);

    for i = 3:n_all_poly
        sj = 1; %only for this example
        pref = (n_order^2 - n_order)/sj;

        Aieq_a(i,i) = 1 * pref;
        Aieq_v(i,i-1) = -2 * pref;
        Aieq_v(i,i-1) = 1 *pref;
        bieq_v(i,1) = a_max;

        Aieq_a(i,i) = -1 * pref;
        Aieq_v(i,i-1) = 2 * pref;
        Aieq_v(i,i-1) = -1 *pref;
        bieq_v(i,1) = -a_max;


    end
    
    %#####################################################
    % combine all components to form Aieq and bieq   
    Aieq = [Aieq_p; Aieq_v; Aieq_a];
    bieq = [bieq_p; bieq_v; bieq_a];
    Aieq = Aieq_p;
    bieq = bieq_p;
end