C = input('Enter C: ');
A = input('Enter A: ');
[m, n] = size(A);
b = input('Enter b: ');
[x_0, z_I] = phase_I(A, b);
if z_I > 0
    disp("The problem is infeasible!");
else
    [w_star, z_B, Dual_Benson_exitflag] = Dual_Benson(C, A, b, x_0);
    if Dual_Benson_exitflag>0
        if ~isequal(z_B,0)
            x_bar = L_P_w_star(w_star, C, A, b);
        else
            x_bar = x_0;
        end
        [ebs_index, ebd] = phase_III(x_bar, C, A, b);
        ebs = [];
        for J_B = ebs_index
            [x_bar, ~, ~, ~, ~] = interpreter(J_B, C, A, b);
            ebs = [ebs x_bar];
        end
        disp("The set of basic efficient solutions is:");
        disp(ebs);
        if isempty(ebd)
            disp("There are no efficient direction!")
        else
            for column = ebd
                disp("The set of basic efficient directions is:");
                disp(ebd(1));
                disp(ebd(2:end));
            end
        end
    else
        disp("The set of efficient solutions is empty!");
    end
end

function [x_0, z_I] = phase_I(A,b)
    [m, n] = size(A);
    e = [zeros(n,1);ones(m,1)];
    Aeq = [A eye(m)];
    beq = b;
    A = [];
    b = [];
    lb = zeros(1,n+m);
    ub = [];
    [x_0, z_I] = linprog(e, A, b, Aeq, beq, lb, ub);
    x_0 = x_0(1:n,1);
end

function [d_star, z_B, Dual_Benson_exitflag] = Dual_Benson(C, A, b, x_0)
    [p, n] = size(C);
    [m, ~] = size(A);
    c = [C*x_0; b].';
    Aeq = [];
    beq = [];
    A = [-C -eye(p); -A zeros(m, p)].';
    b = [zeros(n , 1); -ones(p,1)];
    lb = [];
    ub = [];
    [d_star, z_B, Dual_Benson_exitflag] = linprog(c, A, b, Aeq, beq, lb, ub);
    d_star = d_star(end-p+1 : end).';
end

function x_bar = L_P_w_star(w_star, C, A, b)
    [~, n] = size(A);
    c = w_star*C;
    Aeq = A;
    beq = b;
    A = [];
    b = [];
    lb = zeros(1,n);
    ub = [];
    [x_bar] = linprog(c, A, b, Aeq, beq, lb, ub);
end

function [ebs_index, ebd_index] = phase_III(x_bar, C, A, b)
    [m, n] = size(A);
    J_N = find(x_bar == 0, n-m, 'first');
    J_B = setdiff(1:numel(x_bar), J_N);
    index_check_list = [J_B.'];
    ebs_index = [];
    ebd_index = [];
    iter = 0;
    while ~isempty(index_check_list)
        iter = iter + 1;
        J_B = index_check_list(:, 1).';
        ebs_index = [ebs_index J_B.'];
        [~, B, b_hat, J_N, R] = interpreter(J_B, C, A, b);
        index_check_list = index_check_list(:,2:end);
        nbe = [];
        for j = J_N
            [~, z_r_j] = nbe_finder(j, R, A, J_N);
            if (isequal(z_r_j,0))
                nbe = [nbe j];
            end
        end
        Y = inv(B)*A;
        for j = nbe
            [new_J_B, ~] = pivot(Y( : , j), b_hat, J_B, J_N, j);
            dd = direction_detector(new_J_B, C, A, b);
            if isequal(dd, 1)
                d = direction_creator(J_B, J_N, m, n, j);
                ebd_index = [ebd_index [iter; d]];
            end
            check = ebs_inclusion_check(ebs_index, index_check_list, new_J_B, C, A, b);
            if check == 0
                index_check_list = [index_check_list new_J_B.'];
            end
        end
    end
end

function d = direction_creator(J_B, J_N, m, n, j)
    d = zeros(n,1);
    I = eye(n-m);
    d(J_B) = -Y( : , j);
    d(J_N) = I(:,j);
end

function [new_J_B, new_J_N] = pivot(y_j, b_hat, J_B, J_N, j)
    K = b_hat ./ y_j;
    positive_indices = find(K > 0);
    if ~isempty(positive_indices)
        min_positive_value = min(K(positive_indices));
    end
    r = find( K == min_positive_value, 1, 'first');
    k = find(J_N == j);
    new_J_B = J_B;
    new_J_N = J_N;
    new_J_B(r) = j;
    new_J_N(k) = J_B(r);
end

function check = ebs_inclusion_check(ebs_index, index_check_list, new_J_B, C, A, b)  
    check = 0;
    [x_new, ~, ~, ~, ~] = interpreter(new_J_B, C, A, b);
    for J_B = ebs_index
        [x_bar, ~, ~, ~, ~] = interpreter(J_B, C, A, b);
        if x_bar == x_new
            check = 1;
            return
        end
    end
    for J_B = index_check_list
        [x_bar, ~, ~, ~, ~] = interpreter(J_B, C, A, b);
        if x_bar == x_new
            check = 1;
            return
        end
    end
end

function dd = direction_detector(new_J_B, C, A, b)
    dd = 0;
    [~, new_B, ~, ~, ~] = interpreter(new_J_B, C, A, b);
    new_Y = inv(new_B)*A;
    for y = new_Y
        if y<=0
            dd = 1;
            return
        end
    end
end

function [x_bar, B, b_hat, J_N, R] = interpreter(J_B, C, A, b)
    [~, n] = size(A); 
    B = A(:, J_B);
    b_hat = inv(B)*b;
    x_bar = zeros(n,1);
    x_bar(J_B) = b_hat;
    J_N = (setdiff(1:numel(x_bar), J_B));
    R = C(:, J_N) - C(:, J_B)*inv(B)*A(:, J_N);
end

function [x_r_j, z_r_j] = nbe_finder(j, R, A, J_N)
    [m, n] = size(A);
    k = find(J_N == j);
    r_j = R(:,k);
    [p,~] = size(R);
    c = [zeros(n-m+1,1); -ones(p,1)];
    Aeq = [R -r_j eye(p)];
    beq = zeros(p,1);
    A = [];
    b = [];
    lb = zeros(1, n-m+1+p);
    ub = [];
    [x_r_j, z_r_j] = linprog(c, A, b, Aeq, beq, lb, ub);
end