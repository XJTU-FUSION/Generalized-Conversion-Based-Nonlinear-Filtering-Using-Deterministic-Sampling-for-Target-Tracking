function [Qt n_c] = Q_generate(Z_gen, zk, n_cons, Dev_Matrix, Sum_Matrix, dist_gate)
%% Generating the constraint matrix and calculating Q at current time in generalized conversion based filter (GCF)
% reference:
% "Generalized Conversion Based Nonlinear Filtering Using Deterministic Sampling for Target Tracking" 
% (submitted to IEEE Transactions on Aerospace and Electronic Systems) by Jian Lan, 2023.
% Copyright by Jian Lan, Xi'an Jiaotong University
% Email:lanjian@mail.xjtu.edu.cn

    [n_z, ind] = size(Z_gen);
    
    Z_gen_normal1_1 = Z_gen;
    zk_std_1 = zk;
    Z_dev1 = ([Z_gen, zk]);
    
    C_center0(1, 1) = Z_gen_normal1_1(1, 1);
    count_center0 = 1;
    for tmp3 = 2:ind
        flag = 1;
        for tmp2 = 1:count_center0
            dist_tmp = abs(Z_gen_normal1_1(1, tmp3) - C_center0(1, tmp2));
            if dist_tmp < dist_gate
                'test 1';
                flag = 0;
            end
        end
        if flag == 1
            count_center0 = count_center0 + 1;
            C_center0(1, count_center0) = Z_gen_normal1_1(1, tmp3);
        end
    end

    C_center0_save(1, 1:count_center0) = C_center0;
    count_center0_save(1) = count_center0;

Sel = 1;
zk_std = zk_std_1(1, :);
Z_dev = Z_dev1(1, :);
    
Z_dev_size = size(Z_dev);


C_center = C_center0_save(Sel, 1:count_center0_save(Sel));
count_center = count_center0_save(Sel);

count_center = count_center + 1;
C_center(1, count_center) = Z_dev(1, Z_dev_size(2));



H_tmp = zeros(count_center, Z_dev_size(2));
for tmp1 = 1:Z_dev_size(2)
    dist_cal = abs(C_center - Z_dev(1, tmp1));
    [dev_sort, sort_index] = sort(dist_cal, 'ascend');
    H_tmp(sort_index(1), tmp1) = 1;
end

count_com = count_center;
H_com_final1 = H_tmp;
H_com_est1 = H_com_final1(:, 1:(Z_dev_size(2))); 

Z_dev_new = C_center; %(diag(1./H_com_num)*H_com_final*X_std_dev')';
index_tmp = count_com-1;

[z_std_proj, proj_index2] = sort(Z_dev_new, 'ascend');


Sort_switch2 = zeros(index_tmp+1, index_tmp+1); % Switch the order of the items for sorting
for tmp1 = 1:(index_tmp+1)
    Sort_switch2(tmp1, proj_index2(tmp1)) = 1;
end



dev_gate = 1e-30; 0; 
for tmp1 = 1:n_cons 
	if tmp1 == 1
        Dev = Dev_Matrix(1:(index_tmp), 1:(index_tmp+1));
        Sum_tmp = z_std_proj';
        Z_down1 = Dev*Sum_tmp;
        
        Z_down = abs(Z_down1(:, 1)).^1;
        
      	for tmp2 = 1:(index_tmp - tmp1 + 1)
            if Z_down(tmp2, 1) < dev_gate
                'test';
                Z_down(tmp2, 1) = dev_gate;
            end
        end
        
        Cons_aug = diag(1./Z_down)*Dev*Sort_switch2;
    else
        Dev = Dev_Matrix(1:(index_tmp - tmp1 + 1), 1:(index_tmp - tmp1 + 2));
        Sum_tmp = Sum_Matrix(1:(index_tmp - tmp1 + 2), 1:(index_tmp - tmp1 + 3))*Sum_tmp;
        Z_down1 = Dev*Sum_tmp;
               
        Z_down = abs(Z_down1(:, 1)).^1;
        
      	for tmp2 = 1:(index_tmp - tmp1 + 1)
            if Z_down(tmp2, 1) < dev_gate
                'test';
                Z_down(tmp2, 1) = dev_gate;
            end
        end
        
        Cons_aug = diag(1./Z_down)*Dev*Cons_aug;
    end
    if tmp1 == 2
        Cons_aug_save = Cons_aug;
    end
end
Cons_aug1 = Cons_aug;


Cons_full = Cons_aug1;
H_com_est = H_com_est1;
H_z = H_com_est'; %dimensions: (ind+1)*n

Cons_full_STD = Cons_full(:, 1:(end - 1));

C_final = Cons_full_STD;
Cons_size2 = size(C_final);
[Q_c, R_c] = qr(C_final');
Q2 = Q_c(:, (Cons_size2(1)+1):Cons_size2(2));


n_reduced = 0;
Q2 = Q2(:, 1:(end - n_reduced));

n_c = index_tmp;

Qt = (H_z(1:(end-1), 1:(end - 1))*Q2)';


