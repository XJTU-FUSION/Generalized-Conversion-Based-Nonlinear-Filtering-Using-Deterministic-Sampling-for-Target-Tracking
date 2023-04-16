function [xkk, Pkk, n_c] = GCF(My, xk_1k_1,Pk_1k_1,zk,...
    Qk_1,wmk_1,vmk_1,Fk_1,Gk_1,Rk, GHQ1, n_cons, dist_gate)
%% Generalized conversion based filter (GCF) for target tracking using polar measurements
% reference:
% "Generalized Conversion Based Nonlinear Filtering Using Deterministic Sampling for Target Tracking" 
% (submitted to IEEE Transactions on Aerospace and Electronic Systems) by Jian Lan, 2023.
% Copyright by Jian Lan, Xi'an Jiaotong University
% Email:lanjian@mail.xjtu.edu.cn

[N_GHQ1, tmp] = size(GHQ1.x);
%% Measurement transformed to Cartesian coordinates
z = zk;
zk = [z(1)*cos(z(2)); z(1)*sin(z(2))];

%% Prediction
xkk_1 = Fk_1*xk_1k_1 + Gk_1*wmk_1;
Pkk_1 = Fk_1*Pk_1k_1*Fk_1' + Gk_1*Qk_1*Gk_1';

%% GHQ sample generation
C_GHQ = blkdiag(Pkk_1*1, Rk); % covariances of the predicted state and measurement
C_GHQ_root = (chol(2.*C_GHQ))';

cons_GHQ = 1/(pi^3);
X_GHQ = [xkk_1; vmk_1]; % means of the predicted state and measurement
[N_X, tmp] = size(X_GHQ);
s_aug = zeros(N_X, 1);

ind = 1;
clear z_tmp;
for tmp1 = 1:N_GHQ1
    s_aug(1, 1) = GHQ1.x(tmp1, 1);
    for tmp2 = 1:N_GHQ1
        s_aug(2, 1) = GHQ1.x(tmp2, 1);
        for tmp3 = 1:N_GHQ1
            s_aug(3, 1) = GHQ1.x(tmp3, 1);
            for tmp4 = 1:N_GHQ1
                s_aug(4, 1) = GHQ1.x(tmp4, 1);
                for tmp5 = 1:N_GHQ1
                    s_aug(5, 1) = GHQ1.x(tmp5, 1);
                    for tmp6 = 1:N_GHQ1
                        s_aug(6, 1) = GHQ1.x(tmp6, 1);
                        
                        weight_GHQ = GHQ1.w(tmp1, 1)*GHQ1.w(tmp2, 1)*GHQ1.w(tmp3, 1)*GHQ1.w(tmp4, 1)*GHQ1.w(tmp5, 1)*GHQ1.w(tmp6, 1);
                        x_aug = C_GHQ_root*s_aug + X_GHQ;
                        z_tmp = ...
                            [(sqrt(x_aug(1, 1)^2+x_aug(3, 1)^2) + x_aug(5, 1))*cos(atan2(x_aug(3, 1), x_aug(1, 1)) + x_aug(6, 1));
                             (sqrt(x_aug(1, 1)^2+x_aug(3, 1)^2) + x_aug(5, 1))*sin(atan2(x_aug(3, 1), x_aug(1, 1)) + x_aug(6, 1));];
                         
                        mu(ind, 1) = weight_GHQ;
                        X_mu(:, ind) = x_aug;
                        Z_mu(:, ind) = z_tmp;
                        
                        ind = ind + 1;
                    end
                end
            end
        end
    end
end
ind = ind - 1;

mu = cons_GHQ*mu;

X_mu = (X_mu - repmat(X_GHQ, 1, ind))*diag(mu);

UP_xz = X_mu(1:4, :);
Down_z = My;

index_tmp = ind;
Dev_Matrix = zeros(index_tmp, index_tmp+1);
Sum_Matrix = zeros(index_tmp, index_tmp+1);
for tmp1 = 1:(index_tmp)
    Dev_Matrix(tmp1, tmp1) = -1;
    Dev_Matrix(tmp1, tmp1+1) = 1;
    
    Sum_Matrix(tmp1, tmp1) = 1/2;
    Sum_Matrix(tmp1, tmp1+1) = 1/2;
end

%% Generating the constraint matrix and calculating Q at current time
u1 = eye(2, 2);
for tmp1 = 1:2
    vec_proj = u1(tmp1, :);
    [Qt_tmp{tmp1} n_c(tmp1)] = Q_generate(vec_proj*Z_mu, vec_proj*zk, n_cons, Dev_Matrix, Sum_Matrix, dist_gate);
end
Qt = [Qt_tmp{1}; Qt_tmp{2}];

%% Update
As = UP_xz*Qt';
Bs = Qt*Down_z*Qt';

xkk = xkk_1 - As*inv(Bs)*Qt*mu;
Pkk = Pkk_1 - As*inv(Bs)*As';


    

