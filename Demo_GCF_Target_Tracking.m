%% Generalized conversion based filter (GCF) for target tracking using polar measurements
% reference:
% "Generalized Conversion Based Nonlinear Filtering Using Deterministic Sampling for Target Tracking" 
% (submitted to IEEE Transactions on Aerospace and Electronic Systems) by Jian Lan, 2023.
% Copyright by Jian Lan, Xi'an Jiaotong University
% Email:  lanjian@mail.xjtu.edu.cn

clear;
clc;
warning on;
%% Initialization
runs = 10; % number of Monte Carlo runs

% initial state and covariance
InitPos_x = -20000; % m
InitVel_x = 300; % m/s
InitPos_y = 15000; % m
InitVel_y = -10; % m/s
X = [InitPos_x;
      InitVel_x;
      InitPos_y;
      InitVel_y];
P= diag([10000,5000,10000,5000]);

T=1;
F = [1 T 0 0;
     0 1 0 0;
     0 0 1 T;
     0 0 0 1];
Gama=[0.5*(T^2) 0;
              T 0;
      0 0.5*(T^2);
              0 T];     
w_mu = [0;0];
Q = diag([1,1]).*0.0001;

% mean and covariance of measurement noise
v_mu = [0;0];
r_sigma = 200;
b_sigma = 0.2;            
R = diag([r_sigma^2,b_sigma^2]);

startTime = 0;
endTime = 100;

Time = endTime - startTime;
Step = Time/T;

%% Generation of GHQ weights and samples.
N_GHQ1 = 3;
[GHQ1.x, GHQ1.w] = hermquad(N_GHQ1);

%% Generate matrix "My" a priori
ind = 1;
cons_GHQ = 1/(pi^3);
for tmp1 = 1:N_GHQ1
    for tmp2 = 1:N_GHQ1
        for tmp3 = 1:N_GHQ1
            for tmp4 = 1:N_GHQ1
                for tmp5 = 1:N_GHQ1
                    for tmp6 = 1:N_GHQ1
                        weight_GHQ = GHQ1.w(tmp1, 1)*GHQ1.w(tmp2, 1)*GHQ1.w(tmp3, 1)*GHQ1.w(tmp4, 1)*GHQ1.w(tmp5, 1)*GHQ1.w(tmp6, 1);
                        mu3(ind, 1) = weight_GHQ;
                        ind = ind + 1;
                    end
                end
            end
        end
    end
end
ind = ind - 1;
My = zeros(ind, ind);
mu3 = mu3*cons_GHQ;
for tmp1 = 1:ind
    A3 = -mu3;
    A3(tmp1, 1) = 1-mu3(tmp1, 1);
    
    My = My + mu3(tmp1, 1)*A3*A3';
end
%


'Get Started'

%% Monte Carlo runs
for j = 1:runs
         P0 = P;
         X0 = X;
         
         X = mvnrnd(X0, P0)';   


         X_GCF = X0;
         P_GCF = P0;
         Estimand_GCF = X_GCF;

         X_GCFc = X0;
         P_GCFc = P0;
         Estimand_GCFc = X_GCFc;
         
          for i=1:Step
            %% Measurement generation
            Noise_w = mvnrnd(w_mu,Q)';
            Noise_v = mvnrnd(v_mu,R)';                   
            X = F * X + Gama * Noise_w;
            Z = [sqrt(X(1)^2 + X(3)^2); atan2(X(3), X(1)) ]+Noise_v;       
            
            %% Filtering
            'GCF';
            n_order = 3; % the order of the derivatives constraint
            dist_gate1 = 1e-3;
            [X_GCF, P_GCF, n_c1] = GCF(My, X_GCF, P_GCF, Z, Q, w_mu, v_mu, F, Gama, R, GHQ1, n_order, dist_gate1);
            n_c_GCF(:, i, j) = n_c1;
            
           	'GCF using clustering';
            n_order = 3;
            dist_gate2 = 1e+2;
            [X_GCFc, P_GCFc, n_c2] = GCF(My, X_GCFc, P_GCFc, Z, Q, w_mu, v_mu, F, Gama, R, GHQ1, n_order, dist_gate2);
            n_c_GCFc(:, i, j) = n_c2;
            
            %% Evaluation using RMSE
            if j == 1
                RMSE_GCF(:,i,j) = [(X_GCF(1) - X(1, 1))^2+(X_GCF(3) - X(3, 1))^2; (X_GCF(2) - X(2, 1))^2+(X_GCF(4) - X(4, 1))^2];
                
                RMSE_GCFc(:,i,j) = [(X_GCFc(1) - X(1, 1))^2+(X_GCFc(3) - X(3, 1))^2; (X_GCFc(2) - X(2, 1))^2+(X_GCFc(4) - X(4, 1))^2];
            else
                RMSE_GCF(:,i,j) = RMSE_GCF(:,i,j-1) + ...
                    [(X_GCF(1) - X(1, 1))^2+(X_GCF(3) - X(3, 1))^2; (X_GCF(2) - X(2, 1))^2+(X_GCF(4) - X(4, 1))^2];
                
                RMSE_GCFc(:,i,j) = RMSE_GCFc(:,i,j-1) + ...
                    [(X_GCFc(1) - X(1, 1))^2+(X_GCFc(3) - X(3, 1))^2; (X_GCFc(2) - X(2, 1))^2+(X_GCFc(4) - X(4, 1))^2];
            end
            
          end %i
j
end %j
RMSE_GCF_PRT = (RMSE_GCF(:,:,runs)./runs).^0.5;
RMSE_GCFc_PRT = (RMSE_GCFc(:,:,runs)./runs).^0.5;

%% Calculating the numbers of clusters of GCF and GCFc
M1 = [1/2 1/2];
M2 = ones(Step, 1);
MEAN_n_c_GCF = M1*n_c_GCF(:, :, 1)*M2;
MEAN_n_c_GCFc = M1*n_c_GCFc(:, :, 1)*M2;

for j = 2:runs
    MEAN_n_c_GCF = MEAN_n_c_GCF + M1*n_c_GCF(:, :, j)*M2;
    MEAN_n_c_GCFc = MEAN_n_c_GCFc + M1*n_c_GCFc(:, :, j)*M2;
end
MEAN_n_c_GCF = MEAN_n_c_GCF/runs/Step;
MEAN_n_c_GCFc = MEAN_n_c_GCFc/runs/Step;


%% Showing results
MEAN_n_c_GCF
MEAN_n_c_GCFc

clear i;
i=1:Step;

figure;
plot(i,RMSE_GCF_PRT(1,i), 'b-d', 'MarkerSize', 1.5);
hold on;
plot(i,RMSE_GCFc_PRT(1,i), 'k-d', 'MarkerSize', 1.5);

xlabel('time (s)');ylabel('position RMSE (m)');
legend('GCF ({\it{n}}=3)', 'GCFc ({\it{n}}=3)');
hold off;


figure;
plot(i,RMSE_GCF_PRT(2,i), 'b-d', 'MarkerSize', 1.5);
hold on;
plot(i,RMSE_GCFc_PRT(2,i), 'k-d', 'MarkerSize', 1.5);

xlabel('time (s)');ylabel('velocity RMSE (m/s)');
legend('GCF ({\it{n}}=3)', 'GCFc ({\it{n}}=3)');
hold off;


