function [dp_opt,dt_opt,dcca_opt,misfit] = cca_cv(temp,prox,cca_options)
% Choose an optimal set of CCA truncation parameters.
%
% Three truncation parameters are required for CCA reconstruction:
%     dp:   Truncation parameter of the proxy matrix
%     dt:   Truncation parameter of the temperature matrix
%     dcca: Truncation parameter of the covariance matrix of proxy and temperature
%
%
%     dcca by design is smaller than min(dp,dt). The optimal selection of
%     dt, dp and dcca is based on the k-fold cross-validation RMSE.
%
%     The CCA reconstruction method is inspired by Smerdon et al. [2010a]
%
%
% Input:
%     temp: Temperature matrix
%     prox: Proxy (pseudoproxy) matrix
%           Notice that temperature and proxy matrix should have the same
%           time dimension
%        k: Number of folds of K-fold cross-validation
%     indices: Indices of K-fold cross-validation
%     noise:  Index of noise level, used to save RMSE output
%
% Output:
%     dp_opt:   Optimal choice of dp
%     dt_opt:   Optimal choice of dt
%     dcca_opt: Optimal chocie of dcca
%     misfit:   A struct containing the useful error information as follows:
%     rmse:  RMSE of different folds
%     error: Truncated version of rmse
%
% See also TRUNCATION_CRITERIA, CCA_BP
fopts    = fieldnames(cca_options);
K        = cca_options.K;
indices  = cca_options.indices;
weights  = cca_options.weights;
% noise    = cca_options.noise;
dt_max   = cca_options.dt_max;
dp_max   = cca_options.dp_max;
method   = cca_options.method;

if sum(strcmpi('dp_range',fopts))>0
    dp_range   = cca_options.dp_range;
else
    dp_range   = 1:dp_max;
end

if sum(strcmpi('dt_range',fopts))>0
    dt_range   = cca_options.dt_range;
else
    dt_range   = 1:dt_max;
end

if strcmp('smerdon10',method)
    rmse = nan(min(dp_max,dt_max),dp_max,dt_max,K);
    for k = 1:K
        display(['Fold-',num2str(k)])
        % Define training and test data
        test = (indices == k); train = ~test;
        ntest = length(test);
        
        % Temperature matrix
        temp_train    = temp(train,:);
        temp_test     = temp(test,:);
        [Tds,Tm,Ts]   = standardize(temp_train);
        
        %  PUT OPTION HERE
        % dt_max(k)     = rank(temp_train);
        
        % Proxy matrix
        prox_train    = prox(train,:);
        [Pds,Pm,Ps]   = standardize(prox_train);
        
        [UP,SP,VP]    = svd(Pds');
        [UT,ST,VT]    = svd(Tds');
        
        % dp_max(k)     = rank(Pds);
        
        
        % Cross Covariance matrix
        for dt = dt_range
            disp(['----dt = ',num2str(dt), ', dp_max = ', num2str(max(dp_range)),'----'])
            for dp = dp_range
                % disp(['dp = ',num2str(dp)])
                dcca_max = min(dt,dp);
                for dcca = 1:dcca_max
                    F          = VT(:,1:dt)'*VP(:,1:dp);
                    [UF,SF,VF] = svd(F);
                    % Estimate CV regression matrix
                    B_cv       = UT(:,1:dt)*ST(1:dt,1:dt)*UF(1:dt,1:dcca)...
                        *SF(1:dcca,1:dcca)*VF(1:dp,1:dcca)'*...
                        diag(1.0./diag(SP(1:dp,1:dp)))*UP(:,1:dp)';
                    
                    % Predict temperature over verification sample
                    temp_pred  = (prox-ones(ntest,1)*Pm)*diag(1.0./Ps)*B_cv'*diag(Ts)+ ones(ntest,1)*Tm;
                    
                    % Calculate the misfit (RMSE)
                    mse                = mean((temp_pred(test,:) - temp_test).^2);
                    rmse(dcca,dp,dt,k) = sqrt(mse*weights);
                end
            end
        end
    end
    mean_misfit = sqrt(nmean(rmse.^2,4));
    min_RMSE = min(mean_misfit(:));
    [dcca_opt,dp_opt,dt_opt] = ind2sub(size(mean_misfit),find(mean_misfit == min_RMSE));
    
elseif strcmp('ne08',method)
    rmse  = nan(50,K);
    for k = 1:K
        display(['Fold-',num2str(k)])
        % Define training and test data
        test  = (indices == k); train = ~test;
        ntest = length(test);
        
        % Temperature matrix
        temp_train    = temp(train,:);
        temp_test     = temp(test,:);
        [Tds,Tm,Ts]   = standardize(temp_train);
        
        
        % Proxy matrix
        prox_train    = prox(train,:);
        [Pds,Pm,Ps]   = standardize(prox_train);
        
        [UP,SP,VP]    = svd(Pds');
        [UT,ST,VT]    = svd(Tds');
        
        [~,trunc_t]   = truncation_criteria(temp_train,cca_options);
        dt(k)         = trunc_t.ind;
        
        [~,trunc_p]   = truncation_criteria(prox_train,cca_options);
        dp(k)         = trunc_p.ind;
%         if dp(k) < dt(k)
%             dp(k)     = dt(k);
%         end
        %dp(k) = min(25,dt(k));
        disp(['dt = ',num2str(dt(k)),', dp = ',num2str(dp(k))])
        
        %================================================================
        % Choose truancation parameter for cov(T,P) using k-fold CV
        %================================================================
        F             = VT(:,1:dt(k))'*VP(:,1:dp(k));
        [UF,SF,VF]    = svd(F);
        
        dcca_max      = min(dp(k),dt(k));
        dcca_vec      = [1: dcca_max];
        
        for dcca = dcca_vec
            % Estimate CV regression matrix
            B_cv         = UT(:,1:dt(k))*ST(1:dt(k),1:dt(k))*UF(1:dt(k),1:dcca)...
                *SF(1:dcca,1:dcca)*VF(1:dp(k),1:dcca)'*...
                diag(1.0./diag(SP(1:dp(k),1:dp(k))))*UP(:,1:dp(k))';
            
            % Predict temperature over verification sample
            temp_pred    = (prox-ones(ntest,1)*Pm)*diag(1.0./Ps)*B_cv'*diag(Ts)+ ones(ntest,1)*Tm;
            
            % Calculate the misfit (RMSE)
            mse          = mean((temp_pred(test,:) - temp_test).^2);
            rmse(dcca,k) = sqrt(mse*weights);
            
        end
    end
    mean_misfit          = sqrt(nmean(rmse(1:min(dt),:).^2,2));
    [min_RMSE, dcca_opt] = min(mean_misfit(:));
    % dcca_opt = ind2sub(size(mean_misfit),find(mean_misfit == min_RMSE));
    dt_opt               = min(dt);
    dp_opt               = min(dp);
end

% save(['rmse_cv_' noise '.mat'],'rmse','indices');

display('Computation is done!')


% Save output
misfit.rmse = rmse;
misfit.min  = min_RMSE;


end
