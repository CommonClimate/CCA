function [Tp,B,Tm,Pm,Ts,Ps] = cca_bp(P,T,dp,dt,dcca)
% function to make predictions using CCA-BP
%
% Input parameters:
% P - matrix of predictor vectors (columns)
% T - matrix of predictand vectors (columns)
% Note: P and T should have the same row dimension
% dp and dt are reduction dimensions for P and T respectively, 
% per Barnett and Preisendorfer's CCA version
% dcca is the number of CCA modes to be used for prediction
% Note: dcca cannot exceed min(dp,dt)
% 
% Output parameters:
% Tp - in-sample prediction for T
% B  - regression matrix
% Tm - average of T in row dimension
% Pm - average of P in row dimension
% Ts - std of T in row dimension
% Ps - std of P in row dimension
%
% For explanations see the paper 
% Smerdon, J.E., A. Kaplan, D. Chang, and M.N. Evans, 2010: 
% A pseudoproxy evaluation of the CCA and RegEM methods for reconstructing 
% climate fields of the last millennium, J. Climate, in press.
% http://www.ldeo.columbia.edu/~jsmerdon/papers/2010a_jclim_smerdonetal.pdf
%
%  A.Kaplan, J.E.Smerdon, August 2010 
% Modified by J. Wang, USC, April 2012

nc = size(P,1);
% Time-standardized matrices

[Pds,Pm,Ps] = standardize(P);
[Tds,Tm,Ts] = standardize(T);

[UP,SP,VP] = svd(Pds');
[UT,ST,VT] = svd(Tds');

% Cross Covariance matrix
F = VT(:,1:dt)'*VP(:,1:dp);
[UF,SF,VF] = svd(F);


B=UT(:,1:dt)*ST(1:dt,1:dt)*UF(1:dt,1:dcca)*SF(1:dcca,1:dcca)*VF(1:dp,1:dcca)'*diag(1.0./diag(SP(1:dp,1:dp)))*UP(:,1:dp)';

Tdsp = B*Pds';
% Tp = diag(Ts')*Tdsp + Tm'*ones(1,nc);
Tp = Tdsp'*diag(Ts) + ones(nc,1)*Tm;
