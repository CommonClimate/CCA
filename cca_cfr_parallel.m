function [field_r, diagn] = cca_cfr_parallel(field,proxy,calib,cca_options)
% Function [field_r,diag] = cca_cfr(field,proxy,t,calib,cca_options)
%
%  inputs:  - field, 2D climate field (nt x ns)
%           - proxy, proxy matrix (nt x nr)
%           - t, time axis  (nt)
%           - calib, index of calibration (< nt)
%           - cca_options, structure. in particular:
%              # method: 'smerdon10' (default)
%                        'NE08' 
%              # K: value of k for k-fold CV, K = 2 (default)
%              # early/late: index of period for the 'smerdon10' style CV
%              # block: option for full vs sparse matrices (default = 1,
%              i.e. sparse network, assuming that there are missing values
%              in the input data)
%
%  outputs: - field_r, reconstructed field
%           - diagn, structure of diagnostic outputs
%              # params: calculated CCA parameters (dp, dt, dcca)
%
%
% ====================================================


% Process options
if nargin < 3 || isempty(cca_options)
    fopts      = [];
else
    fopts      = fieldnames(cca_options);
end

%
if sum(strcmp('params', fopts)) ~= 0
    params_given = true;
else
    params_given = false;
end

if sum(strcmp('block', fopts)) ~= 0
    block      = cca_options.block;
else
    block      = false;
end

if sum(strcmpi('dp_max',fopts))>0
    dp_max     = cca_options.dp_max;
else
    dp_max     = 50; 
end

if sum(strcmpi('dt_max',fopts))>0
    dt_max     = cca_options.dt_max;
else
    dt_max     = 50;
    cca_options.dt_max = dt_max;
end

if sum(strcmpi('params_range',fopts))>0
    dp_range         = cca_options.params_range(2,:);
    dt_range         = cca_options.params_range(3,:);
end

[nt np]        = size(proxy);
field_r        = zeros(size(field));

kavlr          = cell(nt,1);
kmisr          = cell(nt,1);
for j=1:nt
    kavlr{j}   = find(~isnan(proxy(j,:)));
    kmisr{j}   = find(isnan(proxy(j,:)));
end

% search for unique patterns
avail          = double(~isnan(proxy));
avail_uniq     = unique(avail,'rows','first');
nps            = size(avail_uniq,1); % number of patterns

fieldc         = field(calib,:);
proxc          = proxy(calib,:);
nc             = length(calib);

if params_given == false
    K = cca_options.K;
    if K == 2
        early             = cca_options.early;
        late              = cca_options.late;
        indices(early)    = 1;
        indices(late)     = 2;
        indices           = indices(calib);
    elseif K ~= 2
        indices           = crossvalind('Kfold',nc,K);
    end
    cca_options.indices   = indices;
else
    dcca(j) = cca_options.params(1);
    dp(j)   = cca_options.params(2);
    dt(j)   = cca_options.params(3);
end
if block
    % cycle over patterns of missing values
    parfor j=1:nps             
        options    = cca_options;
        % Fill in the pattern cell array
        display(['Pattern ' num2str(j) '/' num2str(nps)])
        avail_m    = avail - repmat(avail_uniq(j,:),nt,1);
        pattern    = find(std(avail_m,0,2) == 0)';
        jp         = min(pattern);             % position of this pattern
        pm         = length(kmisr{jp});        % number of missing values in this pattern
        
        mp         = length(pattern);          % number of rows matching this pattern
        pa         = np - pm;                  % number of available values in this pattern
        
        prox_p     = proxc(:,kavlr{jp});
        
        options.dp_max = min(pa,dp_max);       % Max value of truncation parameters
        
        % Compute CCA parameters using cross-validation
        [dp, dt, dcca, misfit] = cca_cv(fieldc, prox_p, options);
        
        % Compute regression coefficients per pattern
        [~,B,Tm,Pm,Ts,Ps]      = cca_bp(prox_p,fieldc,dp,dt,dcca);
        
        
        
        % Perform CCA reconstruction per pattern
        display(['dcca = ', num2str(dcca), ', dp = ', num2str(dp), ', dt = ', num2str(dt)]);
        field_rc = (proxy(pattern,kavlr{jp}) - repmat(Pm,[mp 1]))*diag(1.0./Ps)*B'*diag(Ts);
        field_r  = field_rc  + repmat(Tm,[mp 1]);
        fname = sprintf( 'params_%d', j );
        iSaveX( fname, dp, dt, dcca, field_r, pattern );
        % save('cca_params.mat','dcca','dp','dt');
    end
else
    if params_given == false  % Choose the optimal truncation parameters first
        K = cca_options.K;
        if K == 2
            early          = cca_options.early;
            late           = cca_options.late;
            indices(early) = 1;
            indices(late)  = 2;
            indices        = indices(calib);
        elseif K ~= 2
            indices        = crossvalind('Kfold',nc,K);
        end
        cca_options.indices    = indices;
        [dp, dt, dcca, misfit] = cca_cv(fieldc, proxc, cca_options);
    else
        dcca = cca_options.params(1);
        dp   = cca_options.params(2);
        dt   = cca_options.params(3);
    end
    display(['dcca = ', num2str(dcca), ', dp = ', num2str(dp), ', dt = ', num2str(dt)])
    
    % Compute regression coefficients and perform reconstruction
    [~,B,Tm,Pm,Ts,Ps] = cca_bp(proxc,fieldc,dp,dt,dcca);
    field_rc          = (proxy - repmat(Pm,[nt 1]))*diag(1.0./Ps)*B'*diag(Ts);
    field_r           = field_rc  + repmat(Tm,[nt 1]);
end

for j = 1:nps
    fname = sprintf( 'params_%d', j );
    a = load(fname);
    dp(j)   = a.dp;
    dt(j)   = a.dt;
    dcca(j) = a.dcca;
    pattern{j}            = a.pattern;
    field_r(pattern{j},:) = a.field_r;
end
% diagn.misfit = misfit;
diagn.params = [dcca; dp; dt;];


end


function iSaveX( fname, dp, dt, dcca, field_r, pattern )
save( fname, 'dp','dt','dcca','field_r','pattern' );
end

