% Employ Canonical Correlation Analysis to perform reconstruction. The
% example here is shown using a pseudoproxy network that has missing values
% to mimick proxies' decreasing temporal availability, and the script can
% be easily modified to use on real-world proxy networks.

% If the proxy matrix is full rank, use cca_cfr.m instead of
% cca_cfr_parallel.m (see below line 74)

% Jianghao Wang @USC, May 2012


clear all; clc; close all;
%======================================================
% Define networks, tags for labeling purposes
%======================================================
method = 'cca'; % other choices are cca, ppr, barcast, etc...
ctag   = 'm08';
avail  = 'sparse';  % or full, defines whether the data is temporally invariant
model  = 'csm1p4';% GCMs to be used

dirname = [model '_'  ctag '_' avail '_' method ]; % Experiment name


if ~exist(dirname,'dir') % If it doesn't exist already
    system(['mkdir ', dirname]);% create experiment directory
end
noise = {'local','max'};

if ~exist('figures','dir')
    system(['mkdir ', 'figures']);
end

%======================================================
% Load proxy and temperature data
%======================================================
i = 1; c = 1; n = 1; % Pick one realization to perform the reconstruction
fname = ['csm1p4_pproxy_m08_',avail,'_','realistic','_'];
num = [sprintf('%02d',(i-1)*10+1),'-',num2str(i*10)];
load([fname,num])
nc = size(pproxy,2);
M  = size(pproxy{c},3);

% Define time
T     = ptime;           nt    = length(T);
cbeg  = find(T == 1850); cend  = find(T == 1995);  cmid = find(T == 1922);
calib = [cbeg:cend]';    early = [cbeg:cmid]';     late = [cmid+1:cend]';
nh    = length(calib);   % early and late for leave-half-out CV

% Define temperature field (and weights) subsampled on location grid
temp = temp(1:cend,:);
weights0=cosd(locs_pnas(:,2));
weights=weights0/sum(weights0);
cutoff = 20;

% Tags
stag = [dirname,'_' num];
snr{1} = case_vec{5};
snr{2} = case_vec{6};

% Define CCA options
cca_options.method       = 'smerdon10';   % can be NE08 as well
%cca_options.dp_max      = 50;  % maximum dp, dt
cca_options.K            = 2;        % user-defined, for k-fold CV
cca_options.weights      = weights;
cca_options.early        = early;
cca_options.late         = late;
cca_options.block        = 1;

if strcmpi(method,'cca')== 1
    matlabpool('open','CLIM_72')  
    % needs to be changed to your own parallel configurations
end

proxy = pproxy{c}(:,:,n);
display(['Reconstructing SNR case ',int2str(c),'/', int2str(nc),', noise realization ',int2str(n),'/',int2str(M),',',num]);

%====================================================
% Perform CCA reconstruction and save diagnostics
%====================================================
[temp_r{c,n}, diagn{c,n}] = cca_cfr_parallel(temp,proxy,calib,cca_options);

save([stag, '_step02'],'temp','temp_r','diagn','locs_pnas','T')

% Close matlab pool
if strcmpi(method,'cca')== 1
    matlabpool close
end
%
% %
