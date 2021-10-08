function [sst,esst,varargout] = ssteval(beta,params,p,f,a,varargin)
% This functions takes the set of outputs BETA and PARAMS from lowesstatx.m
% and return the SST estimates
%
% Usage
%
% [SST,ESST] = ssteval(BETA,PARAMS,P,F,A) returns the total SST estimate
% that is the sum of the diurnal and polynomial (non-diurnal) components of
% the estimation model
% 
% BETA is a cell array of matrices and PARAMS is a cell array of structure
% variables as returned by lowesstatx.m
% P is the order of the polynomial of the model
% F is a row vector of cyclic frequencies (in units of cycles per unit
% time) of the model
% Input variable A is the resolution of the input data which is used to correct the 
% error variance estimate for the resolution or quantization error. A must
% be a vector of the same length as beta and params
% A typical SST equation for a drifter SST sensor is SST = A*n +B where n is a bit 
% count. A is the resolution of the input data. If left empty A is set to
% zero.

% [SST1,ESST1] = ssteval(BETA,PARAMS,P,F,'background')
% optionally returns only the (non-diurnal) background polynomial component
% of the model
% [SST,ESST2] = ssteval(BETA,PARAMS,P,F,'diurnal')
% optionally returns only the diurnal component
% of the model
% 
% Note that if you first call this function to obtain the background
% component and second call it again to get the diurnal component you
% should obtain SST1+SST2=SST but not ESST1+ESST2=ESST because of cross
% covariance terms. See Elipot et al. 2021

% initialize outputs and needed variables
sst = cell(size(beta));
esst = sst;
vbeta = sst;

n = length(f); % number of frequencies

if isempty(a) % set quantization error to zero
    a = 0*ones(size(beta));
end

% check for size of a
if length(a) ~= length(beta)
    error('Input parameter A must be of the same length as BETA and PARAMS');
end

% rescale the error variance estimate
resvar = a.^2/12;

for m = 1:length(beta)

    % add the resolution variance
    rss_tmp = params{m}.rss+resvar(m);
    % first rescale the covariance matrix
    params{m}.vbeta = bsxfun(@times,[params{m}.vbeta],rss_tmp./params{m}.rss);
    % second replace the error variance
    params{m}.rss = rss_tmp;
end

for m = 1:length(beta)
    vbeta{m} = params{m}.vbeta;
end

if isempty(varargin) % returns only sst and total error
    
    fac1 = 1;
    fac2 = 1;
    facmatrix = ones(1,n+1,n+1);
    
elseif strcmp(varargin{1},'background') %  only background with sst tendency

    fac1 = 1;
    fac2 = 0;
    facmatrix = zeros(1,n+1,n+1);
    facmatrix(1,1,1) = 1;
        
elseif strcmp(varargin{1},'diurnal') % only diurnal with a and phi
    
    fac1 = 0;
    fac2 = 1;   
    facmatrix = fac2*ones(1,n+1,n+1);
    facmatrix(1,1,:) = fac1;
    facmatrix(1,:,1) = fac1;   
    
end

for m = 1:length(beta)
    sst{m} = fac1*beta{m}(:,1) + fac2*sum(beta{m}(:,p+1+1:p+1+n),2);
    % taking the real part in case of negative result
    esst{m} = real(sum(sum(bsxfun(@times,facmatrix,vbeta{m}(:,[1 p+1+1:p+1+n],[1 p+1+1:p+1+n])),3),2).^0.5);     
end

return
