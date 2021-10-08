function [betax,paramsx,beta,params] = lowesstatx(xi,yi,x,p,f,h,N,varargin)
% "LOWESS at X" applies the LOcally WEighted Scatter plot Smoothing (LOWESS) 
% estimator to time series of sea surface temperature observations as implemented 
% in Elipot et al. 2021. This method is adapted from Cleveland (1979) ﻿doi: 10.2307/2286407.
% 
% The method fits SST observations within a kernel of constant bandwith to a model 
% constituted of two parts: a polynomial function and a sum of cosine and sine 
% function pairs at chosen frequencies.
% 
% Usage
% 
% [BETAX,PARAMSX,BETA,PARAMS] = lowesstatx(XI,YI,X,P,F,H,N)
% XI are observations times 
% YI are the SST observations
% X are the target times at which the estimation is utlimately
% conducted, and can be empty i.e. x = []
% or can coincide to the input times if X is set to XI
% P is order the polynomial s_P = \sum_p=0^P (x-xi)^p
% F is row vector of cyclic frequencies (in units of cycles per unit time)
% input can be empty i.e. F = []
% H is the bandwith of the window kernel (of length 2H)
% N determines the N+1 number of iterations of the estimation
% 
% [BETAX,PARAMSX,BETA,PARAMS] = lowesstatx(XI,YI,X,P,F,H,N,D) 
% optionally specifies the factor within the denominator of the argument for the function
% calculating the robust weights as part of the iteration. Cleveland (1979)
% sets D=6. Default value is D=14 as explained in Elipot et al. 2021.
% Lower values of D result in more data points being rejected as outliers.
% 
% Outputs
% 
% with M being the number of input data points and MX being the number of target times
% the output is as follows:
%
% BETAX is a cell array of arrays of size MX x (P+1+2 x length(F)) containg the estimated
% parameters of the model
% PARAMSX is a cell array of structure variables with the following fields
% PARAMSX{m}.H     MX x 1 vector of final bandwidth of estimation kernel
% PARAMSX{m}.ND  MX x 2  matrix of number of data points per estimation (first dimension) 
 %                         and effective number of degree of freedom
 %                         (second dimension)
% PARAMSX{m}.VBETA = vbeta; % covariance matrix of parameter estimates, of size
% MX x (P+1+ 2 length(F)) x (P+1+2 length(F))
% PARAMSX{m}.RSS ; % MX x 1 vector of normalized weighted residual sum of squares = error variance estimate

% BETA is a cell array of arrays of size M x (P+1+2 x length(F)) containg the estimated
% parameters of the model
% PARAMS is a cell array of structure variable with the following fields 
% PARAMS{m}.H       M x 1 vector of final bandwidth of estimation kernel
% PARAMS{m}.R       M x 1 vector of residuals
% PARAMS{m}.D       M x 1 vector of penultimate weights used for estimation "at x"
% PARAMS{m}.M      (N+1) x 1 vector of median values of absolute residuals per iteration
% PARAMS{m}.ND    M x 2 matrix of number of data points per estimation (first dimension) 
 %                               and effective number of degree of freedom
 %                               (second dimension)
% PARAMS{m}.VBETA = vbeta; % covariance matrix of parameter estimates, of size
% M x (P+1+ 2 length(F)) x (P+1+2 length(F))
% PARAMS{m}.RSS ; % vector of length M of normalized weighted residual sum of squares = error variance estimate

% Shane Elipot, 2021, version 1

% how many optional input argument?
nin = nargin - 7;
if nin == 0
    D = 14; % set the factor of denominator of argument of robust weights calculation
                 % Cleveland 1979 uses D = 6;
elseif nin == 1
    D = varargin{1};
else
    error('Only one optional argument is allowed');
end

% sort input data
[xi,I] = sort(xi);
yi = yi(I);

xi = xi(:);
yi = yi(:);
x = x(:);

[n1,~] = size(xi);
[nx,~] = size(x);

% number of frequencies
n = length(f);

% initialize the outputs
% solution at the estimation times

d = ones(n1,1); % LOWESS weights
M = NaN*ones(N+1,1); % median of residuals

% N+1 is the number of fits, N iterations plus the original one with d = 1 for all
% points

for m = 1:1+N
    
    disp(['Iteration ' num2str(m)]);
    
    r = NaN*ones(n1,1); % residuals
    nd = NaN*ones(n1,2); % number of points used per estimation 
                                        % and effective number of degrees of freedom
    rss = NaN*ones(n1,1); % error variance estimate
    beta = NaN*ones(n1,p+1+2*n); % parameter solution
    hb = NaN*ones(n1,1); % final bandwidth
    vbeta = NaN*ones(n1,(p+1+2*n),(p+1+2*n)); % covariance matrix
    
    for k = 1:length(xi)
        
        hk = h;
        
        if 2*hk < (xi(end)-xi(1)) % if the window does not exceed the size of the time series
            
            g = 0; % a conditional test on the invertible matrix for least squares estimation
            
            while (g < 1) &&  2*hk < 4 % while the matrix conditioning is large
                                                     % and twice the bandwidth does not increase
                                                     % beyond 4 days
                
                w = kernel1(xi-xi(k),hk);% create the kernel window centered on current observation
                w = w.*d;                      % modify that kernel by the robust weights; originally all 1
                
                % need to make sure the problem is not underdetermined:
                % if it is, no estimation possible
                q = find(w~=0);
                
                if length(q)>=p+1+2*n  % there are technically enough points for the estimation
                                            
                    % set up design matrix and weight matrix
                    % reduce input data to points with non-zero weights
                    w = w(q);
                    xi2 = xi(q);
                    yi2 = yi(q);
                    
                    % design matrix
                    X = zeros(length(xi2),p+1);
                    z = xi2-xi(k);
                    z = z(:);
                    for j = 0:p
                        X(:,j+1) = z.^j;
                    end                    
                    if ~isempty(f)
                        X = [X cos(2*pi*z*f) sin(2*pi*z*f)];
                    end
                    
                    % weigh matrix
                    W = diag(w,0);
                    
                    Y = yi2(:);
                    R = transpose(X)*W*X;
                    cR = cond(R);  % condition number of matrix to be inverted
                    
                    if cR> 10^4           % condition number is  poor -> blow up of inversion
                        hk = hk + 1/24; % increase bandwidth by 1/24 (here = 1 hour)
                        clear R
                    else % the matrix condition number is ok
                        g = 1;
                    end
                    
                else % not enough points to do the fit, increase bandwith by 1/24 (= 1 hr)
                    
                    hk = hk + 1/24;
                    %clear R
                    % disp('increase window by 1 hr');
                    
                end % now test again if g<1 and that the bandwidth is not too large
                
            end % condition number of matrix is reasonable and a solution can be obtained
            
            if exist('R','var')
                
                % R = transpose(X)*W*X;
                A = R\transpose(X)*W;
                % LSQ solution
                betah = (A*Y);
                
                % save results at every iteration
                beta(k,:) = betah';
                % solution: not returned see sstestimate.m
                % s_k = sum(beta(k,[1 p+1+1:p+1+n]),2);
                % residuals = data minus solution
                r(k) = yi(k) - sum(beta(k,[1 p+1+1:p+1+n]),2);
                
                % normalized weighted residual sum of squares (nwrss):
                % Eq. 4.8 of Fan and Gijbels
                deno = trace(W-W*X*inv(X'*W*X)*X'*W);
                nwrss = (Y - X*betah)'*W*(Y - X*betah)/deno;
                rss(k) = nwrss;
                
                % variance of the estimates
                % varbetah = R\X'*W*(rss*eye(length(q)))*W*X/R;
                % Eq. 3.6 of Fan and Gijbels:
                S = nwrss*W.^2;
                varbetah = R\X'*S*X/R;
                vbeta(k,:,:) = varbetah;
                
                % final bandwidth
                hb(k) = hk;
                % number of points
                nd(k,1) = length(q);
                nd(k,2) = deno;
                      
                clear R
                
            else % no solution but keep the size of the window and final number of points
                
                hb(k) = hk;
                nd(k,1) = length(q);
                
            end
            
        end % end of condition on length of window
        
        clear R
        
    end % end of loop on each observation
    
    % median of absolute residuals
    %M(m) = nanmedian(abs(r));
    M(m) = median(abs(r),'omitnan');    
    % new choice of weights; Cleveland (1979) chose D=6 in denominator;
    d = kernel2(r/(D*M(m)));
    % when no solution has been obtained, r=NaN and d = NaN so replace by
    % zero to effectively reject the point as "outlier" with no achievable
    % solution
    d(isnan(d)) = 0;
    
    if m == N % save the final weights before the final estimation
        dfinal = d;
    end
    
end

% save structure  output
params.h = hb; % final bandwidth
params.r = r; % residuals
params.d = dfinal; % penultimate weights
params.M = M; % vector of median of absolute residuals
params.nd = nd; % number of data points per estimation and effective number of degree of freedom
params.vbeta = vbeta; % covariance matrix of parameter estimates
params.rss = rss; % normalized weighted residual sum of squares = error variance estimate

if nx ~= 0
    
    % now estimate at specified times ("at x")
    % and get all parameters
    
    betax = NaN*ones(nx,p+1+2*n); % parameter estimates
    hb = NaN*ones(nx,1); % final bandwidth
    nd = NaN*ones(nx,2); % number of points used per estimation
    vbeta = NaN*ones(nx,(p+1+2*n),(p+1+2*n)); % covariance matrix of parameter estimates
    rssx = NaN*ones(nx,1); % error variance estimate
    
    for k = 1:length(x)
        
        hk = h;
        
        if 2*hk < (x(end)-x(1)) % if the window does not exceed the size of the time series
            
            g = 0; % a conditional test on the invertible matrix for least squares estimation
            
            while (g < 1) &&  2*hk < 4 %(xi(end)-xi(1)) % while the matrix conditioning is large
                % and twice the bandwidth does not increase
                % beyond 4 days
                
                w = kernel1(xi-x(k),hk); % create the kernel window centered on estimation time
                w = w.*dfinal;                      % modify that kernel by the robust weights
                
                % need to make sure the problem is not underdetermined:
                
                q = find(w~=0);
                
                if length(q)>=p+1+2*n  % there are technically enough points; set up design matrix and weight matrix
                    
                    w = w(q);
                    xi2 = xi(q);
                    yi2 = yi(q);
                    
                    X = zeros(length(xi2),p+1);
                    z = xi2-x(k);
                    z = z(:);
                    for j = 0:p
                        X(:,j+1) = z.^j;
                    end
                    
                    if ~isempty(f)
                        X = [X cos(2*pi*z*f) sin(2*pi*z*f)];
                    end
                    
                    W = diag(w,0);
                    Y = yi2(:);
                    R = transpose(X)*W*X;
                    cR = cond(R);  % condition number of matrix to be inverted
                    
                    if cR> 10^4 % condition number is too poor -> blow up of inversion
                        hk = hk + 1/24;
                        clear R
                    else % the matrix condition number is ok
                        g = 1;
                    end
                    
                else % not enough points to do the fit, increase bandwith by 1 hr
                    
                    hk = hk + 1/24;
                    % clear R
                    % disp('increase window by 1 hr');
                    
                end % now test again if g<1 and that the bandwidth is not too large
                
            end % condition number of matrix is reasonable and a solution can be obtained
            
            if exist('R','var')
                
                A = R\transpose(X)*W;
                % solution
                betah = (A*Y);
                % solution at estimation times
                %s(k) = sum(betah([1 p+1+1:p+1+n]));
                betax(k,:) = betah';
                % there are no residual here
                
                % normalized weighted residual sum of squares:
                % Eq. 4.8 of Fan and Gijbels
                deno = trace(W-W*X*inv(X'*W*X)*X'*W);
                nwrss = (Y - X*betah)'*W*(Y - X*betah)/deno;
                rssx(k) = nwrss;
                
                % variance of the estimates
                % varbetah = R\X'*W*(rss*eye(length(q)))*W*X/R;
                % Eq. 3.6 of Fan and Gijbels:
                S = nwrss*W.^2;
                varbetah = R\X'*S*X/R;
                vbeta(k,:,:) = varbetah;
                
                nd(k,2) = deno;
               
                clear R
                
            end
            
            % final bandwidth
            hb(k) = hk;
            % number of points
            nd(k,1) = length(q);
            
        end % end of condition on length of window
        
    end
    
    paramsx.vbeta = vbeta; % covariance matrix
    paramsx.h = hb; % final bandwidth
    paramsx.nd = nd; % number of data points and degrees of freedom
    paramsx.rss = rssx; % error variance estimate
    
else % no estimation if no output argument specified
    
    betax = [];
    paramsx = [];

end

function K = kernel1(t,h) % tricube kernel function
% the kernel is zero outside of the normalized bandwith 1 by construction
K = (1-abs(t/h).^3);
qn = K<0;
K(qn) = 0;
K = K.^3; % Cleveland 1979 version
%K = (1/h)*(70/81)*K.^3; % Fan and Gijbels version; factor is not important

function K = kernel2(t) % biweight kernel function
% the kernel is zero outside of the normalized bandwith 1 by construction
K = (1-abs(t).^2);
qn = K<0;
K(qn) = 0;
K = K.^2;
