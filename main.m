% main code

% load the example dataset for drifter with GDP AOML ID 55366 (WMO ID 3100541)
load('level0_sst_data_drifter_id55366.mat','presst','timepresst','numx');
% timepresst is a cell array with observation times in days as returned by Matlab
% datenum.m
% presst is a cell array with the corresponding SST observations
% numx is a cell array the target times for estimation

%% estimation routine

% choice of model parameters:
p = 1; % polynomial order
f = [1 2 3]; % peridoifrequencies of model
% choice of estimation parameters
N = 3; % iteration number is N+1
bw = 1; % starting kernel bandwidth; cannot be equal or larger than 2
D = 14; % factor for robust weight calculation; see Elipot et al. 2021

for m = 1:1
    tic;
    [betax{m},paramsx{m},beta{m},params{m}] = ...
        lowesstatx(timepresst{m},round(1000*presst{m})/1000,numx{m},p,f,bw,N,D);
    toc;
end

% calculate sst estimates at original times
% specify the resolution "a" of the input data
a = 0.05*ones(size(beta));
[sst1,esst1] = ssteval(beta,params,p,f,a,'background');
[sst2,esst2] = ssteval(beta,params,p,f,a,'diurnal');
[sst,esst] = ssteval(beta,params,p,f,a);

% calculate sst estimates at target times (numx)
[sstx1,esstx1] = ssteval(betax,paramsx,p,f,a,'background');
[sstx2,esstx2] = ssteval(betax,paramsx,p,f,a,'diurnal');
[sstx,esstx] = ssteval(betax,paramsx,p,f,a);

%% calculate the statistics used to assess model
R = 10^5; % rounding to 5 decimal places
m = 1;
r = round(R*params{m}.r)/R; % rounded residuals
w = params{m}.d./sum(params{m}.d); % normalized weights
rss = round(R^2*params{m}.rss)/R^2; % rounded error variance
q = isfinite(r) & isfinite(rss);        % successful estimations

% weighted root mean square error
wrmse = sum((r(q).^2).*w(q)).^0.5;

% weighted median of error variance
wmed = weightedMedian(rss(q),w(q));

disp(['weighted root mean square error = ' num2str(wrmse) ' K']);
disp(['square root of weighted median variance = ' num2str(wmed.^0.5) ' K']);

%% figures
cc = lines(4); % colorscale for curves
fac = 1.96; % standard error of estimates times 1.96 will approximately give 95% confidence intervals

m = 1;

q = params{m}.d == 0; % outlier points

% figure of  estimates at original times
figure, hold on
h0 = plot(timepresst{m},round(1000*presst{m})/1000,'+');
ho = plot(timepresst{m}(q),round(1000*presst{m}(q))/1000,'ko');
h1 = plot(timepresst{m},sst1{m},'x');
h2 = plot(timepresst{m},sst2{m},'x');
h12 = plot(timepresst{m},sst{m},'x');
datetick('x');
title('SST estimates at original observation times');
legend([h0 h1 h2 h12 ho],{'Input data','Non-diurnal SST','Diurnal SST anomalies','Total SST','Outliers'});
ylabel('Degree Celsius');

% figure of regular continuous estimates with confidence intervals
figure, hold on
h0 = plot(timepresst{m},round(1000*presst{m})/1000,'.','color',0.5*[1 1 1]);
ho = plot(timepresst{m}(q),round(1000*presst{m}(q))/1000,'ko');
h1 = plot(numx{m},sstx1{m},'linewidth',1.5,'color',cc(1,:));
h11 = plot(numx{m},sstx1{m}+fac*esstx1{m}*[-1 1],'linewidth',0.5,'color',cc(1,:));
h2 = plot(numx{m},sstx2{m},'linewidth',1.5,'color',cc(2,:));
h22 = plot(numx{m},sstx2{m}+fac*esstx2{m}*[-1 1],'linewidth',0.5,'color',cc(2,:));
h3 = plot(numx{m},sstx{m},'linewidth',1.5,'color',cc(3,:));
h33 = plot(numx{m},sstx{m}+fac*esstx{m}*[-1 1],'linewidth',0.5,'color',cc(3,:));
datetick('x');
legend([h0 h1 h2 h3 ho],{'Input data','Non-diurnal SST','Diurnal SST anomalies','Total SST','Outliers'});
title('SST estimates at regular hourly time steps');
ylabel('Degree Celsius');

% It is best to zoom manually on the figure to see the results

return


