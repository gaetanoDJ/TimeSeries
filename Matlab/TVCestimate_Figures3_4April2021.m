% This code computes the time-varying coefficients using the method
% developed in Boivin (JMCB 2006)
% Key modifications:
% a) convert policy persistence into the same quarterly frequency;
% b) allow for three breaks in the size of volatility:
%    pre-1979,
%    1979-1982,
%    post-1982;
% c) do Kozicki and Tinsley (2007) caluclation of the target rate and
% compute the corresponding determinacy rates


% Coefficeints int he Taylor rule are assumed to follow random walk

clear
% close all
randn('seed',1234567)

num_sim_cov=1000; % number of draws to simulated the covariance matrix of converted smoothness in the Taylor rule
its=1000;         % number of draws to compute determinacy


%========================================================================
% Step 1: import data and construct RHS and LHS
%========================================================================

%QUESTION 1. EXPLAIN HOW DATA IS IMPORTED. WHAT IS T, Y, X, inflation?

disp('Step 1: Import data and run OLS regressions ...')
import_data_greenbook;

Y=data_greenbook(:,4);
X=data_greenbook(:,5:10);
inflation=data_greenbook(:,11); % GDP deflator smoothed over current and past 12 quarters
T=length(Y);

trend_inflation=data_greenbook(:,12);
trend_inflation_date1=129;
trend_inflation_date2=157;

% create x-axis labels
labelfigs=[1 data_greenbook(1,2)];
stepfig=2;  % two year step for labels in figures
for i=2:T
    if labelfigs(end,2)+stepfig-1<data_greenbook(i,2)
        labelfigs=[labelfigs; i data_greenbook(i,2)];
    end
end

%========================================================================
% Step 2: estimate the size of the shocks to coefficients
%========================================================================
disp('Step 2: Run OLS regressions and estimate variance of shocks ...')
% The Stock and Watson (JASA 1998) parameter to the size of the shocks to
% coefficients in the Taylor rule. S&W method is used to have a median
% unbiased estimate for the size of the shock. Otherwise, the variances of
% innovations in coefs collapse to zero.

disp(['Quandt Liklihood ratio test statistic: ' num2str(QLR(Y,X))])
lambda=8; % This parameter is taken from Stock and Watson (JASA 1998)

% estimated Taylor rule:
% coef #1 = expected inflation rate, E(PI{t+1}|I_{t}), pi_tp
% coef #2 = current growth rate of output, Y{t}, gry_t
% coef #3 = current output gap, gap{t}, gap{t}, gap_t00
% coef #4 = FFR{t-1}
% coef #5 = FFR{t-2}
% coef #6 = constant terms

% ALLOW FOR THREE PERIODS WITH DIFFERENT VOLATILITIES:
% a) pre-1979
% b) 1979-1982"
% c) post-1982

%QUESTION 2. WHAT IS THE MAIN PURPOSE OF THIS SECOND STEP? WHY THE AUTHORS BREAK TIMELINE INTO 3 PERIODS? 
%WHAT ARE "beta", "resid", "R", "SigmaXX"?

% a) pre-1979
X1=X(1:130,:);
Y1=Y(1:130,1);

beta1=inv(X1'*X1)*(X1'*Y1); % OLS estimate of the policy reaction function across all periods
resid=Y1-X1*beta1;

R=var(resid);   % variance of the error term

SigmaXX=(X1'*X1)/length(X1);

EXeeX=0; % White heteroskedasticity consistent estimate
for i=1:length(X1)
    EXeeX=EXeeX+resid(i)^2*X1(i,:)'*X1(i,:);
end
EXeeX=EXeeX/length(X1);

SigmaVV1=inv(SigmaXX)*EXeeX*inv(SigmaXX);

% estimate of innovations to coefficients
SigmaBB1=(lambda/T)^2*SigmaVV1;
R1=R;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% a) 1979-1982
X1=X(131:159,:);
Y1=Y(131:159,1);

beta2=inv(X1'*X1)*(X1'*Y1); % OLS estimate of the policy reaction function across all periods
resid=(Y1-X1*beta2);

R=var(resid);   % variance of the error term

SigmaXX=(X1'*X1)/length(X1);

EXeeX=0; % White heteroskedasticity consistent estimate
for i=1:length(X1)
    EXeeX=EXeeX+resid(i)^2*X1(i,:)'*X1(i,:);
end
EXeeX=EXeeX/length(X1);

SigmaVV2=inv(SigmaXX)*EXeeX*inv(SigmaXX);

% estimate of innovations to coefficients
SigmaBB2=(lambda/T)^2*SigmaVV2;
R2=R;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% a) post-1982
X1=X(160:end,:);
Y1=Y(160:end,1);

beta3=inv(X1'*X1)*(X1'*Y1); % OLS estimate of the policy reaction function across all periods
resid=Y1-X1*beta3;

R=var(resid);   % variance of the error term

SigmaXX=(X1'*X1)/length(X1);

EXeeX=0; % White heteroskedasticity consistent estimate
for i=1:length(X1)
    EXeeX=EXeeX+resid(i)^2*X1(i,:)'*X1(i,:);
end
EXeeX=EXeeX/length(X1);

SigmaVV3=inv(SigmaXX)*EXeeX*inv(SigmaXX);

% estimate of innovations to coefficients
SigmaBB3=(lambda/T)^2*SigmaVV3;
R3=R;


R2=R1;
R3=R1;
SigmaBB2=SigmaBB1;
SigmaBB3=SigmaBB1;



%========================================================================
% Step 3: HP filter on series
%========================================================================

disp('Step 3: HP filter series ...')

HPlambda=3000000;
rrate_matHP=Y(:,1)-X(:,1)-hpfilter(Y(:,1)-X(:,1),HPlambda,0);
gy_matHP=X(:,2)-hpfilter(X(:,2),HPlambda,0);
gap_matHP=X(:,3)-hpfilter(X(:,3),HPlambda,0);

% ALLOW for breaks in the trend and target
% three breaks
% pre-1979
gy_matHP_br0(1:159,1)=X(1:159,2)-hpfilter(X(1:159,2),HPlambda,0);
gap_matHP_br0(1:159,1)=X(1:159,3)-hpfilter(X(1:159,3),HPlambda,0);
rrate_matHP_br0(1:159,1)=Y(1:159,1)-X(1:159,1)-hpfilter(Y(1:159,1)-X(1:159,1),HPlambda,0);

gy_matHP_br(1:130,1)=gy_matHP_br0(1:130,1);
gap_matHP_br(1:130,1)=gap_matHP_br0(1:130,1);
rrate_matHP_br(1:130,1)=rrate_matHP_br0(1:130,1);

% 1979-1982
gy_matHP_br0(110:179,1)=X(110:179,2)-hpfilter(X(110:179,2),HPlambda,0);
gap_matHP_br0(110:179,1)=X(110:179,3)-hpfilter(X(110:179,3),HPlambda,0);
rrate_matHP_br0(110:179,1)=Y(110:179,1)-X(110:179,1)-hpfilter(Y(110:179,1)-X(110:179,1),HPlambda,0);

gy_matHP_br(131:159,1)=gy_matHP_br0(131:159,1);
gap_matHP_br(131:159,1)=gap_matHP_br0(131:159,1);
rrate_matHP_br(131:159,1)=rrate_matHP_br0(131:159,1);


% post 1982
gy_matHP_br(160:317,1)=X(160:317,2)-hpfilter(X(160:317,2),HPlambda,0);
gap_matHP_br(160:317,1)=X(160:317,3)-hpfilter(X(160:317,3),HPlambda,0);
rrate_matHP_br(160:317,1)=Y(160:317,1)-X(160:317,1)-hpfilter(Y(160:317,1)-X(160:317,1),HPlambda,0);


% % two breaks
% % pre-1982
% rrate_matHP_br(1:159,1)=Y(1:159,1)-X(1:159,1)-hpfilter(Y(1:159,1)-X(1:159,1),HPlambda,0);
% gy_matHP_br(1:159,1)=X(1:159,2)-hpfilter(X(1:159,2),HPlambda,0);
% gap_matHP_br(1:159,1)=X(1:159,3)-hpfilter(X(1:159,3),HPlambda,0);
%
% % post 1982
% rrate_matHP_br(160:317,1)=Y(160:317,1)-X(160:317,1)-hpfilter(Y(160:317,1)-X(160:317,1),HPlambda,0);
% gy_matHP_br(160:317,1)=X(160:317,2)-hpfilter(X(160:317,2),HPlambda,0);
% gap_matHP_br(160:317,1)=X(160:317,3)-hpfilter(X(160:317,3),HPlambda,0);



%========================================================================
% Step 4: run the Kalman filter and the smoother
%========================================================================

disp('Step 4: Run Kalman filter and smoother ...')
[SEstate,SVstate,KFEstate,KFVstate,KFresid,MLE]=KalmanFSbreaks(Y,X,beta1,SigmaBB1,R1,SigmaBB2,R2,SigmaBB3,R3);
% compute standard errors

for i=1:T
    SVstate_se(i,:)=diag(sqrt(squeeze(SVstate(i,:,:))))';
end

for i=1:T
    KFVstate_se(i,:)=diag(sqrt(squeeze(KFVstate(i,:,:))))';
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%compute long-run responses and associated standard errors (use delta
%method)

for i=1:T

    % point estimates for the smoother
    SEstate_LR(i,1)=SEstate(i,1)/(1-SEstate(i,4)-SEstate(i,5));
    SEstate_LR(i,2)=SEstate(i,2)/(1-SEstate(i,4)-SEstate(i,5));
    SEstate_LR(i,3)=SEstate(i,3)/(1-SEstate(i,4)-SEstate(i,5));
    SEstate_LR(i,4)=SEstate(i,4);
    SEstate_LR(i,5)=SEstate(i,5);
    SEstate_LR(i,6)=SEstate(i,6)/(1-SEstate(i,4)-SEstate(i,5));

    Deriv0=[
        1/(1-SEstate(i,4)-SEstate(i,5))  0 0 SEstate(i,1)/(1-SEstate(i,4)-SEstate(i,5))^2 SEstate(i,1)/(1-SEstate(i,4)-SEstate(i,5))^2 0;
        0 1/(1-SEstate(i,4)-SEstate(i,5))  0 SEstate(i,2)/(1-SEstate(i,4)-SEstate(i,5))^2 SEstate(i,2)/(1-SEstate(i,4)-SEstate(i,5))^2 0;
        0  0 1/(1-SEstate(i,4)-SEstate(i,5)) SEstate(i,3)/(1-SEstate(i,4)-SEstate(i,5))^2 SEstate(i,3)/(1-SEstate(i,4)-SEstate(i,5))^2 0;
        0  0 0                               1                                            0                                            0;
        0  0 0                               0                                            1                                            0;
        0  0 0                               SEstate(i,6)/(1-SEstate(i,4)-SEstate(i,5))^2 SEstate(i,6)/(1-SEstate(i,4)-SEstate(i,5))^2 1/(1-SEstate(i,4)-SEstate(i,5));
        ];

    SVstate_LR(i,:,:) =  Deriv0*squeeze(SVstate(i,:,:))*Deriv0';
    SVstate_LR_se(i,:)=diag(sqrt(squeeze(SVstate_LR(i,:,:))))';

end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% compute by simulation the covariance matrix of rho1 and rho2 at quarterly
% frequency (since we do the conversion from monthly and 3/2 monthly
% frequncy into quarterly frequency by simulation)

disp('   convert rhos into quarterly frequency ...')
for i=1:T
    if mod(i,50)==0
        disp(i);
    end

    % Smoother
    param0=SEstate_LR(i,:)';
    Vdraw=squeeze(SVstate_LR(i,:,:));

    store_draws=[];
    for convi=1:num_sim_cov
        rd=randn(size(param0));
        betadraw=param0+(Vdraw^.5)*rd;

        if i<144
            rho_con=averaging([betadraw(4) betadraw(5)],3,2000);
        else
            rho_con=averaging([betadraw(4) betadraw(5)],2,2000);
        end

        betadraw_mod=betadraw;
        betadraw_mod(4)=rho_con(1);
        betadraw_mod(5)=rho_con(2);
        betadraw_mod(6)=betadraw(6)*(1-rho_con(1)-rho_con(2));

        store_draws=[store_draws; betadraw_mod'];
    end

    store_drawsM=nanmean(store_draws);

    SEstate_LR_M(i,:)=SEstate_LR(i,:);
    SEstate_LR_M(i,4)=store_drawsM(4);
    SEstate_LR_M(i,5)=store_drawsM(5);
    SEstate_LR_M(i,6)=param0(6)*(1-store_drawsM(4)-store_drawsM(5));

    SVstate_LR_M(i,:,:)=nancov(store_draws);
    SVstate_LR_M_se(i,:)=diag(sqrt(cov(store_draws)))';

end

% construct sum of autoregressive coefs
sumrho=SEstate_LR_M(:,4)+SEstate_LR_M(:,5);
sumrho_se=[];
for i=1:T
    Deriv0=[1 1];
    MM= Deriv0*squeeze(SVstate_LR_M(i,4:5,4:5))*Deriv0';
    sumrho_se(i,1)=sqrt(squeeze(MM));
end

%========================================================================
% Step 5: Compute and plot implied trend inflation
%========================================================================
Denom=  1-SEstate_LR_M(:,1).*(1-SEstate_LR_M(:,4)-SEstate_LR_M(:,5))...
    -SEstate_LR_M(:,4)-SEstate_LR_M(:,5);

interceptM0=SEstate_LR_M(:,6);

trendpi_br(:,1)=( interceptM0 - rrate_matHP_br(:,1).*(1-SEstate_LR_M(:,4)-SEstate_LR_M(:,5)) ...
    + gy_matHP_br(:,1).*SEstate_LR_M(:,2).*(1-SEstate_LR_M(:,4)-SEstate_LR_M(:,5)) ...
    + gap_matHP_br(:,1).*SEstate_LR_M(:,3).*(1-SEstate_LR_M(:,4)-SEstate_LR_M(:,5))) ...
    ./ Denom;

trendpi(:,1)=( interceptM0 - rrate_matHP(:,1).*(1-SEstate_LR_M(:,4)-SEstate_LR_M(:,5)) ...
    + gy_matHP(:,1).*SEstate_LR_M(:,2).*(1-SEstate_LR_M(:,4)-SEstate_LR_M(:,5)) ...
    + gap_matHP(:,1).*SEstate_LR_M(:,3).*(1-SEstate_LR_M(:,4)-SEstate_LR_M(:,5))) ...
    ./ Denom;


trendpi_br_se0=[];
trendpi_se0=[];

% compute standard errors for trend inflation using delta method
for i=1:T

    % point estimates for the smoother
    Denom000=1-SEstate_LR_M(i,1).*(1-SEstate_LR_M(i,4)-SEstate_LR_M(i,5))-SEstate_LR_M(i,4)-SEstate_LR_M(i,5);
    Denom000=1-SEstate_LR_M(i,1).*(1-SEstate_LR_M(i,4)-SEstate_LR_M(i,5))-SEstate_LR_M(i,4)-SEstate_LR_M(i,5);

    Numerator_br=( interceptM0(i,1) - rrate_matHP_br(i,1).*(1-SEstate_LR_M(i,4)-SEstate_LR_M(i,5)) ...
        + gy_matHP_br(i,1).*SEstate_LR_M(i,2).*(1-SEstate_LR_M(i,4)-SEstate_LR_M(i,5)) ...
        + gap_matHP_br(i,1).*SEstate_LR_M(i,3).*(1-SEstate_LR_M(i,4)-SEstate_LR_M(i,5)));

    Numerator=( interceptM0(i,1) - rrate_matHP(i,1).*(1-SEstate_LR_M(i,4)-SEstate_LR_M(i,5)) ...
        + gy_matHP(i,1).*SEstate_LR_M(i,2).*(1-SEstate_LR_M(i,4)-SEstate_LR_M(i,5)) ...
        + gap_matHP(i,1).*SEstate_LR_M(i,3).*(1-SEstate_LR_M(i,4)-SEstate_LR_M(i,5)));


    Deriv0=[
        Numerator_br*(1-SEstate_LR_M(i,4)-SEstate_LR_M(i,5))/Denom000^2;        % derivative wrt phipi
        gy_matHP_br(i,1)*(1-SEstate_LR_M(i,4)-SEstate_LR_M(i,5))/Denom000;      % derivative wrt phigy
        gap_matHP_br(i,1)*(1-SEstate_LR_M(i,4)-SEstate_LR_M(i,5))/Denom000;     % derivative wrt phix
        interceptM0(i,1)*(1-SEstate_LR_M(i,1))/Denom000^2;                      % derivative wrt rho1
        interceptM0(i,1)*(1-SEstate_LR_M(i,1))/Denom000^2;                      % derivative wrt rho2
        1/Denom000;                                                             % derivative wrt intercept
        ];

    MM =  Deriv0'*squeeze(SVstate_LR_M(i,:,:))*Deriv0;
    trendpi_br_se0(i,1)=diag(sqrt(squeeze(MM)))';


    Deriv0=[
        Numerator*(1-SEstate_LR_M(i,4)-SEstate_LR_M(i,5))/Denom000^2;         % derivative wrt phipi
        gy_matHP(i,1)*(1-SEstate_LR_M(i,4)-SEstate_LR_M(i,5))/Denom000;     % derivative wrt phigy
        gap_matHP(i,1)*(1-SEstate_LR_M(i,4)-SEstate_LR_M(i,5))/Denom000;    % derivative wrt phix
        interceptM0(i,1)*(1-SEstate_LR_M(i,1))/Denom000^2;                      % derivative wrt rho1
        interceptM0(i,1)*(1-SEstate_LR_M(i,1))/Denom000^2;                      % derivative wrt rho2
        1/Denom000;                         % derivative wrt intercept
        ];

    MM =  Deriv0'*squeeze(SVstate_LR_M(i,:,:))*Deriv0;
    trendpi_se0(i,1)=diag(sqrt(squeeze(MM)))';


end



%========================================================================

%   FIGURE 3 in PAPER: SMOOTHED LR RESPONSES AND TREND INFLATION
figure('name','Smoother point long-run response estimates and confidence bands')
for i=1:3
    subplot(3,2,i)
    plot(1:T,SEstate_LR_M(:,i),'blue-','Linewidth',2);                      % here, plot responses to inflation, output growth, and output gap
    hold on
    plot(1:T,SEstate_LR_M(:,i)+SVstate_LR_M_se(:,i),'black:','Linewidth',1);
    plot(1:T,SEstate_LR_M(:,i)-SVstate_LR_M_se(:,i),'black:','Linewidth',1);
    hold off
    if i==1        ylim([0.5 3.5])
    elseif i==2    ylim([-0.5 2.5])
    elseif i==3    ylim([0 1])
    end
    set(gca,'XTick',labelfigs(:,1))
    set(gca,'XTickLabel',labelfigs(:,2))
    xlim([1 T])
end

subplot(3,2,4)
plot(1:T,sumrho(:,1),'blue-','Linewidth',2);                                % here, plot sum of interest smoothing terms
hold on
plot(1:T,sumrho(:,1)+sumrho_se(:,1),'black:','Linewidth',1);
plot(1:T,sumrho(:,1)-sumrho_se(:,1),'black:','Linewidth',1);
hold off
set(gca,'XTick',labelfigs(:,1))
set(gca,'XTickLabel',labelfigs(:,2))
xlim([1 T])
ylim([0.3 1])

subplot(3,2,5)
plot(1:T,SEstate_LR_M(:,6),'blue-','Linewidth',2);                          % here, plot time-varying constant
hold on
plot(1:T,SEstate_LR_M(:,6)+SVstate_LR_M_se(:,6),'black:','Linewidth',1);
plot(1:T,SEstate_LR_M(:,6)-SVstate_LR_M_se(:,6),'black:','Linewidth',1);
hold off
set(gca,'XTick',labelfigs(:,1))
set(gca,'XTickLabel',labelfigs(:,2))
xlim([1 T])
ylim([-2 2])

subplot(3,2,6)
plot(1:T,movmean(trendpi_br(:,1),5),'blue-','Linewidth',2);                  % here, plot time-varying trend inflation (with break)
hold on
plot(1:T,movmean(trendpi_br(:,1),5)+movmean(trendpi_br_se0(:,1),5),'black:','Linewidth',1);
plot(1:T,movmean(trendpi_br(:,1),5)-movmean(trendpi_br_se0(:,1),5),'black:','Linewidth',1);
hold off
%legend('no break trend inflation','actual average inflation rate','90% C.I.')
set(gca,'XTick',labelfigs(:,1))
set(gca,'XTickLabel',labelfigs(:,2))
xlim([1 T])


subplot(3,2,1)
title('\phi_{\pi,t}')

subplot(3,2,2)
title('\phi_{gy,t}')

subplot(3,2,3)
title('\phi_{x,t}')

subplot(3,2,4)
title('\rho_{1,t} + \rho_{2,t}')

subplot(3,2,5)
title('c_t')

subplot(3,2,6)
title('Trend Inflation')



