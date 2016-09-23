% Shell which declares parameters and calls Dynare to solve model using
% first order approximation of aggregate dynamics
%
% Thomas Winberry, July 26th, 2016

clear all
close all
clc

%----------------------------------------------------------------
% Solve for Steady State
%----------------------------------------------------------------
steadyState;

    mAssetsPrime_SS = mAssetsPrime ;
    mConsumption_SS = mConsumption ;
    mHistogram_SS   = mHistogram   ;
    mDistributionFine_SS = mDistributionFine ;
    
save Steady_State.mat mAssetsPrime_SS mConsumption_SS mHistogram_SS mDistributionFine_SS

    


oldFolder = cd('./Auxiliary Functions');

%----------------------------------------------------------------
% Set parameters
%----------------------------------------------------------------
setParameters;

%----------------------------------------------------------------
% Compute approximation tools
%----------------------------------------------------------------

% Grids
computeGrids;

% Polynomials over grids (only if using polynomials to approximate conditional expectation)
if splineOpt == 0
	computePolynomials;
end

%----------------------------------------------------------------
% Save parameters in .mat files to import into Dynare 
%----------------------------------------------------------------

if splineOpt == 0	% if using polynomials to approximate individual decisions

	% Economic parameters
	save economicParameters.mat bbeta ssigma aaBar aalpha ddelta vEpsilonGrid aggEmployment ...
		mmu ttau rrhoTFP ssigmaTFP mEpsilonTransition
		
	% Approximation parameters
	save approximationParameters.mat nEpsilon nAssets nState assetsMin assetsMax nAssetsFine nStateFine nAssetsQuadrature nStateQuadrature ...
		nMeasure nMeasureCoefficients kRepSS maxIterations tolerance dampening
		
	% Grids
	save grids.mat vAssetsGridZeros vAssetsGrid mEpsilonGrid mAssetsGrid mEpsilonPrimeGrid vAssetsGridFine ...
		vAssetsGridFineZeros mEpsilonGridFine mAssetsGridFine mEpsilonPrimeGridFine vQuadratureWeights ...
		vAssetsGridQuadratureZeros vAssetsGridQuadrature mEpsilonGridQuadrature mAssetsGridQuadrature
		
	% Polynomials
	save polynomials.mat vAssetsPoly vAssetsPolySquared vAssetsPolyFine vAssetsPolyQuadrature vAssetsPolyBC
	
else	% if using splines to approximate individual decisions

	% Economic parameters
	save economicParameters.mat bbeta ssigma aaBar aalpha ddelta vEpsilonGrid aggEmployment ...
		mmu ttau rrhoTFP ssigmaTFP mEpsilonTransition
		
	% Approximation parameters
	save approximationParameters.mat nEpsilon nAssets nState assetsMin assetsMax nAssetsFine nStateFine nAssetsQuadrature nStateQuadrature ...
		nMeasure nMeasureCoefficients kRepSS maxIterations tolerance dampening
		
	% Grids
	save grids.mat vAssetsGrid mEpsilonGrid mAssetsGrid mEpsilonPrimeGrid vAssetsGridFine ...
		mEpsilonGridFine mAssetsGridFine mEpsilonPrimeGridFine vQuadratureWeights ...
		vAssetsGridQuadrature mEpsilonGridQuadrature mAssetsGridQuadrature
	
end

%----------------------------------------------------------------
% Run Dynare
%----------------------------------------------------------------

if splineOpt == 0	% if using polynomials to approximate individual decisions

	dynare firstOrderDynamics_polynomials
	
else	% if using splines to approximate individual decisions

	dynare firstOrderDynamics_splines
	
end

    save dynare_results.mat oo_ M_

%% 
%----------------------------------------------------------------
% Get dynamics of distribution and individual variables
%----------------------------------------------------------------

    % Get parameters and grids
        setParameters;
        computeGrids;
        if splineOpt==0
        computePolynomials;
        end 
        load dynare_results.mat
        load Steady_State.mat

    % Number of periods
        T = numel(oo_.irfs.r_aggregateTFPShock);

    % Number of variables
        N_var = M_.endo_nbr ;
        
    % Variable names
        var_name = M_.endo_names ;
    
    % Steady state values
        SS = oo_.steady_state ;
        
    % IRF for all variables
        for i=1:N_var
            eval(['IRF.' , strtrim(var_name(i,:)) , '=' , 'oo_.irfs.' , strtrim(var_name(i,:)) , '_aggregateTFPShock + ' , num2str(SS(i)) , ';' ])
        end 
            
    % Get Distribution
        Dist_IRF = NaN(nEpsilon,nAssetsFine,T) ;
        options = optimoptions(@fminunc,'Algorithm','quasi-newton','Display','notify-detailed',...
                    'MaxFunEvals',50000,'TolFun',1e-12,'GradObj','on','MaxIter',1000);
        for t=1:T
            % Get Moments for that period
            moments      = NaN(nMeasure,1) ;
            vParameters  = NaN(nMeasure,1) ;
            aGridMoments = NaN(nAssetsQuadrature,nMeasure) ;
            aGridMoments_Fine = NaN(nAssetsFine,nMeasure) ;
            for i=1:nEpsilon 
                for j=1:nMeasure 
                   eval([ 'vParameters(',num2str(j),')=IRF.measureCoefficient_',num2str(i),'_',num2str(j),'(',num2str(t),');' ]) 
                   eval([ 'moments(',num2str(j),')=IRF.moment_',num2str(i),'_',num2str(j),'(',num2str(t),');' ]) 
                   

                   % Get moment-asset grid for computing distribution
                   if j==1
                       aGridMoments(:,j) = vAssetsGridQuadrature - moments(j);
                       aGridMoments_Fine(:,j) = vAssetsGridFine - moments(j);
                   else
                       aGridMoments(:,j) = (vAssetsGridQuadrature' - moments(1)) .^ j - moments(j);
                       aGridMoments_Fine(:,j) = (vAssetsGridFine' - moments(1)) .^ j - moments(j);
                   end	
                end
                
                normalization = vQuadratureWeights' * exp(aGridMoments * vParameters);
                mParameters = [1 / normalization; vParameters] ;
                eval([ 'm_hat=IRF.mHat_',num2str(i),'(',num2str(t),');' ]) 
                
                % Get distribution 
                Dist_IRF(i,:,t) = mParameters(1) * exp(aGridMoments_Fine * ...
                        mParameters(2:nMeasure+1));
                % Dist_IRF(i,1,t) = m_hat ;
            end        
            
        end 
        
       
    % Get Policy Functions
        Aprime_IRF = NaN(nEpsilon,nAssetsFine,T) ;
        Cons_IRF   = NaN(nEpsilon,nAssetsFine,T) ;
        for t=1:T
            for i=1:nEpsilon
                for j=1:nAssets
                    eval([ 'mCoefficients(',num2str(i),',',num2str(j),')=IRF.expectationCoefficient_',num2str(i),'_',num2str(j),'(',num2str(t),');' ])          
                end 
            end 
            
            % Get prices
            eval([ 'w=IRF.w(',num2str(t),');' ])  
            eval([ 'r=IRF.r(',num2str(t),');' ])  
            
            % Compute conditional expectation
            mConditionalExpectation = exp(mCoefficients * vAssetsPolyFine');
            
            % Compute savings policy
            mAssetsPrimeStar = w * (mmu * (1 - mEpsilonGridFine) + (1 - ttau) * mEpsilonGridFine) + ...
                (1 + r) * mAssetsGridFine - (mConditionalExpectation .^ (-1 / ssigma));
            Aprime_IRF(:,:,t) = max(mAssetsPrimeStar,aaBar * ones(nEpsilon,nAssetsFine));

            % Compute consumption
            Cons_IRF(:,:,t) = w * (mmu * (1 - mEpsilonGridFine) + (1 - ttau) * mEpsilonGridFine) + ...
                (1 + r) * mAssetsGridFine - Aprime_IRF(:,:,t);
          
        end 
        
%%        
    % Figure for distribution
        %%%%% CLEAN FIGURE TO GET READY TO RECORD%%%%%
        figure(1);
        clf;
        close 1;

        %%%%% Set Parameters%%%%%
        plotlength=T;  % End time of the video
        icut = nAssetsFine;      % Index of maximum wealth being plotted.

        h=figure('Renderer','zbuffer');
        %M(round(plotlength/dt)) = struct('cdata',[],'colormap',[]);
        axis tight manual
        %set(gca,'NextPlot','replaceChildren');
        set(gca,'FontSize',14)
        for k=1:plotlength
            g=Dist_IRF(:,:,k)' ; g = g./repmat(sum(g),nAssetsFine,1) ;
            acut = vAssetsGridFine(1:icut);
            gcut = g(1:icut,:);
            plot(acut,gcut'); 
            xlim([min(acut) max(acut)])
            xlabel('Wealth, a','FontSize',14)
            ylabel('Density','FontSize',14)
            title(sprintf('Evolution of Distribution'));
            legend('Unemployed','Employed','location','northeast');
            set(gcf,'color','w')
            M(k)=getframe(gcf);
        end
        movie2avi(M(1:plotlength),'distribution.avi','fps',4);

        
    % Figure for distribution - Difference
        %%%%% CLEAN FIGURE TO GET READY TO RECORD%%%%%
        figure(1);
        clf;
        close 1;

        %%%%% Set Parameters%%%%%
        plotlength=T;  % End time of the video
        icut = nAssetsFine;      % Index of maximum wealth being plotted.

        h=figure('Renderer','zbuffer');
        %M(round(plotlength/dt)) = struct('cdata',[],'colormap',[]);
        axis tight manual
        %set(gca,'NextPlot','replaceChildren');
        set(gca,'FontSize',14)
        for k=1:plotlength
            g=Dist_IRF(:,:,k)' ;  g = g./repmat(sum(g),nAssetsFine,1) ;
            g_ss = mDistributionFine_SS' ; g_ss = g_ss./repmat(sum(g_ss),nAssetsFine,1) ;
            acut = vAssetsGridFine(1:icut);
            gcut = g(1:icut,:); g_ss_cut = g_ss(1:icut,:) ; 
            plot(acut,gcut'-g_ss_cut',acut,zeros(size(acut)),'--k'); 
            xlim([min(acut) max(acut)])
            ylim([-0.001 0.001])
            xlabel('Wealth, a','FontSize',14)
            ylabel('Difference with SS','FontSize',14)
            title(sprintf('Evolution of Distribution wrt SS'));
            legend('Unemployed','Employed','location','northeast');
            set(gcf,'color','w')
            M(k)=getframe(gcf);
        end
        movie2avi(M(1:plotlength),'distribution_SS.avi','fps',4);
        
%%        
    % Figure for Savings
        %%%%% CLEAN FIGURE TO GET READY TO RECORD%%%%%
        figure(1);
        clf;
        close 1;

        %%%%% Set Parameters%%%%%
        plotlength=T;  % End time of the video
        icut = nAssetsFine;      % Index of maximum wealth being plotted.
        imin = 2 ; %index of minimun wealth being plotted

        h=figure('Renderer','zbuffer');
        %M(round(plotlength/dt)) = struct('cdata',[],'colormap',[]);
        axis tight manual
        %set(gca,'NextPlot','replaceChildren');
        set(gca,'FontSize',14)
        for k=1:plotlength
            g=Aprime_IRF(:,:,k)' ;
            acut = vAssetsGridFine(imin:icut);
            gcut = g(imin:icut,:);
            plot(acut,gcut'); 
            xlim([min(acut) max(acut)])
            xlabel('Wealth, a','FontSize',14)
            ylabel('Savings','FontSize',14)
            title(sprintf('Evolution of Saving Policy Function'));
            legend('Unemployed','Employed','location','southeast');
            set(gcf,'color','w')
            M(k)=getframe(gcf);
        end
        movie2avi(M(1:plotlength),'aprime.avi','fps',4);

        
    % Figure for savings - Difference
        %%%%% CLEAN FIGURE TO GET READY TO RECORD%%%%%
        figure(1);
        clf;
        close 1;

        %%%%% Set Parameters%%%%%
        plotlength=T;  % End time of the video
        icut = nAssetsFine;      % Index of maximum wealth being plotted.
        imin = 3 ; %index of minimun wealth being plotted

        h=figure('Renderer','zbuffer');
        %M(round(plotlength/dt)) = struct('cdata',[],'colormap',[]);
        axis tight manual
        %set(gca,'NextPlot','replaceChildren');
        set(gca,'FontSize',14)
        for k=1:plotlength
            g=Aprime_IRF(:,:,k)' ;
            g_ss = mAssetsPrime_SS' ;
            acut = vAssetsGridFine(imin:icut);
            gcut = g(imin:icut,:);
            g_ss_cut = g_ss(imin:icut,:) ;
            plot(acut,100*(log(gcut)'-log(g_ss_cut)'),acut,zeros(size(acut)),'--k'); 
            xlim([min(acut) max(acut)])
            ylim([-2.0 2.0])
            xlabel('Wealth, a','FontSize',14)
            ylabel('Difference with SS (%)','FontSize',14)
            title(sprintf('Evolution of Savings wrt SS'));
            legend('Unemployed','Employed','location','northeast');
            set(gcf,'color','w')
            M(k)=getframe(gcf);
        end
        movie2avi(M(1:plotlength),'aprime_SS.avi','fps',4);
        
        
%%        
    % Figure for Consumption
        %%%%% CLEAN FIGURE TO GET READY TO RECORD%%%%%
        figure(1);
        clf;
        close 1;

        %%%%% Set Parameters%%%%%
        plotlength=T;  % End time of the video
        icut = nAssetsFine;      % Index of maximum wealth being plotted.

        h=figure('Renderer','zbuffer');
        %M(round(plotlength/dt)) = struct('cdata',[],'colormap',[]);
        axis tight manual
        %set(gca,'NextPlot','replaceChildren');
        set(gca,'FontSize',14)
        for k=1:plotlength
            g=Cons_IRF(:,:,k)' ;
            acut = vAssetsGridFine(1:icut);
            gcut = g(1:icut,:);
            plot(acut,gcut'); 
            xlim([min(acut) max(acut)])
            xlabel('Wealth, a','FontSize',14)
            ylabel('Savings','FontSize',14)
            title(sprintf('Evolution of Consumption Policy Function'));
            legend('Unemployed','Employed','location','southeast');
            set(gcf,'color','w')
            M(k)=getframe(gcf);
        end
        movie2avi(M(1:plotlength),'consumption.avi','fps',4);

        
    % Figure for consumption - Difference
        %%%%% CLEAN FIGURE TO GET READY TO RECORD%%%%%
        figure(1);
        clf;
        close 1;

        %%%%% Set Parameters%%%%%
        plotlength=T;  % End time of the video
        icut = nAssetsFine;      % Index of maximum wealth being plotted.

        h=figure('Renderer','zbuffer');
        %M(round(plotlength/dt)) = struct('cdata',[],'colormap',[]);
        axis tight manual
        %set(gca,'NextPlot','replaceChildren');
        set(gca,'FontSize',14)
        for k=1:plotlength
            g=Cons_IRF(:,:,k)' ;
            acut = vAssetsGridFine(1:icut);
            gcut = g(1:icut,:);
            plot(acut,100*(log(gcut)'-log(mConsumption_SS)),acut,zeros(size(acut)),'--k'); 
            xlim([min(acut) max(acut)])
            ylim([-1.0 1.0])
            xlabel('Wealth, a','FontSize',14)
            ylabel('Difference with SS (%)','FontSize',14)
            title(sprintf('Evolution of Consumption wrt SS'));
            legend('Unemployed','Employed','location','northeast');
            set(gcf,'color','w')
            M(k)=getframe(gcf);
        end
        movie2avi(M(1:plotlength),'consumption_SS.avi','fps',4);
                
%% IRF for aggregate variables

    figure;
    subplot(2,3,1); plot(1:40,100*oo_.irfs.aggregateTFP_aggregateTFPShock,1:40,zeros(1,40),'--k');
    title('TFP'); ylabel('Diff. with SS (%)')
    subplot(2,3,2); plot(1:40,100*oo_.irfs.logAggregateOutput_aggregateTFPShock,1:40,zeros(1,40),'--k');
    title('Output');
    subplot(2,3,3); plot(1:40,100*oo_.irfs.logAggregateConsumption_aggregateTFPShock,1:40,zeros(1,40),'--k');
    title('Consumption');
    subplot(2,3,4); plot(1:40,100*oo_.irfs.logAggregateInvestment_aggregateTFPShock,1:40,zeros(1,40),'--k');
    title('Savings'); ylabel('Diff. with SS (%)')
    subplot(2,3,5); plot(1:40,100*oo_.irfs.logWage_aggregateTFPShock,1:40,zeros(1,40),'--k');
    title('Wage');
    subplot(2,3,6); plot(1:40,100*oo_.irfs.r_aggregateTFPShock,1:40,zeros(1,40),'--k');
    title('Interest Rate');
    set(gcf,'color','w')
    print('-depsc','IRF_Aggregate_Variables.eps')
    

