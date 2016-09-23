% Computes market clearing capital stock and associated distribution and 
% decision rules in steady state
%
% Thomas Winberry, July 26th, 2016

%%
%----------------------------------------------------------------
% Compute approximation tools
%----------------------------------------------------------------

% Grids
computeGrids;

% Polynomials over grids (only if using polynomials to approximate conditional expectation)
if splineOpt == 0
	computePolynomials;
end

%%
%----------------------------------------------------------------
% Compute initial guess of market-clearing capital stock using
% histogram approximation of distribution, from Young (2010)
%----------------------------------------------------------------

t0 = tic;
fprintf('Computing initial guess from histogram...\n')

% Solve for equilibrium prices
x_0 = [0.05 0.12] ;
options = optimoptions('fsolve','Display',displayOpt,'TolFun',1e-4); % In older versions of MATLAB, use: options = optimset('Display',displayOpt); 
[x,err,exitflag] = fsolve(@(x) find_SS_histogram(x),x_0,options);

% Return exitflag if market clearing not solved
if exitflag < 1
    check = 1;
    return; 
end

r_ss = x(1) ;
p_ss = x(2) ;

if strcmp(displayOpt,'iter-detailed') == 1
    fprintf('Done! Time to compute: %2.2f seconds \n\n',toc(t0))
end

%%
%----------------------------------------------------------------
% Compute moments of histogram to use as initial guess for parametric family
%----------------------------------------------------------------

% Compute histogram

[~,mHistogram,mAssetsPrime_Histogram,mConsumption_Histogram] = find_SS_histogram(x) ;
    figure
        subplot(2,2,1); plot(vAssetsGridFine,mHistogram(:,1,1) / sum(mHistogram(:,1,1))); title('Dist: e(low) z(low)')
        subplot(2,2,2); plot(vAssetsGridFine,mHistogram(:,1,2) / sum(mHistogram(:,1,2))); title('Dist: e(low) z(high)')
        subplot(2,2,3); plot(vAssetsGridFine,mHistogram(:,2,1) / sum(mHistogram(:,2,1))); title('Dist: e(high) z(low)')
        subplot(2,2,4); plot(vAssetsGridFine,mHistogram(:,2,2) / sum(mHistogram(:,2,2))); title('Dist: e(high) z(high)')
    figure; 
        plot(vAssetsGridFine,sum(sum(mHistogram,3),2) ); title('Dist: All')
    figure; 
        plot(vAssetsGridFine(1:50),sum(sum(mHistogram(1:50,:,:),3),2) ) ; title('Dist: First 50 nodes')
    figure; 
        plot(log(vAssetsGridFine(1:end-1)),log(1-cumsum(sum(sum(mHistogram(1:end-1,:,:),3),2))) ) ; title('Dist: Pareto Tail')
    figure
        subplot(2,2,1); plot(vAssetsGridFine,mAssetsPrime_Histogram(:,1,1) ,vAssetsGridFine,vAssetsGridFine); title('Ap: e(low) z(low)')
        subplot(2,2,2); plot(vAssetsGridFine,mAssetsPrime_Histogram(:,1,2) ,vAssetsGridFine,vAssetsGridFine); title('Ap: e(low) z(high)')
        subplot(2,2,3); plot(vAssetsGridFine,mAssetsPrime_Histogram(:,2,1) ,vAssetsGridFine,vAssetsGridFine); title('Ap: e(high) z(low)')
        subplot(2,2,4); plot(vAssetsGridFine,mAssetsPrime_Histogram(:,2,2) ,vAssetsGridFine,vAssetsGridFine); title('Ap: e(high) z(high)')

% Compute moments from histogram
mMomentsHistogram = zeros(nEpsilon,nZ,nMeasure);
aGridMoments      = zeros(nAssetsQuadrature,nEpsilon,nZ,nMeasure); % grid for computing PDF

for iZ       = 1 : nZ
for iEpsilon = 1 : nEpsilon
	
	% First moment (uncentered)
	mMomentsHistogram(iEpsilon,iZ,1) = vAssetsGridFine' * mHistogram(:,iEpsilon,iZ) / sum(mHistogram(:,iEpsilon,iZ)) ;
	aGridMoments(:,iEpsilon,iZ,1)    = vAssetsGridQuadrature - mMomentsHistogram(iEpsilon,iZ,1);
		
	% Higher order moments (centered)
	for iMoment = 2 : nMeasure
		mMomentsHistogram(iEpsilon,iZ,iMoment) = ((vAssetsGridFine' - mMomentsHistogram(iEpsilon,iZ,1)) .^ iMoment) * ...
			(mHistogram(:,iEpsilon,iZ) / sum(mHistogram(:,iEpsilon,iZ)));
		aGridMoments(:,iEpsilon,iZ,iMoment) = (vAssetsGridQuadrature - mMomentsHistogram(iEpsilon,iZ,1)) .^ iMoment - ...
            mMomentsHistogram(iEpsilon,iZ,iMoment);
	end	
	
end
end


% Mass at borrowing constraint
mHatHistogram = mHistogram(1,:,:)./sum(mHistogram,1) ; 
mHatHistogram = squeeze(mHatHistogram) ;

%%
%----------------------------------------------------------------
% Compute market-clearing capital stock from parametric family
%----------------------------------------------------------------

t0 = tic; 
fprintf('Compute steady state from parametric family...\n')

% Solve for equilibrium prices
    x_0 = [r_ss p_ss] ;
    f = @(x) find_SS_parametric(x,mMomentsHistogram,mHatHistogram,aGridMoments);
    options = optimoptions('fsolve','Display',displayOpt,'TolFun',1e-4); % In older versions of MATLAB, use: options = optimset('Display',displayOpt); 
    % Run fsolve only if initial guess is not good
    if abs(f(x_0)) > 1e-4
        [x,err,exitflag] = fsolve(@(x) find_SS_histogram(x),x_0,options);
    end 

% Return exitflag if market clearing not solved
    if exitflag < 1
        check = 1;
        return; 
    end
    if abs(f(x)) > 1e-4
        [aggregateCapital,err,exitflag] = fsolve(f,aggregateCapitalInit,options);
    end

r_ss = x(1) ;
p_ss = x(2) ;

if strcmp(displayOpt,'iter-detailed') == 1
    fprintf('Done! Time to compute: %2.2f seconds \n\n',toc(t0))
end

%%
%----------------------------------------------------------------
% Compute other objects from steady state
%----------------------------------------------------------------

[~,mCoefficients,mParameters,mMoments,mHat] = ...
    find_SS_parametric(x,mMomentsHistogram,aGridMoments,mHatHistogram);