%% Data analysis script to accompany Torres et al. 2018 Nature Geoscience
% "Riverine evidence for isotopic mass balance in the Earth's early sulfur cycle"
% The script is provided "as-is" to re-produce portion of the data
% analysis. It is not guaranteed to work for other purposes. Line 104 requires 
% the wmean script from the Matlab file exchange. Without it, comment out
% lines 104 and 127 and uncomment line 128

clear all

%% load Data
data = readtable('canadaRiverData.csv'); %river data
load('roadSalt.mat'); %analyses of road salt from Price et al. 
load('an_precip_chem.mat'); %annual rainfall concentrations from NAtChem
load('bekker2009data.mat'); %rock D33S data from Ref. 19

%% MonteCarlo parameters
n = 1E3; %number of random draws. 1E7 used in paper. 

%% Endmembers
ClRainDist = makedist('normal',mean(an_precip_chem(:,1)),...
    std(an_precip_chem(:,1))); %Cl concentration of rain
SRainDist = makedist('normal',mean(an_precip_chem(:,2)),...
    std(an_precip_chem(:,2))); %SO4 concentration of rain
D33SrainDist = makedist('uniform','lower',-0.1,'upper',0.1); % D33 of rain
SClsaltDist = makedist('uniform','lower',0.001,'upper',0.01); %S/Cl in salt

%% Mixing Analysis
%Pre-allocation
fSsalt = zeros(length(data.Cl),2); %S derived from urban runoff
fSrock = zeros(length(data.Cl),3); %S derived from crust
fSrain = zeros(length(data.Cl),2); %S derived from rain
fSnonrock = zeros(length(data.Cl),2); %non-crustal S (rain + urban runoff)
D33Sac = NaN(length(data.D33S),n); %D33S of Archean Crust (AC)
D33SacStat = NaN(length(D33Sac(:,1)),2); %mean and stdev of D33S_AC

%Generate random values from end-members to perform mixing
SClsalt = random(SClsaltDist,[n,1]); %S/Cl in salt
ClRain = random(ClRainDist,[n,1]); %Cl concentration of rain
SRain = random(SRainDist,[n,1]); %SO4 concentration of rain
D33Srain = random(D33SrainDist,[n,1]); % D33 of rain

for m = 1:length(data.Cl) % for each data point
    %% Analytical Uncertainty
    %10% uncertainty on measured riverine Cl
    ClRivDist = makedist('normal',data.Cl(m),data.Cl(m).*0.1);
    %10% uncertainty on measured riverine SO4
    SRivDist = makedist('normal',data.SO4(m),data.SO4(m).*0.1);
    %0.2 permil uncertainty on measured riverine D33S
    D33StDist = makedist('normal',data.D33S(m),0.2);
    %Generate Random values within measurement uncertainties
    SRiv = random(SRivDist,[n,1]); %random river S concentrations 
    ClRiv = random(ClRivDist,[n,1]); %random river Cl concentrations 
    D33St = random(D33StDist,[n,1]); %random river D33S
    %% Calculate Cl budget
    ClSalt = ClRiv - ClRain; %Cl from roadsalt
    ClInd = find(ClSalt < 0); %find negative numbers
    ClSalt(ClInd) = NaN; %replace negative with NaN
    %% Calculate S budget
    SSalt = ClSalt.*SClsalt; %S dervied from salt from Cl and S/Cl ratio
    Srock = (SRiv - SRain - SSalt); %rock derived S as difference in concentrations.
    %% Determine Mixing Fractions
    f1 = Srock./SRiv; %proportion of S from rock
    f2 = SRain./SRiv; %proportion of S from rain
    f3 = SSalt./SRiv; %proportion of S from urban runoff
    sumTest = abs(f1)+abs(f2)+abs(f3); %are all fractions between 0 & 1?
    SInd = find(sumTest > 1 ); %find non-physical results (f>1 or f<0)
    Srock(SInd) = NaN; %replacing non-physical results with NaN
    SRain(SInd) = NaN;
    SSalt(SInd) = NaN;
    %% Save Mean and Standard Deviation
    fSrock(m,1) = mean( Srock./SRiv,'omitnan'); %mean of fraction rock
    fSrock(m,2) = std( Srock./SRiv,'omitnan'); %stdev of fraction rock
    fSrock(m,3) = min( Srock./SRiv); %minimum fraction rock
    fSsalt(m,1) = mean( SSalt./SRiv,'omitnan'); %mean of fraction urban runoff
    fSsalt(m,2) = std( SSalt./SRiv,'omitnan'); %stdev of fraction urban runoff
    fSrain(m,1) = mean( SRain./SRiv,'omitnan'); %mean of fraction rain
    fSrain(m,2) = std( SRain./SRiv,'omitnan'); %stdev of fraction rain
    fSnonrock(m,1) = mean((SSalt+SRain)./SRiv,'omitnan'); %mean of fraction non-rock
    fSnonrock(m,2) = std((SSalt+SRain)./SRiv,'omitnan'); %stdev of fraction non-rock
    %% D33S Mixing (to calculate D33S archean crust (AC)
    %D33S_AC for all fractions non-rock
    D33Sac(m,:) = (D33St - ((1-f1).*D33Srain))./f1; 
    D33SacStat(m,1) = mean(D33Sac(m,:),'omitnan'); %mean value
    D33SacStat(m,2) =  std(D33Sac(m,:),'omitnan'); %standard deviation
end

%% Area weighting
%sort out samples with with high stdev for D33S_AC
stdInd = find(D33SacStat(:,2) < 1); %stdev less than 1 permil
Area = data.Area_km2_(stdInd); %find catchments areas for the "good" data
D33SacStatHigh = D33SacStat(stdInd,:); % making matrix of good means and stdev
[~,sortInd] = sort(D33SacStatHigh(:,1)); %sorting by mean values
aWeights = repmat(Area,[1,n]); %area weighting
D33Sac_ = D33Sac(stdInd,:); %matrix of good D33S results
D33SacRS = reshape(D33Sac_,[],1); %reshape to one column
aWeights = reshape(aWeights,[],1); %reshape to one column

xInd = find(isfinite(D33SacRS)==1); %index of finite mixing results
D33SacRS = D33SacRS(xInd); %remove NaN values
aWeights = aWeights(xInd); %remove NaN values
aInd = find(isfinite(aWeights)==1);%index of finite areas
aWeights = aWeights(aInd); %remove NaN values
D33SacRS = D33SacRS(aInd); %remove NaN values
AreaMean = wmean(D33SacRS,aWeights); %weighted mean of average values
AreaStd = std(D33SacRS,aWeights); %weighted mean of standard deviations

Aind = find(isfinite(Area)==1); %index of finite areas
D33SacStat_ = D33SacStatHigh(Aind,:); %D33S data with finite areas
riverDist = 1; %pre-allocation
Area_ = Area(Aind); %area data

for m = 1:length(D33SacStat_(:,1)) %for each sample
    %create normal distribution from calculated mean and standard deviation
    Rdist = makedist('normal','mu',D33SacStat_(m,1),'sigma',D33SacStat_(m,2));
    Rdraws = random(Rdist,[10*Area_(m),1]); %re-sample proportional to area
    riverDist = [riverDist;Rdraws]; %append to matrix
end
%fit area-weighted distrbution to t distribution
riverDistFit = fitdist(riverDist,'tLocationScale'); 

%% Plotting D33S AC estimate
figure(1)
clf
hold on
yyaxis right
%x-values for +/- 4 standard deviations
xrange = linspace(AreaMean-4*AreaStd,AreaMean+4*AreaStd);
%xrange = linspace(-2,2);
%plot area-weighted D33S_AC distribution
plot(xrange,pdf(riverDistFit,xrange),'k-','linewidth',1.5) 
axis([-5 6 0 1.3])
yyaxis left
yrange = 1:length(D33SacStatHigh(:,1)); %range of y values
rectangle('position',[1 0 4 1+max(yrange)],...
    'FaceColor',[0.8 0.5 0.5],'EdgeColor',[0.8 0.5 0.5]) % Ref.7 estimate 
plot(0.5.*ones(100,1),linspace(0,36),'r-','linewidth',3) %Ref. 6 estimate
%plot individual river values
pts = ploterr(D33SacStatHigh(sortInd,1),yrange,D33SacStatHigh(sortInd,2),0,'ko','hhxy',0);
set(pts,'markerfacecolor','w','linewidth',1.5,'markersize',12)
axis([-5 8 0 1+max(yrange)])
xlabel(['\Delta^{33}S (',char(8240),')']); ylabel('p(X)')

%% Plotting raw D33S data
figure(2)
clf
subplot(1,2,1)
hold on
dataFit = fit(fSrock(:,1),data.D33S,'poly1'); %linear fit to data
dataFitci = predint(dataFit,linspace(0,1),0.68); %68% confidence interval of fit
plot(feval(dataFit,linspace(0,1)),linspace(0,1),'r-','linewidth',2) %plot best fit
plot(dataFitci(:,1),linspace(0,1),'r--','linewidth',1.5) %plot bound
plot(dataFitci(:,2),linspace(0,1),'r--','linewidth',1.5) %plot bound
pts = ploterr(data.D33S,fSrock(:,1),0.2,fSrock(:,2),'ko','hhxy',0); %plot data
set(pts,'markerfacecolor','w','linewidth',1.5,'markersize',16,...
    'displayname','Canadian Rivers (\pm 1\sigma)')
%plot Ref. 19 data
rocks = plot(bek(:,2),ones(length(bek(:,1)),1),'ws','displayname','Local Sulfides (Ref. 19)');
set(rocks,'markerfacecolor',[0.5 0.5 0.5],'markersize',18)
%plot predicted atmospheric end-member
plot(linspace(-0.1,0.1),zeros(100,1),'k-','linewidth',8,'displayname','Rain / Modern Contamination')
set(gca,'ydir','reverse')
xlabel(['\Delta^{33}S (',char(8240),')']); ylabel('Fraction Crustal Sulfur (f_C)')
legend('location','northeast')

subplot(1,2,2)
hold on
xPro = 1-(data.x_ProI + data.x_ProS); %calculate fraction proterozoic area
pInd = find(isfinite(xPro)==1); %remove data without geologic data
dataFit = fit(xPro(pInd),data.D33S(pInd),'poly1'); %fit data
dataFitci = predint(dataFit,linspace(0,1),0.68); %calculate confidence interval
plot(feval(dataFit,linspace(0,1)),linspace(0,1),'r-','linewidth',2) %plot best fit
plot(dataFitci(:,1),linspace(0,1),'r--','linewidth',1.5) %plot bound
plot(dataFitci(:,2),linspace(0,1),'r--','linewidth',1.5) %plot bound
pts = ploterr(data.D33S,xPro,0.2,0,'ko','hhxy',0); %plot data
set(pts,'markerfacecolor','w','linewidth',1.5,'markersize',16,'displayname','Canadian Rivers')
xlabel(['\Delta^{33}S (',char(8240),')']);ylabel('Fraction Archean Catchment Area')
set(gca,'ydir','reverse')
axis([-0.6 0.6 0 1])

%% Plotting raw d34S data
figure(3)
clf
hold on
%fit fraction of urban runoff vs. d34S
dataFit = fit(fSsalt(:,1),data.d34S,'poly1'); 
dataFitci = predint(dataFit,linspace(0,1),0.68); %calculate confidence interval
plot(feval(dataFit,linspace(0,1)),linspace(0,1),'r-','linewidth',2) %plot best fit
plot(dataFitci(:,1),linspace(0,1),'r--','linewidth',1.5) %plot bound
plot(dataFitci(:,2),linspace(0,1),'r--','linewidth',1.5) %plot bound
pts = ploterr(data.d34S,fSsalt(:,1),0.1,fSsalt(:,2),'ko','hhxy',0); %plot data
set(pts,'markerfacecolor','w','linewidth',1.5,'markersize',16,...
    'displayname','Canadian Rivers (\pm 1\sigma)')

plot(linspace(2,9),zeros(100,1),'k-','linewidth',8) %plot range of d34S rain
%plot range of d34S crust
plot(linspace(min(bek(:,1)),max(bek(:,1))),zeros(100,1),'b-','linewidth',8) 

axis([-1 25 0 0.25])
xlabel(['\delta^{34}S (',char(8240),')']); ylabel('Fraction Urban Runoff')

% Make inset to project fit to 100% urban runoff
axes('Position',[.6 .5 .35 .45])
hold on
yyaxis right
plot(feval(dataFit,linspace(0,1)),linspace(0,1),'r-','linewidth',2)
plot(dataFitci(:,1),linspace(0,1),'r--','linewidth',1.5)
plot(dataFitci(:,2),linspace(0,1),'r--','linewidth',1.5)
pts = ploterr(data.d34S,fSsalt(:,1),0.1,fSsalt(:,2),'ko','hhxy',0);
set(pts,'markerfacecolor','w','linewidth',1.5,'markersize',16,...
    'displayname','Canadian Rivers (\pm 1\sigma)')
box on
ylabel(['\delta^{34}S (',char(8240),')']); ylabel('Fraction Urban Runoff')
