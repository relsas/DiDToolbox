%% Medicaid and Mortality Analysis

clear; close all;
workDir ="d:\Temp\MCP\DID\did-toolbox\";
cd(workDir);

if exist(workDir+"data\"+"Medicaid.mat")==false
    data = readtable(workDir+"Data\"+"Medicaid_own.xlsx");

    save(workDir+"Data\"+"Medicaid.mat","data");
else
    load(workDir+"Data\"+"Medicaid.mat","data");
end
%% Differences in Mean (2014)
fprintf('  ------- Differences in Mean (2014) ------\n')
data.G2014 = data.yaca==2014;
data.control2014 = isnan(data.yaca) | data.yaca>2019;
data.post2014 = data.year>=2014;
% unweighted
sample = data.year>=2013 & data.year<=2014 & (data.G2014|data.control2014);
g2014_noW = groupsummary(data(sample,:),["G2014","year"],"mean","crude_rate_20_64");

% weighted
wmean = @(w,v)mean(w'*v);
g2014_W = groupsummary(data(sample,:),["G2014","year"],wmean,{"population_20_64","crude_rate_20_64"});
sumPop = groupsummary(data(sample,:),["G2014","year"],"sum","population_20_64");
g2014_W(:,4) = g2014_W(:,4)./sumPop.sum_population_20_64;

% calc mean diffs
m2014_noW = [g2014_noW{3,4},g2014_noW{1,4};...
    g2014_noW{4,4},g2014_noW{2,4};...
    g2014_noW{4,4}-g2014_noW{3,4},g2014_noW{2,4}-g2014_noW{1,4}];

m2014_noW = [m2014_noW,[m2014_noW(1,1)-m2014_noW(1,2);...
    m2014_noW(2,1)-m2014_noW(2,2);...
    m2014_noW(3,1)-m2014_noW(3,2)]];


m2014_W = [g2014_W{3,4},g2014_W{1,4};...
    g2014_W{4,4},g2014_W{2,4};...
    g2014_W{4,4}-g2014_W{3,4},g2014_W{2,4}-g2014_W{1,4}];

m2014_W = [m2014_W,[m2014_W(1,1)-m2014_W(1,2);...
    m2014_W(2,1)-m2014_W(2,2);...
    m2014_W(3,1)-m2014_W(3,2)]];


res2014 = array2table([m2014_noW,m2014_W],VariableNames=["Expansion(noW)","NoExpansion(noW)","Gap_nW/DiD",...
    "Expansion(W)","NoExpansion(W)","Gap_W/DiD"],...
    RowNames = ["2013","2014","Trend"])
%% Simple Regression
% Regression - unweighted
fprintf('  ------- Regression unweighted (2014) ------\n')
% same with regression - unweigthed

treatPost2014 = data.G2014.*data.post2014;
data.treatPost2014 = treatPost2014;

data2014=data(sample,:);


y = data2014.crude_rate_20_64;

X =[data2014.G2014, data2014.post2014, data2014.treatPost2014];
reg2014_noW = fitlm(X,y,VarNames={'isTreated','post','treatPost','logDeathRate'})


%% Covariates and IPW
fprintf('  ------- Covariates and IPW (2014) ------\n')
ds2014 = did.Dataset.fromTable(data2014, idVar="county", timeVar="year", dVar="treatPost2014", ...
    yVar="crude_rate_20_64",describe=true);
% Desc statistics Covariates
groupsummary(data2014(data2014.year==2013,:),"expansionstatus",["mean","median"],["perc_female","perc_white","perc_hispanic","unemp_rate_pc","poverty_rate", "median_income_k"])
% Callaway/Sant'Anna OR
resCS = did.fit("cs", ds2014, ...
    Comparison="never", Delta=0, Approach="or", ...
    Covariates=["perc_female","perc_white","perc_hispanic","unemp_rate_pc","poverty_rate", "median_income_k"], ...
    Seed=42, SEMethod="clustered");
% Callaway/Sant'Anna DR
resCS = did.fit("cs", ds2014, ...
    Comparison="never", Delta=0, Approach="dr", ...
    Covariates=["perc_female","perc_white","perc_hispanic","unemp_rate_pc","poverty_rate", "median_income_k"], ...
    Seed=42, SEMethod="clustered");
%% Staggered DID
fprintf('\n  ------- Full sample Staggered DID ------\n')
% prepare data
data.isTreated =  ~isnan(data.yaca);
data.treatPost = data.year >= data.treat_year & data.treat_year~=0;

data.county = string(data.county);

ds = did.Dataset.fromTable(data,idVar="county",timeVar="year",yVar="crude_rate_20_64",dVar="treatPost");
resBacon = did.fit("Bacon",ds);
% No Covariates

% TWFE (biased!)
res = did.fit("twfe", ds,...
    vcov="clustered", clusters=["county","year"]);
% Wooldridge
resWool = did.fit("Wooldridge", ds,clusters=["county"]);
% DID_M
res = did.fit("DID_M", ds);
resCS_Uncond = did.fit("CS", ds, Approach="unconditional",...
    SEMethod="clustered", clusterVar="county", PreTrendBase="varying", Display=true);
disp("=== CS Unconditional Results ===");
disp(resCS_Uncond.summaryTable);
%% With Covariates

resWool = did.fit("Wooldridge", ds,Covariates=["perc_female","perc_white","perc_hispanic","unemp_rate_pc","poverty_rate", "median_income_k"],...
    clusters=["county","year"]);

resCS = did.fit("CS", ds, Approach="dr",Covariates=["perc_female","perc_white","perc_hispanic","unemp_rate_pc","poverty_rate", "median_income_k"],...
    SEMethod="clustered",ClusterVar="county", PreTrendBase="varying");
disp(resCS.summaryTable);
resBJS =did.fit("BJS",ds,Covariates=["perc_female","perc_white","perc_hispanic","unemp_rate_pc","poverty_rate", "median_income_k"],...
    useParallel=8, Display=true)
