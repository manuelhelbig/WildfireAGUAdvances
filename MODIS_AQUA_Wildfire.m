%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% code to analyse MODIS Aqua observations at fire and reference sites %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author: Manuel Helbig (email: helbig@gfz-potsdam.de)
% last update: August 2, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read MODIS LST data
clear all
path1 = '.../MODIS/LST_A/day/FIRE/'; % navigate to data folder
cd('...MODIS/LST_A/day/FIRE/');
filePattern = fullfile(path1, '*.csv');
finfo = dir(filePattern);

pathr = '...MODIS/LST_A/day/REF/';
filePattern = fullfile(pathr, '*.csv');
finfor = dir(filePattern);
%% extract year of fire data from filename
for k=1:length(finfo) % fire sites
    y = str2num(finfo(k).name(34:35));
    
    if y<24
        wYR(k)= str2num(strcat('20',finfo(k).name(34:35))); 
    else
        wYR(k)= str2num(strcat('19',finfo(k).name(34:35)));
    end
    clear y
end

for k=1:length(finfor) % reference sites
    y = str2num(finfor(k).name(32:33));
    if y<24 
        wYRr(k)= str2num(strcat('20',finfor(k).name(32:33)));
    else
        wYRr(k)= str2num(strcat('19',finfor(k).name(32:33)));
    end
    clear y
end

%% read land surface temperature data
for n=1:length(finfo)

    opts = delimitedTextImportOptions("NumVariables", 7);
    filename = [path1 finfo(n).name]; % fire
    % Specify range and delimiter
    opts.DataLines = [1, Inf];
    opts.Delimiter = ",";

    % Specify column names and types
    opts.VariableNames = ["Var1", "MOD", "Var3", "LatLong", "Var5", "Var6", "LST"];
    opts.SelectedVariableNames = ["MOD","LatLong","LST"];
    opts.VariableTypes = ["string", "string", "string", "string", "string", "string", "double"];

    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";

    % Specify variable properties
    opts = setvaropts(opts, ["Var1", "Var3", "Var5", "Var6"], "WhitespaceRule", "preserve");
    opts = setvaropts(opts, ["Var1", "Var3", "Var5", "Var6"], "EmptyFieldRule", "auto");

    % Import the data
    tbl = readtable(filename, opts);

    % Convert to output type
    MAT = tbl.LST;
    
    LatLong = tbl.LatLong;
    MOD{n}=tbl.MOD{1};
       
    newStr = split(LatLong{1},'Lat');
    newStr1 = split(newStr{2},'Lon');
    newStr2 = split(newStr1{2},'Samp');
    Lat(n)=str2num(newStr1{1});
    Lon(n)=str2num(newStr2{1});
    
    % Clear temporary variables
    clear opts tbl
    
    opts = delimitedTextImportOptions("NumVariables", 7);
    
    filename = [pathr finfor(n).name]; % reference
    % Specify range and delimiter
    opts.DataLines = [1, Inf];
    opts.Delimiter = ",";

    % Specify column names and types
    opts.VariableNames = ["Var1", "MOD", "Var3", "LatLong", "Var5", "Var6", "LST"];
    opts.SelectedVariableNames = ["MOD","LatLong","LST"];
    opts.VariableTypes = ["string", "string", "string", "string", "string", "string", "double"];

    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";

    % Specify variable properties
    opts = setvaropts(opts, ["Var1", "Var3", "Var5", "Var6"], "WhitespaceRule", "preserve");
    opts = setvaropts(opts, ["Var1", "Var3", "Var5", "Var6"], "EmptyFieldRule", "auto");

    % Import the data
    tbl = readtable(filename, opts);

    % Convert to output type
    MATr = tbl.LST;
    MODr{n} = tbl.MOD{1};
    LatLong = tbl.LatLong;
       
    newStr = split(LatLong{1},'Lat');
    newStr1 = split(newStr{2},'Lon');
    newStr2 = split(newStr1{2},'Samp');
    Latr(n)=str2num(newStr1{1});
    Lonr(n)=str2num(newStr2{1});
    % Clear temporary variables
    clear opts tbl
 
    % Clear temporary variables
    clear opts tbl
    
    opts = delimitedTextImportOptions("NumVariables", 7);

    % Specify range and delimiter
    opts.DataLines = [1, Inf];
    opts.Delimiter = ",";

    % Specify column names and types
    opts.VariableNames = ["Var1", "Var2", "A2002185", "Var4", "Var5", "Var6", "Var7"];
    opts.SelectedVariableNames = "A2002185";
    opts.VariableTypes = ["string", "string", "double", "string", "string", "string", "string"];

    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";

    % Specify variable properties
    opts = setvaropts(opts, ["Var1", "Var2", "Var4", "Var5", "Var6", "Var7"], "WhitespaceRule", "preserve");
    opts = setvaropts(opts, ["Var1", "Var2", "Var4", "Var5", "Var6", "Var7"], "EmptyFieldRule", "auto");
    opts = setvaropts(opts, "A2002185", "TrimNonNumeric", true);
    opts = setvaropts(opts, "A2002185", "ThousandsSeparator", ",");

    % Import the data
    filename = [path1 finfo(n).name];
    DT = readtable(filename, opts);

    % Convert to output type
    DT = table2cell(DT);
    numIdx = cellfun(@(x) ~isnan(str2double(x)), DT);
    DT(numIdx) = cellfun(@(x) {str2double(x)}, DT(numIdx));

    % Clear temporary variables
    clear opts
    % get year and DOY from datetime vector
    for k=1:length(DT)
        t=num2str(DT{k});
        YEAR(k)=str2num(t(1:4));
        DOY(k)=str2num(t(5:end));
        clearvars t
    end
    % create Matlab datetime, month, and day vector
    DT=doy2date(DOY,YEAR);
    for k=1:length(DT)
        MM(k)=str2num(datestr(DT(k),'mm'));
        DD(k)=str2num(datestr(DT(k),'dd'));
    end

    % get daily mean per day of year
    s=0;
    uYY=unique(YEAR);
    for k=1:length(uYY)
       for r=1:365
           lst=nanmean(MAT(DOY==r & YEAR==uYY(k) & YEAR<2023));
           lst_r=nanmean(MATr(DOY==r & YEAR==uYY(k) & YEAR<2023));
           dlst=lst-lst_r;
           
           if ~isempty(dlst)
              LST{n}(r,k)= dlst;
              LSTf{n}(r,k)= lst;
              LSTr{n}(r,k)= lst_r;
             
           else
              LST{n}(r,k)= NaN; 
              LSTf{n}(r,k)= NaN; 
              LSTr{n}(r,k)= NaN; 
              
           end
           
        clearvars dlst lst lst_r
       end
       
    end
    clearvars -except MOD MODr MODn MODnr LST Lat Latr Lonr Lon LSTr LSTf LSTnf sYY nYY pathr finfor path1 finfonr finfo wYR n path2 finfo2
end
%% read temperature and time since fire
p=0;
uYY=2002:2023;

sYY=[];
lst=[];
lstf=[];
lstr=[];
% loop through sites
for k=1:142
    lst=horzcat(lst,LST{k}); % temperature difference
    lstf=horzcat(lstf,LSTf{k}-273.15); % fire site temperature
    lstr=horzcat(lstr,LSTr{k}-273.15); % reference site temperature
    sz=size(LST{k});
    for r = 1:sz(2)
        sYY=vertcat(sYY,uYY(r)-wYR(k)); % years since fire
    end
    sit(k*22-21:k*22)=ones(1,22).*k;
end
%% only choose days with data
ind = find(sum(isnan(lstf'))<length(sYY));
DOY=ind;
lstf=lstf(ind,:);
lstr=lstr(ind,:);
%% filter for outliers (> 3 x standard deviation)
for k=1:46
    ind=find(abs(lstf(k,:)-nanmedian(lstf(k,:)))>=3.*nanstd(lstf(k,:)));
    if ~isempty(ind)
        lstf(k,ind)=NaN;
    end

    ind=find(abs(lstr(k,:)-nanmedian(lstr(k,:)))>=3.*nanstd(lstr(k,:)));
    if ~isempty(ind)
        lstr(k,ind)=NaN;
    end
end
lst=lstf-lstr;
%% calculate median temperature differences for post-fire decades
ty=10; % for how many years should averages be calculated
LSTmat=NaN(46,90./ty);
LSTe=NaN(46,90./ty);

for k=1:90./ty
    if k==1 % differences pre-fire
        LSTmat(:,k)=nanmedian(lst(:,sYY<0)');
        LSTe(:,k)=nanstd(lst(:,sYY<0)')./sqrt(sum(sYY<0));
    elseif k==90./ty % differences for 
        LSTmat(:,k)=nanmedian(lst(:,sYY>(k-1)*ty-ty)');
        LSTe(:,k)=nanstd(lst(:,sYY>(k-1)*ty-ty)')./sqrt(sum(sYY>(k-1)*ty-ty));
    else
        LSTmat(:,k)=nanmedian(lst(:,sYY>(k-1)*ty-ty & sYY<=(k-1)*ty)');
        LSTe(:,k)=nanstd(lst(:,sYY>(k-1)*ty-ty & sYY<=(k-1)*ty)')./sqrt(sum(sYY>(k-1)*ty-ty & sYY<=(k-1)*ty));
    end
end

%% plot figure 2 (Aqua)
figure,
cmap=colormap(brewermap(90./ty-1,'BrBG'));
clear albx
for k=1:90./ty-1
    n(k)=plotshaded(DOY,([LSTmat(:,k+1)-LSTe(:,k+1) LSTmat(:,k+1)+LSTe(:,k+1)])',[0.6 0.6 0.6]);
    hold on
end
for k=1:90./ty-1
    lstx(k)=plot(DOY,LSTmat(:,k+1),'-','LineWidth',2.5,'Color',cmap(k,:));
    hold on
end
lstx(k+1)=plot(DOY,LSTmat(:,1),':','Color',[0.2 0.2 0.2],'LineWidth',2.5);
set(gca,'FontSize',14)
refline(0,0)
[hleg,att] = legend(lstx,{'1 - 10 yr','11 - 20 yr','21 - 30 yr','31 - 40 yr','41 - 50 yr','51 - 60 yr','61 - 70 yr','> 70 yr','< 0 yr'});
title(hleg,'Years since fire')
xlim([0.5 365.5])
ylabel('\Delta LST (fire - control) [\circC]');
xlabel('Day of year');
text(10,4.7,'Aqua','FontSize',13)
set(gca,'FontSize',14)
text(280,4.7,'Years after fire');
%exportgraphics(gcf,'FigS2.pdf','ContentType','vector')

%% estimate impact of wildfires on surface temperatures %%%%%%%%%%%%%%%%%%%
% read burned area for Canada
data=importdata('...WildfireAreaCan.csv'); % read burned area data
Wf.Year=data.data(:,1);
Wf.FireHist=data.data(:,2); % historical burned area
mFIRE=mean(Wf.FireHist(Wf.Year>2010 & Wf.Year<=2020)); % calculate average burned area between 2011 and 2020
fFIREl=mFIRE+0.36*mFIRE; % low burn scenario (in 2050)
fFIREh=mFIRE+1.5*mFIRE; % high burn scenario (in 2050)
Wf.FutL=Wf.FireHist;
Wf.FutH=Wf.FireHist;

%% estimate future burned area for years 2024 - 2050
for k=1:27
    Wf.FutH(38+k)=mFIRE+((1.5*mFIRE)./30).*(k+3);
end

for k=1:27
    Wf.FutL(38+k)=mFIRE+((0.36*mFIRE)./30).*(k+3);
end

Wf.FirePercL=(Wf.FutL./544163659.1).*100; % calculate burned area in %
Wf.FirePercH=(Wf.FutH./544163659.1).*100;

%% randomly select one year per site to estimate uncertainties
for k=1:1000
    test=37;
    while test<38 % only use time series if all 38 years are included
        ind=randsample(23,142,'true');
        for n=1:142
            % select only one year per site
            sel(n)=ind(n)+n*24-24; 
        end
        LSTsel=lst(:,sel); % get land surface temps per site
        yr_sel=sYY(sel); % get years
        uyr=unique(yr_sel); 
        uyr=uyr(uyr>0); % only choose years after fire
        uyr(uyr>38)=38; % set years past 38 years to 38 years
        uyr=unique(uyr);
        test=length(uyr); % number of years
        LSTselect{k}=LSTsel; % save LST and years
        yr_select{k}=yr_sel;
    end
end
%% calculate median and standard deviation of estimates
for s=1:1000
    yr_sel=yr_select{s};
    LST_sel=LSTselect{s};
    for n=1:46
        for k=1:40
            if k==1 % pre-fire
                i=find(yr_sel<0);
                dT(n,k,s)=nanmedian(LST_sel(n,i));
                sdT(n,k,s)=nanstd(LST_sel(n,i));
            elseif k==40 % > 38 years after fire
                i=find(yr_sel>=k-2);
                dT(n,k,s)=nanmedian(LST_sel(n,i));
                sdT(n,k,s)=nanstd(LST_sel(n,i));
            else
                i=find(yr_sel==k-2);
                dT(n,k,s)=nanmedian(LST_sel(n,i));
                sdT(n,k,s)=nanstd(LST_sel(n,i));
            end
        end
    end
end
%% estimate surface temperature impact per year (high burn scenario annual impact)
for r=1:28 % for year r = 1 (2024) and for years r = 2:29 (2025:2050)
    if r==1
        for s=1:1000
            for k=1:38 % year after fire
                for n=1:46 % day of year
                    if sum(isnan(dT(n,:,s)))==40
                        dT_Fire(k,n,s)=NaN;
                    else
                        % calculate contributions for individual years
                        dT_Fire(k,n,s)=Wf.FirePercH(38-k+1)./100.*dT(n,k+2,s); 
                    end
                end
            end
        end
        % sum contributions for all years to estimate temp impact for 2024
        dT_tot(r,:)=nanmean(squeeze(nansum(dT_Fire,1)));
    else
        % calculate temp contributions for years 2025 - 2050
        for s=1:1000
            for k=1:38 % year after fire
                for n=1:46 % day of year
                    if sum(isnan(dT(n,:,s)))==40
                        dT_Fire(k,n,s)=NaN;
                    else
                        dT_Fire(k,n,s)=Wf.FirePercH((38+r-1)-k+1)./100.*dT(n,k+2,s);
                    end
                end
            end
        end
        % sum contributions for all years to estimate temp impact for 2025 - 2050
        dT_tot(r,:)=nanmean(squeeze(nansum(dT_Fire,1)));
    end
end

%% estimate surface temperature impact per year (high burn scenario summar and winter impacts)
for r=1:28
    if r==1
        for s=1:1000
            for k=1:38 % year after fire
                for n=1:46 % day of year
                    if sum(isnan(dT(n,:,s)))==40
                        dT_Fire(k,n,s)=NaN;
                    else
                        dT_Fire(k,n,s)=Wf.FirePercH(38-k+1)./100.*dT(n,k+2,s);
                    end
                end
            end
        end
        dT_tot_S(r,:)=nanmean(squeeze(nansum(dT_Fire(:,DOY>=182 & DOY<=273,:),1)));
        dT_tot_W(r,:)=nanmean(squeeze(nansum(dT_Fire(:,DOY>=32 & DOY<=120,:),1)));
    else
        for s=1:1000
            for k=1:38 % year after fire
                for n=1:46 % day of year
                    if sum(isnan(dT(n,:,s)))==40
                        dT_Fire(k,n,s)=NaN;
                    else
                        dT_Fire(k,n,s)=Wf.FirePercH((38+r-1)-k+1)./100.*dT(n,k+2,s);
                    end
                end
            end
        end
        dT_tot_S(r,:)=nanmean(squeeze(nansum(dT_Fire(:,DOY>=182 & DOY<=273,:),1)));
        dT_tot_W(r,:)=nanmean(squeeze(nansum(dT_Fire(:,DOY>=32 & DOY<=120,:),1)));
    end
end

%% estimate surface temperature impact per year (low burn scenario annual impact)
for r=1:28 % for year r = 1 (2024) and for years r = 2:29 (2025:2050)
    if r==1
        for s=1:1000
            for k=1:38 % year after fire
                for n=1:46 % day of year
                    if sum(isnan(dT(n,:,s)))==40
                        dT_Fire(k,n,s)=NaN;
                    else
                        % calculate contributions for individual years
                        dT_Fire(k,n,s)=Wf.FirePercL(38-k+1)./100.*dT(n,k+2,s); 
                    end
                end
            end
        end
        % sum contributions for all years to estimate temp impact for 2024
        dT_tot_L(r,:)=nanmean(squeeze(nansum(dT_Fire,1)));
    else
        % calculate temp contributions for years 2025 - 2050
        for s=1:1000
            for k=1:38 % year after fire
                for n=1:46 % day of year
                    if sum(isnan(dT(n,:,s)))==40
                        dT_Fire(k,n,s)=NaN;
                    else
                        dT_Fire(k,n,s)=Wf.FirePercL((38+r-1)-k+1)./100.*dT(n,k+2,s);
                    end
                end
            end
        end
        % sum contributions for all years to estimate temp impact for 2025 - 2050
        dT_tot_L(r,:)=nanmean(squeeze(nansum(dT_Fire,1)));
    end
end

%% estimate surface temperature impact per year (high burn scenario summar and winter impacts)
for r=1:28
    if r==1
        for s=1:1000
            for k=1:38 % year after fire
                for n=1:46 % day of year
                    if sum(isnan(dT(n,:,s)))==40
                        dT_Fire(k,n,s)=NaN;
                    else
                        dT_Fire(k,n,s)=Wf.FirePercL(38-k+1)./100.*dT(n,k+2,s);
                    end
                end
            end
        end
        dT_tot_S_L(r,:)=nanmean(squeeze(nansum(dT_Fire(:,DOY>=182 & DOY<=273,:),1)));
        dT_tot_W_L(r,:)=nanmean(squeeze(nansum(dT_Fire(:,DOY>=32 & DOY<=120,:),1)));
    else
        for s=1:1000
            for k=1:38 % year after fire
                for n=1:46 % day of year
                    if sum(isnan(dT(n,:,s)))==40
                        dT_Fire(k,n,s)=NaN;
                    else
                        dT_Fire(k,n,s)=Wf.FirePercL((38+r-1)-k+1)./100.*dT(n,k+2,s);
                    end
                end
            end
        end
        dT_tot_S_L(r,:)=nanmean(squeeze(nansum(dT_Fire(:,DOY>=182 & DOY<=273,:),1)));
        dT_tot_W_L(r,:)=nanmean(squeeze(nansum(dT_Fire(:,DOY>=32 & DOY<=120,:),1)));
    end
end

%% plot temperature impacts for low and high burn scenarios (Fig. 10)
figure,
subplot(1,2,1)
plotshaded(2024:2050,prctile(dT_tot_S(1:27,:)',[2.5 97.5]),'r')
hold on
l1=plot(2024:2050,nanmedian(dT_tot_S(1:27,:)'));
plotshaded(2024:2050,prctile(dT_tot_W(1:27,:)',[2.5 97.5]),'r')
l2=plot(2024:2050,nanmedian(dT_tot_W(1:27,:)'));

plotshaded(2024:2050,prctile(dT_tot_S_L(1:27,:)',[2.5 97.5]),'r')
hold on
l3=plot(2024:2050,nanmedian(dT_tot_S_L(1:27,:)'));
plotshaded(2024:2050,prctile(dT_tot_W_L(1:27,:)',[2.5 97.5]),'r')
l4=plot(2024:2050,nanmedian(dT_tot_W_L(1:27,:)'));
refline(0,0)
xlim([2022.5 2050.5])
ylim([-0.15 0.5])
ylabel('\Delta LST (fire - no fire) [\circC]')
legend([l1 l2 l3 l4],"Summer [high burn area]","Winter [high burn area]","Summer [low burn area]","Winter [low burn area]")
set(gca,'FontSize',14)
subplot(1,2,2)
plotshaded(2024:2050,prctile(dT_tot(1:27,:)',[2.5 97.5]),'r')
hold on
l1=plot(2024:2050,nanmedian(dT_tot(1:27,:)'));
plotshaded(2024:2050,prctile(dT_tot_L(1:27,:)',[2.5 97.5]),'r')
l2=plot(2024:2050,nanmedian(dT_tot_L(1:27,:)'));
refline(0,0)
xlim([2022.5 2050.5])
legend([l1 l2],"Annual [high burn area]","Annual [low burn area]")
set(gca,'FontSize',14)
ylim([-0.15 0.5])

%% read MODIS leaf area index data
clear all
path1 = '...MODIS/LAI/FIRE/';
filePattern = fullfile(path1, '*.csv');
finfo = dir(filePattern);

pathr = '...MODIS/LAI/REF/';
filePattern = fullfile(pathr, '*.csv');
finfor = dir(filePattern);
%% get years of fire disturbance
for k=1:length(finfo)
    y = str2num(finfo(k).name(29:30));
    if y<24
        wYR(k)= str2num(strcat('20',finfo(k).name(29:30)));
    else
        wYR(k)= str2num(strcat('19',finfo(k).name(29:30)));
    end
    clear y
end

for k=1:length(finfor)
    y = str2num(finfor(k).name(29:30));
    if y<24
        wYRr(k)= str2num(strcat('20',finfor(k).name(29:30)));
    else
        wYRr(k)= str2num(strcat('19',finfor(k).name(29:30)));
    end
    clear y
end

%% read MODIS LAI data
for n=1:length(finfo)
    % reference
    filename = [path1 finfo(n).name];
    opts = delimitedTextImportOptions("NumVariables", 7);
    
    opts.DataLines = [1, Inf];
    opts.Delimiter = ",";

    % Specify column names and types
    opts.VariableNames = ["MYD15A2HA2002185h12v030612020072165944Lai_500m", "MYD15A2H", "A2002185", "Lat5943719508Lon1112102748Samp1Line1", "VarName5", "Lai_500m", "VarName7"];
    opts.VariableTypes = ["string", "categorical", "double", "string", "double", "double", "double"];

    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";

    % Specify variable properties
    opts = setvaropts(opts, "MYD15A2HA2002185h12v030612020072165944Lai_500m", "WhitespaceRule", "preserve");
    opts = setvaropts(opts, ["MYD15A2HA2002185h12v030612020072165944Lai_500m", "MYD15A2H", "Lat5943719508Lon1112102748Samp1Line1"], "EmptyFieldRule", "auto");
    opts = setvaropts(opts, ["A2002185", "Lai_500m"], "TrimNonNumeric", true);
    opts = setvaropts(opts, ["A2002185", "Lai_500m"], "ThousandsSeparator", ",");

    % Import the data
    tbl = readtable(filename, opts);

    DT = tbl.A2002185;
    LatLong = tbl.Lat5943719508Lon1112102748Samp1Line1;
    MAT = tbl.VarName7;
       
    newStr = split(LatLong{1},'Lat');
    newStr1 = split(newStr{2},'Lon');
    newStr2 = split(newStr1{2},'Samp');
    Lat(n)=str2num(newStr1{1});
    Lon(n)=str2num(newStr2{1});
    
    clear opts tbl newStr newStr1 newStr2
        
    % fire site
    filename = [pathr finfor(n).name];
    opts = delimitedTextImportOptions("NumVariables", 7);
    
    opts.DataLines = [1, Inf];
    opts.Delimiter = ",";

    % Specify column names and types
    opts.VariableNames = ["MYD15A2HA2002185h12v030612020072165944Lai_500m", "MYD15A2H", "A2002185", "Lat5943719508Lon1112102748Samp1Line1", "VarName5", "Lai_500m", "VarName7"];
    opts.VariableTypes = ["string", "categorical", "double", "string", "double", "double", "double"];

    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";

    % Specify variable properties
    opts = setvaropts(opts, "MYD15A2HA2002185h12v030612020072165944Lai_500m", "WhitespaceRule", "preserve");
    opts = setvaropts(opts, ["MYD15A2HA2002185h12v030612020072165944Lai_500m", "MYD15A2H", "Lat5943719508Lon1112102748Samp1Line1"], "EmptyFieldRule", "auto");
    opts = setvaropts(opts, ["A2002185", "Lai_500m"], "TrimNonNumeric", true);
    opts = setvaropts(opts, ["A2002185", "Lai_500m"], "ThousandsSeparator", ",");

    % Import the data
    tbl = readtable(filename, opts);

    DTr = tbl.A2002185;
    LatLongr = tbl.Lat5943719508Lon1112102748Samp1Line1;
    MATr = tbl.VarName7;
       
    newStr = split(LatLongr{1},'Lat');
    newStr1 = split(newStr{2},'Lon');
    newStr2 = split(newStr1{2},'Samp');
    Latr(n)=str2num(newStr1{1});
    Lonr(n)=str2num(newStr2{1});
    
    clear opts tbl newStr newStr1 newStr2
    
        for k=1:length(DTr)
            t=num2str(DTr(k));
            YEARr(k)=str2num(t(1:4));
            DOYr(k)=str2num(t(5:end));
            clearvars t
        end
        DTr=doy2date(DOYr,YEARr);
        for k=1:length(DTr)
            MMr(k)=str2num(datestr(DTr(k),'mm'));
            DDr(k)=str2num(datestr(DTr(k),'dd'));
        end
        
        for k=1:length(DT)
            t=num2str(DT(k));
            YEAR(k)=str2num(t(1:4));
            DOY(k)=str2num(t(5:end));
            clearvars t
        end
        DT=doy2date(DOY,YEAR);
        for k=1:length(DT)
            MM(k)=str2num(datestr(DT(k),'mm'));
            DD(k)=str2num(datestr(DT(k),'dd'));
        end

    s=0;
    uYY=unique(YEAR);
    for k=1:length(uYY)
       for r=1:365
           lai=nanmean(MAT(DOY==r & YEAR==uYY(k) & YEAR<2023));
           lai_r=nanmean(MATr(DOYr==r & YEARr==uYY(k) & YEARr<2023));
           dlai=lai-lai_r;
           
           if ~isempty(dlai)
              LAI{n}(r,k)= dlai;
              LAIf{n}(r,k)= lai;
              LAIr{n}(r,k)= lai_r;
           else
              LAI{n}(r,k)= NaN; 
              LAIf{n}(r,k)= NaN; 
              LAIr{n}(r,k)= NaN; 
           end
           
        clearvars dlai lai lai_r
       end
        
    end
    clearvars -except LAI LAIf LAIr sYY nYY path1 finfo wYR n Lat Lon pathr finfor Latr Lonr
end

%% read data per site and remove 0 values
for k=1:142

    lai=LAIf{n};
    ind=find(lai==0);
    lai(ind)=NaN;
    LAIf{n}(ind)=NaN;

    lai=LAIr{n};
    ind=find(lai==0);
    lai(ind)=NaN;
    LAIr{n}(ind)=NaN;
end
%% calculate differences between fire and reference sites
p=0;
uYY=2002:2023;
sYY=[];
lai_1D=[];
for k=1:142

    lai_1D=horzcat(lai_1D,(LAIf{k}-LAIr{k}));
    sz=size(LAIf{k});
    for r = 1:sz(2)
        sYY=vertcat(sYY,uYY(r)-wYR(k));
        syy(r,k)=uYY(r)-wYR(k);
    end
    sit(k*22-21:k*22)=ones(1,22).*k;
end
%% get leaf aread index difference for July to September
clear sYY
uYY=2002:2023;
mLAI=NaN(22,142);
timestamp=doy2date(1:365,ones(1,365).*2017);
MM=str2num(datestr(timestamp,'mm'));
for k=1:142
    lai=LAI{k};
    laif=LAIf{k};
    lair=LAIr{k};
    sz=size(lai);
    for n=1:sz(2)
        mLAI(n,k)=(nanmedian(laif(MM>=7 & MM<=9,n))-nanmedian(lair(MM>=7 & MM<=9,n))); 
        mLAI2(n,k)=(nanmedian(laif(MM>=7 & MM<=9,n))-nanmedian(lair(MM>=7 & MM<=9,n)))./nanmedian(lair(MM>=7 & MM<=9,n));
        sYY(n,k)=uYY(n)-wYR(k);
    end
    clear lai sz laif lair
end

%% calculate median and standard error per year after fire
lng=142;
years = -1:80;
lai = reshape(mLAI2.*100,22*lng,1);%
years_1D = reshape(sYY,22*lng,1);

for k=1:length(years)
    if k==1
        ind = find(years_1D<=years(k));
    elseif k==80
        ind = find(years_1D>=80);
    else
        ind = find(years_1D==years(k));
    end
    if length(ind)>=10
        
        lai_t(k)=nanmedian(lai(ind));
        lai_t_sd(k)=nanstd(lai(ind));
        lai_t_se(k)=nanstd(lai(ind))./sqrt(length(ind));
        sample_n(k)=length(ind);
    else
        lai_t(k)=NaN;
        lai_t_sd(k)=NaN;
        lai_t_se(k)=NaN;
        sample_n(k)=NaN;

    end
    
    clear ind
end
%% remove NaN values
lai_t_sd=lai_t_sd(~isnan(lai_t));
lai_t_se=lai_t_se(~isnan(lai_t));
years=years(~isnan(lai_t));
lai_t=lai_t(~isnan(lai_t));
lai_t(2)=NaN;
%% plot change in LAI since year of fire
figure,
hold on
plotshaded(years(3:end),[(lai_t(3:end)-lai_t_se(3:end))' (lai_t(3:end)+lai_t_se(3:end))'],'r');
plot(years(1)-1,lai_t(1),'ok');
hold on
plot(years(3:end),lai_t(3:end),'k');
refline(0,0)
plot([0 0],[-60 40],'--');

xlim([-4 80]);
xlabel('Year since fire');
ylabel('\Deltaleaf area index (fire - control) [%]');
h=text(-1.5,0.3,'pre-fire');
set(h,'Rotation',90);

%% read MODIS albedo data
clear all
path = '...MODIS/ALB/REF/';
filePattern = fullfile(path, '*.csv');
finfo = dir(filePattern);

path1 = '...MODIS/ALB/FIRE/';
filePattern = fullfile(path1, '*.csv');
finfo1 = dir(filePattern);

%% get years since fire information
for k=1:length(finfo)
    if length(finfo(k).name)==42
        SITE{k}=finfo(k).name(33:37);
    else
        SITE{k}=finfo(k).name(33:39);
    end
    y = str2num(finfo(k).name(36:37));
    if y<24
        wYR(k)= str2num(strcat('20',finfo(k).name(36:37)));
    else
        wYR(k)= str2num(strcat('19',finfo(k).name(36:37)));
    end
    clear y
end

for k=1:length(finfo1)
    if length(finfo1(k).name)==42
        SITE_FIRE{k}=finfo1(k).name(33:37);
    else
        SITE_FIRE{k}=finfo1(k).name(33:39);
    end
    y = str2num(finfo1(k).name(36:37));
    if y<24
        wYR_FIRE(k)= str2num(strcat('20',finfo1(k).name(36:37)));
    else
        wYR_FIRE(k)= str2num(strcat('19',finfo1(k).name(36:37)));
    end
    clear y
end
%% get data from MODIS files

for n=1:length(finfo)
    % reference
    filename = [path finfo(n).name];
    opts = delimitedTextImportOptions("NumVariables", 7);

    opts.DataLines = [1, Inf];
    opts.Delimiter = ",";

    opts.VariableNames = ["MCD43AA2000055h11v030612020038132853shortwave_white", "MCD43A", "A2000055", "Lat59073300000000Lon120666800000000Samp1Line1", "VarName5", "shortwave_white", "F"];
    opts.VariableTypes = ["string", "double", "double", "string", "double", "categorical", "double"];

    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";

    opts = setvaropts(opts, "MCD43AA2000055h11v030612020038132853shortwave_white", "WhitespaceRule", "preserve");
    opts = setvaropts(opts, ["MCD43AA2000055h11v030612020038132853shortwave_white", "Lat59073300000000Lon120666800000000Samp1Line1", "shortwave_white"], "EmptyFieldRule", "auto");
    opts = setvaropts(opts, ["MCD43A", "A2000055"], "TrimNonNumeric", true);
    opts = setvaropts(opts, ["MCD43A", "A2000055"], "ThousandsSeparator", ",");

    tbl = readtable(filename, opts);

    LatLong = tbl.Lat59073300000000Lon120666800000000Samp1Line1;
    DTr = tbl.A2000055;
    MATr = tbl.F;
    
    newStr = split(LatLong{1},'Lat');
    newStr1 = split(newStr{2},'Lon');
    newStr2 = split(newStr1{2},'Samp');
    Latr(n)=str2num(newStr1{1});
    Lonr(n)=str2num(newStr2{1});
    
    clear opts tbl newStr newStr1 newStr2
        
    % fire site
    filename = [path1 finfo1(n).name];
    opts = delimitedTextImportOptions("NumVariables", 7);

    opts.DataLines = [1, Inf];
    opts.Delimiter = ",";

    opts.VariableNames = ["MCD43AA2000055h11v030612020038132853shortwave_white", "MCD43A", "A2000055", "Lat59073300000000Lon120666800000000Samp1Line1", "VarName5", "shortwave_white", "F"];
    opts.VariableTypes = ["string", "double", "double", "string", "double", "categorical", "double"];

    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";

    opts = setvaropts(opts, "MCD43AA2000055h11v030612020038132853shortwave_white", "WhitespaceRule", "preserve");
    opts = setvaropts(opts, ["MCD43AA2000055h11v030612020038132853shortwave_white", "Lat59073300000000Lon120666800000000Samp1Line1", "shortwave_white"], "EmptyFieldRule", "auto");
    opts = setvaropts(opts, ["MCD43A", "A2000055"], "TrimNonNumeric", true);
    opts = setvaropts(opts, ["MCD43A", "A2000055"], "ThousandsSeparator", ",");

    tbl = readtable(filename, opts);

    LatLong = tbl.Lat59073300000000Lon120666800000000Samp1Line1;
    DT = tbl.A2000055;
    MAT = tbl.F;
    
    newStr = split(LatLong{1},'Lat');
    newStr1 = split(newStr{2},'Lon');
    newStr2 = split(newStr1{2},'Samp');
    Lat(n)=str2num(newStr1{1});
    Lon(n)=str2num(newStr2{1});
    
    clear opts tbl newStr newStr1 newStr2
    
        for k=1:length(DTr)
            t=num2str(DTr(k));
            YEARr(k)=str2num(t(1:4));
            DOYr(k)=str2num(t(5:end));
            clearvars t
        end
        DTr=doy2date(DOYr,YEARr);
        for k=1:length(DTr)
            MMr(k)=str2num(datestr(DTr(k),'mm'));
            DDr(k)=str2num(datestr(DTr(k),'dd'));
        end
        
        for k=1:length(DT)
            t=num2str(DT(k));
            YEAR(k)=str2num(t(1:4));
            DOY(k)=str2num(t(5:end));
            clearvars t
        end
        DT=doy2date(DOY,YEAR);
        for k=1:length(DT)
            MM(k)=str2num(datestr(DT(k),'mm'));
            DD(k)=str2num(datestr(DT(k),'dd'));
        end

    s=0;
    uYY=unique(YEAR);
    for k=1:length(uYY)
       for r=1:365
           albedo=nanmedian(MAT(DOY==r & YEAR==uYY(k) & YEAR<2023));
           albedo_r=nanmedian(MATr(DOYr==r & YEARr==uYY(k) & YEARr<2023));
           dalb=albedo-albedo_r;
           
           clearvars albedo albedo_r
       
     
           if ~isempty(dalb)
                ALB{n}(r,k)= dalb;
           else
              ALB{n}(r,k)= NaN; 
           end
           
        clearvars dalb
       end
       
    end
    clearvars -except ALB sYY nYY path1 finfo wYR n Lat Lon path finfo1 Latr Lonr
end

ALB{25}=horzcat(ALB{25},NaN(365,1)); % add missing year to site
%% get year and albedo
p=0;
uYY=2000:2023;
sYY=[];
alb=[];

% loop through sites
for k=1:142
    alb=horzcat(alb,ALB{k});
    sz=size(ALB{k});
    for r = 1:sz(2)
        sYY=vertcat(sYY,uYY(r)-wYR(k));
    end
    sit(k*24-23:k*24)=ones(1,24).*k;
end

%% only use days with data
ind = find(sum(isnan(alb'))<length(sYY));
DOY=ind;
alb=alb(ind,:);
%% get summer and winter data
s=0;
usites=unique(sit);
for k=1:length(usites)
    ind=find(DOY>=182 & DOY<=273);
    ind2=find(DOY>=32 & DOY<=120);
    ind1=find(sit==usites(k));
    for n=1:length(ind1)
        s=s+1;
        alb_sum(s)=nanmedian(alb(ind,ind1(n)));
        alb_win(s)=nanmedian(alb(ind2,ind1(n)));
        sites(s)=usites(k);
        yyear(s)=sYY(ind1(n));
    end
end
%%
ty=10;

ind=find(sum(isnan(ALB{1})')<length(ALB{1}(2,:)));
ALBmat=NaN(length(ind),90./ty);
ALBe=NaN(length(ind),90./ty);
for k=1:90./ty
    if k==1
        sample(k)=sum(sYY<0);
        ALBmat(:,k)=nanmedian(alb(ind,sYY<0)');
        ALBe(:,k)=nanstd(alb(ind,sYY<0)')./sqrt(sum(sYY<0));
    elseif k==90./ty
        sample(k)=sum(sYY>(k-1)*ty-ty);
        ALBmat(:,k)=nanmedian(alb(ind,sYY>(k-1)*ty-ty)');
        ALBe(:,k)=nanstd(alb(ind,sYY>(k-1)*ty-ty)')./sqrt(sum(sYY>(k-1)*ty-ty));
    else
        sample(k)=sum(sYY>(k-1)*ty-ty & sYY<=(k-1)*ty);
        ALBmat(:,k)=nanmedian(alb(ind,sYY>(k-1)*ty-ty & sYY<=(k-1)*ty)');
        ALBe(:,k)=nanstd(alb(ind,sYY>(k-1)*ty-ty & sYY<=(k-1)*ty)')./sqrt(sum(sYY>(k-1)*ty-ty & sYY<=(k-1)*ty));
    end
end
%%
clear YRmat DOYmat
for k=1:90./ty
    DOYmat(1:length(ind),k)=ind;
    if k==1
        YRmat(1:length(ind),k)=-1;
        
    else
        YRmat(1:length(ind),k)=(k-1)*ty-ty;
    end
end

%%
figure,
cmap=colormap(brewermap(90./ty-1,'BrBG'));
clear albx
for k=1:90./ty-1
    n(k)=plotshaded(DOYmat(1:length(ind),1)',([ALBmat(:,k+1)-ALBe(:,k+1) ALBmat(:,k+1)+ALBe(:,k+1)])',[0.6 0.6 0.6]);
    hold on
end
for k=1:90./ty-1
    albx(k)=plot(DOYmat(1:length(ind),1),ALBmat(:,k+1),'-','LineWidth',2.5,'Color',cmap(k,:));
    hold on
end
albx(k+1)=plot(DOYmat(1:length(ind)),ALBmat(1:length(ind)),':','Color',[0.2 0.2 0.2],'LineWidth',2.5);
xlim([0.5 365.5])
refline(0,0)
set(gca,'FontSize',14)

[hleg,att] = legend(albx,{'1 - 10 yr','11 - 20 yr','21 - 30 yr','31 - 40 yr','41 - 50 yr','51 - 60 yr','61 - 70 yr','> 70 yr','< 0 yr'});
title(hleg,'Years since fire')
xlim([0.5 365.5])
ylabel('\Delta albedo (fire - control)');
xlabel('Day of year');
set(gca,'FontSize',14)