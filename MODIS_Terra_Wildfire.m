%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% code to analyse MODIS Terra observations at fire and reference sites %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author: Manuel Helbig (email: helbig@gfz-potsdam.de)
% last update: August 2, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read MODIS LST data
clear all
path1 = '.../MODIS/LST_T/day/FIRE/';
cd('.../MODIS/LST_T/day/FIRE/');
filePattern = fullfile(path1, '*.csv');
finfo = dir(filePattern);

pathr = '.../MODIS/LST_T/day/REF/';
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
uYY=2000:2023;

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
    sit(k*24-23:k*24)=ones(1,24).*k;
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

%% plot figure S2 (Terra)
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
text(10,4.7,'Terra','FontSize',13)
set(gca,'FontSize',14)
text(280,4.7,'Years after fire');
%exportgraphics(gcf,'FigS2_Terra.pdf','ContentType','vector')

%% estimate impact of wildfires on surface temperatures %%%%%%%%%%%%%%%%%%%
% read burned area for Canada
data=importdata('.../WildfireAreaCan.csv'); % read burned area data
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