%% This function calculates the peak heights of 3 baseline scans and 4 standard additions scans for samples run 
%% using NOVA 1.8 software. The peak heights are calculated and then output into an excel file, sorted in the order the
%% samples were run in NOVA
%% to run the code you will need the Matlab curve fitting toolbox
%% You will also need to download export_fig (http://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig)
%% and tight_subplot (http://www.mathworks.com/matlabcentral/fileexchange/27991-tight-subplot-nh--nw--gap--marg-h--marg-w-)
%% and add them to the Matlb search path
%% hint: for export to vector formats (pdf, eps), export_fig relies on other programs (e.g. ghostscript).
%% the easiest way to install these on mac is using macports, or similar, on Mac.
%% updated 4-6-16 RMB, MRC

function COint_final
%function integrates cobalt peak
%last modified by MRC 04/11/2016
tic; %starts counting how long the processing takes

home='/Users/mcape/Documents/WHOI/perso/RB/CO_peakint/matlab/RB_edit'; % set your home directory here

cd /Users/mcape/Dropbox/GT_Arctic/011216; % change to sample path 
% information here to go to the correct folder where you have your cobalt
% data in text files as output from NOVA

fsize=10;

%This is setting up the file format for the sample or "baseline" scans
%get sample file names, file number (sequential), created on date, and save
%to samlist
cosam=dir('*sample*txt');
ns=length(cosam);
log.sam=zeros(ns,1);

for i=1:ns
    temp=sscanf(cosam(i).name,'Co_sample_(%d).txt');
    if ~isempty(temp)
        log.sam(i)=temp;
    end
end

samlist=cellstr(num2str(log.sam));
temp=struct2cell(cosam);
samlist=sortrows(horzcat(samlist,temp(1,:)',temp(2,:)'),1);

%This is setting up the file format for the standard addition scans
%same procedure as above for addition files
coadd=dir('*addition*txt');
nd=length(coadd);
log.add=zeros(nd,1);

for i=1:nd
    temp=sscanf(coadd(i).name,'Co_addition_(%d).txt');
    if ~isempty(temp)
        log.add(i)=temp;
    end
end

addlist=cellstr(num2str(log.add));
temp=struct2cell(coadd);
addlist=sortrows(horzcat(addlist,temp(1,:)',temp(2,:)'),1);

% Calculates the number of sample sets
Nsets=nd/4;

% Setup the output csv file (You can change the file name as necessary)
fid=fopen('GEOTRACES_Arctic.csv','wt');
fprintf(fid,'Date,Filename,Peak height\n');

addst=1;
samst=1;

xh=[0;0;0;25;50;75;100]; % x values are for the peak heights, from the std additions
h=zeros(3,1); % these are the calculated peak heights

% Sets up the figure
figure('renderer','painters')

set(gcf,'color','w','inverthardcopy','off','units','centimeters',...
    'position',[1,1,22,7])

for j=1:Nsets
    ha=tight_subplot(1,3,0.08,[.18,.07],[.1,.05]);
    set(ha,'color',[0.6,0.6,0.6])

    axes(ha(1));hold on;
    for i=1:3
        [h(i),pln]=peakcalc(samlist{samst,2},samlist{samst,3},fid,i);
        samst=samst+1;
    end
    
    hold off
    c=get(gca,'children');
    l=legend(c(2:3:8),samlist(samst-3:samst-1,1),'textcolor','w','location','southeast');
    set(l,'box','off','color','none','fontsize',8)
    
    axes(ha(2));hold on;
    for k=1:4
        [h(i+k),pln]=peakcalc(addlist{addst,2},addlist{addst,3},fid,i+k);
        addst=addst+1;
    end
    
    hold off
    c=get(gca,'children');
    l=legend(c(2:3:11),addlist(addst-4:addst-1,1),'textcolor','w','location','southeast');
    set(l,'box','off','color','none','fontsize',8)
        
    %fit line
    f=fit(xh,h,'poly1','robust','bisquare');
    
    axes(ha(3));
    plot(xh,h,'bo','markerfacecolor','b','markersize',4);
    hold on
    plot(xh,f(xh),'r-','linewidth',1)
    xlabel('Cobalt added [pM]','fontsize',fsize)
    ylabel('Peak height','fontsize',fsize)
    
    clear axes
    
    %choose one of two print options, depending on if you have export_fig
    %installed
    
    %option 1
    export_fig(['Co_additions_' num2str(addlist{addst-4}) '_' num2str(addlist{addst-1})...
        '.pdf']) % Need to have ghostscript and exportfig installed in matlab in order for this to work

    %option 2
%     set(gcf,'PaperPositionMode','auto')
%     print(gcf,'-dpng','-r300',['Co_additions_' num2str(addlist{addst-4}) '_' num2str(addlist{addst-1})...
%         '.png'])%alternative, if export_fig is not installed


    clf;
end

fclose(fid);
cd(home);

disp(['Time elapsed = ' num2str(toc) ' seconds.']);
end

function [height,pln]=peakcalc(fname,dt,fid,ncol)
% set up colors for plot
% separate for sample and additions
% samples are shades of green
% additions follow heatmap -> cool for low addition, warm for high addition
clr=[229,245,249;...
    153,216,201;...
    44,162,95;...
    44,123,182;...
    171,217,233;...
    253,174,97;...
    215,25,28]/255;

% set up parameters for height calculation
% narrows search to area near the known peak
xlims=[-1.25,-1.06];
data=dlmread(fname,'\t',1,0);

ind=data(:,1)>=xlims(1) & data(:,1)<=xlims(2);
x=data(ind,1);
y=data(ind,2);

% smooth the data using SG smooth
ys=smooth(y,5,'sgolay',3);

% calculate derivative
dx=diff(x); % change in x between data points
dy=diff(ys); % change in y between data points

dydx=dy./dx; % first derivative
x=x(1:end-1,1); % adjust size of vector to match derivative (remove 1 with differentiation)
y=y(1:end-1,1); % adjust vector
ys=ys(1:end-1,1); % adjust vector

% smoothed derivative
dydxs=smooth(dydx,7,'sgolay',5);

% plot data with derivative for testing
plot(x,y,'.','markersize',2,'color',clr(ncol,:))
pln(ncol)=plot(x,ys,'k-','color',clr(ncol,:));


% find peak start
% point before negative first derivative
peak=zeros(1,2); % 1.end 2.start

ind=find(dydxs<0,1,'last');
peak(1)=ind+1; % end
ind=find(dydxs>0,1,'first');
peak(2)=ind;

% plot endpoints and baseline
hold on
plot(x(peak),y(peak),'--o','linewidth',1,'markerfacecolor',...
    clr(ncol,:),'markeredgecolor','k','color',clr(ncol,:),'markersize',3)
set(gca,'xlim',xlims)

% interpolate points along baseline (linear)
ind=x>=x(peak(1)) & x<=x(peak(2));
xint=x(ind); % points along baseline
yval=y(ind); % observations between peak start and end
f=fit(x(peak),y(peak),'poly1');
yint=feval(f,xint);

% calculate maximum height relative to baseline
[height,~]=max(-yval+yint);

% write to file
fprintf(fid,'%s,%s,%1.3e\n',dt,fname,height);


end

