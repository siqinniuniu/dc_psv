%downward continuate surface plane wavefield by a station on top of the
%ice sheet.

clear all;close all;clc
%% parameters

%sac
sacdir = 'sac/';
saclst = [sacdir,'evt.lst'];
cmpnm = {'.BHR','.BHZ'};

%time samples
t0 = 0;
t1 = 50;
fs = 20;

%only vary thickness of ice sheet
z = 1:0.1:5;
nmod = length(z);

%Earth model (ice/bedrock)
vp0 = [3.87; 5.8]; % km/s
vs0 = [1.95; 3.46]; % km/s
rho0 = [0.917; 2.72]; % g/cm^3
nlyr = 2;

%model space
vp = repmat(vp0,1,nmod);
vs = repmat(vs0,1,nmod);
rho = repmat(rho0,1,nmod);
thik = zeros(nlyr,nmod); thik(1,:) = z;

% sachd
sachd_rayp = 'user0';

%output
sacout = 'sac_dc/';
figdir = 'figure/';

%% load sac data

sacst = SACST_fread('list',saclst,'prefix',sacdir,'suffix',cmpnm);

%% downward continuation

[REDsu,m1,Esu1,m0,Esu0] = SACST_dc_psv_nmodel(...
    sacst,...                      % sac data
    t0,t1,fs,...                   % time window for downward continuation
    nlyr,nmod,vp,vs,rho,thik,...   % earth models 
    sachd_rayp);                   % sac head field for ray parameter

%% output sac

evnm = textread(saclst,'%s');
nevt = size(REDsu,1);
nt = size(m0,2); dt = 1/fs; 

sacst_dc = sacst(:,ones(4,1));
[sacst_dc.b] = deal(t0);
[sacst_dc.delta] = deal(dt);
[sacst_dc.npts] = deal(nt);
[sacst_dc(:,1).kcmpnm] = deal('Pd');
[sacst_dc(:,2).kcmpnm] = deal('Pu');
[sacst_dc(:,3).kcmpnm] = deal('Sd');
[sacst_dc(:,4).kcmpnm] = deal('Su');

for ievt = 1:nevt
    
    %get the model of maximum REDsu
    [c,imod] = max(REDsu(ievt,:));
    
    %data
    idx0 = 4*(ievt-1);
    sacst_dc(ievt,1).data = m1(idx0+1,:,imod); %Pd
    sacst_dc(ievt,2).data = m1(idx0+2,:,imod); %Pu
    sacst_dc(ievt,3).data = m1(idx0+3,:,imod); %Sd
    sacst_dc(ievt,4).data = m1(idx0+4,:,imod); %Su
    [sacst_dc(ievt,:).user3] = deal(z(imod)); %optimal thickness
end

SACST_fwrite(sacst_dc,evnm,'prefix',sacout,...
    'suffix',{'.Pd','.Pu','.Sd','.Su'})

%% Figure

hf = figure;
% hf = figure('visible','off');

% figure size
set(hf,'paperunits','centimeters','units','centimeters',...
    'PaperPositionMode','manual','papertype','a4',...
    'paperorientation','landscape',...
    'position',[0 0 29.7 21],'paperposition',[0 0 29.7 21])

% paper margins
West = 2;
South = 2;
Width = 29.7-4;
Height = 21-5;

% define axes positions
% Pu
pos_ax1 = [West,South,Width*0.4,Height];
hax1 = axes('Units','centimeters','position',pos_ax1);
% Su
pos_ax2 = [West+Width*0.45,South,Width*0.4,Height];
hax2 = axes('Units','centimeters','position',pos_ax2);
% REDsu
pos_ax3 = [West+Width*0.9,South,Width*0.15,Height];
hax3 = axes('Units','centimeters','position',pos_ax3);
% title
pos_ax4 = [West+Width/2,South+Height+1.5,0.1,0.1];
hax4 = axes('Units','centimeters','position',pos_ax4,...
    'visible','off','Tag','suptitle');

%% plot each event

t = t0+(0:nt-1)*dt;
v_xlim = [t0,t1];

for ievt = 1:nevt
    
    idx0 = 4*(ievt-1);
    
    % plot Pu
    set(hf,'currentaxes',hax1)
    seistrace = squeeze(m1(idx0+2,:,:));
    yoffset = 1:nmod;
    maxamp = mean(max(abs(seistrace)));
    yscale = 1.5/maxamp;
    data = yscale*seistrace+yoffset(ones(nt,1),:);
    v_ylim = [min(data(:))-0.5,max(data(:))+0.5];
    plot(t,data,'k')
    xlim(v_xlim); ylim(v_ylim)
    v_ytick = yoffset;
    str_ytick = cellstr(num2str(z(:),'%3.1f'));
    set(gca,'ytick',v_ytick,'yticklabel',str_ytick)
    title('upgoing P in bedrock')
    xlabel('Time (s)'); ylabel('Thickness of ice sheet (km)')
    
    %plot Su
    set(hf,'currentaxes',hax2)
    seistrace = squeeze(m1(idx0+4,:,:));
    maxamp = mean(max(abs(seistrace)));
    yscale = 1.5/maxamp;
    data = yscale*seistrace+yoffset(ones(nt,1),:);
    plot(t,data,'k')
    xlim(v_xlim); ylim(v_ylim)
    set(gca,'ytick',v_ytick,'yticklabel',str_ytick)
    title('upgoing S in bedrock')
    xlabel('Time (s)')
    
    %plot REDsu
    set(hf,'currentaxes',hax3)
    data = REDsu(ievt,:);
    plot(data,yoffset,'ko')
    ylim(v_ylim); xlim([min(data) 1])
    set(gca,'xdir','reverse')
    set(gca,'ytick',v_ytick,'yticklabel',str_ytick)
    title('1-Esu1/Esu0')
    xlabel('Energy reduction rate')
    
    %title
    set(hf,'currentaxes',hax4)
    str_title = evnm{ievt};
    title(str_title)
    htext = text(0,0,str_title);
    set(htext,'horizontalalignment','center',...
        'verticalalignment','bottom',...
        'interpreter','none');
    
    %font size
    formatfig(hf)
    
    %output figure
    fignm = sprintf('%s/%s.pdf',figdir,evnm{ievt});
    print(hf,'-dpdf','-painters',fignm)
    
    %clear for next plot
    delete(htext)
   
end