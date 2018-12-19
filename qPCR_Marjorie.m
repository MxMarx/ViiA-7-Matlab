folder = 'V:\Russell\qPCR\Marjorie\12.18.18';
fname = fullfile(folder,'2018-12-17 PACAP Cre Mice HM3DQ RiboTag cFos GAPDH.xls');
keyname = fullfile(folder,'2018Dec17_PACAP cre mice cDNA Analysis Plate Layout.xlsx');


% opts = detectImportOptions(fname,'NumHeaderLines',35,'Sheet','Results');
% opts = setvartype(opts,opts.VariableNames([5,6]),'categorical');

% save('ImportOpts','opts');
load('ImportOpts','opts');

t = readtable(fname,opts,'basic',1);
k2 = readtable(keyname,'Sheet','Key','basic',1);
k = readtable(keyname,'Sheet','Plate Template','Range','C4:Z20','basic',1);

for i = 1:width(k)
    k.(i) = categorical(k{:,i});
end

%% Get the key
sampleID = categorical(reshape(k{:,:}',384,1));
t.SampleID = sampleID(t.Well);


% Find samples outside of expected TM
TM_Threshold = inf;
GroupingVariables = {'TargetName'};
InputVariables = {'Tm1'};
OutputVariableNames = {'WithinTM'};
a = rowfun(@(x) abs(x - median(x)) < TM_Threshold,t,'GroupingVariables',GroupingVariables,'InputVariables',InputVariables,'OutputVariableNames',OutputVariableNames);
t.WithinTM = a.WithinTM;

% Find the standards
standard_concentration = cellfun(@(x) strsplit(x,' '),t.SampleName,'UniformOutput',0);
standard_concentration = cellfun(@(x) max([0,str2num(x{1})]),standard_concentration);

% b = find(standard_concentration);
% 
% [~,genes] = findgroups(t{b,'TargetName'});
% %% Calculate efficiency add rtsq to table
% fprintf(1,'\n\n%-20s%-15s%-15s\n','Gene:','Efficiency:','R-Squared:');
% for i = genes'
%     standards = t.TargetName == i & standard_concentration>0 & ~isnan(t.CT) & t.WithinTM;
%     concentration = standard_concentration(standards);
%     CT = t{standards,'CT'};
%     [beta,stats] = robustfit(CT,log(concentration));
%     
%     efficiency = exp( -1 / beta(2));
%     rsquared = corr(CT,beta(1)+beta(2)*log(concentration))^2;
%     
%     fprintf(1,'%-20s%-15.4f%-15.4f \n',i,efficiency,rsquared);
%     
%     % calculate rtsq
%     CT = t{t.TargetName==i,'CT'};
%     rstq = exp((CT - beta(1)) / beta(2)); 
%     t{t.TargetName==i,'rstq'} = rstq;
%     
%     % calculate rtsq, assuming efficiency is 2
%     rstq = log(exp((CT - beta(1)) / 2)); 
%     t{t.TargetName==i,'rstq_efficiency_2'} = rstq;
% end
close all


standards_cutoff = -3;
%% Calculate rstq
% This does the same thing as gramm
t.log2_standard_concentration = log2(standard_concentration);
c = t(t.log2_standard_concentration >= standards_cutoff,:);
mdl = fitlme(c,'log2_standard_concentration ~ CT*TargetName')
t.rstq = 2.^predict(mdl,t);

%% Plot Standards

ip = contains(cellstr(t.SampleID),{'IP'},'IgnoreCase',true);
in = contains(cellstr(t.SampleID),'IN','IgnoreCase',true);
sample = categorical(ip+in*2,[1,2],{'IP','IN'});

figure('Position',[100,100,700,500])
clear g
if any(~isundefined(sample))
    g(1,1) = gram('x',t.CT,'color',t.TargetName,'subset',~isundefined(sample));
    g(1,1).facet_grid(sample,[]);
else
    g(1,1) = gram('x',t.CT,'color',t.TargetName,'subset',standard_concentration==0);
end
g(1,1).stat_bin('geom','overlaid_bar')
g(1,1).set_names('y','Count','x','','color','','row','');
g(1,1).no_legend
g(1,1).set_layout_options('position',[0,.6,.8,.4],'legend',0)

g(2,1) = gramm('y',t.log2_standard_concentration,'x',t.CT,'color',t.TargetName,'subset',t.log2_standard_concentration >= -10);
g(2,1).set_color_options('lightness_range',[30,80]);
g(2,1).geom_point();
g(2,1).set_names('y','Log2 concentration','x','CT','color','Target','lightness','Tm1 out of range');
g(2,1).set_stat_options('alpha',.3173);
g(2,1).set_limit_extra([.1,.1],[.1,.1]);
g(2,1).axe_property('XTick',0:1:100,'YTick',-10:1:10,'XGrid','on','Ygrid','on','GridColor',[0.5 0.5 0.5])
g(2,1).set_layout_options('position',[0,0,1,.6])

g.draw;
g(2,1).update('subset',t.log2_standard_concentration >= standards_cutoff)
g(2,1).stat_glm('fullrange',1);
% g(2,1).axe_property('Ylim',[-4.1,-.9])
g(2,1).no_legend
g.draw

xl = [min(g(1).facet_axes_handles(1).XLim(1),g(2).facet_axes_handles.XLim(1)),...
    max(g(1).facet_axes_handles(1).XLim(2),g(2).facet_axes_handles.XLim(2))];
    
% xl = [min(t.CT)-1, max(t.CT)+1];
g(1).facet_axes_handles(1).XLim = xl;
g(2).facet_axes_handles.XLim = xl;
g(1).facet_axes_handles(1).Position(3) = g(2).facet_axes_handles.Position(3)
g(1).facet_axes_handles(end).Position(3) = g(2).facet_axes_handles.Position(3)

g.export('export_path',folder,'file_name','standards','file_type','png')

%% Heatmaps
rows = 'ABCDEFGHIJKLMNOP';
[I,J] = ind2sub([24,16],t.Well);
t.I = I;
t.J = rows(J)';
for i = {'CT','Tm1'}
figure('position',[100 100 1200 600],'color','w');
heatmap(t,'I','J','ColorVariable',i,'Title',i);
colormap(inferno);
ylabel('Row');
xlabel('Column');
export_fig(fullfile(folder,[i{1} '.png']),'-m3');
end


GroupingVariables = {'SampleID','TargetName'};
InputVariables = {'CT','rstq'};
OutputVariableNames = {'samples','CT','std_CT','rstq','std_rstq'};
func = @(CT,rstq,rstq_efficiency_2) deal(...
    sum(~isnan(CT)),...
    nanmean(CT),...
    nanstd(CT),...
    nanmean(rstq),...
    nanstd(rstq));
a = rowfun(func,t,'GroupingVariables',GroupingVariables,'InputVariables',InputVariables,'OutputVariableNames',OutputVariableNames);

%% Get the other key
k2.Properties.VariableNames(1) = {'SampleID'};
k2.SampleID = categorical(k2.SampleID);
j = innerjoin(k2,a,'Keys','SampleID');

%% Save
j.GroupCount = [];
j = sortrows(j,{'TargetName','SampleID'});
[fname,fpath] = uiputfile(fullfile(folder,'Output.xlsx'));
writetable(j,[fpath,fname]);


figure
g = gram('x',j.type,'y',j.rstq,'group',cellfun(@(x) x(1:6),cellstr(j.SampleID),'UniformOutput',0),'lightness',j.Group,'color',j.Virus)
g.geom_line
g.axe_property('YScale', 'log')
% g.set_limit_extra([.1,.1],[1,1]);

g.facet_grid([],j.TargetName,'scale','free')
g.set_names('x','','y','RSTQ','column','','color','ID','Lightness','Group');
g.set_color_options('lightness_range',[75,25])
g.geom_point
g.draw
g.export('export_path',folder,'file_name','Graphs','file_type','png')





