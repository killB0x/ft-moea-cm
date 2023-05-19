clc, clear all, close all


%% Initial parameters:

case_study = {'monopropellant_propulsion_system_noise_0per'};
save_plot = 1;

for cccc = 1:length(case_study)
    
    save_plots_dir = strcat('',case_study{cccc},'');
    if ~isfolder(save_plots_dir)
        mkdir(save_plots_dir)
    end
    
    
    check_folders = {
        strcat(['data_for_figures/fig8p/',case_study{cccc},'_MCS0TS0Acc1_sfv0']),...
        strcat(['data_for_figures/fig8/',case_study{cccc},'_MCS0TS0Acc1_sfv1']),...
        strcat(['data_for_figures/fig8/',case_study{cccc},'_MCS0TS0Acc1_sfv2']),...
        strcat(['data_for_figures/fig8/',case_study{cccc},'_MCS0TS0Acc1_sfv3']),...
        strcat(['data_for_figures/fig8/',case_study{cccc},'_MCS0TS0Acc1_sfv4']),...
        strcat(['data_for_figures/fig8/',case_study{cccc},'_MCS0TS0Acc1_sfv5']),...
        strcat(['data_for_figures/fig8/',case_study{cccc},'_MCS0TS0Acc1_sfv6']),...
        strcat(['data_for_figures/fig8p/',case_study{cccc},'_MCS1TS1Acc1_sfv0']),...
        strcat(['data_for_figures/fig8/',case_study{cccc},'_MCS1TS1Acc1_sfv1']),...
        strcat(['data_for_figures/fig8/',case_study{cccc},'_MCS1TS1Acc1_sfv2']),...
        strcat(['data_for_figures/fig8/',case_study{cccc},'_MCS1TS1Acc1_sfv3']),...
        strcat(['data_for_figures/fig8/',case_study{cccc},'_MCS1TS1Acc1_sfv4']),...
        strcat(['data_for_figures/fig8/',case_study{cccc},'_MCS1TS1Acc1_sfv5']),...
        strcat(['data_for_figures/fig8/',case_study{cccc},'_MCS1TS1Acc1_sfv6']),...
        };

    n = {'0',...
        '1',...
        '2',...
        '3',...
        '4',...
        '5',...
        '6',...
        '0',...
        '1',...
        '2',...
        '3',...
        '4',...
        '5',...
        '6',...
        };
    
    
    m = {'d','d','d','d','d','d','d','sdc','sdc','sdc','sdc','sdc','sdc','sdc'};
    
    plot_setup.fig_pos    = [ -281   177   210   193];
    plot_setup.axis_pos   = [0.3019    0.2812    0.6462    0.6786];
    plot_setup.text_size  = 16;
    plot_setup.legend_pos = [0.2142    0.7810    0.7230    0.2059];
    plot_setup.font_name  = 'Arial Unicode MS';
    plot_setup.markers = 'oxs^+';
    plot_setup.ylabel = '\rho';
    
    
    compare = {};
    output_results = {};
    compare_all = {};
    for i = 1 : length(check_folders)
        output_results = [output_results;get_all_results_output(get_dir_files_str(check_folders{i}))];
        compare = [compare;get_all_results(get_dir_files_str(check_folders{i}))];
        compare_all = [compare_all;get_all_all_results(get_dir_files_str(check_folders{i}))];
    end
    
    %%
    pop_size_all = [400];
    for pop = 1:length(pop_size_all)
        pop_size = pop_size_all(pop);
        r.tree_size = [];
        r.acc_data = [];
        r.acc_mcs = [];
        r.t = [];
        r.num_super_var = [];
        r.groupByColor = {};
        r.case_study = {};
        for i = 1 : length(check_folders)
            idx = find(output_results{i}.pop_size == pop_size);
            r.tree_size     = [r.tree_size;double(output_results{i}.tree_size(idx))];
            r.acc_data      = [r.acc_data;double(output_results{i}.acc_after_noise(idx))];
            r.acc_mcs       = [r.acc_mcs;double(output_results{i}.acc_mcs_gt(idx))];
            r.num_super_var = [r.num_super_var; double(output_results{i}.num_super_var(idx))];
            r.t = [r.t;double(output_results{i}.time_convergence(idx))/60];
            %     for j = 1 : length(idx)
            %         r.t = [r.t;compare{i}(idx(j)).time(end)/60];
            %     end
            r.groupByColor = [r.groupByColor;repmat({m{i}},length(idx),1,1)];
            r.case_study = [r.case_study;repmat({n{i}},length(idx),1,1)];
        end
        
        % Tree size
        figure('color',[1 1 1],'position',[plot_setup.fig_pos]),hold on
        boxchart(categorical(r.case_study,unique(r.case_study)),r.tree_size,'GroupByColor',r.groupByColor)
        %leg1 = legend('show');
        %set(leg1,'position',plot_setup.legend_pos)
        box on,% grid on
        set(gca,'FontName',plot_setup.font_name,'FontSize',plot_setup.text_size,'position',plot_setup.axis_pos,'FontAngle','normal','YTick',[0:25:100])
        plot([0 max(double(get(gca,'XTick')))+1],[23 23],'r-.')
        leg1 = legend(unique(m));
        set(leg1,'location','best')
        ylabel('\phi_s')
        xlabel(plot_setup.ylabel)
        if save_plot == 1
            h = gcf;
            set(h,'Units','Inches');
            pos = get(h,'Position');
            set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
            print(h,strcat(['',save_plots_dir,'/pop_',num2str(pop_size),'_tree_size']),'-dpdf','-r0')
            savefig(gcf,strcat(['',save_plots_dir,'/pop_',num2str(pop_size),'_tree_size']))
            close(gcf)
        end
        
        % Acc. dataset:
        figure('color',[1 1 1],'position',[plot_setup.fig_pos])
        boxchart(categorical(r.case_study,unique(r.case_study)),1-r.acc_data,'GroupByColor',r.groupByColor)
        leg1 = legend('show');
        set(leg1,'location','best')
        %set(leg1,'position',plot_setup.legend_pos)
        box on%, grid on
        set(gca,'FontName',plot_setup.font_name,'FontSize',plot_setup.text_size,'position',plot_setup.axis_pos,'FontAngle','normal','YTick',[0:0.1:0.4])
        ylabel('\phi_d')
        ylim([-0.01 0.4])
        xlabel(plot_setup.ylabel)
        if save_plot == 1
            h = gcf;
            set(h,'Units','Inches');
            pos = get(h,'Position');
            set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
            print(h,strcat(['',save_plots_dir,'/pop_',num2str(pop_size),'_acc_data']),'-dpdf','-r0')
            savefig(gcf,strcat(['',save_plots_dir,'/pop_',num2str(pop_size),'_acc_data']))
            close(gcf)
        end
        
        % Acc. MCSs:
        figure('color',[1 1 1],'position',[plot_setup.fig_pos])
        boxchart(categorical(r.case_study,unique(r.case_study)),1-r.acc_mcs,'GroupByColor',r.groupByColor)
        leg1 = legend('show');
        set(leg1,'location','best')
        %set(leg1,'position',plot_setup.legend_pos)
        box on%, grid on
        set(gca,'FontName',plot_setup.font_name,'FontSize',plot_setup.text_size,'position',plot_setup.axis_pos,'FontAngle','normal','YTick',[0:0.1:0.4])
        ylabel('\phi_c')
        ylim([-0.01 0.4])
        xlabel(plot_setup.ylabel)
        if save_plot == 1
            h = gcf;
            set(h,'Units','Inches');
            pos = get(h,'Position');
            set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
            print(h,strcat(['',save_plots_dir,'/pop_',num2str(pop_size),'_acc_mcs']),'-dpdf','-r0')
            savefig(gcf,strcat(['',save_plots_dir,'/pop_',num2str(pop_size),'_acc_mcs']))
            close(gcf)
        end
        
        % Time:
        figure('color',[1 1 1],'position',[plot_setup.fig_pos])
        boxchart(categorical(r.case_study,unique(r.case_study)),r.t,'GroupByColor',r.groupByColor)
        %set(leg1,'position',plot_setup.legend_pos)
        box on%, grid on
        set(gca,'FontName',plot_setup.font_name,'FontSize',plot_setup.text_size,'position',plot_setup.axis_pos,'FontAngle','normal','YTick',[0:100:400])
        ylabel('Time (min)')
        xlabel(plot_setup.ylabel)
        leg1 = legend('show');
        set(leg1,'location','best')
        ylim([0 400])
        
        if save_plot == 1
            h = gcf;
            set(h,'Units','Inches');
            pos = get(h,'Position');
            set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
            print(h,strcat(['',save_plots_dir,'/pop_',num2str(pop_size),'_time']),'-dpdf','-r0')
            savefig(gcf,strcat(['',save_plots_dir,'/pop_',num2str(pop_size),'_time']))
            close(gcf)
        end

        % Number of superfluous variables:
        figure('color',[1 1 1],'position',[plot_setup.fig_pos])
        boxchart(categorical(r.case_study,unique(r.case_study)),r.num_super_var,'GroupByColor',r.groupByColor)
        %set(leg1,'position',plot_setup.legend_pos)
        box on%, grid on
        set(gca,'FontName',plot_setup.font_name,'FontSize',plot_setup.text_size,'position',plot_setup.axis_pos,'FontAngle','normal','YTick',[-6:2:6])
        ylabel('\pm BEs')
        xlabel(plot_setup.ylabel)
        ylim([-6 6])
        leg1 = legend('show');
        set(leg1,'location','best')
        
        if save_plot == 1
            h = gcf;
            set(h,'Units','Inches');
            pos = get(h,'Position');
            set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
            print(h,strcat(['',save_plots_dir,'/pop_',num2str(pop_size),'_superfluous_var']),'-dpdf','-r0')
            savefig(gcf,strcat(['',save_plots_dir,'/pop_',num2str(pop_size),'_superfluous_var']))
            close(gcf)
        end
        
        
    end
end
%%
function all_files= get_dir_files_str(folder_results)
%% Get the parent directories:
files = dir(folder_results);
files = {files.name}';
a = [];
for i = 1:length(files)
    if files{i}(1) ~= 'c'
        a = [a;i];
    end
end
files(a) = [];
g = [];
for i = 1 : length(files)
    a = strfind(files{i},'_');
    g = [g;str2num(files{i}(a(end)+1:end))];
end
[~,x] = sort(g,'ascend');
s = g(x);
files = files(x);
%% Working directories:
all_files = [];
for i = 1 : length(files)
    dir_path = strcat(['',folder_results,'/',files{i},'/run_0']);
    subfiles = dir(dir_path);
    subfiles = {subfiles.name}';
    a = [];
    for j = 1:length(subfiles)
        if subfiles{j}(1) ~= 'g'
            a = [a;j];
        end
    end
    subfiles(a) = [];
    g = [];
    for j = 1 : length(subfiles)
        a = strfind(subfiles{j},'_');
        g = [g;str2num(subfiles{j}(a(1)+1:end-4))];
    end
    [~,x] = sort(g,'ascend');
    subfiles = subfiles(x);
    %subfiles
    all_files(i).pop_size = s(i);
    all_files(i).path = dir_path;
    all_files(i).files = subfiles;
    idx = strfind(folder_results,'/');
    all_files(i).case_study = folder_results(idx(end)+1:end);
end

end
function r = get_all_results(c)
for j = 1 : length(c)
    a = c(j);
    pop_size = a.pop_size;
    acc_data = [];
    tree_size =[];
    time = [];
    acc_mcs=[];
    str_fts = {};
    for i = 1 : length(a.files)
        b = load(strcat(['',a.path,'/',a.files{i},'']));
        acc_data = [acc_data;b.parameters_trimmed(end,3)];
        time = [time;b.time];
        tree_size =[tree_size;b.parameters_trimmed(end,2)];
        acc_mcs = [acc_mcs;b.parameters_trimmed(end,1)];
        str_fts = [str_fts;b.fts_sorted_trimmed];
    end

    % Remove part of the convergence time:
    %temp1 = [acc_data,acc_mcs,tree_size];
    %idx = size(temp1,1) - find(flip(sum(diff(temp1)'))~=0,1) + 1;
    idx = length(time);
    %if isempty(idx)
    %   idx = 1; 
    %end

    r(j).time = time(1:idx);
    r(j).tree_size = tree_size(1:idx);
    r(j).acc_data = acc_data(1:idx);
    r(j).acc_mcs = acc_mcs(1:idx);
    r(j).pop_size = pop_size;
    r(j).case_study = a.case_study;
    r(j).str_fts = str_fts(1:idx);
    
end
end
function r = get_all_results_output(c)
r = [];
noise = [];
pop_size = [];
acc_before_noise = [];
acc_after_noise = [];
acc_after_cleaning = [];
tree_size = [];
time_convergence = [];
acc_mcs_gt = [];
num_super_var = [];
for j = 1 : length(c)
    a = c(j);
    b = load(strcat(['',a.path,'/output.mat']));
    
    noise = [noise;b.noise];
    pop_size = [pop_size;c(j).pop_size];
    acc_before_noise = [acc_before_noise;b.acc_data_before_noise];
    acc_after_noise = [acc_after_noise;b.acc_data_after_noise];
    %acc_after_cleaning = [acc_after_cleaning;b.acc_data_after_cleaning];
    time_convergence = [time_convergence;b.conv_time];
    acc_mcs_gt = [acc_mcs_gt;b.acc_mcs];
    tree_size = [tree_size;b.tree_size];
    num_super_var = [num_super_var;length(count_unique_BEs(b.best_ft)) - length(count_unique_BEs(b.original_ft))];
    
end

r.noise = noise;
r.pop_size = pop_size;
r.acc_before_noise = acc_before_noise;
r.acc_after_noise = acc_after_noise;
%r.acc_after_cleaning = acc_after_cleaning;
r.time_convergence = time_convergence;
r.acc_mcs_gt = acc_mcs_gt;
r.tree_size = tree_size;
r.num_super_var = num_super_var;
try
    r.case_study = a.case_study;
catch
    r
end

end
function r = get_all_all_results(c)
for j = 1 : length(c)
    a = c(j);
    pop_size  = a.pop_size;
    acc_data  = nan(a.pop_size,length(a.files));
    tree_size = nan(a.pop_size,length(a.files));
    time      = nan(a.pop_size,length(a.files));
    acc_mcs   = nan(a.pop_size,length(a.files));
    str_fts   = {};
    for i = 1 : length(a.files)
        b = load(strcat(['',a.path,'/',a.files{i},'']));
        acc_data(:,i)  = b.parameters_trimmed(:,3);
        time(:,i)      = b.time;
        tree_size(:,i) = b.parameters_trimmed(:,2);
        acc_mcs(:,i)   = b.parameters_trimmed(:,1);
        str_fts        = [str_fts,num2cell(b.fts_sorted_trimmed,2)];
    end

    r(j).time = time;
    r(j).tree_size = tree_size;
    r(j).acc_data = acc_data;
    r(j).acc_mcs = acc_mcs;
    r(j).pop_size = pop_size;
    r(j).case_study = a.case_study;
    r(j).str_fts = str_fts;
end
end
%%
%% 
function BEs_all = count_unique_BEs(tree_string)
% Identify the total number of BEs in this FT:
BE_idx = strfind(tree_string,'BE');
c_idx = strfind(tree_string,',');
p_idx = strfind(tree_string,')');
BEs_all = {};
for i = 1 : length(BE_idx)
    c1 = p_idx(find(p_idx>BE_idx(i),1));
    c2 = c_idx(find(c_idx>BE_idx(i),1));
    c = min([c1 c2]);
    BEs_all = [BEs_all,tree_string(BE_idx(i):c-1)];
end
BEs_all = unique(BEs_all)';
end
