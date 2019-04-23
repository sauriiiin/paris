%%  Sau MATLAB Colony Analyzer Toolkit
%
%%  TECHNICAL REPLICATE BASED CONTROL DISTRIBUTION
%   Control Distribution of Validation Experiments (4Control)

%   Author: Saurin Parikh, April 2019
%   dr.saurin.parikh@gmail.com

%%  Load Paths to Files and Data

    col_analyzer_path = '/Users/saur1n/Documents/GitHub/Matlab-Colony-Analyzer-Toolkit';
    bean_toolkit_path = '/Users/saur1n/Documents/GitHub/bean-matlab-toolkit';
    sau_toolkit_path = '/Users/saur1n/Documents/GitHub/sau-matlab-toolkit';
    addpath(genpath(col_analyzer_path));
    addpath(genpath(bean_toolkit_path));
    addpath(genpath(sau_toolkit_path));
%     javaaddpath(uigetfile());

%%  Add MCA toolkit to Path
%     add_mca_toolkit_to_path

%%  Initialization

%     Set preferences with setdbprefs.
    setdbprefs('DataReturnFormat', 'structure');
    setdbprefs({'NullStringRead';'NullStringWrite';'NullNumberRead';'NullNumberWrite'},...
                  {'null';'null';'NaN';'NaN'})

    expt_name = '4C3_GA1';
    expt = 'FS1-1';
%     out_path = '/home/sbp29/MATLAB/4C3_Data/GA1/contdist/';
    out_path = '/Users/saur1n/Desktop/4C3/Analysis/GA/S1Analysis/contdist/';
    density = 6144;

%   MySQL Table Details  

%     tablename_norm      = sprintf('%s_%d_NORM',expt_name,density);
    tablename_fit       = sprintf('%s_%d_FITNESS',expt_name,density);
%     tablename_pval       = sprintf('%s_%d_PVALUE',expt_name,density);

    tablename_p2o       = '4C3_pos2orf_name1';
    tablename_bpos      = '4C3_borderpos';
    
%   Reference Strain Name
    cont.name           = 'BF_control';

%   MySQL Connection and fetch initial data
    connectSQL;
    
    hours = fetch(conn, sprintf(['select distinct hours from %s ',...
                 'order by hours asc'], tablename_fit));
    hours = hours.hours;

%%  CREATING DISTRIBUTION PLOTS

    contpos = fetch(conn, sprintf(['select pos from %s ',...
        'where orf_name = ''%s'' and pos < 10000 ',...
        'and pos not in ',...
        '(select pos from %s)'],...
        tablename_p2o,cont.name,tablename_bpos));
    contpos = contpos.pos + [110000,120000,130000,140000,...
        210000,220000,230000,240000];
    yy = 0:0.1:25;

    for iii = 7:length(hours)
%         contfit = [];
%         for ii = 1:length(contpos)
%             temp = fetch(conn, sprintf(['select fitness from %s ',...
%                 'where hours = %d and pos in (%s) ',...
%                 'and fitness is not null'],tablename_fit,hours(iii),...
%                 sprintf('%d,%d,%d,%d,%d,%d,%d,%d',contpos(ii,:))));
%             if nansum(temp.fitness) > 0
%                 contfit = [contfit, nanmean(temp.fitness)];
%             end
%         end
% 
%         contfmean = nanmean(contfit);
%         contfmed = nanmedian(contfit);
%         contfstd = nanstd(contfit);
%         
%         perc25f = prctile(contfit,2.5);
%         perc975f = prctile(contfit,97.5);
% 
%         fig = figure('Renderer', 'painters', 'Position', [10 10 960 600],'visible','off');
% %         figure()
%         [f,xi] = ksdensity(contfit);
%         plot(xi,f,'LineWidth',3)
%         hold on
%         plot(ones(1,length(yy))*perc25f,yy,'--r',...
%             ones(1,length(yy))*perc975f,yy,'--r',...
%             ones(1,length(yy))*contfmed,yy,'--r','LineWidth',1)
%         xlim([0.9,1.15])
%         ylim([0,25])
%         hold off
%         grid on
%         grid minor
%         xlabel('Fitness')
%         ylabel('Density')
%         title(sprintf(['%s | Control Distribution (Fitness) | Time = %d hrs\n',...
%             'Perc 2.5 = %0.4f | Mean = %0.4f | Perc 50 = %0.4f | Perc 97.5 = %0.4f'],...
%             expt,hours(iii),perc25f,contfmean,contfmed,perc975f))
%         saveas(fig,sprintf('%s%s_ContDist_%d.png',...
%                     out_path,expt_name,hours(iii)))

        contavg = [];
        for ii = 1:size(contpos,1)
            temp = fetch(conn, sprintf(['select average from %s ',...
                'where hours = %d and pos in (%s) ',...
                'and fitness is not null'],tablename_fit,hours(iii),...
                sprintf('%d,%d,%d,%d,%d,%d,%d,%d',contpos(ii,:))));
            if nansum(temp.average) > 0
                contavg = [contavg, temp.average];
            end
        end
        
        xmin = round(min(min(contavg)) - 50, -1);
        xmax = round(max(max(contavg)) + 50, -1);

        for i = 1:size(contpos,2)
            contamean = nanmean(contavg(i,:));
            contamed = nanmedian(contavg(i,:));
            contastd = nanstd(contavg(i,:));

            perc25a = prctile(contavg(i,:),2.5);
            perc975a = prctile(contavg(i,:),97.5);

    %         fig = figure('Renderer', 'painters', 'Position', [10 10 960 600],'visible','off');
            figure()
            [f,xi] = ksdensity(contavg(i,:));
            plot(xi,f,'LineWidth',3)
            hold on
            plot(ones(1,length(yy))*perc25a,yy,'--r',...
                ones(1,length(yy))*perc975a,yy,'--r',...
                ones(1,length(yy))*contamed,yy,'--r','LineWidth',1)
            xlim([xmin,xmax])
            ylim([0,0.02])
            hold off
            grid on
            grid minor
            xlabel('Pixel Count')
            ylabel('Density')
            title(sprintf(['%s | Control Distribution (Pixel Count) | Time = %d hrs\n',...
                'Perc 2.5 = %0.4f | Mean = %0.4f | Perc 50 = %0.4f | Perc 97.5 = %0.4f'],...
                expt,hours(iii),perc25a,contamean,contamed,perc975a))
    %         saveas(fig,sprintf('%s%s_ContDistPix_%d.png',...
    %                     out_path,expt_name,hours(iii)))
        end
    end
    
%%  END
