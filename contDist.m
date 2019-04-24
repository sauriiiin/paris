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

%   Set preferences with setdbprefs.
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
    
%     splt = {'1TL','1TR','1BL','1BR','2TL','2TR','2BL','2BR'};
    splt = {'TL','TR','BL','BR'};
    c = categorical({'< mean','> mean'});
    
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
        contavg = [];
        contfit = [];
        for ii = 1:size(contpos,1)
            temp = fetch(conn, sprintf(['select average, fitness from %s ',...
                'where hours = %d and pos in (%s) ',...
                'and fitness is not null'],tablename_fit,hours(iii),...
                sprintf('%d,%d,%d,%d,%d,%d,%d,%d',contpos(ii,:))));
            if nansum(temp.average) > 0
                contavg = [contavg, temp.average];
                contfit = [contfit, temp.fitness];
            end
        end
        
        for i = 1:size(contpos,2)
            contamean{i} = nanmean(contavg(i,:));      contfmean{i} = nanmean(contfit(i,:));
            contamed{i} = nanmedian(contavg(i,:));     contfmed{i} = nanmedian(contfit(i,:));
            contastd{i} = nanstd(contavg(i,:));        contfstd{i} = nanstd(contfit(i,:));
            perc25a{i} = prctile(contavg(i,:),2.5);    perc25f{i} = prctile(contfit(i,:),2.5);
            perc975a{i} = prctile(contavg(i,:),97.5);  perc975f{i} = prctile(contfit(i,:),97.5);
            [fa{i},xia{i}] = ksdensity(contavg(i,:));
            [ff{i},xif{i}] = ksdensity(contfit(i,:));
        end
        
        xmina = round(min(min(contavg))-50,-1);
        xmaxa = round(max(max(contavg))+50,-1);
        xminf = round(min(min(contfit)),3);
        xmaxf = round(max(max(contfit)),3);
        
        for p = 1:2
            fig = figure('Renderer','painters',...
                'Position',[10 10 2000 2000],...
                'visible','off');
%             figure()
            subplot(4,4,1)
            plot(xia{p-1+1},fa{p-1+1},'LineWidth',3)
            hold on
            plot(ones(1,length(yy))*perc25a{p-1+1},yy,'--r',...
                ones(1,length(yy))*perc975a{p-1+1},yy,'--r',...
                ones(1,length(yy))*contamed{p-1+1},yy,'--r','LineWidth',1)
            xlim([xmina,xmaxa])
            ylim([0,0.02])
            hold off
            grid on
            grid minor
            ylabel(sprintf('%s\nDensity',splt{1}))
            title(sprintf('Perc 2.5 = %0.2f | Mean = %0.2f | Median = %0.2f | Perc 97.5 = %0.2f',...
                perc25a{p-1+1},contamean{p-1+1},contamed{p-1+1},perc975a{p-1+1}))
            
            subplot(4,4,2)
            counts = [sum(contavg(p-1+1,:) < contamean{p-1+1}),...
                sum(contavg(p-1+1,:) > contamean{p-1+1})];
            b = bar(c,counts,'FaceColor','flat');
            b.CData(1,:) = [97/255,97/255,97/255];
            b.CData(2,:) = [96/255,125/255,139/255];
            grid on
            grid minor
            ylim([0,300])
            ylabel('Count')

            subplot(4,4,3)
            plot(xif{p-1+1},ff{p-1+1},'LineWidth',3)
            hold on
            plot(ones(1,length(yy))*perc25f{p-1+1},yy,'--r',...
                ones(1,length(yy))*perc975f{p-1+1},yy,'--r',...
                ones(1,length(yy))*contfmed{p-1+1},yy,'--r','LineWidth',1)
            xlim([xminf,xmaxf])
            ylim([0,15])
            hold off
            grid on
            grid minor
            title(sprintf('Perc 2.5 = %0.2f | Mean = %0.2f | Median = %0.2f | Perc 97.5 = %0.2f',...
                perc25f{p-1+1},contfmean{p-1+1},contfmed{p-1+1},perc975f{p-1+1}))
            
            subplot(4,4,4)
            counts = [sum(contfit(p-1+1,:) < contfmean{p-1+1}),...
                sum(contfit(p-1+1,:) > contfmean{p-1+1})];
            b = bar(c,counts,'FaceColor','flat');
            b.CData(1,:) = [97/255,97/255,97/255];
            b.CData(2,:) = [96/255,125/255,139/255];
            grid on
            grid minor
            ylim([0,300])
            ylabel('Count')
            
            subplot(4,4,5)
            plot(xia{p-1+2},fa{p-1+2},'LineWidth',3)
            hold on
            plot(ones(1,length(yy))*perc25a{p-1+2},yy,'--r',...
                ones(1,length(yy))*perc975a{p-1+2},yy,'--r',...
                ones(1,length(yy))*contamed{p-1+2},yy,'--r','LineWidth',1)
            xlim([xmina,xmaxa])
            ylim([0,0.02])
            hold off
            grid on
            grid minor
            ylabel(sprintf('%s\nDensity',splt{2}))
            title(sprintf('Perc 2.5 = %0.2f | Mean = %0.2f | Median = %0.2f | Perc 97.5 = %0.2f',...
                perc25a{p-1+2},contamean{p-1+2},contamed{p-1+2},perc975a{p-1+2}))
            
            subplot(4,4,6)
            counts = [sum(contavg(p-1+2,:) < contamean{p-1+2}),...
                sum(contavg(p-1+2,:) > contamean{p-1+2})];
            b = bar(c,counts,'FaceColor','flat');
            b.CData(1,:) = [97/255,97/255,97/255];
            b.CData(2,:) = [96/255,125/255,139/255];
            grid on
            grid minor
            ylim([0,300])
            ylabel('Count')

            subplot(4,4,7)
            plot(xif{p-1+2},ff{p-1+2},'LineWidth',3)
            hold on
            plot(ones(1,length(yy))*perc25f{p-1+2},yy,'--r',...
                ones(1,length(yy))*perc975f{p-1+2},yy,'--r',...
                ones(1,length(yy))*contfmed{p-1+2},yy,'--r','LineWidth',1)
            xlim([xminf,xmaxf])
            ylim([0,15])
            hold off
            grid on
            grid minor
            title(sprintf('Perc 2.5 = %0.2f | Mean = %0.2f | Median = %0.2f | Perc 97.5 = %0.2f',...
                perc25f{p-1+2},contfmean{p-1+2},contfmed{p-1+2},perc975f{p-1+2}))
            
            subplot(4,4,8)
            counts = [sum(contfit(p-1+2,:) < contfmean{p-1+2}),...
                sum(contfit(p-1+2,:) > contfmean{p-1+2})];
            b = bar(c,counts,'FaceColor','flat');
            b.CData(1,:) = [97/255,97/255,97/255];
            b.CData(2,:) = [96/255,125/255,139/255];
            grid on
            grid minor
            ylim([0,300])
            ylabel('Count')
            
            subplot(4,4,9)
            plot(xia{p-1+3},fa{p-1+3},'LineWidth',3)
            hold on
            plot(ones(1,length(yy))*perc25a{p-1+3},yy,'--r',...
                ones(1,length(yy))*perc975a{p-1+3},yy,'--r',...
                ones(1,length(yy))*contamed{p-1+3},yy,'--r','LineWidth',1)
            xlim([xmina,xmaxa])
            ylim([0,0.02])
            hold off
            grid on
            grid minor
            ylabel(sprintf('%s\nDensity',splt{3}))
            title(sprintf('Perc 2.5 = %0.2f | Mean = %0.2f | Median = %0.2f | Perc 97.5 = %0.2f',...
                perc25a{p-1+3},contamean{p-1+3},contamed{p-1+3},perc975a{p-1+3}))
            
            subplot(4,4,10)
            counts = [sum(contavg(p-1+3,:) < contamean{p-1+3}),...
                sum(contavg(p-1+3,:) > contamean{p-1+3})];
            b = bar(c,counts,'FaceColor','flat');
            b.CData(1,:) = [97/255,97/255,97/255];
            b.CData(2,:) = [96/255,125/255,139/255];
            grid on
            grid minor
            ylim([0,300])
            ylabel('Count')

            subplot(4,4,11)
            plot(xif{p-1+3},ff{p-1+3},'LineWidth',3)
            hold on
            plot(ones(1,length(yy))*perc25f{p-1+3},yy,'--r',...
                ones(1,length(yy))*perc975f{p-1+3},yy,'--r',...
                ones(1,length(yy))*contfmed{p-1+3},yy,'--r','LineWidth',1)
            xlim([xminf,xmaxf])
            ylim([0,15])
            hold off
            grid on
            grid minor
            title(sprintf('Perc 2.5 = %0.2f | Mean = %0.2f | Median = %0.2f | Perc 97.5 = %0.2f',...
                perc25f{p-1+3},contfmean{p-1+3},contfmed{p-1+3},perc975f{p-1+3}))
            
            subplot(4,4,12)
            counts = [sum(contfit(p-1+3,:) < contfmean{p-1+3}),...
                sum(contfit(p-1+3,:) > contfmean{p-1+3})];
            b = bar(c,counts,'FaceColor','flat');
            b.CData(1,:) = [97/255,97/255,97/255];
            b.CData(2,:) = [96/255,125/255,139/255];
            grid on
            grid minor
            ylim([0,300])
            ylabel('Count')
            
            subplot(4,4,13)
            plot(xia{p-1+4},fa{p-1+4},'LineWidth',3)
            hold on
            plot(ones(1,length(yy))*perc25a{p-1+4},yy,'--r',...
                ones(1,length(yy))*perc975a{p-1+4},yy,'--r',...
                ones(1,length(yy))*contamed{p-1+4},yy,'--r','LineWidth',1)
            xlim([xmina,xmaxa])
            ylim([0,0.02])
            hold off
            grid on
            grid minor
            ylabel(sprintf('%s\nDensity',splt{4}))
            xlabel('Pixel Count')
            title(sprintf('Perc 2.5 = %0.2f | Mean = %0.2f | Median = %0.2f | Perc 97.5 = %0.2f',...
                perc25a{p-1+4},contamean{p-1+4},contamed{p-1+4},perc975a{p-1+4}))
            
            subplot(4,4,14)
            counts = [sum(contavg(p-1+4,:) < contamean{p-1+4}),...
                sum(contavg(p-1+4,:) > contamean{p-1+4})];
            b = bar(c,counts,'FaceColor','flat');
            b.CData(1,:) = [97/255,97/255,97/255];
            b.CData(2,:) = [96/255,125/255,139/255];
            grid on
            grid minor
            ylim([0,300])
            ylabel('Count')

            subplot(4,4,15)
            plot(xif{p-1+4},ff{p-1+4},'LineWidth',3)
            hold on
            plot(ones(1,length(yy))*perc25f{p-1+4},yy,'--r',...
                ones(1,length(yy))*perc975f{p-1+4},yy,'--r',...
                ones(1,length(yy))*contfmed{p-1+4},yy,'--r','LineWidth',1)
            xlim([xminf,xmaxf])
            ylim([0,15])
            hold off
            grid on
            grid minor
            xlabel('Fitness')
            title(sprintf('Perc 2.5 = %0.2f | Mean = %0.2f | Median = %0.2f | Perc 97.5 = %0.2f',...
                perc25f{p-1+4},contfmean{p-1+4},contfmed{p-1+4},perc975f{p-1+4}))
            
            subplot(4,4,16)
            counts = [sum(contfit(p-1+4,:) < contfmean{p-1+4}),...
                sum(contfit(p-1+4,:) > contfmean{p-1+4})];
            b = bar(c,counts,'FaceColor','flat');
            b.CData(1,:) = [97/255,97/255,97/255];
            b.CData(2,:) = [96/255,125/255,139/255];
            grid on
            grid minor
            ylim([0,300])
            ylabel('Count')
            
            sgtitle(sprintf('%s | Control Distribution | Time = %d hrs\nPlate - %d',...
                expt,hours(iii),p))
            saveas(fig,sprintf('%s%s_ContDist_%d_P%d.png',...
                        out_path,expt_name,hours(iii),p))
        end
    end
    
%%  IN THE WORKS

    c = {'< mean','> mean'};
            
    x1 = contfit(p-1+4,contfit(p-1+4,:) < contfmean{p-1+4})';
    x2 = contfit(p-1+4,contfit(p-1+4,:) > contfmean{p-1+4})';

    group = [    ones(size(x1));
             2 * ones(size(x2))];

    figure()
    boxplot([x1; x2],group,'Labels',c)
    grid on
    grid minor
    ylim([xminf,xmaxf])
    
%%  END
