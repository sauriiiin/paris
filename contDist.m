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
    out_path = '/Users/saur1n/Desktop/4C3/Analysis/GA/S1Analysis/contdist/rmoutliers/';
    density = 6144;

%   MySQL Table Details  
    tablename_fit       = sprintf('%s_%d_FITNESS',expt_name,density);
    tablename_p2o       = '4C3_pos2orf_name1';
    tablename_bpos      = '4C3_borderpos';
    
    p2c_info(1,:) = '4C3_pos2coor6144';
    p2c_info(2,:) = '6144plate       ';
    p2c_info(3,:) = '6144col         ';
    p2c_info(4,:) = '6144row         ';
    
%   Reference Strain Name
    cont.name           = 'BF_control';

%   MySQL Connection and fetch initial data
    connectSQL;
    
    hours = fetch(conn, sprintf(['select distinct hours from %s ',...
                 'order by hours asc'], tablename_fit));
    hours = hours.hours;
    
    p2c = fetch(conn, sprintf(['select * from %s a ',...
        'order by a.%s, a.%s, a.%s'],...
        p2c_info(1,:),...
        p2c_info(2,:),...
        p2c_info(3,:),...
        p2c_info(4,:)));
    
    n_plates = fetch(conn, sprintf(['select distinct %s from %s a ',...
        'order by %s asc'],...
        p2c_info(2,:),...
        p2c_info(1,:),...
        p2c_info(2,:)));
    
    splt = {'TL','TR','BL','BR'};
    
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
        contbg = [];
        for ii = 1:size(contpos,1)
            temp = fetch(conn, sprintf(['select * from %s ',...
                'where hours = %d and pos in (%s) ',...
                'and fitness is not null'],tablename_fit,hours(iii),...
                sprintf('%d,%d,%d,%d,%d,%d,%d,%d',contpos(ii,:))));
            if nansum(temp.average) > 0
                contavg = [contavg, temp.average];
                contbg = [contbg, temp.bg];
%                 contfit = [contfit, temp.fitness];
                [~,outlier] = rmoutliers(temp.fitness);
                temp.fitness(outlier) = NaN;
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
            conterr{i} = contavg(i,:) - contbg(i,:);   conterr{i} = conterr{i}(~isnan(conterr{i}));
            u{i} = conterr{i}(conterr{i} > 0)';         d{i} = conterr{i}(conterr{i} < 0)';
            xu{i} = [1:length(u{i}), fliplr(1:length(u{i}))];
            inBetweenU{i} = [sort(u{i})', fliplr(zeros(1,length(u{i})))];
            xd{i} = [1:length(d{i}), fliplr(1:length(d{i}))];
            inBetweenD{i} = [sort(d{i})', fliplr(zeros(1,length(d{i})))];
        end
        
        xmina = round(min(min(contavg))-50,-1);
        xmaxa = round(max(max(contavg))+50,-1);
        xminf = round(min(min(contfit)),3);
        xmaxf = round(max(max(contfit)),3);

        fig = figure('Renderer','painters',...
                'Position',[10 10 1300 500],...
                'visible','off');
%         figure()
        subplot(1,2,1)
        hello = nanmean(contfit,1);
        [f,xi] = ksdensity(hello);
        plot(xi,f,'LineWidth',3)
        hold on
        plot(ones(1,length(yy))*prctile(hello,2.5),yy,'--r',...
            ones(1,length(yy))*prctile(hello,97.5),yy,'--r',...
            ones(1,length(yy))*nanmedian(hello),yy,'--r','LineWidth',1)
        ylabel('Density')
        xlabel('Fitness')
        xlim([0.8,1.2])
        ylim([0,25])
        hold off
        grid on
        grid minor
        title(sprintf('Perc 2.5 = %0.3f | Mean = %0.3f | Median = %0.3f | Perc 97.5 = %0.3f',...
            prctile(hello,2.5),nanmean(hello),nanmedian(hello),prctile(hello,97.5)))
        
        subplot(1,2,2)
        hello = rmoutliers(hello);
        [f,xi] = ksdensity(hello);
        plot(xi,f,'LineWidth',3)
        hold on
        plot(ones(1,length(yy))*prctile(hello,2.5),yy,'--r',...
            ones(1,length(yy))*prctile(hello,97.5),yy,'--r',...
            ones(1,length(yy))*nanmedian(hello),yy,'--r','LineWidth',1)
        ylabel('Density')
        xlabel('Fitness')
        xlim([0.8,1.2])
        ylim([0,25])
        hold off
        grid on
        grid minor
        title(sprintf('Perc 2.5 = %0.3f | Mean = %0.3f | Median = %0.3f | Perc 97.5 = %0.3f',...
            prctile(hello,2.5),nanmean(hello),nanmedian(hello),prctile(hello,97.5)))
        
        sgtitle(sprintf('%s | Removing Outliers | Time = %d hrs',...
                expt,hours(iii)))
        saveas(fig,sprintf('%s%s_OUTLIERS_%d.png',...
                    out_path,expt_name,hours(iii)))
                
        
%         for p = 1:2
%             fig = figure('Renderer','painters',...
%                 'Position',[10 10 2000 2000],...
%                 'visible','off');
% %             figure()
%             subplot(4,3,1)
%             plot(xia{p-1+1},fa{p-1+1},'LineWidth',3)
%             hold on
%             plot(ones(1,length(yy))*perc25a{p-1+1},yy,'--r',...
%                 ones(1,length(yy))*perc975a{p-1+1},yy,'--r',...
%                 ones(1,length(yy))*contamed{p-1+1},yy,'--r','LineWidth',1)
%             xlim([xmina,xmaxa])
%             ylim([0,0.02])
%             hold off
%             grid on
%             grid minor
%             ylabel(sprintf('%s\nDensity',splt{1}))
%             title(sprintf('Perc 2.5 = %0.2f | Mean = %0.2f | Median = %0.2f | Perc 97.5 = %0.2f',...
%                 perc25a{p-1+1},contamean{p-1+1},contamed{p-1+1},perc975a{p-1+1}))
% 
%             subplot(4,3,2)
%             plot(xif{p-1+1},ff{p-1+1},'LineWidth',3)
%             hold on
%             plot(ones(1,length(yy))*perc25f{p-1+1},yy,'--r',...
%                 ones(1,length(yy))*perc975f{p-1+1},yy,'--r',...
%                 ones(1,length(yy))*contfmed{p-1+1},yy,'--r','LineWidth',1)
%             xlim([xminf,xmaxf])
%             ylim([0,15])
%             hold off
%             grid on
%             grid minor
%             title(sprintf('Perc 2.5 = %0.2f | Mean = %0.2f | Median = %0.2f | Perc 97.5 = %0.2f',...
%                 perc25f{p-1+1},contfmean{p-1+1},contfmed{p-1+1},perc975f{p-1+1}))
%             
%             subplot(4,3,3)
%             plot(1:max(length(u{p-1+1}),length(d{p-1+1})),...
%                 zeros(1,max(length(u{p-1+1}),length(d{p-1+1}))),...
%                 '--k','LineWidth',2)
%             hold on
%             plot(sort(u{p-1+1}),'.b')
%             plot(sort(d{p-1+1}),'.r')
%             text(max(length(u{p-1+1}),length(d{p-1+1}))/2 - 15,40,...
%                 ['\Sigma = ',sprintf('%.2f',sum(u{p-1+1}))],...
%                 'Color','blue','FontSize',15)
%             text(max(length(u{p-1+1}),length(d{p-1+1}))/2 - 15,-40,...
%                 ['\Sigma = ',sprintf('%.2f',sum(d{p-1+1}))],...
%                 'Color','red','FontSize',15)
%             grid on
%             grid minor
%             xlim([0,max(length(u{p-1+1}),length(d{p-1+1}))])
%             ylim([-80,80])
%             xticks([])
%             ylabel('Standard Error (Pixels)')
%             title('Pos-Wise Standard Error')
%             hold on
%             fill(xu{p-1+1},inBetweenU{p-1+1}, 'b')
%             alpha(0.25)
%             fill(xd{p-1+1}, inBetweenD{p-1+1}, 'r')
%             alpha(0.25)
%             hold off
%             
%             subplot(4,3,4)
%             plot(xia{p-1+2},fa{p-1+2},'LineWidth',3)
%             hold on
%             plot(ones(1,length(yy))*perc25a{p-1+2},yy,'--r',...
%                 ones(1,length(yy))*perc975a{p-1+2},yy,'--r',...
%                 ones(1,length(yy))*contamed{p-1+2},yy,'--r','LineWidth',1)
%             xlim([xmina,xmaxa])
%             ylim([0,0.02])
%             hold off
%             grid on
%             grid minor
%             ylabel(sprintf('%s\nDensity',splt{2}))
%             title(sprintf('Perc 2.5 = %0.2f | Mean = %0.2f | Median = %0.2f | Perc 97.5 = %0.2f',...
%                 perc25a{p-1+2},contamean{p-1+2},contamed{p-1+2},perc975a{p-1+2}))
% 
%             subplot(4,3,5)
%             plot(xif{p-1+2},ff{p-1+2},'LineWidth',3)
%             hold on
%             plot(ones(1,length(yy))*perc25f{p-1+2},yy,'--r',...
%                 ones(1,length(yy))*perc975f{p-1+2},yy,'--r',...
%                 ones(1,length(yy))*contfmed{p-1+2},yy,'--r','LineWidth',1)
%             xlim([xminf,xmaxf])
%             ylim([0,15])
%             hold off
%             grid on
%             grid minor
%             title(sprintf('Perc 2.5 = %0.2f | Mean = %0.2f | Median = %0.2f | Perc 97.5 = %0.2f',...
%                 perc25f{p-1+2},contfmean{p-1+2},contfmed{p-1+2},perc975f{p-1+2}))
%             
%             subplot(4,3,6)
%             plot(1:max(length(u{p-1+2}),length(d{p-1+2})),...
%                 zeros(1,max(length(u{p-1+2}),length(d{p-1+2}))),...
%                 '--k','LineWidth',2)
%             hold on
%             plot(sort(u{p-1+2}),'.b')
%             plot(sort(d{p-1+2}),'.r')
%             text(max(length(u{p-1+2}),length(d{p-1+2}))/2 - 15,40,...
%                 ['\Sigma = ',sprintf('%.2f',sum(u{p-1+2}))],...
%                 'Color','blue','FontSize',15)
%             text(max(length(u{p-1+2}),length(d{p-1+2}))/2 - 15,-40,...
%                 ['\Sigma = ',sprintf('%.2f',sum(d{p-1+2}))],...
%                 'Color','red','FontSize',15)
%             grid on
%             grid minor
%             xlim([0,max(length(u{p-1+2}),length(d{p-1+2}))])
%             ylim([-80,80])
%             xticks([])
%             ylabel('Standard Error (Pixels)')
%             title('Pos-Wise Standard Error')
%             hold on
%             fill(xu{p-1+2},inBetweenU{p-1+2}, 'b')
%             alpha(0.25)
%             fill(xd{p-1+2}, inBetweenD{p-1+2}, 'r')
%             alpha(0.25)
%             hold off
%             
%             subplot(4,3,7)
%             plot(xia{p-1+3},fa{p-1+3},'LineWidth',3)
%             hold on
%             plot(ones(1,length(yy))*perc25a{p-1+3},yy,'--r',...
%                 ones(1,length(yy))*perc975a{p-1+3},yy,'--r',...
%                 ones(1,length(yy))*contamed{p-1+3},yy,'--r','LineWidth',1)
%             xlim([xmina,xmaxa])
%             ylim([0,0.02])
%             hold off
%             grid on
%             grid minor
%             ylabel(sprintf('%s\nDensity',splt{3}))
%             title(sprintf('Perc 2.5 = %0.2f | Mean = %0.2f | Median = %0.2f | Perc 97.5 = %0.2f',...
%                 perc25a{p-1+3},contamean{p-1+3},contamed{p-1+3},perc975a{p-1+3}))
% 
%             subplot(4,3,8)
%             plot(xif{p-1+3},ff{p-1+3},'LineWidth',3)
%             hold on
%             plot(ones(1,length(yy))*perc25f{p-1+3},yy,'--r',...
%                 ones(1,length(yy))*perc975f{p-1+3},yy,'--r',...
%                 ones(1,length(yy))*contfmed{p-1+3},yy,'--r','LineWidth',1)
%             xlim([xminf,xmaxf])
%             ylim([0,15])
%             hold off
%             grid on
%             grid minor
%             title(sprintf('Perc 2.5 = %0.2f | Mean = %0.2f | Median = %0.2f | Perc 97.5 = %0.2f',...
%                 perc25f{p-1+3},contfmean{p-1+3},contfmed{p-1+3},perc975f{p-1+3}))
%             
%             subplot(4,3,9)
%             plot(1:max(length(u{p-1+3}),length(d{p-1+3})),...
%                 zeros(1,max(length(u{p-1+3}),length(d{p-1+3}))),...
%                 '--k','LineWidth',2)
%             hold on
%             plot(sort(u{p-1+3}),'.b')
%             plot(sort(d{p-1+3}),'.r')
%             text(max(length(u{p-1+3}),length(d{p-1+3}))/2 - 15,40,...
%                 ['\Sigma = ',sprintf('%.2f',sum(u{p-1+3}))],...
%                 'Color','blue','FontSize',15)
%             text(max(length(u{p-1+3}),length(d{p-1+3}))/2 - 15,-40,...
%                 ['\Sigma = ',sprintf('%.2f',sum(d{p-1+3}))],...
%                 'Color','red','FontSize',15)
%             grid on
%             grid minor
%             xlim([0,max(length(u{p-1+3}),length(d{p-1+3}))])
%             ylim([-80,80])
%             xticks([])
%             ylabel('Standard Error (Pixels)')
%             title('Pos-Wise Standard Error')
%             hold on
%             fill(xu{p-1+3},inBetweenU{p-1+3}, 'b')
%             alpha(0.25)
%             fill(xd{p-1+3}, inBetweenD{p-1+3}, 'r')
%             alpha(0.25)
%             hold off
%             
%             subplot(4,3,10)
%             plot(xia{p-1+4},fa{p-1+4},'LineWidth',3)
%             hold on
%             plot(ones(1,length(yy))*perc25a{p-1+4},yy,'--r',...
%                 ones(1,length(yy))*perc975a{p-1+4},yy,'--r',...
%                 ones(1,length(yy))*contamed{p-1+4},yy,'--r','LineWidth',1)
%             xlim([xmina,xmaxa])
%             ylim([0,0.02])
%             hold off
%             grid on
%             grid minor
%             ylabel(sprintf('%s\nDensity',splt{4}))
%             xlabel('Pixel Count')
%             title(sprintf('Perc 2.5 = %0.2f | Mean = %0.2f | Median = %0.2f | Perc 97.5 = %0.2f',...
%                 perc25a{p-1+4},contamean{p-1+4},contamed{p-1+4},perc975a{p-1+4}))
%             
%             subplot(4,3,11)
%             plot(xif{p-1+4},ff{p-1+4},'LineWidth',3)
%             hold on
%             plot(ones(1,length(yy))*perc25f{p-1+4},yy,'--r',...
%                 ones(1,length(yy))*perc975f{p-1+4},yy,'--r',...
%                 ones(1,length(yy))*contfmed{p-1+4},yy,'--r','LineWidth',1)
%             xlim([xminf,xmaxf])
%             ylim([0,15])
%             hold off
%             grid on
%             grid minor
%             xlabel('Fitness')
%             title(sprintf('Perc 2.5 = %0.2f | Mean = %0.2f | Median = %0.2f | Perc 97.5 = %0.2f',...
%                 perc25f{p-1+4},contfmean{p-1+4},contfmed{p-1+4},perc975f{p-1+4}))
%             
%             subplot(4,3,12)
%             plot(1:max(length(u{p-1+4}),length(d{p-1+4})),...
%                 zeros(1,max(length(u{p-1+4}),length(d{p-1+4}))),...
%                 '--k','LineWidth',2)
%             hold on
%             plot(sort(u{p-1+4}),'.b')
%             plot(sort(d{p-1+4}),'.r')
%             text(max(length(u{p-1+4}),length(d{p-1+4}))/2 - 15,40,...
%                 ['\Sigma = ',sprintf('%.2f',sum(u{p-1+4}))],...
%                 'Color','blue','FontSize',15)
%             text(max(length(u{p-1+4}),length(d{p-1+4}))/2 - 15,-40,...
%                 ['\Sigma = ',sprintf('%.2f',sum(d{p-1+4}))],...
%                 'Color','red','FontSize',15)
%             grid on
%             grid minor
%             xlim([0,max(length(u{p-1+4}),length(d{p-1+4}))])
%             ylim([-80,80])
%             xticks([])
%             ylabel('Standard Error (Pixels)')
%             title('Pos-Wise Standard Error')
%             hold on
%             fill(xu{p-1+4},inBetweenU{p-1+4}, 'b')
%             alpha(0.25)
%             fill(xd{p-1+4}, inBetweenD{p-1+4}, 'r')
%             alpha(0.25)
%             hold off
% 
%             sgtitle(sprintf('%s | Control Distribution | Time = %d hrs\nPlate - %d',...
%                 expt,hours(iii),p))
%             saveas(fig,sprintf('%s%s_ContDist_%d_P%d.png',...
%                         out_path,expt_name,hours(iii),p))
%         end
    end
    
%%  POSITIONS OF LARGE CONTROLS
%   check contdist.R on paris
 
%%  END
