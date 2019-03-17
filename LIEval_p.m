%%  Sau MATLAB Colony Analyzer Toolkit
%
%%  LIEval.m

%   Author: Saurin Parikh, February 2019
%   dr.saurin.parikh@gmail.com

%   Evaluation of the LI based control normalization

%%  Load Paths to Files and Data

    cd /home/sbp29/MATLAB

    addpath('/home/sbp29/MATLAB/Matlab-Colony-Analyzer-Toolkit')
    addpath('/home/sbp29/MATLAB/bean-matlab-toolkit')
    addpath('/home/sbp29/MATLAB/sau-matlab-toolkit')
    addpath('/home/sbp29/MATLAB/sau-matlab-toolkit/grid-manipulation')
    addpath('/home/sbp29/MATLAB/paris')
    addpath('/home/sbp29/MATLAB/development')

    javaaddpath('/home/sbp29/MATLAB/mysql-connector-java-8.0.12.jar');

%%  Add MCA toolkit to Path
%     add_mca_toolkit_to_path

%%  Initialization

%     Set preferences with setdbprefs.
    setdbprefs('DataReturnFormat', 'structure');
    setdbprefs({'NullStringRead';'NullStringWrite';'NullNumberRead';'NullNumberWrite'},...
                  {'null';'null';'NaN';'NaN'})

    expt_name = '4C3_GA';
    density = 6144;
    
%   MySQL Table Details  
    
    tablename_jpeg      = sprintf('%s_%d_JPEG',expt_name,density);
    tablename_norm      = sprintf('%s_%d_NORM',expt_name,density);
    tablename_fit       = sprintf('%s_%d_FITNESS',expt_name,density);
    tablename_fits      = sprintf('%s_%d_FITNESS_STATS',expt_name,density);
    tablename_es        = sprintf('%s_%d_FITNESS_ES',expt_name,density);
    tablename_pval      = sprintf('%s_%d_PVALUE',expt_name,density);
    tablename_res       = sprintf('%s_%d_RES',expt_name,density);
    
    tablename_p2o       = '4C3_pos2orf_name';
    tablename_bpos      = '4C3_borderpos';
    
%   Reference Strain Name

    cont.name           = 'BF_control';
    
%   MySQL Connection and fetch initial data

    connectSQL;
    
    p2c_info(1,:) = '4C3_pos2coor6144';
    p2c_info(2,:) = '6144plate       ';
    p2c_info(3,:) = '6144col         ';
    p2c_info(4,:) = '6144row         ';

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
    
    hours = fetch(conn, sprintf(['select distinct hours from %s ',...
            'order by hours asc'], tablename_jpeg));
    hours = hours.hours;
    
%%  PLATEWISE RMSE

    for i = 1:length(hours)
        for iii = 1:length(n_plates.x6144plate)
            clear rmse
            bg = fetch(conn, sprintf(['select a.* ',...
                'from %s a, %s b ',...
                'where a.hours = %d ',...
                'and a.pos = b.pos ',...
                'and b.%s = %d ',...
                'order by b.%s, b.%s'],...
                tablename_fit,p2c_info(1,:),hours(i),p2c_info(2,:),...
                iii,p2c_info(3,:),p2c_info(4,:)));

            rmse(:,iii) = (bg.average - bg.bg).^2;

            max_avg = max(bg.average);
            min_avg = min(bg.average);

            fig = figure('Renderer', 'painters', 'Position', [10 10 1920 1200],'visible','off');
%            figure()
            subplot(2,2,1)
            heatmap(col2grid(bg.average),'ColorLimits',[min_avg max_avg]);
            title(sprintf('Observed Pixel Count\n(Plate %d, %d hr)',iii,hours(i)))
            title('Observed Pixel Count')
            subplot(2,2,2)
            heatmap(col2grid(bg.bg),'ColorLimits',[min_avg max_avg]);
%             title(sprintf('Predicted Pixel Count\n(Plate %d, %d hr)',iii,hours(i)))
            title('Predicted Pixel Count')
            subplot(2,2,3)
            heatmap(col2grid(bg.fitness),'ColorLimits',[0.7 1.4]);
            title('Fitness')
            subplot(2,2,4)
            heatmap(col2grid(abs(bg.average - bg.bg)),'ColorLimits',[0 120]);
            title(sprintf('RMSE (%0.3f)',sqrt(nanmean(rmse(:,iii)))))
            colormap parula
            saveas(fig,sprintf('overview%d_%d.png',iii,hours(i)))
        end
    end

%%  POWER, FALSE POSITIVE AND ES
    
    p = 0:0.01:1;
    %hours = 12.5;
    for ii=1:length(hours)
        if conn.isopen == 0
            connectSQL;
        end
        cont_data = fetch(conn, sprintf(['select * from %s ',...
            'where orf_name = ''%s'' ',...
            'and fitness is not NULL and hours = %d'],...
            tablename_fit,cont.name,hours(ii)));

        rest_data = fetch(conn, sprintf(['select * from %s ',...
            'where orf_name != ''%s'' ',...
            'and fitness is not NULL and hours = %d'],...
            tablename_fit,cont.name,hours(ii)));

        ef_size = mean(rest_data.fitness)/mean(cont_data.fitness);

        cont_dist = [];
        cont_means = [];
        for i=1:100000
            cont_dist(i,:) = datasample(cont_data.fitness, 8, 'Replace', false);
            cont_means(i,:) = mean(cont_dist(i,:));
        end
        
        clear restmean reststd stat pvals;
        contmean = nanmean(cont_means);
        contstd = nanstd(cont_means);

        m = cont_means;
        tt = length(m);

        rest_dist =[];
        rest_means = [];
%        ss = 8;
         for ss=1:8
            for i=1:100000
                rest_dist{ss}(i,:) = datasample(rest_data.fitness, ss, 'Replace', false);
                rest_means{ss}(i,:) = mean(rest_dist{ss}(i,:));
            end
            restmean{ss} = nanmean(rest_means{ss});
            reststd{ss} = nanstd(rest_means{ss});
            
            temp_p = [];
            temp_s = [];
            for i = 1:length(rest_means{ss})
                if sum(m<rest_means{ss}(i)) < tt/2
                    if m<rest_means{ss}(i) == 0
                        temp_p = [temp_p; 1/tt];
                        temp_s = [temp_s; (rest_means{ss}(i) - contmean)/contstd];
                    else
                        temp_p = [temp_p; ((sum(m<=rest_means{ss}(i)))/tt)*2];
                        temp_s = [temp_s; (rest_means{ss}(i) - contmean)/contstd];
                    end
                else
                    temp_p = [temp_p; ((sum(m>=rest_means{ss}(i)))/tt)*2];
                    temp_s = [temp_s; (rest_means{ss}(i) - contmean)/contstd];
                end
            end
            pvals{ii}{ss} = temp_p; stat{ii}{ss} = temp_s;
            pow = (sum(temp_p<0.05)/length(rest_means{ss}))*100;

%             figure()
%            fig = figure('Renderer', 'painters', 'Position', [10 10 480 300],'visible','off');
%            [f,xi] = ksdensity(cont_means);
%            plot(xi,f,'LineWidth',3)
%            hold on
%            [f,xi] = ksdensity(rest_means{ss});
%            plot(xi,f,'LineWidth',3)
%            legend('control','rest of plate')
%            title(sprintf('ES = %0.3f',ef_size))
%            xlabel('Fitness')
%            ylabel('Density')
%            grid on
%            hold off
%            saveas(fig,sprintf('es_%d_ss%d.png',hours(ii), ss))
            
            len = length(pvals{ii}{ss});
            fpdat = [];

            for i = 1:length(p)
                fp = sum(pvals{ii}{ss} <= p(i));
                fpdat = [fpdat; [p(i), fp/len]];
            end

%             fig = figure('Renderer', 'painters', 'Position', [10 10 480 300],'visible','off');
%             histogram(pvals{ii}{ss}, 'Normalization', 'cdf')
%             hold on
%             plot(0:0.01:1,0:0.01:1,'--r','LineWidth',3)
%             grid on
%             xlabel('P Value Cut-offs')
%             ylabel('Proportion of Colonies')
%             title(sprintf('Time = %d hrs (SS = %d)',hours(ii),ss))
%             xlim([0,1])
%             ylim([0,1])
%             saveas(fig,sprintf('pval_colonies_%d_ss%d.png',hours(ii),ss))
            
         end

    fprintf('time %d hrs done\n', hours(ii))
    end
    
%%  NULL DISTRIBUTION
%   The fitness distribution of the positions used to create the LI model
    
    m = cont_data.fitness;
    tt = length(m);

    contp = [];
    for i = 1:100000
        temp = mean(datasample(cont_data.fitness, 1, 'Replace', false));
        if sum(m<temp) < tt/2
            if m<temp == 0
                contp = [contp; 1/tt];
            else
                contp = [contp; ((sum(m<=temp))/tt)*2];
            end
        else
            contp = [contp; ((sum(m>=temp))/tt)*2];
        end
    end
    contp(contp>1) = 1;
    
%    figure()
%    histogram(contp, 'NumBins', 20, 'Normalization', 'pdf')
%    grid on
%    xlabel('P Values')
%    ylabel('Probability Density')
%    title('NULL DISTRIBUTION')

%%  DATA UNDER PVAL CUT-OFFS
    
    p = 0:0.01:1;
    
    for ii = 1:length(hours)
        for ss = 1:8
            len = length(pvals{ii}{ss});
            fpdat = [];

            for i = 1:length(p)
                fp = sum(pvals{ii}{ss} <= p(i));
                fpdat = [fpdat; [p(i), fp/len]];
            end

            fig = figure('Renderer', 'painters', 'Position', [10 10 480 300],'visible','off');
            histogram(pvals{ii}{ss}, 'Normalization', 'cdf')
            hold on
            plot(0:0.01:1,0:0.01:1,'--r','LineWidth',3)
            grid on
            xlabel('P Value Cut-offs')
            ylabel('Proportion of Colonies')
            title(sprintf('Time = %d hrs',hours(ii)))
            xlim([0,1])
            ylim([0,1])
            saveas(fig,sprintf('pval_colonies_%d.png',hours(ii)))
        end
    end
    

%%  MISSING DATA WHEN CALCULATING BG
%   Refer gapEVAL_p.m for this.

%%  VIRTUAL PLATE POWER ANALYSIS
 
    for ss=[1:8]
        fprintf('sample size = %d\n',ss)
        cont_hrs = 12.5;
        rest_hrs = [8,9,10,11,14,15.5,17,18];
        data = [];
    %     ss = 1;

        for ii = 1:length(rest_hrs)
    %         if conn.isopen == 0
    %             connectSQL;
    %         end
            plate_fit = [];
            cont_fit = [];
            rest_fit = [];
            for iii = 1:length(n_plates.x6144plate)
                pos.all = fetch(conn, sprintf(['select a.pos ',...
                    'from %s a ',...
                    'where %s = %d ',...
                    'order by %s, %s'],...
                    p2c_info(1,:),...
                    p2c_info(2,:),...
                    n_plates.x6144plate(iii),...
                    p2c_info(3,:),...
                    p2c_info(4,:)));

                pos.cont = fetch(conn, sprintf(['select a.pos ',...
                    'from %s a, %s b ',...
                    'where a.pos = b.pos and %s = %d and a.orf_name = ''%s'' ',...
                    'order by %s, %s'],...
                    tablename_p2o,...
                    p2c_info(1,:),...
                    p2c_info(2,:),...
                    n_plates.x6144plate(iii),...
                    cont.name,...
                    p2c_info(3,:),...
                    p2c_info(4,:)));

                cont_pos = col2grid(ismember(pos.all.pos, pos.cont.pos));
                rest_pos = ~cont_pos;

                cont_data = fetch(conn, sprintf(['select a.* ',...
                    'from %s a, %s b ',...
                    'where a.hours = %d ',...
                    'and a.pos = b.pos and b.%s = %d ',...
                    'order by b.%s, b.%s'],...
                    tablename_fit,p2c_info(1,:),cont_hrs,...
                    p2c_info(2,:),1,p2c_info(3,:),p2c_info(4,:)));

                cont_avg = col2grid(cont_data.average).*cont_pos;
                plate_bg = col2grid(cont_data.bg);

                rest_data = fetch(conn, sprintf(['select a.* ',...
                    'from %s a, %s b ',...
                    'where a.hours = %d ',...
                    'and a.pos = b.pos and b.%s = %d ',...
                    'order by b.%s, b.%s'],...
                    tablename_fit,p2c_info(1,:),rest_hrs(ii),...
                    p2c_info(2,:),1,p2c_info(3,:),p2c_info(4,:)));

                rest_avg = col2grid(rest_data.average).*rest_pos;
                plate_avg = cont_avg + rest_avg;

                plate_fit = plate_avg./plate_bg;

                cont_fit_tmp = plate_fit.*cont_pos;
                cont_fit_tmp = cont_fit_tmp(cont_fit_tmp > 0);
                cont_fit = [cont_fit; cont_fit_tmp];

                rest_fit_tmp = plate_fit.*rest_pos;
                rest_fit_tmp = rest_fit_tmp(rest_fit_tmp>0); 
                rest_fit = [rest_fit; rest_fit_tmp];
            end

        % % % % % % % % % 

            cont_dist = [];
            cont_means = [];
            for i=1:100000
                cont_dist(i,:) = datasample(cont_fit, 8, 'Replace', false);
                cont_means(i,:) = mean(cont_dist(i,:));
            end

            rest_dist =[];
            rest_means = [];
            for i=1:100000
                rest_dist(i,:) = datasample(rest_fit, ss, 'Replace', false);
                rest_means(i,:) = mean(rest_dist(i,:));
            end

            contmean = nanmean(cont_means);
            contstd = nanstd(cont_means);
            restmean = nanmean(rest_means);
            reststd = nanstd(rest_means);

            m = cont_means;
            tt = length(m);

            temp_p = [];
            temp_s = [];
            for i = 1:length(rest_means)
                if sum(m<rest_means(i)) < tt/2
                    if m<rest_means(i) == 0
                        temp_p = [temp_p; 1/tt];
                        temp_s = [temp_s; (rest_means(i) - contmean)/contstd];
                    else
                        temp_p = [temp_p; ((sum(m<=rest_means(i)))/tt)*2];
                        temp_s = [temp_s; (rest_means(i) - contmean)/contstd];
                    end
                else
                    temp_p = [temp_p; ((sum(m>=rest_means(i)))/tt)*2];
                    temp_s = [temp_s; (rest_means(i) - contmean)/contstd];
                end
            end

            ef_size = mean(rest_fit)/mean(cont_fit);
            pow = (sum(temp_p<0.05)/length(rest_means))*100;
            avg_diff = abs(nanmean(nanmean(cont_avg)) - nanmean(nanmean(rest_avg)));

            data = [data; ef_size, pow, avg_diff];

%             figure()
            fig = figure('Renderer', 'painters', 'Position', [10 10 480 300],'visible','off');
            [f,xi] = ksdensity(cont_means);
            plot(xi,f,'LineWidth',3)
            xlim([0.75,1.25])
            ylim([0,30])
            hold on
            [f,xi] = ksdensity(rest_means);
            plot(xi,f,'LineWidth',3)
            legend('control','rest of plate')
            title(sprintf(['TimeC = %d hrs | TimeR = %d hrs \n ',...
                'ES = %0.3f'],...
                cont_hrs,rest_hrs(ii),ef_size))
            xlabel('Fitness')
            ylabel('Density')
            grid on
            grid minor
            hold off
            saveas(fig,sprintf('vp_powes_%d_%d.png',cont_hrs,rest_hrs(ii)))
        fprintf('Virtual plate %d hrs V/S %d hrs done!\n', cont_hrs,rest_hrs(ii))
        end

    %%  POWER vs ES

        [~, i] = sort(data(:,1));
        es_pow = data(i, :);

        x{cont_hrs}{ss}   = es_pow(:,1);
        y{cont_hrs}{ss}   = es_pow(:,2);
        xx{cont_hrs}{ss}  = min(es_pow(:,1)):.001:max(es_pow(:,1));
        yy{cont_hrs}{ss}  = interp1(x{cont_hrs}{ss},y{cont_hrs}{ss},...
            xx{cont_hrs}{ss},'pchip');

    %     figure()
        fig = figure('Renderer', 'painters', 'Position', [10 10 960 800],'visible','off');
        plot(xx{cont_hrs}{ss},yy{cont_hrs}{ss},'r','LineWidth',2)
        hold on
        scatter(x{cont_hrs}{ss}, y{cont_hrs}{ss},'MarkerEdgeColor',[0 .5 .5],...
                  'MarkerFaceColor',[0 .7 .7],...
                  'LineWidth',2);
        grid on
        grid minor
        xlim([0,3])
        ylim([0,101])
        xlabel('Effect Size')
        ylabel('Power')
        title(sprintf('Power V/S ES (SS = %d)',ss))
        hold off
        saveas(fig,sprintf('vp_powes_%d_%d.png',cont_hrs,ss))
    end
    
%%  EFFECT OF SCHEMES AND INTERLEAVING

%    expt1 = "4C2_R1";
%    expt2 = "4C2_R2";
%    expt3 = "4C2_R1_NIL";
%    expt4 = "4C2_R2_NIL";
%    
%%     hours = 24;
%    better12 = [];
%    better1n = [];
%    better2n = [];
%    
%    for i = 1:length(hours)
%        for iii = 1:length(n_plates.x6144plate)
%            bg1 = fetch(conn, sprintf(['select a.* ',...
%                'from %s_6144_FITNESS a, %s b ',...
%                'where a.hours = %d ',...
%                'and a.pos = b.pos ',...
%                'and b.%s = %d ',...
%                'order by b.%s, b.%s'],...
%                expt1,p2c_info(1,:),hours(i),p2c_info(2,:),...
%                iii,p2c_info(3,:),p2c_info(4,:)));
%
%            bg2 = fetch(conn, sprintf(['select a.* ',...
%                'from %s_6144_FITNESS a, %s b ',...
%                'where a.hours = %d ',...
%                'and a.pos = b.pos ',...
%                'and b.%s = %d ',...
%                'order by b.%s, b.%s'],...
%                expt2,p2c_info(1,:),hours(i),p2c_info(2,:),...
%                iii,p2c_info(3,:),p2c_info(4,:)));
%
%            bg3 = fetch(conn, sprintf(['select a.* ',...
%                'from %s_6144_FITNESS a, %s b ',...
%                'where a.hours = %d ',...
%                'and a.pos = b.pos ',...
%                'and b.%s = %d ',...
%                'order by b.%s, b.%s'],...
%                expt3,p2c_info(1,:),hours(i),p2c_info(2,:),...
%                iii,p2c_info(3,:),p2c_info(4,:)));
%
%            bg4 = fetch(conn, sprintf(['select a.* ',...
%                'from %s_6144_FITNESS a, %s b ',...
%                'where a.hours = %d ',...
%                'and a.pos = b.pos ',...
%                'and b.%s = %d ',...
%                'order by b.%s, b.%s'],...
%                expt4,p2c_info(1,:),hours(i),p2c_info(2,:),...
%                iii,p2c_info(3,:),p2c_info(4,:)));
%
%            for ii = 1:length(bg1.average)
%                e1((iii-1)*6144+ii,1) = (bg1.average(ii) - bg1.bg(ii)).^2;
%                e2((iii-1)*6144+ii,1) = (bg2.average(ii) - bg2.bg(ii)).^2;
%                e3((iii-1)*6144+ii,1) = (bg3.average(ii) - bg3.bg(ii)).^2;
%                e4((iii-1)*6144+ii,1) = (bg4.average(ii) - bg4.bg(ii)).^2;
%            end
%        end
%
%        e1 = e1(~isnan(e1));
%        e2 = e2(~isnan(e2));
%        e3 = e3(~isnan(e3));
%        e4 = e4(~isnan(e4));
%
%        better12 = [better12;[hours(i),sqrt(mean(e1.^2)),sqrt(mean(e2.^2)),...
%            ranksum(e1,e2,'tail','right'),...
%            ranksum(e1,e2,'tail','left')]];
%        better1n = [better1n;[hours(i),sqrt(mean(e1.^2)),sqrt(mean(e3.^2)),...
%            ranksum(e1,e3,'tail','right'),...
%            ranksum(e1,e3,'tail','left')]];
%        better2n = [better2n;[hours(i),sqrt(mean(e2.^2)),sqrt(mean(e4.^2)),...
%            ranksum(e2,e4,'tail','right'),...
%            ranksum(e2,e4,'tail','left')]];
%    end
%
%    
%    sum(better12(:,4)<0.05)/length(better12)
%    sum(better1n(:,4)<0.05)/length(better1n)
%    sum(better2n(:,4)<0.05)/length(better2n)
    
    
%     exec(conn, 'drop table 4C2_R2R2N_RMSE');
%     exec(conn, ['create table 4C2_R2R2N_RMSE( ',...
%                 'hours int(11) not NULL, ',...
%                 'R2_rmse double default NULL, ',...
%                 'R2N_rmse double default NULL, ',...
%                 'p_right double default NULL, ',...
%                 'p_left double default NULL)']);
%             
%     datainsert(conn, '4C2_R2R2N_RMSE',...
%         {'hours','R2_rmse','R2N_rmse','p_right','p_left'},better2n);

