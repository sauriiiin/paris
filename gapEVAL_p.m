

    try
        
    %%  INITIALIZATION

        cd /home/sbp29/MATLAB

        addpath('/home/sbp29/MATLAB/Matlab-Colony-Analyzer-Toolkit')
        addpath('/home/sbp29/MATLAB/bean-matlab-toolkit')
        addpath('/home/sbp29/MATLAB/sau-matlab-toolkit')
        addpath('/home/sbp29/MATLAB/sau-matlab-toolkit/grid-manipulation')
        addpath('/home/sbp29/MATLAB/paris')
        addpath('/home/sbp29/MATLAB/development')

        javaaddpath('/home/sbp29/MATLAB/mysql-connector-java-8.0.12.jar');

    %%  DATA GATHER

        connectSQL;

        ss = 0;
        IL = 1; % 0,1 : no, yes
        up = 2; % 0,1,2 : random gaps, upscale pattern gaps, random 96 source gaps

        expt_name = '4C2_GAP';
        density = 6144;

    %   MySQL Table Details  

        tablename_jpeg      = sprintf('%s_%d_JPEG',expt_name,density);
        tablename_norm      = sprintf('%s_%d_NORM',expt_name,density);
        tablename_fit       = sprintf('%s_%d_FITNESS',expt_name,density);
        tablename_p2o       = 'VP_pos2orf_name1';
        tablename_bpos      = 'VP_borderpos';

    %   Reference Strain Name

        cont.name           = 'BF_control';

    %   MySQL Connection and fetch initial data

        p2c_info(1,:) = 'VP_pos2coor6144';
        p2c_info(2,:) = '6144plate      ';
        p2c_info(3,:) = '6144col        ';
        p2c_info(4,:) = '6144row        ';

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

        pos.border = fetch(conn, sprintf('select pos from %s',...
             tablename_bpos));

        for iii = 1:length(ss)

            fprintf('Starting analysis for %d missing references.\n',ss(iii))

        %%  TEMPORARY FITNESS AND NORM TABLES

            [data_fit, pos_miss] = LIGap_p(hours,n_plates,p2c_info,cont.name,...
                    tablename_p2o,tablename_jpeg,pos.border,ss(iii),IL,up);

            exec(conn, sprintf('drop table %s',tablename_norm));
            exec(conn, sprintf(['create table %s ( ',...
                        'pos int(11) not NULL, ',...
                        'hours int(11) not NULL, ',...
                        'bg double default NULL, ',...
                        'average double default NULL, ',...
                        'fitness double default NULL ',...
                        ')'],tablename_norm));
            for i=1:length(hours)
                datainsert(conn, tablename_norm,...
                    {'pos','hours','bg','average','fitness'},data_fit{i});
            end

            exec(conn, sprintf('drop table %s',tablename_fit)); 
            exec(conn, sprintf(['create table %s ',...
                '(select b.orf_name, a.pos, a.hours, a.bg, a.average, a.fitness ',...
                'from %s a, %s b ',...
                'where a.pos = b.pos ',...
                'order by a.pos asc)'],tablename_fit,tablename_norm,tablename_p2o));

            if up == 0
                pos_miss = sprintf('%d,',[pos_miss;pos_miss + 100000]);
            else
                pos_miss = sprintf('%d,',pos_miss);
            end

            exec(conn, sprintf(['update %s ',...
                'set fitness = NULL ',...
                'where pos in ',...
                '(%s)'],tablename_fit,pos_miss(1:end-1)));

        %%  FALSE POSITIVE RATE

            p = 0:0.01:1;

            for ii = 1:length(hours)
                cont_data = fetch(conn, sprintf(['select * from %s ',...
                    'where orf_name = ''%s'' ',...
                    'and fitness is not NULL and hours = %d'],...
                    tablename_fit,cont.name,hours(ii)));

                rest_data = fetch(conn, sprintf(['select * from %s ',...
                    'where orf_name != ''%s'' ',...
                    'and fitness is not NULL and hours = %d'],...
                    tablename_fit,cont.name,hours(ii)));

                cont_dist = [];
                cont_means = [];
                for i=1:100000
                    cont_dist(i,:) = datasample(cont_data.fitness, 8, 'Replace', false);
                    cont_means(i,:) = mean(cont_dist(i,:));
                end

                contmean = nanmean(cont_means);
                contstd = nanstd(cont_means);

                m = cont_means;
                tt = length(m);

                rest_dist =[];
                rest_means = [];

                for i=1:100000
                    rest_dist(i,:) = datasample(rest_data.fitness, 8, 'Replace', false);
                    rest_means(i,:) = mean(rest_dist(i,:));
                end
                restmean = nanmean(rest_means);
                reststd = nanstd(rest_means);

                pvals = [];
                stat = [];
                for i = 1:length(rest_means)
                    if sum(m<rest_means(i)) < tt/2
                        if m<rest_means(i) == 0
                            pvals = [pvals; 1/tt];
                            stat = [stat; (rest_means(i) - contmean)/contstd];
                        else
                            pvals = [pvals; ((sum(m<=rest_means(i)))/tt)*2];
                            stat = [stat; (rest_means(i) - contmean)/contstd];
                        end
                    else
                        pvals = [pvals; ((sum(m>=rest_means(i)))/tt)*2];
                        stat = [stat; (rest_means(i) - contmean)/contstd];
                    end
                end

                len = length(pvals);
                fpdat = [];

                for i = 1:length(p)
                    fp = sum(pvals <= p(i));
                    fpdat = [fpdat; [p(i), fp/len]];
                end

                fig = figure('Renderer', 'painters', 'Position', [10 10 960 600],'visible','off');
                histogram(pvals, 'Normalization', 'cdf')
                hold on
                plot(0:0.01:1,0:0.01:1,'--r','LineWidth',3)
                grid on
                xlabel('P Value Cut-offs')
                ylabel('Proportion of Colonies')
                title(sprintf('Time = %d hrs | Missing %d References',hours(ii),ss(iii)))
                xlim([0,1])
                ylim([0,1])
                saveas(fig,sprintf('fpr_%d_%d.png',ss(iii),hours(ii)))
            end

            fprintf('FPR analysis for %d missing references complete.\n',ss(iii))

        %%  POWER ANALYSIS        

            cont_hrs = 19;
            rest_hrs = [14:18,20:21,24:31];
            fpr = fpr4c(tablename_fit, cont.name, cont_hrs, 8);

            data = [];

            for ii = 1:length(rest_hrs)
                plate_fit = [];
                cont_fit = [];
                rest_fit = [];
                for pp = 1:length(n_plates.x6144plate)
                    pos.all = fetch(conn, sprintf(['select a.pos ',...
                        'from %s a ',...
                        'where %s = %d ',...
                        'order by %s, %s'],...
                        p2c_info(1,:),...
                        p2c_info(2,:),...
                        n_plates.x6144plate(pp),...
                        p2c_info(3,:),...
                        p2c_info(4,:)));

                    pos.cont = fetch(conn, sprintf(['select a.pos ',...
                        'from %s a, %s b ',...
                        'where a.pos = b.pos and %s = %d and a.orf_name = ''%s'' ',...
                        'order by %s, %s'],...
                        tablename_p2o,...
                        p2c_info(1,:),...
                        p2c_info(2,:),...
                        n_plates.x6144plate(pp),...
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
                    rest_dist(i,:) = datasample(rest_fit, 8, 'Replace', false);
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

                all = length(rest_means);
                pp = sum(temp_p<0.05);
                pow = (pp/all)*100;
                fdr = 0;
                avg_diff = abs(nanmean(nanmean(cont_avg)) - nanmean(nanmean(rest_avg)));

                data = [data; ef_size, pow, fdr, avg_diff];

        %         figure()
        %         fig = figure('Renderer', 'painters', 'Position', [10 10 480 300],'visible','off');
        %         [f,xi] = ksdensity(cont_means);
        %         plot(xi,f,'LineWidth',3)
        %         xlim([0.75,1.25])
        %         ylim([0,30])
        %         hold on
        %         [f,xi] = ksdensity(rest_means);
        %         plot(xi,f,'LineWidth',3)
        %         legend('control','rest of plate')
        %         title(sprintf(['TimeC = %d hrs | TimeR = %d hrs \n ',...
        %             'ES = %0.3f | Power = %0.3f'],...
        %             cont_hrs,rest_hrs(ii),ef_size,pow))
        %         xlabel('Fitness')
        %         ylabel('Density')
        %         grid on
        %         grid minor
        %         hold off
        %         saveas(fig,sprintf('vp_powes_%d_%d.png',cont_hrs,rest_hrs(ii)))
            fprintf('Virtual plate %d hrs V/S %d hrs for %d missing references done.\n',...
                cont_hrs,rest_hrs(ii),ss(iii))
            fprintf('%d more comparisons to go!\n',length(rest_hrs)-ii)
            end

        %%  POWER vs ES

            [~, i] = sort(data(:,1));
            es_pow = data(i, :);

            x   = es_pow(:,1);
            y   = es_pow(:,2);
            xx  = min(es_pow(:,1)):.001:max(es_pow(:,1));
            yy  = interp1(x,y,xx,'pchip');

            fig = figure('Renderer', 'painters', 'Position', [10 10 960 800],'visible','off');
            plot(xx,yy,'Color',[0.5 0.75 1],'LineWidth',2)
            grid on
            grid minor
            xlim([0.6,1.4])
            ylim([-1,101])
            xlabel('Effect Size (Relative Fitness)')
            ylabel('Power')
            hold on
            scatter(x, y,'MarkerEdgeColor',[0 .5 .5],...
                      'MarkerFaceColor',[0 .7 .7],...
                      'LineWidth',2);
            hold on
            title(sprintf('ES V/S Power\nTime = %dhrs | FPR = %.2f%% | Missing %d',cont_hrs,fpr,ss(iii)))
            hold off
            saveas(fig,sprintf('powes_%d_%d.png',cont_hrs,ss(iii)))

            fprintf('Power analysis done for %d missing references.\n',ss(iii))
            send_message(4124992194,'fi','gapEVAL',...
                sprintf('Power analysis done for %d missing references.',ss(iii)))
        end

        send_message(4124992194,'fi','gapEVAL','Task Complete!')
        
    catch me
        
        warning(me.message)
        send_message(4124992194,'fi','Error: gapEVAL',me.message)

    end


