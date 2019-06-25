%%  Sau MATLAB Colony Analyzer Toolkit
%
%%  POWER ANALYSIS
%   Power Analysis of Validation Experiments (4Control)

%   Author: Saurin Parikh, March 2019
%   dr.saurin.parikh@gmail.com

    %%  Load Paths to Files and Data
    try

        cd /home/sbp29/MATLAB

        addpath(genpath('/home/sbp29/MATLAB/Matlab-Colony-Analyzer-Toolkit'))
        addpath(genpath('/home/sbp29/MATLAB/bean-matlab-toolkit'))
        addpath(genpath('/home/sbp29/MATLAB/sau-matlab-toolkit'))
        addpath(genpath('/home/sbp29/MATLAB/sau-matlab-toolkit/grid-manipulation'))
        addpath(genpath('/home/sbp29/MATLAB/paris'))
        addpath(genpath('/home/sbp29/MATLAB/development'))

        javaaddpath('/home/sbp29/MATLAB/mysql-connector-java-8.0.16.jar');

    %%  Initialization

    %     Set preferences with setdbprefs.
        setdbprefs('DataReturnFormat', 'structure');
        setdbprefs({'NullStringRead';'NullStringWrite';'NullNumberRead';'NullNumberWrite'},...
                      {'null';'null';'NaN';'NaN'})
        setdbprefs('NullNumberWrite','999999')

        expt_name = '4C3_GA1_TRBLBR';
        expt = 'FS1-1-TRBLBR';
        out_path = '/home/sbp29/MATLAB/4C3_Data/GA/S1Analysis/power/';
%         out_path = '/Users/saur1n/Desktop/4C3/Analysis/GA/S1Analysis/fnfp/';
        density = 6144;

    %   MySQL Table Details  

        tablename_norm      = sprintf('%s_%d_NORM',expt_name,density);
        tablename_fit       = sprintf('%s_%d_FITNESS',expt_name,density);

        tablename_p2o       = '4C3_TRBLBR_pos2orf_name';
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
                 'order by hours asc'], tablename_fit));
        hours = hours.hours;
        
        fpr = 0.05;
        cdata = [];

        for t = 6:length(hours)
            cont_hrs = hours(t);
            rest_hrs = hours(6:length(hours));%hours(hours~=cont_hrs);
            fprintf("Analysis for Control Hour = %0.1f Started.\n",cont_hrs);
            
            for ss = 2:2:20%4:8
                fprintf('Control Hour = %0.1f | Sample Size = %d\n',cont_hrs,ss)

%                 fpr = fpr4c(tablename_fit, cont.name, cont_hrs, ss);
                p5 = pforfpr(tablename_fit, cont.name, cont_hrs, fpr, ss);
                
                data = [];

                for ii = 1:length(rest_hrs)
                    plate_fit = [];
                    cont_fit = [];
                    rest_fit = [];
                    for iii = 1:length(n_plates.x6144plate_1)
                        pos.all = fetch(conn, sprintf(['select a.pos ',...
                            'from %s a ',...
                            'where %s = %d ',...
                            'order by %s, %s'],...
                            p2c_info(1,:),...
                            p2c_info(2,:),...
                            n_plates.x6144plate_1(iii),...
                            p2c_info(3,:),...
                            p2c_info(4,:)));

                        pos.cont = fetch(conn, sprintf(['select a.pos ',...
                            'from %s a, %s b ',...
                            'where a.pos = b.pos and %s = %d and a.orf_name = ''%s'' ',...
                            'order by %s, %s'],...
                            tablename_p2o,...
                            p2c_info(1,:),...
                            p2c_info(2,:),...
                            n_plates.x6144plate_1(iii),...
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
                            p2c_info(2,:),n_plates.x6144plate_1(iii),...
                            p2c_info(3,:),p2c_info(4,:)));

                        cont_avg = col2grid(cont_data.average).*cont_pos;
                        plate_bg = col2grid(cont_data.bg);

                        rest_data = fetch(conn, sprintf(['select a.* ',...
                            'from %s a, %s b ',...
                            'where a.hours = %d ',...
                            'and a.pos = b.pos and b.%s = %d ',...
                            'order by b.%s, b.%s'],...
                            tablename_fit,p2c_info(1,:),rest_hrs(ii),...
                            p2c_info(2,:),n_plates.x6144plate_1(iii),...
                            p2c_info(3,:),p2c_info(4,:)));

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
                        cont_dist(i,:) = datasample(cont_fit, ss, 'Replace', false);
                        cont_means(i,:) = mean(cont_dist(i,~isoutlier(cont_dist(i,:))));
                    end

                    rest_dist =[];
                    rest_means = [];
                    for i=1:100000
                        rest_dist(i,:) = datasample(rest_fit, ss, 'Replace', false);
                        rest_means(i,:) = mean(rest_dist(i,~isoutlier(rest_dist(i,:))));
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

                    pp = sum(temp_p<p5);
                    pow = (pp/length(rest_means))*100;
%                     avg_diff = abs(nanmean(nanmean(cont_avg)) - nanmean(nanmean(rest_avg)));
                    
                    neu = length(rest_means) - pp;
                    del = neu + sum(temp_s(temp_p < p5) < 0);
                    ben = del + sum(temp_s(temp_p < p5) > 0);
                    
                    neu_p = length(rest_means) - sum(temp_p < 0.05);
                    del_p = neu_p + sum(temp_s(temp_p < 0.05) < 0);
                    ben_p = del_p + sum(temp_s(temp_p < 0.05) > 0);
                    
                    cdata = [cdata; ss, cont_hrs, rest_hrs(ii), ef_size, pow,...
                        neu, del, ben,...
                        neu_p, del_p, ben_p];
%                     data = [data; ef_size, pow, avg_diff];

% %                     figure()
%                     fig = figure('Renderer', 'painters', 'Position', [10 10 480 300],'visible','off');
%                     [f,xi] = ksdensity(cont_means);
%                     plot(xi,f,'LineWidth',3)
%                     xlim([0.75,1.25])
%                     ylim([0,30])
%                     hold on
%                     [f,xi] = ksdensity(rest_means);
%                     plot(xi,f,'LineWidth',3)
%                     legend('control','rest of plate')
%                     title(sprintf(['TimeC = %d hrs | TimeR = %d hrs \n ',...
%                         'ES = %0.3f | Power = %0.3f'],...
%                         cont_hrs,rest_hrs(ii),ef_size,pow))
%                     xlabel('Fitness')
%                     ylabel('Density')
%                     grid on
%                     grid minor
%                     hold off
%                     saveas(fig,sprintf('vp_powes_%d_%d.png',cont_hrs,rest_hrs(ii)))
%                     saveas(fig,sprintf('%s%s_ContRest_%d%d_%d.png',...
%                         out_path,expt_name,cont_hrs,rest_hrs(ii),ss))
                    fprintf('%0.1f hrs V/S %0.1f hrs and %d replicates done!\n',...
                        cont_hrs,rest_hrs(ii),ss)
                end

            %%  POWER vs ES
% 
%                 [~, i] = sort(data(:,1));
%                 es_pow = data(i, :);
% 
%                 x{cont_hrs}{ss}   = es_pow(:,1);
%                 y{cont_hrs}{ss}   = es_pow(:,2);
%                 xx{cont_hrs}{ss}  = min(es_pow(:,1)):.001:max(es_pow(:,1));
%                 yy{cont_hrs}{ss}  = interp1(x{cont_hrs}{ss},y{cont_hrs}{ss},...
%                     xx{cont_hrs}{ss},'pchip');
% 
%         %         figure()
%                 fig = figure('Renderer', 'painters', 'Position', [10 10 960 800],'visible','off');
%                 plot(xx{cont_hrs}{ss},yy{cont_hrs}{ss},'Color',[0.5 0.75 1],'LineWidth',2)
%                 grid on
%                 grid minor
%                 xlim([0.6,1.4])
%                 ylim([-1,101])
%                 xlabel('Effect Size (Relative Fitness)')
%                 ylabel('Power')
%                 hold on
%                 scatter(x{cont_hrs}{ss}, y{cont_hrs}{ss},'MarkerEdgeColor',[0 .5 .5],...
%                           'MarkerFaceColor',[0 .7 .7],...
%                           'LineWidth',2);
%                 hold on
%                 title(sprintf('ES V/S Power\nTime = %dhrs | SS = %d | FPR = %.2f%%',cont_hrs, ss, fpr))
%                 hold off
%                 saveas(fig,sprintf('%s%s_powES_%d_%d.png',...
%                     out_path,expt_name,cont_hrs,ss))

                fprintf('powAnalysis for %s at cont_hrs %0.1f and %d replicates is done.\n',...
                    expt_name,cont_hrs,ss)
            
            end
                    
            fprintf('powAnalysis for %s at %d hrs is done.\n',...
                expt_name,cont_hrs)
            send_message(4124992194,'fi','powAnalysis Update',...
                sprintf('powAnalysis for %s at %d hrs is done.\n',...
                    expt_name,cont_hrs))
        end
%%  COMPOSITE ES AND POW RELATIONSHIP
            
%         [~, i] = sort(cdata(:,1));
%         es_pow = cdata(i, :);
% 
%         x   = es_pow(:,1);
%         y   = es_pow(:,2);
%         xx  = min(es_pow(:,1)):.001:max(es_pow(:,1));
%         yy  = interp1(x,y,xx,'pchip');
% 
% %         figure()
%         fig = figure('Renderer', 'painters', 'Position', [10 10 960 800],'visible','off');
%         plot(xx,yy,'Color',[0.5 0.75 1],'LineWidth',2)
%         grid on
%         grid minor
%         xlim([0.6,1.4])
%         ylim([-1,101])
%         xlabel('Effect Size (Relative Fitness)')
%         ylabel('Power')
%         hold on
%         scatter(x, y,'MarkerEdgeColor',[0 .5 .5],...
%                   'MarkerFaceColor',[0 .7 .7],...
%                   'LineWidth',2);
%         hold on
%         title(sprintf('%s\nES V/S Power',expt))
%         hold off
%         saveas(fig,sprintf('%s%s_powES_%d.png',out_path,expt_name,ss))

        csvwrite(sprintf('%s_POWANA.csv',expt_name),cdata)

        fprintf("Power V/S Effect Size Analysis For %s Complete!\n",expt_name);
        send_message(4124992194,'fi','powAnalysis Complete',...
            sprintf("Power V/S Effect Size Analysis For %s Complete!",expt_name))

    catch me

        warning(me.message)
        send_message(4124992194,'fi','powAnalysis Error',me.message)

    end   
