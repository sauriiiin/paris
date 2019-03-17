
%%  VIRTUAL PLATE POWER ANALYSIS AND FDR
%   Benchmark FPR and then use it to calculate FDR for each virtual plate

    %%  Load Paths to Files and Data

    cd /home/sbp29/MATLAB

    addpath('/home/sbp29/MATLAB/Matlab-Colony-Analyzer-Toolkit')
    addpath('/home/sbp29/MATLAB/bean-matlab-toolkit')
    addpath('/home/sbp29/MATLAB/sau-matlab-toolkit')
    addpath('/home/sbp29/MATLAB/sau-matlab-toolkit/grid-manipulation')
    addpath('/home/sbp29/MATLAB/paris')
    addpath('/home/sbp29/MATLAB/development')

    javaaddpath('/home/sbp29/MATLAB/mysql-connector-java-8.0.12.jar');

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

    for ss=1:8
        
        fprintf('sample size = %d\n',ss)
        cont_hrs = 14;
        rest_hrs = hours(hours~=cont_hrs);
        fpr = fpr4c(tablename_fit, cont.name, cont_hrs, ss);
        
        data = [];

        for ii = 1:length(rest_hrs)
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

%             s = (((length(cont_fit)*(std(cont_fit))^2 +...
%                 length(rest_fit)*(std(rest_fit))^2)/...
%                 (length(cont_fit) +...
%                 length(rest_fit) - 2))^(0.5));
% 
%             ef_size = abs(mean(cont_fit) - mean(rest_fit))/s;
            ef_size = mean(rest_fit)/mean(cont_fit);
            
            all = length(rest_means);
            pp = sum(temp_p<0.05);
%             fp = length(rest_means)*fpr/100;
            pow = (pp/all)*100;
            fdr = 0;
            
%             if fp >= pp
%                 tp = 0;
%                 pow = 0;
%                 fdr = 100;
%             else
%                 tp = pp - fp;
%                 pow = (tp/all)*100;
%                 fdr = (fp/pp)*100;
%             end

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
        
%         yr{cont_hrs}{ss}  = es_pow(:,3);
%         yyr{cont_hrs}{ss}  = interp1(x{cont_hrs}{ss},yr{cont_hrs}{ss},...
%             xx{cont_hrs}{ss},'pchip');

%         figure()
        fig = figure('Renderer', 'painters', 'Position', [10 10 960 800],'visible','off');
%         yyaxis left
        plot(xx{cont_hrs}{ss},yy{cont_hrs}{ss},'Color',[0.5 0.75 1],'LineWidth',2)
        grid on
        grid minor
        xlim([0.6,1.4])
        ylim([-1,101])
        xlabel('Effect Size (Relative Fitness)')
        ylabel('Power')
        hold on
        scatter(x{cont_hrs}{ss}, y{cont_hrs}{ss},'MarkerEdgeColor',[0 .5 .5],...
                  'MarkerFaceColor',[0 .7 .7],...
                  'LineWidth',2);
        hold on
%         yyaxis right
%         plot(xx{cont_hrs}{ss},yyr{cont_hrs}{ss},'Color',[0.91 0.41 0.17],'LineWidth',2)
%         ylim([-1,101])
%         ylabel('False Discovery')
%         hold on
%         scatter(x{cont_hrs}{ss}, yr{cont_hrs}{ss},'MarkerEdgeColor',[0 .6 .6],...
%                   'MarkerFaceColor',[0 1 1],...
%                   'LineWidth',2);
        title(sprintf('ES V/S Power\nTime = %dhrs | SS = %d | FPR = %.2f%%',cont_hrs, ss, fpr))
        hold off
        saveas(fig,sprintf('vp_powes(rf)_%d_%d.png',cont_hrs,ss))
    end
    
    
