
%%  INITIALIZATION

    cd /home/sbp29/MATLAB

    addpath('/home/sbp29/MATLAB/Matlab-Colony-Analyzer-Toolkit')
    addpath('/home/sbp29/MATLAB/bean-matlab-toolkit')
    addpath('/home/sbp29/MATLAB/sau-matlab-toolkit')
    addpath('/home/sbp29/MATLAB/sau-matlab-toolkit/grid-manipulation')
    addpath('/home/sbp29/MATLAB/paris')

    javaaddpath('/home/sbp29/MATLAB/mysql-connector-java-8.0.12.jar');
    
%%  DATA GATHER

    connectSQL;
    
    ss = [100,150,200,250,300];
    IL = 1;
    
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
 
    %%  TEMPORARY FITNESS AND NORM TABLES

        [data_fit, pos_miss] = LIGap(hours,n_plates,p2c_info,cont.name,...
                tablename_p2o,tablename_jpeg,pos.border,ss(iii),IL);

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
        
        pos_miss = sprintf('%d,',[pos_miss;pos_miss + 100000]);
        
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
            pvals{ii} = temp_p;
            stat{ii} = temp_s;

            len = length(pvals{ii});
            fpdat = [];

            for i = 1:length(p)
                fp = sum(pvals{ii} <= p(i));
                fpdat = [fpdat; [p(i), fp/len]];
            end

            fig = figure('Renderer', 'painters', 'Position', [10 10 960 600],'visible','off');
            histogram(pvals{ii}, 'Normalization', 'cdf')
            hold on
            plot(0:0.01:1,0:0.01:1,'--r','LineWidth',3)
            grid on
            xlabel('P Value Cut-offs')
            ylabel('Proportion of Colonies')
            title(sprintf('Time = %d hrs',hours(ii)))
            xlim([0,1])
            ylim([0,1])
            saveas(fig,sprintf('fpr_%d_%d.png',ss(iii),hours(ii)))
        end
        
    end



