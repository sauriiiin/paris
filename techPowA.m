%%  Sau MATLAB Colony Analyzer Toolkit
%
%%  TECHNICAL REPLICATES BASED POWER ANALYSIS
%   Power Analysis of Validation Experiments (4Control)

%   Author: Saurin Parikh, March 2019
%   dr.saurin.parikh@gmail.com

    %%  Load Paths to Files and Data
    try

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

        expt_name = '4C3_GA1';
        expt = 'FS1-1';
        out_path = '/home/sbp29/MATLAB/4C3_Data/GA/S1Analysis/fnfp/';
%         out_path = '/Users/saur1n/Desktop/4C3/Analysis/GA/S1Analysis/fnfp/';
        density = 6144;

    %   MySQL Table Details  

        tablename_norm      = sprintf('%s_%d_NORM',expt_name,density);
        tablename_fit       = sprintf('%s_%d_FITNESS',expt_name,density);
        tablename_pval       = sprintf('%s_%d_PVALUE',expt_name,density);

        tablename_p2o       = '4C3_pos2orf_name1';
        tablename_bpos      = '4C3_borderpos';

        temp_norm      = sprintf('%s_TEMP_%d_NORM',expt_name,density);
        temp_fit       = sprintf('%s_TEMP_%d_FITNESS',expt_name,density);
        temp_fits      = sprintf('%s_TEMP_%d_FITNESS_STATS',expt_name,density);
        temp_pval      = sprintf('%s_TEMP_%d_PVALUE',expt_name,density);
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
        
        cdata = [];

        for t = 5:length(hours)
%%  GENERATE FITNESS DATA

            exec(conn, sprintf('drop table %s',temp_norm));
            exec(conn, sprintf(['create table %s ( ',...
                'pos int(11) not NULL, ',...
                'hours int(11) not NULL, ',...
                'bg double default NULL, ',...
                'average double default NULL, ',...
                'fitness double default NULL ',...
                ')'],temp_norm));
                
            cont_hrs = hours(t);
            rest_hrs = hours;%hours(hours~=cont_hrs);
            fprintf("Analysis for Control Hour = %0.1f Started.\n",cont_hrs);
            fprintf('Control Hour = %0.1f\n',cont_hrs)

            fpr = fetch(conn, sprintf(['select count(*) from %s ',...
                'where hours = %d and p < 0.05'],tablename_pval,cont_hrs));
            orfs = fetch(conn, sprintf(['select count(*) from %s ',...
                'where hours = %d'],tablename_pval,cont_hrs));
            fpr = (fpr.count___/orfs.count___)*100;
            
            data = [];
            fdata = [];
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
                        p2c_info(2,:),n_plates.x6144plate(iii),...
                        p2c_info(3,:),p2c_info(4,:)));

                    cont_avg = col2grid(cont_data.average).*cont_pos;
                    plate_bg = col2grid(cont_data.bg);

                    rest_data = fetch(conn, sprintf(['select a.* ',...
                        'from %s a, %s b ',...
                        'where a.hours = %d ',...
                        'and a.pos = b.pos and b.%s = %d ',...
                        'order by b.%s, b.%s'],...
                        tablename_fit,p2c_info(1,:),rest_hrs(ii),...
                        p2c_info(2,:),n_plates.x6144plate(iii),...
                        p2c_info(3,:),p2c_info(4,:)));

                    rest_avg = col2grid(rest_data.average).*rest_pos;
                    plate_avg = cont_avg + rest_avg;

                    plate_fit = plate_avg./plate_bg;
                    fdata = [fdata; pos.all.pos, ones(length(pos.all.pos),1).*rest_hrs(ii),...
                        grid2row(plate_bg)', grid2row(plate_avg)', grid2row(plate_fit)'];
                end
                fprintf('%0.1f hrs V/S %0.1f hrs done!\n', cont_hrs,rest_hrs(ii))
            end
            datainsert(conn, temp_norm,...
                    {'pos','hours','bg','average','fitness'},fdata);
                
            exec(conn, sprintf('drop table %s',temp_fit)); 
            exec(conn, sprintf(['create table %s ',...
                '(select b.orf_name, a.pos, a.hours, a.bg, a.average, a.fitness ',...
                'from %s a, %s b ',...
                'where a.pos = b.pos ',...
                'order by a.pos asc)'],temp_fit,temp_norm,tablename_p2o));
%%  FITNESS STATS            
            exec(conn, sprintf('drop table %s', temp_fits));
            exec(conn, sprintf(['create table %s (orf_name varchar(255) null, ',...
                'hours int not null, N int not null, cs_mean double null, ',...
                'cs_median double null, cs_std double null)'],temp_fits));

            colnames_fits = {'orf_name','hours','N','cs_mean','cs_median','cs_std'};

            stat_data = fit_stats(temp_fit);
            tic
            datainsert(conn,temp_fits,colnames_fits,stat_data)
            toc
%%  STATS TO PVALUES            
            exec(conn, sprintf('drop table %s',temp_pval));
            exec(conn, sprintf(['create table %s (orf_name varchar(255) null,'...
                'hours int not null, p double null, stat double null)'],temp_pval));
            colnames_pval = {'orf_name','hours','p','stat'};

            for iii = 1:length(rest_hrs)
                contpos = fetch(conn, sprintf(['select pos from %s ',...
                    'where orf_name = ''%s'' and pos < 10000 ',...
                    'and pos not in ',...
                    '(select pos from %s)'],...
                    tablename_p2o,cont.name,tablename_bpos));
                contpos = contpos.pos + [110000,120000,130000,140000,...
                    210000,220000,230000,240000];

                contfit = [];
                for ii = 1:length(contpos)
                    temp = fetch(conn, sprintf(['select fitness from %s ',...
                        'where hours = %d and pos in (%s) ',...
                        'and fitness is not null'],temp_fit,rest_hrs(iii),...
                        sprintf('%d,%d,%d,%d,%d,%d,%d,%d',contpos(ii,:))));
                    if isstruct(temp) ~= 0
                        if nansum(temp.fitness) > 0
                            contfit = [contfit, nanmean(temp.fitness)];
                        end
                    end
                end

                contmean = nanmean(contfit);
                contstd = nanstd(contfit);

                orffit = fetch(conn, sprintf(['select orf_name, cs_median, ',...
                    'cs_mean, cs_std from %s ',...
                    'where hours = %d and orf_name != ''%s'' ',...
                    'order by orf_name asc'],temp_fits,rest_hrs(iii),cont.name));

                m = contfit';
                tt = length(m);

                pvals = [];
                stat = [];
                for i = 1:length(orffit.orf_name)
                    if sum(m<orffit.cs_mean(i)) < tt/2
                        if m<orffit.cs_mean(i) == 0
                            pvals = [pvals; 1/tt];
                            stat = [stat; (orffit.cs_mean(i) - contmean)/contstd];
                        else
                            pvals = [pvals; ((sum(m<=orffit.cs_mean(i)))/tt)*2];
                            stat = [stat; (orffit.cs_mean(i) - contmean)/contstd];
                        end
                    else
                        pvals = [pvals; ((sum(m>=orffit.cs_mean(i)))/tt)*2];
                        stat = [stat; (orffit.cs_mean(i) - contmean)/contstd];
                    end
                end

                pdata{iii}.orf_name = orffit.orf_name;
                pdata{iii}.hours = ones(length(pdata{iii}.orf_name),1)*rest_hrs(iii);
                pdata{iii}.p = num2cell(pvals);
                pdata{iii}.p(cellfun(@isnan,pdata{iii}.p)) = {[]};
                pdata{iii}.stat = num2cell(stat);
                pdata{iii}.stat(cellfun(@isnan,pdata{iii}.stat)) = {[]};

                datainsert(conn,temp_pval,colnames_pval,pdata{iii})
            end
%%  CALCULATE POWER & ES
            cfit = fetch(conn, sprintf(['select fitness from %s ',...
                'where hours = %d and orf_name = ''%s'' ',...
                'and fitness is not NULL'],...
                temp_fit,cont_hrs,cont.name));
            cont_fit = cfit.fitness;

            for iii = 1:length(rest_hrs)
                rfit = fetch(conn, sprintf(['select fitness from %s ',...
                    'where hours = %d and orf_name != ''%s'' ',...
                    'and fitness is not NULL'],...
                    temp_fit,rest_hrs(iii),cont.name));
                rest_fit = rfit.fitness;
                ef_size = mean(rest_fit)/mean(cont_fit);
                
                pp = fetch(conn, sprintf(['select p from %s ',...
                    'where hours = %d '],...
                    temp_pval,rest_hrs(iii)));
                
                pow = (sum(pp.p<=0.05)/length(pp.p))*100;
                cdata = [cdata; ef_size, pow, abs(1-ef_size), 100-pow];
                data = [data; ef_size, pow, abs(1-ef_size), 100-pow]; 
            end

        %%  POWER vs ES

%             [~, i] = sort(data(:,1));
%             es_pow = data(i, :);
% 
%             x   = es_pow(:,1);
%             y   = es_pow(:,2);
%             xx  = min(es_pow(:,1)):.001:max(es_pow(:,1));
%             yy  = interp1(x,y,...
%                 xx,'pchip');
% 
%     %         figure()
%             fig = figure('Renderer', 'painters', 'Position', [10 10 960 800],'visible','off');
%             plot(xx,yy,'Color',[0.5 0.75 1],'LineWidth',2)
%             grid on
%             grid minor
%             xlim([0.6,1.4])
%             ylim([-1,101])
%             xlabel('Effect Size (Relative Fitness)')
%             ylabel('Power')
%             hold on
%             scatter(x, y,'MarkerEdgeColor',[0 .5 .5],...
%                       'MarkerFaceColor',[0 .7 .7],...
%                       'LineWidth',2);
%             hold on
%             title(sprintf('ES V/S Power\nTime = %dhrs | SS = %d | FPR = %.2f%%',cont_hrs, ss, fpr))
%             hold off
%             saveas(fig,sprintf('%s%s_TpowES_%d_%d.png',...
%                 out_path,expt_name,cont_hrs,ss))
            
%             
%             [~, i] = sort(data(:,3));
%             es_fn = data(i, :);
%             x   = log10(es_fn(:,3));
%             y   = es_fn(:,4);
% 
%     %         figure()
%             fig = figure('Renderer', 'painters', 'Position', [10 10 960 800],'visible','off');
%             grid on
%             grid minor
%             ylim([-1,101])
%             xlabel('Log10(Effect Size)')
%             ylabel('False Negative Rate')
%             hold on
%             scatter(x, y,'MarkerEdgeColor',[0 .5 .5],...
%                       'MarkerFaceColor',[0 .7 .7],...
%                       'LineWidth',2);
%             hold on
%             title(sprintf('ES V/S FN\nTime = %dhrs | SS = %d | FPR = %.2f%%',cont_hrs, ss, fpr))
%             hold off
%             saveas(fig,sprintf('%s%s_TFNES_%d_%d.png',...
%                 out_path,expt_name,cont_hrs,ss))
                    
            fprintf('techPowA for %s at %d hrs is done.\n',...
                expt_name,cont_hrs)
            send_message(4124992194,'fi','techPowA Update',...
                sprintf('techPowA for %s at %d hrs is done.\n',...
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
%         saveas(fig,sprintf('%s%s_TpowES_%d.png',out_path,expt_name,ss))
        
        
        [~, i] = sort(cdata(:,3));
        es_fn = cdata(i, :);
        x   = log10(es_fn(:,3));
        y   = es_fn(:,4);

%         figure()
        fig = figure('Renderer', 'painters', 'Position', [10 10 960 800],'visible','off');
        grid on
        grid minor
        ylim([-1,101])
        xlabel('Log10(| 1- Effect Size |)')
        ylabel('False Negative Rate')
        hold on
        scatter(x, y,'MarkerEdgeColor',[0 .5 .5],...
                  'MarkerFaceColor',[0 .7 .7],...
                  'LineWidth',2);
        hold on
        title(sprintf('ES V/S FN\n%s',expt))
        hold off
        saveas(fig,sprintf('%s%s_TFNES_%d.png',out_path,expt_name,ss))

        fprintf("TechRep Based Power V/S Effect Size Analysis For %s Complete!\n",expt_name);
        send_message(4124992194,'fi','techPowA Complete',...
            sprintf("TechRep Based Power V/S Effect Size Analysis For %s Complete!",expt_name))

    catch me

        warning(me.message)
        send_message(4124992194,'fi','techPowA Error',me.message)

    end   
