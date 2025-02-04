%%  Sau MATLAB Colony Analyzer Toolkit
%
%%  TECHNICAL REPLICATES BASED POWER ANALYSIS USING BEAN METHOD
%   Power Analysis of Validation Experiments (4Control)
%   With technical replicates

%   Author: Saurin Parikh, June 2019
%   dr.saurin.parikh@gmail.com

    %%  Load Paths to Files and Data
    try

    load_toolkit;

    %%  Initialization

    %     Set preferences with setdbprefs.
        setdbprefs('DataReturnFormat', 'structure');
        setdbprefs({'NullStringRead';'NullStringWrite';'NullNumberRead';'NullNumberWrite'},...
                      {'null';'null';'NaN';'NaN'})
%         setdbprefs('NullNumberWrite','999999')

        expt_name = '4C4_FS_BEAN';
        expt = '4C4 FS BEAN';
        out_path = '/home/sbp29/MATLAB/4C3_Data/GA/S1Analysis/power/';
%         out_path = '/Users/saur1n/Desktop/4C3/Analysis/GA/S1Analysis/fnfp/';
        density = 6144;

    %   MySQL Table Details  
        tablename_norm      = sprintf('%s_%d_NORM',expt_name,density);
        tablename_fit       = sprintf('%s_%d_FITNESS',expt_name,density);
        tablename_pval       = sprintf('%s_%d_PVALUE',expt_name,density);

        tablename_p2o       = '4C4_pos2orf_name';
        tablename_bpos      = '4C4_borderpos';

        temp_norm      = sprintf('%s_TEMP_%d_NORM',expt_name,density);
        temp_fit       = sprintf('%s_TEMP_%d_FITNESS',expt_name,density);
        temp_fits      = sprintf('%s_TEMP_%d_FITNESS_STATS',expt_name,density);
        temp_pval      = sprintf('%s_TEMP_%d_PVALUE',expt_name,density);
        
    %   Reference Strain Name
        cont.name           = 'BF_control';

    %   MySQL Connection and fetch initial data
        connectSQL;

        p2c_info(1,:) = '4C4_pos2coor';
        p2c_info(2,:) = 'plate       ';
        p2c_info(3,:) = 'col         ';
        p2c_info(4,:) = 'row         ';

         p2c = fetch(conn, sprintf(['select * from %s a ',...
            'where density = %d ',...
            'order by a.%s, a.%s, a.%s'],...
            p2c_info(1,:),...
            density,...
            p2c_info(2,:),...
            p2c_info(3,:),...
            p2c_info(4,:)));

        n_plates = fetch(conn, sprintf(['select distinct %s from %s a ',...
            'where density = %d ',...
            'order by %s asc'],...
            p2c_info(2,:),...
            p2c_info(1,:),...
            density,...
            p2c_info(2,:)));

        hours = fetch(conn, sprintf(['select distinct hours from %s ',...
                 'order by hours asc'], tablename_fit));
        hours = hours.hours;
        
        cdata = [];
        contfitall = [];

        for atmpt = 1:3
            for t = 1:length(hours)
    %%  GENERATE FITNESS DATA

                cont_hrs = hours(t);
                rest_hrs = hours;%hours(hours~=cont_hrs);

                exec(conn, sprintf('drop table %s',temp_norm));
                exec(conn, sprintf(['create table %s ( ',...
                    'pos int(11) not NULL, ',...
                    'hours double not NULL, ',...
                    'bg double default NULL, ',...
                    'average double default NULL, ',...
                    'fitness double default NULL ',...
                    ')'],temp_norm));

                fprintf("Analysis for Control Hour = %0.1f Started.\n",cont_hrs);
                fprintf('Control Hour = %0.1f\n',cont_hrs)
    % 
    %             fpr = fetch(conn, sprintf(['select count(*) from %s ',...
    %                 'where hours = %d and p < 0.05'],tablename_pval,cont_hrs));
    %             orfs = fetch(conn, sprintf(['select count(*) from %s ',...
    %                 'where hours = %d'],tablename_pval,cont_hrs));
    %             fpr = (fpr.count____1/orfs.count____1)*100;

                data = [];
                fdata = [];
                for ii = 1:length(rest_hrs)
                    plate_fit = [];
                    for iii = 1:length(n_plates.plate)
                        pos.all = fetch(conn, sprintf(['select a.pos ',...
                            'from %s a ',...
                            'where %s = %d ',...
                            'and density = %d ',...
                            'order by %s, %s'],...
                            p2c_info(1,:),...
                            p2c_info(2,:),...
                            n_plates.plate(iii),...
                            density,...
                            p2c_info(3,:),...
                            p2c_info(4,:)));

                        pos.cont = fetch(conn, sprintf(['select a.pos ',...
                                'from %s a, %s b ',...
                                'where a.pos = b.pos and density = %d ',...
                                'and %s = %d and a.orf_name = ''%s'' ',...
                                'order by %s, %s'],...
                                tablename_p2o,...
                                p2c_info(1,:),...
                                density,...
                                p2c_info(2,:),...
                                n_plates.plate(iii),...
                                cont.name,...
                                p2c_info(3,:),...
                                p2c_info(4,:)));

                        cont_pos = col2grid(ismember(pos.all.pos, pos.cont.pos));
                        rest_pos = ~cont_pos;

                        cont_data = fetch(conn, sprintf(['select a.* ',...
                            'from %s a, %s b ',...
                            'where a.hours = %0.2f ',...
                            'and a.pos = b.pos and b.%s = %d ',...
                            'order by b.%s, b.%s'],...
                            tablename_fit,p2c_info(1,:),cont_hrs,...
                            p2c_info(2,:),n_plates.plate(iii),...
                            p2c_info(3,:),p2c_info(4,:)));
                        cont_avg = col2grid(cont_data.average).*cont_pos;

                        rest_data = fetch(conn, sprintf(['select a.* ',...
                            'from %s a, %s b ',...
                            'where a.hours = %0.2f ',...
                            'and a.pos = b.pos and b.%s = %d ',...
                            'order by b.%s, b.%s'],...
                            tablename_fit,p2c_info(1,:),rest_hrs(ii),...
                            p2c_info(2,:),n_plates.plate(iii),...
                            p2c_info(3,:),p2c_info(4,:)));
                        rest_avg = col2grid(rest_data.average).*rest_pos;
                        plate_avg = cont_avg + rest_avg;

                        plate_fit = col2grid(apply_correction( ...
                            grid2row(plate_avg), 'dim', 2, ...
                            InterleaveFilter(SpatialBorderMedian('SpatialFilter', ...
                            SpatialMedian('windowSize', 9))),..., 'windowShape', 'square'))), ...
                            PlateMode() ));

                        fdata = [fdata; pos.all.pos, ones(length(pos.all.pos),1).*rest_hrs(ii),...
                            grid2row(plate_avg)', grid2row(plate_avg)', grid2row(plate_fit)'];
                    end
                    fprintf('%0.1f hrs V/S %0.1f hrs done!\n', cont_hrs,rest_hrs(ii))
                end
                datainsert(conn, temp_norm,...
                        {'pos','hours','bg','average','fitness'},fdata);

    %             exec(conn, sprintf(['update %s ',...
    %                 'set fitness = NULL ',...
    %                 'where pos in ',...
    %                 '(select pos from %s)'],temp_norm, tablename_bpos));

                exec(conn, sprintf('drop table %s',temp_fit)); 
                exec(conn, sprintf(['create table %s ',...
                    '(select b.orf_name, a.pos, a.hours, a.bg, a.average, a.fitness ',...
                    'from %s a, %s b ',...
                    'where a.pos = b.pos ',...
                    'order by a.pos asc)'],temp_fit,temp_norm,tablename_p2o));

    %%  FITNESS STATS
                for rep = 2:2:16
                    rep_pos = combnk(1:16,rep);

                    if rep < 16 %16 for entire plate
                        N = datasample(1:length(rep_pos),1);
                        c = rep_pos(N,:);
                    else
                        c = rep_pos;
                    end

                    exec(conn, sprintf('drop table %s', temp_fits));
                    exec(conn, sprintf(['create table %s (orf_name varchar(255) null, ',...
                        'hours double not null, N int not null, cs_mean double null, ',...
                        'cs_median double null, cs_std double null)'],temp_fits));

                    colnames_fits = {'orf_name','hours','N','cs_mean','cs_median','cs_std'};

%                     stat_data = fit_stats(temp_fit);
% 
%                     tic
%                     datainsert(conn,temp_fits,colnames_fits,stat_data)
%                     toc
% 
%     %                 exec(conn, sprintf(['update %s ',...
%     %                     'set cs_mean = NULL, cs_median = NULL, cs_std = NULL ',...
%     %                     'where cs_mean = 999999'],temp_fits));
    
                    stat_data = fit_statsN(temp_fit,c);
                    stat_data.cs_mean(isnan(stat_data.cs_mean)) = 999999;
                    stat_data.cs_median(isnan(stat_data.cs_median)) = 999999;
                    stat_data.cs_std(isnan(stat_data.cs_std)) = 999999;

                    tic
                    datainsert(conn,temp_fits,colnames_fits,stat_data)
                    toc

                    exec(conn, sprintf(['delete from %s ',...
                        'where cs_mean = 999999'],temp_fits));


        %%  STATS TO PVALUES            
                    exec(conn, sprintf('drop table %s',temp_pval));
                    exec(conn, sprintf(['create table %s (orf_name varchar(255) null,'...
                        'hours double not null, p double null, stat double null)'],temp_pval));
                    colnames_pval = {'orf_name','hours','p','stat'};

                    contpos = fetch(conn, sprintf(['select pos from %s ',...
                        'where orf_name = ''%s'' and pos < 10000 ',...
                        'and pos not in ',...
                        '(select pos from %s)'],...
                        tablename_p2o,cont.name,tablename_bpos));
                    contpos = contpos.pos + [110000,120000,130000,140000,...
                        210000,220000,230000,240000,...
                        310000,320000,330000,340000,...
                        410000,420000,430000,440000];

    %                 contfitall = [contfitall, [ones(1,length(contfit))*cont_hrs;...
    %                         ones(1,length(contfit))*rest_hrs(iii);...
    %                         contfit]];

                    for iii = 1:length(rest_hrs)
                        contfit = [];
                        for ii = 1:length(contpos)
                            temp = fetch(conn, sprintf(['select fitness from %s ',...
                                'where hours = %0.2f and pos in (%s)'],temp_fit,rest_hrs(iii),...
                                sprintf('%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d',contpos(ii,:))));
                            if nansum(temp.fitness) > 0
                                outlier = isoutlier(temp.fitness);
                                temp.fitness(outlier) = NaN;
                                contfit = [contfit, nanmean(temp.fitness(c))];
                            end
                        end
                        contmean = nanmean(contfit);
                        contstd = nanstd(contfit);

                        orffit = fetch(conn, sprintf(['select orf_name, cs_median, ',...
                            'cs_mean, cs_std from %s ',...
                            'where hours = %0.2f and orf_name != ''%s'' ',...
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
                        tic
    %                         datainsert(conn,temp_pval,colnames_pval,pdata{iii})
                        sqlwrite(conn,temp_pval,struct2table(pdata{iii}))
                        toc
                    end

        %%  SAVING ALL DATA

                    temp_stat_p = fetch(conn, sprintf(['select a.*, b.p ',...
                            'from %s a, %s b ',...
                            'where a.orf_name = b.orf_name ',...
                            'and a.orf_name != ''BFC100'' ',...
                            'and a.hours = b.hours ',...
                            'order by a.hours, a.orf_name'], temp_fits, temp_pval));
                    writetable(temp_stat_p,...
                            sprintf('%s_%d_%d_%d_STATS_P.csv',expt_name,rep,...
                            t,atmpt),...
                            'Delimiter',',',...
                            'QuoteStrings',true)

                    temp_fitness = fetch(conn, sprintf(['select a.hours, a.pos, b.plate, b.col, b.row, ',...
                            'a.orf_name, a.bg, a.average, a.fitness ',...
                            'from %s a, %s b ',...
                            'where a.pos = b.pos ',...
                            'order by a.hours, b.%s, b.%s, b.%s'],...
                            temp_fit,p2c_info(1,:),...
                                        p2c_info(2,:),...
                                        p2c_info(3,:),...
                                        p2c_info(4,:)));        
                     writetable(temp_fitness,...
                            sprintf('%s_%d_%d_%d_FITNESS.csv',expt_name,rep,...
                            t,atmpt),...
                            'Delimiter',',',...
                            'QuoteStrings',true)

    %                 fprintf('techPowA for %s and %d (%s) replicates at %d hrs is done.\n',...
    %                     expt_name,rep,join(string(c),''),cont_hrs)
    %                     send_message(4124992194,'fi','techPowA Update',...
    %                         sprintf('techPowA for %s and %d (%s) replicates at %d hrs is done.\n',...
    %                             expt_name,rep,string(c),cont_hrs))

                end

                fprintf('techPowA for %s at %d hrs is done.\n',...
                    expt_name,cont_hrs)
    %             send_message(4124992194,'fi','techPowA Update',...
    %                 sprintf('techPowA for %s at %d hrs is done.\n',...
    %                     expt_name,cont_hrs))

            end
        end
        
%         writematrix(contfitall',...
%             sprintf('%s_CONTFITALL_%d.csv',expt_name,rep),...
%             'Delimiter',',',...
%             'QuoteStrings',true)
        
        fprintf("TechRep Based Power V/S Effect Size Analysis For %s Complete!\n",expt_name);    
%         send_message(4124992194,'fi','techPowA Complete',...
%             sprintf("TechRep Based Power V/S Effect Size Analysis For %s Complete!",expt_name))

    catch me

        warning(me.message)
        send_message(4124992194,'fi','techPowA Error',me.message)

    end   
