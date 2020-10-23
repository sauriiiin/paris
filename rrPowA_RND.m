%%  Sau MATLAB Colony Analyzer Toolkit
%
%%  POWER ANALYSIS VARRYING TECH REP AND REF PROP - RND VIRTUAL PLATES
%   Power Analysis of Validation Experiments (4Control)
%   With change in Reference Proportion & technical replicates
%   specifically designed for the virtual plates with random cs
%   distribution

%   Author: Saurin Parikh, July 2020
%   dr.saurin.parikh@gmail.com

    %%  Load Paths to Files and Data
    try

        load_toolkit;

    %%  Initialization

        expt_name = '4C4_FS_RND2';
        expt      = '4C4 FS RND2';
        density   = 6144;

    %   MySQL Table Details  
        tablename_norm = sprintf('%s_%d_NORM',expt_name,density);
        tablename_fit  = sprintf('%s_%d_FITNESS',expt_name,density);
        tablename_pval = sprintf('%s_%d_PVALUE',expt_name,density);

        tablename_p2o  = '4C4_pos2orf_name';
        tablename_bpos = '4C4_borderpos';

        temp_fit       = tablename_fit;
        temp_fits      = sprintf('%s_TEMP_%d_FITNESS_STATS',expt_name,density);
        temp_pval      = sprintf('%s_TEMP_%d_PVALUE',expt_name,density);
        
    %   Reference Strain Name
        cont.name = 'BF_control';

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
        
        for atmpt = 1:5
    %%  FITNESS STATS
            max_rep = 16;
            for rep = 2:2:16
                rep_pos = combnk(1:max_rep,rep);

                if rep < max_rep %max_rep for entire plate
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
%                     contpos = contpos.pos + [110000,120000,130000,140000,...
%                         210000,220000,230000,240000];
                contfit = [];
                for iii = 1:length(hours)
                    for ii = 1:length(contpos)
                        temp = fetch(conn, sprintf(['select fitness from %s ',...
                            'where hours = %0.2f and pos in (%s)'],temp_fit,hours(iii),...
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
                        'order by orf_name asc'],temp_fits,hours(iii),cont.name));

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
                    pdata{iii}.hours = ones(length(pdata{iii}.orf_name),1)*hours(iii);
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

                t = 99;
                temp_stat_p = fetch(conn, sprintf(['select c.*, d.rnd_hrs from  ',...
                    '(select a.*, b.p from %s a, %s b ',...
                    'where a.orf_name = b.orf_name and a.orf_name != ''BFC100'' ',...
                    'and a.hours = b.hours order by a.hours, a.orf_name) c ',...
                    'left join ',...
                    '(select orf_name, hours, round(avg(rnd_hrs),2) rnd_hrs ',...
                    'from 4C4_FS_RND2_6144_DATA ',...
                    'group by orf_name, hours ',...
                    'order by hours, orf_name) d ',...
                    'on c.orf_name = d.orf_name and c.hours = d.hours'], temp_fits, temp_pval));
                
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

%                 fprintf('techPowA for %s and %d (%s) replicates at %0.2f hrs is done.\n',...
%                     expt_name,rep,join(string(c),''),cont_hrs)
%                     send_message(4124992194,'fi','techPowA Update',...
%                         sprintf('techPowA for %s and %d (%s) replicates at %d hrs is done.\n',...
%                             expt_name,rep,string(c),cont_hrs))
            end
            
        end
        
        fprintf("TechRep Based Power V/S Effect Size Analysis For %s Complete!\n",expt_name);    
        send_message(4124992194,'fi','techPowA Complete',...
            sprintf("TechRep Based Power V/S Effect Size Analysis For %s Complete!",expt_name))

    catch me

        warning(me.message)
        send_message(4124992194,'fi','techPowA Error',me.message)

    end   
