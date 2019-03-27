%%  Sau MATLAB Colony Analyzer Toolkit
%
%%  nullNPVAL.m

%   Author: Saurin Parikh, March 2019
%   dr.saurin.parikh@gmail.com

%   Evaluation of the LI based control normalization
%   Null Distribution and pval cutoff proportions

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

        expt_name = '4C3_GA4';
        out_path = '/home/sbp29/MATLAB/4C3_Data/GA4/';
        density = 6144;

    %   MySQL Table Details  
        tablename_norm      = sprintf('%s_%d_NORM',expt_name,density);
        tablename_fit       = sprintf('%s_%d_FITNESS',expt_name,density);

        tablename_p2o       = '4C3_pos2orf_name4';
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

        p = 0:0.01:1;

    %%  NULL DISTRIBUTION
    %   The fitness distribution of the positions used to create the LI model

        for ii=1:length(hours)
            cont_data = fetch(conn, sprintf(['select * from %s ',...
                'where orf_name = ''%s'' ',...
                'and fitness is not NULL and hours = %d'],...
                tablename_fit,cont.name,hours(ii)));
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

            fig = figure('Renderer', 'painters', 'Position', [10 10 960 800],'visible','off');
            histogram(contp, 'NumBins', 20, 'Normalization', 'pdf');
            grid on
            xlabel('P Values')
            ylabel('Probability Density')
            title(sprintf('NULL DISTRIBUTION | Time = %d hrs',hours(ii)))
            saveas(fig,sprintf('%snulldist/%s_nullDIST_%d.png',...
                out_path,expt_name,hours(ii)))

            %%  DATA UNDER PVAL CUT-OFFS

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

            contmean = nanmean(cont_means);
            contstd = nanstd(cont_means);
            m = cont_means;
            tt = length(m);
            rest_dist =[];
            rest_means = [];

            for ss=4:8
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

                % proportion of colonies below a pval cutoff
                len = length(pvals{ii}{ss});
                fpdat = [];
                for i = 1:length(p)
                    fp = sum(pvals{ii}{ss} <= p(i));
                    fpdat = [fpdat; [p(i), fp/len]];
                end

                fig = figure('Renderer', 'painters', 'Position', [10 10 960 600],'visible','off');
                histogram(pvals{ii}{ss}, 'Normalization', 'cdf');
                hold on
                plot(0:0.01:1,0:0.01:1,'--r','LineWidth',3)
                grid on
                xlabel('P Value Cut-offs')
                ylabel('Proportion of Colonies')
                title(sprintf('Time = %d hrs | FP = %0.2f%%',hours(ii),fpdat(6,2)*100))
                xlim([0,1])
                ylim([0,1])
                saveas(fig,sprintf('%sprop/%s_cutoffPROP_%d_%d.png',...
                    out_path,expt_name,hours(ii),ss))
            end
            
            fprintf('nullNPVAL Calculation for %s at %d hrs is done.\n',...
                expt_name,hours(ii))
            %send_message(4124992194,'fi','nullNPVAL Update',...
            %sprintf("nullNPVAL Calculation for %s at %d hrs is done.",...
            %    expt_name,hours(ii)))
        
        end

        send_message(4124992194,'fi','nullNPVAL Complete',...
        sprintf("Null Distribution & Colony Proportion Measurements for %s Complete!",expt_name))

    catch me

        warning(me.message)
        send_message(4124992194,'fi','nullNPVAL Error',me.message)

    end     