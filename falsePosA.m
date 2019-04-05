%%  Sau MATLAB Colony Analyzer Toolkit
%
%%  TECHNICAL REPLICATES BASED POWER ANALYSIS
%   FALSE POSITIVE of Validation Experiments (4Control)

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

        javaaddpath('/home/sbp29/MATLAB/mysql-connector-java-8.0.15.jar');

    %%  Initialization

    %     Set preferences with setdbprefs.
        setdbprefs('DataReturnFormat', 'structure');
        setdbprefs({'NullStringRead';'NullStringWrite';'NullNumberRead';'NullNumberWrite'},...
                      {'null';'null';'NaN';'NaN'})

        expt_name = '4C3_GA1';
        expt = 'FS1-1';
%         out_path = '/home/sbp29/MATLAB/4C3_Data/GA/S1Analysis/power/';
        out_path = '/Users/saur1n/Desktop/4C3/Analysis/GA/S1Analysis/fnfp/';
        density = 6144;
        
        fprintf("Analysis for %s Started.\n",expt_name);

    %   MySQL Table Details  

        tablename_pval      = sprintf('%s_%d_PVALUE',expt_name,density);
        tablename_es        = sprintf('%s_%d_FITNESS_ES',expt_name,density);

    %   MySQL Connection and fetch initial data
        connectSQL;
        hours = fetch(conn, sprintf(['select distinct hours from %s ',...
                 'order by hours asc'], tablename_pval));
        hours = hours.hours;
        
        ul = 0.005:0.005:0.120;

        for t = 5:length(hours)
%%  GET EFFECT SIZE DATA

            fprintf("Analysis for Hour = %0.1f Started.\n",hours(t));

            alldat = fetch(conn, sprintf(['select a.orf_name,a.effect_size,b.p ',...
                'from %s a, %s b ',...
                'where a.hours = %d ',...
                'and a.hours = b.hours and a.orf_name = b.orf_name'],...
                tablename_es,tablename_pval,hours(t)));
            alldat.es_perc = abs(1-alldat.effect_size);
            
            data = [];
            for i = 1:length(ul)
                temp = alldat.p(alldat.es_perc <= ul(i));
                N = length(temp);
                n = sum(temp <= 0.05);
                fp = n/N*100;
                data = [data,fp];
            end
       
            N = length(alldat.p);
            n = sum(alldat.p <= 0.05);
            max = n/N*100;
            
            x   = ul;
            y   = data;
            yy = smooth(x,y,'rloess');

%             figure()
            fig = figure('Renderer', 'painters', 'Position', [10 10 960 800],'visible','off');
            plot(x,yy,'r--','LineWidth',3)
            grid on
            grid minor
            ylim([-max*0.01,max*1.01])
            xlabel('| 1 - Effect Size |')
            ylabel('False Positive Rate')
            hold on
            scatter(x, y,'MarkerEdgeColor',[0 .5 .5],...
                      'MarkerFaceColor',[0 .7 .7],...
                      'LineWidth',2);
%             legend('fitted curve','real data')
            hold on
            title(sprintf('ES V/S FP\n%s | Time = %dhrs | FPR = %.2f%%',...
                expt,hours(t), max))
            hold off
            saveas(fig,sprintf('%s%s_TFPES_%d.png',...
                out_path,expt_name,hours(t)))
                    
            fprintf('falsePosA for %s at %d hrs is done.\n',...
                expt_name,hours(t))
        end   
%% 
        fprintf("TechRep Based ES V/S FP Analysis For %s Complete!\n",expt_name);
        send_message(4124992194,'fi','falsePosA Complete',...
            sprintf("TechRep Based Power V/S Effect Size Analysis For %s Complete!",expt_name))

    catch me

        warning(me.message)
        send_message(4124992194,'fi','falsePosA Error',me.message)

    end
    
    