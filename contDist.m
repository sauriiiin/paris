%%  Sau MATLAB Colony Analyzer Toolkit
%
%%  TECHNICAL REPLICATE BASED CONTROL DISTRIBUTION
%   Control Distribution of Validation Experiments (4Control)

%   Author: Saurin Parikh, April 2019
%   dr.saurin.parikh@gmail.com

%%  Load Paths to Files and Data

%     cd /home/sbp29/MATLAB
% 
%     addpath('/home/sbp29/MATLAB/Matlab-Colony-Analyzer-Toolkit')
%     addpath('/home/sbp29/MATLAB/bean-matlab-toolkit')
%     addpath('/home/sbp29/MATLAB/sau-matlab-toolkit')
%     addpath('/home/sbp29/MATLAB/sau-matlab-toolkit/grid-manipulation')
%     addpath('/home/sbp29/MATLAB/paris')
%     addpath('/home/sbp29/MATLAB/development')
% 
%     javaaddpath('/home/sbp29/MATLAB/mysql-connector-java-8.0.15.jar');

%%  Initialization

%     Set preferences with setdbprefs.
    setdbprefs('DataReturnFormat', 'structure');
    setdbprefs({'NullStringRead';'NullStringWrite';'NullNumberRead';'NullNumberWrite'},...
                  {'null';'null';'NaN';'NaN'})

    expt_name = '4C3_GA1';
    expt = 'FS1-1';
%     out_path = '/home/sbp29/MATLAB/4C3_Data/GA1/contdist/';
    out_path = '/Users/saur1n/Desktop/4C3/Analysis/GA/S1Analysis/contdist/';
    density = 6144;

%   MySQL Table Details  

%     tablename_norm      = sprintf('%s_%d_NORM',expt_name,density);
    tablename_fit       = sprintf('%s_%d_FITNESS',expt_name,density);
%     tablename_pval       = sprintf('%s_%d_PVALUE',expt_name,density);

    tablename_p2o       = '4C3_pos2orf_name1';
    tablename_bpos      = '4C3_borderpos';
    
%   Reference Strain Name
    cont.name           = 'BF_control';

%   MySQL Connection and fetch initial data
    connectSQL;
    
    hours = fetch(conn, sprintf(['select distinct hours from %s ',...
                 'order by hours asc'], tablename_fit));
    hours = hours.hours;

%%

    contpos = fetch(conn, sprintf(['select pos from %s ',...
        'where orf_name = ''%s'' and pos < 10000 ',...
        'and pos not in ',...
        '(select pos from %s)'],...
        tablename_p2o,cont.name,tablename_bpos));
    contpos = contpos.pos + [110000,120000,130000,140000,...
        210000,220000,230000,240000];

    for iii = 1:length(hours)
        contfit = [];
        for ii = 1:length(contpos)
            temp = fetch(conn, sprintf(['select fitness from %s ',...
                'where hours = %d and pos in (%s) ',...
                'and fitness is not null'],tablename_fit,hours(iii),...
                sprintf('%d,%d,%d,%d,%d,%d,%d,%d',contpos(ii,:))));
            if nansum(temp.fitness) > 0
                contfit = [contfit, nanmean(temp.fitness)];
            end
        end

        contmean = nanmean(contfit);
        contstd = nanstd(contfit);

%         fig = figure('Renderer', 'painters', 'Position', [10 10 480 300],'visible','off');
        figure()
        [f,xi] = ksdensity(contfit);
        plot(xi,f,'LineWidth',3)
        grid on
        grid minor
        xlabel('Fitness')
        ylabel('Density')
        title(sprintf('%s | Control Distribution\nMean = %0.4f',...
            expt,contmean))
%         saveas(fig,sprintf('%s%s_ContDist_%d.png',...
%                     out_path,expt_name,cont_hrs))
    end