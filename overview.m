%%  Sau MATLAB Colony Analyzer Toolkit
%
%%  overview.m

%   Author: Saurin Parikh, March 2019
%   dr.saurin.parikh@gmail.com

%   Hourwise-platewise RMSE overview of the results

%%  Load Paths to Files and Data

    cd /home/sbp29/MATLAB

    addpath('/home/sbp29/MATLAB/Matlab-Colony-Analyzer-Toolkit')
    addpath('/home/sbp29/MATLAB/bean-matlab-toolkit')
    addpath('/home/sbp29/MATLAB/sau-matlab-toolkit')
    addpath('/home/sbp29/MATLAB/sau-matlab-toolkit/grid-manipulation')
    addpath('/home/sbp29/MATLAB/paris')
    addpath('/home/sbp29/MATLAB/development')

    javaaddpath('/home/sbp29/MATLAB/mysql-connector-java-8.0.12.jar');

%%  Add MCA toolkit to Path
%     add_mca_toolkit_to_path

%%  Initialization

%     Set preferences with setdbprefs.
    setdbprefs('DataReturnFormat', 'structure');
    setdbprefs({'NullStringRead';'NullStringWrite';'NullNumberRead';'NullNumberWrite'},...
                  {'null';'null';'NaN';'NaN'})

    expt_name = '4C3_GA_NIL';
    out_path = '/home/sbp29/MATLAB/4C3_Data/GA/S1NILAnalysis/overview/';
    density = 6144;
    
%   MySQL Table Details  
    tablename_fit       = sprintf('%s_%d_FITNESS',expt_name,density);

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
    
%%  PLATEWISE RMSE

    for i = 1:length(hours)
        for iii = 1:length(n_plates.x6144plate)
            clear rmse
            bg = fetch(conn, sprintf(['select a.* ',...
                'from %s a, %s b ',...
                'where a.hours = %d ',...
                'and a.pos = b.pos ',...
                'and b.%s = %d ',...
                'order by b.%s, b.%s'],...
                tablename_fit,p2c_info(1,:),hours(i),p2c_info(2,:),...
                iii,p2c_info(3,:),p2c_info(4,:)));

            rmse(:,iii) = (bg.average - bg.bg).^2;

            max_avg = max(bg.average);
            min_avg = min(bg.average);

            fig = figure('Renderer', 'painters', 'Position', [10 10 1920 1200],'visible','off');
%            figure()
            subplot(2,2,1)
            heatmap(col2grid(bg.average),'ColorLimits',[min_avg max_avg]);
            title(sprintf('Observed Pixel Count\n(Plate %d, %d hr)',iii,hours(i)))
            title('Observed Pixel Count')
            subplot(2,2,2)
            heatmap(col2grid(bg.bg),'ColorLimits',[min_avg max_avg]);
%             title(sprintf('Predicted Pixel Count\n(Plate %d, %d hr)',iii,hours(i)))
            title('Predicted Pixel Count')
            subplot(2,2,3)
            heatmap(col2grid(bg.fitness),'ColorLimits',[0.7 1.4]);
            title('Fitness')
            subplot(2,2,4)
            heatmap(col2grid(abs(bg.average - bg.bg)),'ColorLimits',[0 120]);
            title(sprintf('RMSE (%0.3f)',sqrt(nanmean(rmse(:,iii)))))
            colormap parula
            saveas(fig,sprintf('%s%s_OVW_%d_%d.png',...
                out_path,expt_name,hours(i),iii))
        end
    end
