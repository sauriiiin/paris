%%  Sau MATLAB Colony Analyzer Toolkit
%
%%  plateeffectP.m

%   Author: Saurin Parikh, February 2019
%   dr.saurin.parikh@gmail.com

%   Observe Plate Effect

%%  INITIALIZATION

    cd /home/sbp29/MATLAB

    addpath('/home/sbp29/MATLAB/Matlab-Colony-Analyzer-Toolkit')
    addpath('/home/sbp29/MATLAB/bean-matlab-toolkit')
    addpath('/home/sbp29/MATLAB/sau-matlab-toolkit')
    addpath('/home/sbp29/MATLAB/sau-matlab-toolkit/grid-manipulation')
    addpath('/home/sbp29/MATLAB/sau-matlab-toolkit/other-sources')
    addpath('/home/sbp29/MATLAB/paris')
    addpath('/home/sbp29/MATLAB/development')

    javaaddpath('/home/sbp29/MATLAB/mysql-connector-java-8.0.12.jar');

%     Set preferences with setdbprefs.
    setdbprefs('DataReturnFormat', 'structure');
    setdbprefs({'NullStringRead';'NullStringWrite';'NullNumberRead';'NullNumberWrite'},...
                  {'null';'null';'NaN';'NaN'})

    expt_name = '4C3_GA';
    out_path = '/home/sbp29/MATLAB/4C3_Data/GA/S1Analysis/plateffect/';
    density = 6144;
    
%   MySQL Table Details  
    
    tablename_jpeg      = sprintf('%s_%d_JPEG',expt_name,density);
    
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
    
%%  GET PLATE DATA

    for ii = 1:length(hours)
        for iii = 1:length(n_plates.x6144plate)
            platedat{ii}{iii} = fetch(conn, sprintf(['select average from %s ',...
                'where hours = %d and pos in ',...
                '(select pos from %s ',...
                'where %s = %d ',...
                'order by %s, %s, %s)'],...
                tablename_jpeg,...
                hours(ii),...
                p2c_info(1,:),...
                p2c_info(2,:),...
                n_plates.x6144plate(iii),...
                p2c_info(2,:),...
                p2c_info(3,:),...
                p2c_info(4,:)));
       
%%  SOURCE-PLATE-WISE DIVISION
            
            [sourcedat{ii}{iii}{1},sourcedat{ii}{iii}{2},...
                sourcedat{ii}{iii}{3},sourcedat{ii}{iii}{4}] = ...
                downscale(col2grid(platedat{ii}{iii}.average));
            
            sourcedat{ii}{iii}{1}(sourcedat{ii}{iii}{1} == 0) = NaN;
            sourcedat{ii}{iii}{2}(sourcedat{ii}{iii}{2} == 0) = NaN;
            sourcedat{ii}{iii}{3}(sourcedat{ii}{iii}{3} == 0) = NaN;
            sourcedat{ii}{iii}{4}(sourcedat{ii}{iii}{4} == 0) = NaN;
            
%%  COMPARE SOURCE PLATES
            
            for t = 1:4
                for tt = 1:4
                    p{ii}{iii}(t,tt) = ranksum(grid2row(sourcedat{ii}{iii}{t}),grid2row(sourcedat{ii}{iii}{tt}));
                end
            end
           
%%  PLOT THE DISTRIBUTION
            
            fig = figure('Renderer', 'painters', 'Position', [10 10 1200 1000],'visible','off');
%             figure()
            subplot(2,2,1)
            bplot(grid2row(sourcedat{ii}{iii}{1}),hours(ii));
            xlim([hours(ii)-1,hours(ii)+1])
            ylim([0,600])
            title('Top Left Plate')
            xlabel('Hours')
            ylabel('Raw Pixel Count')

            subplot(2,2,2)
            bplot(grid2row(sourcedat{ii}{iii}{2}),hours(ii));
            xlim([hours(ii)-1,hours(ii)+1])
            ylim([0,600])
            title('Top Right Plate')
            xlabel('Hours')
            ylabel('Raw Pixel Count')

            subplot(2,2,3)
            bplot(grid2row(sourcedat{ii}{iii}{3}),hours(ii));
            xlim([hours(ii)-1,hours(ii)+1])
            ylim([0,600])
            title('Bottom Left Plate')
            xlabel('Hours')
            ylabel('Raw Pixel Count')

            subplot(2,2,4)
            bplot(grid2row(sourcedat{ii}{iii}{4}),hours(ii));
            xlim([hours(ii)-1,hours(ii)+1])
            ylim([0,600])
            title('Bottom Right Plate')
            xlabel('Hours')
            ylabel('Raw Pixel Count')
            
            %sgtitle(sprintf('%s\n%d hrs | Plate %d',expt_name,hours(ii),iii))
            saveas(fig,sprintf('%s%s_PE_%d_%d.png',...
                out_path,expt_name,hours(ii),iii))
            
        end
    end
    
%%  ENDs
     