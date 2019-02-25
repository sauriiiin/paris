%%  Sau MATLAB Colony Analyzer Toolkit
%
%%  LinearInNorm6144.m

%   Author: Saurin Parikh, February 2019
%   dr.saurin.parikh@gmail.com

%   Linear Interpolation to predict plate background using control colony
%   pixel counts.
%   Without interleaving.

%%

    function data_fit = LinearInNorm6144(hours,IL,n_plates,p2c_info,cont_name,...
        tablename_p2o,tablename_jpeg)
    
        connectSQL;
        
        for ii = 1:length(hours)
            temp = [];
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
                    cont_name,...
                    p2c_info(3,:),...
                    p2c_info(4,:)));

                avg_data = fetch(conn, sprintf(['select a.pos, a.hours, a.average ',...
                    'from %s a, %s b ',...
                    'where a.pos = b.pos and hours = %d and %s = %d ',...
                    'order by %s, %s'],...
                    tablename_jpeg,...
                    p2c_info(1,:),...
                    hours(ii),...
                    p2c_info(2,:),...
                    n_plates.x6144plate(iii),...
                    p2c_info(3,:),...
                    p2c_info(4,:)));

        %%  CALCULATE BACKGROUND

%                 all_avg  = col2grid(avg_data.average);
%                 nonzero  = all_avg > 0;
                cont_pos = col2grid(ismember(pos.all.pos, pos.cont.pos));
                cont_avg = col2grid(avg_data.average).*cont_pos;

                cont_avg(cont_avg == 0) = NaN;
                
                if IL == 0
                    [a,b,c,d] = downscale(cont_avg);
                    plates = {a,b,c,d};

                    for i=1:4
                        [p,q,r,s] = downscale(plates{i});
                        [xq,yq] = ndgrid(1:32,1:48);

                        if nansum(nansum(p)) ~= 0 %Top Left
                            P = contBG(p);
                            p = (fillmissing(fillmissing(p, 'linear',2),'linear',1) +...
                                (fillmissing(fillmissing(p, 'linear',1),'linear',2)))/2;
                            [x,y] = ndgrid(1:2:32,1:2:48);
                            f = griddedInterpolant(x,y,p,'linear');
                            plates{i} = f(xq,yq);
                            [~,x,y,z] = downscale(plates{i});
                            bground{i} = plategen(P,x,y,z);

                        elseif nansum(nansum(q)) ~= 0 % Top Right
                            Q = contBG(q);
                            q = (fillmissing(fillmissing(q, 'linear',2),'linear',1) +...
                                (fillmissing(fillmissing(q, 'linear',1),'linear',2)))/2;
                            [x,y] = ndgrid(1:2:32,2:2:48); 
                            f = griddedInterpolant(x,y,q,'linear');
                            plates{i} = f(xq,yq);
                            [x,~,y,z] = downscale(plates{i});
                            bground{i} = plategen(x,Q,y,z);

                        elseif nansum(nansum(r)) ~= 0 % Bottom Left
                            R = contBG(r);
                            r = (fillmissing(fillmissing(r, 'linear',2),'linear',1) +...
                                (fillmissing(fillmissing(r, 'linear',1),'linear',2)))/2;
                            [x,y] = ndgrid(2:2:32,1:2:48); 
                            f = griddedInterpolant(x,y,r,'linear');
                            plates{i} = f(xq,yq);
                            [x,y,~,z] = downscale(plates{i});
                            bground{i} = plategen(x,y,R,z);

                        else % Bottom Right
                            S = contBG(s);
                            s = (fillmissing(fillmissing(s, 'linear',2),'linear',1) +...
                                (fillmissing(fillmissing(s, 'linear',1),'linear',2)))/2;
                            [x,y] = ndgrid(2:2:32,2:2:48); 
                            f = griddedInterpolant(x,y,s,'linear');
                            plates{i} = f(xq,yq);
                            [x,y,z,~] = downscale(plates{i});
                            bground{i} = plategen(x,y,z,S);

                        end
                    end
                    bg{iii} = grid2row(plategen(bground{1},bground{2},bground{3},bground{4}))';%.*nonzero)';
                    bg{iii}(bg{iii} == 0) = NaN;
                    temp = abs([temp; [pos.all.pos, ones(length(pos.all.pos),1)*hours(ii),...
                        bg{iii}, avg_data.average, avg_data.average./bg{iii}]]);
                else
                    print('working on it')
                end
            end
            data_fit{ii} = temp;
        end
    end
    
%%  END    
    
