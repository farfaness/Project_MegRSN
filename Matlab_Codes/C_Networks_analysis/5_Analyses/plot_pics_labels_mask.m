%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Creates mask to plot specifs ROI
%%%%%%%%% or agregated ROI
%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialisation
% load('/pclxserver/raid11/data/SELFI/archive_sourceplot_surface/pval_plot_int_spm8template.mat')
% load('/pclxserver/raid11/data/SELFI/scripts/Networks_analysis/5_Analyses/MyColormap.mat')

file_ROI = 'cluster' %or 'pics'

%% ROI to plots pics
%i = [14, 20, 18] % Prefrontal medial R 
%i = [25, 27, 24] % Occipital R : sup, mid, cuneus
%i = [37, 38, 39, 41, 42, 44] % Temporal R
%i = [2, 3, 5, 6, 7, 8, 9] % Frontal Lateral R

%% ROI to plots cluster
%i = [41, 43, 45, 37, 39, 35] % Occipital R: sup, mid, inf, cuneus, lingual, calcarine R
%i = [68, 70, 72, 74, 76, 78, 47, 32] %Temporal : Heschl, sup, pole sup, mid, pole mid, inf, fusiform, parahippocampal (limbique?)
%i = [68, 70, 72, 74, 76, 78] %Temporal : Heschl, sup, pole sup, mid, pole mid, inf R
%i = [68, 70, 72] %Temporal sup : Heschl, sup, pole sup R
%i = [74, 76, 78, 47] % Temporal mid et inf: mid, pole mid, inf, fusiform
%i = [68, 70, 72, 74, 76] %Temporal sup and mid: mid 74, pole mid 76
%i = [78, 47] %Temporal inf : inf, fusiform
%i = [4, 5, 7, 8, 10, 11, 12] %Frontal : sup, sup orb, mid, mid orb, inf oper, inf tri, inf orb R (+rolandic oper: 14?)
i = [19, 25, 27] %Front med R, cing ant, cing mid
%% --------------------for the points of the cluster -----------------------
if strcmp(file_ROI, 'cluster') %or pics
    load('/pclxserver/raid11/data/SELFI/archive_sourceplot_surface/cluster_label_025.mat')
    
    % --------------------for regrouping several ROI -----------------------
    if length(i)>1
        initialisation = cluster.(cell2mat(cluster.List_ROI(i(1)))).mask;
        mask_all_Roi(:, 1) = reshape(initialisation, [256*256*256 1]);
        
        for nbROI = 1:length(i)
            name_ROI{nbROI} = cell2mat(cluster.List_ROI(i(nbROI)));
            display(['***************************' name_ROI(nbROI) '*******************************'])
            
            mask_Roi{nbROI} = reshape(cluster.(name_ROI{nbROI}).mask, [256*256*256 1]);
            
            mask_all_Roi(:, 1) = (mask_all_Roi(:, 1)==1 | mask_Roi{nbROI}(:, 1)==1);
        end
        
        pval_plot_intp.mask(:, 1) = ((pval_plot_intp.pow(:, 1)>2.4450) & (mask_all_Roi(:, 1)==1));
        
        %--------------------for only one ROI------------------------------------
        elseif length(i)==1
        %Choose the Roi to plot
        name_ROI = cell2mat(cluster.List_ROI(i));
        display(['***************************' name_ROI '*******************************'])
        
        %Mask to plot
        mask_Roi = reshape(cluster.(name_ROI).mask, [256*256*256 1]);
        pval_plot_intp.mask(:, 1) = ((pval_plot_intp.pow(:, 1)>2.4450) & (mask_Roi(:, 1)==1));
        
    end

%%--------------------for the pics of the cluster -----------------------  
elseif strcmp(file_ROI, 'pics')
    load('/pclxserver/raid11/data/SELFI/archive_sourceplot_surface/pics_labels_005.mat')
    
    % --------------------for regrouping several ROI -----------------------
    if length(i)>1
        initialisation = pics.(cell2mat(pics.List_ROI(i(1)))).mask;
        mask_all_Roi(:, 1) = reshape(initialisation, [256*256*256 1]);
        
        for nbROI = 1:length(i)
            name_ROI{nbROI} = cell2mat(pics.List_ROI(i(nbROI)));
            display(['***************************' name_ROI(nbROI) '*******************************'])
            
            mask_Roi{nbROI} = reshape(pics.(name_ROI{nbROI}).mask, [256*256*256 1]);
            
            mask_all_Roi(:, 1) = (mask_all_Roi(:, 1)==1 | mask_Roi{nbROI}(:, 1)==1);
        end
        
        pval_plot_intp.mask(:, 1) = ((pval_plot_intp.pow(:, 1)>2.4450) & (mask_all_Roi(:, 1)==1));
        
        %--------------------for only one ROI------------------------------------
        elseif length(i)==1
        %Choose the Roi to plot
        name_ROI = cell2mat(pics.List_ROI(i));
        display(['***************************' name_ROI '*******************************'])
        
        %Mask to plot
        mask_Roi = reshape(pics.(name_ROI).mask, [256*256*256 1]);
        pval_plot_intp.mask(:, 1) = ((pval_plot_intp.pow(:, 1)>2.4450) & (mask_Roi(:, 1)==1));
        
    end
end
%% Surface plots

%Colormap to use
mymap = mymap_5hotColor;
mymap.color = mymap_red;
path2ft = '/lena13/home_users/users/hazem/fieldtrip-20161231'% 
%path2ft = 'C:\Users\Ness\Documents\Toolboxes\fieldtrip-20170116'

%----------------Surface Plot of the mask --------------------------

            % Left hemisphere internal side
            cfg              = [];
            cfg.method       = 'surface';
            cfg.funcolorlim  = 'zeromax' %'minzero'%'maxabs'
            cfg.funparameter = 'pow';
            cfg.camlight     = 'yes';
            cfg.colorbar = 'no'
            cfg.funcolormap   =  'hot';            
            cfg.maskparameter = 'mask'
            cfg.surffile     = [path2ft filesep 'template' filesep 'anatomy' filesep 'surface_white_left.mat']   
            ft_sourceplot(cfg, pval_plot_intp);
            view ([90 0]);
            set(gcf,'color','black')
            title(['Left hemisphere ' name_ROI] ,'Color', 'w')
            set(gca, 'CLim', [2 4.7805]);
            colormap(mymap.color);
            hc = colorbar;
%             set(hc, 'Color', [1 1 1], 'FontName', 'cambria', 'box', 'off','Location', 'southoutside', 'FontWeight', 'bold', 'Limits', [2.1  4.8], 'Ticks', [2.1, 2.9, 3.2, 3.9, 4.8], 'Ticklabels', {'p=0.05','p=0.01','p=0.005', 'p=0.001', 'p=0.00016'})   %color bar in white
            set(hc, 'Color', [1 1 1], 'FontName', 'cambria', 'box', 'off','Location', 'southoutside', 'FontWeight', 'bold', 'Limits', [2.4450  4.7805], 'Ticks', [mymap.Ticks(1, :), mymap.Ticks(2, :), mymap.Ticks(3, :), mymap.Ticks(4, :), mymap.Ticks(5, :), mymap.Ticks(6, :)], 'Ticklabels', {mymap.Index(1, :), mymap.Index(2, :), mymap.Index(3, :), mymap.Index(4, :), mymap.Index(5, :), mymap.Index(6, :)})   %color bar in white
                          
            % Left hemisphere external side
            cfg              = [];
            cfg.method       = 'surface';
            cfg.funcolorlim  = 'zeromax' %'minzero'%'maxabs'
            cfg.funparameter = 'pow';
            cfg.camlight     = 'yes';
            cfg.colorbar = 'no'
            cfg.funcolormap   =  'hot';              
            cfg.maskparameter = 'mask'
            cfg.surffile     = [path2ft filesep 'template' filesep 'anatomy' filesep 'surface_white_left.mat']
            ft_sourceplot(cfg, pval_plot_intp);
            view ([-90 0]);
            set(gcf,'color','black')
            title(['Left hemisphere ' name_ROI],'Color', 'w')
            set(gca, 'CLim', [2 4.7805]);
            colormap(mymap.color);
            hc = colorbar;
%             set(hc, 'Color', [1 1 1], 'FontName', 'cambria', 'box', 'off','Location', 'southoutside', 'FontWeight', 'bold', 'Limits', [2.1  4.8], 'Ticks', [2.1, 2.9, 3.2, 3.9, 4.8], 'Ticklabels', {'p=0.05','p=0.01','p=0.005', 'p=0.001', 'p=0.00016'})   %color bar in white
            set(hc, 'Color', [1 1 1], 'FontName', 'cambria', 'box', 'off','Location', 'southoutside', 'FontWeight', 'bold', 'Limits', [2.4450  4.7805], 'Ticks', [mymap.Ticks(1, :), mymap.Ticks(2, :), mymap.Ticks(3, :), mymap.Ticks(4, :), mymap.Ticks(5, :), mymap.Ticks(6, :)], 'Ticklabels', {mymap.Index(1, :), mymap.Index(2, :), mymap.Index(3, :), mymap.Index(4, :), mymap.Index(5, :), mymap.Index(6, :)})   %color bar in white
              
%-----------------------------------------------------------------------------------------------------------------------------            
            
            % Right hemisphere internal side
            cfg              = [];
            cfg.method       = 'surface';
            cfg.funcolorlim  = 'zeromax' %'minzero'%'maxabs'
            cfg.funparameter = 'pow';
            cfg.camlight     = 'yes';
            cfg.colorbar = 'no';
            cfg.funcolormap   =  'hot';             
            cfg.maskparameter = 'mask';
            cfg.surffile     = [path2ft filesep 'template' filesep 'anatomy' filesep 'surface_white_right.mat']
            ft_sourceplot(cfg, pval_plot_intp);
            view ([-90 0]);
            set(gcf,'color','black');
            title(['Right hemisphere ' name_ROI],'Color', 'w')
            set(gca, 'CLim', [2 4.7805]);
            colormap(mymap.color);
            hc = colorbar;
%             set(hc, 'Color', [1 1 1], 'FontName', 'cambria', 'box', 'off','Location', 'southoutside', 'FontWeight', 'bold', 'Limits', [2.1  4.8], 'Ticks', [2.1, 2.9, 3.2, 3.9, 4.8], 'Ticklabels', {'p=0.05','p=0.01','p=0.005', 'p=0.001', 'p=0.00016'})   %color bar in white
            set(hc, 'Color', [1 1 1], 'FontName', 'cambria', 'box', 'off','Location', 'southoutside', 'FontWeight', 'bold', 'Limits', [2.4450  4.7805], 'Ticks', [mymap.Ticks(1, :), mymap.Ticks(2, :), mymap.Ticks(3, :), mymap.Ticks(4, :), mymap.Ticks(5, :), mymap.Ticks(6, :)], 'Ticklabels', {mymap.Index(1, :), mymap.Index(2, :), mymap.Index(3, :), mymap.Index(4, :), mymap.Index(5, :), mymap.Index(6, :)})   %color bar in white             
            
            % Right hemisphere external side
            cfg              = [];
            cfg.method       = 'surface';
            cfg.funcolorlim  = 'zeromax' %'minzero'%'maxabs'
            cfg.funparameter = 'pow';
            cfg.camlight     = 'yes';
            cfg.colorbar = 'no';
            cfg.funcolormap   =  'hot';             
            cfg.maskparameter = 'mask';
            cfg.surffile     = [path2ft filesep 'template' filesep 'anatomy' filesep 'surface_white_right.mat']
            ft_sourceplot(cfg, pval_plot_intp);
            view ([90 0]);
            set(gcf,'color','black')
            title(['Right hemisphere ' name_ROI],'Color', 'w')
            set(gca, 'CLim', [2 4.7805]);
            colormap(mymap.color);
            hc = colorbar;
%             set(hc, 'Color', [1 1 1], 'FontName', 'cambria', 'box', 'off','Location', 'southoutside', 'FontWeight', 'bold', 'Limits', [2.1  4.8], 'Ticks', [2.1, 2.9, 3.2, 3.9, 4.8], 'Ticklabels', {'p=0.05','p=0.01','p=0.005', 'p=0.001', 'p=0.00016'})   %color bar in white
            set(hc, 'Color', [1 1 1], 'FontName', 'cambria', 'box', 'off','Location', 'southoutside', 'FontWeight', 'bold', 'Limits', [2.4450  4.7805], 'Ticks', [mymap.Ticks(1, :), mymap.Ticks(2, :), mymap.Ticks(3, :), mymap.Ticks(4, :), mymap.Ticks(5, :), mymap.Ticks(6, :)], 'Ticklabels', {mymap.Index(1, :), mymap.Index(2, :), mymap.Index(3, :), mymap.Index(4, :), mymap.Index(5, :), mymap.Index(6, :)})   %color bar in white
                       
 %-----------------------------------------------------------------------------------------------------------------------------            
 
            % Both hemisphere 
            cfg              = [];
            cfg.method       = 'surface';
            cfg.funcolorlim  = 'zeromax'; %'minzero'%'maxabs'
            cfg.funparameter = 'pow';
            cfg.camlight     = 'yes';
            cfg.funcolormap   =  'hot';     
            cfg.maskparameter = 'mask';
%             cfg.funcolorlim  = [0 5]
%             cfg.opacitylim  = [0 2]
%             cfg.opacitymap  = 'rampup'
            cfg.surffile     = [path2ft filesep 'template' filesep 'anatomy' filesep 'surface_white_both.mat']
            ft_sourceplot(cfg, pval_plot_intp);
            view ([0 90]);
            %view ([0 20])
            set(gcf,'color','black');
            title(['Both hemisphere' name_ROI] ,'Color', 'w');
            set(gca, 'CLim', [2 4.7805]);
            colormap(mymap.color);
            hc = colorbar;
%             set(hc, 'Color', [1 1 1], 'FontName', 'cambria', 'box', 'off','Location', 'southoutside', 'FontWeight', 'bold', 'Limits', [2.1  4.8], 'Ticks', [2.1, 2.9, 3.2, 3.9, 4.8], 'Ticklabels', {'p=0.05','p=0.01','p=0.005', 'p=0.001', 'p=0.00016'})   %color bar in white
            set(hc, 'Color', [1 1 1], 'FontName', 'cambria', 'box', 'off','Location', 'southoutside', 'FontWeight', 'bold', 'Limits', [2.4450  4.7805], 'Ticks', [mymap.Ticks(1, :), mymap.Ticks(2, :), mymap.Ticks(3, :), mymap.Ticks(4, :), mymap.Ticks(5, :), mymap.Ticks(6, :)], 'Ticklabels', {mymap.Index(1, :), mymap.Index(2, :), mymap.Index(3, :), mymap.Index(4, :), mymap.Index(5, :), mymap.Index(6, :)})   %color bar in white
 
            display(['***************************' name_ROI '*******************************'])         
            