% Copy right: Jia Guo, Assistant Professor in Neurology, Columbia University
% Contributor: Nanyan "Rosalie" Zhu, Columbia University


function [] = DRC_modeling(varargin)
    % Handle the optional input parameters.
    % =====================================================================
    p = inputParser;
    p.KeepUnmatched = true;
    % Specify the input datapath.
    addParamValue(p,'data_path','./Input/',@isstring);
    % Specify the 4d scan suffix.
    addParamValue(p,'scan_4d_suffix','4d.nii.gz',@isstring);
    % Specify the brain mask suffix.
    addParamValue(p,'BM_suffix','brainmask.nii.gz',@isstring);
    % Specify the background image suffix.
    addParamValue(p,'median_suffix','median.nii.gz',@isstring);
    % Specify the peudoinverse tolerance.
    addParamValue(p,'pinv_tolerance',8e-8, @isnumeric);
    
    % Parse the optional parameters.
    parse(p,varargin{:});

    % Handle the unexpected input parameter.
    UnmatchedParam = fieldnames(p.Unmatched);
    if ~isempty(UnmatchedParam)
        error(['"',UnmatchedParam{1},'" is not a valid parameter.']);
    end
    
    % glob all the subfolders
    data_path = p.Results.data_path;
    file_list = dir(data_path);
    file_flag = [file_list.isdir];
    folder_list = file_list(file_flag);
    folder_name_list = cell(1, length(folder_list) - 2);
    for i = 1:length(folder_list) - 2
        folder_name = folder_list(i + 2).name;
        folder_name_list{i} = folder_name;
    end
    for subject = 1: length(folder_name_list)
        close all
        clearvars -except p data_path folder_name_list subject
        % Tune parameters here.
        pinv_tolerance = p.Results.pinv_tolerance;
        coarse_cluster_K = 20;
        refine_cluster_K = 8;
        gaussian_sigma = [0.5 0.5 0.001];
        SplineSmoothingParam = 0.004;
        NM_discrete_datalength = 200;
        smoothing_param2 = 0.004;
        % selected CTC Cluster on line 298 or so

        % Ploting parameters
        start_slice = 2;
        plot_step = 1;
        selected_slice_idx = 6;

        % Setup Pathes 
        Input_Path = strcat(data_path, folder_name_list{subject}, '/');
        Output_Path = strcat('./Output/', folder_name_list{subject}, '/', 'pinv_tolerance_', num2str(pinv_tolerance), '_guassian_sigma05/');
        Output_Figure_Path = strcat(Output_Path, '/', 'Figures/');

        % Mask sure do not rerun the scans that have already been run
        if exist(Output_Path, 'dir')
            path_file_list = dir(Output_Path);
            if ~isempty(strfind(path_file_list(3).name, '.nii.gz'))
                ;
            else
                error(strcat(Output_Path, ' has not complete, please delete the folder.'))
            end

        else

            fprintf(strcat(num2str(subject), ' / ', num2str(length(folder_name_list)), ': ', folder_name_list{subject}, ...
                 '   Pinv Tolerance: ', num2str(pinv_tolerance)))

            %%
            % Import the matlab nifti toolbox
            addpath(genpath('./matlab_toolbox/'))

            % If no Output_Path, make a new folder
            if ~exist(Output_Path, 'dir')
                mkdir(Output_Path);
            end

            if ~exist(Output_Figure_Path, 'dir')
                mkdir(Output_Figure_Path);
            end

            backgroud_str = load_nifti(strcat(Input_Path,'*', p.Results.median_suffix));
            mask_str = load_nifti(strcat(Input_Path,'*', p.Results.BM_suffix));
            DSC4DRaw_str = load_nifti(strcat(Input_Path,'*', p.Results.scan_4d_suffix));
            %%
            mask = double(mask_str.vol);
            DSC_4D_img = double(DSC4DRaw_str.vol);


            scan_num = size(DSC_4D_img,4);


            NR_select = size(DSC_4D_img,4);

            for scan_i = 1:NR_select
                DSC3D = DSC_4D_img(:, :, :, scan_i);

                DSC_3D_smoothed = gauss3filter(DSC3D, gaussian_sigma);

                DSC_4D_orig_mat(:, :, :, scan_i) = DSC3D;
                DSC_4D_orig_all_array(:, scan_i) = DSC3D(:);
                DSC_4D_orig_masked_mat(:, :, :, scan_i) = DSC3D.*mask;
                DSC_4D_orig_masked_array(:, scan_i) = DSC3D(mask==1);

                DSC_4D_smoothed_mat(:, :, :, scan_i) = DSC_3D_smoothed;

            end


            %% DSC K-means
            start_ind = 5;
            TE = 45 * 10^(-3);
            TR = 3500 * 10^(-3);
            DeltaT = 1.46; % in mins

            pre = mean(DSC_4D_orig_all_array(:,1:start_ind),2);
            post = DSC_4D_orig_all_array(:,start_ind:end);

            % sum all the repititions 
            delta_R2_array_all = double(log(repmat(pre, 1, size(post, 2)) ./ post) ./ TE);
            delta_R2_array_all(isinf(delta_R2_array_all)) = 0;
            delta_R2_array_all(isnan(delta_R2_array_all)) = 0;
            delta_R2_array_AUC_all = sum(delta_R2_array_all(:, 10:end), 2);

            DSC_AUC_img = zeros(size(DSC3D));
            DSC_AUC_img(~isinf(mask)) = delta_R2_array_AUC_all;

            %%%%%%%%%%%%%%% Derive mask with positive AUC voxels in masked brain %%%%%%%%%%%%%%%% 
            DSC_AUC_masked_img=DSC_AUC_img .* mask;

            DSC_AUC_masked_img(DSC_AUC_masked_img < 0) = 0;
            DSC_AUC_masked_positive_mask = zeros(size(DSC_AUC_masked_img));
            DSC_AUC_masked_positive_mask(DSC_AUC_masked_img > 0) = 1;


            DSC_positive_mask = DSC_AUC_masked_positive_mask; % mask or DSC_AUC_masked_positive_mask (no CSF region)


            %%%%%%%%%%%%%%%%%%%%%%%%% whole brain cluestering %%%%%%%%%%%%%%%%%%%%

            figure();
            subplot(2,1,1)
            imagesc(imrotate(mask(:,:,9), 90))
            axis off;
            title('Brain Mask')
            subplot(2,1,2)
            imagesc(imrotate(DSC_positive_mask(:,:,9), 90))
            axis off;
            title('Positive Brain Mask')


            for scan_i=1:NR_select
                DSC3D = DSC_4D_img(:, :, :, scan_i);

                DSC_3D_smoothed = gauss3filter(DSC3D, gaussian_sigma);

                DSC_4D_smoothed_array(:, scan_i) = DSC_3D_smoothed(:);
                DSC_4D_smoothed_masked_brain_array(:, scan_i) = DSC_3D_smoothed(DSC_positive_mask == 1);
            end
            brain_ind = find(DSC_positive_mask == 1);
            nonebrain_ind = find(mask == 1 & DSC_positive_mask == 0);


            %% Perfusion Modeling for CBF, CBV, MTT and AUC, T_max, CTC_max 
            % ************** step 1: geneate the CTC **************
            % generate the CTC(as y) and Time(as x)
            DSC_4D_vol = DSC_4D_smoothed_mat;
            DSC_mean_vol = mean(DSC_4D_smoothed_mat, 4);
            DSC_pre_vol = DSC_4D_vol(:, :, :, start_ind);
            %%%%% pre scan: from ? to ? %%%%

            %%%%%%%%%%%%%%%% pre %%%%%%%%%%%%%%%%%%%%%%%%
            DSC_post_4D_vol = DSC_4D_vol(:, :, :, start_ind:end);

            %%%%%%%%%%%%%%%% post %%%%%%%%%%%%%%%%%%%%%%%
            DSC_pre_4D_vol = repmat(DSC_pre_vol, 1, 1, 1, size(DSC_post_4D_vol, 4));
            DSC_CTC_4D_vol = double(log(DSC_pre_4D_vol ./ DSC_post_4D_vol) ./ TE);
            Time = linspace(0, DeltaT * (size(DSC_CTC_4D_vol, 4)-1), size(DSC_CTC_4D_vol, 4));
            NM_discrete_Time = linspace(0, DeltaT * (size(DSC_CTC_4D_vol, 4) - 1), NM_discrete_datalength);

            %% ************** step 2: non-linear fitting of ICA CTC using smooth spline **************
            % Load the brain mask for generating the tissue CTC
            TissueMask = mask;
            TissueMask_4d = repmat(TissueMask, 1, 1, 1, size(DSC_CTC_4D_vol,4));
            DSC_CTC_Tissue_4D_vol = DSC_CTC_4D_vol .* TissueMask_4d;
            Size_4D_vol = size(DSC_CTC_Tissue_4D_vol);
            DSC_CTC_Tissue_4D_array = reshape(DSC_CTC_Tissue_4D_vol, Size_4D_vol(1) * Size_4D_vol(2) * Size_4D_vol(3), Size_4D_vol(4));
            DSC_CTC_Tissue_4D_masked_array = DSC_CTC_Tissue_4D_array(TissueMask == 1, :);

            % Option 1: automatic ICA finding method
            % coarse ICA finding
            DSC_cluster = kmeans(DSC_4D_smoothed_masked_brain_array, coarse_cluster_K);
            DSC_cluster_img = zeros(size(mask));
            DSC_cluster_img(DSC_positive_mask == 1) = DSC_cluster;

            for i = 1:coarse_cluster_K
                delta_R2_cluster_timecourse(i, :) = mean(delta_R2_array_all(DSC_cluster_img == i, :));
                delta_R2_cluster_timecourse_AUC(i) = sum(mean(delta_R2_array_all(DSC_cluster_img == i, :)));
            end
            [~, Cluster_Order] = sort(delta_R2_cluster_timecourse_AUC, 'descend');

            ICAMask = zeros(size(mask));
            ICAMask(DSC_cluster_img == Cluster_Order(1)) = 1;

            ICAMask_4d=repmat(ICAMask, 1, 1, 1, size(DSC_CTC_4D_vol, 4));
            DSC_CTC_ICA_4D_vol = DSC_CTC_4D_vol .* ICAMask_4d;
            DSC_CTC_ICA_4D_array = reshape(DSC_CTC_ICA_4D_vol, Size_4D_vol(1) * Size_4D_vol(2) * Size_4D_vol(3), Size_4D_vol(4));
            DSC_CTC_ICA_4D_masked_array = DSC_CTC_ICA_4D_array(ICAMask == 1, :);

            % ICA K-means to find the real ICA
            ICA_cluster = kmeans(DSC_CTC_ICA_4D_masked_array, refine_cluster_K);
            ICA_cluster_img = zeros(size(mask));
            ICA_cluster_img(ICAMask == 1) = ICA_cluster;

            for i = 1:refine_cluster_K
                DSC_CTC_ICA_4D_Kmeans_array(i, :) = mean(DSC_CTC_ICA_4D_array(ICA_cluster_img == i, :));
                DSC_CTC_ICA_4D_Kmeans_array_AUC(i) = sum(DSC_CTC_ICA_4D_Kmeans_array(i, 1:end));
            end
            % smooth spline fit of the real ICA CTC

            figure()
            [~, ICA_plot_order] = sort(DSC_CTC_ICA_4D_Kmeans_array_AUC, 'descend');
            for i = 1:refine_cluster_K
                hold on
                plot(Time, DSC_CTC_ICA_4D_Kmeans_array(ICA_plot_order(i), :), '.');
                xlabel('min');
                ylabel('CTC ICA');
            end
            legend(strcat('ICA Cluster=', num2str(ICA_plot_order')), 'FontSize', 8)
            CTC_ICA = DSC_CTC_ICA_4D_Kmeans_array(ICA_plot_order(1), :);
            CTC_ICA_fit = fit(Time', CTC_ICA', 'smoothingspline', 'SmoothingParam', SplineSmoothingParam);
            hold on;
            plot(CTC_ICA_fit, Time, CTC_ICA);
            xlim([0 max(Time)]);
            title('Smooth Spline Fit of ICA CTC');
            legend('Location', 'Southeast')
            saveas(gcf, strcat(Output_Figure_Path, 'CTC_fitting.png'))


            for i = 1:size(DSC_CTC_ICA_4D_masked_array, 1)
                DSC_CTC_ICA_4D_temp = DSC_CTC_ICA_4D_masked_array(i, :);
                CTC_ICA_fit_temp = fit(Time', DSC_CTC_ICA_4D_temp', 'smoothingspline', 'SmoothingParam', smoothing_param2);
                NM_discrete_CTC_ICA(i, :) = CTC_ICA_fit_temp(NM_discrete_Time);
            end

            % Option 2: load the manually modified ICA mask for generating the AIF
            % ICAMask_str=load_nii('.\Output\DSC\Test00\DSC_cluster_ICA_label.nii.gz');
            % ICAMask=double(ICAMask_str.vol);


            %% **** step 3: generate the AIF matrix, AIF_mat *******
            NM_discrete_CTC_ICA_mean = CTC_ICA_fit(NM_discrete_Time);
            AIF = NM_discrete_CTC_ICA_mean;
            AIF = AIF-min(AIF);
            AIF_mat = zeros(length(NM_discrete_Time), length(NM_discrete_Time));

            tic
            DeltaT = NM_discrete_Time(2)-NM_discrete_Time(1);
            for i = 1:length(NM_discrete_Time)
                for j = 1:length(NM_discrete_Time)
                    if i == j
                        AIF_mat(i, j) = (DeltaT/5) * (4 * AIF(1) + AIF(2));
                    elseif i > j && (i - j) < length(NM_discrete_Time) - 1
                        AIF_mat(i, j) = (DeltaT/6) * (AIF(i-j) + 4 * AIF(i-j+1) + AIF(i-j+2));
                    elseif (i-j) == length(NM_discrete_Time) - 1
                        AIF_mat(i, j) = (DeltaT/5) * (AIF(i-j) + 4 * AIF(i-j+1));
                    else
                        AIF_mat(i, j) = 0;
                    end
                end
            end
            toc

            %% ************** step 4: non-linear fitting of Tissue CTC using smooth spline and calculate the CBF using decovolution **************
            % Initialize Parameters
            NM_discrete_CTC_Tissue = zeros(size(DSC_CTC_Tissue_4D_masked_array, 1), length(NM_discrete_Time));
            CBF_mul_r = zeros(size(DSC_CTC_Tissue_4D_masked_array, 1), length(NM_discrete_Time));

            % fit all orig CTC in ICA mask
            tic 
            parfor i = 1:size(DSC_CTC_Tissue_4D_masked_array, 1)
                DSC_CTC_Tissue_4D_temp = DSC_CTC_Tissue_4D_masked_array(i,:);
                CTC_Tissue_fit_temp = fit(Time', DSC_CTC_Tissue_4D_temp', 'smoothingspline', 'SmoothingParam', SplineSmoothingParam);
                NM_discrete_CTC_Tissue_temp = CTC_Tissue_fit_temp(NM_discrete_Time);
                NM_discrete_CTC_Tissue(i, :) = NM_discrete_CTC_Tissue_temp;

                % generate CBF using deconvolution: 
                CBF_mul_r(i,:) = pinv(AIF_mat, pinv_tolerance) * (NM_discrete_CTC_Tissue_temp - min(NM_discrete_CTC_Tissue_temp)); 
            end
            toc

            % CBF = max[CBF * R] 
            CBF_Tissue_masked_array = max((CBF_mul_r), [], 2);
            CBF_vol = zeros(size(mask));
            CBF_vol(TissueMask == 1) = CBF_Tissue_masked_array;

            %% Save figure for paper 
            close all

            % selected CTC cluster
            Selected_Cluster = Cluster_Order(end-2);

            % AIF generation
            AIF_timecourse_selected = DSC_CTC_ICA_4D_Kmeans_array(ICA_plot_order(1),:);
            delta_R2_cluster_timecourse_selected = delta_R2_cluster_timecourse(Selected_Cluster,:);
            delta_R2_cluster_timecourse_selected_fit = fit(Time',delta_R2_cluster_timecourse_selected(1:length(Time))',...
                'smoothingspline', 'SmoothingParam', SplineSmoothingParam);
            AIF_timecourse_selected_fit = fit(Time', AIF_timecourse_selected', 'smoothingspline', 'SmoothingParam', SplineSmoothingParam);

            AIF_timecourse_selected_sp = AIF_timecourse_selected_fit(NM_discrete_Time);
            AIF_timecourse_selected_sp_min = min(AIF_timecourse_selected_sp(:));
            AIF_timecourse_selected_sp = AIF_timecourse_selected_sp-AIF_timecourse_selected_sp_min;
            delta_R2_cluster_timecourse_selected_sp = delta_R2_cluster_timecourse_selected_fit(NM_discrete_Time);
            delta_R2_cluster_timecourse_selected_sp_min = min(delta_R2_cluster_timecourse_selected_sp(:));
            delta_R2_cluster_timecourse_selected_sp = delta_R2_cluster_timecourse_selected_sp - delta_R2_cluster_timecourse_selected_sp_min;

            AIF = AIF_timecourse_selected_sp;
            AIF_mat = zeros(length(NM_discrete_Time), length(NM_discrete_Time));

            tic
            NM_DeltaT = NM_discrete_Time(2) - NM_discrete_Time(1);
            for i = 1:length(NM_discrete_Time)
                for j = 1:length(NM_discrete_Time)
                    if i == j
                        AIF_mat(i, j) = (NM_DeltaT/6) * (0 + 4 * AIF(1) + AIF(2));
                    elseif i > j && (i-j) < length(NM_discrete_Time) - 1
                        AIF_mat(i, j) = (NM_DeltaT/6) * (AIF(i-j) + 4*  AIF(i-j+1) + AIF(i-j+2));
                    elseif (i-j) == length(NM_discrete_Time) - 1
                        AIF_mat(i, j) = (NM_DeltaT/6) * (AIF(i-j) + 5 * AIF(i-j+1));
                    else
                        AIF_mat(i, j) = 0;
                    end
                end
            end
            toc

            figure(); 
            imagesc(AIF_mat); 
            axis image; 
            colormap('jet'); colorbar; 
            title('The AIF Matrix')
            saveas(gcf, strcat(Output_Figure_Path, 'The_AIF_Matrix.png'))

            figure(); 
            plot(AIF_mat(end,:)); 
            title('Last Row of the AIF Matirx')
            saveas(gcf, strcat(Output_Figure_Path, 'Last_Row_AIF_Matrix.png'))


            % Tissue_Residual is CBF * R (Tissue Residual Function)
            Tissue_Residual = pinv(AIF_mat, pinv_tolerance) * delta_R2_cluster_timecourse_selected_sp;


            figure('Renderer', 'painters', 'Position', [10 10 1000 800]);
            subplot(2,3,1)
            plot(Time, AIF_timecourse_selected - AIF_timecourse_selected_sp_min, '.');
            hold on;
            plot(NM_discrete_Time, AIF_timecourse_selected_sp, 'r');
            xlabel('time (min)'); ylabel('AIF (a.u.)');
            xlim([0 max(NM_discrete_Time)]); ylim([0 1.1 * max(AIF_timecourse_selected_sp)])
            legend({'Raw', 'S-Spline Fit'} ,'FontSize', 10, 'location', 'Southeast')

            subplot(2,3,2)
            plot(Time, delta_R2_cluster_timecourse_selected - delta_R2_cluster_timecourse_selected_sp_min,'.');
            hold on;
            plot(NM_discrete_Time, delta_R2_cluster_timecourse_selected_sp,'b');
            xlabel('time (min)'); ylabel('CTC (a.u.)');
            xlim([0 max(NM_discrete_Time)]); ylim([0 1.1 * max(delta_R2_cluster_timecourse_selected_sp)])
            legend({'Raw','S-Spline Fit'}, 'FontSize', 10, 'location', 'Southeast')

            subplot(2,3,3)
            plot(NM_discrete_Time, Tissue_Residual,'M');
            xlabel('time (min)');ylabel('CBF x R (a.u.)');
            xlim([0 max(NM_discrete_Time)]);ylim([-0.2 * max(Tissue_Residual) 1.1 * max(Tissue_Residual)])

            subplot(2,3,4:6)
            plot(NM_discrete_Time, AIF_timecourse_selected_sp, 'r');
            xlabel('time (min)');xlim([0 max(NM_discrete_Time)]);
            ylim([0 1.1 * max(AIF_timecourse_selected_sp)])
            hold on
            plot(NM_discrete_Time, delta_R2_cluster_timecourse_selected_sp, 'b');
            set(gcf,'color','white')
            hold on;
            delta_R2_cluster_timecourse_estimated = conv(Tissue_Residual, NM_DeltaT .* AIF_timecourse_selected_sp);
            plot(NM_discrete_Time, delta_R2_cluster_timecourse_estimated(1:length(NM_discrete_Time)), 'go');
            legend({'AIF', 'CTC', 'CTC Derived'}, 'FontSize', 10, 'location', 'Northwest')

            saveas(gcf, strcat(Output_Figure_Path, 'Modeling_perfusion_haemodynamics.png'))
            % Tissue_Residual=ifft(fft(delta_R2_cluster_timecourse_selected_2000)./fft(AIF_timecourse_selected_2000));

            % Derived Tissue Residual, estimated AIF and estimated CTC normalization
            Tissue_Residual_N = Tissue_Residual ./ max(abs(Tissue_Residual(2:end)));
            AIF_timecourse_selected_sp_N = AIF_timecourse_selected_sp ./ max(AIF_timecourse_selected_sp);
            delta_R2_cluster_timecourse_selected_sp_N = delta_R2_cluster_timecourse_selected_sp ./ max(delta_R2_cluster_timecourse_selected_sp);

            figure();
            plot(AIF_timecourse_selected_sp_N, 'linew', 2)
            hold on
            plot(delta_R2_cluster_timecourse_selected_sp_N, 'linew', 2);
            hold on;
            plot(Tissue_Residual_N(1:length(NM_discrete_Time)), 'linew', 2);
            legend('AIF', 'CTC', 'CBF * R')
            title('Normalized Derived Tissue Residual, estimated AIF and CTC')
            saveas(gcf, strcat(Output_Figure_Path, 'Normalized_haemodynamics.png'))


            figure()
            plot(delta_R2_cluster_timecourse_estimated(1:length(NM_discrete_Time)), 'or');
            hold on;
            plot(delta_R2_cluster_timecourse_selected_sp, 'g', 'linew', 2);
            legend('CTC estimated', 'CTC')

            %% Generate CBF, CBV, MTT, TPP, WIR

            % CBV Generation
            CBF_Tissue_masked_array = max((CBF_mul_r),[],2);
            CBF_vol = zeros(size(mask));
            CBF_vol(TissueMask == 1) = CBF_Tissue_masked_array;

            rho_tissue = 1;% g/ml
            k_H = 0.733;

            CBF_vol = 100 .* CBF_vol .* k_H ./ NM_DeltaT; % ml/min/100g
            CBF_vol(DSC_positive_mask ~= 1) = 0;

            figure()
            suptitle('Cerebral Blood Flow (CBF)')
            for i=1:9
                subplot(3,3,i)
                slice_idx = start_slice + plot_step .* (i-1);
                imagesc(imrotate(squeeze(CBF_vol(:, :, slice_idx)), -90)); colorbar
                title(strcat('Slice ', num2str(slice_idx)))
                axis image;axis off;colormap('jet');
                caxis([0 250]);
            end
            set(gcf, 'color', 'white')
            saveas(gcf, strcat(Output_Figure_Path, 'sclices_Cerebral_Blood_Flow.png'))

            %%
            % ************** step 5: calculate the CBV using the ratio between Integral(CTC) and Integral(AIF) **************
            % CTC_Tissue_Integral=sum(NM_discrete_CTC_Tissue,2);
            % AIF_Tissue_Integral=sum(NM_discrete_CTC_ICA_mean);
            CTC_Tissue_Max = max(NM_discrete_CTC_Tissue, [], 2);
            %CTC_Tissue_Max = NM_discrete_CTC_Tissue(:, round(size(NM_discrete_CTC_Tissue, 2) .* 37.5 / 64));
            AIF_Tissue_Max = max(max(NM_discrete_CTC_ICA, [], 2));

            % generate the Blood Volume Fraction (BVf)
            BVf_Tissue_array = CTC_Tissue_Max ./ AIF_Tissue_Max;

            % save the CBV as 3D volume
            CBV_vol = zeros(size(mask));
            CBV_vol(TissueMask == 1) = 50 .* BVf_Tissue_array .* k_H ./ rho_tissue;% ml/100g
            % why 50 * CBV ???????????????

            CBV_vol(DSC_positive_mask ~= 1) = 0;
            % CBV_vol=medfilt3(CBV_vol,[3 3 1]);
            % view the CBV

            figure()
            suptitle('Cerebral Blood Volume (CBV)_(ml/100g)')
            for i = 1:9
                subplot(3,3,i)
                slice_idx = start_slice + plot_step .* (i-1);
                imagesc(imrotate(squeeze(CBV_vol(:, :, slice_idx)), -90));
                title(strcat('Slice ', num2str(slice_idx)))
                axis image; axis off;
                colormap('jet'); caxis([0 10]); colorbar
            end
            set(gcf, 'color', 'white')
            saveas(gcf, strcat(Output_Figure_Path, 'sclices_Cerebral_Blood_Volume.png'))

            % ************** step 6: derive the MTT using MTT=CBV/CBF **************
            MTT_vol = 60 .* CBV_vol ./ CBF_vol;% in sec
            MTT_vol(isinf(MTT_vol)) = 0;
            MTT_vol(isnan(MTT_vol)) = 0;
            MTT_vol(mask ~= 1) = 0;
            MTT_vol_Th = 10;
            MTT_vol(MTT_vol > MTT_vol_Th) = MTT_vol_Th;


            figure()
            suptitle('Mean Transit Time (MTT)')
            for i = 1:9
                subplot(3,3,i)
                slice_idx = start_slice + plot_step .* (i-1);
                imagesc(imrotate(squeeze(MTT_vol(:, :, slice_idx)), -90));
                title(strcat('Slice ', num2str(slice_idx)))
                axis image; axis off;
                colormap('jet'); colorbar
                caxis([0 9])
            end
            set(gcf,'color','white')
            saveas(gcf, strcat(Output_Figure_Path, 'sclices_Mean_Transit_Time.png'))

            % ************** step 7: derive the TTP, CAP, WIR=d(CTC)/d(t)|t in [0, TTP] for brain tissue **************
            [CAP_array, TTP_idx_array] = max(NM_discrete_CTC_Tissue, [], 2);
            TTP_array = NM_discrete_Time(1, TTP_idx_array);

            % Initialize WIR_array
            WIR_array = [];

            % averageg wash-in rate
            for i = 1:length(NM_discrete_CTC_Tissue)
                NM_discrete_CTC_Tissue_temp = NM_discrete_CTC_Tissue(i, :);
                TTP_idx_temp=TTP_idx_array(i);
                if TTP_idx_temp > 100
                   TTP_idx_temp = 100;
                end
                WIR_array(i) = mean(gradient(NM_discrete_CTC_Tissue_temp(1:TTP_idx_temp)));
            end

            % concentration at peak; time to peak
            CAP_vol = zeros(size(mask));
            TTP_vol = zeros(size(mask));
            WIR_vol = zeros(size(mask));
            CAP_vol(TissueMask == 1) = CAP_array;
            TTP_vol(TissueMask == 1) = TTP_array;
            WIR_vol(TissueMask == 1) = WIR_array;

            figure()
            suptitle('Concentration At Peak (CAP)')
            for i = 1:9
                subplot(3,3,i)
                slice_idx = start_slice + plot_step .* (i-1);
                imagesc(imrotate(squeeze(CAP_vol(:, :, slice_idx)), -90));
                title(strcat('Slice ', num2str(slice_idx)))
                axis image; axis off;
                colormap('jet'); colorbar
                caxis([0 max(CAP_vol(:))/2]);
            end
            set(gcf,'color','white')
            saveas(gcf, strcat(Output_Figure_Path, 'sclices_Concentration_At_Peak.png'))


            figure()
            suptitle('Time To Peak (TTP)')
            for i = 1:9
                subplot(3,3,i)
                slice_idx = start_slice + plot_step .* (i-1);
                imagesc(imrotate(squeeze(TTP_vol(:, :, slice_idx)), -90));
                title(strcat('Slice ', num2str(slice_idx)))
                axis image; axis off;
                colormap('jet'); 
                caxis([0 64])
            end
            set(gcf, 'color', 'white')
            saveas(gcf, strcat(Output_Figure_Path, 'sclices_Time_To_Peak.png'))

            figure()
            suptitle('Wash-in Rate (WIR)')
            for i = 1:9
                subplot(3,3,i)
                slice_idx = start_slice + plot_step .* (i-1);
                imagesc(imrotate(squeeze(WIR_vol(:, :, slice_idx)), -90));
                title(strcat('Slice ', num2str(slice_idx)))
                axis image; axis off;
                colormap('jet'); 
                caxis([0 0.5])
            end
            set(gcf, 'color', 'white')
            saveas(gcf, strcat(Output_Figure_Path, 'sclices_Washin_Rate.png'))

            % GEMRI
            GEMRI_vol = 2000 - 1 .* DSC_mean_vol;
            GEMRI_vol = GEMRI_vol .* mask;

            figure()
            suptitle('GEMRI')
            for i = 1:9
                subplot(3,3,i)
                slice_idx = start_slice + plot_step .* (i-1);
                imagesc(imrotate(squeeze(GEMRI_vol(:, :, slice_idx)), -90));
                title(strcat('Slice ', num2str(slice_idx)))
                axis image;axis off; 
                colormap('jet'); 
                caxis([0 max(GEMRI_vol(:).*0.8)])
            end
            set(gcf,'color','white')
            saveas(gcf, strcat(Output_Figure_Path, 'sclices_MRI.png'))

            figure();

            suptitle('Histograms in DSC positive brain region')

            subplot(6,1,1)
            histogram(CBF_vol(CBF_vol ~= 0), 255)
            title('CBF')

            subplot(6,1,2)
            histogram(CBV_vol(DSC_positive_mask == 1))
            title('CBV');

            subplot(6,1,3);
            histogram(TTP_vol(mask==1&TTP_vol~=64))
            title('TTP')

            subplot(6,1,4)
            histogram(WIR_vol(mask==1&WIR_vol>0))
            title('WIR')

            subplot(6,1,5)
            histogram(GEMRI_vol(mask == 1 & GEMRI_vol > 0));
            xlim([-500 max(GEMRI_vol(mask == 1 & GEMRI_vol ~= 0))])
            title('MRI')

            subplot(6,1,6)
            histogram(MTT_vol(mask == 1));
            title('MTT');
            xlim([0 10])
            saveas(gcf, strcat(Output_Figure_Path, 'sclices_Histograms_Positive.png'))

            %% plot selected slice
            close all
            figure()
            suptitle('Cerebral Blood Flow (CBF)')
            imagesc(imrotate(squeeze(CBF_vol(:, :, selected_slice_idx)), -90));
            axis image; axis off;
            colormap('jet'); colorbar
            caxis([0 2500]);
            set(gcf,'color','white')
            saveas(gcf, strcat(Output_Figure_Path, 'Slice', num2str(selected_slice_idx), '_Cerebral_Blood_Flow.png'))

            figure()
            suptitle('Cerebral Blood Volume (CBV)')
            imagesc(imrotate(squeeze(CBV_vol(:, :, selected_slice_idx)), -90));
            axis image; axis off;
            colormap('jet'); colorbar
            caxis([0 30]);
            set(gcf, 'color', 'white')
            saveas(gcf, strcat(Output_Figure_Path, 'Slice', num2str(selected_slice_idx), '_Cerebral_Blood_Volume.png'))

            figure()
            suptitle('Mean Transit Time (MTT)')
            imagesc(imrotate(squeeze(MTT_vol(:, :, selected_slice_idx)), -90));
            axis image; axis off;
            colormap('jet'); colorbar
            caxis([0 1]);
            set(gcf, 'color', 'white')
            saveas(gcf, strcat(Output_Figure_Path, 'Slice', num2str(selected_slice_idx), '_Mean_Transit_Time.png'))

            figure()
            suptitle('Time To Peak (TTP)')
            imagesc(imrotate(squeeze(TTP_vol(:,:,selected_slice_idx)),-90));
            axis image; axis off;
            colormap('jet'); colorbar
            caxis([0 35]);
            set(gcf, 'color', 'white')
            saveas(gcf, strcat(Output_Figure_Path, 'Slice', num2str(selected_slice_idx), '_Time_To_Peak.png'))

            figure()
            suptitle('Wash-in Rate (WIR)')
            imagesc(imrotate(squeeze(WIR_vol(:, :, selected_slice_idx)), -90));
            axis image; axis off;
            colormap('jet'); colorbar
            caxis([0 0.5]);
            set(gcf, 'color', 'white')
            saveas(gcf, strcat(Output_Figure_Path, 'Slice', num2str(selected_slice_idx), '_sclices_Washin_Rate.png'))
            %% save nii

            backgroud_str.vol = DSC_cluster_img;
            save_nifti(backgroud_str,strcat(Output_Path,'ADRC_cluster_x',num2str(refine_cluster_K),'_label.nii.gz'));

            backgroud_str.vol=GEMRI_vol;
            save_nifti(backgroud_str,strcat(Output_Path,'ADRC_GEMRI.nii.gz'));

            backgroud_str.vol=CBF_vol;
            save_nifti(backgroud_str,strcat(Output_Path,'ADRC_CBF.nii.gz'));

            backgroud_str.vol=CBV_vol;
            save_nifti(backgroud_str,strcat(Output_Path,'ADRC_CBV.nii.gz'));

            backgroud_str.vol=MTT_vol;
            save_nifti(backgroud_str,strcat(Output_Path,'ADRC_MTT.nii.gz'));

            backgroud_str.vol=TTP_vol;
            save_nifti(backgroud_str,strcat(Output_Path,'ADRC_TTP.nii.gz'));

            backgroud_str.vol=WIR_vol;
            save_nifti(backgroud_str,strcat(Output_Path,'ADRC_WIR.nii.gz'));

            backgroud_str.vol=WIR_vol;
            save_nifti(backgroud_str,strcat(Output_Path,'ADRC_WIR.nii.gz'));
        end
    end
end
    
