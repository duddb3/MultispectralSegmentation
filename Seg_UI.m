function Seg_UI(varargin)

    persistent RIM
    if nargin == 0
        selector = 0;
        S  = get(0,'ScreenSize');
    else
        selector = cell2mat(varargin);
    end

    % clr1 = [1 170 199]./255;
    clr1 = [0 0 0];
    clr2 = [132 198 66]./255;
    clr3 = [1 85 99]./255;
    switch selector
        case 0 % Initialize variables and objects
            RIM.location = pwd;
            RIM.F = figure('IntegerHandle', 'off', ...
                'Name',sprintf('Seg_UI'), ...
                'NumberTitle','off',...
                'Tag','Froot',...
                'Units','Normalized',...
                'Position',[.1 .1 .8 .76],...
                'Pointer','arrow',...
                'Color',clr1,... 
                'Resize','on',...
                'MenuBar','none',...
                'DockControls','off',...
                'HandleVisibility','On',...
                'WindowScrollWheelFcn',@scroll,...
                'WindowButtonDownFcn',@AddPoint,...
                'WindowKeyPressFcn',@keyitup,...
                'WindowButtonMotionFcn',@mouseMove,...
                'CloseRequestFcn','Seg_UI(''exit'')',...
                'Visible','On');
            %==============================================================
            % Menu Bars
            %==============================================================
            p1 = uimenu('Label', 'File',...
                'Parent', RIM.F);
            uimenu(p1,'Label','Start New Case', ...
                'Callback', 'Seg_UI(''New'')');
            uimenu(p1,'Label','Save Analysis', ...
                'Callback', 'Seg_UI(''Save'')');
            uimenu(p1,'Label','Load Analyis',...
                'Callback', 'Seg_UI(''Load'')');
            uimenu(p1,'Label','Exit', ...
                'Callback', 'Seg_UI(''exit'')') ;

            %==============================================================
            % Display Axes
            %==============================================================
            RIM.v1 = axes(...
                'HandleVisibility', 'on',...
                'Visible', 'on',...
                'Tag','ax1',...
                'Position',[0.26 0.05 .48 0.9],...
                'NextPlot', 'replacechildren',...
                'ButtonDownFcn',@AddPoint,...
                'Parent', RIM.F);
            RIM.v1.XAxis.Visible = 'off';
            RIM.v1.YAxis.Visible = 'off';

            RIM.v2 = axes(...
                'HandleVisibility','on',...
                'Visible','on',...
                'Tag','ax2',...
                'Position',[0.75 0.525 .24 0.425],...
                'NextPlot', 'replacechildren',...
                'Parent', RIM.F);
            RIM.v2.XAxis.Visible = 'off';
            RIM.v2.YAxis.Visible = 'off';

            RIM.v3 = axes(...
                'HandleVisibility','on',...
                'Visible','on',...
                'Tag','ax3',...
                'Position',[0.75 0.05 .24 0.425],...
                'NextPlot', 'replacechildren',...
                'Parent', RIM.F);
            RIM.v3.XAxis.Visible = 'off';
            RIM.v3.YAxis.Visible = 'off';


            uicontrol("Style",'text',...
                'String','Image Display Type',...
                'Units','normalized',...
                'Position',[0.01 0.90 0.24 0.05],...
                'ForegroundColor',clr2,...
                'BackgroundColor',clr1,...
                'FontSize',12,...
                'HorizontalAlignment','left',...
                'Parent',RIM.F);
            uicontrol('Style','popupmenu',...
                'Tag','ViewType',...
                'String',{'Water signal','Fat signal','OP','Multispectral'},...
                'Value',4,...
                'Units','normalized',...
                'Position',[0.01 0.85 .24 0.05],...
                'Callback',@IMDisp,...
                'FontSize',12,...
                'Parent',RIM.F);


            RIM.bg1 = uibuttongroup('Visible','on',...
                'Position',[0.01 0.8 0.24 0.05],...
                'Title','Drawing tool',...
                'FontSize',12,...
                'ForegroundColor',clr2,...
                'BackgroundColor',clr1,...
                Parent=RIM.F);
            uicontrol('Style','radiobutton',...
                'String','MCKmeans (k)',...
                'Min',0,'Max',1,'Value',1,...
                'Units','normalized',...
                'Position',[0.025 0.05 0.45 0.9],...
                'ForegroundColor',clr2,...
                'BackgroundColor',clr1,...
                'FontSize',10,...
                'HorizontalAlignment','left',...
                'Visible','on',...
                Parent=RIM.bg1);
            uicontrol('Style','radiobutton',...
                'String','Manual contour (d)',...
                'Min',0,'Max',1,'Value',0,...
                'Units','normalized',...
                'Position',[0.525 0.05 0.45 0.9],...
                'ForegroundColor',clr2,...
                'BackgroundColor',clr1,...
                'FontSize',10,...
                'HorizontalAlignment','left',...
                'Visible','on',...
                Parent=RIM.bg1);

            RIM.bg2 = uibuttongroup('Visible','on',...
                'Position',[0.01 0.6 0.24 0.2],...
                'Title','Segmentation operations',...
                'FontSize',12,...
                'ForegroundColor',clr2,...
                'BackgroundColor',clr1,...
                Parent=RIM.F);
            uicontrol('Style','pushbutton',...
                'String','Erode (e)',...
                'Units','normalized',...
                'ForegroundColor',clr2,...
                'BackgroundColor',clr3,...
                'Position', [0.05 0.7625 0.425 0.1875],...
                'Callback',@contour_shrink,...
                'Parent',RIM.bg2);
            uicontrol('Style','pushbutton',...
                'String','Grow (g)',...
                'Units','normalized',...
                'ForegroundColor',clr2,...
                'BackgroundColor',clr3,...
                'Position', [0.525 0.7625 0.425 0.1875],...
                'Callback',@contour_grow,...
                'Parent',RIM.bg2);

            uicontrol('Style','pushbutton',...
                'String','Smooth 2D (s)',...
                'Units','normalized',...
                'ForegroundColor',clr2,...
                'BackgroundColor',clr3,...
                'Position', [0.05 0.525 0.425 0.1875],...
                'Callback',@Smooth2D,...
                'Parent',RIM.bg2);
            uicontrol('Style','pushbutton',...
                'String','Smooth 3D',...
                'Units','normalized',...
                'ForegroundColor',clr2,...
                'BackgroundColor',clr3,...
                'Position', [0.525 0.525 0.425 0.1875],...
                'Callback',@Smooth3D,...
                'Parent',RIM.bg2);

            uicontrol('Style','pushbutton',...
                'String','Remove Islands (i)',...
                'Units','normalized',...
                'Position', [0.05 0.2875 0.425 0.1875],...
                'ForegroundColor',clr2,...
                'BackgroundColor',clr3,...
                'Callback',@NoIslands3D,...
                'Parent',RIM.bg2);
            uicontrol('Style','pushbutton',...
                'String','View 3D (v)',...
                'Units','normalized',...
                'Position', [0.525 0.2875 0.425 0.1875],...
                'ForegroundColor',clr2,...
                'BackgroundColor',clr3,...
                'Callback',@ThreeDView,...
                'Parent',RIM.bg2);
            uicontrol('Style','pushbutton',...
                'String','Interpolate Segmentation',...
                'Units','normalize',...
                'ForegroundColor',clr2,...
                'BackgroundColor',clr3,...
                'Position',[0.275 0.05 0.425 0.1875],...
                'Callback',@interpolate,...
                Parent=RIM.bg2);
            % uicontrol('Style','pushbutton',...
            %     'String','Undo (ctrl-z)',...
            %     'Units','normalize',...
            %     'ForegroundColor',clr2,...
            %     'BackgroundColor',clr3,...
            %     'Position',[0.05 0.05 0.425 0.1875],...
            %     'Callback',@UndoMask,...
            %     Parent=RIM.bg2);
            % uicontrol('Style','pushbutton',...
            %     'String','Redo (ctrl-y)',...
            %     'Units','normalize',...
            %     'ForegroundColor',clr2,...
            %     'BackgroundColor',clr3,...
            %     'Position',[0.525 0.05 0.425 0.1875],...
            %     'Callback',@RedoMask,...
            %     Parent=RIM.bg2);

            

        %-----------------------------------------------------------------%
        %             Start new case, load DICOM information              %
        %-----------------------------------------------------------------%
        case 'New'
            % Disable clickable axes
            RIM.F.ButtonDownFcn = '';
            RIM.isdrawing = true;

            % Clear analysis fields if they exist
            analfields = {'H2O','Fat','OP','I','BW','prev','next'};
            for sf=1:length(analfields)
                if isfield(RIM,analfields{sf})
                    RIM = rmfield(RIM,analfields{sf});
                end
            end

            % Read in the images
            % Get user to select reference and dynamic images
            path = uigetdir(RIM.location,'Select folder containing 3D VANE NII images');
            if isequal(path,0)
                return
            else
                % store selected path as location for (probably) quicker
                % navigation on saving results/loading next sessions
                RIM.location = fileparts(path);
            end
            RIM.starttime = datetime;
            RIM.H2O.name = fullfile(path,'3D_VANE_W.nii');
            RIM.Fat.name = fullfile(path,'3D_VANE_F.nii');
            RIM.OP.name = fullfile(path,'3D_VANE_OP.nii');
            RIM.H2O.info = niftiinfo(RIM.H2O.name);

            LoadImagingData;


            % instantiate segmentation mask
            RIM.BW = false(size(RIM.I.H2O));

            % instantiate previous/future BW states
            RIM.prev = cell(0);
            RIM.next = cell(0);

            % Enable clickable axes, populate fields
            RIM.F.ButtonDownFcn = '@Seg_UI/AddPoint';
            RIM.isdrawing = false;
            RIM.onax = false;
            IMDisp;            


        case 'Save'
            % Stop the clock
            RIM.savetime(length(RIM.starttime)) = datetime;
            RIM.totalanalysistime = sum(seconds(RIM.savetime-RIM.starttime));
            % Select the file name and location to save
            defname = fullfile(fileparts(RIM.H2O.info.Filename),'MCseg');
            [savefile,savepath] = uiputfile('*.mat','Save file as...',defname);
            if isequal(savefile,0)
                return
            end
            % Create a structure with necessary fields
            savefields = {'starttime','savetime','totalanalysistime','H2O','Fat','OP','BW','prev','next'};
            SaveStruct = struct();
            for sf=1:length(savefields)
                SaveStruct.(savefields{sf}) = RIM.(savefields{sf});
            end
            save(fullfile(savepath,savefile),'SaveStruct')
            % Write the mask as a nifti
            minfo = RIM.H2O.info;
            minfo.Datatype = 'uint8';
            niftiwrite(rot90(flip(uint8(RIM.BW)),-1),fullfile(savepath,strrep(savefile,'.mat','.nii')),minfo);

        case 'Load'
            % Select the analysis file to load
            if isfield(RIM,'location')
                defname = RIM.location;
            else
                defname = pwd;
            end
            [loadfile,loadpath] = uigetfile('*.mat','Select file to load',defname);
            if isequal(loadfile,0)
                return
            end
            % Disable clickable axes
            RIM.F.ButtonDownFcn = '';
            RIM.isdrawing = true;
            load(fullfile(loadpath,loadfile),'SaveStruct')
            fldnms = fieldnames(SaveStruct);
            for fn=1:length(fldnms)
                RIM.(fldnms{fn}) = SaveStruct.(fldnms{fn});
            end
            % Start the clock on a new session
            RIM.starttime(end+1) = datetime;
            LoadImagingData;
            % Enable clickable axes
            RIM.F.ButtonDownFcn = '@Seg_UI/AddPoint';
            RIM.isdrawing = false;
            IMDisp;    



        case 'exit' 
            disp('Thanks for Viewing!')
            delete(RIM.F)
            clear RIM
            return
    end


    function LoadImagingData(~,~)
        w = waitbar(0,'Loading Images');
        niifiles = {RIM.H2O.name,RIM.Fat.name,RIM.OP.name};
        for nf=1:length(niifiles)
            futarr(nf) = parfeval(backgroundPool,@niftiread,1,niifiles{nf});
        end
        wait(futarr)
        RIM.I.H2O = mat2gray(flip(rot90(fetchOutputs(futarr(1)))));
        RIM.I.Fat = mat2gray(flip(rot90(fetchOutputs(futarr(2)))));
        RIM.I.OP = mat2gray(flip(rot90(fetchOutputs(futarr(3)))));
        
        % Create a background mask from the input images
        waitbar(0.5,w,'Detecting background');
        BKG = ~(imbinarize(RIM.I.H2O) | imbinarize(RIM.I.Fat) | imbinarize(RIM.I.OP));
        FatFrac = RIM.I.Fat./(RIM.I.Fat+RIM.I.H2O);
        FatFrac(BKG) = 0;
        BKG = BKG | FatFrac>.7;

        % Create Multichannel Image
        waitbar(0.75,w,'Decorrelating channels');
        % First decorrelate image - must operate on 2D n-channel image
        tempI = cat(4,RIM.I.H2O,RIM.I.Fat,RIM.I.OP);
        tempI = reshape(tempI,[],size(BKG,3),3);
        tempB = reshape(BKG,[],size(BKG,3));
        [row,col] = find(~tempB);
        tempI = decorrstretch(tempI,'Tol',0.0001,'SampleSubs',{row,col});
        RIM.I.RGB = reshape(tempI,[size(BKG) 3]);
        RIM.I.RGB = RIM.I.RGB.*~(repmat(BKG,1,1,1,3));


        % instantiate image display properties
        RIM.DataAspectRatio = RIM.H2O.info.PixelDimensions./RIM.H2O.info.PixelDimensions(1);
        RIM.zsli = round(size(RIM.I.RGB,3)/2);
        RIM.zslipre = 0;
        RIM.xsli = round(size(RIM.I.RGB,2)/2);
        RIM.xslipre = 0;
        RIM.ysli = round(size(RIM.I.RGB,1)/2);
        RIM.yslipre = 0;
        i = ismember({RIM.bg1.Children.String},'MCKmeans (k)');
        RIM.bg1.Children(i).Value = 1;
        RIM.ToolPointer = 'Crosshair';
        RIM.F.Pointer = RIM.ToolPointer;
        close(w);



    end


    %---------------------------------------------------------------------%
    %      Subroutine to display axial, sagittal, and coronal images      %
    %---------------------------------------------------------------------%
    function IMDisp(~,~)
        imtype = get(findobj(RIM.F,'Tag','ViewType'),'String');
        sel = get(findobj(RIM.F,'Tag','ViewType'),'Value');
        switch imtype{sel}
            case 'Water signal'
                img = RIM.I.H2O(:,:,RIM.zsli);
                cimg = rot90(squeeze(RIM.I.H2O(RIM.ysli,:,:)));
                simg = rot90(squeeze(RIM.I.H2O(:,RIM.xsli,:)));
                cmap = [1 1 0; 0 0 0];
            case 'Fat signal'
                img = RIM.I.Fat(:,:,RIM.zsli);
                cimg = rot90(squeeze(RIM.I.Fat(RIM.ysli,:,:)));
                simg = rot90(squeeze(RIM.I.Fat(:,RIM.xsli,:)));
                cmap = [1 1 0; 0 0 0];
            case 'OP'
                img = RIM.I.OP(:,:,RIM.zsli);
                cimg = rot90(squeeze(RIM.I.OP(RIM.ysli,:,:)));
                simg = rot90(squeeze(RIM.I.OP(:,RIM.xsli,:)));
                cmap = [1 1 0; 0 0 0];
            case 'Multispectral'
                img = squeeze(RIM.I.RGB(:,:,RIM.zsli,:));
                cimg = rot90(squeeze(RIM.I.RGB(RIM.ysli,:,:,:)));
                simg = rot90(squeeze(RIM.I.RGB(:,RIM.xsli,:,:)));
                cmap = [1 1 1; 0 0 0];
        end
        if RIM.zsli~=RIM.zslipre
            imshow(labeloverlay(img,bwperim(RIM.BW(:,:,RIM.zsli)),'Colormap',cmap),'Parent',RIM.v1)
            axis image
        end
        if RIM.ysli~=RIM.yslipre
            imshow(labeloverlay(cimg,rot90(bwperim(squeeze(RIM.BW(RIM.ysli,:,:)))),'Colormap',cmap),'Parent',RIM.v2)
            set(RIM.v2,'DataAspectRatio',RIM.DataAspectRatio([3 1 2]));
        end
        if RIM.xsli~=RIM.xslipre
            imshow(labeloverlay(simg,rot90(bwperim(squeeze(RIM.BW(:,RIM.xsli,:)))),'Colormap',cmap),'Parent',RIM.v3)
            set(RIM.v3,'DataAspectRatio',RIM.DataAspectRatio([3 1 2]));
        end
        drawnow
    end

    %---------------------------------------------------------------------%
    %        Subroutines for manual editing of segmentation result        %
    %---------------------------------------------------------------------%
    function AddPoint(src,cbd)
        if ~isfield(RIM,'I')
            % no imaging data loaded
            return
        end
        if ~RIM.onax
            return
        end
        
        if RIM.isdrawing
            return
        end
        cp = cbd.Source.CurrentAxes.CurrentPoint;
        c1 = round(cp(1,2));
        c2 = round(cp(1,1));
        switch RIM.axtag
            case 'ax1'
                sp = RIM.I.H2O(:,:,RIM.zsli);
            case 'ax2'
                sp = rot90(squeeze(RIM.I.H2O(RIM.ysli,:,:)));
            case 'ax3'
                sp = rot90(squeeze(RIM.I.H2O(:,RIM.xsli,:)));
        end
        
        
        try
            pixid = sub2ind(size(sp),c1,c2);
            switch RIM.axtag
                case 'ax2'
                    RIM.zslipre = RIM.zsli;
                    RIM.zsli = size(sp,1)-c1+1;
                    RIM.xslipre = RIM.xsli;
                    RIM.xsli = c2;
                    IMDisp;
                    return
                case 'ax3'
                    RIM.zslipre = RIM.zsli;
                    RIM.zsli = size(sp,1)-c1+1;
                    RIM.yslipre = RIM.ysli;
                    RIM.ysli = c2;
                    IMDisp;
                    return
            end
        catch
            % click is out of bounds
            fprintf('oob\n')
            return
        end

        switch RIM.bg1.SelectedObject.String
            case 'MCKmeans (k)'
                switch src.SelectionType
                    case 'normal'
                        % Left click: select and add
                        h = images.roi.Freehand('Color','y','Parent',gca,'Position',[c2 c1]);
                        h.beginDrawingFromPoint([c2 c1])
                        if isvalid(h)
                            roi = createMask(h);
                        else
                            return
                        end
                        delete(h)
                        if ~any(roi(:))
                            RIM.isdrawing = false;
                            return
                        end
                        backup;
                        altL = squeeze(RIM.I.RGB(:,:,RIM.zsli,:));
                        altL = altL.*repmat(roi,1,1,3);
                        [L,centers] = imsegkmeans(single(altL),4);
                        % detect and remove background
                        [~,i] = min(pdist2(centers,single([0 0 0])));
                        L = double(L);
                        L(L==i) = NaN;
                        if ~any(RIM.BW(:))
                            cdists = pdist2(centers,single([1 .1 .75]),'cosine');
                        else
                            target_color = mean(reshape(RIM.I.RGB(repmat(RIM.BW,1,1,1,3)),[],3));
                            cdists = pdist2(centers,single(target_color),'euclidean');
                        end
                        i = find(cdists==min(cdists(~ismember(1:length(cdists),i))));
                        L = L==i;
                        
                        
                        % i = mode(L(:));
                        % L = L==i;
                        L = contour_shrink(L);
                        L = NoIsland2D(L);
                        L = contour_grow(L);
                        L(~roi) = 0;
                        L = Smooth2D(L);
                        RIM.BW(:,:,RIM.zsli) = RIM.BW(:,:,RIM.zsli) | L;

                    case 'alt'
                        % Right click: select and remove
                        h = images.roi.Freehand('Color','r','Parent',gca,'Position',[c2 c1]);
                        h.beginDrawingFromPoint([c2 c1])
                        if isvalid(h)
                            rmm = createMask(h);
                        else
                            return
                        end
                        delete(h)
                        if ~any(rmm(:))
                            RIM.isdrawing = false;
                            return
                        end
                        backup;
                        RIM.BW(:,:,RIM.zsli) = min(RIM.BW(:,:,RIM.zsli),~rmm);
                end

            case 'Manual contour (d)'
                switch src.SelectionType
                    case 'normal'
                        % Left click: select and add
                        h = images.roi.Freehand('Color','y','Parent',gca,'Position',[c2 c1]);
                        h.beginDrawingFromPoint([c2 c1])
                        if isvalid(h)
                            roi = createMask(h);
                        else
                            return
                        end
                        delete(h)
                        if ~any(roi(:))
                            RIM.isdrawing = false;
                            return
                        end
                        backup;
                        RIM.BW(:,:,RIM.zsli) = max(RIM.BW(:,:,RIM.zsli),roi);
                    case 'alt'
                        % Right click: select and remove
                        h = images.roi.Freehand('Color','r','Parent',gca,'Position',[c2 c1]);
                        h.beginDrawingFromPoint([c2 c1])
                        if isvalid(h)
                            rmm = createMask(h);
                        else
                            return
                        end
                        delete(h)
                        if ~any(rmm(:))
                            RIM.isdrawing = false;
                            return
                        end
                        backup;
                        RIM.BW(:,:,RIM.zsli) = min(RIM.BW(:,:,RIM.zsli),~rmm);
                        

                end
                RIM.isdrawing = false;
                RIM.F.Pointer = RIM.ToolPointer;

        end
        % Update coronal and sagital images to center of ROI on current
        % slice
        centroid = regionprops(RIM.BW(:,:,RIM.zsli),{'Area','Centroid'});
        if ~isempty(centroid)
            [~,biggest] = max([centroid(:).Area]);
            RIM.xsli = round(centroid(biggest).Centroid(1));
            RIM.ysli = round(centroid(biggest).Centroid(2));
        end
        IMDisp;
    end

    function backup(~,~)
        % Add current BW state to history
        RIM.prev{end+1} = RIM.BW;
        % Only keep 10 in memory
        RIM.prev(1:end-10) = [];
    end

    function UndoMask(~,~)
        if isempty(RIM.prev)
            return
        end
        % Add current BW state to futures
        RIM.next{end+1} = RIM.BW;
        % Only keep 10 in memory
        RIM.next(1:end-10) = [];
        
        % Set current BW state to most recent history
        RIM.BW = RIM.prev{end};
        % Remove new state from history
        RIM.prev(end) = [];

        % Add code to make this jump to the slice where the change is
        % taking place?
        IMDisp;
    end

    function RedoMask(~,~)
        if isempty(RIM.next)
            return
        end
        % Add current BW state to history
        RIM.prev{end+1} = RIM.BW;
        % Only keep 10 in memory
        RIM.prev(1:end-10) = [];
        
        % Set current BW state to most recent future
        RIM.BW = RIM.next{end};
        % Remove new state from history
        RIM.next(end) = [];

        % Add code to make this jump to the slice where the change is
        % taking place?
        IMDisp;
    end

    function ResetSlice(~,~)
        RIM.BW(:,:,RIM.zsli) = 0;
        IMDisp;
    end

    function ResetStudy(~,~)
        RIM.BW = false(size(RIM.I.H2O));
        IMDisp;
    end





    function sbw = contour_grow(sbw,~)
        if nargin==0
            sbw = 'keyboardpress';
        end
        if ~islogical(sbw)
            backup;
            sbw = RIM.BW(:,:,RIM.zsli);
            updateRIM = true;
        else
            updateRIM = false;
        end
        imtype = get(findobj(RIM.F,'Tag','ViewType'),'String');
        sel = get(findobj(RIM.F,'Tag','ViewType'),'Value');
        switch imtype{sel}
            case 'Water signal'
                imcon = RIM.I.H2O(:,:,RIM.zsli);
            case 'Fat signal'
                imcon = RIM.I.Fat(:,:,RIM.zsli);
            case 'OP'
                imcon = RIM.I.OP(:,:,RIM.zsli);
            case 'Multispectral'
                YIQ = rgb2ntsc(squeeze(RIM.I.RGB(:,:,RIM.zsli,:)));
                imcon = YIQ(:,:,3);
        end
        sbw = activecontour(imcon,sbw,5,'Chan-Vese','ContractionBias',0,'SmoothFactor',.1);
        if updateRIM
            RIM.BW(:,:,RIM.zsli) = sbw;
            IMDisp;
        end
    end

    function sbw = contour_shrink(sbw,~)
        if nargin==0
            sbw = 'keyboardpress';
        end
        if ~islogical(sbw)
            backup;
            sbw = RIM.BW(:,:,RIM.zsli);
            updateRIM = true;
        else
            updateRIM = false;
        end
        imtype = get(findobj(RIM.F,'Tag','ViewType'),'String');
        sel = get(findobj(RIM.F,'Tag','ViewType'),'Value');
        switch imtype{sel}
            case 'Water signal'
                imcon = RIM.I.H2O(:,:,RIM.zsli);
            case 'Fat signal'
                imcon = RIM.I.Fat(:,:,RIM.zsli);
            case 'OP'
                imcon = RIM.I.OP(:,:,RIM.zsli);
            case 'Multispectral'
                YIQ = rgb2ntsc(squeeze(RIM.I.RGB(:,:,RIM.zsli,:)));
                imcon = YIQ(:,:,3);
        end
        sbw = activecontour(imcon,sbw,5,'Chan-Vese','ContractionBias',.75,'SmoothFactor',0.1);
        if updateRIM
            RIM.BW(:,:,RIM.zsli) = sbw;
            IMDisp;
        end
    end

    function sbw = Smooth2D(sbw,~)
        if nargin==0
            sbw = 'keyboardpress';
        end
        if ~islogical(sbw)
            backup;
            updateRIM = true;
            sbw = RIM.BW(:,:,RIM.zsli);
        else
            updateRIM = false;
        end
        fwhm = 3;
        sigma = fwhm/(2*sqrt(2*log(2)))/RIM.H2O.info.PixelDimensions(1);
        sbw = round(imgaussfilt(double(sbw),sigma))==1;
        if updateRIM
            RIM.BW(:,:,RIM.zsli) = sbw;
            IMDisp;
        end
    end

    function Smooth3D(~,~)
        backup;
        fwhm = 3;
        sigma = repmat(fwhm/(2*sqrt(2*log(2))),1,3)./RIM.H2O.info.PixelDimensions;
        RIM.BW = round(imgaussfilt3(double(RIM.BW),sigma))==1;
        IMDisp;
    end





    function sbw = NoIsland2D(sbw,~)
        % Remove islands from the input slice
        if nargin==0
            sbw = 'keyboardpress';
        end
        if ~islogical(sbw)
            backup;
            sbw = RIM.BW(:,:,RIM.zsli);
            updateRIM = true;
        else
            updateRIM = false;
        end
        [L,in] = bwlabel(sbw,4);
        hc = histcounts(L,1:1:in+1);
        [~,in] = max(hc);
        sbw = L==in;
        if updateRIM
            RIM.BW(:,:,RIM.zsli) = sbw;
            IMDisp;
        end
    end

    function NoIslands3D(~,~)
        backup;
        % Remove islands
        [L,in] = bwlabeln(RIM.BW);
        hc = histcounts(L,1:1:in+1);
        [~,in] = max(hc);
        RIM.BW = L==in;
        IMDisp;
    end

    function sbw = fill2D(sbw,~)
        if nargin==0
            sbw = 'keyboardpress';
        end
        if ~islogical(sbw)
            backup;
            sbw = RIM.BW(:,:,RIM.zsli);
            updateRIM = true;
        else
            updateRIM = false;
        end
        sbw = imfill(sbw,'holes');
        if updateRIM
            RIM.BW(:,:,RIM.zsli) = sbw;
            IMDisp;
        end
    end

    function fill3D(~,~)
        backup;
        RIM.BW = imfill(RIM.BW,'holes');
    end

    function interpolate(~,~)
        backup;
        sidx = find(squeeze(sum(sum(RIM.BW))));
        for n=sidx(1):sidx(end)
            if ismember(n,sidx)
                % original label is on this slice
                continue
            end
            % Get the coordinates along the perimeter of lower and higher
            % thresholds
            lwr = RIM.BW(:,:,n-1);
            upr = RIM.BW(:,:,n+1);
            if ~any(lwr(:)) || ~any(upr(:))
                continue
            end
            [x1,y1] = find(bwperim(lwr));
            [x2,y2] = find(bwperim(upr));
            % get pairwise distances
            dist = pdist2([x1,y1],[x2,y2]);
            [~,i] = min(dist,[],2);
            tgx = x2(i);
            tgy = y2(i);
            newx = (mean([x1 tgx],2));
            newy = (mean([y1 tgy],2));
    
            % Re-order points to make ROI curve
            dist = squareform(pdist([newx newy]));
            dist(eye(size(dist))==1) = NaN;
            o = 1;
            for p=1:length(dist)
                [pd(p),i] = min(dist(:,o(end)));
                dist(o(end),:) = NaN;
                dist(:,o(end)) = NaN;
                if pd(p)<5
                    o = cat(1,o,i);
                end
            end
            slice = poly2mask(newy(o),newx(o),size(lwr,1),size(lwr,2));
            slice(upr) = true;
            slice(~lwr) = false;
            RIM.BW(:,:,n) = slice;
            
        end
        IMDisp;
    end




    function ThreeDView
        RIM.vol = volshow(RIM.I.RGB, ...
            'RenderingStyle','GradientOpacity',...
            'Transformation',RIM.H2O.info.Transform.T',...
            'OverlayData',RIM.BW,'OverlayAlphamap',1);
        RIM.vol.Alphamap = RIM.vol.Alphamap.^2;
        % % resample to 1mm isotropic
        % V = imresize3(RIM.BW,round(size(RIM.BW).*RIM.DataAspectRatio),'nearest');
        % 
        % fv = isosurface(V,0);
        % figure
        % patch('faces',fv.faces,'vertices',fv.vertices,...
        %     'FaceColor',[1 .5 .5],'EdgeColor','none',...
        %     'FaceLighting','gouraud');
        % axis equal off
        % rotate3d on
        % 
        % view(0,0)
        % 
        % l1 = light;
        % l1.Color = [.25 .25 .25];
        % l1.Position = [1 1 1];
        % l2 = light;
        % l2.Color = [.25 .25 .25];
        % l2.Position = [-1 -1 -1];
        % l3 = light;
        % l3.Color = [.25 .25 .25];
        % l3.Position = [-1 -1 1];
        % l4 = light;
        % l4.Color = [.25 .25 .25];
        % l4.Position = [1 1 -1];
        % 
        % l5 = light;
        % l5.Color = [.25 .25 .25];
        % l5.Position = [1 -1 1];
        % l6 = light;
        % l6.Color = [.25 .25 .25];
        % l6.Position = [-1 1 -1];
        % l7 = light;
        % l7.Color = [.25 .25 .25];
        % l7.Position = [-1 1 1];
        % l8 = light;
        % l8.Color = [.25 .25 .25];
        % l8.Position = [1 -1 -1];
    end

    %---------------------------------------------------------------------%
    %        Sub-routine to scroll through slices with mouse wheel        %
    %---------------------------------------------------------------------%
    function scroll(~,callbackdata)
        if ~isfield(RIM,'isdrawing')             % if no data is loaded
            return
        end
        if RIM.isdrawing || ~RIM.onax
            return
        end
        cp = callbackdata.Source.CurrentAxes.CurrentPoint;
        c1 = round(cp(1,2));
        c2 = round(cp(1,1));
        sp = RIM.I.H2O(:,:,RIM.zsli);
        % try
        %     pixid = sub2ind(size(sp),c1,c2);
        % catch
        %     % click is out of bounds
        %     return
        % end
        switch RIM.axtag
            case 'ax1'
                if callbackdata.VerticalScrollCount > 0
                    if RIM.zsli>1
                        RIM.zslipre = RIM.zsli;
                        RIM.zsli = RIM.zsli-1;
                        IMDisp;
                    end
                elseif callbackdata.VerticalScrollCount < 0
                    if RIM.zsli<size(RIM.I.RGB,3)
                        RIM.zslipre = RIM.zsli;
                        RIM.zsli  = RIM.zsli+1;
                        IMDisp;
                    end
                end
            case 'ax2'
                if callbackdata.VerticalScrollCount > 0
                    if RIM.ysli>1
                        RIM.yslipre = RIM.ysli;
                        RIM.ysli = RIM.ysli-1;
                        IMDisp;
                    end
                elseif callbackdata.VerticalScrollCount < 0
                    if RIM.ysli<size(RIM.I.RGB,1)
                        RIM.yslipre = RIM.ysli;
                        RIM.ysli  = RIM.ysli+1;
                        IMDisp;
                    end
                end
            case 'ax3'
                if callbackdata.VerticalScrollCount > 0
                    if RIM.xsli>1
                        RIM.xslipre = RIM.xsli;
                        RIM.xsli = RIM.xsli-1;
                        IMDisp;
                    end
                elseif callbackdata.VerticalScrollCount < 0
                    if RIM.xsli<size(RIM.I.RGB,2)
                        RIM.xslipre = RIM.xsli;
                        RIM.xsli  = RIM.xsli+1;
                        IMDisp;
                    end
                end
        end
    end

    %---------------------------------------------------------------------%
    %   Update current axis based on cursor position when mouse moves     %
    %---------------------------------------------------------------------%
    function mouseMove(~,cbd)
        if ~isfield(RIM,'isdrawing')             % if no data is loaded
            return
        end
        
        if RIM.isdrawing
            return
        end
        % try ax1 first
        set(RIM.F,'CurrentAxes',RIM.v1)
        sp = RIM.I.H2O(:,:,RIM.zsli);
        
        C = get(gca,'CurrentPoint');
        c1 = round(C(1,2));
        c2 = round(C(1,1));     
        try
            pixid = sub2ind(size(sp),c1,c2);
            RIM.axtag = 'ax1';
        catch
            % cursor is out of bounds for current axis, try swapping
            set(RIM.F,'CurrentAxes',RIM.v2)
            sp = rot90(squeeze(RIM.I.H2O(RIM.ysli,:,:)));
            C = get(gca,'CurrentPoint');
            c1 = round(C(1,2));
            c2 = round(C(1,1));     
            try
                pixid = sub2ind(size(sp),c1,c2);
                RIM.axtag = 'ax2';
            catch
                % cursor is out of bounds for current axis, try axis 3
                set(RIM.F,'CurrentAxes',RIM.v3)
                sp = rot90(squeeze(RIM.I.H2O(:,RIM.xsli,:)));
                C = get(gca,'CurrentPoint');
                c1 = round(C(1,2));
                c2 = round(C(1,1));     
                try
                    pixid = sub2ind(size(sp),c1,c2);
                    RIM.axtag = 'ax3';
                catch
                    % Cursor not on any axis
                    RIM.F.Pointer = "arrow";
                    RIM.onax = false;
                    return
                end
            end
        end
        RIM.onax = true;
        if isfield(RIM,'ToolPointer')
            RIM.F.Pointer = RIM.ToolPointer;
        else
            RIM.F.Pointer = "arrow";
        end

    end

    %---------------------------------------------------------------------%
    %                         Keyboard shortcuts!                         %
    %---------------------------------------------------------------------%
    function keyitup(~,callbackdata)
        
        if ~isfield(RIM,'H2O')             % if no data is loaded
            return
        end
        if ~isempty(callbackdata.Modifier)
            if strcmp(callbackdata.Key,'z') && strcmp(callbackdata.Modifier,'control')
                UndoMask
            end
            if strcmp(callbackdata.Key,'y') && strcmp(callbackdata.Modifier,'control')
                RedoMask
            end
            if strcmp(callbackdata.Key,'s') && strcmp(callbackdata.Modifier,'control')
                Seg_UI('Save')
            end
            if strcmp(callbackdata.Key,'l') && strcmp(callbackdata.Modifier,'control')
                Seg_UI('Load')
            end
            if strcmp(callbackdata.Key,'r') && strcmp(callbackdata.Modifier,'control')
                ResetStudy
            end
            if strcmp(callbackdata.Key,'s') && strcmp(callbackdata.Modifier,'shift')
                Smooth3D;
            end
        end
        
        switch callbackdata.Character
            case {'d'}
                chname = {RIM.bg1.Children.String};
                i = ismember(chname,'Manual contour (d)');
                RIM.bg1.Children(i).Value = 1;
            case {'k'}
                chname = {RIM.bg1.Children.String};
                i = ismember(chname,'MCKmeans (k)');
                RIM.bg1.Children(i).Value = 1;
            case {'e'}
                contour_shrink;
            case {'g'}
                contour_grow;
            case {'s'}
                Smooth2D;
            case {'i'}
                NoIslands3D;
            case {'v'}
                ThreeDView;
            case {'r'}
                ResetSlice;

            otherwise
                return
        end
    end


end