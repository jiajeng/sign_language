function o = checkGUI(app,local_flag,o,connF,ISCF)
    steps = dictionary('proces_use_procD',struct("Setup",1, ...
                                              "Preprocessing",0, ...
                                              "Denoising",1, ...
                                              "fst_Analysis",1, ...
                                              "snd_Analysis",1, ...
                                              "Add_Roi",0) ...
                      ,'Denois_use_procD',struct("Setup",1, ...
                                              "Preprocessing",0, ...
                                              "Denoising",1, ...
                                              "fst_Analysis",0, ...
                                              "snd_Analysis",0, ...
                                              "Add_Roi",0) ...
                      ,'proces_new',struct("Setup",1, ...
                                          "Preprocessing",1, ...
                                          "Denoising",1, ...
                                          "fst_Analysis",1, ...
                                          "snd_Analysis",1, ...
                                          "Add_Roi",0) ...
                      ,'Denois_new',struct("Setup",1, ...
                                          "Preprocessing",1, ...
                                          "Denoising",1, ...
                                          "fst_Analysis",0, ...
                                          "snd_Analysis",0, ...
                                          "Add_Roi",0) ...
                     ,'Add_roi',struct("Setup",0, ...
                                       "Preprocessing",0, ...
                                       "Denoising",0, ...
                                       "fst_Analysis",0, ...
                                       "snd_Analysis",0, ...
                                       "Add_Roi",1));
    BasePixelX = 2560;
    BasePixelY = 1440;
    CurrentPixel = get(0,'ScreenSize');
    CurrentPixelX = CurrentPixel(3);
    CurrentPixelY = CurrentPixel(4);
    ResizeX = CurrentPixelX/BasePixelX;
    ResizeY = CurrentPixelY/BasePixelY;
    resizeF = ResizeX;
    % Create UIFigure and hide until all components are created
    app.UIFigure = uifigure('Visible', 'off');
    app.UIFigure.Position = [100*ResizeX 100*ResizeY 1415*ResizeX 1100*ResizeY];
    app.UIFigure.Name = 'MATLAB App';
    app.UIFigure.UserData = o;
    %% -----------------preprocess File panel------------------------
    % Create preprocessfiledefinePanel
    app.preprocessfiledefinePanel = uipanel(app.UIFigure);
    app.preprocessfiledefinePanel.BorderColor = [0.302 0.7451 0.9333];
    app.preprocessfiledefinePanel.BorderWidth = 3;
    app.preprocessfiledefinePanel.Title = 'preprocess file define';
    app.preprocessfiledefinePanel.BackgroundColor = [1 1 1];
    app.preprocessfiledefinePanel.FontAngle = 'italic';
    app.preprocessfiledefinePanel.FontWeight = 'bold';
    app.preprocessfiledefinePanel.FontSize = 18*resizeF;
    app.preprocessfiledefinePanel.Position = [10*ResizeX 542*ResizeY 674*ResizeX 550*ResizeY];

    % Create subjectTable
    if class(o.sub) ~= "cell", o.sub = convertStringsToChars(o.sub); end
    tab = o.sub;
    if size(tab,2) > 1,tab = tab'; end
    tab = cat(1,tab,repmat({""},10,1));
    tab = table(tab); tab.Properties.VariableNames = {'subject'};
    app.subjectTable = uitable(app.preprocessfiledefinePanel,"Data",tab,"FontSize",12*resizeF);
    app.subjectTable.ColumnName = {'subject'};
    app.subjectTable.RowName = {};
    app.subjectTable.Position = [26*ResizeX 6*ResizeY 100*ResizeX 253*ResizeY];
    app.subjectTable.ColumnEditable = true;

    % Create sessionTable
    tab = o.sess;
    if size(tab,2) > 1,tab = tab'; end
    tab = cat(1,tab,repmat({""},10,1));
    tab = table(tab); tab.Properties.VariableNames = {'sess'};
    app.sessionTable = uitable(app.preprocessfiledefinePanel,"Data",tab,"FontSize",12*resizeF);
    app.sessionTable.ColumnName = {'sess'};
    app.sessionTable.RowName = {};
    app.sessionTable.Position = [152*ResizeX 6*ResizeY 100*ResizeX 253*ResizeY];
    app.sessionTable.ColumnEditable = true;


    % Create roundTable
    tab = o.round;
    if size(tab,2) > 1,tab = tab'; end
    tab = cat(1,tab,repmat({""},10,1));
    tab = table(tab); tab.Properties.VariableNames = {'round'};
    app.roundTable = uitable(app.preprocessfiledefinePanel,"Data",tab,"FontSize",12*resizeF);
    app.roundTable.ColumnName = {'round'};
    app.roundTable.RowName = {};
    app.roundTable.Position = [271*ResizeX 6*ResizeY 100*ResizeX 253*ResizeY];
    app.roundTable.ColumnEditable = true;

    if local_flag
        % Create localsubjectpathLabel
        app.localsubjectpathLabel = uilabel(app.preprocessfiledefinePanel,"FontSize",12*resizeF);
        app.localsubjectpathLabel.Interpreter = 'latex';
        app.localsubjectpathLabel.HorizontalAlignment = 'right';
        app.localsubjectpathLabel.Position = [31*ResizeX 486*ResizeY 106*ResizeX 22*ResizeY];
        app.localsubjectpathLabel.Text = 'local subject path';

        % Create localsubjectpath
        app.localsubjectpath = uieditfield(app.preprocessfiledefinePanel, "Value",o.subpath,'HorizontalAlignment',"right","FontSize",12*ResizeY);
        app.localsubjectpath.Position = [152*ResizeX 486*ResizeY 491*ResizeX 22*ResizeY];

        % Create localoutputpathLabel
        app.localoutputpathLabel = uilabel(app.preprocessfiledefinePanel,"FontSize",12*resizeF);
        app.localoutputpathLabel.Interpreter = 'latex';
        app.localoutputpathLabel.HorizontalAlignment = 'right';
        app.localoutputpathLabel.Position = [40*ResizeX 453*ResizeY 94*ResizeX 22*ResizeY];
        app.localoutputpathLabel.Text = 'local outputpath';

        % Create localoutputpath
        app.localoutputpath = uieditfield(app.preprocessfiledefinePanel,"Value",o.outsubpath,'HorizontalAlignment',"right","FontSize",12*ResizeY);
        app.localoutputpath.Position = [152*ResizeX 453*ResizeY 491*ResizeX 22*ResizeY];
    else
        % Create ftpnasipLabel
        app.ftpnasipLabel = uilabel(app.preprocessfiledefinePanel,"FontSize",12*resizeF);
        app.ftpnasipLabel.Interpreter = 'latex';
        app.ftpnasipLabel.HorizontalAlignment = 'right';
        app.ftpnasipLabel.Position = [71*ResizeX 486*ResizeY 66*ResizeX 22*ResizeY];
        app.ftpnasipLabel.Text = 'ftp(nas) ip';
        
        % Create ftpnasip
        app.ftpnasip = uieditfield(app.preprocessfiledefinePanel,"Value", o.ftpServer.ip,'HorizontalAlignment',"right","FontSize",12*resizeF);
        app.ftpnasip.Position = [152*ResizeX 486*ResizeY 491*ResizeX 22*ResizeY];
        
        % Create ftpnasaccountLabel
        app.ftpnasaccountLabel = uilabel(app.preprocessfiledefinePanel,"FontSize",12*resizeF);
        app.ftpnasaccountLabel.Interpreter = 'latex';
        app.ftpnasaccountLabel.HorizontalAlignment = 'right';
        app.ftpnasaccountLabel.Position = [38*ResizeX 453*ResizeY 99*ResizeX 22*ResizeY];
        app.ftpnasaccountLabel.Text = 'ftp(nas) account';
        
        % Create ftpnasaccount
        app.ftpnasaccount = uieditfield(app.preprocessfiledefinePanel, "Value",o.ftpServer.account,'HorizontalAlignment',"right","FontSize",12*resizeF);
        app.ftpnasaccount.Position = [152*ResizeX 453*ResizeY 491*ResizeX 22*ResizeY];
        
        % Create ftpnaspasswordLabel
        app.ftpnaspasswordLabel = uilabel(app.preprocessfiledefinePanel,"FontSize",12*resizeF);
        app.ftpnaspasswordLabel.Interpreter = 'latex';
        app.ftpnaspasswordLabel.HorizontalAlignment = 'right';
        app.ftpnaspasswordLabel.Position = [30*ResizeX 417*ResizeY 107*ResizeX 22*ResizeY];
        app.ftpnaspasswordLabel.Text = 'ftp(nas) password';
        
        % Create ftpnaspassword
        app.ftpnaspassword = uieditfield(app.preprocessfiledefinePanel, "Value",o.ftpServer.password,'HorizontalAlignment',"right","FontSize",12*resizeF);
        app.ftpnaspassword.Position = [152*ResizeX 417*ResizeY 491*ResizeX 22*ResizeY];
        
        % Create ftpnasinputfolderLabel
        app.ftpnasinputfolderLabel = uilabel(app.preprocessfiledefinePanel,"FontSize",12*resizeF);
        app.ftpnasinputfolderLabel.Interpreter = 'latex';
        app.ftpnasinputfolderLabel.HorizontalAlignment = 'right';
        app.ftpnasinputfolderLabel.Position = [20*ResizeX 379*ResizeY 117*ResizeX 22*ResizeY];
        app.ftpnasinputfolderLabel.Text = 'ftp(nas) inputfolder';
        
        % Create ftpnasinputfolder
        app.ftpnasinputfolder = uieditfield(app.preprocessfiledefinePanel, "Value",o.ftpServer.infolder,'HorizontalAlignment',"right","FontSize",12*resizeF);
        app.ftpnasinputfolder.Position = [152*ResizeX 379*ResizeY 491*ResizeX 22*ResizeY];
        
        % Create ftpnasoutputfolderLabel
        app.ftpnasoutputfolderLabel = uilabel(app.preprocessfiledefinePanel,"FontSize",12*resizeF);
        app.ftpnasoutputfolderLabel.Interpreter = 'latex';
        app.ftpnasoutputfolderLabel.HorizontalAlignment = 'right';
        app.ftpnasoutputfolderLabel.Position = [13*ResizeX 344*ResizeY 124*ResizeX 22*ResizeY];
        app.ftpnasoutputfolderLabel.Text = 'ftp(nas) outputfolder';
        
        % Create ftpnasoutputfolder
        app.ftpnasoutputfolder = uieditfield(app.preprocessfiledefinePanel, "Value",o.ftpServer.outfolder,'HorizontalAlignment',"right","FontSize",12*resizeF);
        app.ftpnasoutputfolder.Position = [152*ResizeX 344*ResizeY 491*ResizeX 22*ResizeY];
        
        % Create localsubjectpathEditFieldLabel
        app.localoutputpathLabel = uilabel(app.preprocessfiledefinePanel,"FontSize",12*resizeF);
        app.localoutputpathLabel.Interpreter = 'latex';
        app.localoutputpathLabel.HorizontalAlignment = 'right';
        app.localoutputpathLabel.Position = [30*ResizeX 310*ResizeY 106*ResizeX 22*ResizeY];
        app.localoutputpathLabel.Text = 'local subject path';
        
        % Create localsubjectpathEditField
        app.localoutputpath = uieditfield(app.preprocessfiledefinePanel,"Value",o.outsubpath,'HorizontalAlignment',"right","FontSize",12*resizeF);
        app.localoutputpath.Position = [151*ResizeX 310*ResizeY 492*ResizeX 22*ResizeY];

    end
    % Create foldernestEditFieldLabel_2
    app.foldernestEditFieldLabel = uilabel(app.preprocessfiledefinePanel,"FontSize",12*resizeF);
    app.foldernestEditFieldLabel.Interpreter = 'latex';
    app.foldernestEditFieldLabel.Position = [412*ResizeX 196*ResizeY 235*ResizeX 90*ResizeY];
    app.foldernestEditFieldLabel.Text = {'folder nest '; 'define the folder nest  under suject folder'; 'sess and round are define by the table'; ''; 'e.x. sess/mri/round'};

    % Create Panel
    app.Panel = uipanel(app.preprocessfiledefinePanel,"FontSize",12*resizeF);
    app.Panel.BorderColor = [1 1 1];
    app.Panel.BackgroundColor = [1 1 1];
    app.Panel.Position = [396*ResizeX 6*ResizeY 268*ResizeX 150*ResizeY];

    % Create NestfolderEditFieldLabel
    app.NestfolderEditFieldLabel = uilabel(app.Panel,"FontSize",12*resizeF);
    app.NestfolderEditFieldLabel.HorizontalAlignment = 'right';
    app.NestfolderEditFieldLabel.Position = [5*ResizeX 115*ResizeY 70*ResizeX 22*ResizeY];
    app.NestfolderEditFieldLabel.Text = 'Nest folder';

    % Create NestfolderEditField
    app.NestfolderEditField = uieditfield(app.Panel,"text", "Value",strjoin(o.folder_nest,filesep),"FontSize",12*resizeF);
    app.NestfolderEditField.Position = [83*ResizeX 115*ResizeY 164*ResizeX 22*ResizeY];

    % % Create Tree
    % nest = o.folder_nest;
    % nest = split(nest,filesep);
    % app = createFolderNest(app,nest);
    %
    % % Create ButtonGroup
    % app.ButtonGroup = uibuttongroup(app.preprocessfiledefinePanel,"SelectionChangedFcn", @(src,evt) ButtonGroupSelectionChanged(src,evt,o,app));
    % app.ButtonGroup.BorderColor = [1 1 1];
    % app.ButtonGroup.BackgroundColor = [1 1 1];
    % app.ButtonGroup.Position = [525 159 123 67];
    % 
    % % Create selectfolderButton
    % app.selectfolderButton = uiradiobutton(app.ButtonGroup);
    % app.selectfolderButton.Text = 'select folder';
    % app.selectfolderButton.Position = [11 21 87 22];
    % 
    % % Create keyfilepathButton
    % app.keyfilepathButton = uiradiobutton(app.ButtonGroup);
    % app.keyfilepathButton.Text = 'key file path';
    % app.keyfilepathButton.Position = [11 -1 86 22];
    % app.keyfilepathButton.Value = true;
    % -----------------------------------------------------------------------

    %% -----------------CONN Project panel-------------------------
    o.round(contains(o.round,"T1")) = [];
    % Create CONNvariablePanel
    app.CONNvariablePanel = uipanel(app.UIFigure);
    app.CONNvariablePanel.BorderColor = [0.302 0.7451 0.9333];
    app.CONNvariablePanel.BorderWidth = 3;
    app.CONNvariablePanel.Title = 'CONN variable';
    app.CONNvariablePanel.BackgroundColor = [1 1 1];
    app.CONNvariablePanel.FontAngle = 'italic';
    app.CONNvariablePanel.FontWeight = 'bold';
    app.CONNvariablePanel.FontSize = 18*resizeF;
    app.CONNvariablePanel.Position = [12*ResizeX 12*ResizeY 673*ResizeX 515*ResizeY];
    
    % Create projectNameEditFieldLabel
    app.projectNameEditFieldLabel = uilabel(app.CONNvariablePanel,"FontSize",12*resizeF);
    app.projectNameEditFieldLabel.Interpreter = 'latex';
    app.projectNameEditFieldLabel.HorizontalAlignment = 'right';
    app.projectNameEditFieldLabel.Position = [8*ResizeX 446*ResizeY 83*ResizeX 22*ResizeY];
    app.projectNameEditFieldLabel.Text = 'project Name';
    
    % Create projectNameEditField
    app.projectNameEditField = uieditfield(app.CONNvariablePanel, 'text',"Value",o.conn_prjName,"FontSize",12*resizeF);
    app.projectNameEditField.Position = [111*ResizeX 446*ResizeY 543*ResizeX 22*ResizeY];
    
    % Create projectpathEditFieldLabel
    app.projectpathEditFieldLabel = uilabel(app.CONNvariablePanel,"FontSize",12*ResizeY);
    app.projectpathEditFieldLabel.Interpreter = 'latex';
    app.projectpathEditFieldLabel.HorizontalAlignment = 'right';
    app.projectpathEditFieldLabel.Position = [17*ResizeX 414*ResizeY 74*ResizeX 22*ResizeY];
    app.projectpathEditFieldLabel.Text = 'project path';
    
    % Create projectpathEditField
    prjPath = fullfile(o.conn_prjPath,o.conn_prjName);
    app.projectpathEditField = uieditfield(app.CONNvariablePanel, 'text',"Value",prjPath,"FontSize",12*resizeF);
    app.projectpathEditField.Position = [111*ResizeX 414*ResizeY 543*ResizeX 22*ResizeY];

    % ----------------------------------------define ROI
    % Create userdefineROIpathEditFieldLabel
    app.userdefineROIpathEditFieldLabel = uilabel(app.CONNvariablePanel,"FontSize",12*resizeF);
    app.userdefineROIpathEditFieldLabel.Interpreter = 'latex';
    app.userdefineROIpathEditFieldLabel.HorizontalAlignment = 'right';
    app.userdefineROIpathEditFieldLabel.Position = [205*ResizeX 107*ResizeY 124*ResizeX 22*ResizeY];
    app.userdefineROIpathEditFieldLabel.Text = 'user define ROI path';
    app.userdefineROIpathEditFieldLabel.Visible = false;


    % Create userdefineROIpathEditField
    app.userdefineROIpathEditField = uieditfield(app.CONNvariablePanel, 'text',"Value",o.conn_ROIpath,"FontSize",12*resizeF);
    app.userdefineROIpathEditField.Position = [349*ResizeX 107*ResizeY 165*ResizeX 22*ResizeY];
    app.userdefineROIpathEditField.Visible = false;
    % --------------------------------------------------
        
    % ----------------------------------preprocessing
    % Create smoothkernelEditFieldLabel
    app.smoothkernelEditFieldLabel = uilabel(app.CONNvariablePanel,"FontSize",12*resizeF);
    app.smoothkernelEditFieldLabel.Interpreter = 'latex';
    app.smoothkernelEditFieldLabel.HorizontalAlignment = 'right';
    app.smoothkernelEditFieldLabel.Position = [283*ResizeX 261*ResizeY 88*ResizeX 22*ResizeY];
    app.smoothkernelEditFieldLabel.Text = 'smooth kernel';
    app.smoothkernelEditFieldLabel.Visible = false;

    % Create smoothkernelEditField
    app.smoothkernelEditField = uieditfield(app.CONNvariablePanel, 'numeric','Value',o.conn_fwhm,"FontSize",12*resizeF);
    app.smoothkernelEditField.Position = [391*ResizeX 261*ResizeY 121*ResizeX 22*ResizeY];
    app.smoothkernelEditField.Visible = false;
    
    % Create StructurevoxelresolutionEditFieldLabel
    app.StructurevoxelresolutionEditFieldLabel = uilabel(app.CONNvariablePanel,"FontSize",12*resizeF);
    app.StructurevoxelresolutionEditFieldLabel.Interpreter = 'latex';
    app.StructurevoxelresolutionEditFieldLabel.HorizontalAlignment = 'right';
    app.StructurevoxelresolutionEditFieldLabel.Position = [219*ResizeX 229*ResizeY 152*ResizeX 22*ResizeY];
    app.StructurevoxelresolutionEditFieldLabel.Text = 'Structure voxel resolution';
    app.StructurevoxelresolutionEditFieldLabel.Visible = false;
    
    % Create StructurevoxelresolutionEditField
    app.StructurevoxelresolutionEditField = uieditfield(app.CONNvariablePanel, 'numeric','Value',o.conn_strVres,"FontSize",12*resizeF);
    app.StructurevoxelresolutionEditField.Position = [391*ResizeX 229*ResizeY 121*ResizeX 22*ResizeY];
    app.StructurevoxelresolutionEditField.Visible = false;
    
    % Create FunctionalvoxelresolutionEditFieldLabel
    app.FunctionalvoxelresolutionEditFieldLabel = uilabel(app.CONNvariablePanel,"FontSize",12*resizeF);
    app.FunctionalvoxelresolutionEditFieldLabel.Interpreter = 'latex';
    app.FunctionalvoxelresolutionEditFieldLabel.HorizontalAlignment = 'right';
    app.FunctionalvoxelresolutionEditFieldLabel.Position = [210*ResizeX 197*ResizeY 161*ResizeX 22*ResizeY];
    app.FunctionalvoxelresolutionEditFieldLabel.Text = 'Functional voxel resolution';
    app.FunctionalvoxelresolutionEditFieldLabel.Visible = false;
    
    % Create FunctionalvoxelresolutionEditField
    app.FunctionalvoxelresolutionEditField = uieditfield(app.CONNvariablePanel, 'numeric','Value',o.conn_funcVres,"FontSize",12*resizeF);
    app.FunctionalvoxelresolutionEditField.Position = [391*ResizeX 197*ResizeY 121*ResizeX 22*ResizeY];
    app.FunctionalvoxelresolutionEditField.Visible = false;
    

    % ---------------------------------------------------
   
    % ------------------------------Denoising

    % Create filterBandEditFieldLabel
    app.filterBandEditFieldLabel = uilabel(app.CONNvariablePanel,"FontSize",12*resizeF);
    app.filterBandEditFieldLabel.Interpreter = 'latex';
    app.filterBandEditFieldLabel.HorizontalAlignment = 'right';
    app.filterBandEditFieldLabel.Position = [304*ResizeX 167*ResizeY 67*ResizeX 22*ResizeY];
    app.filterBandEditFieldLabel.Text = 'filter Band';
    app.filterBandEditFieldLabel.Visible = false;

    % Create filterBandEditField
    app.filterBandEditField = uieditfield(app.CONNvariablePanel, 'text', 'Value',num2str(o.conn_filtBand),"FontSize",12*resizeF);
    app.filterBandEditField.Position = [391*ResizeX 167*ResizeY 121*ResizeX 22*ResizeY];
    app.filterBandEditField.Visible = false;

    % ----------------------------------------


    % ---------------------------------1st level Analysis
    % Create AnalysisNameEditFieldLabel
    app.AnalysisNameEditFieldLabel = uilabel(app.CONNvariablePanel,"FontSize",12*resizeF);
    app.AnalysisNameEditFieldLabel.Interpreter = 'latex';
    app.AnalysisNameEditFieldLabel.HorizontalAlignment = 'right';
    app.AnalysisNameEditFieldLabel.Position = [278*ResizeX 137*ResizeY 93*ResizeX 22*ResizeY];
    app.AnalysisNameEditFieldLabel.Text = 'Analysis Name';
    app.AnalysisNameEditFieldLabel.Visible = false;

    % Create AnalysisNameEditField
    app.AnalysisNameEditField = uieditfield(app.CONNvariablePanel, 'text','Value',o.conn_AnalysisName,'HorizontalAlignment','right',"FontSize",12*resizeF);
    app.AnalysisNameEditField.Position = [391*ResizeX 137*ResizeY 121*ResizeX 22*ResizeY];
    app.AnalysisNameEditField.Visible = false;

    % ---------------------------------------------------



    % -------------------------------2nd level
    % Create UITable
    tab = o.conn_contrast;
    if size(tab,1) > 1,tab = tab';end
    tab = cat(1,tab,repmat({[]},10,1));
    tab = table(tab);
    tab.Properties.VariableNames = {'contrast'};
    app.UITable = uitable(app.CONNvariablePanel,"Data",tab,"FontSize",12*resizeF);
    app.UITable.ColumnName = {'contrast'};
    app.UITable.RowName = {};
    app.UITable.Position = [551*ResizeX 15*ResizeY 98*ResizeY 123*ResizeX];
    app.UITable.Visible = false;
    app.UITable.ColumnEditable = true;


    % --------------------------------------------

    % ---------------------------------------setup
   
    % Create TREditFieldLabel
    app.TREditFieldLabel = uilabel(app.CONNvariablePanel,"FontSize",12*resizeF);
    app.TREditFieldLabel.Interpreter = 'latex';
    app.TREditFieldLabel.HorizontalAlignment = 'right';
    app.TREditFieldLabel.Position = [345*ResizeX 353*ResizeY 26*ResizeX 22*ResizeY];
    app.TREditFieldLabel.Text = 'TR';
    app.TREditFieldLabel.Visible = false;

    
    % Create TREditField
    app.TREditField = uieditfield(app.CONNvariablePanel,'numeric',"Value",o.conn_TR,"FontSize",12*resizeF);
    app.TREditField.Position = [391*ResizeX 353*ResizeY 121*ResizeX 22*ResizeY];
    app.TREditField.Visible = false;

    
    % Create sliceorderEditFieldLabel
    app.sliceorderEditFieldLabel = uilabel(app.CONNvariablePanel,"FontSize",12*resizeF);
    app.sliceorderEditFieldLabel.Interpreter = 'latex';
    app.sliceorderEditFieldLabel.HorizontalAlignment = 'right';
    app.sliceorderEditFieldLabel.Position = [305*ResizeX 291*ResizeY 66*ResizeX 22*ResizeY];
    app.sliceorderEditFieldLabel.Text = 'slice order';
    app.sliceorderEditFieldLabel.Visible = false;

    
    % Create sliceorderEditField
    app.sliceorderEditField = uieditfield(app.CONNvariablePanel, 'text',"Value",num2str(o.conn_slice_order),"FontSize",12*resizeF);
    app.sliceorderEditField.Position = [391*ResizeX 291*ResizeY 121*ResizeX 22*ResizeY];
    app.sliceorderEditField.Visible = false;

    % Create RestconditionListBoxLabel
    app.RestconditionListBoxLabel = uilabel(app.CONNvariablePanel,"FontSize",12*resizeF);
    app.RestconditionListBoxLabel.HorizontalAlignment = 'right';
    app.RestconditionListBoxLabel.Position = [569*ResizeX 247*ResizeY 53*ResizeX 22*ResizeY];
    app.RestconditionListBoxLabel.Text = 'condition';
    app.RestconditionListBoxLabel.Visible = false;

    % Create RestconditionListBox
    app.RestconditionListBox = uilistbox(app.CONNvariablePanel,"Items",{'REST'});
    app.RestconditionListBox.Multiselect = 'on';
    app.RestconditionListBox.Position = [549*ResizeX 158*ResizeY 100*ResizeX 90*ResizeY];
    app.RestconditionListBox.Value = {'REST'};
    app.RestconditionListBox.Visible = false;

    % Create TaskconditionListBoxLabel
    app.TaskconditionListBoxLabel = uilabel(app.CONNvariablePanel,"FontSize",12*resizeF);
    app.TaskconditionListBoxLabel.HorizontalAlignment = 'right';
    app.TaskconditionListBoxLabel.Position = [569*ResizeX 247*ResizeY 53*ResizeX 22*ResizeY];
    app.TaskconditionListBoxLabel.Text = 'condition';
    app.TaskconditionListBoxLabel.Visible = false;

    % Create TaskconditionListBox
    app.TaskconditionListBox = uilistbox(app.CONNvariablePanel,"Items",o.taskconditon,"FontSize",12*resizeF);
    app.TaskconditionListBox.Multiselect = 'on';
    app.TaskconditionListBox.Position = [549*ResizeX 158*ResizeY 100*ResizeX 90*ResizeY];
    if ~isempty(o.conn_conditon)
        app.TaskconditionListBox.Value = o.conn_conditon;
    else
        app.TaskconditionListBox.Value = o.taskconditon;
    end
    app.TaskconditionListBox.Visible = false;

    % Create roundListBoxLabel
    app.roundListBoxLabel = uilabel(app.CONNvariablePanel,"FontSize",12*resizeF);
    app.roundListBoxLabel.HorizontalAlignment = 'right';
    app.roundListBoxLabel.Position = [582*ResizeX 378*ResizeY 36*ResizeX 22*ResizeY];
    app.roundListBoxLabel.Text = 'round';
    app.roundListBoxLabel.Visible = false;

    % Create roundListBox
    app.roundListBox = uilistbox(app.CONNvariablePanel,"Items",o.round,"FontSize",12*resizeF);
    app.roundListBox.Multiselect = 'on';
    app.roundListBox.Position = [548*ResizeX 291*ResizeY 100*ResizeX 88*ResizeY];
    app.roundListBox.Value = o.conn_round;
    app.roundListBox.Visible = false;
    app.roundListBox.ValueChangedFcn = @(src,event,x,y) roundListValueChanged(src,event,app,o) ;
    if any(contains(o.conn_round,'Rest') | contains(o.conn_round,'REST') | contains(o.conn_round,'rest'))
        app.roundListBox.UserData = 'Rest';
    else
        app.roundListBox.UserData = 'Task';
    end

    % Create sessionEditFieldLabel
    app.sessionEditFieldLabel = uilabel(app.CONNvariablePanel,"FontSize",12*resizeF);
    app.sessionEditFieldLabel.Interpreter = 'latex';
    app.sessionEditFieldLabel.HorizontalAlignment = 'right';
    app.sessionEditFieldLabel.Position = [341*ResizeX 322*ResizeY 30*ResizeX 22*ResizeY];
    app.sessionEditFieldLabel.Text = 'sess';
    app.sessionEditFieldLabel.Visible = false;

    % Create sessionEditField
    app.sessionEditField = uieditfield(app.CONNvariablePanel, 'text','Value',o.conn_ses,'HorizontalAlignment','right',"FontSize",12*resizeF);
    app.sessionEditField.Position = [391*ResizeX 322*ResizeY 121*ResizeX 22*ResizeY];
    app.sessionEditField.Visible = false;

    switch class(o.conn_steps)
    case "char" 
        step = steps(o.conn_steps);
    case "string"
        step = steps(o.conn_steps);
    case "struct"
        step = o.conn_steps;
    end
    % Create Tree2
    app.Tree2 = uitree(app.CONNvariablePanel, 'checkbox',"FontSize",12*resizeF);
    app.Tree2.Position = [12*ResizeX 15*ResizeY 155*ResizeX 205*ResizeY];
    app.Tree2.CheckedNodesChangedFcn = @(src,event,x) Tree2checkchange(src,event,app);
    
    stepCheck = [];

    % Create SetupNodeconditionListBox
    app.SetupNode = uitreenode(app.Tree2);
    app.SetupNode.Text = 'Setup';
    if step.Setup
        stepCheck = cat(2,stepCheck,app.SetupNode); 
        app.TREditFieldLabel.Visible = true;
        app.TREditField.Visible = true;
        app.sliceorderEditFieldLabel.Visible = true;
        app.sliceorderEditField.Visible = true;
        app.roundListBoxLabel.Visible = true;
        app.roundListBox.Visible = true;
        if app.roundListBox.UserData == "Rest"
            app.RestconditionListBoxLabel.Visible = true;
            app.RestconditionListBox.Visible = true;
            app.TaskconditionListBoxLabel.Visible = false;
            app.TaskconditionListBox.Visible = false;
        else
            app.RestconditionListBoxLabel.Visible = false;
            app.RestconditionListBox.Visible = false;
            app.TaskconditionListBoxLabel.Visible = true;
            app.TaskconditionListBox.Visible = true;
        end
        app.sessionEditFieldLabel.Visible = true;
        app.sessionEditField.Visible = true;
    end
    
    % Create PreprocessingNode
    app.PreprocessingNode = uitreenode(app.Tree2);
    app.PreprocessingNode.Text = 'Preprocessing';
    if step.Preprocessing
        stepCheck = cat(2,stepCheck,app.PreprocessingNode); 
        app.FunctionalvoxelresolutionEditField.Visible = true;
        app.FunctionalvoxelresolutionEditFieldLabel.Visible = true;
        app.StructurevoxelresolutionEditField.Visible = true;
        app.StructurevoxelresolutionEditFieldLabel.Visible = true;
        app.smoothkernelEditField.Visible = true;
        app.smoothkernelEditFieldLabel.Visible = true;
    end
    
    % Create DenoisingNode
    app.DenoisingNode = uitreenode(app.Tree2);
    app.DenoisingNode.Text = 'Denoising';
    if step.Denoising 
        stepCheck = cat(2,stepCheck,app.DenoisingNode); 
        app.filterBandEditFieldLabel.Visible = true;
        app.filterBandEditField.Visible = true;
    end
    
    % Create stLevelNode
    app.stLevelNode = uitreenode(app.Tree2);
    app.stLevelNode.Text = '1st Level';
    if step.fst_Analysis
        stepCheck = cat(2,stepCheck,app.stLevelNode); 
        app.AnalysisNameEditFieldLabel.Visible = true;
        app.AnalysisNameEditField.Visible = true;
    end
    
    % Create ndLevelNode
    app.ndLevelNode = uitreenode(app.Tree2);
    app.ndLevelNode.Text = '2nd Level';
    if step.snd_Analysis
        stepCheck = cat(2,stepCheck,app.ndLevelNode); 
        app.UITable.Visible = true;
    end
    
    % Create AddROINode
    app.AddROINode = uitreenode(app.Tree2);
    app.AddROINode.Text = 'Add ROI';
    if step.Add_Roi
        stepCheck = cat(2,stepCheck,app.AddROINode);
        app.userdefineROIpathEditFieldLabel.Visible = true;
        app.userdefineROIpathEditField.Visible = true;
    end
    
    app.Tree2.CheckedNodes = stepCheck;

    % Create newprojectCheckBox
    app.newprojectCheckBox = uicheckbox(app.CONNvariablePanel,"FontSize",12*ResizeY);
    app.newprojectCheckBox.Text = 'raw data';
    app.newprojectCheckBox.Position = [432*ResizeX 74*ResizeY 83*ResizeX 22*ResizeY];
    app.newprojectCheckBox.Value = o.conn_rawdata;
    
    % ------------------------------------------------------------
    % -------------------------------------------------------------

    % --------------------------------------------------------ISC 
    % Create ISCvariablePanel
    app.ISCvariablePanel = uipanel(app.UIFigure);
    app.ISCvariablePanel.BorderColor = [0.302 0.7451 0.9333];
    app.ISCvariablePanel.BorderWidth = 3;
    app.ISCvariablePanel.Title = 'ISC variable';
    app.ISCvariablePanel.BackgroundColor = [1 1 1];
    app.ISCvariablePanel.FontAngle = 'italic';
    app.ISCvariablePanel.FontWeight = 'bold';
    app.ISCvariablePanel.FontSize = 18*resizeF;
    app.ISCvariablePanel.Position = [711*ResizeX 649*ResizeY 679*ResizeX 443*ResizeY];

    % Create inputpathEditFieldLabel
    app.inputpathEditFieldLabel = uilabel(app.ISCvariablePanel,"FontSize",12*resizeF);
    app.inputpathEditFieldLabel.Interpreter = 'latex';
    app.inputpathEditFieldLabel.HorizontalAlignment = 'right';
    app.inputpathEditFieldLabel.Position = [27*ResizeX 374*ResizeY 64*ResizeX 22*ResizeY];
    app.inputpathEditFieldLabel.Text = 'input path';

    % Create inputpathEditField
    app.inputpathEditField = uieditfield(app.ISCvariablePanel, 'text','Value', o.ISCinpath,"FontSize",12*resizeF);
    app.inputpathEditField.Position = [111*ResizeX 374*ResizeY 543*ResizeX 22*ResizeY];

    % Create outputpathEditFieldLabel
    app.outputpathEditFieldLabel = uilabel(app.ISCvariablePanel,"FontSize",12*resizeF);
    app.outputpathEditFieldLabel.Interpreter = 'latex';
    app.outputpathEditFieldLabel.HorizontalAlignment = 'right';
    app.outputpathEditFieldLabel.Position = [20*ResizeX 342*ResizeY 71*ResizeX 22*ResizeY];
    app.outputpathEditFieldLabel.Text = 'output path';

    % Create outputpathEditField
    app.outputpathEditField = uieditfield(app.ISCvariablePanel, 'text','Value', o.ISCoutpath,"FontSize",12*resizeF);
    app.outputpathEditField.Position = [111*ResizeX 342*ResizeY 543*ResizeX 22*ResizeY];

    % Create UITable
    app.UITable2 = uitable(app.ISCvariablePanel,'ColumnEditable',true,"FontSize",12*resizeF);
    if size(o.ISCcondition,1) == 1, o.ISCcondition = o.ISCcondition'; end
    group = repmat(o.ISCgroup,length(o.ISCcondition),1);
    cn = [];
    for i = 1:size(group,2)
        cn = cat(1,cn,strcat(o.ISCcondition,'_',group(:,i)));
    end
    app.UITable2.ColumnName = cn;
    L = unique(cellfun(@length,o.ISCcon));
    if length(L)==1 && L == length(cn)
        tab = array2table(num2cell(cell2mat(o.ISCcon)));
        tab = cat(1,tab,cell(10,size(tab,2)));
        app.UITable2.Data = tab;
    else
        tab = cell(length(o.ISCcon)+10,length(cn));
        tab = cell2table(tab);
        app.UITable2.Data = tab;
    end
    app.UITable2.RowName = {};
    app.UITable2.Position = [28*ResizeX 48*ResizeY 625*ResizeX 122*ResizeY];

    % Create conditionListBoxLabel
    app.conditionListBoxLabel = uilabel(app.ISCvariablePanel,"FontSize",12*resizeF);
    app.conditionListBoxLabel.HorizontalAlignment = 'right';
    app.conditionListBoxLabel.Position = [442*ResizeX 308*ResizeY 53*ResizeX 22*ResizeY];
    app.conditionListBoxLabel.Text = 'condition';

    % Create conditionListBox
    item = o.round;
    item = item(~contains(item,'T1'));
    app.conditionListBox = uilistbox(app.ISCvariablePanel,'Items',item,"FontSize",12*resizeF);
    app.conditionListBox.Multiselect = 'on';
    app.conditionListBox.Position = [422*ResizeX 183*ResizeY 100*ResizeX 126*ResizeY];
    app.conditionListBox.Value = o.ISCcondition;
    

    % Create conditionListBoxLabel
    app.groupListBoxLabel = uilabel(app.ISCvariablePanel,"FontSize",12*resizeF);
    app.groupListBoxLabel.HorizontalAlignment = 'right';
    app.groupListBoxLabel.Position = [175*ResizeX 308*ResizeY 53*ResizeX 22*ResizeY];
    app.groupListBoxLabel.Text = 'group';

    % Create conditionListBox
    app.groupListBox = uilistbox(app.ISCvariablePanel,'Items',o.ISCgroup,"FontSize",12*resizeF);
    app.groupListBox.Multiselect = 'on';
    app.groupListBox.Position = [155*ResizeX 183*ResizeY 100*ResizeX 126*ResizeY];
    app.groupListBox.Value =  o.ISCgroup;

    app.conditionListBox.ValueChangedFcn = @(src,evt,x) ISCconValueChanged(src,evt,app);
    app.groupListBox.ValueChangedFcn = @(src,evt,x) ISCconValueChanged(src,evt,app);
    % ------------------------------------------------------------

    % Create DoneButton
    app.DoneButton = uibutton(app.preprocessfiledefinePanel, 'push','ButtonPushedFcn',@(scr,event) DonebutomPush(scr,event,app,local_flag),"FontSize",12*resizeF);
    app.DoneButton.BackgroundColor = [1 1 1];
    app.DoneButton.FontWeight = 'bold';
    app.DoneButton.FontAngle = 'italic';
    app.DoneButton.Position = [577*ResizeX 521*ResizeY 88*ResizeX 23*ResizeY];
    app.DoneButton.Text = 'Done';

    % Show the figure after all components are created
    app.UIFigure.Visible = 'on';
    if ~connF
        app.CONNvariablePanel.Visible = false;
    end
    if ~ISCF
        app.UIFigure.Position = [100*ResizeX 100*ResizeY 700*ResizeX 1100*ResizeY];
        app.ISCvariablePanel.Visible = false;
    end
    
    while 1
        o = app.UIFigure.UserData;
        if o.done
            break;
        else
            app.UIFigure.UserData = [];
            waitfor(app.UIFigure,'UserData');
        end
    end
    delete(app.UIFigure);
end

function DonebutomPush(~,~,app,local_flag)
    Ef = 0;
    % define folder nest
    fNest = app.NestfolderEditField.Value;
    fNest = split(fNest,'/');
    fNest = split(fNest,'\');
    if size(fNest,1) > 1, fNest = fNest'; end

    app.UIFigure.UserData.folder_nest = fNest;
    % define path and nas info
    if local_flag
        app.UIFigure.UserData.subpath = app.localsubjectpath.Value;
        app.UIFigure.UserData.outsubpath = app.localoutputpath.Value;
        app.UIFigure.UserData.ftpServer = struct();
    else 
        app.UIFigure.UserData.subpath = [];
        app.UIFigure.UserData.ftpServer.ip = app.ftpnasip.Value;
        app.UIFigure.UserData.ftpServer.account = app.ftpnasaccount.Value;
        app.UIFigure.UserData.ftpServer.password = app.ftpnaspassword.Value;
        app.UIFigure.UserData.ftpServer.infolder = app.ftpnasinputfolder.Value;
        app.UIFigure.UserData.ftpServer.outfolder = app.ftpnasoutputfolder.Value;
        app.UIFigure.UserData.outsubpath = app.localoutputpath.Value;
    end
    % define subject, session, round
    tmp = app.subjectTable.Data.subject;
    tmp = getuiTableData(tmp);
    app.UIFigure.UserData.sub = tmp;

    tmp = app.sessionTable.Data.sess;
    tmp = getuiTableData(tmp);
    app.UIFigure.UserData.sess = tmp;

    tmp = app.roundTable.Data.round;
    tmp = getuiTableData(tmp);
    app.UIFigure.UserData.round = tmp;

    app.UIFigure.UserData.conn_prjName = app.projectNameEditField.Value;
    app.UIFigure.UserData.conn_prjPath = app.projectpathEditField.Value;
    app.UIFigure.UserData.conn_TR = app.TREditField.Value;
    
    tmp = app.sliceorderEditField.Value;

    if any(isstrprop(tmp,'digit'))
        if contains(tmp,',')
            tmp(tmp == ' ') = [];
            tmp = split(tmp,', ');
        elseif contains(tmp,' ')
            tmp = split(tmp,' ');
        end
        tmp(cellfun(@isempty,tmp)) = [];
        tmp = cellfun(@str2double,tmp);
        if size(tmp,1) > 1,tmp = tmp'; end
    end 
    
    app.UIFigure.UserData.conn_slice_order = tmp;

    app.UIFigure.UserData.conn_fwhm = app.smoothkernelEditField.Value;
    app.UIFigure.UserData.conn_strVres = app.StructurevoxelresolutionEditField.Value;
    app.UIFigure.UserData.conn_funcVres = app.FunctionalvoxelresolutionEditField.Value;
    tmp = app.filterBandEditField.Value;
    if contains(tmp,',')
        tmp(tmp == ' ') = [];
        tmp = split(tmp,',');
    elseif contains(tmp,' ')
        tmp = split(tmp,' ');
    end
    tmp(cellfun(@isempty,tmp)) = [];
    tmp = cellfun(@str2double,tmp);
    if size(tmp,1) > 1,tmp = tmp'; end
    app.UIFigure.UserData.conn_filtBand = tmp;

    app.UIFigure.UserData.conn_AnalysisName = app.AnalysisNameEditField.Value;

    if app.roundListBox.UserData == "Rest"
        app.UIFigure.UserData.conn_conditon = app.RestconditionListBox.Value;
    else
        app.UIFigure.UserData.conn_conditon = app.TaskconditionListBox.Value;
    end
    
    tmp = app.UITable.Data.contrast;
    tmp(cellfun(@isempty,tmp)) = [];
    app.UIFigure.UserData.conn_contrast = tmp;

    app.UIFigure.UserData.conn_ses = app.sessionEditField.Value;

    app.UIFigure.UserData.conn_round = app.roundListBox.Value;

    app.UIFigure.UserData.conn_ROIpath = app.userdefineROIpathEditField.Value;

    conn_steps = struct("Setup",0, ...
                       "Preprocessing",0, ...
                       "Denoising",0, ...
                       "fst_Analysis",0, ...
                       "snd_Analysis",0, ...
                       "Add_Roi",0);
    checkedNode = {app.Tree2.CheckedNodes.Text};

    app.UIFigure.UserData.conn_steps = app.Tree2.CheckedNodes;
    if any(contains(checkedNode,"Setup")), conn_steps.Setup = 1;end
    if any(contains(checkedNode,"Preprocessing")), conn_steps.Preprocessing = 1;end
    if any(contains(checkedNode,"Denoising")), conn_steps.Denoising = 1;end
    if any(contains(checkedNode,"1st Level")), conn_steps.fst_Analysis = 1;end
    if any(contains(checkedNode,"2nd Level")), conn_steps.snd_Analysis = 1;end
    if any(contains(checkedNode,"Add ROI")), conn_steps.Add_Roi = 1;end
    app.UIFigure.UserData.conn_steps = conn_steps;
    app.UIFigure.UserData.conn_rawdata = app.newprojectCheckBox.Value;

    tmp = app.inputpathEditField.Value;
    group = app.groupListBox.Value;
    for i = 1:length(group)
        if contains(tmp,group{i})
            tmp = erase(tmp,group{i});
        end
    end
    
    app.UIFigure.UserData.ISCinpath = tmp;
    app.UIFigure.UserData.ISCoutpath = app.outputpathEditField.Value;
    app.UIFigure.UserData.ISCcondition = app.conditionListBox.Value;
    app.UIFigure.UserData.ISCgroup = app.groupListBox.Value;

    tmp = cell2mat(table2array(app.UITable2.Data));
    app.UIFigure.UserData.ISCcon = mat2cell(tmp,ones(1,size(tmp,1)),size(tmp,2));

    if isempty(app.UIFigure.UserData.conn_conditon)
        warndlg('need select more than one condition');
        Ef = true;
    end

    if ~Ef
        app.UIFigure.UserData.done = true;
    else
        app.UIFigure.UserData.done = false;
    end

    function data = getuiTableData(data)  
        data(cellfun(@(x) x=="",data)) = [];
        idx = find(cellfun(@(x) class(x) == "string",data));
        for I = 1:length(idx)
            data{idx(I)} = char(data{idx(I)});
        end
    end
end

function app = createFolderNest(app,nest)
    % Create Tree for folder nest
    app.Tree = uitree(app.Panel);
    app.Tree.Position = [8 10 250 131];


    ses = app.sessionTable.Value;
    round = app.roundTable.Value;
    % find session and round
    for i = 1:length(nest)
        d = string(nest{i});
        if any(contains(ses,d)), nest{i} = 'sess'; end
        if any(contains(round,d)), nest{i} = 'round'; end
    end
    for i = 1:length(nest)
        if i == 1, parent = 'Tree';else, parent = ['Node',num2str(i-1)]; end
        tmp = fieldnames(app);
        if any(contains(tmp,['Node',num2str(i)]))
            app.(['Node',num2str(i)]).delete;
        end
        app.(['Node',num2str(i)]) = uitreenode(app.(parent));
        app.(['Node',num2str(i)]).Text = nest{i};
        app.(['Node',num2str(i)]).Icon = fullfile(pwd,'icon','folder.jpg');
        expand(app.(parent))
    end
end


% steps checked box change function 
function Tree2checkchange(src,~,app)
    if ~isempty(src.CheckedNodes)
        steps = {src.CheckedNodes.Text};
    else
        steps = {' '};
    end
    if any(contains(steps,"Setup"))
        app.TREditFieldLabel.Visible = true;
        app.TREditField.Visible = true;
        app.sliceorderEditFieldLabel.Visible = true;
        app.sliceorderEditField.Visible = true;
        app.roundListBoxLabel.Visible = true;
        app.roundListBox.Visible = true;
        if app.roundListBox.UserData == "Rest"
            app.RestconditionListBoxLabel.Visible = true;
            app.RestconditionListBox.Visible = true;
            app.TaskconditionListBoxLabel.Visible = false;
            app.TaskconditionListBox.Visible = false;
        else
            app.TaskconditionListBoxLabel.Visible = true;
            app.TaskconditionListBox.Visible = true;
            app.RestconditionListBoxLabel.Visible = false;
            app.RestconditionListBox.Visible = false;
        end
        app.sessionEditFieldLabel.Visible = true;
        app.sessionEditField.Visible = true;
    else
        app.TREditFieldLabel.Visible = false;
        app.TREditField.Visible = false;
        app.sliceorderEditFieldLabel.Visible = false;
        app.sliceorderEditField.Visible = false;
        app.roundListBoxLabel.Visible = false;
        app.roundListBox.Visible = false;
        app.RestconditionListBoxLabel.Visible = false;
        app.RestconditionListBox.Visible = false;
        app.TaskconditionListBoxLabel.Visible = false;
        app.TaskconditionListBox.Visible = false;
        app.sessionEditFieldLabel.Visible = false;
        app.sessionEditField.Visible = false;
    end
    if any(contains(steps,"Preprocessing"))
        app.FunctionalvoxelresolutionEditField.Visible = true;
        app.FunctionalvoxelresolutionEditFieldLabel.Visible = true;
        app.StructurevoxelresolutionEditField.Visible = true;
        app.StructurevoxelresolutionEditFieldLabel.Visible = true;
        app.smoothkernelEditField.Visible = true;
        app.smoothkernelEditFieldLabel.Visible = true;
    else
        app.FunctionalvoxelresolutionEditField.Visible = false;
        app.FunctionalvoxelresolutionEditFieldLabel.Visible = false;
        app.StructurevoxelresolutionEditField.Visible = false;
        app.StructurevoxelresolutionEditFieldLabel.Visible = false;
        app.smoothkernelEditField.Visible = false;
        app.smoothkernelEditFieldLabel.Visible = false;
    end
    if any(contains(steps,"Denoising"))
        app.filterBandEditFieldLabel.Visible = true;
        app.filterBandEditField.Visible = true;
    else
        app.filterBandEditFieldLabel.Visible = false;
        app.filterBandEditField.Visible = false;
    end
    if any(contains(steps,"1st Level"))
        app.AnalysisNameEditFieldLabel.Visible = true;
        app.AnalysisNameEditField.Visible = true;
    else
        app.AnalysisNameEditFieldLabel.Visible = false;
        app.AnalysisNameEditField.Visible = false;
    end
    if any(contains(steps,"2nd Level"))
        app.UITable.Visible = true;
    else
        app.UITable.Visible = false;
    end
    if any(contains(steps,"Add ROI"))
        app.userdefineROIpathEditFieldLabel.Visible = true;
        app.userdefineROIpathEditField.Visible = true;
    else
        app.userdefineROIpathEditFieldLabel.Visible = false;
        app.userdefineROIpathEditField.Visible = false;
    end
end
function roundListValueChanged(src,~,app,o)
    selectValue = src.Value;
    if any(contains(selectValue,'REST') | contains(selectValue,'Rest') | contains(selectValue,'rest'))
        src.UserData = 'Rest';
        app.RestconditionListBoxLabel.Visible = true;
        app.RestconditionListBox.Visible = true;
        app.TaskconditionListBoxLabel.Visible = false;
        app.TaskconditionListBox.Visible = false;
    else
        src.UserData = 'Task';
        app.TaskconditionListBoxLabel.Visible = true;
        app.TaskconditionListBox.Visible = true;
        app.RestconditionListBoxLabel.Visible = false;
        app.RestconditionListBox.Visible = false;
    end
    app.projectNameEditField.Value = ['conn_',cell2mat(selectValue)];
    app.projectpathEditField.Value = fullfile(o.conn_prjPath,['conn_',cell2mat(selectValue)]);
    
end

function ISCconValueChanged(src,~,app)
    selectValue = src.Value;
    oldData = app.UITable2.Data;
    oldcn = app.UITable2.ColumnName;
    group = app.groupListBox.Value;
    condition = app.conditionListBox.Value;
    group = repmat(group,length(condition),1);
    if size(condition,1)==1,condition = condition';end
    cn = [];
    for i = 1:size(group,2)
        cn = cat(1,cn,strcat(condition,'_',group(:,i)));
    end
    chk_i = cellfun(@(x) find(string(x)==string(oldcn)),cn,'UniformOutput',false);
    chk_i = cell2mat(chk_i(~cellfun(@isempty,chk_i)));
    chk = length(chk_i)==length(cn) && length(cn) < size(oldData,2);
    if chk
       app.UITable2.Data = oldData(:,chk_i); 
    else
        app.UITable2.ColumnName = cn;
        tmp = cell(3^length(cn),length(cn));
        app.UITable2.Data = cell2table(tmp);
    end
end