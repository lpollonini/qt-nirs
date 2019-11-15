classdef bpGoodQuality < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure              matlab.ui.Figure
        barPlot               matlab.ui.control.UIAxes
        ThresholdSliderLabel  matlab.ui.control.Label
        ThresholdSlider       matlab.ui.control.Slider
        ExporttonirsButton    matlab.ui.control.Button
    end

    
    properties (Access = private)
        qMats % Quality matrices from NIRSPlot quality computation
        qThld % Threshold for marking good-quality channels
        bp %boxplot
        thldLn % threshold line
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, qualityMats, qualityThld)
            app.qMats = qualityMats;
            app.qThld = qualityThld;            
            app.bp = bar(app.barPlot,app.qMats.gcl,'Horizontal','on');
            app.bp.FaceColor = 'flat';
            app.barPlot.Box = 'on';
            app.barPlot.YLabel.String  = 'Channel';
            app.barPlot.XLim = [0,1];
            app.thldLn = xline(app.barPlot,app.qThld,'--r');
            app.ThresholdSlider.Value = round(app.qThld*100);
            idxThl = qualityMats.gcl<qualityThld;
            app.bp.CData(idxThl,:) = repmat([0.6, 0.6, 0.6],sum(idxThl),1);
            app.bp.CData(~idxThl,:) = repmat([0 1 0],sum(~idxThl),1);
            
        end

        % Value changed function: ThresholdSlider
        function ThresholdSliderValueChanged(app, event)
%             app.qThld = (app.ThresholdSlider.Value)/100;
%             b = bar(app.barPlot,app.qMats.gcl,'Horizontal','on');
%             idxThl = app.qMats.gcl<app.qThld;
%             b.CData(idxThl,:) = repmat([0.6, 0.6, 0.6],sum(idxThl),1);
        end

        % Value changing function: ThresholdSlider
        function ThresholdSliderValueChanging(app, event)
            app.qThld  = (event.Value)/100;
            app.thldLn.Value = app.qThld;
            %bar(app.barPlot,app.qMats.gcl,'Horizontal','on');
            %app.bp = bar(app.barPlot,app.qMats.gcl,'Horizontal','on');
            idxThl = app.qMats.gcl<app.qThld;
            app.bp.CData(idxThl,:) = repmat([0.6, 0.6, 0.6],sum(idxThl),1);
            app.bp.CData(~idxThl,:) = repmat([0 1 0],sum(~idxThl),1);
            %xline(app.barPlot,app.qThld,'--r');           
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 640 480];
            app.UIFigure.Name = 'UI Figure';

            % Create barPlot
            app.barPlot = uiaxes(app.UIFigure);
            title(app.barPlot, 'Good-quality Channels')
            xlabel(app.barPlot, 'X')
            ylabel(app.barPlot, 'Y')
            app.barPlot.Position = [57 112 402 319];

            % Create ThresholdSliderLabel
            app.ThresholdSliderLabel = uilabel(app.UIFigure);
            app.ThresholdSliderLabel.HorizontalAlignment = 'right';
            app.ThresholdSliderLabel.Position = [486 112 59 22];
            app.ThresholdSliderLabel.Text = 'Threshold';

            % Create ThresholdSlider
            app.ThresholdSlider = uislider(app.UIFigure);
            app.ThresholdSlider.Orientation = 'vertical';
            app.ThresholdSlider.ValueChangedFcn = createCallbackFcn(app, @ThresholdSliderValueChanged, true);
            app.ThresholdSlider.ValueChangingFcn = createCallbackFcn(app, @ThresholdSliderValueChanging, true);
            app.ThresholdSlider.Position = [566 121 3 303];

            % Create ExporttonirsButton
            app.ExporttonirsButton = uibutton(app.UIFigure, 'push');
            app.ExporttonirsButton.Position = [509 57 100 22];
            app.ExporttonirsButton.Text = 'Export to .nirs';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = bpGoodQuality(varargin)

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @(app)startupFcn(app, varargin{:}))

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end