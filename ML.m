classdef ML < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                       matlab.ui.Figure
        AMPIEZZACICLOLIMITESTABILELabel  matlab.ui.control.Label
        AMPIEZZACICLOLIMITEINSTABILELabel  matlab.ui.control.Label
        AMPIEZZACICLOLIMITESTABILELabel_2  matlab.ui.control.Label
        EQUILIBRIINSTABILILabel        matlab.ui.control.Label
        EQUILIBRISTABILILabel          matlab.ui.control.Label
        bulletEQUILIBRIINSTABILILabel  matlab.ui.control.Label
        bulletEQUILIBRISTABILILabel    matlab.ui.control.Label
        corrente                       matlab.ui.control.NumericEditField
        CORRENTEAPPLICATAEditFieldLabel  matlab.ui.control.Label
        bf                             matlab.ui.control.UIAxes
        rf                             matlab.ui.control.UIAxes
    end

    

    % Callbacks that handle component events
    methods (Access = private)

        % Value changed function: corrente
        function correnteValueChanged(app, event)
            
            cla(app.rf)
            cla(app.bf)
            
            value = app.corrente.Value;
            gca = 1;
            gk = 2;
            gl = 0.5;
            Eca = 1;
            Ek = -0.7;
            El = -0.5;
            phi = 1/3;
            Iapp = value;
            V1 = -0.01;
            V2 = 0.15;
            V3 = 0.1;
            V4 = 0.145;
            
            minf = @(V) 1/2 * (1 + tanh((V-V1)/(V2)));
            winf = @(V) 1/2 * (1 + tanh((V-V3)/(V4)));
            tau = @(V) (cosh((V-V3)/(2*V4)))^(-1);
            Ica = @(V) gca*minf(V).*(V-Eca);
            Il = @(V) gl*(V-El);
            
            f = @(t,Y) [-gca*minf(Y(1))*(Y(1)-Eca)-gk*Y(2)*(Y(1)-Ek)-gl*(Y(1)-El)+Iapp;...
                phi*(winf(Y(1))-Y(2))/tau(Y(1))];
            
            y1 = linspace(-0.6,0.3,20);
            y2 = linspace(-0.1,0.6,20);
            [x,y] = meshgrid(y1,y2);
            
            u = zeros(size(x));
            v = zeros(size(x));
            
            t=0; % we want the derivatives at each point at t=0, i.e. the starting time
            for i = 1:numel(x)
                Yprime = f(t,[x(i); y(i)]);
                u(i) = Yprime(1);
                v(i) = Yprime(2);
            end
            u = u./(sqrt(u.^2+v.^2));
            v = v./(sqrt(u.^2+v.^2));
            W_nullocline = winf(y1);
            V_nullocline = (Iapp-Ica(y1)-Il(y1))./(gk.*(y1-Ek));
            
            %% fixed point with Newton
            f = @(Y) [-gca*minf(Y(1))*(Y(1)-Eca)-gk*Y(2)*(Y(1)-Ek)-gl*(Y(1)-El)+Iapp;...
                phi*(winf(Y(1))-Y(2))/tau(Y(1))];
            
            h1 = @(x) tanh((x-V1)./(V2));
            h2 = @(x) (1./V2) * (1-(tanh((x-V1)./(V2))).^2);
            h3 = @(x) tanh((x-V3)./(V4));
            h4 = @(x) (1./V4) * (1-(tanh((x-V3)./(V4))).^2);
            
            c = @(x) cosh( (x-V3)./(2*V4) );
            s = @(x) (1./(2*V4)) * sinh( (x-V3)./(2*V4) );
            
            J1 = @(Y) -gca*(1/2*h2(Y(1)))*(Y(1)-Eca)-gca*(1/2+1/2*h1(Y(1)))-gk*Y(2)-gl;
            J2 = @(Y) -gk*(Y(1)-Ek);
            J3 = @(Y) phi.*(1/2*h4(Y(1))*c(Y(1))+(1/2+1/2*h3(Y(1))-Y(2))*s(Y(1)));
            J4 = @(Y) -phi./tau(Y(1));
            Jf = @(Y) [J1(Y),J2(Y);J3(Y),J4(Y)];
            
            
            tol = 1e-10;
            y0 = [-0.3535;0.002];
            e1 = y0;
            delta = -Jf(e1)\f(e1);
            while norm(delta,inf)>tol
                e1 = e1+delta;
                delta = -Jf(e1)\f(e1);
            end 
            e1 = e1+delta;
            eig_e1 = eig(Jf(e1));
            
            y0 = [-0.16;0.02];
            e2 = y0;
            delta = -Jf(e2)\f(e2);
            while norm(delta,inf)>tol
                e2 = e2+delta;
                delta = -Jf(e2)\f(e2);
            end 
            e2 = e2+delta;
            eig_e2 = eig(Jf(e2));
            
            y0 = [0.03;0.27];
            e3 = y0;
            delta = -Jf(e3)\f(e3);
            while norm(delta,inf)>tol
                e3 = e3+delta;
                delta = -Jf(e3)\f(e3);
            end 
            e3 = e3+delta;
            eig_e3 = eig(Jf(e3));
            
            hold(app.rf, 'on')
            quiver(app.rf, x,y,u,v,0.3, 'color', [0.6 0.6 0.6]);
            plot(app.rf, y1,V_nullocline, 'color', [0.4940 0.1840 0.5560]);
            plot(app.rf,y1,W_nullocline, 'Color', [0.9290 0.6940 0.1250]);
            if real(eig_e1(1)) < 0 && real(eig_e1(2)) < 0
                scatter(app.rf,e1(1),e1(2),'r', 'filled')
            end
            if real(eig_e1(1)) > 0 || real(eig_e1(2)) > 0
                scatter(app.rf,e1(1),e1(2),'k', 'filled')
            end
            
            if real(eig_e2(1)) < 0 && real(eig_e2(2)) < 0
                scatter(app.rf,e2(1),e2(2),'r', 'filled')
            end
            if real(eig_e2(1)) > 0 || real(eig_e2(2)) > 0
                scatter(app.rf,e2(1),e2(2),'k', 'filled')
            end
                        
            if real(eig_e3(1)) < 0 && real(eig_e3(2)) < 0
                scatter(app.rf,e3(1),e3(2),'r', 'filled')
            end
            if real(eig_e3(1)) > 0 || real(eig_e3(2)) > 0
                scatter(app.rf,e3(1),e3(2),'k', 'filled')
            end
            
            % ciclo limite 
            if value > 0.0832 && value < 0.2412
                f = @(t,Y) [-gca*minf(Y(1))*(Y(1)-Eca)-gk*Y(2)*(Y(1)-Ek)-gl*(Y(1)-El)+Iapp;...
                phi*(winf(Y(1))-Y(2))/tau(Y(1))];
                
                y0 = [0.005;0.23];
                [~,yout] = ode45(f,[0,200],y0);
                plot(app.rf, yout(end-100:end,1),yout(end-100:end,2), 'g')
            end
            hold(app.rf, 'off')
            
            
            
            
            %biforcazione
            
            A = importdata('diagram.dat');
            n = length(A(:,1));
            % 1 = eq stabili
            % 2 = eq instabili
            % 3 = ciclo limite stabile
            % 4 = ciclo limite instabile
            
            % stable
            b = zeros(1,n);
            counter = 1;
            for i = 1:n
                if A(i,4) == 1
                   b(counter) = i;
                   counter = counter + 1;
                end 
            end 
            b = b(1:36);
            
            stabile = zeros(36,2);
            for i = 1:36
                stabile(i,:) = [A(b(i),1),A(b(i),2)];
            end
            
            % instable
            c = zeros(1,n);
            counter = 1;
            for i = 1:n
                if A(i,4) == 2
                   c(counter) = i;
                   counter = counter + 1;
                end 
            end 
            c = c(1:36);
            
            instabile = zeros(36,2);
            for i = 1:36
                instabile(i,:) = [A(c(i),1),A(c(i),2)];
            end 
            instabile = instabile(2:end,:);
            
            
            % congiungo nera e rossa
            xi = [max(instabile(:,1)), min(stabile(19:end,1))];
            xf = [ max(instabile(end,2)), min(stabile(19:end,2))];
           
            
            
            % cl stabile
            d = zeros(1,n);
            counter = 1;
            for i = 1:n
                if A(i,4) == 3
                   d(counter) = i;
                   counter = counter + 1;
                end 
            end 
            d = d(1:43);
            
            cls = zeros(43,2);
            for i = 1:43
                cls(i,:) = [A(d(i),1),A(d(i),2)];
            end 
            
            % cl instabile
            g = zeros(1,n);
            counter = 1;
            for i = 1:n
                if A(i,4) == 4
                   g(counter) = i;
                   counter = counter + 1;
                end 
            end 
            g = g(1:18);
            
            cl = zeros(18,2);
            for i = 1:18
                cl(i,:) = [A(g(i),1),A(g(i),2)];
            end
            
            
            %%
            B = A(74:end,:);
            n1 = length(B(:,1));
            % cl stabile
            d1 = zeros(1,n1);
            counter = 1;
            for i = 1:n1
                if B(i,4) == 3
                   d1(counter) = i;
                   counter = counter + 1;
                end 
            end 
            d1 = d1(1:43);
            
            cls1 = zeros(43,2);
            for i = 1:43
                cls1(i,:) = [B(d1(i),1),B(d1(i),3)];
            end 
            
            
            % cl instabile
            g1 = zeros(1,n1);
            counter = 1;
            for i = 1:n1
                if B(i,4) == 4
                   g1(counter) = i;
                   counter = counter + 1;
                end 
            end 
            g1 = g1(1:17);
            
            cl1 = zeros(17,2);
            for i = 1:17
                cl1(i,:) = [B(g1(i),1),B(g1(i),3)];
            end 
            
            % congiungo blu e blu
            xi1 = [max(instabile(:,1)), min(cl1(1:11,1))];
            xf1 = [max(instabile(:,2)), max(cl1(1:11,2))];
            
            
            % congiungo blu e verde
            xi2 = [cls(1,1), max(cl(1:12,1))];
            xf2 = [cls(1,2), max(cl(1:12,2))];
            
            xi3 = [max(cls1(:,1)), max(cl1(1:11,1))];
            xf3 = [max(cls1(:,2)), min(cl1(1:11,2))];
            
            hold(app.bf,'on')
            plot(app.bf,stabile(1:18,1),stabile(1:18,2),'r')
            plot(app.bf,stabile(19:end,1),stabile(19:end,2),'r')
            plot(app.bf,instabile(:,1),instabile(:,2),'k')
            plot(app.bf,xi,xf,'r')
            plot(app.bf,cls(1:end-4,1),cls(1:end-4,2),'g')
            plot(app.bf,cl(1:12,1),cl(1:12,2),'b')
            plot(app.bf,cls1(1:end-4,1),cls1(1:end-4,2),'g')
            plot(app.bf,cl1(1:11,1),cl1(1:11,2),'b')
            plot(app.bf,xi1,xf1,'b')
            plot(app.bf,xi2,xf2,'b')
            plot(app.bf,xi3,xf3,'b')
            l = linspace(-0.6,0.3,100);
            l1 = ones(1,100)*0.2038;
            plot(app.bf, l1,l, 'color', [0.9 0.9 0.9])
            l2 = ones(1,100)*0.2420;
            plot(app.bf, l2,l, 'color', [0.9 0.9 0.9])
            l3 = ones(1,100)*0.0832;
            plot(app.bf, l3,l, 'color', [0.9 0.9 0.9])
            ly = linspace(-0.6,0.3,100);
            lx = ones(1,100)*value;
            plot(app.bf, lx,ly,'--', 'color', [0.3 0.7 0.9], 'LineWidth', 1.2)
            
            text(app.bf, 0.18, -0.56, 'subHB', 'color', [0.6 0.6 0.6])
            text(app.bf, 0.06, -0.56, 'superHB', 'color', [0.6 0.6 0.6])
            text(app.bf, 0.23, -0.56, 'SNP','color', [0.6 0.6 0.6])
            
            hold(app.bf,'off')
        end

        % Callback function
        function RUNButtonPushed(app, event)
            cla(app.rf)
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1112 626];
            app.UIFigure.Name = 'MATLAB App';

            % Create rf
            app.rf = uiaxes(app.UIFigure);
            title(app.rf, 'Ritratto in fase')
            xlabel(app.rf, 'V')
            ylabel(app.rf, 'W')
            zlabel(app.rf, 'Z')
            app.rf.XLim = [-0.6 0.3];
            app.rf.YLim = [-0.1 0.6];
            app.rf.XTick = [-0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3];
            app.rf.YTick = [-0.1 0 0.1 0.2 0.3 0.4 0.5 0.6];
            app.rf.YTickLabel = {'-0.1'; '0'; '0.1'; '0.2'; '0.3'; '0.4'; '0.5'; '0.6'};
            app.rf.HandleVisibility = 'off';
            app.rf.BusyAction = 'cancel';
            app.rf.PickableParts = 'none';
            app.rf.Position = [1 208 505 371];

            % Create bf
            app.bf = uiaxes(app.UIFigure);
            title(app.bf, 'Diagramma di biforcazione')
            xlabel(app.bf, 'I_{app}')
            ylabel(app.bf, 'V')
            zlabel(app.bf, 'Z')
            app.bf.XLim = [-0.05 0.4];
            app.bf.YLim = [-0.6 0.3];
            app.bf.XTick = [-0.2 -0.15 -0.1 -0.05 0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4];
            app.bf.YTick = [-0.6 -0.5 -0.4 -0.3 -0.2 -0.1 1.11022302462516e-16 0.1 0.2 0.3];
            app.bf.YTickLabel = {'-0.6'; '-0.5'; '-0.4'; '-0.3'; '-0.2'; '-0.1'; '0'; '0.1'; '0.2'; '0.3'};
            app.bf.Position = [566 208 500 371];

            % Create CORRENTEAPPLICATAEditFieldLabel
            app.CORRENTEAPPLICATAEditFieldLabel = uilabel(app.UIFigure);
            app.CORRENTEAPPLICATAEditFieldLabel.HorizontalAlignment = 'right';
            app.CORRENTEAPPLICATAEditFieldLabel.FontSize = 13;
            app.CORRENTEAPPLICATAEditFieldLabel.FontWeight = 'bold';
            app.CORRENTEAPPLICATAEditFieldLabel.Position = [476 27 78 32];
            app.CORRENTEAPPLICATAEditFieldLabel.Text = {'CORRENTE'; 'APPLICATA'; ''};

            % Create corrente
            app.corrente = uieditfield(app.UIFigure, 'numeric');
            app.corrente.Limits = [-0.02 0.4];
            app.corrente.ValueChangedFcn = createCallbackFcn(app, @correnteValueChanged, true);
            app.corrente.FontSize = 13;
            app.corrente.FontWeight = 'bold';
            app.corrente.Position = [566 32 37 22];

            % Create bulletEQUILIBRISTABILILabel
            app.bulletEQUILIBRISTABILILabel = uilabel(app.UIFigure);
            app.bulletEQUILIBRISTABILILabel.Interpreter = 'latex';
            app.bulletEQUILIBRISTABILILabel.BackgroundColor = [0.902 0.902 0.902];
            app.bulletEQUILIBRISTABILILabel.FontSize = 10;
            app.bulletEQUILIBRISTABILILabel.FontColor = [1 0 0];
            app.bulletEQUILIBRISTABILILabel.Position = [115 148 257 22];
            app.bulletEQUILIBRISTABILILabel.Text = '$\bullet$ EQUILIBRI STABILI';

            % Create bulletEQUILIBRIINSTABILILabel
            app.bulletEQUILIBRIINSTABILILabel = uilabel(app.UIFigure);
            app.bulletEQUILIBRIINSTABILILabel.Interpreter = 'latex';
            app.bulletEQUILIBRIINSTABILILabel.BackgroundColor = [0.902 0.902 0.902];
            app.bulletEQUILIBRIINSTABILILabel.FontSize = 10;
            app.bulletEQUILIBRIINSTABILILabel.Position = [115 127 257 22];
            app.bulletEQUILIBRIINSTABILILabel.Text = '$\bullet$ EQUILIBRI INSTABILI';

            % Create EQUILIBRISTABILILabel
            app.EQUILIBRISTABILILabel = uilabel(app.UIFigure);
            app.EQUILIBRISTABILILabel.Interpreter = 'latex';
            app.EQUILIBRISTABILILabel.BackgroundColor = [0.902 0.902 0.902];
            app.EQUILIBRISTABILILabel.FontSize = 10;
            app.EQUILIBRISTABILILabel.FontColor = [1 0 0];
            app.EQUILIBRISTABILILabel.Position = [757 148 263 22];
            app.EQUILIBRISTABILILabel.Text = '$-$ EQUILIBRI STABILI';

            % Create EQUILIBRIINSTABILILabel
            app.EQUILIBRIINSTABILILabel = uilabel(app.UIFigure);
            app.EQUILIBRIINSTABILILabel.Interpreter = 'latex';
            app.EQUILIBRIINSTABILILabel.BackgroundColor = [0.902 0.902 0.902];
            app.EQUILIBRIINSTABILILabel.FontSize = 10;
            app.EQUILIBRIINSTABILILabel.Position = [757 127 263 22];
            app.EQUILIBRIINSTABILILabel.Text = '$-$ EQUILIBRI INSTABILI';

            % Create AMPIEZZACICLOLIMITESTABILELabel_2
            app.AMPIEZZACICLOLIMITESTABILELabel_2 = uilabel(app.UIFigure);
            app.AMPIEZZACICLOLIMITESTABILELabel_2.Interpreter = 'latex';
            app.AMPIEZZACICLOLIMITESTABILELabel_2.BackgroundColor = [0.902 0.902 0.902];
            app.AMPIEZZACICLOLIMITESTABILELabel_2.FontSize = 10;
            app.AMPIEZZACICLOLIMITESTABILELabel_2.FontColor = [0 0.6784 0];
            app.AMPIEZZACICLOLIMITESTABILELabel_2.Position = [757 106 263 22];
            app.AMPIEZZACICLOLIMITESTABILELabel_2.Text = '$-$ AMPIEZZA CICLO LIMITE STABILE';

            % Create AMPIEZZACICLOLIMITEINSTABILELabel
            app.AMPIEZZACICLOLIMITEINSTABILELabel = uilabel(app.UIFigure);
            app.AMPIEZZACICLOLIMITEINSTABILELabel.Interpreter = 'latex';
            app.AMPIEZZACICLOLIMITEINSTABILELabel.BackgroundColor = [0.902 0.902 0.902];
            app.AMPIEZZACICLOLIMITEINSTABILELabel.FontSize = 10;
            app.AMPIEZZACICLOLIMITEINSTABILELabel.FontColor = [0 0 1];
            app.AMPIEZZACICLOLIMITEINSTABILELabel.Position = [757 85 263 22];
            app.AMPIEZZACICLOLIMITEINSTABILELabel.Text = '$-$ AMPIEZZA CICLO LIMITE INSTABILE';

            % Create AMPIEZZACICLOLIMITESTABILELabel
            app.AMPIEZZACICLOLIMITESTABILELabel = uilabel(app.UIFigure);
            app.AMPIEZZACICLOLIMITESTABILELabel.Interpreter = 'latex';
            app.AMPIEZZACICLOLIMITESTABILELabel.BackgroundColor = [0.902 0.902 0.902];
            app.AMPIEZZACICLOLIMITESTABILELabel.FontSize = 10;
            app.AMPIEZZACICLOLIMITESTABILELabel.FontColor = [0 0.6784 0];
            app.AMPIEZZACICLOLIMITESTABILELabel.Position = [115 106 257 22];
            app.AMPIEZZACICLOLIMITESTABILELabel.Text = '$-$ AMPIEZZA CICLO LIMITE STABILE';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = ML

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

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