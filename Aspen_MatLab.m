
clear

pause(5)

rmatdatestr = 'Base_';


Recirculation = 1.5*10^6;
Reboiler_Duty = 100000; % [kW]
Strip_pres = 1.8; %
Distilate_rate = 120000;
HX_Temp = 5;
Inter_Temp = 5; 

% Recirculation = 1.1*10^6:2.0*10^5:2.1*10^6;
% Reboiler_Duty = 80000:10000:130000; % [kW]
% Strip_pres = 1.0:0.2:2.4; %
% Distilate_rate = 100000:20000:160000;
% HX_Temp = [5 10 15];
% Inter_Temp = 5:5:15; 
mat = combvec(Recirculation,Reboiler_Duty,Strip_pres,Distilate_rate,HX_Temp,Inter_Temp);

[M, N] = size(mat);  % Matrix size
nSub = 1;           % Number of submatrices
cMat = mat2cell(mat, M, diff(round(linspace(0, N, nSub+1))));

Cell_row = cellfun('size',cMat,1); % rows

Cell_columns = max(cellfun('size',cMat,2)); % columns

h = zeros(1,Cell_columns);
for t = 1:nSub
        Captured_CO2{1,t} = h;
        Fluegas_CO2{1,t} = h;
        Removal{1,t} = h;
        Respec{1,t} = h;
        Flow_GAS4_CO2{1,t} = h;
        IMP_MFCO2{1,t} = h;
        IMP_MFH2O{1,t} = h;
        IMP_MFO2{1,t} = h;
        IMP_MFN2{1,t} = h;
        IMP_MFNO{1,t} = h;
        IMP_MFNO2{1,t} = h;
        IMP_MFN2O{1,t} = h;
        IMP_MFCO{1,t} = h;
        IMP_MFSO2{1,t} = h;
        IMP_MFNH3{1,t} = h;
        IMP_MFAR{1,t} = h;
        IMP_MFH2{1,t} = h;
        DEOX_MFCO2{1,t} = h;
        DEOX_MFH2O{1,t} = h;
        DEOX_MFO2{1,t} = h;
        DEOX_MFN2{1,t} = h;
        DEOX_MFNO{1,t} = h;
        DEOX_MFNO2{1,t} = h;
        DEOX_MFN2O{1,t} = h;
        DEOX_MFCO{1,t} = h;
        DEOX_MFSO2{1,t} = h;
        DEOX_MFNH3{1,t} = h;
        DEOX_MFAR{1,t} = h;
        DEOX_MFH2{1,t} = h;
        MFCO2{1,t} = h;
        MFH2O{1,t} = h;
        MFO2{1,t} = h;
        MFN2{1,t} = h;
        MFNO{1,t} = h;
        MFNO2{1,t} = h;
        MFN2O{1,t} = h;
        MFCO{1,t} = h;
        MFSO2{1,t} = h;
        MFNH3{1,t} = h;
        MFAR{1,t} = h;
        MFH2{1,t} = h;
        Cooler_1{1,t} = h;
        Cooler_2{1,t} = h;
        Cooler_3{1,t} = h;
        Inter_Cooler1{1,t} = h;
        Inter_Cooler2{1,t} = h;
        Inter_Cooler3{1,t} = h;
        Inter_Cooler4{1,t} = h;
        Compressor1{1,t} = h;
        Compressor2{1,t} = h;
        Compressor3{1,t} = h;
        Compressor4{1,t} = h;
        Water_removed1{1,t} = h;
        Water_removed2{1,t} = h;
        Water_removed3{1,t} = h;
        Water_removed4{1,t} = h;
        Comp1_CO2{1,t} = h;
        Comp2_CO2{1,t} = h;
        Comp3_CO2{1,t} = h;
        Comp1_H2O{1,t} = h;
        Comp2_H2O{1,t} = h;
        Comp3_H2O{1,t} = h;
        HX{1,t} = h;
        Temp_cold{1,t} = h;
        Temp_hot{1,t} = h;
        Abs_flood{1,t} = h;
        Water_flood{1,t} = h;
        Strip_flood{1,t} = h;
        Abs_dia{1,t} = h;
        Water_dia{1,t} = h;
        Strip_dia{1,t} = h;
        Pump{1,t} = h;
        Reactor{1,t} = h;
        Cooler_Reactor{1,t} = h;
        Stripper_Condenser{1,t} = h;
        Stripper_Boilup{1,t} = h;
        Stripper_Refluxrate{1,t} = h;
        Stripper_ReboilerTemp{1,t} = h;
        Rich{1,t} = h;
        Lean{1,t} = h;
        Water_out_stripper{1,t} = h;
        Lean_out{1,t} = h;
        MEA_MU{1,t} = h;
        Simulation_Convergency{1,t} = h;
        timen{1,t} = h;
        error{1,t} = string(h);
        text1{1,t} = string(h);
        text2{1,t} = string(h);
        text3{1,t} = string(h);
        text4{1,t} = string(h);
        text5{1,t} = string(h);
        text6{1,t} = string(h);
        text7{1,t} = string(h);
        text8{1,t} = string(h);
        text9{1,t} = string(h);
        fra{1,t} = h;
end


rmat = cell(1,nSub);

for i = 1:nSub
rmat{i} = zeros(88,length(h));
%Log_File{i} = zeros(1,length(h));
end

text = cell(1,nSub);

Log_File = cell(1,nSub);

for i = 1:nSub
    if exist(['Log_File_',num2str(i),'.mat']) == 2
        disp(['Log_File',num2str(i),'.mat found: Loading it'])
        LLog_File = load(['Log_File_',num2str(i),'.mat'],['Log_File']);
        Log_File{1,i} = LLog_File.Log_File;

    else
        disp(['Log_File',num2str(i),'.mat ?? Hælp: Making empty ones'])

        for i = 1:nSub
            disp(['Saving Log_File',num2str(i),'.mat'])
            Log_File = zeros(1,length(h));
            save(['Log_File_',num2str(i),'.mat'],'Log_File')
        end
        disp('Stitching the log files together into -> Log_File')
        Log_File = cell(1,nSub);
        for i = 1:nSub
        LLog_File = load(['Log_File_',num2str(i),'.mat'],'Log_File');
        Log_File{1,i} = Log_File;
        end
    end
end
disp('Done loading log files...')

for i = 1:nSub
        if  exist([rmatdatestr,num2str(i),'.mat']) == 2
            disp([rmatdatestr,num2str(i),'.mat found: Loading it...'])
            bmat{1,i} = load([rmatdatestr,num2str(i),'.mat']);
            rmat{1,i} = bmat{i}.rmat;
        else
            warning(['rmat ',num2str(i),' was not found, continuing...']);

        end
end
disp('Done loading rmat files and stiching it together...')

%%


Aspen = cell(1,nSub);

tic
%p = parpool(nSub,"SpmdEnabled",false);

disp('Parfor running going to the parfor loop...')
try
    for i = 1:nSub

        disp(['Starting Aspen ',num2str(i),'...'])
        Aspen{i} = actxserver('Apwn.Document.39.0');
        pause(0.5);
        [~ , mess]=fileattrib; % get attributes of folder (Necessary to establish the location of the simulation)
        Simulation_Name{i} = append('Aspen_Model_',num2str(i));% Aspen{i} Plus Simulation Name
        pause(0.5);
        Aspen{i}.invoke('InitFromArchive2',[mess.Name '\' Simulation_Name{i} '.apw']);
        pause(0.5);
        Aspen{i}.Visible = 1; % 1 ---> Aspen{i} is Visible; 0 ---> Aspen{i} is open but not visible
        pause(0.5);
        Aspen{i}.SuppressDialogs = 1; % Suppress windows dialogs.
        pause(15);

        for k = 1:Cell_columns

            if Log_File{1,i}(:,k) == 0

                disp(['Worker ',num2str(i),' starting simulation ',num2str(k),' out of ',num2str(Cell_columns)])

                Aspen{i}.Tree.FindNode("\Data\Streams\LEANIN\Input\TOTFLOW\MIXED").Value = cMat{1,i}(1,k); % Recirculation rate
                Aspen{i}.Tree.FindNode("\Data\Blocks\STRIPPER\Input\QN").Value = cMat{1,i}(2,k); % Reboiler
                Aspen{i}.Tree.FindNode("\Data\Blocks\STRIPPER\Input\PRES1").Value = cMat{1,i}(3,k); % Pressure in stripper
                Aspen{i}.Tree.FindNode("\Data\Blocks\STRIPPER\Input\BASIS_D").Value = cMat{1,i}(4,k); % Distilate flow
                %Aspen{i}.Tree.FindNode("\Data\Blocks\HX1\Input\TEMP").Value = cMat{1,i}(5,k); % Heatexchanger temperature
                Aspen{i}.Tree.FindNode("\Data\Blocks\INTER1\Input\TEMP").Value = cMat{1,i}(6,k); % Inter cooler temperatures

                if cMat{1,i}(5,k) == 5
                    Aspen{i}.Tree.FindNode("\Data\Flowsheeting Options\Calculator\HX\Input\FORTRAN_EXEC\#2").Value = "      IF (HOTOUT .LT. (COLDIN+5)) THEN";
                    Aspen{i}.Tree.FindNode("\Data\Flowsheeting Options\Calculator\HX\Input\FORTRAN_EXEC\#11").Value = "      IF (HOTOUT .GT. (COLDIN+7)) THEN";
                end

                if cMat{1,i}(5,k) == 10
                    Aspen{i}.Tree.FindNode("\Data\Flowsheeting Options\Calculator\HX\Input\FORTRAN_EXEC\#2").Value = "      IF (HOTOUT .LT. (COLDIN+10)) THEN";
                    Aspen{i}.Tree.FindNode("\Data\Flowsheeting Options\Calculator\HX\Input\FORTRAN_EXEC\#11").Value = "      IF (HOTOUT .GT. (COLDIN+12)) THEN";
                end

                if cMat{1,i}(5,k) == 15
                    Aspen{i}.Tree.FindNode("\Data\Flowsheeting Options\Calculator\HX\Input\FORTRAN_EXEC\#2").Value = "      IF (HOTOUT .LT. (COLDIN+15)) THEN";
                    Aspen{i}.Tree.FindNode("\Data\Flowsheeting Options\Calculator\HX\Input\FORTRAN_EXEC\#11").Value = "      IF (HOTOUT .GT. (COLDIN+17)) THEN";
                end

                pause(5);

                Aspen{i}.Reinit; % Reinit simulation
                pause(5);

                Aspen{i}.Engine.Run2(1); %Run the simulation. (1) ---> Matlab isnt busy; (0) Matlab is Busy;


                timen{1,i}(1,k) = 1;
                while Aspen{i}.Engine.IsRunning == 1 % 1 --> If Aspen{i} is running; 0 ---> If Aspen{i} stop.
                    pause(0.5)
                    timen{1,i}(1,k) = timen{1,i}(1,k) + 1;
                    if timen{1,i}(1,k) == 600 % Control of simulation time.
                        Aspen{i}.Engine.Stop;
                        pause(5)
                        Aspen{i}.Reinit;
                        pause(5)
                        Simulation_Convergency{1,i}(1,k) = 0;

                    end
                end


                if timen{1,i}(1,k) < 600
                    pause(5)
                    Simulation_Convergency{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Results Summary\Run-Status\Output\PER_ERROR").Value;

                    if Simulation_Convergency{1,i}(1,k) == 1
                        error{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Results Summary\Run-Status\Output\PER_ERROR\2").Value;
                    end
                end

                pause(3)



                % ----------------------------- RESULTS ------------------------------

                % Convergency



                if (Simulation_Convergency{1,i}(1,k) == 1) && (error{1,i}(1,k) ~= "completed with errors:")

                    % Flow rate
                    Captured_CO2{1,i}(1,k) =        Aspen{i}.Tree.FindNode("\Data\Streams\IMPCO2\Output\MASSFLOW\MIXED\CO2").Value; % [kg/h]
                    Fluegas_CO2{1,i}(1,k) =         Aspen{i}.Tree.FindNode("\Data\Streams\FLUEGAS2\Output\MASSFLOW\MIXED\CO2").Value; % [kg/h]
                    Flow_GAS4_CO2{1,i}(1,k) =       Aspen{i}.Tree.FindNode("\Data\Streams\GAS4\Output\MASSFLOW\MIXED\CO2").Value; % [kg/h]
                    Water_out_stripper{1,i}(1,k) =  Aspen{i}.Tree.FindNode("\Data\Streams\COND\Output\MASSFLMX\MIXED").Value; % [kg/h]
                    Lean_out{1,i}(1,k) =            Aspen{i}.Tree.FindNode("\Data\Streams\LEAN1\Output\MASSFLMX\MIXED").Value; % [kg/h]
                    MEA_MU{1,i}(1,k) =              Aspen{i}.Tree.FindNode("\Data\Streams\MEAMU\Output\MASSFLMX\MIXED").Value; % [kg/h]

                    % Impurities in CO2 stream after stripper
                    IMP_MFCO2{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\IMPCO2\Output\MOLEFRAC\MIXED\CO2").Value; % Molefraction of CO2 in IMPCO2
                    IMP_MFH2O{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\IMPCO2\Output\MOLEFRAC\MIXED\H2O").Value; % Molefraction of water in IMPCO2
                    IMP_MFO2{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\IMPCO2\Output\MOLEFRAC\MIXED\O2").Value; % Molefraction of in IMPCO2
                    IMP_MFN2{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\IMPCO2\Output\MOLEFRAC\MIXED\N2").Value; % Molefraction in IMPCO2
                    IMP_MFNO{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\IMPCO2\Output\MOLEFRAC\MIXED\NO").Value; % Molefraction in IMPCO2
                    IMP_MFNO2{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\IMPCO2\Output\MOLEFRAC\MIXED\NO2").Value; % Molefraction in IMPCO2
                    IMP_MFN2O{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\IMPCO2\Output\MOLEFRAC\MIXED\N2O").Value; % Molefraction in IMPCO2
                    IMP_MFCO{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\IMPCO2\Output\MOLEFRAC\MIXED\CO").Value; % Molefraction in IMPCO2
                    IMP_MFSO2{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\IMPCO2\Output\MOLEFRAC\MIXED\SO2").Value; % Molefraction in IMPCO2
                    IMP_MFNH3{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\IMPCO2\Output\MOLEFRAC\MIXED\NH3").Value; % Molefraction in IMPCO2
                    IMP_MFAR{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\IMPCO2\Output\MOLEFRAC\MIXED\AR").Value; % Molefraction in IMPCO2
                    IMP_MFH2{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\IMPCO2\Output\MOLEFRAC\MIXED\H2").Value; % Molefraction in IMPCO2

                    % Impurities in CO2 stream after reactor
                    DEOX_MFCO2{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\DEOX\Output\MOLEFRAC\MIXED\CO2").Value; % Molefraction of CO2 in DEOX
                    DEOX_MFH2O{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\DEOX\Output\MOLEFRAC\MIXED\H2O").Value; % Molefraction of water in DEOX
                    DEOX_MFO2{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\DEOX\Output\MOLEFRAC\MIXED\O2").Value; % Molefraction of in DEOX
                    DEOX_MFN2{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\DEOX\Output\MOLEFRAC\MIXED\N2").Value; % Molefraction in DEOX
                    DEOX_MFNO{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\DEOX\Output\MOLEFRAC\MIXED\NO").Value; % Molefraction in DEOX
                    DEOX_MFNO2{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\DEOX\Output\MOLEFRAC\MIXED\NO2").Value; % Molefraction in DEOX
                    DEOX_MFN2O{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\DEOX\Output\MOLEFRAC\MIXED\N2O").Value; % Molefraction in DEOX
                    DEOX_MFCO{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\DEOX\Output\MOLEFRAC\MIXED\CO").Value; % Molefraction in DEOX
                    DEOX_MFSO2{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\DEOX\Output\MOLEFRAC\MIXED\SO2").Value; % Molefraction in DEOX
                    DEOX_MFNH3{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\DEOX\Output\MOLEFRAC\MIXED\NH3").Value; % Molefraction in DEOX
                    DEOX_MFAR{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\DEOX\Output\MOLEFRAC\MIXED\AR").Value; % Molefraction in DEOX
                    DEOX_MFH2{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\DEOX\Output\MOLEFRAC\MIXED\H2").Value; % Molefraction in DEOX

                    % Impurities in final CO2 stream
                    MFCO2{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\GAS4\Output\MOLEFRAC\MIXED\CO2").Value; % Molefraction of CO2 in GAS4
                    MFH2O{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\GAS4\Output\MOLEFRAC\MIXED\H2O").Value; % Molefraction of water in GAS4
                    MFO2{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\GAS4\Output\MOLEFRAC\MIXED\O2").Value; % Molefraction of in GAS4
                    MFN2{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\GAS4\Output\MOLEFRAC\MIXED\N2").Value; % Molefraction in GAS4
                    MFNO{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\GAS4\Output\MOLEFRAC\MIXED\NO").Value; % Molefraction in GAS4
                    MFNO2{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\GAS4\Output\MOLEFRAC\MIXED\NO2").Value; % Molefraction in GAS4
                    MFN2O{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\GAS4\Output\MOLEFRAC\MIXED\N2O").Value; % Molefraction in GAS4
                    MFCO{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\GAS4\Output\MOLEFRAC\MIXED\CO").Value; % Molefraction in GAS4
                    MFSO2{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\GAS4\Output\MOLEFRAC\MIXED\SO2").Value; % Molefraction in GAS4
                    MFNH3{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\GAS4\Output\MOLEFRAC\MIXED\NH3").Value; % Molefraction in GAS4
                    MFAR{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\GAS4\Output\MOLEFRAC\MIXED\AR").Value; % Molefraction in GAS4
                    MFH2{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\GAS4\Output\MOLEFRAC\MIXED\H2").Value; % Molefraction in GAS4

                    % All coolers besides compression train
                    Cooler_1{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Blocks\COOLER\Output\QCALC").Value; % Cooler duty [kW]
                    Cooler_2{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Blocks\COOLER2\Output\QCALC").Value; % Cooler duty [kW]
                    Cooler_3{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Blocks\COOLER3\Output\QCALC").Value; % Cooler duty [kW]

                    % HX
                    HX{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\Q\Output\QCALC").Value; % HX duty [kW]
                    Temp_cold{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\RICHSOL2\Output\TEMP_OUT\MIXED").Value; % [C]
                    Temp_hot{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\LEAN2\Output\TEMP_OUT\MIXED").Value; % [C]

                    % Flodding
                    Abs_flood{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Blocks\ABSORBER\Output\CA_FLD_FAC1\INT-1\CS-2").Value; % [-]
                    Water_flood{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Blocks\ABSORBER\Output\CA_FLD_FAC1\INT-1\CS-1").Value; % [-]
                    Strip_flood{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Blocks\STRIPPER\Output\CA_FLD_FAC1\INT-1\CS-1").Value; % [-]

                    % Tower diameter
                    Abs_dia{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Blocks\ABSORBER\Output\CA_DIAM1\INT-1\CS-2").Value; % [m]
                    Water_dia{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Blocks\ABSORBER\Output\CA_DIAM1\INT-1\CS-1").Value; % [m]
                    Strip_dia{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Blocks\STRIPPER\Output\CA_DIAM1\INT-1\CS-2").Value; % [m]

                    % Intercoolers in compression train
                    Inter_Cooler1{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Blocks\INTER1\Output\QCALC").Value; % Cooler duty [kW]
                    Inter_Cooler2{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Blocks\INTER2\Output\QCALC").Value; % Cooler duty [kW]
                    Inter_Cooler3{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Blocks\INTER3\Output\QCALC").Value; % Cooler duty [kW]
                    Inter_Cooler4{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Blocks\INTER4\Output\QCALC").Value; % Cooler duty [kW]

                    % Compressors in compression train
                    Compressor1{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Blocks\B2\Output\IND_POWER").Value; % Compressor duty [kW]
                    Compressor2{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Blocks\B3\Output\IND_POWER").Value; % Compressor duty [kW]
                    Compressor3{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Blocks\B4\Output\IND_POWER").Value; % Compressor duty [kW]
                    Compressor4{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Blocks\B5\Output\IND_POWER").Value; % Compressor duty [kW]

                    % Water removed after every step in compression train
                    Water_removed1{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\W1\Output\MASSFLOW\MIXED\H2O").Value; % [kg/h]
                    Water_removed2{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\W2\Output\MASSFLOW\MIXED\H2O").Value; % [kg/h]
                    Water_removed3{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\W3\Output\MASSFLOW\MIXED\H2O").Value; % [kg/h]
                    Water_removed4{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\W4\Output\MASSFLOW\MIXED\H2O").Value; % [kg/h]

                    % Purity of CO2 after every step in compression train NOTE: the last
                    % step in compression train is MFCO2 !
                    Comp1_CO2{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\GAS1\Output\MOLEFRAC\MIXED\CO2").Value; % Molefraction of CO2 in GAS1
                    Comp2_CO2{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\GAS2\Output\MOLEFRAC\MIXED\CO2").Value; % Molefraction of CO2 in GAS2
                    Comp3_CO2{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\GAS3\Output\MOLEFRAC\MIXED\CO2").Value; % Molefraction of CO2 in GAS3

                    Comp1_H2O{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\GAS1\Output\MOLEFRAC\MIXED\H2O").Value; % Molefraction of H2O in GAS1
                    Comp2_H2O{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\GAS2\Output\MOLEFRAC\MIXED\H2O").Value; % Molefraction of H2O in GAS2
                    Comp3_H2O{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\GAS3\Output\MOLEFRAC\MIXED\H2O").Value; % Molefraction of H2O in GAS3


                    % Pump before stripper
                    Pump{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Blocks\PUMP\Output\ELEC_POWER").Value; % Pump duty [kW]

                    % Oxygen reactor after stripper
                    Reactor{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Blocks\ROXYGEN\Output\QCALC").Value; % Reactor heating duty [kW]
                    Cooler_Reactor{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Blocks\COOLERY\Output\QCALC").Value; % Cooler duty after reactor [kW]

                    % Stripper
                    Stripper_Condenser{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Blocks\STRIPPER\Output\COND_DUTY").Value; % Condenser heat duty of stripper [kW]
                    Stripper_Boilup{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Blocks\STRIPPER\Output\MASS_BR").Value; % Boil up ratio [-] mass basis
                    Stripper_Refluxrate{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Blocks\STRIPPER\Output\MASS_RR").Value; % Reflux ratio [-] mass basis
                    Stripper_ReboilerTemp{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Blocks\STRIPPER\Output\REB_TOUT").Value; % Stripper reboiler temperature out [°C]

                    Lean{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\LEAN1\Output\STRM_UPP\ML-LOAD\MIXED\LIQUID").Value;
                    Rich{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Streams\RICHSOL1\Output\STRM_UPP\ML-LOAD\MIXED\LIQUID").Value;


                    text1{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Results Summary\Run-Status\Output\PER_ERROR\1").Value;
                    text2{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Results Summary\Run-Status\Output\PER_ERROR\2").Value;
                    text3{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Results Summary\Run-Status\Output\PER_ERROR\3").Value;
                    text4{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Results Summary\Run-Status\Output\PER_ERROR\4").Value;
                    text5{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Results Summary\Run-Status\Output\PER_ERROR\5").Value;
                    text6{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Results Summary\Run-Status\Output\PER_ERROR\6").Value;
                    text7{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Results Summary\Run-Status\Output\PER_ERROR\7").Value;
                    text8{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Results Summary\Run-Status\Output\PER_ERROR\8").Value;
                    text9{1,i}(1,k) = Aspen{i}.Tree.FindNode("\Data\Results Summary\Run-Status\Output\PER_ERROR\9").Value;
                else
                    Captured_CO2{1,i}(1,k) = 0;
                    Fluegas_CO2{1,i}(1,k) = 0;
                    Flow_GAS4_CO2{1,i}(1,k) = 0;
                    IMP_MFCO2{1,i}(1,k) = 0;
                    IMP_MFH2O{1,i}(1,k) = 0;
                    IMP_MFO2{1,i}(1,k) = 0;
                    IMP_MFN2{1,i}(1,k) = 0;
                    IMP_MFNO{1,i}(1,k) = 0;
                    IMP_MFNO2{1,i}(1,k) = 0;
                    IMP_MFN2O{1,i}(1,k) = 0;
                    IMP_MFCO{1,i}(1,k) = 0;
                    IMP_MFSO2{1,i}(1,k) = 0;
                    IMP_MFNH3{1,i}(1,k) = 0;
                    IMP_MFAR{1,i}(1,k) = 0;
                    IMP_MFH2{1,i}(1,k) = 0;
                    DEOX_MFCO2{1,i}(1,k) = 0;
                    DEOX_MFH2O{1,i}(1,k) = 0;
                    DEOX_MFO2{1,i}(1,k) = 0;
                    DEOX_MFN2{1,i}(1,k) = 0;
                    DEOX_MFNO{1,i}(1,k) = 0;
                    DEOX_MFNO2{1,i}(1,k) = 0;
                    DEOX_MFN2O{1,i}(1,k) = 0;
                    DEOX_MFCO{1,i}(1,k) = 0;
                    DEOX_MFSO2{1,i}(1,k) = 0;
                    DEOX_MFNH3{1,i}(1,k) = 0;
                    DEOX_MFAR{1,i}(1,k) = 0;
                    DEOX_MFH2{1,i}(1,k) = 0;
                    MFCO2{1,i}(1,k) = 0;
                    MFH2O{1,i}(1,k) = 0;
                    MFO2{1,i}(1,k) = 0;
                    MFN2{1,i}(1,k) = 0;
                    MFNO{1,i}(1,k) = 0;
                    MFNO2{1,i}(1,k) = 0;
                    MFN2O{1,i}(1,k) = 0;
                    MFCO{1,i}(1,k) = 0;
                    MFSO2{1,i}(1,k) = 0;
                    MFNH3{1,i}(1,k) = 0;
                    MFAR{1,i}(1,k) = 0;
                    MFH2{1,i}(1,k) = 0;
                    Cooler_1{1,i}(1,k) = 0;
                    Cooler_2{1,i}(1,k) = 0;
                    Cooler_3{1,i}(1,k) = 0;
                    Inter_Cooler1{1,i}(1,k) = 0;
                    Inter_Cooler2{1,i}(1,k) = 0;
                    Inter_Cooler3{1,i}(1,k) = 0;
                    Inter_Cooler4{1,i}(1,k) = 0;
                    Compressor1{1,i}(1,k) = 0;
                    Compressor2{1,i}(1,k) = 0;
                    Compressor3{1,i}(1,k) = 0;
                    Compressor4{1,i}(1,k) = 0;
                    Water_removed1{1,i}(1,k) = 0;
                    Water_removed2{1,i}(1,k) = 0;
                    Water_removed3{1,i}(1,k) = 0;
                    Water_removed4{1,i}(1,k) = 0;
                    Comp1_CO2{1,i}(1,k) = 0;
                    Comp2_CO2{1,i}(1,k) = 0;
                    Comp3_CO2{1,i}(1,k) = 0;
                    Comp1_H2O{1,i}(1,k) = 0;
                    Comp2_H2O{1,i}(1,k) = 0;
                    Comp3_H2O{1,i}(1,k) = 0;
                    Pump{1,i}(1,k) = 0;
                    Reactor{1,i}(1,k) = 0;
                    HX{1,i}(1,k) = 0;
                    Temp_hot{1,i}(1,k) = 0;
                    Temp_cold{1,i}(1,k) = 0;
                    Abs_flood{1,i}(1,k) = 0;
                    Water_flood{1,i}(1,k) = 0;
                    Strip_flood{1,i}(1,k) = 0;
                    Abs_dia{1,i}(1,k) = 0;
                    Water_dia{1,i}(1,k) = 0;
                    Strip_dia{1,i}(1,k) = 0;
                    Cooler_Reactor{1,i}(1,k) = 0;
                    Stripper_Condenser{1,i}(1,k) = 0;
                    Stripper_Boilup{1,i}(1,k) = 0;
                    Stripper_Refluxrate{1,i}(1,k) = 0;
                    Stripper_ReboilerTemp{1,i}(1,k) = 0;
                    Lean{1,i}(1,k) = 0;
                    Rich{1,i}(1,k) = 0;
                    Water_out_stripper{1,i}(1,k) = 0;
                    Lean_out{1,i}(1,k) = 0;
                    MEA_MU{1,i}(1,k) = 0;
                    Simulation_Convergency{1,i}(1,k) = 0;

                    text1{1,i}(1,k) = 0;
                    text2{1,i}(1,k) = 0;
                    text3{1,i}(1,k) = 0;
                    text4{1,i}(1,k) = 0;
                    text5{1,i}(1,k) = 0;
                    text6{1,i}(1,k) = 0;
                    text7{1,i}(1,k) = 0;
                    text8{1,i}(1,k) = 0;
                    text9{1,i}(1,k) = 0;

                end


                rmat{1,i}(:,k) = [...
                    cMat{1,i}(1,k), ...
                    cMat{1,i}(2,k), ...
                    cMat{1,i}(3,k), ...
                    cMat{1,i}(4,k), ...
                    cMat{1,i}(5,k), ...
                    cMat{1,i}(6,k), ...
                    Captured_CO2{1,i}(1,k),...
                    Fluegas_CO2{1,i}(1,k), ...
                    Flow_GAS4_CO2{1,i}(1,k),...
                    IMP_MFCO2{1,i}(1,k),...
                    IMP_MFH2O{1,i}(1,k),...
                    IMP_MFO2{1,i}(1,k),...
                    IMP_MFN2{1,i}(1,k),...
                    IMP_MFNO{1,i}(1,k),...
                    IMP_MFNO2{1,i}(1,k),...
                    IMP_MFN2O{1,i}(1,k),...
                    IMP_MFCO{1,i}(1,k),...
                    IMP_MFSO2{1,i}(1,k),...
                    IMP_MFNH3{1,i}(1,k),...
                    IMP_MFAR{1,i}(1,k),...
                    IMP_MFH2{1,i}(1,k),...
                    DEOX_MFCO2{1,i}(1,k),...
                    DEOX_MFH2O{1,i}(1,k),...
                    DEOX_MFO2{1,i}(1,k),...
                    DEOX_MFN2{1,i}(1,k),...
                    DEOX_MFNO{1,i}(1,k),...
                    DEOX_MFNO2{1,i}(1,k),...
                    DEOX_MFN2O{1,i}(1,k),...
                    DEOX_MFCO{1,i}(1,k),...
                    DEOX_MFSO2{1,i}(1,k),...
                    DEOX_MFNH3{1,i}(1,k),...
                    DEOX_MFAR{1,i}(1,k),...
                    DEOX_MFH2{1,i}(1,k),...
                    MFCO2{1,i}(1,k),...
                    MFH2O{1,i}(1,k),...
                    MFO2{1,i}(1,k),...
                    MFN2{1,i}(1,k),...
                    MFNO{1,i}(1,k),...
                    MFNO2{1,i}(1,k), ...
                    MFN2O{1,i}(1,k),...
                    MFCO{1,i}(1,k),...
                    MFSO2{1,i}(1,k),...
                    MFNH3{1,i}(1,k),...
                    MFAR{1,i}(1,k),...
                    MFH2{1,i}(1,k),...
                    HX{1,i}(1,k),...
                    Temp_cold{1,i}(1,k),...
                    Temp_hot{1,i}(1,k),...
                    Cooler_1{1,i}(1,k),...
                    Cooler_2{1,i}(1,k),...
                    Cooler_3{1,i}(1,k),...
                    Inter_Cooler1{1,i}(1,k),...
                    Inter_Cooler2{1,i}(1,k),...
                    Inter_Cooler3{1,i}(1,k),...
                    Inter_Cooler4{1,i}(1,k),...
                    Compressor1{1,i}(1,k),...
                    Compressor2{1,i}(1,k),...
                    Compressor3{1,i}(1,k),...
                    Compressor4{1,i}(1,k),...
                    Water_removed1{1,i}(1,k),...
                    Water_removed2{1,i}(1,k),...
                    Water_removed3{1,i}(1,k),...
                    Water_removed4{1,i}(1,k),...
                    Comp1_CO2{1,i}(1,k),...
                    Comp2_CO2{1,i}(1,k),...
                    Comp3_CO2{1,i}(1,k),...
                    Comp1_H2O{1,i}(1,k),...
                    Comp2_H2O{1,i}(1,k),...
                    Comp3_H2O{1,i}(1,k),...
                    Pump{1,i}(1,k),...
                    Reactor{1,i}(1,k),...
                    Cooler_Reactor{1,i}(1,k),...
                    Stripper_Condenser{1,i}(1,k),...
                    Stripper_Boilup{1,i}(1,k),  ...
                    Stripper_Refluxrate{1,i}(1,k), ...
                    Stripper_ReboilerTemp{1,i}(1,k), ...
                    Lean{1,i}(1,k),...
                    Rich{1,i}(1,k),...
                    Water_out_stripper{1,i}(1,k),...
                    Lean_out{1,i}(1,k),...
                    MEA_MU{1,i}(1,k),...
                    Abs_flood{1,i}(1,k),...
                    Water_flood{1,i}(1,k),...
                    Strip_flood{1,i}(1,k),...
                    Abs_dia{1,i}(1,k),...
                    Water_dia{1,i}(1,k),...
                    Strip_dia{1,i}(1,k),...
                    Simulation_Convergency{1,i}(1,k)]';

                Log_File{1,i}(:,k) = 1;


                logsave(['Log_File_',num2str(i),'.mat'],Log_File{1,i});
                parsave([rmatdatestr,num2str(i),'.mat'],rmat{1,i});
                
                pause(5)

            else %Log file says it is already done, skip

            end
        end

    end


catch
    save('Log_File.mat','Log_File');
    save([rmatdatestr,'.mat'],'rmat');
    disp('Error occured in parfor, going to "catch" closing all aspen...');

    pause(5)
    system('taskkill /F /IM ASPENPLUS.EXE')
    system('taskkill /F /IM apmain.EXE')
    pause(5)

    disp('Ending parfor before restarting...');
    delete(p);
    delete(gcp('nocreate'))
    disp('Restarting the script in 15 seconds...');
    pause(15)
    Aspen_MatLab
end


pause(5)
system('taskkill /F /IM ASPENPLUS.EXE')
system('taskkill /F /IM apmain.EXE')
pause(5)

disp('Ending parfor...');
delete(p);
delete(gcp('nocreate'))

toc
save('Log_File.mat','Log_File');
save([rmatdatestr,'.mat'],'rmat')
time = toc


function parsave(fname, rmat)
  save(fname, 'rmat')
end

function logsave(fname,Log_File)
  save(fname,'Log_File')
end

