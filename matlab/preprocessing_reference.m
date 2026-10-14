% PREPROCESAMIENTO EEG: revisar las etapas manuales antes de ejecutar.
% Estructura esperada: <proyecto>/Code/matlab y <proyecto>/EEG.

%% CONFIGURACION DE RUTAS
clearvars;
clc;
scriptPath = fileparts(mfilename('fullpath'));
projectPath = fileparts(fileparts(scriptPath));
eegPath = fullfile(projectPath, 'EEG');
controlPath = fullfile(eegPath, 'CONTROLES EEGLAB S4 corrected hpf 0p05Hz');
patientPath = fullfile(eegPath, 'PACIENTES EEGLAB S4 corrected hpf 0p05Hz');

%% PROTECCION PARA LAS ETAPAS MANUALES
% Para matrices ya promediadas, ejecutar CONFIGURACION DE RUTAS y despues
% el apartado 8, sin repetir el preprocesamiento de los registros.
runPreprocessing = false;
if ~runPreprocessing
    error('ReferenceOnly:ReadDocumentation', ...
        'Revise las etapas manuales antes de habilitar runPreprocessing.');
end

%% CONTROLES

% AVISO: antes hay que tener los archivos en .set 

% 1. INSPECCIÓN VISUAL: lo primero que hay que hacer es inspeccionar cada 
% sujeto individualmente de manera visual para ver si hay algún canal 
% ruidoso e interpolarlo en ese caso

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2. CARGAR REGISTROS

% Se conservan las rutas comunes definidas al inicio.
eeglab;  % Inicia EEGLAB

cd(controlPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3. REFERENCIAR PROMEDIO Y FILTRAR

data_path = controlPath;

channel_file = resolve_channel_file();

S = dir(fullfile(data_path, '*adjusted.set'));

% Iniciar pool paralelo con 8 trabajadores

pool = gcp('nocreate');

if isempty(pool)
    parpool('local', 8);
else
    fprintf('Pool activo con %d trabajadores.\n', ...
        pool.NumWorkers);
end

% Procesar sujetos en paralelo

parfor x = 1:length(S)

    name = S(x).name;

    fprintf('Procesando %d/%d: %s\n', ...
        x, length(S), name);

    try
        % Cargar dataset
        EEG = pop_loadset( ...
            'filename', name, ...
            'filepath', data_path);

        EEG = eeg_checkset(EEG);

        % Añadir localizaciones
        EEG = pop_chanedit( ...
            EEG, ...
            'lookup', channel_file);

        EEG = eeg_checkset(EEG);

        % Referencia promedio
        EEG = pop_reref( ...
            EEG, [], ...
            'keepref', 'on');

        EEG = eeg_checkset(EEG);

        % High-pass 0.05 Hz, orden 66000 y filtrado FFT
        EEG = pop_eegfiltnew( ...
            EEG, 0.05, [], 66000, ...
            0, 0, 0, 0, 1);

        % Low-pass 30 Hz, orden 440
        EEG = pop_eegfiltnew( ...
            EEG, [], 30, 440, ...
            0, 0, 0, 0, 0);

        EEG = eeg_checkset(EEG);

        % Nombre de salida
        [~, base_name] = fileparts(name);

        output_name = ...
            [base_name '_filtros.set'];

        EEG.setname = ...
            [base_name '_filtros'];

        % Guardar en la misma carpeta
        EEG = pop_saveset( ...
            EEG, ...
            'filename', output_name, ...
            'filepath', data_path);

        fprintf('Completado: %s\n', name);

    catch ME
        fprintf( ...
            'ERROR en %s: %s\n', ...
            name, ME.message);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 4. ICA

data_path = controlPath;

S = dir(fullfile(data_path, '*filtros.set'));

% Iniciar pool paralelo si no está activo

pool = gcp('nocreate');

if isempty(pool)
    parpool('local', 8); % 8 trabajadores
else
    fprintf('Pool activo con %d trabajadores.\n', ...
        pool.NumWorkers);
end

% ICA en paralelo

parfor x = 1:length(S)

    name = S(x).name;

    fprintf('ICA %d/%d: %s\n', ...
        x, length(S), name);

    try
        % Cargar datos 0.05–30 Hz

        EEG = pop_loadset( ...
            'filename', name, ...
            'filepath', data_path);

        EEG = eeg_checkset(EEG);

        % Crear copia temporal para entrenar ICA

        EEGica = EEG;

        % Reducir de 1000 a 250 Hz
        EEGica = pop_resample(EEGica, 250);

        EEGica = eeg_checkset(EEGica);

        % High-pass 1 Hz para mejorar ICA
        % A 250 Hz, el orden recomendado es aproximadamente 826
        EEGica = pop_eegfiltnew( ...
            EEGica, 1, [], 826, ...
            0, 0, 0, 0, 0);

        EEGica = eeg_checkset(EEGica);

        % Semilla reproducible por participante

        rng(20260720 + x, 'twister');

        % Ejecutar Extended Infomax

        EEGica = pop_runica( ...
            EEGica, ...
            'icatype', 'runica', ...
            'options', {'extended', 1});

        EEGica = eeg_checkset(EEGica);

        % Transferir ICA a los datos 0.05–30 Hz y 1000 Hz

        EEG.icaweights  = EEGica.icaweights;
        EEG.icasphere   = EEGica.icasphere;
        EEG.icawinv     = EEGica.icawinv;
        EEG.icachansind = EEGica.icachansind;
        EEG.icaact      = [];

        % Documentar cómo se entrenó ICA
        EEG.etc.ica_training.algorithm = 'runica';
        EEG.etc.ica_training.extended = 1;
        EEG.etc.ica_training.highpass = 1;
        EEG.etc.ica_training.lowpass = 30;
        EEG.etc.ica_training.srate = 250;
        EEG.etc.ica_training.seed = 20260720 + x;

        EEG = eeg_checkset(EEG);

        % Guardar datos originales con ICA transferida

        [~, base_name] = fileparts(name);

        output_name = [base_name '_ICA.set'];

        EEG.setname = [base_name '_ICA'];

        EEG = pop_saveset( ...
            EEG, ...
            'filename', output_name, ...
            'filepath', data_path);

        fprintf('ICA completada: %s\n', name);

    catch ME
        fprintf( ...
            'ERROR ICA en %s: %s\n', ...
            name, ME.message);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 5. RECHAZAR COMPONENTES ARTEFACTUALES ICA AUTOMATICAMENTE

dataPath = controlPath;

S = dir(fullfile(dataPath, '*_ICA.set'));

logData = { ...
    'Archivo', ...
    'N_comps_eliminados', ...
    'Comps_eliminados'};

for x = 1:length(S)

    name = S(x).name;
    [~, stem] = fileparts(name);

    EEG = pop_loadset( ...
        'filename', name, ...
        'filepath', dataPath);

    EEG = eeg_checkset(EEG);

    EEG = pop_iclabel(EEG, 'default');

    EEG = pop_icflag(EEG, ...
        [NaN NaN; ...
         0.8 1; ...
         0.8 1; ...
         NaN NaN; ...
         NaN NaN; ...
         NaN NaN; ...
         NaN NaN]);

    rejected = find(EEG.reject.gcompreject);

    % Guardar probabilidades antes de eliminar los componentes
    probabilities = ...
        EEG.etc.ic_classification.ICLabel.classifications;

    EEG.etc.ica_qc.iclabel_probabilities = probabilities;
    EEG.etc.ica_qc.rejected_components   = rejected;
    EEG.etc.ica_qc.threshold_eye         = 0.80;
    EEG.etc.ica_qc.threshold_muscle      = 0.80;
    EEG.etc.ica_qc.iclabel_path          = which('pop_iclabel');

    logData(end+1, :) = { ...
        name, ...
        numel(rejected), ...
        mat2str(rejected)};

    if ~isempty(rejected)
        EEG = pop_subcomp(EEG, rejected, 0);
    end

    EEG = eeg_checkset(EEG);

    EEG.setname = [stem '_removed'];

    EEG = pop_saveset(EEG, ...
        'filename', [stem '_removed.set'], ...
        'filepath', dataPath);
end

% Guardar log en archivo .csv
writecell(logData, ...
    fullfile(dataPath, ...
    'ICAs_rejected_EEGLAB_Exp_PEP_2_CON_Automatico_0p05Hz.csv'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%. 6. CREAR EPOCAS

% Época maestra: desde 100 ms antes de S1 hasta 1000 ms después de S4

dataPath = controlPath;

eventFolder = fullfile(projectPath, 'Eventlist_Exp_PEP_2');

% BinLister nuevo: solamente códigos 10, 11, 20, 21, 30 y 40
bdfFile = fullfile(scriptPath, 'binlist_S4.txt');

% Carpeta para guardar un EventList independiente por sujeto
eventLogPath = fullfile( ...
    eventFolder, ...
    'EventLists_hpf_0p05Hz_CONTROLES');

if ~exist(eventLogPath, 'dir')
    mkdir(eventLogPath);
end

% Comprobar que existe el BinLister
if ~isfile(bdfFile)
    error('No se encontró el archivo BinLister: %s', bdfFile);
end

% Buscar los archivos generados después del rechazo ICA
S = dir(fullfile(dataPath, '*_removed.set'));

% Log de procesamiento
epochLog = { ...
    'Archivo', ...
    'Estado', ...
    'N_epocas', ...
    'Puntos_por_epoca', ...
    'Xmin_ms', ...
    'Xmax_ms', ...
    'Mensaje'};

for x = 1:length(S)

    name = S(x).name;
    [~, stem] = fileparts(name);

    fprintf('\nProcesando %d/%d: %s\n', ...
        x, length(S), name);

    try
        
        % Cargar dataset

        EEG = pop_loadset( ...
            'filename', name, ...
            'filepath', dataPath);

        EEG = eeg_checkset(EEG);

        % Eliminar metadatos antiguos reservados por ERPLAB

        erplabReservedFields = { ...
            'bepoch', ...
            'bini', ...
            'binlabel', ...
            'code', ...
            'codelabel', ...
            'enable', ...
            'flag', ...
            'item'};

        fieldsToRemove = intersect( ...
            erplabReservedFields, ...
            fieldnames(EEG.event), ...
            'stable');

        if ~isempty(fieldsToRemove)
            EEG.event = rmfield(EEG.event, fieldsToRemove);
        end

        % El EventList se generará de nuevo en esta ejecución

        if isfield(EEG, 'EVENTLIST')
            EEG.EVENTLIST = [];
        end

        EEG = eeg_checkset(EEG, 'eventconsistency');

        % Definir archivos EventList específicos del sujeto

        eventListInitial = fullfile( ...
            eventLogPath, ...
            [stem '_EventList_inicial.txt']);

        eventListBinned = fullfile( ...
            eventLogPath, ...
            [stem '_EventList_bins.txt']);

        % Crear EventList

        EEG = pop_creabasiceventlist( ...
            EEG, ...
            'AlphanumericCleaning', 'on', ...
            'BoundaryNumeric', {-99}, ...
            'BoundaryString', {'boundary'}, ...
            'Eventlist', eventListInitial, ...
            'Warning', 'on');

        % Asignar bins

        EEG = pop_binlister( ...
            EEG, ...
            'BDF', bdfFile, ...
            'ExportEL', eventListBinned, ...
            'IndexEL', 1, ...
            'SendEL2', 'EEG&Text', ...
            'UpdateEEG', 'on', ...
            'Voutput', 'EEG', ...
            'Warning', 'on');

        EEG = eeg_checkset(EEG);

        % Crear épocas

        EEG = pop_epochbin( ...
            EEG, ...
            [-1000 1000], ...
            'none', ...
            'Warning', 'on');

        EEG = eeg_checkset(EEG);

        % Registrar control de calidad

        epochLog(end+1, :) = { ...
            name, ...
            'OK', ...
            EEG.trials, ...
            EEG.pnts, ...
            EEG.xmin * 1000, ...
            EEG.xmax * 1000, ...
            ''};

        % Guardar dataset epocado

        EEG.setname = [stem '_EPOCH'];

        EEG = pop_saveset( ...
            EEG, ...
            'filename', [stem '_EPOCH.set'], ...
            'filepath', dataPath);

    catch ME

        % Registrar el error y continuar con el siguiente sujeto
        epochLog(end+1, :) = { ...
            name, ...
            'ERROR', ...
            NaN, ...
            NaN, ...
            NaN, ...
            NaN, ...
            ME.message};

        warning( ...
            'No se pudo procesar %s: %s', ...
            name, ME.message);
    end
end

% Guardar log

logFile = fullfile( ...
    dataPath, ...
    'Epoching_QC_hpf_0p05Hz_CONTROLES.csv');

writecell(epochLog, logFile);

fprintf('\nCreación de épocas finalizada.\n');
fprintf('Log guardado en:\n%s\n', logFile);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 7. RECHAZAR ARTEFACTOS Y PROMEDIAR POR SEPARADO CNV/PINV Y MMN

% Rutas

inputPath = controlPath;

outputPath = eegPath;

csvOutput = fullfile(projectPath, 'EEGLAB_accep_rej_trials_CON_0p05Hz_componentes.csv');

% Parámetros de rechazo por amplitud absoluta

channelsToCheck = 1:16;

absoluteThreshold = [-100 100];

% Ventana CNV/PINV respecto a S4

cnvPinvStartMs = -1000;
cnvPinvEndMs = 1000;

cnvPinvBaseline = [-1000 -900];

% Ventana MMN respecto a S4

mmnStartMs = -50;
mmnEndMs = 300;

mmnBaseline = [-50 0];

% Localizar datasets epocados

S = dir(fullfile(inputPath, '*_EPOCH.set'));

% Ordenar los archivos alfabéticamente

[~, sortIndex] = sort({S.name});
S = S(sortIndex);

nSubjects = length(S);

if nSubjects == 0
    error('No se encontraron archivos *_EPOCH.set en la carpeta de entrada.');
end

% Inicializar nombres y estado de procesamiento

subject_names = strings(nSubjects, 1);
processing_status = strings(nSubjects, 1);

initial_trials = nan(nSubjects, 1);

% Inicializar resultados CNV/PINV

cnvPinv_std_accepted = nan(nSubjects, 1);
cnvPinv_std_rejected = nan(nSubjects, 1);

cnvPinv_dev_accepted = nan(nSubjects, 1);
cnvPinv_dev_rejected = nan(nSubjects, 1);

cnvPinv_bin2_accepted = nan(nSubjects, 1);
cnvPinv_bin4_accepted = nan(nSubjects, 1);

cnvPinv_total_rejection_pct = nan(nSubjects, 1);
cnvPinv_dev_rejection_pct = nan(nSubjects, 1);

% Inicializar resultados MMN

mmn_std_accepted = nan(nSubjects, 1);
mmn_std_rejected = nan(nSubjects, 1);

mmn_dev_accepted = nan(nSubjects, 1);
mmn_dev_rejected = nan(nSubjects, 1);

mmn_bin2_accepted = nan(nSubjects, 1);
mmn_bin4_accepted = nan(nSubjects, 1);

mmn_total_rejection_pct = nan(nSubjects, 1);
mmn_dev_rejection_pct = nan(nSubjects, 1);

% Inicializar matrices grupales

ERPtotal_CON_CNV_PINV = [];
ERRORtotal_CON_CNV_PINV = [];

ERPtotal_CON_MMN = [];
ERRORtotal_CON_MMN = [];

% Procesar cada sujeto

for x = 1:nSubjects

    name = S(x).name;

    subject_names(x) = string(name);

    fprintf('\nProcesando %d/%d: %s\n', x, nSubjects, name);

    try

        % Cargar el dataset epocado original

        EEG = pop_loadset( ...
            'filename', name, ...
            'filepath', inputPath);

        EEG = eeg_checkset(EEG);

        initial_trials(x) = EEG.trials;

        % Crear dos copias independientes del dataset

        EEG_CNV_PINV = EEG;
        EEG_MMN = EEG;

        % Procesamiento específico para CNV/PINV

        EEG_CNV_PINV = pop_select( ...
            EEG_CNV_PINV, ...
            'time', [cnvPinvStartMs cnvPinvEndMs] / 1000);

        EEG_CNV_PINV = eeg_checkset(EEG_CNV_PINV);

        % Corrección de línea base pre-S1

        EEG_CNV_PINV = pop_rmbase( ...
            EEG_CNV_PINV, ...
            cnvPinvBaseline);

        EEG_CNV_PINV = eeg_checkset(EEG_CNV_PINV);

        % Umbral absoluto para CNV/PINV

        EEG_CNV_PINV = pop_artextval( ...
            EEG_CNV_PINV, ...
            'Channel', channelsToCheck, ...
            'Flag', 1, ...
            'Threshold', absoluteThreshold, ...
            'Twindow', [cnvPinvStartMs cnvPinvEndMs]);

        EEG_CNV_PINV = eeg_checkset(EEG_CNV_PINV);

        % Obtener los ensayos aceptados y rechazados para CNV/PINV

        [EEG_CNV_PINV, tprej_CNV_PINV, acce_CNV_PINV, ...
            rej_CNV_PINV, histoflags_CNV_PINV] = ...
            pop_summary_AR_eeg_detection(EEG_CNV_PINV, 'none');

        % Combinar los indicadores de rechazo

        EEG_CNV_PINV = eeg_rejsuperpose( ...
            EEG_CNV_PINV, ...
            1, 1, 0, 0, 0, 0, 0, 0);

        cnvPinvRejectMask = logical(EEG_CNV_PINV.reject.rejglobal);

        cnvPinv_total_rejection_pct(x) = ...
            100 * sum(cnvPinvRejectMask) / EEG_CNV_PINV.trials;

        % Guardar recuentos por condición

        cnvPinv_std_accepted(x) = ...
            sum(acce_CNV_PINV([1 3 5 6]));

        cnvPinv_std_rejected(x) = ...
            sum(rej_CNV_PINV([1 3 5 6]));

        cnvPinv_dev_accepted(x) = ...
            sum(acce_CNV_PINV([2 4]));

        cnvPinv_dev_rejected(x) = ...
            sum(rej_CNV_PINV([2 4]));

        cnvPinv_bin2_accepted(x) = acce_CNV_PINV(2);
        cnvPinv_bin4_accepted(x) = acce_CNV_PINV(4);

        cnvPinvTotalDev = ...
            cnvPinv_dev_accepted(x) + cnvPinv_dev_rejected(x);

        if cnvPinvTotalDev > 0
            cnvPinv_dev_rejection_pct(x) = ...
                100 * cnvPinv_dev_rejected(x) / cnvPinvTotalDev;
        end

        % Eliminar los ensayos rechazados para CNV/PINV

        EEG_CNV_PINV = pop_rejepoch( ...
            EEG_CNV_PINV, ...
            cnvPinvRejectMask, ...
            0);

        EEG_CNV_PINV = eeg_checkset(EEG_CNV_PINV);

        % Promediar los ensayos aceptados para CNV/PINV

        ERP_CNV_PINV = pop_averager( ...
            EEG_CNV_PINV, ...
            'Criterion', 'good', ...
            'ExcludeBoundary', 'on', ...
            'SEM', 'on');

        % Inicializar matrices grupales después del primer sujeto válido

        if isempty(ERPtotal_CON_CNV_PINV)

            cnvPinvSize = size(ERP_CNV_PINV.bindata);

            if length(cnvPinvSize) < 3
                cnvPinvSize(3) = 1;
            end

            ERPtotal_CON_CNV_PINV = nan( ...
                cnvPinvSize(1), ...
                cnvPinvSize(2), ...
                cnvPinvSize(3), ...
                nSubjects);

            ERRORtotal_CON_CNV_PINV = nan( ...
                cnvPinvSize(1), ...
                cnvPinvSize(2), ...
                cnvPinvSize(3), ...
                nSubjects);
        end

        ERPtotal_CON_CNV_PINV(:, :, :, x) = ...
            ERP_CNV_PINV.bindata;

        ERRORtotal_CON_CNV_PINV(:, :, :, x) = ...
            ERP_CNV_PINV.binerror;

        % Procesamiento específico para MMN

        EEG_MMN = pop_select( ...
            EEG_MMN, ...
            'time', [mmnStartMs mmnEndMs] / 1000);

        EEG_MMN = eeg_checkset(EEG_MMN);

        % Corrección de línea base pre-S4

        EEG_MMN = pop_rmbase( ...
            EEG_MMN, ...
            mmnBaseline);

        EEG_MMN = eeg_checkset(EEG_MMN);

        % Umbral absoluto para MMN

        EEG_MMN = pop_artextval( ...
            EEG_MMN, ...
            'Channel', channelsToCheck, ...
            'Flag', 1, ...
            'Threshold', absoluteThreshold, ...
            'Twindow', [mmnStartMs mmnEndMs]);

        EEG_MMN = eeg_checkset(EEG_MMN);

        % Obtener los ensayos aceptados y rechazados para MMN

        [EEG_MMN, tprej_MMN, acce_MMN, ...
            rej_MMN, histoflags_MMN] = ...
            pop_summary_AR_eeg_detection(EEG_MMN, 'none');

        % Combinar los indicadores de rechazo

        EEG_MMN = eeg_rejsuperpose( ...
            EEG_MMN, ...
            1, 1, 0, 0, 0, 0, 0, 0);

        mmnRejectMask = logical(EEG_MMN.reject.rejglobal);

        mmn_total_rejection_pct(x) = ...
            100 * sum(mmnRejectMask) / EEG_MMN.trials;

        % Guardar recuentos por condición

        mmn_std_accepted(x) = ...
            sum(acce_MMN([1 3 5 6]));

        mmn_std_rejected(x) = ...
            sum(rej_MMN([1 3 5 6]));

        mmn_dev_accepted(x) = ...
            sum(acce_MMN([2 4]));

        mmn_dev_rejected(x) = ...
            sum(rej_MMN([2 4]));

        mmn_bin2_accepted(x) = acce_MMN(2);
        mmn_bin4_accepted(x) = acce_MMN(4);

        mmnTotalDev = ...
            mmn_dev_accepted(x) + mmn_dev_rejected(x);

        if mmnTotalDev > 0
            mmn_dev_rejection_pct(x) = ...
                100 * mmn_dev_rejected(x) / mmnTotalDev;
        end

        % Eliminar los ensayos rechazados para MMN

        EEG_MMN = pop_rejepoch( ...
            EEG_MMN, ...
            mmnRejectMask, ...
            0);

        EEG_MMN = eeg_checkset(EEG_MMN);

        % Promediar los ensayos aceptados para MMN

        ERP_MMN = pop_averager( ...
            EEG_MMN, ...
            'Criterion', 'good', ...
            'ExcludeBoundary', 'on', ...
            'SEM', 'on');

        % Inicializar matrices grupales después del primer sujeto válido

        if isempty(ERPtotal_CON_MMN)

            mmnSize = size(ERP_MMN.bindata);

            if length(mmnSize) < 3
                mmnSize(3) = 1;
            end

            ERPtotal_CON_MMN = nan( ...
                mmnSize(1), ...
                mmnSize(2), ...
                mmnSize(3), ...
                nSubjects);

            ERRORtotal_CON_MMN = nan( ...
                mmnSize(1), ...
                mmnSize(2), ...
                mmnSize(3), ...
                nSubjects);
        end

        ERPtotal_CON_MMN(:, :, :, x) = ...
            ERP_MMN.bindata;

        ERRORtotal_CON_MMN(:, :, :, x) = ...
            ERP_MMN.binerror;

        processing_status(x) = "OK";

        fprintf(['CNV/PINV: %d desviantes aceptados; ' ...
            'MMN: %d desviantes aceptados\n'], ...
            cnvPinv_dev_accepted(x), ...
            mmn_dev_accepted(x));

    catch ME

        processing_status(x) = string(ME.message);

        warning( ...
            'No se pudo procesar %s: %s', ...
            name, ...
            ME.message);
    end
end

% Crear tabla conjunta de control de calidad

results_table = table( ...
    subject_names, ...
    processing_status, ...
    initial_trials, ...
    cnvPinv_std_accepted, ...
    cnvPinv_std_rejected, ...
    cnvPinv_dev_accepted, ...
    cnvPinv_dev_rejected, ...
    cnvPinv_bin2_accepted, ...
    cnvPinv_bin4_accepted, ...
    cnvPinv_total_rejection_pct, ...
    cnvPinv_dev_rejection_pct, ...
    mmn_std_accepted, ...
    mmn_std_rejected, ...
    mmn_dev_accepted, ...
    mmn_dev_rejected, ...
    mmn_bin2_accepted, ...
    mmn_bin4_accepted, ...
    mmn_total_rejection_pct, ...
    mmn_dev_rejection_pct, ...
    'VariableNames', { ...
    'Subject', ...
    'Status', ...
    'InitialTrials', ...
    'CNV_PINV_StdAccepted', ...
    'CNV_PINV_StdRejected', ...
    'CNV_PINV_DevAccepted', ...
    'CNV_PINV_DevRejected', ...
    'CNV_PINV_Bin2Accepted', ...
    'CNV_PINV_Bin4Accepted', ...
    'CNV_PINV_TotalRejectionPct', ...
    'CNV_PINV_DevRejectionPct', ...
    'MMN_StdAccepted', ...
    'MMN_StdRejected', ...
    'MMN_DevAccepted', ...
    'MMN_DevRejected', ...
    'MMN_Bin2Accepted', ...
    'MMN_Bin4Accepted', ...
    'MMN_TotalRejectionPct', ...
    'MMN_DevRejectionPct'});

% Guardar tabla de control de calidad

writetable(results_table, csvOutput);

% Guardar resultados CNV/PINV

save( ...
    fullfile(outputPath, ...
    'ERPtotal_CON_CNV_PINV_0p05Hz.mat'), ...
    'ERPtotal_CON_CNV_PINV', ...
    'ERRORtotal_CON_CNV_PINV', ...
    'subject_names', ...
    'results_table', ...
    'cnvPinvStartMs', ...
    'cnvPinvEndMs', ...
    'cnvPinvBaseline', ...
    '-v7.3');

% Guardar resultados MMN

save( ...
    fullfile(outputPath, ...
    'ERPtotal_CON_MMN_0p05Hz.mat'), ...
    'ERPtotal_CON_MMN', ...
    'ERRORtotal_CON_MMN', ...
    'subject_names', ...
    'results_table', ...
    'mmnStartMs', ...
    'mmnEndMs', ...
    'mmnBaseline', ...
    '-v7.3');

fprintf('\nProcesamiento completado.\n');
fprintf('Tabla de control de calidad: %s\n', csvOutput);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PACIENTES

% AVISO: antes hay que tener los archivos en .set 

% 1. INSPECCIÓN VISUAL: lo primero que hay que hacer es inspeccionar cada 
% sujeto individualmente de manera visual para ver si hay algún canal 
% ruidoso e interpolarlo en ese caso

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2. CARGAR REGISTROS

% Se conservan las rutas comunes definidas al inicio.
eeglab;  % Inicia EEGLAB

cd(patientPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3. REFERENCIAR PROMEDIO Y FILTRAR

data_path = patientPath;

channel_file = resolve_channel_file();

S = dir(fullfile(data_path, '*adjusted.set'));

% Iniciar pool paralelo con 8 trabajadores

pool = gcp('nocreate');

if isempty(pool)
    parpool('local', 8);
else
    fprintf('Pool activo con %d trabajadores.\n', ...
        pool.NumWorkers);
end

% Procesar sujetos en paralelo

parfor x = 1:length(S)

    name = S(x).name;

    fprintf('Procesando %d/%d: %s\n', ...
        x, length(S), name);

    try
        % Cargar dataset
        EEG = pop_loadset( ...
            'filename', name, ...
            'filepath', data_path);

        EEG = eeg_checkset(EEG);

        % Añadir localizaciones
        EEG = pop_chanedit( ...
            EEG, ...
            'lookup', channel_file);

        EEG = eeg_checkset(EEG);

        % Referencia promedio
        EEG = pop_reref( ...
            EEG, [], ...
            'keepref', 'on');

        EEG = eeg_checkset(EEG);

        % High-pass 0.05 Hz, orden 66000 y filtrado FFT
        EEG = pop_eegfiltnew( ...
            EEG, 0.05, [], 66000, ...
            0, 0, 0, 0, 1);

        % Low-pass 30 Hz, orden 440
        EEG = pop_eegfiltnew( ...
            EEG, [], 30, 440, ...
            0, 0, 0, 0, 0);

        EEG = eeg_checkset(EEG);

        % Nombre de salida
        [~, base_name] = fileparts(name);

        output_name = ...
            [base_name '_filtros.set'];

        EEG.setname = ...
            [base_name '_filtros'];

        % Guardar en la misma carpeta
        EEG = pop_saveset( ...
            EEG, ...
            'filename', output_name, ...
            'filepath', data_path);

        fprintf('Completado: %s\n', name);

    catch ME
        fprintf( ...
            'ERROR en %s: %s\n', ...
            name, ME.message);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 4. ICA

data_path = patientPath;

S = dir(fullfile(data_path, '*filtros.set'));

% Iniciar pool paralelo si no está activo

pool = gcp('nocreate');

if isempty(pool)
    parpool('local', 8); % 8 trabajadores
else
    fprintf('Pool activo con %d trabajadores.\n', ...
        pool.NumWorkers);
end

% ICA en paralelo

parfor x = 1:length(S)

    name = S(x).name;

    fprintf('ICA %d/%d: %s\n', ...
        x, length(S), name);

    try
        % Cargar datos 0.05–30 Hz

        EEG = pop_loadset( ...
            'filename', name, ...
            'filepath', data_path);

        EEG = eeg_checkset(EEG);

        % Crear copia temporal para entrenar ICA

        EEGica = EEG;

        % Reducir de 1000 a 250 Hz
        EEGica = pop_resample(EEGica, 250);

        EEGica = eeg_checkset(EEGica);

        % High-pass 1 Hz para mejorar ICA
        % A 250 Hz, el orden recomendado es aproximadamente 826
        EEGica = pop_eegfiltnew( ...
            EEGica, 1, [], 826, ...
            0, 0, 0, 0, 0);

        EEGica = eeg_checkset(EEGica);

        % Semilla reproducible por participante

        rng(20260720 + x, 'twister');

        % Ejecutar Extended Infomax

        EEGica = pop_runica( ...
            EEGica, ...
            'icatype', 'runica', ...
            'options', {'extended', 1});

        EEGica = eeg_checkset(EEGica);

        % Transferir ICA a los datos 0.05–30 Hz y 1000 Hz

        EEG.icaweights  = EEGica.icaweights;
        EEG.icasphere   = EEGica.icasphere;
        EEG.icawinv     = EEGica.icawinv;
        EEG.icachansind = EEGica.icachansind;
        EEG.icaact      = [];

        % Documentar cómo se entrenó ICA
        EEG.etc.ica_training.algorithm = 'runica';
        EEG.etc.ica_training.extended = 1;
        EEG.etc.ica_training.highpass = 1;
        EEG.etc.ica_training.lowpass = 30;
        EEG.etc.ica_training.srate = 250;
        EEG.etc.ica_training.seed = 20260720 + x;

        EEG = eeg_checkset(EEG);

        % Guardar datos originales con ICA transferida

        [~, base_name] = fileparts(name);

        output_name = [base_name '_ICA.set'];

        EEG.setname = [base_name '_ICA'];

        EEG = pop_saveset( ...
            EEG, ...
            'filename', output_name, ...
            'filepath', data_path);

        fprintf('ICA completada: %s\n', name);

    catch ME
        fprintf( ...
            'ERROR ICA en %s: %s\n', ...
            name, ME.message);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 5. RECHAZAR COMPONENTES ARTEFACTUALES ICA AUTOMATICAMENTE

dataPath = patientPath;

S = dir(fullfile(dataPath, '*_ICA.set'));

logData = { ...
    'Archivo', ...
    'N_comps_eliminados', ...
    'Comps_eliminados'};

for x = 1:length(S)

    name = S(x).name;
    [~, stem] = fileparts(name);

    EEG = pop_loadset( ...
        'filename', name, ...
        'filepath', dataPath);

    EEG = eeg_checkset(EEG);

    EEG = pop_iclabel(EEG, 'default');

    EEG = pop_icflag(EEG, ...
        [NaN NaN; ...
         0.8 1; ...
         0.8 1; ...
         NaN NaN; ...
         NaN NaN; ...
         NaN NaN; ...
         NaN NaN]);

    rejected = find(EEG.reject.gcompreject);

    % Guardar probabilidades antes de eliminar los componentes
    probabilities = ...
        EEG.etc.ic_classification.ICLabel.classifications;

    EEG.etc.ica_qc.iclabel_probabilities = probabilities;
    EEG.etc.ica_qc.rejected_components   = rejected;
    EEG.etc.ica_qc.threshold_eye         = 0.80;
    EEG.etc.ica_qc.threshold_muscle      = 0.80;
    EEG.etc.ica_qc.iclabel_path          = which('pop_iclabel');

    logData(end+1, :) = { ...
        name, ...
        numel(rejected), ...
        mat2str(rejected)};

    if ~isempty(rejected)
        EEG = pop_subcomp(EEG, rejected, 0);
    end

    EEG = eeg_checkset(EEG);

    EEG.setname = [stem '_removed'];

    EEG = pop_saveset(EEG, ...
        'filename', [stem '_removed.set'], ...
        'filepath', dataPath);
end

% Guardar log en archivo .csv
writecell(logData, ...
    fullfile(dataPath, ...
    'ICAs_rejected_EEGLAB_Exp_PEP_2_PEP_Automatico_0p05Hz.csv'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%. 6. CREAR EPOCAS

% Época maestra: desde 100 ms antes de S1 hasta 1000 ms después de S4

dataPath = patientPath;

eventFolder = fullfile(projectPath, 'Eventlist_Exp_PEP_2');

% BinLister nuevo: solamente códigos 10, 11, 20, 21, 30 y 40
bdfFile = fullfile(scriptPath, 'binlist_S4.txt');

% Carpeta para guardar un EventList independiente por sujeto
eventLogPath = fullfile( ...
    eventFolder, ...
    'EventLists_hpf_0p05Hz_PACIENTES');

if ~exist(eventLogPath, 'dir')
    mkdir(eventLogPath);
end

% Comprobar que existe el BinLister
if ~isfile(bdfFile)
    error('No se encontró el archivo BinLister: %s', bdfFile);
end

% Buscar los archivos generados después del rechazo ICA
S = dir(fullfile(dataPath, '*_removed.set'));

% Log de procesamiento
epochLog = { ...
    'Archivo', ...
    'Estado', ...
    'N_epocas', ...
    'Puntos_por_epoca', ...
    'Xmin_ms', ...
    'Xmax_ms', ...
    'Mensaje'};

for x = 1:length(S)

    name = S(x).name;
    [~, stem] = fileparts(name);

    fprintf('\nProcesando %d/%d: %s\n', ...
        x, length(S), name);

    try
        % Cargar dataset

        EEG = pop_loadset( ...
            'filename', name, ...
            'filepath', dataPath);

        EEG = eeg_checkset(EEG);

  % Eliminar metadatos antiguos reservados por ERPLAB

        erplabReservedFields = { ...
            'bepoch', ...
            'bini', ...
            'binlabel', ...
            'code', ...
            'codelabel', ...
            'enable', ...
            'flag', ...
            'item'};

        fieldsToRemove = intersect( ...
            erplabReservedFields, ...
            fieldnames(EEG.event), ...
            'stable');

        if ~isempty(fieldsToRemove)
            EEG.event = rmfield(EEG.event, fieldsToRemove);
        end

        % El EventList se generará de nuevo en esta ejecución

        if isfield(EEG, 'EVENTLIST')
            EEG.EVENTLIST = [];
        end

        EEG = eeg_checkset(EEG, 'eventconsistency');

        % Definir archivos EventList específicos del sujeto

        eventListInitial = fullfile( ...
            eventLogPath, ...
            [stem '_EventList_inicial.txt']);

        eventListBinned = fullfile( ...
            eventLogPath, ...
            [stem '_EventList_bins.txt']);

        % Crear EventList

        EEG = pop_creabasiceventlist( ...
            EEG, ...
            'AlphanumericCleaning', 'on', ...
            'BoundaryNumeric', {-99}, ...
            'BoundaryString', {'boundary'}, ...
            'Eventlist', eventListInitial, ...
            'Warning', 'on');

        % Asignar bins

        EEG = pop_binlister( ...
            EEG, ...
            'BDF', bdfFile, ...
            'ExportEL', eventListBinned, ...
            'IndexEL', 1, ...
            'SendEL2', 'EEG&Text', ...
            'UpdateEEG', 'on', ...
            'Voutput', 'EEG', ...
            'Warning', 'on');

        EEG = eeg_checkset(EEG);

        % Crear épocas

        EEG = pop_epochbin( ...
            EEG, ...
            [-1000 1000], ...
            'none', ...
            'Warning', 'on');

        EEG = eeg_checkset(EEG);

        % Registrar control de calidad

        epochLog(end+1, :) = { ...
            name, ...
            'OK', ...
            EEG.trials, ...
            EEG.pnts, ...
            EEG.xmin * 1000, ...
            EEG.xmax * 1000, ...
            ''};

        % Guardar dataset epocado

        EEG.setname = [stem '_EPOCH'];

        EEG = pop_saveset( ...
            EEG, ...
            'filename', [stem '_EPOCH.set'], ...
            'filepath', dataPath);

    catch ME

        % Registrar el error y continuar con el siguiente sujeto
        epochLog(end+1, :) = { ...
            name, ...
            'ERROR', ...
            NaN, ...
            NaN, ...
            NaN, ...
            NaN, ...
            ME.message};

        warning( ...
            'No se pudo procesar %s: %s', ...
            name, ME.message);
    end
end

% Guardar log

logFile = fullfile( ...
    dataPath, ...
    'Epoching_QC_hpf_0p05Hz_PACIENTES.csv');

writecell(epochLog, logFile);

fprintf('\nCreación de épocas finalizada.\n');
fprintf('Log guardado en:\n%s\n', logFile);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 7. RECHAZAR ARTEFACTOS Y PROMEDIAR POR SEPARADO CNV/PINV Y MMN

% Rutas

inputPath = patientPath;

outputPath = eegPath;

csvOutput = fullfile(projectPath, 'EEGLAB_accep_rej_trials_PEP_0p05Hz_componentes.csv');

% Parámetros de rechazo por amplitud absoluta

channelsToCheck = 1:16;

absoluteThreshold = [-100 100];

% Ventana CNV/PINV respecto a S4

cnvPinvStartMs = -1000;
cnvPinvEndMs = 1000;

cnvPinvBaseline = [-1000 -900];

% Ventana MMN respecto a S4

mmnStartMs = -50;
mmnEndMs = 300;

mmnBaseline = [-50 0];

% Localizar datasets epocados

S = dir(fullfile(inputPath, '*_EPOCH.set'));

% Ordenar los archivos alfabéticamente

[~, sortIndex] = sort({S.name});
S = S(sortIndex);

nSubjects = length(S);

if nSubjects == 0
    error('No se encontraron archivos *_EPOCH.set en la carpeta de entrada.');
end

% Inicializar nombres y estado de procesamiento

subject_names = strings(nSubjects, 1);
processing_status = strings(nSubjects, 1);

initial_trials = nan(nSubjects, 1);

% Inicializar resultados CNV/PINV

cnvPinv_std_accepted = nan(nSubjects, 1);
cnvPinv_std_rejected = nan(nSubjects, 1);

cnvPinv_dev_accepted = nan(nSubjects, 1);
cnvPinv_dev_rejected = nan(nSubjects, 1);

cnvPinv_bin2_accepted = nan(nSubjects, 1);
cnvPinv_bin4_accepted = nan(nSubjects, 1);

cnvPinv_total_rejection_pct = nan(nSubjects, 1);
cnvPinv_dev_rejection_pct = nan(nSubjects, 1);

% Inicializar resultados MMN

mmn_std_accepted = nan(nSubjects, 1);
mmn_std_rejected = nan(nSubjects, 1);

mmn_dev_accepted = nan(nSubjects, 1);
mmn_dev_rejected = nan(nSubjects, 1);

mmn_bin2_accepted = nan(nSubjects, 1);
mmn_bin4_accepted = nan(nSubjects, 1);

mmn_total_rejection_pct = nan(nSubjects, 1);
mmn_dev_rejection_pct = nan(nSubjects, 1);

% Inicializar matrices grupales

ERPtotal_PEP_CNV_PINV = [];
ERRORtotal_PEP_CNV_PINV = [];

ERPtotal_PEP_MMN = [];
ERRORtotal_PEP_MMN = [];

% Procesar cada sujeto

for x = 1:nSubjects

    name = S(x).name;

    subject_names(x) = string(name);

    fprintf('\nProcesando %d/%d: %s\n', x, nSubjects, name);

    try

        % Cargar el dataset epocado original

        EEG = pop_loadset( ...
            'filename', name, ...
            'filepath', inputPath);

        EEG = eeg_checkset(EEG);

        initial_trials(x) = EEG.trials;

        % Crear dos copias independientes del dataset

        EEG_CNV_PINV = EEG;
        EEG_MMN = EEG;

        % Procesamiento específico para CNV/PINV

        EEG_CNV_PINV = pop_select( ...
            EEG_CNV_PINV, ...
            'time', [cnvPinvStartMs cnvPinvEndMs] / 1000);

        EEG_CNV_PINV = eeg_checkset(EEG_CNV_PINV);

        % Corrección de línea base pre-S1

        EEG_CNV_PINV = pop_rmbase( ...
            EEG_CNV_PINV, ...
            cnvPinvBaseline);

        EEG_CNV_PINV = eeg_checkset(EEG_CNV_PINV);

        % Umbral absoluto para CNV/PINV

        EEG_CNV_PINV = pop_artextval( ...
            EEG_CNV_PINV, ...
            'Channel', channelsToCheck, ...
            'Flag', 1, ...
            'Threshold', absoluteThreshold, ...
            'Twindow', [cnvPinvStartMs cnvPinvEndMs]);

        EEG_CNV_PINV = eeg_checkset(EEG_CNV_PINV);

        % Obtener los ensayos aceptados y rechazados para CNV/PINV

        [EEG_CNV_PINV, tprej_CNV_PINV, acce_CNV_PINV, ...
            rej_CNV_PINV, histoflags_CNV_PINV] = ...
            pop_summary_AR_eeg_detection(EEG_CNV_PINV, 'none');

        % Combinar los indicadores de rechazo

        EEG_CNV_PINV = eeg_rejsuperpose( ...
            EEG_CNV_PINV, ...
            1, 1, 0, 0, 0, 0, 0, 0);

        cnvPinvRejectMask = logical(EEG_CNV_PINV.reject.rejglobal);

        cnvPinv_total_rejection_pct(x) = ...
            100 * sum(cnvPinvRejectMask) / EEG_CNV_PINV.trials;

        % Guardar recuentos por condición

        cnvPinv_std_accepted(x) = ...
            sum(acce_CNV_PINV([1 3 5 6]));

        cnvPinv_std_rejected(x) = ...
            sum(rej_CNV_PINV([1 3 5 6]));

        cnvPinv_dev_accepted(x) = ...
            sum(acce_CNV_PINV([2 4]));

        cnvPinv_dev_rejected(x) = ...
            sum(rej_CNV_PINV([2 4]));

        cnvPinv_bin2_accepted(x) = acce_CNV_PINV(2);
        cnvPinv_bin4_accepted(x) = acce_CNV_PINV(4);

        cnvPinvTotalDev = ...
            cnvPinv_dev_accepted(x) + cnvPinv_dev_rejected(x);

        if cnvPinvTotalDev > 0
            cnvPinv_dev_rejection_pct(x) = ...
                100 * cnvPinv_dev_rejected(x) / cnvPinvTotalDev;
        end

        % Eliminar los ensayos rechazados para CNV/PINV

        EEG_CNV_PINV = pop_rejepoch( ...
            EEG_CNV_PINV, ...
            cnvPinvRejectMask, ...
            0);

        EEG_CNV_PINV = eeg_checkset(EEG_CNV_PINV);

        % Promediar los ensayos aceptados para CNV/PINV

        ERP_CNV_PINV = pop_averager( ...
            EEG_CNV_PINV, ...
            'Criterion', 'good', ...
            'ExcludeBoundary', 'on', ...
            'SEM', 'on');

        % Inicializar matrices grupales después del primer sujeto válido

        if isempty(ERPtotal_PEP_CNV_PINV)

            cnvPinvSize = size(ERP_CNV_PINV.bindata);

            if length(cnvPinvSize) < 3
                cnvPinvSize(3) = 1;
            end

            ERPtotal_PEP_CNV_PINV = nan( ...
                cnvPinvSize(1), ...
                cnvPinvSize(2), ...
                cnvPinvSize(3), ...
                nSubjects);

            ERRORtotal_PEP_CNV_PINV = nan( ...
                cnvPinvSize(1), ...
                cnvPinvSize(2), ...
                cnvPinvSize(3), ...
                nSubjects);
        end

        ERPtotal_PEP_CNV_PINV(:, :, :, x) = ...
            ERP_CNV_PINV.bindata;

        ERRORtotal_PEP_CNV_PINV(:, :, :, x) = ...
            ERP_CNV_PINV.binerror;

        % Procesamiento específico para MMN

        EEG_MMN = pop_select( ...
            EEG_MMN, ...
            'time', [mmnStartMs mmnEndMs] / 1000);

        EEG_MMN = eeg_checkset(EEG_MMN);

        % Corrección de línea base pre-S4

        EEG_MMN = pop_rmbase( ...
            EEG_MMN, ...
            mmnBaseline);

        EEG_MMN = eeg_checkset(EEG_MMN);

        % Umbral absoluto para MMN

        EEG_MMN = pop_artextval( ...
            EEG_MMN, ...
            'Channel', channelsToCheck, ...
            'Flag', 1, ...
            'Threshold', absoluteThreshold, ...
            'Twindow', [mmnStartMs mmnEndMs]);

        EEG_MMN = eeg_checkset(EEG_MMN);

        % Obtener los ensayos aceptados y rechazados para MMN

        [EEG_MMN, tprej_MMN, acce_MMN, ...
            rej_MMN, histoflags_MMN] = ...
            pop_summary_AR_eeg_detection(EEG_MMN, 'none');

        % Combinar los indicadores de rechazo

        EEG_MMN = eeg_rejsuperpose( ...
            EEG_MMN, ...
            1, 1, 0, 0, 0, 0, 0, 0);

        mmnRejectMask = logical(EEG_MMN.reject.rejglobal);

        mmn_total_rejection_pct(x) = ...
            100 * sum(mmnRejectMask) / EEG_MMN.trials;

        % Guardar recuentos por condición

        mmn_std_accepted(x) = ...
            sum(acce_MMN([1 3 5 6]));

        mmn_std_rejected(x) = ...
            sum(rej_MMN([1 3 5 6]));

        mmn_dev_accepted(x) = ...
            sum(acce_MMN([2 4]));

        mmn_dev_rejected(x) = ...
            sum(rej_MMN([2 4]));

        mmn_bin2_accepted(x) = acce_MMN(2);
        mmn_bin4_accepted(x) = acce_MMN(4);

        mmnTotalDev = ...
            mmn_dev_accepted(x) + mmn_dev_rejected(x);

        if mmnTotalDev > 0
            mmn_dev_rejection_pct(x) = ...
                100 * mmn_dev_rejected(x) / mmnTotalDev;
        end

        % Eliminar los ensayos rechazados para MMN

        EEG_MMN = pop_rejepoch( ...
            EEG_MMN, ...
            mmnRejectMask, ...
            0);

        EEG_MMN = eeg_checkset(EEG_MMN);

        % Promediar los ensayos aceptados para MMN

        ERP_MMN = pop_averager( ...
            EEG_MMN, ...
            'Criterion', 'good', ...
            'ExcludeBoundary', 'on', ...
            'SEM', 'on');

        % Inicializar matrices grupales después del primer sujeto válido

        if isempty(ERPtotal_PEP_MMN)

            mmnSize = size(ERP_MMN.bindata);

            if length(mmnSize) < 3
                mmnSize(3) = 1;
            end

            ERPtotal_PEP_MMN = nan( ...
                mmnSize(1), ...
                mmnSize(2), ...
                mmnSize(3), ...
                nSubjects);

            ERRORtotal_PEP_MMN = nan( ...
                mmnSize(1), ...
                mmnSize(2), ...
                mmnSize(3), ...
                nSubjects);
        end

        ERPtotal_PEP_MMN(:, :, :, x) = ...
            ERP_MMN.bindata;

        ERRORtotal_PEP_MMN(:, :, :, x) = ...
            ERP_MMN.binerror;

        processing_status(x) = "OK";

        fprintf(['CNV/PINV: %d desviantes aceptados; ' ...
            'MMN: %d desviantes aceptados\n'], ...
            cnvPinv_dev_accepted(x), ...
            mmn_dev_accepted(x));

    catch ME

        processing_status(x) = string(ME.message);

        warning( ...
            'No se pudo procesar %s: %s', ...
            name, ...
            ME.message);
    end
end

% Crear tabla conjunta de control de calidad

results_table = table( ...
    subject_names, ...
    processing_status, ...
    initial_trials, ...
    cnvPinv_std_accepted, ...
    cnvPinv_std_rejected, ...
    cnvPinv_dev_accepted, ...
    cnvPinv_dev_rejected, ...
    cnvPinv_bin2_accepted, ...
    cnvPinv_bin4_accepted, ...
    cnvPinv_total_rejection_pct, ...
    cnvPinv_dev_rejection_pct, ...
    mmn_std_accepted, ...
    mmn_std_rejected, ...
    mmn_dev_accepted, ...
    mmn_dev_rejected, ...
    mmn_bin2_accepted, ...
    mmn_bin4_accepted, ...
    mmn_total_rejection_pct, ...
    mmn_dev_rejection_pct, ...
    'VariableNames', { ...
    'Subject', ...
    'Status', ...
    'InitialTrials', ...
    'CNV_PINV_StdAccepted', ...
    'CNV_PINV_StdRejected', ...
    'CNV_PINV_DevAccepted', ...
    'CNV_PINV_DevRejected', ...
    'CNV_PINV_Bin2Accepted', ...
    'CNV_PINV_Bin4Accepted', ...
    'CNV_PINV_TotalRejectionPct', ...
    'CNV_PINV_DevRejectionPct', ...
    'MMN_StdAccepted', ...
    'MMN_StdRejected', ...
    'MMN_DevAccepted', ...
    'MMN_DevRejected', ...
    'MMN_Bin2Accepted', ...
    'MMN_Bin4Accepted', ...
    'MMN_TotalRejectionPct', ...
    'MMN_DevRejectionPct'});

% Guardar tabla de control de calidad

writetable(results_table, csvOutput);

% Guardar resultados CNV/PINV

save( ...
    fullfile(outputPath, ...
    'ERPtotal_PEP_CNV_PINV_0p05Hz.mat'), ...
    'ERPtotal_PEP_CNV_PINV', ...
    'ERRORtotal_PEP_CNV_PINV', ...
    'subject_names', ...
    'results_table', ...
    'cnvPinvStartMs', ...
    'cnvPinvEndMs', ...
    'cnvPinvBaseline', ...
    '-v7.3');

% Guardar resultados MMN

save( ...
    fullfile(outputPath, ...
    'ERPtotal_PEP_MMN_0p05Hz.mat'), ...
    'ERPtotal_PEP_MMN', ...
    'ERRORtotal_PEP_MMN', ...
    'subject_names', ...
    'results_table', ...
    'mmnStartMs', ...
    'mmnEndMs', ...
    'mmnBaseline', ...
    '-v7.3');

fprintf('\nProcesamiento completado.\n');
fprintf('Tabla de control de calidad: %s\n', csvOutput);

%% 8. ANALISIS DE PINV CON LINEA BASE PRE-S4: AMBOS GRUPOS

% Aplicar despues de completar el preprocesamiento, el rechazo de artefactos
% y el promedio de las epocas largas de controles y pacientes (pasos 1-7).
% La linea base principal de PINV es de -100 a 0 ms respecto a S4.
% Se utilizan las matrices guardadas, sin repetir filtrado, ICA o rechazo
% de ensayos y sin modificar las matrices de origen utilizadas para CNV.

experimentPathPinv = projectPath;

eegPathPinv = fullfile( ...
    experimentPathPinv, ...
    'EEG');

conSourceFile = fullfile( ...
    eegPathPinv, ...
    'ERPtotal_CON_CNV_PINV_0p05Hz.mat');

pepSourceFile = fullfile( ...
    eegPathPinv, ...
    'ERPtotal_PEP_CNV_PINV_0p05Hz.mat');

conPinvFile = fullfile( ...
    eegPathPinv, ...
    ['ERPtotal_CON_PINV_preS4BL_' ...
     '0p05Hz.mat']);

pepPinvFile = fullfile( ...
    eegPathPinv, ...
    ['ERPtotal_PEP_PINV_preS4BL_' ...
     '0p05Hz.mat']);

requiredSourceFiles = {conSourceFile; pepSourceFile};

for fileIndex = 1:numel(requiredSourceFiles)
    if ~isfile(requiredSourceFiles{fileIndex})
        error('No se encuentra la matriz de origen: %s', ...
            requiredSourceFiles{fileIndex});
    end
end

if isfile(conPinvFile) || isfile(pepPinvFile)
    error([ ...
        'Ya existe al menos una matriz del analisis principal de PINV. ' ...
        'No se ha sobrescrito ningun archivo.']);
end

CON_SOURCE = load(conSourceFile);
PEP_SOURCE = load(pepSourceFile);

if ~isfield(CON_SOURCE, 'ERPtotal_CON_CNV_PINV') || ...
        ~isfield(CON_SOURCE, 'subject_names')
    error('La matriz de origen de controles no contiene lo esperado.');
end

if ~isfield(PEP_SOURCE, 'ERPtotal_PEP_CNV_PINV') || ...
        ~isfield(PEP_SOURCE, 'subject_names')
    error('La matriz de origen de pacientes no contiene lo esperado.');
end

ERP_CON_SOURCE_LONG = CON_SOURCE.ERPtotal_CON_CNV_PINV;
ERP_PEP_SOURCE_LONG = PEP_SOURCE.ERPtotal_PEP_CNV_PINV;

if size(ERP_CON_SOURCE_LONG, 2) ~= ...
        size(ERP_PEP_SOURCE_LONG, 2)
    error('Las matrices de origen CON y PEP difieren en longitud.');
end

samplingRateHz = 1000;
timePinvS4 = ...
    -940 + (0:(size(ERP_CON_SOURCE_LONG, 2) - 1)) * ...
    (1000 / samplingRateHz);

pinvBaselineMs = [-100 0];

pinvBaselineMask = ...
    timePinvS4 >= pinvBaselineMs(1) & ...
    timePinvS4 <= pinvBaselineMs(2);

if sum(pinvBaselineMask) ~= 101
    error(['La linea base pre-S4 no contiene las 101 muestras ' ...
        'esperadas entre -100 y 0 ms, ambas incluidas.']);
end

% La transformacion es lineal y conserva exactamente los mismos
% participantes, bins y ensayos aceptados de las matrices de origen.

conBaselineMean = mean( ...
    ERP_CON_SOURCE_LONG(:, pinvBaselineMask, :, :), ...
    2);

pepBaselineMean = mean( ...
    ERP_PEP_SOURCE_LONG(:, pinvBaselineMask, :, :), ...
    2);

ERPtotal_CON_PINV_preS4BL = ...
    ERP_CON_SOURCE_LONG - conBaselineMean;

ERPtotal_PEP_PINV_preS4BL = ...
    ERP_PEP_SOURCE_LONG - pepBaselineMean;

conResidualBaseline = max(abs(mean( ...
    ERPtotal_CON_PINV_preS4BL( ...
    :, pinvBaselineMask, :, :), 2)), [], 'all');

pepResidualBaseline = max(abs(mean( ...
    ERPtotal_PEP_PINV_preS4BL( ...
    :, pinvBaselineMask, :, :), 2)), [], 'all');

if conResidualBaseline > 1e-5 || pepResidualBaseline > 1e-5
    error('La comprobacion de la nueva linea base no fue satisfactoria.');
end

pinvDefinition = [ ...
    'Primary PINV analysis: each subject, electrode and bin was ' ...
    'baseline-corrected using its mean from -100 to 0 ms relative to S4. ' ...
    'Participants, bins and accepted trials are unchanged from the ' ...
    'source long-epoch matrices used for CNV.'];

sourceBaselineMsRelativeToS4 = [-1000 -900];
pinvBaselineMsRelativeToS4 = ...
    pinvBaselineMs;

sourceFile = conSourceFile;
subject_names = CON_SOURCE.subject_names;

save( ...
    conPinvFile, ...
    'ERPtotal_CON_PINV_preS4BL', ...
    'subject_names', ...
    'timePinvS4', ...
    'sourceBaselineMsRelativeToS4', ...
    'pinvBaselineMsRelativeToS4', ...
    'pinvDefinition', ...
    'sourceFile', ...
    '-v7.3');

sourceFile = pepSourceFile;
subject_names = PEP_SOURCE.subject_names;

save( ...
    pepPinvFile, ...
    'ERPtotal_PEP_PINV_preS4BL', ...
    'subject_names', ...
    'timePinvS4', ...
    'sourceBaselineMsRelativeToS4', ...
    'pinvBaselineMsRelativeToS4', ...
    'pinvDefinition', ...
    'sourceFile', ...
    '-v7.3');

fprintf('\nAnalisis principal de PINV pre-S4 generado sin sobrescribir datos.\n');
fprintf('Controles: %s\n', conPinvFile);
fprintf('Pacientes: %s\n', pepPinvFile);
fprintf('Linea base: -100 a 0 ms respecto a S4.\n');
fprintf(['Los participantes, bins y ensayos aceptados son identicos ' ...
    'a los de las matrices de origen.\n']);

%% 9. PLOTEAR PROMEDIOS ENTRE TODOS LOS SUJETOS

% Para generar las figuras, ejecutar plot_figures.m por separado.