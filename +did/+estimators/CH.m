% ============================================================
% File: +did/+estimators/CH.m
% Purpose: OOP wrapper for de Chaisemartin & d'Haultfoeuille (2020) DID_M
%          compatible with did.Model (subclass of did.estimators.Estimator).
%
% Usage (factory route):
%   est = did.estimators.CH(B=399, ComputePlacebo=true, Seed=42);
%   ds  = did.Dataset.fromTable(T, "id","time","y","D");
%   mdl = did.Model(ds, est);
%   out = mdl.fit();
%
% Notes:
%   • Heavy lifting is delegated to did.ch_estimator (functional core).
%   • We provide 'summaryTable' (rich: SE, tStat, pValue) so Model can print it.
%   • We also provide a Model-friendly 'coef' table.
% ============================================================
classdef CH < did.estimators.Estimator

    properties
        % ---- Estimator options ----
        B (1,1) double {mustBeInteger, mustBeNonnegative} = 100
        ComputePlacebo (1,1) logical = true
        Seed (1,1) double = randi([1,1e7],1,1)
        WeightVar (1,1) string = ""     % optional weight variable in the input table
        Print (1,1) logical = true
        Details (1,1) logical = false
        Covariates string = string.empty(1,0)
        CovarSample (1,1) string {mustBeMember(CovarSample,["D0","never","all"])} = "D0"
    
    end

    methods
        function obj = CH(varargin)
            % Constructor: accepts Name–Value pairs. Common synonyms mapped.
            if nargin > 0
                if mod(nargin,2) ~= 0
                    error('did:estimators:CH:InvalidArgs', ...
                        'Constructor expects Name–Value pairs.');
                end
                for k = 1:2:nargin
                    name = string(varargin{k});
                    val  = varargin{k+1};
                    switch lower(name)
                        case 'display',        name = "Print";
                        % case 'rngseed',        name = "Seed";
                        case 'computeplacebo', name = "ComputePlacebo";
                        case 'weight',         name = "WeightVar";
                        case 'covariates',  name = "Covariates";
                        case 'x',           name = "Covariates";
                        case 'covarsample', name = "CovarSample";
                    end
                    if isprop(obj, name)
                        obj.(name) = val;
                    end
                end
            end
            % Set a friendly name on the inherited property (if present)
            try
                obj.Name = "CH2020";
            catch
                % superclass may not expose Name; ignore
            end
        end

        function out = fit(obj, ds)
            % FIT  Run DID_M via the functional core on the dataset provided by did.Model

            % -- extract data + var names from Dataset
            T = obj.getTable_(ds);
            [idVar, timeVar, yVar, dVar] = obj.getVarNames_(ds);

            % -- build NV args for functional core
            nv = {'idVar', idVar, 'timeVar', timeVar, 'yVar', yVar, 'dVar', dVar, ...
                'B', obj.B, 'ComputePlacebo', obj.ComputePlacebo, 'Seed', obj.Seed, ...
                'Print', obj.Print, 'Details', obj.Details, ...
                'Covariates', obj.Covariates, 'CovarSample', obj.CovarSample};

            if strlength(obj.WeightVar) > 0 && any(strcmp(obj.WeightVar, T.Properties.VariableNames))
                nv = [nv, {'WeightVar', obj.WeightVar}]; 
            end

            % -- delegate to functional core
            out = did.ch_estimator(T, nv{:});

            % -- standardize metadata
            if ~isfield(out, 'Estimator') || isempty(out.Estimator)
                try
                    out.Estimator = obj.Name;  % inherited from Estimator base
                catch
                    out.Estimator = "CH2020";
                end
            end
            if ~isfield(out, 'Vcov') || isempty(out.Vcov)
                out.Vcov = struct('Type','cluster-bootstrap', ...
                    'Clusters', idVar, ...
                    'B', obj.B, ...
                    'Seed', obj.Seed);
            end
            
            out.Method = "didM/CH";

            % -- provide a rich summaryTable the Model can use directly
            % out.Summary has: Effect, Estimate, SE, t, p
            sumT = out.summaryTable;
            if ~ismember('Name', sumT.Properties.VariableNames)
                sumT.Name = string(sumT.Effect);
            end
            if ismember('t', sumT.Properties.VariableNames) && ~ismember('tStat', sumT.Properties.VariableNames)
                sumT.tStat = sumT.t;
            end
            if ismember('p', sumT.Properties.VariableNames) && ~ismember('pValue', sumT.Properties.VariableNames)
                sumT.pValue = sumT.p;
            end
            sumT = movevars(sumT, 'Name', 'Before', 1);
            wantSummary = {'Name','Effect','Estimate','SE','tStat','pValue'};
            haveSummary = wantSummary(ismember(wantSummary, sumT.Properties.VariableNames));
            out.summaryTable = sumT(:, haveSummary);

            % -- also provide coef for utilities that expect it (Model fallback)
            co = out.summaryTable;
            wantCoef = {'Name','Estimate','SE','tStat','pValue'};
            haveCoef = wantCoef(ismember(wantCoef, co.Properties.VariableNames));
            out.coef = co(:, haveCoef);
        end
    end

    methods (Access = private)
        function T = getTable_(~, ds)
            % Try common Dataset representations
            if isprop(ds, 'Tbl')
                T = ds.Tbl;
            elseif isprop(ds, 'T')
                T = ds.T;
            elseif isprop(ds, 'Table')
                T = ds.Table;
            else
                error('did:estimators:CH:NoTable', ...
                    'Cannot locate table inside Dataset object.');
            end
        end

        function [idVar, timeVar, yVar, dVar] = getVarNames_(~, ds)
        % Use Dataset metadata directly (no guessing, no defaults)
        if isprop(ds, "idVar"),   idVar   = ds.idVar;   else, error("Dataset missing idVar");   end
        if isprop(ds, "timeVar"), timeVar = ds.timeVar; else, error("Dataset missing timeVar"); end
        if isprop(ds, "yVar"),    yVar    = ds.yVar;    else, error("Dataset missing yVar");    end
        if isprop(ds, "dVar"),    dVar    = ds.dVar;    else, error("Dataset missing dVar");    end
    end
    end
end
% ============================================================
