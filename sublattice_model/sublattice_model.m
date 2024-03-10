classdef sublattice_model
    %SUBLATTICE_MODEL Class for calculation of sublattice molar fractions
    %for sublattice model
    %   (todo) Detailed explanation goes here
    
    properties
        components
        sublattices
        capacities
    end
    
    methods
        function obj = sublattice_model(varargin)
            %SUBLATTICE_MODEL Initializes a sublattice description of a
            %phase
            %   Accepts a cell array of strings representing components of
            %   the system
            obj.components = varargin;
            obj.sublattices = {};
            obj.capacities = [];
        end
        
        function obj = addSublattice(obj, capacity, varargin)
            %ADDSUBLATTICE Adds a sublattice to the phase
            %   Accepts capacity and pairs of components and their
            %   charges (vacancy denoted by an empty string), e.g.
            %   'A', 2, 'B', 3, '', 0
            comps = varargin;
            sublattice = {};
            lastid = 0;
            if size(obj.sublattices, 2) > 0
                lastsubl = obj.sublattices{size(obj.sublattices, 2)};
                lastid = lastsubl{size(lastsubl, 2)}.y_id;
            end
            if rem(size(comps), 2) == 1
                error('Odd number of component arguments, possibly wrong format');
            end
            for i = 1:(size(comps,2)/2)
                i1 = (i - 1)*2 + 1;
                i2 = (i - 1)*2 + 2;
                id = 0;
                comp_name = comps{i1};
                %b = size(find(strcmp(obj.components, comp_name)),2)
                if size(find(strcmp(obj.components, comp_name)),2)==0 && not(isempty(comp_name))
                    error('Component %s not found in components list', comp_name);
                end
                if not(isempty(comp_name))
                    id = find(strcmp(obj.components, comp_name));
                end
                charge = comps{i2};
                lastid = lastid + 1;
                comp_i = struct('id', id, 'charge', charge, 'y_id', lastid);
                comp_i.name = sprintf('%s_(%d)', obj.makeName(comp_i), size(obj.sublattices, 2) + 1);
                sublattice = [sublattice comp_i];
            end
            obj.sublattices = [obj.sublattices {sublattice}];
            obj.capacities = [obj.capacities capacity];
        end
        
        function [comps, y_funcs] = calculateModel(obj, varargin)
            matr = sym([]);
            rhs = sym([]);
            row = 1;
            % conditions on sublattices
            for s = 1:size(obj.sublattices, 2)
                subl = obj.sublattices{s};
                for sp = 1:size(subl, 2)
                    matr(row, subl{sp}.y_id) = 1;
                end
                rhs(row, 1) = 1;
                row = row + 1;
            end
            % electroneutrality
            for s = 1:size(obj.sublattices, 2)
                subl = obj.sublattices{s};
                for sp = 1:size(subl, 2)
                    matr(row, subl{sp}.y_id) = subl{sp}.charge * obj.capacities(s);
                end
            end
            rhs(row, 1) = 0;
            row = row + 1;
            % parsing input free parameter args
            if rem(size(varargin, 2), 2) == 1
                error('Odd number of elements in input, expected pairs of component name and symbolic variable');
            end
            for i = 1:size(varargin, 2) / 2
                i_name = (i - 1)*2 + 1;
                i_sym = (i - 1)*2 + 2;
                name = varargin{i_name};
                x = varargin{i_sym};
                in_comps = find(strcmp(obj.components, name));
                in_subl = find(strcmp(obj.sublatticeNames()', name));
                if size(in_comps, 2) ~= 0
                    % a component
                    id = in_comps(1);
                    for s = 1:size(obj.sublattices, 2)
                        subl = obj.sublattices{s};
                        for sp = 1:size(subl, 2)
                            if subl{sp}.id == id
                                matr(row, subl{sp}.y_id) = (x - 1) * obj.capacities(s);
                            elseif subl{sp}.id ~= 0
                                matr(row, subl{sp}.y_id) = x * obj.capacities(s);
                            end
                        end
                    end
                    rhs(row, 1) = 0;
                elseif size(in_subl, 2) ~= 0
                    % a component in sublattice
                    y_id = in_subl(1);
                    matr(row, y_id) = 1;
                    rhs(row, 1) = x;
                else
                    error('Component "%s" not found among components or sublattices', name);
                end
                row = row + 1;
            end
            comps = obj.sublatticeNames();
            y_funcs = matr\rhs;
        end
        
        function names = sublatticeNames(obj)
            names = string([]);
            for s = 1:size(obj.sublattices, 2)
                subl = obj.sublattices{s};
                for sp = 1:size(subl, 2)
                    names = [names; sprintf('%s_(%d)', obj.makeName(subl{sp}), s)];
                end
            end
        end
        
        function formula = printFormula(obj)
            formula = '';
            for s = 1:size(obj.sublattices, 2)
                if s ~= 1
                    formula = [formula ':'];
                end
                formula = [formula '('];
                subl = obj.sublattices{s};
                for sp = 1:size(subl, 2)
                    if sp ~= 1
                        formula = [formula ','];
                    end
                    formula = [formula obj.makeName(subl{sp})];
                end
                formula = [formula sprintf(')%d', obj.capacities(s))];
            end
        end
        
        function name = makeName(obj, spec)
            if spec.id == 0
                name = 'Va';
            else
                name = obj.components{spec.id};
            end
            if spec.charge ~= 0
                sign = '+';
                if spec.charge < 0
                    sign = '-';
                end
                name = sprintf('%s^%d%s', name, abs(spec.charge), sign);
            end
        end
    end
end