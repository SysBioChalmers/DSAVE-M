classdef ProgrBar < handle
    properties (Access = private)
        id % id of this progress bar, used to identify children in the parent progress bar
        progress % number between 0 and 1
        parent % another ProgrBar or []
        text % The text printed out before the progress
        lastPrintLength % number of chars that will be backspaced before new print
        children % map of id and fraction
        done % avoid printing done more than once if called twice for some reason
        silent % true if no printouts should be made
    end
    methods(Static)
      function ctxt = GetSilentContext
         ctxt = ProgrBarContext(true, [], []);
      end
    end
    methods %public methods
        function obj = ProgrBar(text_, context)
            persistent nextId;
            if isempty(nextId)
                nextId = 1;
            end
            obj@handle();%call super class, not sure if that is needed
            obj.id = nextId;
            nextId = nextId + 1;
            if nargin < 2 || isempty(context)
                obj.parent = [];
                obj.silent = false;
            else
                %Check that the fraction is ok
                %need to take care of the case where a context with no
                %parent is sent in as well
                if ~isempty(context.parent) && (isempty(context.fraction) || context.fraction < 0 || context.fraction > 1.1)
                    error('ProgrBar: If a parent is supplied, a valid fraction must be supplied as well')
                end
                obj.parent = context.parent;
                if ~isempty(obj.parent)
                    obj.parent.RegisterChild(obj.id, context.fraction);
                end
                
                if isempty(context.silent)
                    obj.silent = false;
                else
                    obj.silent = context.silent;
                end
            end
            if nargin == 0
                obj.text = 'Progress';
            else
                obj.text = text_;
            end
            obj.lastPrintLength = 0;
            obj.children = containers.Map('KeyType','uint32','ValueType','any');%will use ProgrBarChild
            obj.progress = 0;
            obj.done = false;
            
            %print first phrase (before the progress in itself) in case of
            %no parent
            if isempty(obj.parent) && ~obj.silent
                fprintf('%s: ',obj.text);
            end
        end
        
        function ctxt = GetSubContext(this, fraction)
            ctxt = ProgrBarContext(false, this, fraction);
        end
        
        function Progress(this, val) % 0 <= val <= 1
            this.progress = val;
            this.ReportProgress();
        end
        
        function Done(this)
            if ~this.done % make sure done is only called once; looks messy if it is
                if ~isempty(this.parent)
                    %if we have a parent, don't print any output but instead
                    %send the progress to the parent progress bar
                    this.parent.ProgressFromChild(this.id, 1);
                elseif ~this.silent
                    this.PrintProgress(1);
                    fprintf(' Done\n');
                end

                this.done = true;
            end
        end

    end
    methods (Access = private)
        
        function ProgressFromChild(this, id, progress) %logical vector
            child = this.children(id);
            child.progress = progress;
            this.children(id) = child;
            this.ReportProgress();
        end
        
        function RegisterChild(this, id, fraction)
            s = ProgrBarChild();
            s.progress = 0;
            s.fraction = fraction;
            this.children(id) = s;
            
            %make a sanity check that the sum of all children fractions do not exceed 1
            ss = 0;
            values = this.children.values;
            for i = 1:size(this.children.values,2)
                ss = ss + values{1,i}.fraction;
            end
            
            if (ss > 1.1) % leave some slack due to roundoff errors
                error(strcat('ProgrBar/RegisterChild - The sum of the children fractions exceed 1. Value=', num2str(ss)));
            end
            
        end
        
        %don't call this directly, use Progress
        function ReportProgress(this) % 0 <= val <= 1
            %calculate total progress from children and progress val
            progr = 0;
            fractionSum = 0;
            values = this.children.values;
            for i = 1:size(this.children.values,2)
                progr = progr + values{1,i}.progress * values{1,i}.fraction;
                fractionSum = fractionSum + values{1,i}.fraction;
            end
            progr = progr + this.progress * (1-fractionSum); %fills up the space not taken by children
            
            if ~isempty(this.parent)
                %if we have a parent, don't print any output but instead
                %send the progress to the parent progress bar
                this.parent.ProgressFromChild(this.id, progr);
            else
                this.PrintProgress(progr);
            end
        end
        
        function PrintProgress(this, progr) % 0 <= progr <= 1
            if ~this.silent
                %constants for customizing the printouts
                percentageLength = 10;   %   Length of percentage string
                dotLength        = 10;   %   The total number of dots in a progress bar

                %first remove old progress:
                backspaces = '';
                if (this.lastPrintLength > 0)
                    backspaces = repmat('\b', 1, this.lastPrintLength);
                end
                perc = [num2str(round(progr*100)) '%%'];%percentage
                perc = [perc repmat(' ', 1, percentageLength - length(perc) -1 )];%always tab to percentageLength before adding bar, -1 is for %% -> %
                nDots = floor(progr*dotLength);
                dots = ['[' repmat('.',1,nDots) repmat(' ', 1, dotLength - nDots) ']'];
                text = [' ' perc dots];

                fprintf([backspaces text]);

                this.lastPrintLength = length(text) - 1; % -1 for %% -> %
            end            
        end        
    end
end
    