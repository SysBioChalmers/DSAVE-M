classdef ProgrBarContext
    properties
        silent % if true no output will be made
        parent % parent progress bar
        fraction % fraction of the parent's progress bar
    end
    methods
        function obj = ProgrBarContext(silent_, parent_, fraction_)
            if nargin < 1
                silent_ = false;
            end
            if nargin < 2
                parent_ = [];
            end
            if nargin < 3
                fraction_ = 0;
            end
            
            obj.silent = silent_;
            obj.parent = parent_;
            obj.fraction = fraction_;
        end
    end
end
    