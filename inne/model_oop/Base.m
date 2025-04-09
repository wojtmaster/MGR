classdef Base
    methods(Abstract)
        obj = fuzzyfication(obj);
        Y = fuzzy_system(obj, params, arg);
        model(obj, a, b, u, y, K, kk, tau);
        Y = find_value(obj, arg, index)
        show_fuzzy_system(obj)
        show_static_characteristic(obj)
    end
end