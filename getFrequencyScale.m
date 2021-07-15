function [k, scale] = getFrequencyScale(f)

    exponent = ceil(log10(f(end)+1));
    index = ceil(exponent/3);

    switch index
        
        case 1
            k = 1;
            scale = "[Hz]";
        case 2
            k = 1e3;
            scale = "[kHz]";
        case 3
            k = 1e6;
            scale = "[MHz]";
        case 4
            k = 1e9;
            scale = "[GHz]";
        case 5
            k = 1e12;
            scale = "[THz]";
            
    end

end