function [freq, S, Y] = extractSY(filename, varargin)
    %varargin represents the source Z. If not specified is automatically
    %set to 50 Ohm.

    if isempty(varargin) == 1
    
        Z0 = 50;
    
    else
        
        Z0 = varargin;
        
    end

    data = read(rfdata.data, filename);
    freq = data.freq;
    S = extract(data, 'S_PARAMETERS', Z0);
    Y = s2y(S, Z0);
    
%     %Extract S-parameters
%     s11 = squeeze(S(1,1,:));
%     s12 = squeeze(S(1,2,:));
%     s21 = squeeze(S(2,1,:));
%     s22 = squeeze(S(2,2,:));
% 
%     %Convert to Y-parameters
%     y11 = squeeze(Y(1,1,:));
%     y12 = squeeze(Y(1,2,:));
%     y21 = squeeze(Y(2,1,:));
%     y22 = squeeze(Y(2,2,:));
    
    
end