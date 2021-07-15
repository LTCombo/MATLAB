function fig = linplot(type, f, X, varargin)
    %Plots in linear
    %Requires type, f, X as inputs
    %X can be 'Y', 'Z', or 'G' (all base 20) - If incorrect, 'Y' is automatically
    %selected
    %Varargin can be either 'Magnitude' or 'Phase', if not selected both
    %are plotted

    [k, index] = getFrequencyScale(f(end)); 
    XPlot = abs(X);
    fPlot = f/k;
    frequencyLabel = strcat('Frequency, {\itf}', " ",index);
   
    switch type
        case 'Y'
            variableLabel = 'Admittance, {\itY} [S]';
            deltaXp = 10;
            deltaXm = 0;
        case 'Z'
            variableLabel = 'Impedance, {\itZ} [\Omega]';
            deltaXp = 10;
            deltaXm = 0;
        case 'G'
            variableLabel = 'Gain, {\itG}_v [V/V]';
            deltaXp = 10;
            deltaXm = 0;
        otherwise
            variableLabel = 'Admittance, {\itY} [S]';
            deltaXp = 10;
            deltaXm = 0;
    end
        
    
    if isempty(varargin) == 0
       
        switch char(varargin)
            
            case 'Magnitude'
            
                fig = figure;
                plot(fPlot, XPlot, 'LineWidth', 3)
                xlabel(frequencyLabel)
                ylabel(variableLabel)
                plotPreferences
                axis([fPlot(1) fPlot(end) min(XPlot, [], 'all')-deltaXm max(XPlot, [], 'all')+deltaXp])

            case 'Phase'
                
                fig = figure;
                plot(fPlot, angle(X)*180/pi, 'LineWidth', 3)
                xlabel(frequencyLabel)
                ylabel('Phase, \theta [\circ]')
                plotPreferences           
                axis([fPlot(1) fPlot(end) -100 100])
               
        end
    
    else
        
            fig = figure;
            subplot(2,1,1)
            plot(fPlot, 20*log10(abs(X)), 'LineWidth', 3)
            xlabel(frequencyLabel)
            ylabel(variableLabel)
            plotPreferences
            axis([fPlot(1) fPlot(end) min(XPlot, [], 'all')-deltaXm max(XPlot, [], 'all')+deltaXp])

            subplot(2,1,2)
            plot(fPlot, angle(X)*180/pi, 'LineWidth', 3)
            xlabel(frequencyLabel)
            ylabel('Phase, \theta  [\circ]')
            plotPreferences           
            axis([fPlot(1) fPlot(end) -100 100])
        
    end
    
end