function [Q3dB, BW] = extractQ3dB(YdB, freq)

    df = freq(2) - freq(1);

    [YdBMax, posMax] = max(YdB);
    Y3dB = YdBMax + 10*log10(0.5);
    fs = freq(posMax);

    posBW = find(YdB>Y3dB);
    posBW2 = [(posBW(1)-1), posBW, (posBW(end)+1)];

    pRise = polyfit([0; 1], YdB(posBW2(1:2)), 1);
    pFall = polyfit([0; 1], YdB(posBW2(end-1:end)), 1);

    Yrise = @(f) pRise(1)*f + pRise(2);
    Yfall = @(f) pFall(1)*f + pFall(2);

    Yrise0 = @(f) pRise(1)*f + pRise(2) - Y3dB;
    Yfall0 = @(f) pFall(1)*f + pFall(2) - Y3dB;

    xRise = fzero(Yrise0, 0.5);
    xFall = fzero(Yfall0, 0.5);

    BWn = length(posBW2)-1+xFall-xRise;
    BW = BWn*df;
    Q3dB = fs/BW;

    space = linspace(0,1,100);
    frise = linspace(1,2,100);
    ffall = linspace(length(posBW2)-1,length(posBW2),100);

    dBv = Y3dB*ones(length(posBW2),1);

    %Figure testing
%     figure(1)
%     plot(YdB(posBW2),'b')
%     hold on
%     plot(YdB(posBW2),'ko')
% 
%     plot(dBv)
% 
%     plot(frise,Yrise(space),'r')
%     plot(ffall,Yfall(space),'r')
% 
%     plot(1+xRise,Yrise(xRise),'ro')
%     plot(length(posBW2)-1+xFall,Yfall(xFall),'ro')
%     legend('Y_{dB}','Y_{dB} (points)','3 dB threshold','Rise (fitting)','Fall (fitting)','Rise (intpl)','Fall (intpl)')
% 
%     xlabel('Points')
%     ylabel('Admittance, Y_{12} [dB]')
%     set(gcf,'color','white')

end