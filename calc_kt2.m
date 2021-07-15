function [kt2_1, kt2_2, kt2_3] = calckt2(fs, fp)
    %Calculates kt2 according to different definitions
    %kt2_1 is the most conservative
    
    keff2 = (fp.^2-fs.^2)./fp.^2;
    K2 = (fp.^2-fs.^2)./fs.^2;
    kt2_1 = pi/2*fs./fp*1./tan(pi/2*fs./fp);
    kt2_2 = pi^2/4*(1-fs./fp);
    kt2_3 = pi^2/8.*K2;

end