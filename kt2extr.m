function [kt2] = kt2extr(fs, fp)

    kt2 = pi/2*fs/fp*1/tan(pi/2*fs/fp);

end