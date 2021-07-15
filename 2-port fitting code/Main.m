clear all
close all

[filename1, pathname] = uigetfile('Select file.s2p','Select the S-parameter file');

listing = dir(pathname);
mkdir(pathname,'Figures\');

list_length = size(listing,1);
j=1;

for i=2:list_length

names_tmp(i-1) = string(listing(i).name);

    if(contains(names_tmp(i-1),'s2p')==1)
        names(j) = names_tmp(i-1);
        j=j+1;
    end

end

for i=1:length(names)

    filename = strcat(pathname, names(i));
    [Q_3dB(i),Q_slope(i),kt2_a(i),fs(i),Rm(i),Cm(i),Lm(i),C0(i)]=fit_resonator_cld(pathname,filename,names(i));
    
end

namesWoExt=char(erase(names,'.s1p'));

processData.names=names';
processData.fs=fs';
processData.kt=kt2_a';
processData.Q_3dB=Q_3dB';
processData.Q_slope=Q_slope';
processData.Rm=Rm';
processData.Cm=Cm';
processData.Lm=Lm';
processData.C0=C0';

processDataTable = struct2table(processData);
writetable(processDataTable,'ResonatorData.csv');