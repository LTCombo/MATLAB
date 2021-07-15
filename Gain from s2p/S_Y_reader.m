%Preconditioning
%clear all
%close all
clc

%Flags
multiple=0;         %0 for single extraction, 1 for extraction of all the files from the same folder

[filename1, pathname] = uigetfile('*.s1p; *.s2p','Select the S-parameter file');

switch multiple

    %SINGLE READ
    case 0
        
        filename = strcat(pathname, filename1);
        data = read(rfdata.data, filename);
        extension=extractAfter(filename1, length(filename1)-3);
          
        if(extension=='s1p')
            
            % Read freq and S-parameters
            freq = data.freq;
            om=(2*pi).*freq;
            s_params = extract(data, 'S_PARAMETERS',50);
            s11 = squeeze(s_params(1,1,:));
            
            % Read Y-parameters
            y_params = s2y(s_params, 50);
            y11 = squeeze(y_params(1,1,:));
            
        else
           
            % Read freq and S-parameters
            freq = data.freq;
            om=(2*pi).*freq;

            s_params = extract(data, 'S_PARAMETERS',50);
            s11 = squeeze(s_params(1,1,:));
            s12 = squeeze(s_params(1,2,:));
            s21 = squeeze(s_params(2,1,:));
            s22 = squeeze(s_params(2,2,:));

            % Read Y-parameters
            y_params = s2y(s_params, 50);
            y11 = squeeze(y_params(1,1,:));
            y12 = squeeze(y_params(1,2,:));
            y21 = squeeze(y_params(2,1,:));
            y22 = squeeze(y_params(2,2,:));
            
        end

    %MULTIPLE READ
    case 1
        
        %Extension extraction
        extension=extractAfter(filename1, length(filename1)-3);
        
        %Listing
        
        listing = dir(pathname);
        %mkdir(pathname,'Figures\');
        list_length = size(listing,1);
        
        j=1;

        for i=2:list_length

            names_tmp(i-1) = string(listing(i).name);

                if(contains(names_tmp(i-1),extension)==1)
                    names(j) = names_tmp(i-1);
                    j=j+1;
                end

        end
        
        %Reading
        
        if(extension=='s1p')
             for i=1:length(names)

                filename = strcat(pathname, names(i));
                data = read(rfdata.data, char(filename));

                % read freq and s-parameters
                freq = data.freq;
                om=(2*pi).*freq;
                df = freq(2) - freq(1);

                s_params = extract(data, 'S_PARAMETERS',50);
                s11 = squeeze(s_params(1,1,:));

                % for reading y11, y12, y21, y22 from y-parameter
                y_params = s2y(s_params, 50);
                y11 = squeeze(y_params(1,1,:));
            
            end
        else
            for i=1:length(names)

                filename = strcat(pathname, names(i));
                data = read(rfdata.data, char(filename));

                % read freq and s-parameters
                freq = data.freq;
                om=(2*pi).*freq;
                df = freq(2) - freq(1);

                s_params = extract(data, 'S_PARAMETERS',50);
                s11 = squeeze(s_params(1,1,:));
                s12 = squeeze(s_params(1,2,:));
                s21 = squeeze(s_params(2,1,:));
                s22 = squeeze(s_params(2,2,:));

                % for reading y11, y12, y21, y22 from y-parameter
                y_params = s2y(s_params, 50);
                y11 = squeeze(y_params(1,1,:));
                y12 = squeeze(y_params(1,2,:));
                y21 = squeeze(y_params(2,1,:));
                y22 = squeeze(y_params(2,2,:));
            
            end
            
        end
        
        
end

