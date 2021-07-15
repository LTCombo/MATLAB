function [Y,freq]=COMSOL_reader()

[filename1, pathname] = uigetfile('*.txt','Select the txt file for wide-span COMSOL simulation');
filename = strcat(pathname, filename1);

fid = fopen(filename);                      %Open txt files

tline = fgetl(fid);                         %Get first line
tmp_v=[];                                   %Pre-allocate tmp variable
i=1;                                        %Set indexes 
j=1;

while ischar(tline)                         %Checks if there is a character

    tmp_cell=cellstr(tline);                
    a(i)=tmp_cell;
    
    if(tline(1)~='%')
    b(j)=cellstr(tline);
    j=j+1;
    end
    
    tline = fgetl(fid);
    i=i+1;
    
end

fclose(fid);

k=0;
j=0;

for i=1:length(b)
   
    if(mod(i,2)==0)
        k=k+1;
        Y(k)=str2double(b(i));
    else
        j=j+1;
        freq(j)=str2double(b(i));
    end
    
end

end