
function init(obj)
  
    answer = questdlg('Choose .mph loading approach','Selection','Automatic','Manual','Automatic');

    switch answer

        case 'Automatic'

            files=dir;

            j=0;

            goodindex=[];

            for i=1:length(files)

                if endsWith(files(i).name,'.mph')

                    j=j+1;

                    goodindex(j)=i;

                end

            end

            if length(goodindex)==1

                goodfile = files(goodindex);

                obj.model = mphload([goodfile.folder,filesep,goodfile.name]);
                obj.tag = strrep(goodfile.name,'.mph','');
                obj.save_folder = goodfile.folder;

                fprintf(sprintf('Model %s initialized...\n',goodfile.name));
                fprintf(sprintf('Tag %s initialized...\n',obj.tag));
                fprintf(sprintf('Save folder initialized...\n',obj.save_folder));

            else

                fprintf('No model initialized...\n');
                fprintf('No tag initialized...\n');
                fprintf('No folder initialized...\n');

            end

        case 'Manual'
            
            [filename, path] = uigetfile('*.mph');
            
                if isequal(filename,0)
                   
                  disp('User selected Cancel');
                  fprintf('No model initialized...\n');
                  fprintf('No tag initialized...\n');
                  fprintf('No folder initialized...\n');

                else
                   
                  obj.model = mphload(strcat(path,filename));
                  obj.tag = strrep(filename,'.mph','');
                  obj.save_folder = path;
                   
                  fprintf(sprintf('Model %s initialized...\n',filename));
                  fprintf(sprintf('Tag %s initialized...\n',obj.tag));
                  fprintf(sprintf('Save folder initialized...\n',obj.save_folder));
                   
                end

     end
    
end
        
        