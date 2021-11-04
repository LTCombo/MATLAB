function filename = getFileGUI(type)
%Returns the path of the selected file
%Requires the file extension in the format '.format' (i.e. '.txt')

    dialog = strcat('Select file',type);
    [filename1, pathname] = uigetfile(dialog,'Select the file to load');
    filename = strcat(pathname, filename1);

end