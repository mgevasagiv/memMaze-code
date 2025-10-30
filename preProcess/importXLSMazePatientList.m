function List = importXLSMazePatientList(excel_sheet,sheet,cell_rows)

if strcmp(sheet,'Subjects')
p_in.numeric_fields = {'exp','X','N','D','Z','G'}; % numeric fields in excel
else
    p_in.numeric_fields = {'Time_onset','Time_offset','Block','Trial','Maze','Repetition',...
        'Test_contextual_success','Test_non_contextual_success','Test_contextual_distance',...
        'Test_non_contextual_distance'}; % numeric fields in excel
end

if ~exist('cell_rows') || isempty(cell_rows)
    cell_rows = 1:9999;
else
    cell_rows = 1:cell_rows;
end

List = read_excel_sheet(excel_sheet,sheet,cell_rows,p_in.numeric_fields,p_in);

end % func