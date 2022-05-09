%%to compile all these stats mat into the Wildtype_stats_cell.mat
for i = [2,3,4,5,7,8,9,11,12,13,14,15,16]
    path_to_file = "/Users/mikahlbk/Google Drive File Stream/Shared drives/Division_Paper_Final/Experimental_Images/Wildtype/Stats/Stats_With_Layer/";
    filename = "stats_" + num2str(i) + ".mat";
    stats_temp = load(path_to_file + filename);
    temp_field_cell = fields(stats_temp);
    temp_field = temp_field_cell{1};
    temp_table = getfield(stats_temp,temp_field);
    Stats_cell{1,i} = temp_table;
end
