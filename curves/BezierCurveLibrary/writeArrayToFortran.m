function fid = writeArrayToFortran(fid, array, arrayName)

fprintf(fid,'\n      %s =\\ \n',arrayName);

if size(array,1) == 1 && size(array,2) > 1
    row=1;
    for col = 1:1:size(array,2)
        if(col==1)
            fprintf(fid,'      (%1.12e,',array(row,col));
        elseif(col < size(array,2))
            fprintf(fid,'%1.12e,',array(row,col));
        else
            fprintf(fid,'%1.12e)\n',array(row,col));
        end
    end
elseif size(array,2) == 1 && size(array,1) > 1
    col=1;
    for row = 1:1:size(array,1)
        if(row==1)
            fprintf(fid,'      (%1.12e,',array(row,col));
        elseif(row < size(array,1))
            fprintf(fid,'%1.12e,',array(row,col));
        else
            fprintf(fid,'%1.12e)\n',array(row,col));
        end
    end
else
    for col = 1:1:size(array,2)
        for row=1:1:size(array,1)
            if( row == 1 && col == 1)
                fprintf(fid,'      (%1.12e,',array(row,col));
            elseif(row == 1 && col >= 1)
                fprintf(fid,'       %1.12e,',array(row,col));                
            elseif(row > 1 && row < size(array,1) && col <= size(array,2))
                fprintf(fid,'%1.12e,',array(row,col));
            elseif(row == size(array,1) && col < size(array,2))
                fprintf(fid,'%1.12e,\\ \n',array(row,col));
            elseif( row == size(array,1) && col == size(array,2) )
                fprintf(fid,'%1.12e)\n',array(row,col));  
            else
                assert(0);
            end
        end
    end
end