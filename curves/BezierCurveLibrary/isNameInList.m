function flag = isNameInList(name,cellList)

flag=0;
i=1;
while i <= length(cellList) && flag==0
    if strcmp(cellList{i},name)==1
        flag=1;
    end
    i=i+1;
end
