function tf = isEmptyTable(T)
tf = isempty(T) || (istable(T) && height(T)==0);
end
