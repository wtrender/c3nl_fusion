function mdl2Table(mdl,fn,output,caption,label,resize)
if ~exist('resize','var');resize=0;end

if isfield(mdl.Properties,'VarNames')
    fin =mdl.Properties.VarNames;
else 
    fin = mdl.Properties.VariableNames;
end
data = table();
for ii=1:numel(fin)
    tmp  = mdl.(fin{ii});
    data.(fin{ii}) =  mdl.(fin{ii});
    switch class(tmp)
        case 'categorical';format{ii} = '%s';
        case 'cell';format{ii} = '%s';
        case 'char';format{ii} = '%s';
        case 'double' 
            if any(mod(tmp,1))
               format{ii} = '%0.4f';
            else 
               format{ii} = '%i';
            end
        case {'uint16','uint8'};format{ii} = '%i';
    end
    
end
rows = {size(data,1)};
for ii=1:size(data,1)
row ='';
for k=1:numel(format)
    tmp = data{ii,k};
    if iscell(tmp)
        tmp=tmp{1};
        if ~isempty(tmp)
            tmp(strfind(tmp,'_'))=' ';
        else
            tmp = ' ';
        end
    end
    row = [row,sprintf([format{k} ' & '],tmp)];
end
if ii~=size(data,1);row=[row(1:end-2), '\\\midrule'];else; row= [row(1:end-2), '\\\bottomrule[0.2em]'];end
rows{ii,1} = row;
end
toprow = [strjoin(data.Properties.VariableNames,' & ') '\\\toprule[0.2em]'];
toprow(strfind(toprow,'_'))=' ';
rows = [toprow;rows];
rows = strrep(rows,'NaN',' ');
rows{end+1} = '\end{tabular}';
switch  output
    case 'latex'
        top = {
            '\begin{table}';...
            '\centering'};
        if resize
            top{end+1}='\resizebox{\textwidth}{!}{%';
            rows{end}(end+1) = '}';
        end
        bottom = {'\end{table}'};
        caption = {['\caption{' caption '\label{tabel:' label '}}']};
        header = {['\begin{tabular}[0.2em]','{@{}l',repmat('l',1,size(data,2)),'@{}}\toprule']};
        outputtable = [top;header;rows;caption;bottom];
        fid=fopen(fn,'w');
        [nrows,~] = size(outputtable);
        for row = 1:nrows
            fprintf(fid,'%s\n',outputtable{row});
        end
        fclose(fid);
    case 'csv'
        writetable(data,fn);
end
        
end