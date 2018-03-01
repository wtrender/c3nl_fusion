function [varargout]=deal(M,type)
if nargin<2;type='col';end
if strcmpi(type,'row'); M = M';end
varargout = cell(size(M,2));   
for k = 1:size(M,2)
    varargout{k} = M(:,k);
end

end