function [actx_excel,excel_handle]=StartExcel(file_p,Visible)
% Start an ActiveX session with Excel:
actx_excel = actxserver('Excel.Application');
actx_excel.Visible = Visible;
trace(actx_excel.Visible);
if ~exist(file_p,'file')
    % Create new document:
    %excel_handle = invoke(actx_excel.Documents,'Add');
    
else
    % Open existing document:
    excel_handle=actx_excel.Workbooks.Open(file_p);
end
