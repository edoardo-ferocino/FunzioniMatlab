function [SH]=GoToPosExcel(EH,SheetNum)
SH = EH.Sheets.Item(SheetNum);
SH.Select;
end
