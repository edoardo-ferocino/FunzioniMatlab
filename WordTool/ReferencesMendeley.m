%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Parameters                                %
clear all
clc
RefParagID = 1;
RefParagPOS = 2;
RefParagStyle = 3;
RefParagCheckRef = 4;
RefStyle = 5;
RefMultiStyle = 6;
RefMultipleStyle = 7;
FForwardT = true;
FForwardF = false;
Eof = 0;
Bof = 1;

NumRef = 20;
Ref(RefParagID).V = 'References and links';
Ref(RefParagPOS).V = Bof;
Ref(RefParagStyle).V = '%d.';
Ref(RefParagCheckRef).V = '';
Ref(RefStyle).Bra = '[';
Ref(RefStyle).Ket = ']';
Ref(RefStyle).V = '[%d]';
Ref(RefStyle).MultiSeparator = '–';
Ref(RefStyle).MultipleSeparator = ',';
Ref(RefMultiStyle).Style = strcat('[%d]',char(8211),'[%d]'); %'[%d]-[%d]'
Ref(RefMultipleStyle).StyleRightUnit = '%d]'; %'[%d],[%d]' %',[%d]'
Ref(RefMultipleStyle).StyleLeftUnit = '[%d,';
CheckOverRefIDPos_Condition = 1;
RefSelectionMode = 'wh.Selection.MoveRight(wdWord,3,wdExtend)'; % wh.Selection.MoveRight(wdWord,1,wdExtend);




wdStory = 6;
wdFindStop = 0;
wdCharacter = 1;
wdExtend = 1;
wdWord = 2;
wdActiveEndPageNumber = 3;
wdHorizontalPositionRelativeToPage = 5;
wdVerticalPositionRelativeToPage = 6;
wdSortByLocation = 1;

[FileName,PathName,FileIndex] = uigetfile({'*.doc;*.docx'},'Load doc file');
if FileIndex == 0, return, end
Path = [PathName FileName];
%Path = 'C:\Users\Edo\OneDrive - Politecnico di Milano\Articoli\Articoli Miei\In progress\Articolo Mammografo\2.0\Word\test.docm';
wh = actxserver('Word.Application');
wh.Visible = true;
trace(wh.Visible);
WH = invoke(wh.Documents,'Open',Path);

if Ref(RefParagPOS).V == Eof
    wh.Selection.EndKey(wdStory);
    is_forward = false;
else
    wh.Selection.HomeKey(wdStory);
    is_forward = true;
end

is_found=MWFind(wh,Ref(RefParagID).V,'',[is_forward,wdFindStop,false,false,false,false,false,false]);
if is_found == false, return, end

RefParagIDPos = wh.Selection.Start;
wh.Selection.End = RefParagIDPos;
RefParagIDPage = wh.Selection.Information(wdActiveEndPageNumber);
RefParagIDVp = wh.Selection.Information(wdVerticalPositionRelativeToPage);
RefParagIDHp = wh.Selection.Information(wdHorizontalPositionRelativeToPage);

iRef = 1;
while iRef<=NumRef
    switch Ref(RefParagStyle).V
        case '[%d]'
            Ref2Find = strcat('[',num2str(iRef),']');
        case '%d.'
            Ref2Find = strcat(num2str(iRef),'.');
    end
    
    is_found=MWFind(wh,Ref2Find,'',[FForwardT,wdFindStop,false,false,false,false,false,false]);
    if is_found
        wh.Selection.Start = wh.Selection.End;
        %perform check
        Valid = true;
        if Valid
            NameBookMark = strcat('Ref',num2str(iRef));
            wh.ActiveDocument.Bookmarks.Add(NameBookMark,wh.Selection.Range);
            wh.ActiveDocument.Bookmarks.DefaultSorting = wdSortByLocation;
            iRef = iRef +1;
        end
    end
end

if Ref(RefParagPOS).V == Eof
    wh.Selection.Start = RefParagIDPos;
    wh.Selection.End = wh.Selection.Start;
    wh.Selection.HomeKey(wdStory);
    is_forward = true;
else
    %wh.Selection.HomeKey(wdStory);
    is_forward = true;
end

OUT_WHILE = false;
while ~OUT_WHILE
    switch Ref(RefStyle).V
        case '[%d]'
            Ref2Find = Ref(RefStyle).Bra;
    end
    
    is_found=MWFind(wh,Ref2Find,'',[true,wdFindStop,false,false,false,false,false,false]);
    switch Ref(RefParagPOS).V
        case Eof
            if(wh.Selection.Information(wdActiveEndPageNumber) == RefParagIDPage)
                if(wh.Selection.Information(wdVerticalPositionRelativeToPage) >= RefParagIDVp)
                    if(wh.Selection.PageSetup.TextColumns.Count >=2)
                        if(wh.Selection.Information(wdHorizontalPositionRelativeToPage) >= RefParagIDHp)
                            OUT_WHILE  = true;
                            if OUT_WHILE, is_found = false; end
                        end
                    else
                        OUT_WHILE  = true;
                        if OUT_WHILE, is_found = false; end
                    end
                end
            end
        case Bof
            OUT_WHILE = is_found == false;
    end
    
    if is_found
        % Perform Check
        
        
        is_complete_ref = false;
        MoveDir = 'wh.Selection.MoveRight(wdCharacter,1,wdExtend);';
        while ~is_complete_ref
            eval(MoveDir);
            Buffer=wh.Selection.Text;
            is_complete_ref = strfind(Buffer,Ref(RefStyle).Ket);
            if ~isempty(is_complete_ref)
                if is_complete_ref(end)
                    %is_complete_ref = is_complete_ref == length(Buffer);
                    if is_complete_ref(end) < length(Buffer)
                        MoveDir = 'wh.Selection.MoveLeft(wdCharacter,1,wdExtend);';
                        is_complete_ref = 0;
                    else
                        if(is_complete_ref(end) == length(Buffer))
                            is_complete_ref = true;
                        end
                    end
                end
            else
                is_complete_ref = false;
            end
        end
        
        RefString = wh.Selection.Text;
        RefString = RefString(~isspace(RefString));
        ReadRefNum = str2num(RefString(~isspace(RefString)));
        MultiRef = strfind(RefString,Ref(RefStyle).MultiSeparator);
        MultipleRef = strfind(RefString,Ref(RefStyle).MultipleSeparator);
        clear('BufferBM','BufferToolTip')
        is_single_ref = false;
        is_multiref = false;
        is_multipleref = false;
        if ~isempty(ReadRefNum)
            BMRef = strcat('Ref',num2str(ReadRefNum));
            ToolTip = strcat('Ref:',num2str(ReadRefNum));
            is_single_ref = true;
        end
        if ~isempty(MultiRef)
            ReadRefNum = sscanf(RefString,Ref(RefMultiStyle).Style);
            BMRef = strcat('Ref',num2str(ReadRefNum(1)));
            ToolTip = strcat('Ref:',num2str(ReadRefNum(1)),'-',num2str(ReadRefNum(2)));
            is_multiref = true;
        end
        if ~isempty(MultipleRef)
%             Format = Ref(RefMultipleStyle).StyleLeftUnit;
%             for iMR = 1:length(MultipleRef)
%                 Format = strcat(Format,Ref(RefMultipleStyle).StyleRightUnit);
%             end
            ReadRefNum = str2num(RefString);
            is_multipleref = true;
            for iMR = 1 : length(ReadRefNum)
                BufferBM{iMR} = strcat('Ref',num2str(ReadRefNum(iMR))); %#ok<SAGROW>
                BufferToolTip{iMR} = strcat('Ref:',num2str(ReadRefNum(iMR))); %#ok<SAGROW>
            end
%             ToolTip = strcat('Ref:',num2str(ReadRefNum,['%d',Ref(RefStyle).MultipleSeparator]));
%             ToolTip(end) = [];
        end
        if ~isempty(MultipleRef) && ~isempty(MultipleRef)
            
        end
        if is_multipleref
            StringSelection = wh.Selection.Text;
            LenghtSelection = length(wh.Selection.Text);
            wh.Selection.Delete
            wh.Selection.TypeText(StringSelection)
            wh.Selection.MoveLeft(wdCharacter,LenghtSelection);
            TempPos = wh.Selection.End;
             wh.Selection.Start = wh.Selection.End;
             %while (wh.Selection.Characters.First(1).Text) ~= Ref(RefStyle).Bra
             for iMR = 1:length(ReadRefNum)
                 iMove = 0;
                 if iMR == length(ReadRefNum)
                     EndSearch = Ref(RefStyle).Ket;
                     %EndSearch = Ref(RefStyle).Bra;
                 else
                     EndSearch = Ref(RefStyle).MultipleSeparator;
                 end
                 while (wh.Selection.Characters.Last(1).Text) ~= EndSearch
                    wh.Selection.MoveRight(wdCharacter,1,wdExtend);
                    iMove = iMove +1;
                 end
                 
%                  iMove = 0;
%                  while (wh.Selection.Characters.Last(1).Text) ~= EndSearch
%                     wh.Selection.MoveRight(wdCharacter,1,wdExtend);
%                     iMove = iMove +1;
%                  end
                 BMRef = BufferBM{iMR};
                 ToolTip =  BufferToolTip{iMR};
                 wh.ActiveDocument.Hyperlinks.Add(wh.Selection.Range,'',BMRef,ToolTip,wh.Selection.Text);
                 %wh.Selection.MoveRight(wdCharacter,1);
             end
             %wh.Selection.MoveRight(wdCharacter,1);
%              while (wh.Selection.Characters.Last(1).Text) == Ref(RefStyle).Ket
%                 wh.Selection.MoveRight(wdCharacter,1);
%                 iMove = iMove +1;
%              end
             uiwait(msgbox(strcat({'Check:',StringSelection}),'Exception','modal'));
        else
            wh.ActiveDocument.Hyperlinks.Add(wh.Selection.Range,'',BMRef,ToolTip,wh.Selection.Text);
        end
        wh.Selection.Start=wh.Selection.End;
        HLSet = true;
    end
end






%%%%% Functions
function [is_found] = MWFind(wh,text2find,text2replace,FOptions)
wh.Selection.Find.ClearFormatting;
wh.Selection.Find.Replacement.ClearFormatting;
wh.Selection.Find.Text = text2find;
wh.Selection.Find.Replacement.Text = text2replace;
wh.Selection.Find.Forward = FOptions(1);
wh.Selection.Find.Wrap = FOptions(2);
wh.Selection.Find.Format = FOptions(3);
wh.Selection.Find.MatchCase = FOptions(4);
wh.Selection.Find.MatchWholeWord = FOptions(5);
wh.Selection.Find.MatchWildcards = FOptions(6);
wh.Selection.Find.MatchSoundsLike = FOptions(7);
wh.Selection.Find.MatchAllWordForms = FOptions(8);
wh.Selection.Find.Execute;
is_found = wh.Selection.Find.Found;
end