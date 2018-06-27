function WordText(actx_word_p,text_p,style_p,enters_p,color_p)
	%VB Macro
	%Selection.TypeText Text:="Test!"
	%in Matlab
	%set(word.Selection,'Text','test');
	%this also works
	%word.Selection.TypeText('This is a test'); 
    if ~isempty(text_p)
    if(enters_p(1))
        WordInsertParagraph(actx_word_p,enters_p(1));
    end
	actx_word_p.Selection.Style = style_p;
    if(nargin == 5)%check to see if color_p is defined
        actx_word_p.Selection.Font.ColorIndex=color_p;     
    end
    if strcmpi(style_p,'Normal')
        actx_word_p.Selection.Font.Name = 'Times New Roman';
        actx_word_p.Selection.ParagraphFormat.SpaceAfter = 0;
        actx_word_p.Selection.ParagraphFormat.LineSpacingRule = 0;
    end
    [nr,nc]=size(text_p);
    if nr>1
    for inr=1:nr
        actx_word_p.Selection.TypeText(text_p(inr,:));
        actx_word_p.Selection.TypeParagraph; %enter
    end
    else
        actx_word_p.Selection.TypeText(text_p);
    end
    actx_word_p.Selection.Font.ColorIndex='wdAuto';%set back to default color
    WordInsertParagraph(actx_word_p,enters_p(2));
    end
return