function [Found] = WordFind(actx_word,text,options)
actx_word.Selection.Find.ClearFormatting;
actx_word.Selection.Find.Replacement.ClearFormatting;
actx_word.Selection.Find.Text = text;

wdFindContinue = 1;
actx_word.Selection.Find.Format = wdFindContinue;
actx_word.Selection.Find.MatchCase = false;
actx_word.Selection.Find.MatchWholeWord = false;
actx_word.Selection.Find.MatchWildcards = false;
actx_word.Selection.Find.MatchSoundsLike = false;
actx_word.Selection.Find.MatchAllWordForms = false;
actx_word.Selection.Find.Forward = false;
actx_word.Selection.Find.Execute(text,options(1),options(2),options(3),options(4),options(5),options(6),options(7));
Found=actx_word.Selection.Find.Found;
