function [actx_word,word_handle]=StartWord(word_file_p,Visible)
% Start an ActiveX session with Word:
actx_word = actxserver('Word.Application');
actx_word.Visible = Visible;
trace(actx_word.Visible);
if ~exist(word_file_p,'file')
    % Create new document:
    word_handle = invoke(actx_word.Documents,'Add');
    
else
    % Open existing document:
    word_handle = invoke(actx_word.Documents,'Open',word_file_p);
    actx_word.Selection.WholeStory;
    wdCharacter = 1;
    actx_word.Selection.Delete(wdCharacter,1);
end
return