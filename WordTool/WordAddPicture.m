function WordAddPicture(actx_word,file_name,link_to_file,save_with_doc)
    actx_word.Selection.InlineShapes.AddPicture(file_name,link_to_file,save_with_doc);
    WordInsertParagraph(actx_word,1);
end