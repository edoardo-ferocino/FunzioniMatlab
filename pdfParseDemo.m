% Purpose : Demonstrate extracting text from a PDF file using PDFBox Java library
% Usage   : Modify file paths
%           Enable cell mode and step through the code  
% Example : none (Oh, the FEX code metrics..)
% Author  : Dimitri Shvorob, dimitri.shvorob@gmail.com, 5/1/08  

%% 
clear java
javaaddpath('C:\Users\Edo\OneDrive - Politecnico di Milano\LabWorks\Funzioni Matlab\FontBox-0.1.0\lib\FontBox-0.1.0.jar')

javaaddpath(...
'C:\Users\Edo\OneDrive - Politecnico di Milano\LabWorks\Funzioni Matlab\PDFBox-0.7.3\lib\PDFBox-0.7.3.jar')

%%
pdfdoc = org.pdfbox.pdmodel.PDDocument;
reader = org.pdfbox.util.PDFTextStripper;
pdfdoc.close;

%%
pdfdoc = pdfdoc.load('C:\Users\Edo\Documents\Personale\Estratti conto e cedolini\Estratti conto\Estratto_Conto_17_9.pdf');
pdfdoc.isEncrypted;

%% text, with planty of padding
pdfstr = reader.getText(pdfdoc);                  %#ok

%%
class(pdfstr);

%%
pdfstr = char(pdfstr);                            %#ok

%%
class(pdfstr);

%% text 'unpadded'
pdfstr = deblank(pdfstr);                         %#ok

% write pdfstr in the file text.txt 
fid = fopen('text.txt','w'); 
fprintf(fid,'%s',pdfstr);
fclose(fid);
fid = fopen('text.txt','r'); 
il = 1;
while feof(fid)==false
    readtext = {fgetl(fid)};
    ReadText(il,1) = readtext;
    switch il
    case 1
        formtext.IBAN = readtext;
    case 5
        formtext.Data = readtext;
    end
    il = il+1;
end
fclose(fid);
[starts,ends]=regexp(ReadText,'\d[/]');
indexs = find(~cellfun(@custom,starts));
Data = find(cellfun(@custom2,ReadText(indexs)));



function vect = custom(vect)
    if isempty(vect)
        vect = 1;
    else
        if numel(vect)<4
            vect = 1;
        else
            vect = 0;
        end
    end
end

function vect = custom2(vect)
buffer = regexp(vect,'[ ]','split');
Data1 = buffer(1);
Data2 = buffer(2);
% value = buffer(Testo = vect(18:end);
if isempty(vect)
        vect = 1;
    else
        if numel(vect)<4
            vect = 1;
        else
            vect = 0;
        end
    end
end