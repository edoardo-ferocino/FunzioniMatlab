function [ Data, varargout ] = DatRead(FileName,NumRep,NumChan,NumBin,varargin)
HeadLen=764;
SubLen=204;
FilePath=[FileName '.DAT'];
fid=fopen(FilePath,'r');
if fid<0, error('File not found'); return; end
Head=fread(fid,HeadLen,'char');
A=zeros(NumRep,NumChan,NumBin);
Sub=zeros(NumRep,NumChan,SubLen);
info=dir(FilePath);
% if(info.bytes*8==((NumRep*NumChan)*(16*NumBin+8*SubLen)+HeadLen*8))
%     datatype = 'ushort';
% else
%     datatype = 'uint32';
% end
if nargin==5 && strcmp(varargin{1},'uint32')
    datatype = 'uint32';
else
    datatype = 'ushort';
end
for ir=1:NumRep
   for ich=1:NumChan
        Sub(ir,ich,:)=fread(fid,SubLen,'char');
        A(ir,ich,:)=fread(fid,NumBin,datatype);
   end
end
fclose(fid);

switch nargin
    case 4
        Data=(sum(A,1));
    case 5
        if strcmp(varargin{1},'Squeeze')
            Data=squeeze(sum(A,1));
        end
        if strcmp(varargin{1},'All')
                Data=A;
        end
        if strcmp(varargin{1},'single')
                Data=sum(A,1);
        end
end

switch nargout
    case 2
        varargout = {Head};
    case 3
        varargout(1) = {Head};
        varargout(2) = {Sub};

end

