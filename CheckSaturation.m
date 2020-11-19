function SaturatedCurves = CheckSaturation(Data,DNLcoeff,OPTODE,FileName)
Ndims = ndims(Data);
if Ndims == 4
    [NumRep,NumSeq,~,~]=size(Data);
else
    [NumSeq,~,~]=size(Data);
    NumRep = 1;
    Data = permute(Data,[4 1 2 3]);
end
Data = Data(:,:,OPTODE,:);
NumOptode = 1;
RoiV = 1:120;
% IntensityCh = [128 127];
% Intensity = Data(:,:,:,IntensityCh);
% Intensity = Intensity(:,:,:,1)*power(2,16)+Intensity(:,:,:,2);
% CW = sum(Data(:,:,:,RoiV),4);

Data = Data(:,:,:,RoiV);
Data = reshape(Data,NumRep*NumSeq*NumOptode,size(Data,4));
DNL_matrix = repmat(permute(DNLcoeff(RoiV),[2 1]),[NumRep*NumSeq*NumOptode,1]);
Data = Data.*DNL_matrix;
Logical = find(max(Data,[],2) > 0.8 * 4096);
Deriv = diff(Data,1,2);
[~,MaxIndex]=max(Deriv,[],2);
zcd = dsp.ZeroCrossingDetector;
Zc = zeros(numel(Logical),2);
iz = 0;
for il = Logical'
    iz = iz +1;
    Zc(iz,:)= [zcd(Deriv(il,MaxIndex(il)-3:MaxIndex(il)+5)'),il];
    zcd.reset
end

% DiffCounts = (Intensity-CW)./Intensity * 100;

warning('off','verbose')
warning('off','backtrace')
% if any(DiffCounts(:)>10)
%     warning(strcat('Could be saturation in file: ',FileName))
% end
if any(Zc(:,1)>1)
    warning(strcat('Most probably saturation in file: ',FileName))
end
SaturatedCurves = Data(Zc(Zc(:,1)>1,2),:);
warning('off','verbose')
warning('off','backtrace')
end