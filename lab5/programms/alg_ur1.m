pl=[2   0  1];
p2=[1   0   0   -1];
%��������� ���������:   (2�^2+1) (�^3-1) = 2�^5+�^3-2�^2-1
p=conv(p2,pl)
%������� ���������:   (2�^5+�^3-2�^2-1) / (�^3-1) = (2�^2+1) 
deconv(p,p2) 