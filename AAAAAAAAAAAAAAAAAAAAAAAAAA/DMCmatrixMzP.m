function MzP=DMCmatrixMzP(Sz,N,Nu)
% script "DMCmatrixMzP" macierze MzP regulatora DMC,
% DANE WEJŚCIOWE (wymagane):
% Sz(ny,nz,Dz) - macierz odpowiedzi na skoki zakłóceń, gdzie ny = dim y,
% nz = dim z; Dz - horyzont dynamiki zakłóceń od k=1 do k=Dz (uwzględniany w Sz)
% N - horyzont predykcji,
% Nu - horyzont sterowania.
%
    [ny,nz,Dz]=size(Sz);
% Mz1 - pierwsza kolumna macierz dynamicznej Mz:
    Mz1=zeros(N*ny,nz);
    for i=1:N
        Mz1((i-1)*ny+1:(i-1)*ny+ny,1:nz)=Sz(:,:,min(i,Dz)); 
    end
% macierz MP:
    MzP0=zeros(ny*N,nz*(Dz-1));
    for i=1:Dz-1
        for j=1:N
            MzP0((j-1)*ny+1:j*ny,(i-1)*nz+1:i*nz)=Sz(:,:,min(i+j,Dz))-Sz(:,:,i);
        end
    end
    MzP=[Mz1 MzP0];
end