%%--------------------------------
function Sz=DMC_zaklocenia_odp_skok_Sz(Tp,Gz,timefinal)
% konstrukcja macierzowej odpowiedzi skokowej zakłóceniowej dla DMC,
% zakładając timefinal= Tp*Dz, Dz -- liczba elementów odp. skokowej

% przykład obiektu: transmitancje ciągłe zakłóceniowe:
Gzdtf=c2d(Gz,Tp);   
% Dyskretna odpowiedź skokowa wielowymiarowa wg konwencji Matlaba:
Ydstepz=step(Gzdtf,timefinal);% macierz wymiaru (timefinal/Tp)+1 x ny x nz na odcinku
% czasu od t=0 do t=timefinal z krokiem Tp
% Dyskretna odpowiedź skokowa macierzowa Sz o wymiarze ny x nz x timefinal/Tp; pierwsza
% macierz dla czasu dyskr. k=1, ostatnia dla k=Dz=timefinal/Tp (Dz macierzy):
    [nt,ny,nz]=size(Ydstepz); % nt=D+1 (Ydstepz zawiera wartość wyjść w chwili 0, Sz nie)
    Sz=zeros(ny,nz,nt-1);
    for i=1:ny
        for j=1:nz
            Sz(i,j,:)=Ydstepz(2:nt,i,j); 
        end
    end
end