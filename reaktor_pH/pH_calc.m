clear all

Ka1=4.47*10^-7;
Ka2=5.62*10^-11;
pK1=-log10(Ka1);
pK2=-log10(Ka2);

a4=10^-(pK2+14);
a0=-10^pK1;

Wa4=-4.315622133383698e-004;
Wb4=5.277571688626062e-004;

a3=10^-14+Wa4*10^-pK2+2*Wb4*10^-pK2;
a2=10^(pK1-14)+Wa4+Wb4-10^-pK2;
a1=Wa4*10^pK1-1;

p_kand=roots([a4 a3 a2 a1 a0]);
p=0;
for j=1:4
    if isreal(p_kand(j))&&(p<p_kand(j))
        p=p_kand(j);
    end
end
pH=log10(p)
