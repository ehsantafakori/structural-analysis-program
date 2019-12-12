 

scale=1;

scaled=1;

 



XOfNode=[4,2,2,0,0];

YOfNode=[0,2,0,0,2];

 mass=[1,1,1,1,1,1];

 



COMF=[4,3,1,3,2,2];

COMS=[3,1,2,2,4,5];

DgrOfFreDm=[1,1,1,1,1,1,0,0,0,0];

A=0.0015*[1,1,1,1,1,1];

E=2000*[1,1,1,1,1,1];

force=[0,-30,0,0,0,-30,0,0,0,0];

SbsdbcOfsprt=[0,0,0,0];

ErrInBld=[0,0,0,0,0,0];

DltaT=[0,0,42.55,0,0,0];

Alfa=0.000011*[1,1,1,1,1,1];





%     Matrix analysis





NmbrOfEle=length(E);

Asmblelemat=zeros(NmbrOfEle,4);

teta=zeros(NmbrOfEle);

tele=zeros(4,4);

l=zeros(NmbrOfEle);

for i=1:1:NmbrOfEle;

l(i)=((XOfNode(COMS(i))-XOfNode(COMF(i)))^2+(YOfNode(COMS(i))-YOfNode(COMF(i)))^2)^0.5;

teta(i)=atan((YOfNode(COMS(i))-YOfNode(COMF(i)))/(XOfNode(COMS(i))-XOfNode(COMF(i))));

Tele(:,:,i)=[cos(teta(i)),sin(teta(i)),0,0;-1*sin(teta(i)),cos(teta(i)),0,0;0,0,cos(teta(i)),sin(teta(i));0,0,-1*sin(teta(i)),cos(teta(i))];

tele(:,:)=Tele(:,:,i);

kelelcl(:,:,i)=[A(i)*E(i)/l(i),0,-A(i)*E(i)/l(i),0;0,0,0,0;-A(i)*E(i)/l(i),0,A(i)*E(i)/l(i),0;0,0,0,0];

kele(:,:,i)=transpose(tele)*kelelcl(:,:,i)*tele;

Asmblelemat(i,:)=[2*COMF(i)-1,2*COMF(i),2*COMS(i)-1,2*COMS(i)];

end

ktalsz= length(DgrOfFreDm);

ktotal=zeros(ktalsz,ktalsz);

fdltta=zeros(4,1);

ferr=zeros(4,1);

for i=1:1:NmbrOfEle;

tele=Tele(:,:,i);

fdltta=transpose(tele)*[-Alfa(i)*A(i)*E(i)*DltaT(i);0;Alfa(i)*A(i)*E(i)*DltaT(i);0];

ferr=transpose(tele)*[-ErrInBld(i)*A(i)*E(i)/l(i);0;ErrInBld(i)*A(i)*E(i)/l(i);0];

for j=1:1:4

 force(Asmblelemat(i,j))=force(Asmblelemat(i,j))+fdltta(j);

force(Asmblelemat(i,j))=force(Asmblelemat(i,j))+ferr(j);

for k=1:1:4

ktotal(Asmblelemat(i,j),Asmblelemat(i,k))=ktotal(Asmblelemat(i,j),Asmblelemat(i,k))+kele(j,k,i);

end

end

end

keliminateff=ktotal;

keliminatefs=ktotal;

keliminatess=ktotal;

forceeleminate=force;

for i=ktalsz:-1:1

if DgrOfFreDm(i)==0

keliminateff(i,:)=[];

keliminateff(:,i)=[];

keliminatefs(i,:)=[];

forceeleminate(i)=[];

else

keliminatess(i,:)=[];

keliminatess(:,i)=[];

keliminatefs(:,i)=[];

end

end

degOffr=length(forceeleminate);

keliminateForce=keliminateff;

keliminatefree=keliminateff;

keliminateforcefree=keliminateff;

forceeleminatefective=forceeleminate;

for i=degOffr:-1:1

if forceeleminate(i)==0



forceeleminatefective(i)=[];

keliminateForce(i,:)=[];

keliminateForce(:,i)=[];

keliminateforcefree(i,:)=[];

else



keliminatefree(i,:)=[];

keliminatefree(:,i)=[];

keliminateforcefree(:,i)=[];

end

end

displforce=zeros(length(forceeleminatefective),1);

displfree=zeros(length(forceeleminate)-length(forceeleminatefective),1);

flaggg=det(keliminatefree);

if flaggg<1*10^-10

disp truss is unstable

else 

flaggg=det(keliminateForce-keliminateforcefree*keliminatefree^-1*transpose(keliminateforcefree));

end

if flaggg<1*10^-10

disp truss is unstable

else

displforce=(keliminateForce-keliminateforcefree*keliminatefree^-1*transpose(keliminateforcefree))^-1*transpose(forceeleminatefective);

displfree=-1*keliminatefree^-1*transpose(keliminateforcefree)*displforce;

end



assmbldisplforcedisplfree=zeros(length(forceeleminate),1);

zz=0;

kk=0;

for i=1:1:length(forceeleminate);

if forceeleminate(i)==0

zz=zz+1;

assmbldisplforcedisplfree(i)=displfree(zz);

else

kk=kk+1;

assmbldisplforcedisplfree(i)=displforce(kk);

end

end



displ=zeros(degOffr,1);

rActn=zeros(ktalsz-degOffr,1);

flagg=det( keliminateff);

if flagg<1*10^-10

disp truss is unstable

else

% displ= keliminateff^-1*transpose(forceeleminate)-keliminateff^-1*keliminatefs*transpose(SbsdbcOfsprt);

displ= assmbldisplforcedisplfree-keliminateff^-1*keliminatefs*transpose(SbsdbcOfsprt);



rActn=transpose(keliminatefs)*displ+keliminatess*transpose(SbsdbcOfsprt);

 disp (displ)

disp (rActn)

end

dispfinal=zeros(ktalsz,1);

zz=0;

kk=0;

for i=1:1:ktalsz;

if DgrOfFreDm(i)==0

zz=zz+1;

dispfinal(i)=SbsdbcOfsprt(zz);

else

kk=kk+1;

dispfinal(i)= displ(kk);

end

end

felefinalle=zeros(NmbrOfEle,1);

 for i=1:1:NmbrOfEle;

tele=Tele(:,:,i);

felefinalle(i)=[-A(i)*E(i)/l(i),0,A(i)*E(i)/l(i),0]*tele*[dispfinal(2*COMF(i)-1);dispfinal(2*COMF(i));dispfinal(2*COMS(i)-1);dispfinal(2*COMS(i))]-Alfa(i)*A(i)*E(i)*DltaT(i)-ErrInBld(i)*A(i)*E(i)/l(i);

end

massmatrix=zeros(degOffr,degOffr);

for i=1:1:degOffr

 massmatrix(i,i)=mass(i);

end

[mode,perid]=eig(massmatrix^-1*keliminateff);





%    Graphical output





XGrph=zeros(2*NmbrOfEle,1);

YGrph=zeros(2*NmbrOfEle,1);

for i=1:1:NmbrOfEle

XGrph(2*i-1,1)=XOfNode(COMF(i));

XGrph(2*i,1)=XOfNode(COMS(i));

YGrph(2*i-1,1)=YOfNode(COMF(i));

YGrph(2*i,1)=YOfNode(COMS(i));

end

figure (1)

 plot (XGrph,YGrph)

 axis([min(XGrph)-0.1*max(XGrph),max(XGrph)+0.1*max(XGrph),min(YGrph)-0.1*max(XGrph),max(YGrph)+0.1*max(XGrph)])

xlabel('X')

ylabel('Y')

 XGrphAftr=zeros(2*NmbrOfEle,1);

YGrphAftr=zeros(2*NmbrOfEle,1);

XOfNodeAftr=zeros(length(XOfNode),1);

YOfNodeAftr=zeros(length(XOfNode),1);

j=0;

zz=0;

 for i=1:1:length(XOfNode)

if   DgrOfFreDm(2*i-1)==1

j=j+1;

XOfNodeAftr(i)=XOfNode(i)+displ(j)*scale;

else

zz=zz+1;

XOfNodeAftr(i)=XOfNode(i)+SbsdbcOfsprt(zz)*scaled;

end

if   DgrOfFreDm(2*i)==1

j=j+1;

YOfNodeAftr(i)=YOfNode(i)+displ(j)*scale;

else

zz=zz+1;

YOfNodeAftr(i)=YOfNode(i)+SbsdbcOfsprt(zz)*scaled;

 end

end

 for i=1:1:NmbrOfEle

XGrphAftr(2*i-1,1)=XOfNodeAftr(COMF(i));

XGrphAftr(2*i,1)=XOfNodeAftr(COMS(i));

YGrphAftr(2*i-1,1)=YOfNodeAftr(COMF(i));

YGrphAftr(2*i,1)=YOfNodeAftr(COMS(i));

end

figure (2)

 plot (XGrphAftr,YGrphAftr)

axis([min(XGrphAftr)-0.1*max(XGrph),max(XGrphAftr)+0.1*max(XGrph),min(YGrphAftr)-0.1*max(XGrph),max(YGrphAftr)+0.1*max(XGrph)])

 xlabel('X')

ylabel('Y')

figure (3)

plot (XGrphAftr,YGrphAftr,'r',XGrph,YGrph,'b')

 axis([min(min(XGrphAftr),min(XGrph))-0.1*max(XGrph),max(max(XGrphAftr),max(XGrph))+0.1*max(XGrph),min(min(YGrphAftr),min(YGrph))-0.1*max(XGrph),max(max(YGrphAftr),max(YGrph))+0.1*max(XGrph)])

xlabel('X')

ylabel('Y')

  



