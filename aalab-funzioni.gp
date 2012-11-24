{\\lettera -> codice Ascii,  trasforma anche una frase: le lettere devono essere tra "".
lett2num(a)=Vec(Vecsmall(a));
}



{\\numero -> lettera; trasforma anche un vettore di numeri
num2lett(n)=Strchr(n);
}



{\\inverte l'ordine in un vettore
contr(v)=Vec(Polrev(v));
}



{\\restituisce la lista delle posizioni di x nella lista lista
posizioni(lista,x)=local(nn,ris);
nn=length(lista);ris=[];
for(kk=1,nn,if(lista[kk]==x,ris=concat(ris,kk)));
ris;
}



{\\scambia nel vettore v gli elementi di posto i e j
scambia(v,i,j)=local(a,provv);
if(i<=0||j<=0,error("gli indici devono essere positivi!"));
if(i>length(v)||j>length(v),error("gli indici non devono essere maggiori della lunghezza del vettore!"));
a=v;
provv=a[i];
a[i]=a[j];
a[j]=provv;
a;
}


{\\trova il periodo di elem (con operazione oper e neutro operid, vedi pot) se ordine è un suo multiplo
cercaperiodo(elem,ordine)=local();
fordiv(ordine,n,if(pot(elem,n)==operid,return(n)));
0;
}


{\\converte n alla base b, restituisce il vettore delle cifre
converti(n,b)=local(nn,ris,gg);
if(n<0,return([]));
if(n==0,return([0]));
if(n==1,return(1));
ris=[];
nn=n;
while(nn>0, gg=divrem(nn,b); ris=concat(ris,gg[2]); nn=gg[1]);
contr(ris);
}


{\\converte n alla base b, restituisce un vettore di lunghezza L che contiene le cifra ed eventuali zeri iniziali
converti2(n,b,L)=local(nn,ris,gg);
if(n<0,return([]));
if(n==0,return(vector(L)));
if(n==1,return(concat(vector(L-1),[1])));
ris=[];
nn=n;
while(nn>0, gg=divrem(nn,b); ris=concat(ris,gg[2]); nn=gg[1]);
if(length(ris)>L, error("L troppo piccolo!"));
if(length(ris)<L,vv=vector(L-length(ris));ris=concat(ris,vv));
contr(ris);
}


{\\rimescola un vettore v mediante m scambi
rimescola(v,m)=local(n,w);
n=length(v);
w=v;
for(i=1,m,g=random(n);h=random(n);w=scambia(w,g+1,h+1));
w;
}



{\\permutazione casuale di lunghezza n, ottenuta com rimescola 1..n, m volte
randomperm(n,m)=rimescola(vector(n,k,k),m);
}




{\\applica la permutazione p al vettore v; serve anche a comporre due permutazioni p e v
\\se la lunghezza di v è n, p deve essere una permtazione dei numeri da 1 a n
applica(p,v)=local(n,pc,w);
n=length(v);
if(length(p)!=length(v),error("le lunghezze devono essere uguali!"));
pc=vecsort(p);
if(pc-vector(n,k,k)!=vector(n),error("v non è una permutazione corretta"));
w=v;
for(i=1,n,w[i]=v[p[i]]);
w;
}



{\\permutazione inversa
inversa(p)=local(n,pc,w);
n=length(p);
pc=vecsort(p);
if(pc-vector(n,k,k)!=vector(n),error("v non è una permutazione corretta"));
w=p;
for(i=1,n,g=posizioni(p,i);w[i]=g[1]);
w;
}



{\\lista 1...n
identica(n)=vector(n,k,k);
}



{\\schema dell'algoritmo potenza veloce; operid è l'elemento neutro, oper è l'operazione
pot(a,m)=local(u,t,v);
u = m; t = operid; v = a;
while(u != 0, if(u%2 == 1, t = oper(t,v)); v = oper(v,v);
u = floor(u/2));
t;
}




{\\dati due vettori a,b restituisce il vettore (senza ripetizioni) formato dagli elementi di a che non sono in b
complemento(a,b)=local();
eval(setminus(Set(a),Set(b)));
}



{\\composto p1 o p2, si applica prima p2 e poi p1
permcomp(p1,p2)=local(ris);
nn=length(p1);
ris=vector(nn);
for(i=1,nn,ris[i]=p1[p2[i]]);
ris;
}



{\\restituisce la decomposizione di una permutazione nei suoi cicli disgiunti
permdecomp(perm)=local(nn,ris,ss,oo,kk);
nn=length(perm);
ris=[];
ss=vector(nn,k,k);
while(length(ss)!=0,ele=ss[1];oo=[ele];kk=ele;while(perm[kk]!=ele,kk=perm[kk];oo=concat(oo,kk));ris=concat(ris,[oo]);
ss=complemento(ss,oo));
ris;
}



{\\trova il periodo di una permutazione
permperiod(perm)=local(dec);
dec=permdecomp(perm);
lcm(vector(length(dec),k,length(dec[k])));
}



{\\restituisce la parità di una permutazione
permparity(perm)=local(uu,cont);
uu=permdecomp(perm);
cont=0;
nn=length(uu);
for(k=1,nn,if(length(uu[k])%2==0,cont++));
if(cont%2==1,return(1),return(0));
}



{\\crea la matrice permutazionale associata alla permutazione perm
permat(perm)=local(nn,mm);
nn=length(perm);
mm=matrix(nn,nn);
for(i=1,nn,mm[perm[i],i]=1);
mm;
}



{\\inversa di permat
matperm(a)=local(ris,n);
ris=vector(matsize(a)[1]);
n=length(ris);
for(i=1,n,ris[i]=posizioni(a[,i],1)[1]);
ris;
}



{\\media aritmetica delle componenti di un vettore
media(v)=local();
nn=length(v);
cc=sum(i=1,nn,v[i]);
cc/nn;
}



{\\dice quanti elementi y ci sono nel vettore v
quanti(v,y)=length(posizioni(v,y));
}



{\\funzione di Moebius
moeb(n)=local(ee);
if(n==1,return(1));
ee=factor(n)[,2];
if(vecmax(ee)>1,return(0));
return((-1)^length(ee));
}


{\\numirr, calcola il numero dei polinomi monici irriducibili di grado m con coefficienti in GF(q)
numirr(q,m)=local();
ris=0;
fordiv(m,d,ris=ris+moeb(m/d)*q^d);
1/m*ris;
}


{\\prodirr, calcola il prodotto dei polinomi monici irriducibili di grado m con coefficienti in GF(q)
prodirr(q,m)=local(ris);
ris=1;
fordiv(m,d,ris=ris*(x^(q^d)-x)^(moeb(m/d)));
simplify(ris);
}


{\\restituisce l'insieme dei coprimi con n
coprimi(n)=local();
ris=[];
for(i=1,n-1,if(gcd(i,n)==1,ris=concat(ris,i)));
ris;
}



{\\Restituisce in polinomio primitivo di grado n in Z_p[var]
\\L'algoritmo elimina via via gli elementi il cui ordine è un divisore proprio di p^n-1
primpoly(p, n, var) = local(s,f);
f = var^(p^n-1) -1; s = divisors(p^n -1); 
for(k=1,length(s)-1, f=f/gcd(f,var^s[k] - 1));
factormod(f, p)[1, 1];
}



{\\costruisce il campo GF(p^n)
\\utilizza un polinomio primitivo opportuno, ff=primpoly(p,n,var), nella variabile desiderata var,  e costruisce una \\radice generica al 
\\restituisce l'insieme delle P^n-1 potenze distinte di al
\\unito con 0, al primo posto. All'ultimo posto c'è 1
campop(p,n,var)=local(ff,al,kk);
ff=primpoly(p,n,var);
al=Mod(Mod(1,p)*var,ff);
kk=vector(p^n-1,h,al^h);
kk=concat([0],kk);
}


{\\sovrappone due matrici quadrate dello stesso ordine
sovrapponi(a,b)=local(s,ris);
s=length(a);
ris=matrix(s,s);
forvec(x=[[1,s],[1,s]],ris[x[1],x[2]]=[a[x[1],x[2]],b[x[1],x[2]]]);
ris;
}


{\\verifica che gli elementi di una matrice quadrata siano tutti distinti
verifica(m)=local(s,w);
s=length(m);
w=[];
forvec(x=[[1,s],[1,s]],w=concat(w,[m[x[1],x[2]]]));
if(length(Set(w))==s^2,return(1));
0;
}


{\\Dati p primo > 2 ed e > 0 calcola i p^e-1 quadrati ortogonali risultanti
\\dal campo finito
fflatin(p,e)=local();
ff=campop(p,e,x);
mm=matrix(p^e,p^e);
ris=vector(p^e-1);
for(f=2,p^e,forvec(x=[[1,p^e],[1,p^e]],mm[x[1],x[2]]=posizioni(ff,ff[f]*ff[x[1]]+ff[x[2]])[1]);ris[f-1]=mm);
lift(lift(ris));
}