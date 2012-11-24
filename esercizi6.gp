/**
 *  Esercizi per il Laboratorio di Applicazioni dell'Algebra, 20/11/12
 */


/**
 * 1. Se si pone y = Mod( Mod(1,p) * x, x^2 - x - 1 ) con p primo, si ha
 * ovviamente, in Pari, y^2 - y - 1 = 0
 *
 * Spiegare perché, nelle stesse ipotesi, si ha che y^(2p) - y^p - 1 = 0.
 *
 *   In caratteristica p abbiamo (a+b)^p = a^p + b^p, quindi
 *
 *   y^(2p) - y^p - 1  =  (y^2 - y - 1)^p  =  0^p  =  0
 *
 * Trovate quali interi p, non divisibili per 5, tra 3 e 10000, soddisfano
 * la  y^(2p) - y^p - 1 = 0  e non sono primi.
 */

polinomio(n,x) = x^(2*n) - x^n - 1;

campo(n) = Mod( Mod(1,n) * x, x^2 - x - 1 );

test(n) = lift(lift(polinomio(n,campo(n)))) == 0;

es1() =
{
	local(ns);

	ns = List();
	for(n = 3, 10000, if((n%5!=0) && !isprime(n) && test(n), listput(ns, n)));
	ns = Vec(ns);

	return(ns);
}

/* [2737, 4181, 5777, 6721]
 *
 * usando allocatemem() tre volte e alzando il limite a 200000 ho trovato che;
 * (1) tutti questi numeri sono liberi da quadrati
 * (2) il primo della lista con più di 3 fattori è 179697
 */


/**
 * 2. generatori(p,e): restituisce l'elenco dei generatori di GF(p^e)
 * per rendere unico il risultato si ponga GF(p^e) = Z_p[x] / (f_p) dove
 * f_p(x) = primpoly(p,e)
 */

generatori(p, e) =
{
	local(f, y, coprimi, generatori);

	f = primpoly(p, e, x);
	y = Mod(Mod(1,p)*x, f);

	coprimi=List();
	for(i=2, 80, if(gcd(80,i)==1, listput(coprimi,i)));
	coprimi = vecsort(Vec(coprimi));

	generatori = apply((i)->y^i, coprimi);

	return(generatori);
}

