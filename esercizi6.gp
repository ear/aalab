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

es1(limite = 10000) =
{
    my(
        numeri = List(),
        polinomio(n,x) = x^(2*n)-x^n-1,
        radice(n) = Mod(Mod(1,n)*x, x^2-x-1),
        test(n) = lift(lift(polinomio(n,radice(n)))) == 0
    );

    for(n = 3, limite,
        if((n%5!=0) && !isprime(n) && test(n), listput(numeri, n))
    );

    return(Vec(numeri));
}

/* [2737, 4181, 5777, 6721]
 *
 * due osservazioni
 *
 *  (1) tutti i numeri trovati sono liberi da quadrati, è un caso?
 *      ho controllato fino a 2 milioni e lo sono tutti
 *  (2) il primo della lista con più di 3 fattori è 179697 come si può
 *      vedere con factor(es1(180000))
 */


/**
 * 2. generatori(p,e): restituisce l'elenco dei generatori di GF(p^e)
 * per rendere unico il risultato si ponga GF(p^e) = Z_p[x] / (f_p) dove
 * f_p(x) = primpoly(p,e)
 */

generatori(p, e) =
{
    my(f, y, ordine, coprimi, generatori);

    f = primpoly(p, e, x);
    y = Mod(Mod(1,p)*x, f);

    ordine = p^e - 1;
    coprimi = List();
    for(i=1, ordine, if(gcd(ordine,i) == 1, listput(coprimi,i)));
    coprimi = vecsort(Vec(coprimi));

    generatori = apply(i->y^i, coprimi);

    return(generatori);
}


/**
 * 3. fglatin(p,e): restituisce i  p^e - 1  quadrati latini ortogonali
 * costruiti a partire dal piano affine finito su GF(p,e)
 */

fglatin(p, e) =
{
    my(
        f = primpoly(p, e, x),
        g = Mod(Mod(1,p)*x, f),
        gf = concat(0, vector(p^e-1, i, g^i)),
        lt = apply(x->Str(lift(lift(x))), gf),
        l(g) = for(i=1, #lt, if(Str(lift(lift(g))) == lt[i], return(i))),
        a(k, i, j) = l(gf[k]*gf[i] + gf[j])
    );
    return(vector(p^e-1,k,matrix(p^e, p^e, i, j, a(1+k, i, j))));
}

