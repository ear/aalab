/**
 *  Esercizi per il Laboratorio di Applicazioni dell'Algebra, 18/12/12
 */


/**
 * 1. Scrivere una funzione kesava(n,q) che riceve q ed n coprimi e calcola
 * il numero dei fattori irriducibili di x^n - 1 in F_q[x], utilizzando
 * il Lemma di Burnside.
 */

polcyclofactorsnum(n, q) =
{
    my( p = 0 );
    if( isprime( q ), p = q,
    if( !ispower( q, , &p ) || !isprime( p ), error("q is not a prime power") );
    );

    \\ the smallest i such that q^i - 1 = 0 (mod n) is the size of every orbit
    my( i = 1 );
    while( (q^i - 1) % n, i++ );

    \\ the degree of the n-th cyclotomic polynomial is eulerphi(n)
    my( d = eulerphi(n) );

    \\ thus the number of its irreducible factors is
    return( d \ i );
}
addhelp(polcyclofactorsnum, "polcyclofactorsnum(n,p): number of irreducible factors of the n-th cyclotomic polynomial over F_p[x].");

vecsum(v) =
{
    sum( i=1, #v, v[i] );
}

kesava(n, q, verbose=0) =
{
    my( p = 0 );
    if( isprime( q ), p = q,
    if( !ispower( q, , &p ) || !isprime( p ), error("q is not a prime power") );
    );

    if( gcd( n, q ) != 1 , error("n and q must be coprime") );

    /*  x^n - 1 factors into a product of d-th cyclotomic polynomials where  *
     *  d runs on the divisors of n, so here is a vector of the number of    *
     *  the irreducible factors of said cyclotomic polynomials over F_q[x]   */

    my( numirr = apply( ((d) -> polcyclofactorsnum(d, q)), divisors(n) ) );

    if( verbose, print( divisors(n), "\n", numirr ) );

    return( vecsum( numirr ) );
}
addhelp(kesava, "kesava(n,q): number of irreducible factors of x^n - 1 over F_q[x] for coprime n, q.");

/* Some tests */

kesava_table(p, max_n=30, max_e=10) =
{
    matrix(max_e, max_n, i, j, if( gcd(j, p^i) == 1, kesava( j, p^i ) ));
}
addhelp(kesava_table, "kesava_table(p,max_n=30,max_e=10): the number of irreducible factors of the n-th cyclotomic polynomial (1 ≤ n ≤ max_n, rows) over F_(p^e)[x] (1 ≤ e ≤ max_e, columns).")

numirrpol(n, q) =
{
    return(
        vecsum( apply( ( (d) -> moebius(n/d) * q^d ), divisors(n) ) ) / n
    );
}
addhelp(numirrpol, "numirrpol(n,q): number of monic irreducible polynomials of degree n over F_q[x].")


/**
 * 2. Scrivere una funzione laterali(n,q), che calcola i laterali ciclotomici,
 * corrispondenti alle classi ciclotomiche delle radici n-esime dell'unità su
 * F_q.
 */

laterali(n, q, verbose=0) =
{
    my( c1 = vector(n, i, -1) );
    if( verbose, print("i -> cyclotomic cosets") );

    /* As shown by calling the function with the verbose parameter 1, it does
     * do some more work than strictly needed, in exchange for the minimal
     * complete set of representatives of cyclotomic cosets of q modulo n.
     * (An optimization is possible for n = q^m - 1 for some m.)
     */

    forstep( i = n-1, 0, -1,
        c1[ (i*q^0 % n) + 1 ] = i;
        my( j = 1, iq0 = i*q^0 );
        while( (i*q^j % n) != iq0,
            c1[ (i*q^j % n) + 1 ] = i;
            j++;
        );
        if( verbose, print(i, " -> ", c1) );
    );

    return( c1 );
}

laterale(n, q, l) =
{
    my(
        powers = vector(n, k, 2^(k-1)),
        bitmask(bits) = bits*powers~,
        cosets = laterali(n, q),
        i = cosets[ l + 1 ],
        coset(i) = vecextract(
            vector(n,i,i-1),
            bitmask(
                apply( ((x) -> x == i), cosets )
            )
        )
    );
    return( coset(i) );
}


/**
 * 3. Scrivere una funzione irrnql(n,p,l) che, dati n,p,l con 0 ≤ l ≤ n-1,
 * calcola il fattore irriducibile di x^n - 1 in F_p[x] corrispondente all'unico
 * laterale ciclotomico a cui l appartiene.
 */

vecprod(v) =
{
    prod( i=1, #v, v[i] );
}

irrnql(n, p, l) =
{
    my( h = 1 );
    while( (p^h - 1) % n, h++ );
    my( f = primpoly(p, h, t),
        m = (p^h - 1) / n,
        alpha = Mod(t^m, f)
    );
    my(
        coset = laterale(n, p, l),
        factors = apply( ((i) -> x - alpha^i), coset )
    );
    return( vecprod( factors ) );
}


/**
 * 4. Scrivere una funzione circolgrorder(n,q) che restituisce l'ordine del
 * gruppo moltiplicativo di R_{n,q} = F_q[x] / ( x^n - 1 ).
 */

circolgrorder(n, q) =
{
    /* R a ring w/ 1, A_1,…,A_n pairwise coprime ideals
     * (i.e. i =/= j => (A_i) + (A_j) = R) then
     *
     *     R / ( A_1 ∩ … ∩ A_n )  ≃  (R / A_1) ⨯ … ⨯ (R / A_n)
     *
     * in particular for R = F_p[x] and A_i = (m_i(x))
     * where m_i runs through the irreducible factors of x^n - 1 over F_p[x]
     * i.e. i runs through a complete system of cyclotomic coset representatives
     * each factor of the product of rings is a field
     * and the number of invertibles of R_{n,q} is
     *
     *     \prod_{i} deg( m_i(x) ) - 1  =  \prod_{i} |C_i| - 1
     *
     * where C_i is the unique cyclotomic coset of q modulo n containing i
     */

    my(
        cosets = vecsort(laterali(n,q)), \\ sort the system of representatives
        p = -1, \\ fake representative to skip 0
        t = 0 \\ accumulate (for each representative) sum of repetitions - 1
    );
    for( i=1, #cosets,
        if( cosets[i] == p, t++ );
        p = cosets[i];
    );
    return( t );
}


