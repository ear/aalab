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


