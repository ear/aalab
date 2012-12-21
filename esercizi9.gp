/**
 *  Esercizi per il Laboratorio di Applicazioni dell'Algebra, 18/12/12
 */


/**
 * 1. Scrivere una funzione kesava(n,q) che riceve q ed n coprimi e calcola
 * il numero dei fattori irriducibili di x^n - 1 in F_q[x], utilizzando
 * il Lemma di Burnside.
 */

f(n, q) =
{
    if( gcd(q,n) != 1, error("q and n must be coprime"));
    my( p = divisors(q)[2], \\ XXX: fragile way to get p from q = p^k
        i = 1
    );
    while( (p^i - 1) % n, i++ );
    return( i );
}

vecsum(v) =
{
    sum( i=1, #v, v[i] );
}

kesava(n, q) =
{
    my( p = divisors(q)[2], \\ XXX: only used in the factormod check (1)
        ds = divisors(n), \\ divisors
        os = apply( ((d) -> floor(eulerphi(d)/f(d, q))), ds )     \\ orders
    );
    print(length(factormod( x^n - 1, p )[,2])); \\ XXX: (1)
    return( vecsum(os) );
}


