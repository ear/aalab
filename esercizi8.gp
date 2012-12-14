/**
 *  Esercizi per il Laboratorio di Applicazioni dell'Algebra, 11/12/12
 */


/**
 * 1. Scrivere una funzione kron(a, b), che riceve due matrici a, b e
 * restituisce il prodotto di Kronecker, o tensoriale, delle medesime.
 * Ricordiamo che se c = kron(a, b), si deve avere c_{i,j} = a_{i,j}*b.
 */

kron(a, b) =
{
    my(
        ra = #a[,1],
        ca = #a[1,],
        rb = #b[,1],
        cb = #b[1,],
        n = ra * rb,
        m = ca * cb,
        M(i,j) = a[ floor( (i-1) / rb ) + 1, floor( (j-1) / cb ) + 1 ] *
                 b[ ((i-1) % rb) + 1, ((j-1) % cb) + 1 ]
    );
    return( matrix( n, m, i, j, M(i,j)) );
}


/**
 * 2. Si usi kron(•, •) per scrivere una funzione che crei H_m, con m = 2^k,
 * come in [1].
 *
 * [1] http://www.dm.unito.it/personalpages/cerruti/aalab/Materiali/Hadamard.pdf
 */

h(m) =
{
    if( m == 1, return( Mat(1) ),
        if( m == 2, return( [1,1;1,-1] ),
                    return( kron( h(2), h(m-1) ) )
        )
    );
}


/**
 * 3. Si scriva hcodice(m) che genera il codice binario Ham(m) formato dalle
 * righe di H_m e -H_m trasformate in binario.
 */

hcodice(m) =
{
    my(
        h_m = h(m),
        a = h_m % 3 % 2,
        b = h_m * 2 % 3 % 2
    );
    return(concat(
        vector(#a, i, a[i,]),
        vector(#b, i, b[i,])
    ));
}


