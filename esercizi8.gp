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
        M(i,j) = a[ floor( (i-1) / ra ) + 1, floor( (j-1) / ca ) + 1 ] *
                 b[ ((i-1) % rb) + 1, ((j-1) % cb) + 1 ]
    );
    return( matrix( n, m, i, j, M(i,j)) );
}

