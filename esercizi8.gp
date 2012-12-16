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

h(n) =
{
    if(     n == 1 , return( Mat(1)               ),
    if(     n == 2 , return( [1,1;1,-1]           ),
    if( n % 2 == 0 , return( kron( h(2), h(n/2) ) ), \\ sub-optimal
        error("n in h(n) must be a power of 2.");
    )));
}


/**
 * 3. Si scriva hcodice(m) che genera il codice binario Ham(m) formato dalle
 * righe di H_m e -H_m trasformate in binario.
 */

hcodice(m) =
{
    Ham(2*m);
}

/* In [1] il parametro è la lunghezza del codice: il doppio della dimensione
 * della matrice.
 */

Ham(n, verbose=0) =
{
    my(
        n = n / 2,
        h_n = h(n),
        a = h_n % 3 % 2,
        b = h_n * 2 % 3 % 2
    );
    if( verbose,
      print("\n", "Hadamard code of order ", 2*n, " (block size ", n, ")", "\n",
            "Minimum distance: ", n/2, "\n",
            "Can correct: ", n/4 - 1, " errors", "\n");
    );
    return(concat(
        vector(#a, i, a[i,]),
        vector(#b, i, b[i,])
    ));
}


/**
 * Si scriva hdecod(w) che riceve una stringa binaria di lunghezza 2^k e tenta
 * di decodificarla, associandole una stringa di H_m, o dichiarando che ci sono
 * stati troppi errori. Si confrontino i risultati con gli esempi dati per H_8
 * in [1].
 */

hdecod(w_, verbose=0) =
{
    my(
        n = #w_,
        w = w_ * 2 + vector(n, i, -1),
        s = w * h(n),
        t = n/4 - 1,
        k = n - 2*t, \\ error treshold
        v
    );
    if( verbose,
      print("If there are no more than ", t, " errors:", "\n",
            "- the zero entries become of magnitude no larger than ", 2*t, "\n",
            "- the entry of magnitude n becomes no smaller than ", n-2*t, "\n");
    );
    for( i = 1, n,
        if(
            abs(s[i]) >= k
          ,
            v = hcodice(n)[( sign(s[i]) * 2 % 3 % 2 ) * n + i];
            if( verbose,
              print("|s[", i, "]| = |", s[i], "| >= ", k, "\n", "\n",
                    "Decode sucessful (", hamming(w_, v), " errors)", "\n");
            );
            return( v );
          ,
            if( verbose, print("|s[", i, "]| = |", s[i], "| <= ", k) );
        );
    );
    print("\n", "Decode failed: more than ", t, " errors have occurred.", "\n");
}


hamming(a, b) =
{
  my( s = 0 );
  if( #a != #b, error("different vector lengths in hamming(a,b)") );
  for( i = 1, #a, s += a[i] != b[i] );
  return( s );
}


