/**
 *  Esercizi per il Laboratorio di Applicazioni dell'Algebra, 27/11/12
 */


/**
 * 1. da esercizi6.gp:
 *   3. fglatin(p,e): restituisce i  p^e - 1  quadrati latini ortogonali
 *   costruiti a partire dal piano affine finito su GF(p,e).
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
addhelp(fglatin, "fglatin(p,e): restituisce i p^e - 1 quadrati latini ortogonali costruiti a partire dal piano affine finito su GF(p,e).")


/**
 * 2. aincidenza(p,e): restituisce la matrice di incidenza del piano affine
 * generato dal campo F_q, con q = p^e.
 */

aincidenza(p, e) =
{
    my(
        m = pincidenza(p, e),

        q = p^e,
        n = q^2 + q + 1,

        rows = vector(n-1, i, i),

        l = m[n,],
        powers = vector(n, k, 2^(k-1)),
        bitmask(bits) = bits*powers~,
        columns = bitmask(apply((x) -> !x, l)),

        paq = vecextract(m, rows, columns)
    );
    return(paq);
}
addhelp(aincidenza, "aincidenza(p,e): restituisce la matrice di incidenza del piano affine finito generato dal campo F_q con q=p^e.")

/**
 * 3. pincidenza(p,e): restituisce la matrice di incidenza del piano proiettivo
 * generato dal campo F_q, con q = p^e.
 */

pincidenza(p, e) =
{
    my(
        q = p^e,
        fq = campop(p, e, X),

        pairs(xs) = concat(vector(#xs, i, vector(#xs, j, [xs[i],xs[j]]))),

        points = concat([
            apply((yz -> concat(fq[q],yz)), pairs(fq)),  \\ [1,y,z]
            apply(( z -> [0,fq[q],z]     ),        fq),  \\ [0,1,z]
            [            [0,0,fq[q]]                 ]   \\ [0,0,1]
        ]),
        point(i) = points[i],
        line(i) = points[i],

        scalar(v,w) = v*w~,
        pg2q = matrix(q^2+q+1, q^2+q+1, i, j, scalar(line(i), point(j)) == 0)
    );
    return(pg2q);
}
addhelp(pincidenza, "pincidenza(p,e): restituisce la matrice di incidenza del piano proiettivo finito generato dal campo F_q con q=p^e.")


/**
 * 5. In un esame vengono assegnati agli studenti complessivamente 36 problemi.
 * Ad ogni studente sono dati 6 problemi, e ogni coppia di problemi va
 * esattamente a due studenti. Quanti sono gli studenti?
 *
 *   v=36, k=6, l=2: r(6-1) = 2(36-1)
 *
 *   r = 2*35/5 = 14 studenti
 */


/**
 * 6. Ci sono 13 persone. È possibile creare 13 gruppi in modo tale che:
 *  (1) Ogni gruppo contiene 4 persone.
 *  (2) Due gruppi hanno esattamente una persona in comune.
 *  (3) Due persone hanno un solo gruppo in comune.
 * Motivare la risposta teoricamente.
 * Nel caso di risposta sì:
 *  4) A quanti gruppi appartiene una persona?
 *  5) Si dia una risposta effettiva, identificando le persone con
 *     i numeri 1, …, 13.
 *
 *  Siamo alla ricerca di un (13,4,1)-design (v = 13, k = 4, l = 1).
 *
 *  Questo corrisponde al piano proiettivo finito di ordine 3 PG(2,3)
 *  perché in esso abbiamo:
 *
 *    (13 persone) 3^2 + 3 + 1 = 13 rette
 *    (1) 3 + 1 = 4 punti su ogni retta
 *    (2) ogni coppia di rette si interseca in un punto
 *    (3) ogni coppia di punti individua un'unica retta
 *
 *  4) Ogni persona appartiene a 4 gruppi, perché per ogni punto di PG(2,3)
 *  passano 4 rette. (Oppure perché ogni elemento dell'insieme base compare
 *  in r blocchi, con r(k-1) = l(v-1), r = (13-1)/(4-1) = 4.)
 *
 *  5) ? \r aalab-funzioni
 *     ? \r esercizi7
 *     ? tredici()
 *     %37 = [[10, 11, 12, 13], [3, 6, 9, 10], [2, 5, 8, 10], [7, 8, 9, 13],
 *               [3, 5, 7, 11], [2, 6, 7, 12], [4, 5, 6, 13], [3, 4, 8, 12], 
 *               [2, 4, 9, 11], [1, 2, 3, 13], [1, 5, 9, 12], [1, 6, 8, 11],
 *               [1, 4, 7, 10]]
 *
 */

tredici() =
{
    my(
        pg23    = pincidenza(3, 1),
        numbers = vector(13, i, i),
        powers  = vector(13, k, 2^(k-1)),
        bitmask(bits) = bits*powers~
    );
    return(vector(13, i, vecextract(numbers, bitmask(pg23[i,]))));
}


