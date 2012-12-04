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
 * 3. pincidenza(p,e): restituisce la matrice di incidenza del piano proiettivo
 * generato dal campo F_q, con q = p^e.
 */

pincidenza(p, e) =
{
    my(
        q = p^e,
        fq = campop(p, e, X),

        pairs(xs) = concat(vector(length(xs),i,vector(length(xs),j,[xs[i],xs[j]]))),

        points = concat([
            apply((yz -> concat(fq[q],yz)), pairs(fq)),  \\ [1,y,z]
            apply(( z -> [0,fq[q],z]     ),        fq),  \\ [0,1,z]
            [            [0,0,fq[q]]                 ]   \\ [0,0,1]
        ]),
        point(i) = points[i],
        line(i) = points[i],

        scalar(v,w) = v*w~,
        pg2q = matrix(q^2+q+1, q^2+q+1, i, j, scalar(line(i),point(j)) == 0)
    );
    return(pg2q);
}
addhelp(pincidenza, "pincidenza(p,e): restituisce la matrice di incidenza del piano proiettivo generato dal campo F_q con q=p^e.")


