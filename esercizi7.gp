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
