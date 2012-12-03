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

/* vector(2^3,k,vector(3,i,((k-1)\2^(3-i))%2)) */

\\fqn(q, n) =
\\{
\\    return(vector(q^n,k,vector(n,i,((k-1)\q^(n-i))%q)));
\\}

fqn(p, e, n) =
{
    my(
        fq = campop(p, e, X),
        q = p^e,
        numbers = vector(q^n, k, vector(n, i, ((k-1)\q^(n-i)) % q)),
        v(indices) = vector(n, i, fq[indices[i] + 1]),
        fqn = apply(v, numbers)
    );
    return(fqn);
}
addhelp(fqn, "fqn(p,e,n): restituisce lo spazio vettoriale (F_q)^n con q=p^e.")

pincidenza(p, e) =
{
    my(
        fq3 = fqn(p, e, 3),
        scalar(v,w) = v*w~,
        q = p^e,
        point(i) = fq3[i+1],
        line(i) = fq3[i+1],
        pg2q = matrix(q^2+q+1, q^2+q+1, i, j, scalar(line(i),point(j)) == 0)
    );
    return(pg2q);
}
addhelp(pincidenza, "pincidenza(p,e): restituisce la matrice di incidenza del piano proiettivo generato dal campo F_q con q=p^e.")

/* Piano di Fano
 *
 * ? pincidenza(2,1)
 * [0 1 0 1 0 1 0]
 * [1 0 0 1 1 0 0]
 * [0 0 1 1 0 0 1]
 * [1 1 1 0 0 0 0]
 * [0 1 0 0 1 0 1]
 * [1 0 0 0 0 1 1]
 * [0 0 1 0 1 1 0]
 *
 * Come a pagina 15 del pdf (pagina 7 del documento)
 * http://www.dm.unito.it/personalpages/cerruti/aalab/Materiali/OnProjectivePlanes.pdf
 * con numerazione 1->4 3->6 e i rimanenti fissati.
 */


