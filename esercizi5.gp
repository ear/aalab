/**
 *  Esercizi per il Laboratorio di Applicazioni dell'Algebra, 13/11/12
 */


/**
 *  sn(n): matrice m di 3 colonne che da la struttura di S_n.
 *  m[,1] sono i periodi delle permutazioni,
 *  m[,2] sono il numero di permutazioni di ogni periodo,
 *  m[,3] sono le permutazioni.
 */

sn(n) =
{
	local(periodi, occorrenze, permutazioni);

	ordine = n!;

	periodi = Set();
	occorrenze = vector(sigma(ordine, 0));
	permutazioni = vector(sigma(ordine, 0), i, []);

	for(i = 1, ordine,
		p = Vecsmall(numtoperm(n, i));
		id = Vecsmall(vector(n, i, i));

		k = 1;
		until(p^k == id, k++);

		periodi = setunion(periodi, Set(k));
		occorrenze[k] += 1;
		permutazioni[k] = concat(permutazioni[k], [Vec(p)]);
	);

	periodi = vecsort(eval(periodi));
	occorrenze = vecextract(occorrenze, periodi);
	permutazioni = vecextract(permutazioni, periodi);
	return(concat(concat(Mat(periodi)~, Mat(occorrenze)~), Mat(permutazioni)~));
}


/**
 * 1. Distribuzione dei periodi degli elementi di S_7
 */

es1() =
{
	local(s7);
	s7 = sn(7);

	return(Mat([s7[,1], s7[,2]]));
}


/**
 * 2. Determinare quanti e quali sono gli elementi di periodo dispari in S_7
 */

es2() =
{
	local(s7, n);
	s7 = sn(7);

	periodi = s7[,1]~;
	occorrenze = s7[,2]~;

	n = 0;
	for(i = 1, length(periodi), if(periodi[i]%2 == 1, n += occorrenze[i]));

	return(n);
}


/**
 * 3. Generatori di F_(3^4)[x]/(x^4+x^3+x^2+x+1) =~ GF(81)^* =~ C_80  (54 = eulerphi(80)).
 */

es3() =
{
	local(ps, gs);

	ps = Set();
	forvec(X = [[0,2], [0,2], [0,2], [0,2], [0,2]],
		p = Mod(X, 3)*[1, yy, yy^2, yy^3, yy^4]~;
		ps = setunion(ps, Set(p))
	);
	ps = eval(vecsort(ps));

	gs = Set();
	for(i = 1, length(ps), if(cercaperiodo(ps[i], 80) == 80, gs=setunion(gs, Set(ps[i]))));
	gs = eval(vecsort(gs));

	print("trovati ", length(gs), " polinomi che generano GF(81)^*");
	for(i = 1, length(gs), print(lift(lift(gs[i]))));
}

