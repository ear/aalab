/**
 *  Esercizi per il Laboratorio di Applicazioni dell'Algebra, 13/11/12
 */


/**
 * 1. Distribuzione dei periodi degli elementi di S_7
 */

es1(n, silenzioso) =
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
		until(p^k == id, k++); \\ appare più veloce che fordiv
		\\fordiv(ordine, d, k=d; if(p^k == id, break));
		if(!silenzioso, print(Vec(p), " ha ordine ", k));

		periodi = setunion(periodi, Set(k));
		occorrenze[k] += 1;
		permutazioni[k] = concat(permutazioni[k], [Vec(p)]);
	);

	periodi = vecsort(eval(periodi));
	occorrenze = vecextract(occorrenze, periodi);
	permutazioni = vecextract(permutazioni, periodi);
	if(!silenzioso, print("\nordini/occorrenze:\n"));
	return(concat(Mat(periodi)~, Mat(occorrenze)~));
}

\\ es1(7)

