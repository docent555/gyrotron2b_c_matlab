/*
 * main.c
 *
 * Code generation for function 'main'
 *
 */

 /*************************************************************************/
 /* This automatically generated example C main file shows how to call    */
 /* entry-point functions that MATLAB Coder generated. You must customize */
 /* this file for your application. Do not modify this file directly.     */
 /* Instead, make a copy of this file, modify it, and integrate it into   */
 /* your development environment.                                         */
 /*                                                                       */
 /* This file initializes entry-point function arguments to a default     */
 /* size and value before calling the entry-point functions. It does      */
 /* not store or use any values returned from the entry-point functions.  */
 /* If necessary, it does pre-allocate memory for returned values.        */
 /* You can use this file as a starting point for a main function that    */
 /* you can deploy in your application.                                   */
 /*                                                                       */
 /* After you copy the file, and before you deploy it, you must make the  */
 /* following changes:                                                    */
 /* * For variable-size function arguments, change the example sizes to   */
 /* the sizes that your application requires.                             */
 /* * Change the example values of function arguments to the values that  */
 /* your application requires.                                            */
 /* * If the entry-point functions return values, store these values or   */
 /* otherwise use them as required by your application.                   */
 /*                                                                       */
 /*************************************************************************/

 /* Include files */
#include <stdio.h>
#include "main.h"
#include "gyrotron.h"
#include "gyrotron_emxAPI.h"
#include "gyrotron_terminate.h"
#include "rt_nonfinite.h"

static void main_gyrotron(int Ne, double z0, double zex, double zin, double Tend, double Delta,
	double I0, double dz, double dt, double tol, int INTT, int INTZ, double a0, int nexl, int nexr, int mode, double akp2)
{
	emxArray_real_T* ConLow;
	emxArray_real_T* OUTFre;
	emxArray_real_T* OUTFim;
	emxArray_real_T* OUTJre;
	emxArray_real_T* OUTJim;
	emxArray_real_T* OUTZAxis;
	emxArray_real_T* OUTTAxis;
	emxArray_real_T* Eff;
	emxArray_real_T* Omega;
	emxArray_real_T* pim;
	emxArray_real_T* pre;
	double jout;
	emxInitArray_real_T(&OUTFre, 2);
	emxInitArray_real_T(&OUTFim, 2);
	emxInitArray_real_T(&OUTJre, 2);
	emxInitArray_real_T(&OUTJim, 2);
	emxInitArray_real_T(&pre, 2);
	emxInitArray_real_T(&pim, 2);
	emxInitArray_real_T(&ConLow, 1);
	emxInitArray_real_T(&OUTZAxis, 1);
	emxInitArray_real_T(&OUTTAxis, 1);
	emxInitArray_real_T(&Eff, 1);
	emxInitArray_real_T(&Omega, 1);

	FILE* fre, * fim, * ire, * iim, * feff, * fomega, * fileID, * fpre, * fpim, *fcl;
	errno_t err;	
	double Ic, TauEnd, Lz, Lzi;	

	int Nz1, Nz2, Nt, OUTNz, OUTNt, idx0, idx1;
	static int iv2[2], piv[2];

	Lz = zex - z0;
	Lzi = zin - z0;

	Ic = I0;
	Nz1 = (int)(Lz / dz) + 1;
	Nz2 = (int)(Lzi / dz) + 1;
	Nt = (int)(Tend / dt) + 1;	
	TauEnd = Tend;		

	if (INTT < 1) {
		printf("Too small Tend\n");
		exit(-1);
	}
	if (INTZ < 1) {
		printf("Too small Zend\n");
		exit(-1);
	}

	if (INTT > 1 && INTZ > 1) {
		OUTNt = (int)((Nt - 1) / INTT) + 1;
		OUTNz = (int)((Nz1 - 1) / INTZ) + 1;
	}
	else if (INTT == 1 && INTZ > 1) {
		OUTNt = Nt;
		OUTNz = (int)((Nz1 - 1) / INTZ) + 1;
	}
	else if (INTT > 1 && INTZ == 1) {
		OUTNt = (int)((Nt - 1) / INTT) + 1;
		OUTNz = Nz1;
	}
	else {
		OUTNt = Nt;
		OUTNz = Nz1;
	}

	iv2[0] = OUTNz;
	iv2[1] = OUTNt;
	piv[0] = Nz2;
	piv[1] = Ne;

	OUTFre = emxCreateND_real_T(2, iv2);
	OUTFim = emxCreateND_real_T(2, iv2);
	OUTJre = emxCreateND_real_T(2, iv2);
	OUTJim = emxCreateND_real_T(2, iv2);
	pre = emxCreateND_real_T(2, piv);
	pim = emxCreateND_real_T(2, piv);
	OUTZAxis = emxCreateND_real_T(1, &OUTNz);
	OUTTAxis = emxCreateND_real_T(1, &OUTNt);
	Eff = emxCreateND_real_T(1, &OUTNt);
	Omega = emxCreateND_real_T(1, &OUTNt);
	
	err = fopen_s(&fileID, "input_c.txt", "w");
	fprintf(fileID, "Nz1 = % i\n", Nz1);
	fprintf(fileID, "Nz2 = % i\n", Nz2);
	fprintf(fileID, "Nt = % i\n", Nt);
	fprintf(fileID, "Ne = % i\n", (int)Ne);
	fprintf(fileID, "Lz = % f\n", Lz);
	fprintf(fileID, "Lzi = % f\n", Lzi);
	fprintf(fileID, "TauEnd = % f\n", TauEnd);
	fprintf(fileID, "Delta = % f\n", Delta);
	fprintf(fileID, "I0 = % f\n", I0);
	fprintf(fileID, "dz = % f\n", dz);
	fprintf(fileID, "dt = % f\n", dt);
	fprintf(fileID, "tol = % g\n", tol);
	fprintf(fileID, "INTT = % i\n", (int)INTT);
	fprintf(fileID, "INTZ = % i\n", (int)INTZ);
	fprintf(fileID, "a0 = % f\n", a0);
	fprintf(fileID, "nexl = % i\n", nexl);
	fprintf(fileID, "nexr = % i\n", nexr);
	fprintf(fileID, "mode = % i\n", mode);
	fprintf(fileID, "akp2 = % f\n", akp2);
	fclose(fileID);

	/* Call the entry-point 'orotron'. */
	gyrotron(Ne, z0, zex, zin, Tend, Delta, I0, dz, dt, tol, INTT, INTZ, a0, nexl, nexr, mode, akp2,
		OUTFre, OUTFim, OUTJre, OUTJim, pre, pim, ConLow, OUTZAxis, OUTTAxis, Eff, Omega, &jout);

	err = fopen_s(&fpre, "pre.dat", "w");
	if (err == 0)
	{
		printf("The file 'pre.dat' was opened\n");
	}
	else
	{
		printf("The file 'pre.dat' was not opened\n");
	}

	err = fopen_s(&fpim, "pim.dat", "w");
	if (err == 0)
	{
		printf("The file 'pim.dat' was opened\n");
	}
	else
	{
		printf("The file 'pim.dat' was not opened\n");
	}

	err = fopen_s(&fre, "fre.dat", "w");
	if (err == 0)
	{
		printf("The file 'fre.dat' was opened\n");
	}
	else
	{
		printf("The file 'fre.dat' was not opened\n");
	}

	err = fopen_s(&fim, "fim.dat", "w");
	if (err == 0)
	{
		printf("The file 'fim.dat' was opened\n");
	}
	else
	{
		printf("The file 'fim.dat' was not opened\n");
	}

	err = fopen_s(&ire, "ire.dat", "w");
	if (err == 0)
	{
		printf("The file 'ire.dat' was opened\n");
	}
	else
	{
		printf("The file 'ire.dat' was not opened\n");
	}

	err = fopen_s(&iim, "iim.dat", "w");
	if (err == 0)
	{
		printf("The file 'iim.dat' was opened\n");
	}
	else
	{
		printf("The file 'iim.dat' was not opened\n");
	}

	err = fopen_s(&feff, "e.dat", "w");
	if (err == 0)
	{
		printf("The file 'e.dat' was opened\n");
	}
	else
	{
		printf("The file 'e.dat' was not opened\n");
	}

	err = fopen_s(&fomega, "w.dat", "w");
	if (err == 0)
	{
		printf("The file 'w.dat' was opened\n");
	}
	else
	{
		printf("The file 'w.dat' was not opened\n");
	}

	err = fopen_s(&fcl, "cl.dat", "w");
	if (err == 0)
	{
		printf("The file 'cl.dat' was opened\n");
	}
	else
	{
		printf("The file 'cl.dat' was not opened\n");
	}

	for (idx0 = 0; idx0 < Nz2; idx0++) {
		fprintf(fpre, "%e\t", idx0 * dz);
		for (int idx1 = 0; idx1 < Ne; idx1++) {
			fprintf(fpre, "%e\t", pre->data[idx0 + pre->size[0] * idx1]);
		}
		fprintf(fpre, "\n");
	}

	for (idx0 = 0; idx0 < Nz2; idx0++) {
		fprintf(fpim, "%e\t", idx0 * dz);
		for (idx1 = 0; idx1 < Ne; idx1++) {
			fprintf(fpim, "%e\t", pim->data[idx0 + pim->size[0] * idx1]);
		}
		fprintf(fpim, "\n");
	}

	for (idx0 = 0; idx0 < OUTNz; idx0++) {
		fprintf(fre, "%e\t", OUTZAxis->data[idx0]);
		for (idx1 = 0; idx1 < OUTNt; idx1++) {
			fprintf(fre, "%e\t", OUTFre->data[idx0 + OUTFre->size[0] * idx1]);
		}
		fprintf(fre, "\n");
	}

	for (idx0 = 0; idx0 < OUTNz; idx0++) {
		fprintf(fim, "%e\t", OUTZAxis->data[idx0]);
		for (idx1 = 0; idx1 < OUTNt; idx1++) {
			fprintf(fim, "%e\t", OUTFim->data[idx0 + OUTFim->size[0] * idx1]);
		}
		fprintf(fim, "\n");
	}

	for (idx0 = 0; idx0 < OUTNz; idx0++) {
		fprintf(ire, "%e\t", OUTZAxis->data[idx0]);
		for (idx1 = 0; idx1 < OUTNt; idx1++) {
			fprintf(ire, "%e\t", OUTJre->data[idx0 + OUTJre->size[0] * idx1]);
		}
		fprintf(ire, "\n");
	}

	for (idx0 = 0; idx0 < OUTNz; idx0++) {
		fprintf(iim, "%e\t", OUTZAxis->data[idx0]);
		for (idx1 = 0; idx1 < OUTNt; idx1++) {
			fprintf(iim, "%e\t", OUTJim->data[idx0 + OUTJim->size[0] * idx1]);
		}
		fprintf(iim, "\n");
	}

	for (int idx1 = 0; idx1 < Nt; idx1++) {
		fprintf(feff, "%e\t%e\n", idx1*dt, Eff->data[idx1]);
	}

	for (int idx1 = 0; idx1 < Nt; idx1++) {
		fprintf(fomega, "%e\t%e\n", idx1*dt, Omega->data[idx1]);
	}

	for (int idx1 = 0; idx1 < Nt; idx1++) {
		fprintf(fcl, "%e\t%e\n", idx1 * dt, ConLow->data[idx1]);
	}


	if (ire)
	{
		err = fclose(ire);
		if (err == 0)
		{
			printf("The file 'ire.dat' was closed\n");
		}
		else
		{
			printf("The file 'ire.dat' was not closed\n");
		}
	}

	if (iim)
	{
		err = fclose(iim);
		if (err == 0)
		{
			printf("The file 'iim.dat' was closed\n");
		}
		else
		{
			printf("The file 'iim.dat' was not closed\n");
		}
	}

	if (fre)
	{
		err = fclose(fre);
		if (err == 0)
		{
			printf("The file 'fre.dat' was closed\n");
		}
		else
		{
			printf("The file 'fre.dat' was not closed\n");
		}
	}

	if (fim)
	{
		err = fclose(fim);
		if (err == 0)
		{
			printf("The file 'fim.dat' was closed\n");
		}
		else
		{
			printf("The file 'fim.dat' was not closed\n");
		}
	}

	if (feff)
	{
		err = fclose(feff);
		if (err == 0)
		{
			printf("The file 'e.dat' was closed\n");
		}
		else
		{
			printf("The file 'e.dat' was not closed\n");
		}
	}

	if (fomega)
	{
		err = fclose(fomega);
		if (err == 0)
		{
			printf("The file 'w.dat' was closed\n");
		}
		else
		{
			printf("The file 'w.dat' was not closed\n");
		}
	}

	if (fcl)
	{
		err = fclose(fcl);
		if (err == 0)
		{
			printf("The file 'cl.dat' was closed\n");
		}
		else
		{
			printf("The file 'cl.dat' was not closed\n");
		}
	}

	emxDestroyArray_real_T(Omega);
	emxDestroyArray_real_T(Eff);
	emxDestroyArray_real_T(OUTTAxis);
	emxDestroyArray_real_T(OUTZAxis);
	emxDestroyArray_real_T(pre);
	emxDestroyArray_real_T(pim);
	emxDestroyArray_real_T(ConLow);
	emxDestroyArray_real_T(OUTJim);
	emxDestroyArray_real_T(OUTJre);
	emxDestroyArray_real_T(OUTFim);
	emxDestroyArray_real_T(OUTFre);
}

int main(int argc, const char* const argv[])
{
	int Ne, INTT, INTZ, mode, nexl, nexr;
	double z0, zex, zin, Tend, Delta, I0, dz, dt, tol, a0, akp2;

	if (argc != 18) {
		printf("Expected 17 arguments: Ne, z0, zex, zin, Tend, Delta, I0, dz, dt, tol, INTT, INTZ, a0, nexl, nexr, MODE, akp2\n");
		exit(-1);
	}

	Ne = atoi(argv[1]);	
    z0 = atof(argv[2]);
	zex = atof(argv[3]);
	zin = atof(argv[4]);
	Tend = atof(argv[5]);
	Delta = atof(argv[6]);
	I0 = atof(argv[7]);	
	dz = atof(argv[8]);
	dt = atof(argv[9]);
	tol = atof(argv[10]);
	INTT = atoi(argv[11]);
	INTZ = atoi(argv[12]);		
	a0 = atof(argv[13]);	
	nexl = atoi(argv[14]);
	nexr = atoi(argv[15]);
	mode = atoi(argv[16]);
	akp2 = atof(argv[17]);

	/* The initialize function is being called automatically from your entry-point function. So, a call to initialize is not included here. */
	/* Invoke the entry-point functions.
	   You can call entry-point functions multiple times. */
	main_gyrotron(Ne, z0, zex, zin, Tend, Delta, I0, dz, dt, tol, INTT, INTZ, a0, nexl, nexr, mode, akp2);

	/* Terminate the application.
	   You do not need to do this more than one time. */
	gyrotron_terminate();
	return 0;
}

/* End of code generation (main.c) */
