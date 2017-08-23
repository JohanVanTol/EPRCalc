#include <stdio.h>
#include <values.h>

int read_int(FILE *fichier)
{
	char line[50];
	int n, ok=0;

	while (!ok)
	if (fichier != NULL)
	{
		if (fscanf(fichier,"%s", line) != EOF)
		{
			if (((line[0] > 47) && (line[0] < 58))          // first char 0..9
				 || (line[0] == '+') || (line[0] == '-')) // or + or - or .
				ok += sscanf(line, "%d", &n);      // read int from string
		}
			else
		return MAXINT;
	}
	else
	{
		if (scanf("%s", line) != EOF)
		{
			if (((line[0] > 47) && (line[0] < 58))          // first char 0..9
				 || (line[0] == '+') || (line[0] == '-')) // or + or - or .
				ok += sscanf(line, "%d", &n);      // read int from string
		}
			else
		return MAXINT;
	}
	return n;
}

unsigned int read_uint(FILE *fichier)
{
	char line[50];
	unsigned int n;
	int ok=0;

	while (!ok)
	if (fichier != NULL)
	{
		if (fscanf(fichier,"%s", line) != EOF)
		{
			if (((line[0] > 47) && (line[0] < 58))          // first char 0..9
				 || (line[0] == '+') || (line[0] == '-')) // or + or - or .
				ok += sscanf(line, "%ud", &n);      // read unsigned int from string
		}
			else
		return MAXINT;
	}
	else
	{
		if (scanf("%s", line) != EOF)
		{
			if (((line[0] > 47) && (line[0] < 58))          // first char 0..9
				 || (line[0] == '+') || (line[0] == '-')) // or + or - or .
				ok += sscanf(line, "%ud", &n);      // read unsigned int from string
		}
			else
		return MAXINT;
	}
	return n;
}

double read_double(FILE *fichier)
{
	char line[50];
	int ok=0;
	double d;
	while (!ok)
	if (fichier != 0)
	{
		if (fscanf(fichier,"%s", line) != EOF)
		{
			if (((line[0] > 47) && (line[0] < 58))          // first char 0..9
				 || (line[0] == '+') || (line[0] == '-')
				 || (line[0] == '.'))                     // or + or - or .
				ok += sscanf(line, "%lf", &d);      // read double from string
		}
		 else return MAXDOUBLE;
	}
	else
	{
		if (scanf("%s", line) != EOF)
		{
			if (((line[0] > 47) && (line[0] < 58))          // first char 0..9
				 || (line[0] == '+') || (line[0] == '-')
				 || (line[0] == '.'))                     // or + or - or .
				ok += sscanf(line, "%lf", &d);      // read double from string
		}
		 else return MAXDOUBLE;
	}
	return d;
}
