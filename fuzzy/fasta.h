
typedef struct
{
	char *name;
	char *name2;
	char *path;
	char *SCFname;
	int  length;
	char *data;
	char *qual;
	int  finished;
} fasta;

typedef struct
{
	int  length;
	char *data;
} fast;

fasta *decodeFastq (char *fname, int *nContigs, long *tB, char* pdata, long Size_pdata,fasta *segg);
fast *readFastq (char **argv, int argc, int *nContigs);
fasta *splitFastq ( fasta *iseg, int inSeg, int tB, int *nReads, int length, int step);
void fastaLC (fasta *seg, int nSeg);
void fastaUC (fasta *seg, int nSeg);

