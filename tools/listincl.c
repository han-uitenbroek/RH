#include <stdio.h>
#include <string.h>

#define MAX_LINE_SIZE  80
#define MAX_INCL_SIZE 400
#define MAX_NAME_SIZE  20

int main( int argc, char *argv[] )
{
  register int n;

  int    nFound, n0;
  char   inclStr[] = {"#include \""}, line[MAX_LINE_SIZE], *result,
         theIncl[MAX_INCL_SIZE], *ext, inclFileName[MAX_NAME_SIZE],
         dummy[MAX_LINE_SIZE], fileName[MAX_NAME_SIZE],
         libName[MAX_NAME_SIZE];
  FILE  *fp;

  if ((argc >= 3)  &&  strstr(argv[1], "-l")) {
    strcpy(libName, argv[2]);
    n0 = 3;
  } else
    n0 = 1;

  for (n = n0;  n < argc;  n++) {
    nFound = 0;
    strcpy(fileName, argv[n]);
    if (((ext = strstr(fileName, ".c")) == NULL)) {
      fprintf( stderr, "Not a .c file: %s\n", fileName );
      continue;
    }
    fp = fopen( fileName, "r" );

    *theIncl = '\0';
    while( (result = fgets(line, MAX_LINE_SIZE, fp)) != NULL ) {
      if (strstr(line, inclStr)) {
	nFound++;
        sscanf(line, "%s \"%s\"", dummy, inclFileName);
        strncat(theIncl, inclFileName, strlen(inclFileName) - 1);
	strcat(theIncl, "  ");
      }
    }
    fclose( fp );
    if (nFound) {
      if (n0 == 1) {
	strcpy(ext, ".o:");
	printf("%-18s%s\n", fileName, theIncl);
      } else if (n0 == 3) {
	strcpy(ext, ".o):");
	printf("%s(%-18s%s\n", libName, fileName, theIncl);
      }
    }
  }
}
