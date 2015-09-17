
#include <R.h>

void hello(int *n){
	int i;
	for(i=0; i < *n; i++) {
	Rprintf("Hello, world!\n");
	}
}


void helloname(char **name){
	int i;
	Rprintf("Hello, %s!\n",name[0]);
}
