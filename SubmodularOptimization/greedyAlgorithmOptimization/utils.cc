#include <string.h>
#include "utils.h"
#include "error.h"

int
checkline(char *line,unsigned num) {

  if (strlen(line) == MAXLEN && line[strlen(line)-1] != '\n') {
    error("Line %li too long\n",num);
  }
  else {
    line[strlen(line)-1] = '\0';
  }
  return strlen(line);
}




  

  
  

