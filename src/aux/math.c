#include <stdio.h>
#include <stdlib.h>

int trivial_log2_int(int x){
  int log2x = 0;

  while(x != 1){
    x /= 2;
    ++log2x;
  }

  return log2x;
}

int trivial_pow_int(int base, int power){
  int answer = base;

  if(power == 0)
    return 1;
  else
    for(int i=1; i < power; ++i)
      answer *= base;

  return answer;
}
