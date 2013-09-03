#ifndef _GETOPT_H_
#define _GETOPT_H_

int checkparam_dic(char *param_dic);
void getparam_dic(char *param_dic,
		char *option_str, int *option_num);
int getopt(int argc, char *argv[],
		char *param_dic);

#endif
