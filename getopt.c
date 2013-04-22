/*
 written by JCh
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char *opt_str;
char *opt_arr[16];

int checkparam_dic(char *param_dic)
{
	int param_num=0,i=-1;
	if ( param_dic[0]<'A' || \
			(param_dic[0]>'Z' && param_dic[0]<'a') || \
			param_dic[0]>'z' ){
		fprintf(stderr,"checkparam: parameter string ");
		fprintf(stderr,"should start with 'a-z' or 'A-Z'!\n");
		exit(-2);
	}
	while ( param_dic[++i] != '\0' ){
		if ( (param_dic[i]>='A' && param_dic[i]<='Z') || \
				(param_dic[i]>='a' && param_dic[i]<='z') ){
			param_num++;
		} else if ( param_dic[i]<'0' || param_dic[i]>'9' ){
			fprintf(stderr,"checkparam: illegal parameter strings!\n");
			exit(-2);
		}
	}
	return param_num;
}

void getparam_dic(char *param_dic,
		char *option_str, int *option_num)
{
	int i=-1,j=-1;
	
	while ( param_dic[++i] != '\0' ){
		if ( (param_dic[i]>='A' && param_dic[i]<='Z') || \
				(param_dic[i]>='a' && param_dic[i]<='z') ){
			option_str[++j] = param_dic[i];
		} else {
			option_num[j] = atoi(param_dic+i);
			if ( option_num[j] > 15 ){
				fprintf(stderr,"every option should have no more than 16 parameters!\n");
				exit(-2);
			}
			while ( param_dic[i]>='0' && param_dic[i]<='9' ){
				i++;
			}
			i--;
		}
	}
	return;
}


int getopt(int argc, char *argv[],
		char *param_dic)
{
	static char *option_str=NULL;
	static int  *option_num=NULL;
	int i;
	static int firstrun=1;
	if ( firstrun ){
	int param_num=checkparam_dic(param_dic);
	option_str = (char*) malloc ( param_num*sizeof(char) );
	option_num = (int*)  malloc ( param_num*sizeof(int)  );
	for ( i=0;i<param_num;++i ) option_num[i]=0;
	getparam_dic(param_dic,option_str,option_num);
	firstrun=0;
	}
	
	static int arg_count=1;
//	arg_count++;
		// argument counter, skip 0 which is the command name
	char *arg_name;
		// argument(s) name, should be start with "-"
	static int arg_offset=1;
		// offset of arg_name[]
		// skip 1 which should be "-"
	int param_offset=0;
		// offset of param_dic[]
//	int param_flag=false;	// there is no bool type in standard C!
	int arg_flag=0;
		// whether the argument is allowed
	static int narg=0;
		// number of parameter(s) after argument
	int narg_temp;
	static int param_arg=0;
		// number of argument(s) with parameter(s)
		// within one "-"
		// it should be no large than 1

	int arg_len;	

	if ( arg_count >= argc ){
		free(option_str);
		free(option_num);
		return -1; // end of the arguments
	}

	arg_name = argv[arg_count];
	arg_len = strlen(arg_name);


	if ( arg_name[0] != '-' ){
		fprintf(stderr,"getopt: %s should start with \"-\"\n",arg_name);
		exit(1);
	} else {
	}

	if ( arg_name[0] != '-' ){
		fprintf(stderr,"getopt: %s should start with \"-\"\n",arg_name);
		exit(1);
	} else {
		param_offset = -1;
		if ( arg_name[arg_offset] != '\0' ){
			// decide the behavior of this argument which starts with "-"
//			while ( param_dic[++param_offset] != '\0' ){
			while ( option_str[++param_offset] != '\0' ){
				// whether this argument can be found in the param_dic
//				if ( arg_name[arg_offset] == param_dic[param_offset] ){
				if ( arg_name[arg_offset] == option_str[param_offset] ){
					arg_flag=1;
//					narg_temp = 0;
//					while ( param_dic[++param_offset] == ':' ){
//						narg_temp++;
//					}
					narg_temp = option_num[param_offset];
					// how many parameters/options should this parameter have
					if ( narg_temp>=1 ){
						narg = narg_temp;
						// check if the argument has reach its end
						if ( arg_count+narg >= argc ){
							fprintf(stderr,"not enough parameter(s) after \"-%c\"!\n",
									arg_name[arg_offset]);
							fprintf(stderr,
									"there should be %d parameters after \"-%c\"!\n",
								narg,arg_name[arg_offset]);
							exit(-2);
						}
						for ( i=0;i<narg;++i ){
							if ( argv[arg_count+i+1][0] == '-' ){
								fprintf(stderr,
										"getopt: parameter number %d error! ",i);
								fprintf(stderr,
										"there should be %d parameters after \"-%c\"!\n",
									narg,arg_name[arg_offset]);
								exit(-2);
							}
							opt_arr[i] = argv[arg_count+i+1];
						}
						opt_str = opt_arr[0];
						param_arg++;
					}
		//			opt_str = opt_arr[0];
					break;
				}
			}
			if ( !arg_flag ){
				fprintf(stderr,"getopt: illegal argument \"-%c\"\n",
						arg_name[arg_offset]);
				exit(-2);
			}
			if ( param_arg>1 ){
				fprintf(stderr,"getopt: too many arguments in \"%s\"!\n",
						arg_name);
				exit(-2);
			}
			if ( arg_offset == (arg_len-1) ){
				// reach the end of this argument
				arg_count += narg+1;
					// jump narg arguments (narg parameters/options)

				// re-initialize these static values
				narg=0;
				arg_offset=1;
				param_arg=0;
				
				return arg_name[arg_len-1];
			}
			return arg_name[arg_offset++];
		} else {
			fprintf(stderr,"getopt: no argument after \"-\"!\n");
			exit(-2);
		}
	}
//	return '?';
}
